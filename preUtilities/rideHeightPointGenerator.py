"""
Aero map DOE + 6x6 table builder for a VI-grade style FRH x RRH lookup.

Step 1 - choose CFD / tunnel runs from a VD trace:
    python rideHeightPointGenerator.py design  vd_export.csv  [--n 16] [--ref 50] [--ls 12,15]
    --ls: GP length scales (FRH,RRH mm) - how far a run's information reaches when choosing runs
    -> runs.csv (run list incl. validation runs), breakpoints.csv,
       design_points.png (runs over the trace), design_steps.png (how the runs were picked)

Step 1b - roll / yaw / steer correction sweeps:
    python rideHeightPointGenerator.py sweeps  vd_export.csv  [--nsweep 6] [--roll-bps 0,0.5,...] [--yaw-bps ...]
                                               [--steer-bps ...] [--rh1 dF,dR] [--rh2 dF,dR]
    One-sided breakpoints (car assumed left/right symmetric: mirror results for negative angles),
    spaced by V^2*(ay/ay_max)^2*dt (lateral-g weight; --sweep-weight grip for combined G) from 0 to the
    max seen. Sweeps run at the cornering ride height
    (rh1); mid + top levels repeated at a second, lower ride height (rh2) to test RH dependence;
    two combined corner runs (all angles with the sign of the corner) test additivity.
    -> sweeps.csv (run list; pass 1 = do first, pass 2 = fill in if the curve bends), sweep_breakpoints.csv,
       sweeps.png

Step 2 - after the runs, build the table:
    python rideHeightPointGenerator.py table   vd_export.csv  results.csv  [--nbp 6] [--ref 50] [--no-opt]
    results.csv columns: dFRH_mm, dRRH_mm, then one column per quantity
    (e.g. CzF, CzR, Cx). Rows with role == 'validation' (optional column) are held out.
    -> table_<quantity>.csv for each quantity (rows = RRH, cols = FRH, aero RH in mm),
       breakpoints.csv, table_<quantity>.png (fitted map + table error map), and a printed report.
    Add --no-plot to skip the PNGs.

Weighting (per quantity):
  downforce / balance -> V^2 * dt * (G/Gmax)^2, G = combined |ax,ay| in g: load x time x how
                         grip-limited the car is (corners, braking, traction), straights count little
  drag                -> V^2 * dt * V on straights (|ay| < 0.3 g, throttle > 90%): drag power
  Columns whose name starts with Cx / CD or contains 'drag' get the drag weight (or use --drag-cols).
  Run selection uses the grip weight; the lowest-RH (top speed) point is always included for drag.
Ride heights: aero_forces.height_front/rear, as deltas from --ref (mm).
"""
import sys, argparse, numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
from scipy.spatial import cKDTree

def warn(msg):
    print('WARNING: ' + msg, file=sys.stderr)

# ---------------------------------------------------------------- VD trace
def load_vd(path, ref):
    raw = pd.read_csv(path, sep=',', dtype=str).iloc[1:]          # row 0 = units
    g = lambda c: pd.to_numeric(raw[c], errors='coerce').values
    d = pd.DataFrame(dict(t=g('time.TIME'), hf=g('aero_forces.height_front'),
                          hr=g('aero_forces.height_rear'), V=g('aero_forces.longitudinal_velocity'),
                          ax=g('chassis_accelerations.longitudinal'), ay=g('chassis_accelerations.lateral'),
                          thr=g('driving_machine_monitor.driver_throttle')))
    d = d.dropna().drop_duplicates('t').reset_index(drop=True)
    d['dF'] = d.hf * 1000 - ref
    d['dR'] = d.hr * 1000 - ref
    q = d.V ** 2 * np.gradient(d.t)
    G = np.hypot(d.ax, d.ay)
    wg = q * (G / G.max()) ** 2
    straight = (d.ay.abs() < 0.3) & (d.thr > 90)
    wd = q * d.V * straight
    if wd.sum() <= 0:   # e.g. skidpad: no full-throttle straights
        warn('no straights (|ay| < 0.3 g and throttle > 90%) in the trace; '
             'drag weight falls back to V^3*dt over the whole trace')
        wd = q * d.V
    wl =q * (d.ay.abs() / d.ay.abs().max()) ** 2
    d['w_q'] = q / q.sum(); d['w'] = wg / wg.sum(); d['w_drag'] = wd / wd.sum(); d['w_lat'] = wl / wl.sum()
    return d

# ---------------------------------------------------------------- breakpoints
def weighted_breakpoints(x, w, n=6, power=1/3):
    """Interior spacing ~ (V^2-weighted density)^power; ends at the trace extremes."""
    h, e = np.histogram(x, bins=60, range=(x.min(), x.max()), weights=w)
    h = gaussian_filter1d(h, 2) ** power
    c = np.r_[0, np.cumsum(h)] / h.sum()
    return np.interp(np.linspace(0, 1, n), c, e)

# ---------------------------------------------------------------- step 1: run selection
def design(d, n=16, ls=(12., 15.), alpha=0.6, dmin=8.0):
    X = d[['dF', 'dR']].values; w = d.w.values
    anchors = [(0., 0.), tuple(X[np.argmin(X[:, 0] + X[:, 1])])]   # reference + deepest heave
    C = X[::2]; wc = w[::2] / w[::2].sum()
    m = alpha * wc + (1 - alpha) / len(C)
    L = np.array(ls)
    K = np.exp(-0.5 * (((C[:, None] - C[None]) / L) ** 2).sum(-1))
    S = [int(np.argmin(((C - np.array(a)) ** 2).sum(1))) for a in anchors]
    while len(S) < n:
        Kss = K[np.ix_(S, S)] + 1e-4 * np.eye(len(S)); Ksc = K[S]
        cov = K - Ksc.T @ np.linalg.solve(Kss, Ksc)
        gain = (m[:, None] * cov ** 2).sum(0) / np.maximum(np.diag(cov), 1e-9)
        dist = np.hypot(*(C[:, None] - C[S][None]).transpose(2, 0, 1)).min(1)
        gain[dist < dmin] = -1
        if gain.max() < 0:   # trace fully covered (e.g. skidpad): more runs would duplicate existing ones
            warn('only %d of %d runs placed: every trace point is within %g mm of a chosen run '
                 '(ride-height range too small for more runs)' % (len(S), n, dmin))
            break
        S.append(int(np.argmax(gain)))
    return np.round(C[S] * 2) / 2


def uncertainty(P, C, ls=(12., 15.)):
    """GP posterior std (0 = known, 1 = unknown) at C given runs at P. Distance only."""
    L = np.array(ls)
    k = lambda A, B: np.exp(-0.5 * (((A[:, None] - B[None]) / L) ** 2).sum(-1))
    Kpc = k(P, C)
    var = 1 - np.einsum('ij,ij->j', Kpc, np.linalg.solve(k(P, P) + 1e-4 * np.eye(len(P)), Kpc))
    return np.sqrt(np.clip(var, 0, 1))

def pick_validation(d, P, nval=2, ls=(12., 15.), dmin=8.0):
    """Hold-out runs where the car spends weighted time but the fit runs cover least:
    first on drag weight (straights), then on grip weight (corners)."""
    X = d[['dF', 'dR']].values; C = X[::2]; tree = cKDTree(X)
    sig = uncertainty(P, C, ls); V = []
    for wname in (['w_drag', 'w'] * nval)[:nval]:
        w = d[wname].values
        dens = np.array([w[tree.query_ball_point(c, 3)].sum() for c in C])
        score = dens * sig
        for j in np.argsort(-score):
            ok = np.min(np.hypot(*(P - C[j]).T)) > 4 and all(np.hypot(*(C[j] - v)) > dmin for v in V)
            if ok:
                V.append(C[j]); break
    if len(V) < nval:
        warn('only %d of %d validation runs placed: no trace point is > 4 mm from the fit runs '
             'and > %g mm from other validation runs' % (len(V), nval, dmin))
    return np.round(np.array(V).reshape(-1, 2) * 2) / 2


# ---------------------------------------------------------------- step 1b: roll / yaw / steer sweeps
def load_attitude(path):
    raw = pd.read_csv(path, sep=',', dtype=str).iloc[1:]
    g = lambda c: pd.to_numeric(raw[c], errors='coerce').values
    a = pd.DataFrame(dict(t=g('time.TIME'), roll=np.degrees(g('aero_forces.roll')),
                          yaw=np.degrees(g('aero_forces.side_slip')), steer=np.degrees(g('aero_forces.steer'))))
    off = {}
    for ax, sc, ac in (('F', 'rhl_front_setup_sensor_data.ride_height', 'aero_forces.height_front'),
                       ('R', 'rhs_rear_setup_sensor_data.ride_height', 'aero_forces.height_rear')):
        if sc in raw.columns:
            off[ax] = float(np.nanmedian(g(sc) - g(ac)) * 1000)       # setup sensor - aero height, mm
    a = a.dropna().drop_duplicates('t').reset_index(drop=True)
    return a, off

def _wq(x, w, p):
    o = np.argsort(x); c = np.cumsum(w[o]) / w.sum()
    return x[o][min(np.searchsorted(c, p), len(x) - 1)]

def sweep_breakpoints(x, w, n=6):
    """One-sided: 0 .. max|x|, interior spaced by grip-weighted density^(1/3) (tail beyond p99 clipped)."""
    x = np.abs(x); top = x.max(); p99 = _wq(x, w, .99)
    step = 0.1 if top <= 5 else (0.5 if top <= 20 else 1.0)
    inner = weighted_breakpoints(np.r_[0, np.minimum(x, p99)], np.r_[0, w], n)[:-1]
    bp = np.r_[np.round(inner / step) * step, np.ceil(top / step) * step]
    bp[0] = 0.0
    for i in range(1, len(bp)):                       # keep strictly increasing after rounding
        bp[i] = max(bp[i], bp[i - 1] + step)
    return np.round(bp, 3)

def density_peak(X, w, exclude=None, rmin=8.0, h=4.0):
    tree = cKDTree(X); C = X[::3]
    dens = np.array([w[tree.query_ball_point(c, h)].sum() for c in C])
    if exclude is not None:
        dist = np.hypot(*(C - exclude).T)
        if (dist >= rmin).any():
            dens[dist < rmin] = -1
        else:   # e.g. skidpad: whole cloud within rmin -> take the farthest point
            j = int(np.argmax(dist))
            warn('no cornering ride height is >= %g mm from rh1; rh2 set to the farthest one '
                 '(%.1f mm away), the RH check will be weak. Pass --rh2 to choose it.' % (rmin, dist[j]))
            return C[j]
    return C[np.argmax(dens)]

def sweeps(d, att, nsweep=6, bps=None, rh1=None, rh2=None, wcol='w_lat'):
    w = d[wcol].values; ay = d.ay.abs().values
    corner = ay > 1.0
    X = d[['dF', 'dR']].values
    if rh1 is None:
        rh1 = np.average(X[corner], axis=0, weights=w[corner])
    rh1 = np.round(np.asarray(rh1, float) * 2) / 2
    if rh2 is None:   # busiest cornering location at least 8 mm away from rh1
        rh2 = density_peak(X[corner], w[corner], exclude=rh1)
    rh2 = np.round(np.asarray(rh2, float) * 2) / 2
    B = {}
    for k in ('roll', 'yaw', 'steer'):
        B[k] = np.array(bps[k]) if bps and bps.get(k) is not None else sweep_breakpoints(att[k].values, w, nsweep)
    # combined corner states: median attitude in corners, and in the hardest quarter of corners
    hard = corner & (ay > _wq(ay[corner], w[corner], .75))
    comb = {}
    for lab, m in (('moderate corner', corner), ('hard corner', hard)):
        comb[lab] = {k: round(float(_wq(np.abs(att[k].values[m]), w[m], .5)), 1) for k in ('roll', 'yaw', 'steer')}
    rows = []
    add = lambda role, rh, r, y, s, p: rows.append(dict(role=role, dFRH_mm=rh[0], dRRH_mm=rh[1], roll_deg=r,
                                                        yaw_deg=y, steer_deg=s, pass_=p))
    add('base (0 deg)', rh1, 0, 0, 0, 1)
    for k in ('roll', 'yaw', 'steer'):
        lv = B[k][1:]; first = {len(lv) - 1, len(lv) // 2}        # top + middle level first
        for i, v in enumerate(lv):
            add('sweep %s' % k, rh1, *(v if k == q else 0 for q in ('roll', 'yaw', 'steer')), 1 if i in first else 2)
    add('base (0 deg)', rh2, 0, 0, 0, 1)
    for k in ('roll', 'yaw', 'steer'):
        lv = B[k][1:]
        for v in (lv[len(lv) // 2], lv[-1]):
            add('RH check %s' % k, rh2, *(v if k == q else 0 for q in ('roll', 'yaw', 'steer')), 1)
    for lab, c in comb.items():
        add('combined, %s' % lab, rh1, c['roll'], c['yaw'], c['steer'], 1)
    S = pd.DataFrame(rows).rename(columns={'pass_': 'pass'})
    S.insert(0, 'run', range(1, len(S) + 1))
    return S, B, rh1, rh2, comb

def plot_sweeps(d, att, B, rh1, rh2, ref, fname='sweeps.png', wcol='w_lat'):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    w = d[wcol].values; corner = d.ay.abs().values > 1.0
    wlab = {'w_lat': 'lateral-g-weighted', 'w': 'grip-weighted'}[wcol]
    fig, axs = plt.subplots(1, 4, figsize=(19, 4.6), dpi=140, gridspec_kw=dict(width_ratios=[1.15, 1, 1, 1]))
    ax = axs[0]
    o = np.argsort(w)
    ax.scatter(d.dF.values[o], d.dR.values[o], c=w[o] / w.max(), cmap='Blues', s=5, lw=0, vmin=0,
               vmax=np.quantile(w / w.max(), .99))
    for p, lab in ((rh1, 'Sweeps (cornering RH)'), (rh2, 'RH check (lower)')):
        ax.scatter(*p, s=150, c='#c0392b' if p is rh1 else 'white', ec='#c0392b', lw=2, zorder=5,
                   marker='o' if p is rh1 else 's', label=lab)
    ax.set_title('Sweep ride heights', loc='left', fontsize=11); _ax_style(ax, ref); _limits(ax, d, 4)
    ax.legend(loc='upper right', frameon=False, fontsize=8.5)
    units = dict(roll='Roll', yaw='Yaw / sideslip', steer='Steer')
    for ax, k in zip(axs[1:], ('roll', 'yaw', 'steer')):
        x = np.abs(att[k].values); top = max(B[k][-1], x.max())
        ax.hist(x, bins=40, range=(0, top), weights=w / w.sum() * 100, color='#6baed6',
                edgecolor='white', lw=.5)
        for i, b in enumerate(B[k]):
            ax.axvline(b, c='#c0392b', lw=1.4)
            ax.annotate('%g' % b, (b, 1), xycoords=('data', 'axes fraction'), xytext=(2, -12),
                        textcoords='offset points', fontsize=8.5, color='#c0392b')
        ax.set_title('%s: %d one-sided breakpoints' % (units[k], len(B[k])), loc='left', fontsize=11)
        ax.set_xlabel('|%s| (deg)' % k)
        ax.set_ylabel('%% of %s time' % wlab)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    fig.tight_layout(); fig.savefig(fname); plt.close(fig)

# ---------------------------------------------------------------- step 2: surrogate
class Surrogate:
    """Quadratic trend + anisotropic GP on the residuals (trend keeps extrapolation to
    the empty table corners sensible instead of decaying to the mean)."""
    def __init__(self, X, y):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
        self.mu, self.sd = X.mean(0), X.std(0)
        Z = (X - self.mu) / self.sd
        A = self._basis(Z)
        self.beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A @ self.beta
        self.rs = r.std() if r.std() > 0 else 1.0
        k = ConstantKernel(1.0, (1e-3, 1e2)) * RBF([1.0, 1.0], (0.05, 20)) + WhiteKernel(1e-3, (1e-8, 1e-1))
        self.gp = GaussianProcessRegressor(k, normalize_y=False, n_restarts_optimizer=5, random_state=0)
        self.gp.fit(Z, r / self.rs)
    @staticmethod
    def _basis(Z):
        f, r = Z[:, 0], Z[:, 1]
        return np.c_[np.ones(len(Z)), f, r, f * f, f * r, r * r]
    def __call__(self, X):
        Z = (np.atleast_2d(X) - self.mu) / self.sd
        return self._basis(Z) @ self.beta + self.gp.predict(Z) * self.rs

def loo_error(X, y):
    e = []
    for i in range(len(X)):
        m = np.arange(len(X)) != i
        e.append(Surrogate(X[m], y[m])(X[i:i + 1])[0] - y[i])
    return np.array(e)

# ---------------------------------------------------------------- step 2: table
def fill_table(models, bf, br):
    G = np.array(np.meshgrid(bf, br, indexing='ij')).reshape(2, -1).T
    return {k: m(G).reshape(len(bf), len(br)) for k, m in models.items()}

def table_error(models, scales, bf, br, Xq, wq, truth):
    # wq: dict quantity -> weights
    """Weighted RMS of (bilinear table - surrogate) at the VD samples, summed over
    quantities after normalising each by its spread."""
    T = fill_table(models, bf, br); tot = 0; per = {}
    for k in models:
        f = RegularGridInterpolator((bf, br), T[k], method='linear', bounds_error=False, fill_value=None)
        e = f(Xq) - truth[k]
        per[k] = np.sqrt((wq[k] * e ** 2).sum())
        tot += (per[k] / scales[k]) ** 2
    return np.sqrt(tot), per

def optimise_bps(models, scales, bf0, br0, Xq, wq, truth):
    lo_f, hi_f, lo_r, hi_r = bf0[0], bf0[-1], br0[0], br0[-1]
    nf, nr = len(bf0) - 2, len(br0) - 2
    # parametrise interior breakpoints by positive gaps -> ordering always holds
    def unpack(p):
        gf = np.exp(p[:nf + 1]); gr = np.exp(p[nf + 1:])
        bf = lo_f + (hi_f - lo_f) * np.r_[0, np.cumsum(gf)] / gf.sum()
        br = lo_r + (hi_r - lo_r) * np.r_[0, np.cumsum(gr)] / gr.sum()
        return bf, br
    p0 = np.r_[np.log(np.diff(bf0)), np.log(np.diff(br0))]
    obj = lambda p: table_error(models, scales, *unpack(p), Xq, wq, truth)[0]
    res = minimize(obj, p0, method='Nelder-Mead', options=dict(maxiter=3000, xatol=1e-3, fatol=1e-7))
    bf, br = unpack(res.x)
    return np.round(bf * 2) / 2, np.round(br * 2) / 2


# ---------------------------------------------------------------- plots
def _ax_style(ax, ref):
    ax.set_xlabel('dFRH (mm, from %g mm aero RH)' % ref); ax.set_ylabel('dRRH (mm, from %g mm aero RH)' % ref)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    ax.set_aspect('equal')

def _limits(ax, d, pad=6):
    ax.set_xlim(d.dF.min() - pad, d.dF.max() + pad); ax.set_ylim(d.dR.min() - pad, d.dR.max() + pad)

def _grid(ax, bf, br, label=True):
    for x in bf: ax.axvline(x, c='#999', lw=.8, ls=':')
    for y in br: ax.axhline(y, c='#999', lw=.8, ls=':')
    if label: ax.plot([], [], c='#999', lw=.8, ls=':', label='%dx%d breakpoints' % (len(bf), len(br)))

def _runs(ax, fit, val, numbers=True):
    ax.scatter(fit[:, 0], fit[:, 1], s=120, c='#c0392b', ec='white', lw=2, zorder=5, label='Fit runs (run order)')
    if len(val):
        ax.scatter(val[:, 0], val[:, 1], s=105, marker='s', facecolor='white', ec='#c0392b', lw=2,
                   zorder=5, label='Validation (hold out)')
    if numbers:
        for i, p in enumerate(np.vstack([fit, val]) if len(val) else fit):
            ax.annotate(str(i + 1), p, xytext=(6, 5), textcoords='offset points', fontsize=8.5, color='#222', zorder=6)

def plot_design(d, fit, val, bf, br, ref, fname='design_points.png'):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
    o = np.argsort(d.w.values); wr = d.w.values / d.w.max()
    sc = ax.scatter(d.dF.values[o], d.dR.values[o], c=wr[o], cmap='Blues', s=9, alpha=.7, lw=0, vmin=0,
                    vmax=np.quantile(wr, .99))
    cb = fig.colorbar(sc, ax=ax, pad=.01, shrink=.85); cb.set_label('Grip weight  V^2 x (G/Gmax)^2  (relative)')
    _grid(ax, bf, br); _runs(ax, fit, val)
    _ax_style(ax, ref); _limits(ax, d)
    ax.set_title('Aero map runs: %d fit + %d validation' % (len(fit), len(val)), loc='left', fontsize=13)
    ax.legend(loc='upper right', frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(fname); plt.close(fig)

def plot_steps(d, fit, ref, ls=(12., 15.), fname='design_steps.png'):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    C = d[['dF', 'dR']].values[::2]; n = len(fit)
    stages = sorted({2, min(4, n), min(8, n), n})
    fig, axs = plt.subplots(1, len(stages), figsize=(4 * len(stages), 4.6), dpi=140, sharey=True)
    axs = np.atleast_1d(axs)
    for ax, k in zip(axs, stages):
        sig = uncertainty(fit[:k], C, ls)
        sc = ax.scatter(C[:, 0], C[:, 1], c=sig, cmap='Blues', vmin=0, vmax=1, s=5, lw=0)
        ax.scatter(fit[:k, 0], fit[:k, 1], s=60, c='#c0392b', ec='white', lw=1.5, zorder=5)
        title = '%d runs chosen' % k
        if k < n:
            ax.scatter(*fit[k], s=150, facecolor='none', ec='#222', lw=1.8, zorder=6); title += '  (o = next pick)'
        ax.set_title(title, loc='left', fontsize=11); _ax_style(ax, ref); _limits(ax, d, 4)
        if ax is not axs[0]: ax.set_ylabel('')
    cb = fig.colorbar(sc, ax=list(axs), pad=.01, fraction=.015); cb.set_label('Uncertainty (1 = unknown)')
    fig.suptitle('Run selection: uncertainty over the trace shrinks as runs are added', x=.06, ha='left', fontsize=13)
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)

def plot_table(k, model, table, d, fit, val, bf, br, w, ref, fname):
    """Left: fitted map with trace, runs and breakpoints. Right: table - fitted map at every point."""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fx = np.linspace(bf[0], bf[-1], 120); fy = np.linspace(br[0], br[-1], 120)
    G = np.array(np.meshgrid(fx, fy, indexing='ij')).reshape(2, -1).T
    Z = model(G).reshape(len(fx), len(fy))
    Zt = RegularGridInterpolator((bf, br), table, method='linear')(G).reshape(len(fx), len(fy))
    E = Zt - Z; lim = np.abs(E).max() or 1.0
    Xq = d[['dF', 'dR']].values; eq = RegularGridInterpolator((bf, br), table, method='linear',
                                                            bounds_error=False, fill_value=None)(Xq) - model(Xq)
    rms = np.sqrt((w * eq ** 2).sum() / w.sum())
    fig, axs = plt.subplots(1, 2, figsize=(16, 7), dpi=140)
    cs = axs[0].contourf(fx, fy, Z.T, 20, cmap='Blues')
    axs[0].contour(fx, fy, Z.T, 10, colors='white', linewidths=.5, alpha=.6)
    fig.colorbar(cs, ax=axs[0], pad=.01, shrink=.85).set_label(k)
    axs[0].scatter(d.dF, d.dR, s=2, c='#333', alpha=.15, lw=0, label='VD trace')
    _grid(axs[0], bf, br); _runs(axs[0], fit, val)
    axs[0].set_title('%s: fitted map' % k, loc='left', fontsize=12)
    cs = axs[1].pcolormesh(fx, fy, E.T, cmap='RdBu_r', vmin=-lim, vmax=lim, shading='auto')
    fig.colorbar(cs, ax=axs[1], pad=.01, shrink=.85).set_label('%s table - fitted map' % k)
    axs[1].scatter(d.dF, d.dR, s=2, c='#333', alpha=.15, lw=0)
    _grid(axs[1], bf, br, label=False); _runs(axs[1], fit, val, numbers=False)
    axs[1].set_title('%s: %dx%d table error (weighted RMS on trace %.3g)' % (k, len(bf), len(br), rms),
                     loc='left', fontsize=12)
    for ax in axs:
        _ax_style(ax, ref); ax.set_xlim(bf[0] - 3, bf[-1] + 3); ax.set_ylim(br[0] - 3, br[-1] + 3)
    leg = axs[0].legend(loc='upper right', frameon=False, fontsize=8.5)
    for h in leg.legend_handles[:1]: h.set_sizes([20]); h.set_alpha(.6)
    fig.tight_layout(); fig.savefig(fname); plt.close(fig)

# ---------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['design', 'sweeps', 'table'])
    ap.add_argument('vd'); ap.add_argument('results', nargs='?')
    ap.add_argument('--n', type=int, default=16); ap.add_argument('--nbp', type=int, default=6)
    ap.add_argument('--ref', type=float, default=50.0); ap.add_argument('--no-opt', action='store_true')
    ap.add_argument('--drag-cols', default='')
    ap.add_argument('--nval', type=int, default=2); ap.add_argument('--no-plot', action='store_true')
    ap.add_argument('--nsweep', type=int, default=6)
    for k in ('roll', 'yaw', 'steer'): ap.add_argument('--%s-bps' % k, default='')
    ap.add_argument('--rh1', default=''); ap.add_argument('--rh2', default='')
    ap.add_argument('--sweep-weight', choices=['lat', 'grip'], default='lat',
                    help="sweeps: lat = V^2*(ay/ay_max)^2*dt (default), grip = V^2*(G/Gmax)^2*dt")
    ap.add_argument('--ls', default='12,15',
                    help='design: GP length scales FRH,RRH in mm for run / validation selection (default 12,15)')
    a = ap.parse_args()
    ls = tuple(float(v) for v in a.ls.split(','))
    if len(ls) != 2 or min(ls) <= 0:
        raise SystemExit('--ls needs two positive values: FRH,RRH (mm)')
    d = load_vd(a.vd, a.ref)
    bf0 = weighted_breakpoints(d.dF.values, d.w.values, a.nbp)
    br0 = weighted_breakpoints(d.dR.values, d.w.values, a.nbp)
    if a.mode != 'sweeps':   # rounded breakpoints must stay strictly increasing
        for lab, b in (('FRH', bf0), ('RRH', br0)):
            if np.unique(np.round(b * 2) / 2).size < a.nbp:
                raise SystemExit('%s range in the trace (%.1f mm) is too small for %d breakpoints at 0.5 mm '
                                 'resolution; lower --nbp' % (lab, b[-1] - b[0], a.nbp))

    if a.mode == 'sweeps':
        att, off = load_attitude(a.vd)
        if len(att) != len(d):
            raise SystemExit('attitude and ride-height channels have different lengths')
        num = lambda t: [float(v) for v in t.split(',')] if t else None
        bps = {k: num(getattr(a, '%s_bps' % k)) for k in ('roll', 'yaw', 'steer')}
        wcol = {'lat': 'w_lat', 'grip': 'w'}[a.sweep_weight]
        S, B, rh1, rh2, comb = sweeps(d, att, a.nsweep, bps, num(a.rh1), num(a.rh2), wcol)
        S['FRH_aero_mm'] = S.dFRH_mm + a.ref; S['RRH_aero_mm'] = S.dRRH_mm + a.ref
        if off:
            S['FRH_setup_mm'] = (S.FRH_aero_mm + off.get('F', np.nan)).round(1)
            S['RRH_setup_mm'] = (S.RRH_aero_mm + off.get('R', np.nan)).round(1)
        S.to_csv('sweeps.csv', index=False); print(S.to_string(index=False))
        L = max(len(v) for v in B.values())
        pd.DataFrame({k + '_deg': np.r_[v, [np.nan] * (L - len(v))] for k, v in B.items()}).to_csv(
            'sweep_breakpoints.csv', index=False)
        print('\nOne-sided breakpoints (deg, mirror for negative angles):')
        for k, v in B.items(): print('  %-5s %s' % (k, ', '.join('%g' % x for x in v)))
        print('Sweep ride height (dF, dR) = (%g, %g); RH check at (%g, %g)' % (*rh1, *rh2))
        print('Signs: give roll, yaw and steer the sign they have in the same corner (map to your CFD convention).')
        print('%d runs: %d in pass 1, %d in pass 2' % (len(S), (S['pass'] == 1).sum(), (S['pass'] == 2).sum()))
        if not a.no_plot:
            plot_sweeps(d, att, B, rh1, rh2, a.ref, wcol=wcol); print('Wrote sweeps.png')
        return

    if a.mode == 'design':
        P = design(d, a.n, ls=ls)
        Vp = pick_validation(d, P, a.nval, ls=ls) if a.nval > 0 else np.zeros((0, 2))
        A = np.vstack([P, Vp]) if len(Vp) else P
        out = pd.DataFrame(dict(run=range(1, len(A) + 1), role=['fit'] * len(P) + ['validation'] * len(Vp),
                                dFRH_mm=A[:, 0], dRRH_mm=A[:, 1],
                                FRH_aero_mm=A[:, 0] + a.ref, RRH_aero_mm=A[:, 1] + a.ref))
        out.to_csv('runs.csv', index=False); print(out.to_string(index=False))
        bf, br = np.round(bf0 * 2) / 2, np.round(br0 * 2) / 2
        if not a.no_plot:
            plot_design(d, P, Vp, bf, br, a.ref); plot_steps(d, P, a.ref, ls=ls)
            print('Wrote design_points.png, design_steps.png')
    else:
        R = pd.read_csv(a.results)
        val = R.get('role', pd.Series([''] * len(R))).astype(str).str.contains('valid')
        qcols = [c for c in R.columns if c not in ('dFRH_mm', 'dRRH_mm', 'role', 'run')
                 and np.issubdtype(R[c].dtype, np.number)]
        Xf = R.loc[~val, ['dFRH_mm', 'dRRH_mm']].values
        models, scales = {}, {}
        print('Surrogate check (fit on %d runs):' % len(Xf))
        for k in qcols:
            y = R.loc[~val, k].values
            models[k] = Surrogate(Xf, y); scales[k] = np.ptp(y) or 1.0
            loo = loo_error(Xf, y)
            msg = '  %-10s LOO RMS %.4g (%.1f%% of range)' % (k, np.sqrt((loo ** 2).mean()), 100 * np.sqrt((loo ** 2).mean()) / scales[k])
            if val.any():
                ev = models[k](R.loc[val, ['dFRH_mm', 'dRRH_mm']].values) - R.loc[val, k].values
                msg += ' | validation errors ' + ', '.join('%.4g' % x for x in ev)
            print(msg)
        Xq = d[['dF', 'dR']].values[::2]
        dragc = set(a.drag_cols.split(',')) if a.drag_cols else {k for k in qcols
                  if k.lower().startswith(('cx', 'cd')) or 'drag' in k.lower()}
        wq = {}
        for k in qcols:
            ww = (d.w_drag if k in dragc else d.w).values[::2]; wq[k] = ww / ww.sum()
        print('Weighting: drag-weighted ->', sorted(dragc) or 'none', '| grip-weighted -> the rest')
        truth = {k: models[k](Xq) for k in models}
        bf, br = np.round(bf0 * 2) / 2, np.round(br0 * 2) / 2
        e0, per0 = table_error(models, scales, bf, br, Xq, wq, truth)
        print('Table vs surrogate, weighted RMS, weighted breakpoints: ' +
              ', '.join('%s %.4g' % kv for kv in per0.items()))
        if not a.no_opt:
            bfo, bro = optimise_bps(models, scales, bf, br, Xq, wq, truth)
            if np.all(np.diff(bfo) > 0) and np.all(np.diff(bro) > 0):
                bf, br = bfo, bro
            else:   # optimiser squeezed a gap below the 0.5 mm rounding
                warn('optimised breakpoints collapse after 0.5 mm rounding; keeping the weighted breakpoints')
            e1, per1 = table_error(models, scales, bf, br, Xq, wq, truth)
            print('Table vs surrogate, weighted RMS, optimised breakpoints: ' +
                  ', '.join('%s %.4g' % kv for kv in per1.items()) + '  (%.0f%% lower)' % (100 * (1 - e1 / e0)))
        T = fill_table(models, bf, br)
        for k, M in T.items():
            df = pd.DataFrame(M.T, index=br + a.ref, columns=bf + a.ref)   # rows RRH, cols FRH (aero mm)
            df.index.name = 'RRH_mm \\ FRH_mm'
            df.round(5).to_csv('table_%s.csv' % k)
        if not a.no_plot:
            fitP = R.loc[~val, ['dFRH_mm', 'dRRH_mm']].values; valP = R.loc[val, ['dFRH_mm', 'dRRH_mm']].values
            for k in models:
                wk = (d.w_drag if k in dragc else d.w).values
                plot_table(k, models[k], T[k], d, fitP, valP, bf, br, wk, a.ref, 'table_%s.png' % k)
            print('Wrote', ', '.join('table_%s.png' % k for k in models))
    pd.DataFrame({'FRH_aero_mm': bf + a.ref, 'RRH_aero_mm': br + a.ref,
                  'dFRH_mm': bf, 'dRRH_mm': br}).to_csv('breakpoints.csv', index=False)
    print('Breakpoints FRH (aero mm):', bf + a.ref)
    print('Breakpoints RRH (aero mm):', br + a.ref)

if __name__ == '__main__':
    main()