"""
Figures for the "Aero Map Run Selection and Table Building" theory section.

    python makeFigures.py [vd_export.csv]

Uses preUtilities/rideHeightPointGenerator.py on a VI-grade lap trace (default: preUtilities/F27_mxpDamped.csv)
and writes the PNGs next to this file. With NORMALISE = True every axis that comes from the trace is scaled to
0-1 (ride heights over the trace range, angles over their maximum), so no real values can be read off the figures.
Figures 7 and 8 use a synthetic aero response, not CFD results.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NORMALISE = True

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'preUtilities'))
import rideHeightPointGenerator as rh  # noqa: E402

TRACE = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'preUtilities' / 'F27_mxpDamped.csv'
REF = 50.0
LS = (12., 15.)
RED, BLUE = '#c0392b', '#2c7fb8'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 200, 'savefig.bbox': 'tight'})

d = rh.load_vd(str(TRACE), REF)
LO = np.array([d.dF.min(), d.dR.min()]); SPAN = np.array([np.ptp(d.dF), np.ptp(d.dR)])


def n(X):
    """Ride heights (..., 2) in mm -> plotted units."""
    return (np.asarray(X, float) - LO) / SPAN if NORMALISE else np.asarray(X, float)


def nf(x): return n(np.c_[x, np.zeros_like(x)])[:, 0]
def nr(x): return n(np.c_[np.zeros_like(x), x])[:, 1]


def rh_axes(ax, pad=0.04):
    if NORMALISE:
        ax.set_xlabel('normalised FRH'); ax.set_ylabel('normalised RRH')
        ax.set_xlim(-pad, 1 + pad); ax.set_ylim(-pad, 1 + pad)
    else:
        ax.set_xlabel('dFRH (mm)'); ax.set_ylabel('dRRH (mm)')
    ax.set_aspect('equal')


def cloud(ax, w, cmap='Blues', s=3):
    o = np.argsort(w); P = n(d[['dF', 'dR']].values)
    wr = w / w.max()
    return ax.scatter(P[o, 0], P[o, 1], c=wr[o], cmap=cmap, s=s, lw=0, vmin=0, vmax=np.quantile(wr, .99))


def save(fig, name):
    fig.savefig(HERE / name); plt.close(fig); print('wrote', name)


# ------------------------------------------------------------------ 1: trace weighting
fig, axs = plt.subplots(1, 2, figsize=(7.5, 3.6))
for ax, col, title in ((axs[0], 'w', 'Grip weight  $V^2\\Delta t\\,(G/G_{max})^2$'),
                       (axs[1], 'w_drag', 'Drag weight  $V^3\\Delta t$ on straights')):
    sc = cloud(ax, d[col].values)
    ax.set_title(title, loc='left'); rh_axes(ax)
fig.colorbar(sc, ax=axs, shrink=.8, pad=.02).set_label('relative weight')
save(fig, 'rhWeights.png')

# ------------------------------------------------------------------ 2: breakpoint spacing exponent
x = nf(d.dF.values); w = d.w.values
fig, (a0, a1) = plt.subplots(2, 1, figsize=(6.5, 3.8), sharex=True, gridspec_kw=dict(height_ratios=[2, 1.2]))
a0.hist(x, bins=60, weights=w * 100, color='#9ecae1', edgecolor='white', lw=.4)
a0.set_ylabel('% of grip weight'); a0.set_title('Front ride-height breakpoints for different exponents', loc='left')
for i, (p, lab) in enumerate(((0, 'power 0 (uniform)'), (1 / 3, 'power 1/3 (used)'), (1, 'power 1 (follow density)'))):
    bp = rh.weighted_breakpoints(x, w, 6, power=p) if p > 0 else np.linspace(x.min(), x.max(), 6)
    a1.plot(bp, np.full_like(bp, -i), '|', ms=14, mew=2, color=RED if p == 1 / 3 else '#555')
    a1.text(1.02, -i, lab, va='center', fontsize=8, transform=a1.get_yaxis_transform())
a1.set_ylim(-2.6, .6); a1.set_yticks([]); a1.spines['left'].set_visible(False)
a1.set_xlabel('normalised FRH' if NORMALISE else 'dFRH (mm)')
save(fig, 'rhBreakpointPower.png')

# ------------------------------------------------------------------ 3: GP length scale in 1D
xs = np.linspace(0, 1, 400); xp = np.array([.15, .4, .55, .85])
f = lambda t: .5 * np.sin(8 * t) + t
yp = f(xp) - f(xp).mean()
fig, axs = plt.subplots(1, 2, figsize=(7.5, 2.8), sharey=True)
for ax, l, lab in ((axs[0], .1, 'short length scale'), (axs[1], .25, 'long length scale')):
    k = lambda A, B: np.exp(-.5 * ((A[:, None] - B[None]) / l) ** 2)
    Kpp = k(xp, xp) + 1e-6 * np.eye(len(xp)); Ksp = k(xs, xp)
    mu = Ksp @ np.linalg.solve(Kpp, yp)
    sd = np.sqrt(np.clip(1 - np.einsum('ij,ji->i', Ksp, np.linalg.solve(Kpp, Ksp.T)), 0, 1))
    ax.fill_between(xs, mu - 2 * sd, mu + 2 * sd, color='#c6dbef', label='$\\pm 2\\sigma$')
    ax.plot(xs, mu, c=BLUE, lw=1.5, label='posterior mean')
    ax.plot(xp, yp, 'o', c=RED, ms=6, label='runs')
    ax.annotate('', xy=(.4 + l, -1.6), xytext=(.4, -1.6), arrowprops=dict(arrowstyle='<->', color='#333'))
    ax.text(.4 + l / 2, -1.5, '$\\ell$', ha='center', fontsize=10)
    ax.set_title(lab, loc='left'); ax.set_xlabel('ride height (normalised)')
axs[0].set_ylabel('quantity'); axs[0].set_ylim(-1.9, 1.6); axs[1].legend(frameon=False, fontsize=8, loc='upper left')
save(fig, 'rhLengthScale.png')

# ------------------------------------------------------------------ 4: greedy run selection stages
P = rh.design(d, 16, ls=LS)
C = d[['dF', 'dR']].values[::2]; Cn = n(C)
stages = [2, 4, 8, 16]
fig, axs = plt.subplots(1, 4, figsize=(11, 3.2), sharey=True)
for ax, kk in zip(axs, stages):
    sig = rh.uncertainty(P[:kk], C, LS)
    sc = ax.scatter(Cn[:, 0], Cn[:, 1], c=sig, cmap='Blues', vmin=0, vmax=1, s=3, lw=0)
    Pn = n(P[:kk]); ax.scatter(Pn[:, 0], Pn[:, 1], s=40, c=RED, ec='white', lw=1, zorder=5)
    t = '%d runs' % kk
    if kk < len(P):
        nx = n(P[kk]); ax.scatter(*nx, s=110, facecolor='none', ec='#222', lw=1.5, zorder=6); t += ' (o = next)'
    ax.set_title(t, loc='left'); rh_axes(ax)
    if ax is not axs[0]: ax.set_ylabel('')
fig.colorbar(sc, ax=list(axs), pad=.01, fraction=.015).set_label('$\\sigma$ (1 = unknown)')
save(fig, 'rhSelectionSteps.png')

# ------------------------------------------------------------------ 5: effect of --ls on the layout
fig, axs = plt.subplots(1, 3, figsize=(10, 3.4), sharey=True)
for ax, m in zip(axs, (.5, 1., 2.)):
    ls = (LS[0] * m, LS[1] * m)
    Pm = rh.design(d, 16, ls=ls)
    sig = rh.uncertainty(Pm, C, ls)
    sc = ax.scatter(Cn[:, 0], Cn[:, 1], c=sig, cmap='Blues', vmin=0, vmax=1, s=3, lw=0)
    Pn = n(Pm); ax.scatter(Pn[:, 0], Pn[:, 1], s=40, c=RED, ec='white', lw=1, zorder=5)
    ax.set_title('$\\ell$ = %g x default' % m, loc='left'); rh_axes(ax)
    if ax is not axs[0]: ax.set_ylabel('')
fig.colorbar(sc, ax=list(axs), pad=.01, fraction=.015).set_label('$\\sigma$ after 16 runs')
save(fig, 'rhLengthScaleLayout.png')

# ------------------------------------------------------------------ 6: sweeps
att, _ = rh.load_attitude(str(TRACE))
if len(att) == len(d):
    S, B, rh1, rh2, comb = rh.sweeps(d, att, 6)
    wl = d.w_lat.values
    fig, axs = plt.subplots(1, 4, figsize=(12, 3), gridspec_kw=dict(width_ratios=[1.1, 1, 1, 1], wspace=.3))
    cloud(axs[0], wl)
    for p, mk, fc, lab in ((rh1, 'o', RED, 'rh1 (sweeps)'), (rh2, 's', 'white', 'rh2 (RH check)')):
        axs[0].scatter(*n(p), s=90, marker=mk, c=fc, ec=RED, lw=1.8, zorder=5, label=lab)
    axs[0].legend(frameon=False, fontsize=7.5, loc='upper left'); axs[0].set_title('Sweep ride heights', loc='left')
    rh_axes(axs[0])
    for ax, k, lab in zip(axs[1:], ('roll', 'yaw', 'steer'), ('roll', 'yaw (side slip)', 'steer')):
        xa = np.abs(att[k].values); top = max(B[k][-1], xa.max())
        sc_ = top if NORMALISE else 1.0
        ax.hist(xa / sc_, bins=40, range=(0, top / sc_), weights=wl * 100, color='#9ecae1', edgecolor='white', lw=.4)
        for b in B[k] / sc_:
            ax.axvline(b, c=RED, lw=1.2)
        ax.set_title('|%s|' % lab, loc='left')
        ax.set_xlabel('|%s| / max' % k if NORMALISE else '|%s| (deg)' % k)
        if ax is axs[1]: ax.set_ylabel('% of lateral weight')
    save(fig, 'rhSweeps.png')
else:
    print('skipped rhSweeps.png: attitude and ride-height rows differ')

# ------------------------------------------------------------------ 7: quadratic trend + GP vs plain GP (synthetic, 1D)
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
g = lambda t: 1.2 - 2.5 * (t - .45) ** 2 + .03 * np.sin(14 * t)   # mostly quadratic, small local detail
xr = np.linspace(.3, .7, 7); yr = g(xr); xs = np.linspace(0, 1, 300)
kern = lambda: ConstantKernel(1., (1e-3, 1e2)) * RBF(1., (.05, 20)) + WhiteKernel(1e-3, (1e-8, 1e-1))
mu, sd = xr.mean(), xr.std(); z, zs = (xr - mu) / sd, (xs - mu) / sd
plain = GaussianProcessRegressor(kern(), normalize_y=True, n_restarts_optimizer=5, random_state=0).fit(z[:, None], yr)
A = np.c_[np.ones_like(z), z, z * z]; beta = np.linalg.lstsq(A, yr, rcond=None)[0]; res = yr - A @ beta
rs = res.std() or 1.0
gp = GaussianProcessRegressor(kern(), n_restarts_optimizer=5, random_state=0).fit(z[:, None], res / rs)
trend = np.c_[np.ones_like(zs), zs, zs * zs] @ beta
fig, ax = plt.subplots(figsize=(6.5, 3))
ax.axvspan(xr.min(), xr.max(), color='#f0f0f0', zorder=0, label='range covered by runs')
ax.plot(xs, g(xs), c='#999', lw=1, ls='--', label='true response (synthetic)')
ax.plot(xs, plain.predict(zs[:, None]), c='#555', lw=1.5, label='GP only')
ax.plot(xs, trend + gp.predict(zs[:, None]) * rs, c=BLUE, lw=1.8, label='quadratic trend + GP (used)')
ax.plot(xr, yr, 'o', c=RED, ms=6, zorder=5, label='runs')
ax.set_xlabel('ride height (normalised)'); ax.set_ylabel('quantity')
ax.set_title('Extrapolating to the table corners', loc='left'); ax.legend(frameon=False, fontsize=7.5, loc='lower center')
save(fig, 'rhSurrogateTrend.png')

# ------------------------------------------------------------------ 8: table error, weighted vs optimised breakpoints (synthetic CzR)
from scipy.interpolate import RegularGridInterpolator
czr = lambda X: -1.6 + .003 * X[:, 0] - .010 * X[:, 1] - .0002 * X[:, 1] ** 2 + .08 * np.tanh((X[:, 1] + 30) / 6)
model = rh.Surrogate(P, czr(P)); scale = np.ptp(czr(P))
bf = np.round(rh.weighted_breakpoints(d.dF.values, d.w.values, 6) * 2) / 2
br = np.round(rh.weighted_breakpoints(d.dR.values, d.w.values, 6) * 2) / 2
Xq = d[['dF', 'dR']].values[::2]; wq = d.w.values[::2] / d.w.values[::2].sum()
args = ({'CzR': model}, {'CzR': scale})
bfo, bro = rh.optimise_bps(*args, bf, br, Xq, {'CzR': wq}, {'CzR': model(Xq)})
fx = np.linspace(bf[0], bf[-1], 150); fy = np.linspace(br[0], br[-1], 150)
G = np.array(np.meshgrid(fx, fy, indexing='ij')).reshape(2, -1).T; Z = model(G)
E, rms = [], []
for b1, b2 in ((bf, br), (bfo, bro)):
    T = rh.fill_table(*args[:1], b1, b2)['CzR']
    E.append((RegularGridInterpolator((b1, b2), T)(G) - Z).reshape(len(fx), len(fy)) / scale * 100)
    eq = RegularGridInterpolator((b1, b2), T, bounds_error=False, fill_value=None)(Xq) - model(Xq)
    rms.append(np.sqrt((wq * eq ** 2).sum()) / scale * 100)
lim = max(np.abs(e).max() for e in E)
fig, axs = plt.subplots(1, 2, figsize=(9, 3.9), sharey=True)
Fn, Rn = nf(fx), nr(fy); Tn = n(d[['dF', 'dR']].values)
for ax, e, (b1, b2), r, t in zip(axs, E, ((bf, br), (bfo, bro)), rms, ('Weighted', 'Optimised')):
    pc = ax.pcolormesh(Fn, Rn, e.T, cmap='RdBu_r', vmin=-lim, vmax=lim, shading='auto')
    ax.scatter(Tn[:, 0], Tn[:, 1], s=1, c='#333', alpha=.15, lw=0)
    for v in nf(b1): ax.axvline(v, c='#555', lw=.7, ls=':')
    for v in nr(b2): ax.axhline(v, c='#555', lw=.7, ls=':')
    ax.set_title('%s breakpoints, RMS %.2f%%' % (t, r), loc='left'); rh_axes(ax, pad=.01)
    if ax is not axs[0]: ax.set_ylabel('')
fig.colorbar(pc, ax=list(axs), pad=.01, fraction=.02).set_label('table - fitted map (% of range, synthetic CzR)')
save(fig, 'rhTableError.png')
