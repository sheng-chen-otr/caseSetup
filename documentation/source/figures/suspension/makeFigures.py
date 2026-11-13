"""
Figures for the "Suspension Kinematics" theory section.

    python makeFigures.py [hardpoints.cfg] [scale]

Runs the real solver (rideHeightUtils/doubleWishbonePushrodKinematics.py) on a hardpoint file (default: the
template, in mm) and writes the PNGs next to this file. With NORMALISE = True no real value can be read off:
linkage drawings have no ticks and are scaled by the wheel-centre height, travel is scaled by MAX_TRAVEL, and
camber, toe, damper stroke and steer are shown as a fraction of their largest value in the plot.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch

NORMALISE = True
MAX_TRAVEL = 0.03          # m, bump/droop range shown
MAX_STEER = 20.0           # deg, front-left road-wheel angle range shown

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT))
import utilities as u  # noqa: E402
from rideHeightUtils.doubleWishbonePushrodKinematics import (  # noqa: E402
    DoubleWishbonePushrodSolver, resolve_pushrod_mount, _rodrigues, RESIDUAL_TOL)

HP_FILE = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'rideHeightUtils' / 'suspensionHardpoints_template.cfg'
SCALE = float(sys.argv[2]) if len(sys.argv) > 2 else 0.001

RED, BLUE, GREY, GREEN, ORANGE = '#c0392b', '#2c7fb8', '#555555', '#2e8b57', '#e67e22'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 200, 'savefig.bbox': 'tight'})

corners, _ = u._load_suspension_hardpoints_cfg(str(HP_FILE))
corners = u._scale_suspension_hardpoints(corners, SCALE)
for c in corners.values():
    c['pushrod_mount'] = resolve_pushrod_mount(c)
FL = DoubleWishbonePushrodSolver(corners['fl'])
FR = DoubleWishbonePushrodSolver(corners['fr'])
HP = corners['fl']
WC0 = np.array(HP['wheel_center_static'])
LSCALE = WC0[2]            # drawing length scale: wheel-centre height


def save(fig, name):
    fig.savefig(HERE / name); plt.close(fig); print('wrote', name)


def frac(v):
    """Scale an array to a fraction of its largest magnitude (when normalising)."""
    v = np.asarray(v, float)
    if not NORMALISE:
        return v
    m = np.abs(v).max()
    return v / m if m > 0 else v


def travel_axis(ax):
    ax.set_xlabel('wheel travel' + (' (fraction of max)' if NORMALISE else ' (mm)'))
    ax.axvline(0, color='#bbb', lw=.6); ax.axhline(0, color='#bbb', lw=.6)


def yz(p):
    """Front view coordinates (y, z) relative to the static wheel centre, scaled."""
    p = np.asarray(p, float)
    return np.array([p[1] - WC0[1], p[2] - WC0[2]]) / LSCALE


def bare(ax):
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def draw_corner(ax, sol=None, color=None, lw=1.6, alpha=1.0, labels=False, wheel=True):
    """Front view of the FL linkage, static (sol=None) or at a solved pose."""
    i = FL.inner
    if sol is None:
        up = {k: v for k, v in FL.outer0_global.items()}
        prod = FL.pushrod_outer_static
        rk = FL._rocker_points(0.0)
        R = np.eye(3)
        tie_in = i['tie']
    else:
        up = sol['upright_points']; prod = up['pushrod']; rk = sol['rocker_points']
        R = _rodrigues(sol['state'].rvec)
        tie_in = FL.tie_inner_static + sol['rack_travel'] * FL.rack_axis
    col = lambda c: color or c   # noqa: E731
    seg = lambda a, b, c, **k: ax.plot(*np.c_[yz(a), yz(b)], color=col(c), lw=lw, alpha=alpha,  # noqa: E731
                                        solid_capstyle='round', **k)
    for arm, c in (('uca', BLUE), ('lca', BLUE)):
        seg(i['%s_f' % arm], up[arm], c); seg(i['%s_r' % arm], up[arm], c)
    seg(tie_in, up['tie'], GREEN)
    seg(up['uca'], up['lca'], GREY, )                         # upright
    seg(up['uca'], up['wheel_center'], GREY); seg(up['lca'], up['wheel_center'], GREY)
    seg(prod, rk['pushrod'], ORANGE)
    piv = np.array(HP['rocker']['pivot'])
    ax.add_patch(Polygon([yz(piv), yz(rk['pushrod']), yz(rk['damper'])], closed=True,
                         fc=col(ORANGE), ec=col(ORANGE), alpha=.35 * alpha))
    seg(rk['damper'], HP['rocker']['damper_chassis'], RED)
    if wheel:
        r, hw = WC0[2], 0.1                                    # tyre radius ~ wheel-centre height
        box = np.array([[0, -hw, -r], [0, hw, -r], [0, hw, r], [0, -hw, r]])
        wc = up['wheel_center']
        pts = [yz(wc + R @ b) for b in box]
        ax.add_patch(Polygon(pts, closed=True, fc='none', ec=col('#333'), lw=lw * .8, alpha=alpha))
    pts = [i['uca_f'], i['uca_r'], i['lca_f'], i['lca_r'], tie_in, piv, HP['rocker']['damper_chassis']]
    P = np.array([yz(p) for p in pts])
    ax.plot(P[:, 0], P[:, 1], 's', ms=3.5, color=col('k'), alpha=alpha)
    P = np.array([yz(up[k]) for k in ('uca', 'lca', 'tie')] + [yz(prod)])
    ax.plot(P[:, 0], P[:, 1], 'o', ms=3.5, mfc='white', mec=col('k'), alpha=alpha)
    if labels:
        lab = [('UCA', up['uca'], (.12, .10)), ('LCA', up['lca'], (.12, -.12)),
               ('tie rod', up['tie'], (.18, -.05)), ('pushrod', .5 * (prod + rk['pushrod']), (.18, .02)),
               ('rocker', piv, (.05, -.35)), ('damper', .5 * (rk['damper'] + np.array(HP['rocker']['damper_chassis'])),
                                             (0, .2)),
               ('upright', .5 * (up['uca'] + up['lca']), (-.25, .05)), ('wheel', wc + [0, -.1, .1], (-.3, .25))]
        for t, p, off in lab:
            q = yz(p)
            ax.annotate(t, q, q + off, fontsize=8, arrowprops=dict(arrowstyle='-', color='#888', lw=.6))


def ground(ax, x0, x1):
    z = -WC0[2] / LSCALE
    ax.plot([x0, x1], [z, z], color='#999', lw=1); ax.fill_between([x0, x1], z - .06, z, color='#eee', lw=0)


# ------------------------------------------------------------------ 1: pitch, roll, heave
fig, (a0, a1) = plt.subplots(1, 2, figsize=(7.5, 2.8))
for ax, title, (zl, zr), lab, (cl, cr) in (
        (a0, 'Side view: front up, rear down = pitch', (.18, -.08), ('front', 'rear'), ('$\\Delta z_F$', '$\\Delta z_R$')),
        (a1, 'Rear view: right up, left down = roll', (-.06, .16), ('left', 'right'), ('$\\Delta z_L$', '$\\Delta z_R$'))):
    L = 2.0 if ax is a0 else 1.2
    x = np.array([0, L]); z = np.array([zl, zr]) + .45
    ax.plot([0, L], [.45, .45], '--', color='#aaa', lw=1)
    ax.add_patch(Polygon([[0, .45 + .0], [L, .45], [L, .65], [0, .65]], fc='#eee', ec='#bbb', ls='--'))
    ang = np.arctan2(zr - zl, L)
    body = np.array([[0, 0], [L, 0], [L, .2], [0, .2]])
    rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    ax.add_patch(Polygon(body @ rot.T + [0, z[0]], fc='#cfe3f3', ec=BLUE, lw=1.3))
    for xi, zi, ci in zip(x, z, (cl, cr)):
        ax.add_patch(FancyArrowPatch((xi, .45), (xi, zi), arrowstyle='-|>', mutation_scale=9, color=RED, lw=1.2))
        ax.text(xi + (-.08 if xi == 0 else .06), (.45 + zi) / 2, ci, ha='right' if xi == 0 else 'left',
                va='center', color=RED)
    mid = z.mean()
    ax.plot([L / 2], [mid], 'k+', ms=8)
    ax.text(L / 2, mid - .1, 'heave = average', ha='center', va='top', fontsize=8)
    ax.text(0, .2, lab[0], ha='center'); ax.text(L, .2, lab[1], ha='center')
    ax.set_xlim(-.45, L + .45); ax.set_ylim(.1, 1.05); bare(ax); ax.set_title(title, loc='left')
save(fig, 'suspPitchRoll.png')

# ------------------------------------------------------------------ 2: the corner layout
fig, ax = plt.subplots(figsize=(6.5, 4.0))
draw_corner(ax, labels=True)
ground(ax, -.6, 2.3)
ax.plot([], [], 's', color='k', ms=3.5, label='fixed to the chassis')
ax.plot([], [], 'o', mfc='white', mec='k', ms=3.5, label='moves with the wheel')
ax.legend(loc='lower right', frameon=False, fontsize=8)
bare(ax); ax.set_title('Front-left corner seen from behind (car centreline off to the right)', loc='left', pad=14)
save(fig, 'suspCornerLayout.png')

# ------------------------------------------------------------------ 3: solver convergence
fig, ax = plt.subplots(figsize=(5.5, 2.8))
cases = ((0.03, 0.0, 'none', 'bump, no steer'), (-0.03, 0.0, 'none', 'droop, no steer'),
         (0.0, 15.0, 'angle', 'static, 15 deg steer'))
for (dz, st, mode, lab), c in zip(cases, (BLUE, GREEN, RED)):
    err = [FL.solve_for_wheel_travel(dz, steer=st, steer_mode=mode, max_iter=k)['residual_norm'] * 1000
           for k in range(0, 9)]
    stop = next((k for k, e in enumerate(err) if e < 1e-6), len(err) - 1)   # solver stops below 1e-9 m
    err = err[:stop + 1]
    ax.semilogy(range(len(err)), np.maximum(err, 1e-12), 'o-', color=c, ms=3.5, label=lab)
ax.axhline(RESIDUAL_TOL * 1000, color='#888', ls='--', lw=.8)
ax.text(4.2, RESIDUAL_TOL * 1000 * 1.6, 'pass line (0.001 mm)', ha='right', fontsize=8, color='#666')
ax.set_xlabel('solver step'); ax.set_xlim(-.3, 4.3); ax.set_ylabel('total error (mm)')
ax.set_title('How fast the error falls (each solve starts from the static pose)', loc='left')
ax.legend(frameon=False, fontsize=8)
save(fig, 'suspSolverConvergence.png')

# ------------------------------------------------------------------ 4: camber, toe, damper vs travel
s = np.linspace(-MAX_TRAVEL, MAX_TRAVEL, 41)
res, state = [], None
for dz in s:          # warm start in a sweep, as caseSetup does
    r = FL.solve_for_wheel_travel(dz, initial_state=state); res.append(r); state = r['state']
order = np.argsort(s)
cam = np.array([r['camber_deg'] for r in res]); toe = np.array([r['toe_deg'] for r in res])
dmp = np.array([r['damper_delta'] for r in res]) * 1000
ms = s / MAX_TRAVEL if NORMALISE else s * 1000
fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.6))
for ax, v, t, unit in ((axs[0], cam, 'Camber change', 'deg'), (axs[1], toe, 'Toe change', 'deg'),
                       (axs[2], dmp, 'Damper length change', 'mm')):
    ax.plot(ms, frac(v), color=BLUE, lw=1.6); travel_axis(ax)
    ax.set_title(t, loc='left'); ax.set_ylabel('fraction of max' if NORMALISE else unit)
mr = -np.gradient(dmp, s * 1000)
axs[2].text(.05, .08, 'slope = motion ratio', transform=axs[2].transAxes, fontsize=8, color='#444')
fig.tight_layout(w_pad=1.5)
fig.text(.5, -.04, 'positive travel is bump (wheel moving up towards the body)', ha='center', fontsize=8, color='#444')
save(fig, 'suspCurves.png')
print('  FL motion ratio min/max: %.3f %.3f' % (mr.min(), mr.max()))

# ------------------------------------------------------------------ 5: steering / Ackermann
steer = np.linspace(-MAX_STEER, MAX_STEER, 41)
fl_ang, fr_ang, rack = [], [], []
sl = sr = None
for a in steer:
    r1 = FL.solve_for_wheel_travel(0.0, steer=a, steer_mode='angle', initial_state=sl); sl = r1['state']
    r2 = FR.solve_for_wheel_travel(0.0, steer=r1['rack_travel'], steer_mode='rack', initial_state=sr); sr = r2['state']
    fl_ang.append(r1['toe_deg']); fr_ang.append(r2['toe_deg']); rack.append(r1['rack_travel'])
fl_ang, fr_ang = np.array(fl_ang), np.array(fr_ang)
m = np.abs(np.r_[fl_ang, fr_ang]).max() if NORMALISE else 1.0
fig, (a0, a1) = plt.subplots(1, 2, figsize=(7.5, 2.9))
a0.plot(fl_ang / m, fl_ang / m, '--', color='#999', lw=1, label='both wheels the same')
a0.plot(fl_ang / m, fr_ang / m, color=RED, lw=1.6, label='front-right from the solver')
a0.set_xlabel('front-left angle' + (' (fraction of max)' if NORMALISE else ' (deg)'))
a0.set_ylabel('front-right angle' + (' (fraction of max)' if NORMALISE else ' (deg)'))
a0.set_title('Same rack travel, different angles', loc='left'); a0.legend(frameon=False, fontsize=8)
a0.set_aspect('equal')
diff = fr_ang - fl_ang
a1.plot(fl_ang / m, frac(diff), color=RED, lw=1.6); a1.axhline(0, color='#bbb', lw=.6)
a1.axvline(0, color='#bbb', lw=.6)
a1.set_xlabel(a0.get_xlabel()); a1.set_ylabel('FR minus FL' + (' (fraction of max)' if NORMALISE else ' (deg)'))
a1.set_title('The difference between the two wheels', loc='left')
fig.tight_layout(w_pad=2)
save(fig, 'suspAckermann.png')

# ------------------------------------------------------------------ 6: poses
fig, ax = plt.subplots(figsize=(6.5, 4.0))
draw_corner(ax, FL.solve_for_wheel_travel(-MAX_TRAVEL), color=GREEN, lw=1.1, alpha=.8)
draw_corner(ax, FL.solve_for_wheel_travel(MAX_TRAVEL), color=RED, lw=1.1, alpha=.8)
draw_corner(ax, None, color='#222', lw=1.4)
for c, t in ((RED, 'bump'), ('#222', 'static'), (GREEN, 'droop')):
    ax.plot([], [], color=c, label=t)
ax.legend(loc='lower right', frameon=False, fontsize=8)
bare(ax); ax.set_title('The same corner at droop, static and bump (seen from behind)', loc='left', pad=14)
save(fig, 'suspPoses.png')
