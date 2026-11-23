"""
Copyright (c) 2026 Zhi Sheng Chen

Licensed under the PolyForm Noncommercial License 1.0.0 (the "License").
You may not use this file except in compliance with the License. Any
commercial use requires a separate commercial license from the copyright
holder.

You may obtain a copy of the License at:
    https://polyformproject.org/licenses/noncommercial/1.0.0

Or see the LICENSE file in the root of this repository.

This software is provided "as is", without warranty of any kind, express or
implied. See the License for the specific language governing permissions
and limitations.
"""

#gas spring results per cycle, as in Kornhauser & Smith / Willich & White:
#  p/p_min vs crank angle (0 = TDC), loss Psi = -loop p dV / loop |p dV| (net work into the gas, positive),
#  Peclet number Pe = omega Dh^2/(4 alpha)
#  (Dh = 4 V0/Aw0 at mid stroke, alpha = (alpha_max + alpha_min)/2 over the cycle from the mean gas state),
#  mean T at TDC per cycle (periodic when it changes by less than 0.1 K)
#reads postProcessing gasAvg (volAverage p, T, rho) and gasMass, gas volume = mass/mean rho
#writes gasSpring.csv (last cycle) and gasSpringCycles.csv (per cycle) in the case

import os
import sys
import math
import argparse

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, 'postUtilities'))
import postInternal as pi
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
import writeInternal


def gasProperties(config):
    #cp, R, mu(T) (constant or Sutherland), Pr from the case's GLOBAL_MATERIAL as the setup resolves them
    full = {'GLOBAL_MATERIAL': dict((k, [v]) for k, v in config['GLOBAL_MATERIAL'].items())}
    props = writeInternal.getGasProperties(full)
    sutherland = pi.getSetting(config, 'GLOBAL_MATERIAL', 'TRANSPORT', 'const').lower() == 'sutherland'
    def mu(T):
        if sutherland:
            return props['AS']*T**1.5/(T + props['TS'])
        return props['MU']*np.ones_like(T)
    return props['CP'], 8314.46/props['MOL_WEIGHT'], mu, props['PR']


def main():
    parser = argparse.ArgumentParser(description='gas spring validation results')
    parser.add_argument('case', nargs='?', default='.', help='case directory')
    parser.add_argument('--bore', type=float, default=0.0508, help='bore (m)')
    parser.add_argument('--stroke', type=float, default=0.0762, help='stroke (m)')
    parser.add_argument('--clearance', type=float, default=0.0762, help='head to crown at TDC (m)')
    args = parser.parse_args()

    case = os.path.abspath(args.case)
    config = pi.readCaseSetup(case)
    post = os.path.join(case, 'postProcessing')
    freq = float(pi.getSetting(config, 'GLOBAL_CONTROL', 'FREQ'))
    start = float(pi.getSetting(config, 'SETUP_INFO', 'START_ANGLE', '0'))
    names, g = pi.readDat(os.path.join(post, 'gasAvg'), 'volFieldValue')
    mNames, m = pi.readDat(os.path.join(post, 'gasMass'), 'volFieldValue')
    t = g[:, 0]
    p = g[:, names.index('volAverage(p)')]
    T = g[:, names.index('volAverage(T)')]
    rho = g[:, names.index('volAverage(rho)')]
    mass = np.interp(t, m[:, 0], m[:, 1])
    V = mass/rho
    #crank angle from the start angle (prescribed speed), 0 = TDC
    theta = start + 360.0*freq*t
    cycle = np.floor((theta - start)/360.0 + 1e-9).astype(int)

    cp, R, mu, pr = gasProperties(config)
    omega = 2.0*math.pi*freq
    area = 0.25*math.pi*args.bore**2
    v0 = area*(args.clearance + 0.5*args.stroke)
    aw0 = 2.0*area + math.pi*args.bore*(args.clearance + 0.5*args.stroke)
    dh = 4.0*v0/aw0
    alpha = mu(T)*cp/pr/(rho*cp)

    rows = []
    print('\tcase %s: %g Hz (%g rpm), Dh %.5g m, gas mass drift %.2e' % (os.path.basename(case), freq, 60.0*freq, dh,
                                                                     mass[-1]/mass[0] - 1.0))
    print('\tcycle  Psi        Pe        p_max/p_min  T_TDC (K)   dT_TDC (K)')
    tdcPrev = None
    for c in range(cycle.max() + 1):
        sel = cycle == c
        if sel.sum() < 10 or (theta[sel][-1] - theta[sel][0]) < 359.0:
            continue
        #closed loop (piecewise linear p over V), loop of |p dV|; loss positive = net work put into the gas
        pc, vc = p[sel], V[sel]
        dv = np.diff(np.append(vc, vc[0]))
        pm = 0.5*(pc + np.append(pc[1:], pc[0]))
        psi = -np.sum(pm*dv)/np.sum(np.abs(pm*dv))
        pe = omega*dh**2/(4.0*0.5*(alpha[sel].max() + alpha[sel].min()))
        #mean T at TDC (crank angle a multiple of 360)
        th = np.mod(theta[sel], 360.0)
        iTdc = int(np.argmin(np.minimum(th, 360.0 - th)))
        tdc = T[sel][iTdc]
        rows.append((c + 1, psi, pe, pc.max()/pc.min(), tdc, tdc - tdcPrev if tdcPrev is not None else float('nan')))
        print('\t%5d  %.5f  %8.4g  %11.5f  %10.3f  %+10.3f' % rows[-1])
        tdcPrev = tdc
    with open(os.path.join(case, 'gasSpringCycles.csv'), 'w', newline='\n') as f:
        f.write('cycle,psi,pe,pRatio,T_TDC,dT_TDC\n')
        for row in rows:
            f.write('%d,%.9g,%.9g,%.9g,%.9g,%.9g\n' % row)
    #last complete cycle: p/p_min vs crank angle (-180..180, 0 = TDC)
    if rows:
        sel = cycle == rows[-1][0] - 1
        th = np.mod(theta[sel] + 180.0, 360.0) - 180.0
        order = np.argsort(th)
        with open(os.path.join(case, 'gasSpring.csv'), 'w', newline='\n') as f:
            f.write('theta,pRatio,volumeRatio,T\n')
            for k in order:
                f.write('%.6g,%.9g,%.9g,%.9g\n' % (th[k], p[sel][k]/p[sel].min(), V[sel][k]/V[sel].min(), T[sel][k]))
        print('\twritten gasSpring.csv (last cycle) and gasSpringCycles.csv')


if __name__ == '__main__':
    main()
