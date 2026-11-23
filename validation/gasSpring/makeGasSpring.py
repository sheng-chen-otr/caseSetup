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

#gas spring validation case (piston in a closed cylinder, Kornhauser & Smith geometry by default):
#WALL-gasSpring.stl with PIDs head / liner / MOVE-piston (crown at BDC, separate shell, slightly larger than the bore)
#and the crank motion table (theta, piston z of the crown), axis z, head at the top

import os
import sys
import math
import argparse

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
import motionUtilities


def disc(radius, z, n, up):
    #fan of triangles, normal +z (up) or -z
    a = np.linspace(0.0, 2.0*math.pi, n + 1)
    tris = []
    for a0, a1 in zip(a[:-1], a[1:]):
        p0 = [radius*math.cos(a0), radius*math.sin(a0), z]
        p1 = [radius*math.cos(a1), radius*math.sin(a1), z]
        tris.append([[0.0, 0.0, z], p0, p1] if up else [[0.0, 0.0, z], p1, p0])
    return np.array(tris)


def tube(radius, z0, z1, n, dz):
    #cylinder wall, normals pointing inwards (into the gas)
    a = np.linspace(0.0, 2.0*math.pi, n + 1)
    zs = np.linspace(z0, z1, max(int(math.ceil((z1 - z0)/dz)), 1) + 1)
    tris = []
    for za, zb in zip(zs[:-1], zs[1:]):
        for a0, a1 in zip(a[:-1], a[1:]):
            p00 = [radius*math.cos(a0), radius*math.sin(a0), za]
            p01 = [radius*math.cos(a1), radius*math.sin(a1), za]
            p10 = [radius*math.cos(a0), radius*math.sin(a0), zb]
            p11 = [radius*math.cos(a1), radius*math.sin(a1), zb]
            tris.append([p00, p10, p01])
            tris.append([p01, p10, p11])
    return np.array(tris)


def crankPosition(theta, crank, conrod):
    #piston travel below TDC (m) at crank angle theta (deg, 0 = TDC)
    t = np.radians(theta)
    return crank*(1.0 - np.cos(t)) + conrod - np.sqrt(conrod**2 - (crank*np.sin(t))**2)


def main():
    parser = argparse.ArgumentParser(description='gas spring geometry and motion table for internalSetup')
    parser.add_argument('mshDir', help='geometry folder (<job>/02_reference/MSH)')
    parser.add_argument('motionFile', help='motion table to write (csv)')
    parser.add_argument('--bore', type=float, default=0.0508, help='bore (m)')
    parser.add_argument('--stroke', type=float, default=0.0762, help='stroke (m)')
    parser.add_argument('--clearance', type=float, default=0.0762, help='head to crown at TDC (m)')
    parser.add_argument('--conrod', type=float, default=0.183, help='connecting rod length (m)')
    parser.add_argument('--segments', type=int, default=180, help='facets around the circumference')
    parser.add_argument('--dz', type=float, default=0.005, help='axial facet length of the liner (m)')
    args = parser.parse_args()

    n = args.segments
    #polygon radius with the circle's area (volume exact)
    radius = 0.5*args.bore*math.sqrt(2.0*math.pi/(n*math.sin(2.0*math.pi/n)))
    zBdc = 0.0
    zTdc = args.stroke
    zHead = zTdc + args.clearance
    print('\tbore %g m, stroke %g m, clearance %g m, conrod %g m: volume ratio %.4f'
          % (args.bore, args.stroke, args.clearance, args.conrod, (args.stroke + args.clearance)/args.clearance))

    #crown drawn at BDC (REF_ANGLE 180), liner runs past it so the crown seals at every angle
    pidTris = {'head': disc(radius, zHead, n, False),
               'liner': tube(radius, zBdc - 0.05*args.stroke, zHead, n, args.dz),
               'MOVE-piston': disc(1.02*radius, zBdc, n, True)}
    os.makedirs(args.mshDir, exist_ok=True)
    stlPath = os.path.join(args.mshDir, 'WALL-gasSpring.stl')
    motionUtilities.writeStlTriangles(stlPath, pidTris)
    print('\twritten %s (head z %g, crown z %g at BDC)' % (stlPath, zHead, zBdc))

    #crown z vs crank angle, 0 = TDC
    theta = np.arange(0.0, 360.0, 0.5)
    z = zTdc - crankPosition(theta, 0.5*args.stroke, args.conrod)
    os.makedirs(os.path.dirname(os.path.abspath(args.motionFile)), exist_ok=True)
    with open(args.motionFile, 'w', newline='\n') as f:
        f.write('theta,piston\n')
        for t, p in zip(theta, z):
            f.write('%.9g,%.12g\n' % (t, p))
    print('\twritten %s (crown z %.6g..%.6g m)' % (args.motionFile, z.min(), z.max()))


if __name__ == '__main__':
    main()
