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

#writes the multi-solid ASCII STLs for the v2412 check case (all dims in m)
import os
import numpy as np

os.makedirs('constant/triSurface', exist_ok=True)

def tri(p0, p1, p2, n):
    #orients the triangle so its normal points along n
    p0, p1, p2 = np.array(p0, float), np.array(p1, float), np.array(p2, float)
    if np.dot(np.cross(p1-p0, p2-p0), n) < 0:
        p1, p2 = p2, p1
    return (p0, p1, p2)

def quad(p0, p1, p2, p3, n):
    return [tri(p0, p1, p2, n), tri(p0, p2, p3, n)]

def rectZ(x0, x1, y0, y1, z, n):
    return quad((x0,y0,z), (x1,y0,z), (x1,y1,z), (x0,y1,z), (0,0,n))

def frameZ(X0, X1, Y0, Y1, x0, x1, y0, y1, z, n):
    #rectangle at z with a rectangular opening
    return (rectZ(X0,x0,Y0,Y1,z,n) + rectZ(x1,X1,Y0,Y1,z,n) +
            rectZ(x0,x1,Y0,y0,z,n) + rectZ(x0,x1,y1,Y1,z,n))

def sides(x0, x1, y0, y1, z0, z1, s):
    #4 vertical walls, s=+1 normals out of the box, s=-1 into the box
    return (quad((x0,y0,z0),(x0,y1,z0),(x0,y1,z1),(x0,y0,z1),(-s,0,0)) +
            quad((x1,y0,z0),(x1,y1,z0),(x1,y1,z1),(x1,y0,z1),(s,0,0)) +
            quad((x0,y0,z0),(x1,y0,z0),(x1,y0,z1),(x0,y0,z1),(0,-s,0)) +
            quad((x0,y1,z0),(x1,y1,z0),(x1,y1,z1),(x0,y1,z1),(0,s,0)))

def writeStl(fileName, solids):
    with open('constant/triSurface/%s' % fileName, 'w') as f:
        for name, tris in solids:
            f.write('solid %s\n' % name)
            for p0, p1, p2 in tris:
                n = np.cross(p1-p0, p2-p0)
                n = n/np.linalg.norm(n)
                f.write(' facet normal %g %g %g\n  outer loop\n' % tuple(n))
                for p in (p0, p1, p2):
                    f.write('   vertex %.9g %.9g %.9g\n' % tuple(p))
                f.write('  endloop\n endfacet\n')
            f.write('endsolid %s\n' % name)

#chamber 100 x 100 x 60 mm, piston bore 20 x 20 mm offset in y, crown at z = 0.08
C, H = 0.05, 0.06
bx0, bx1, by0, by1, bTop = -0.01, 0.01, 0.015, 0.035, 0.10
writeStl('WALL-chamber.stl', [
    ('hot',   rectZ(-C, C, -C, C, 0.0, -1)),
    ('cold',  frameZ(-C, C, -C, C, bx0, bx1, by0, by1, H, 1)),
    ('side',  sides(-C, C, -C, C, 0.0, H, 1)),
    ('liner', sides(bx0, bx1, by0, by1, H, bTop, 1)),
    #crown is a separate shell, slightly larger than the bore
    ('crown', rectZ(bx0-0.001, bx1+0.001, by0-0.001, by1+0.001, 0.08, 1)),
])

#displacer 96 x 96 x 34 mm (2 mm side gap), 10 x 10 mm hole through the middle
D, z0, z1, h = 0.048, 0.013, 0.047, 0.005
writeStl('MOVE-displacer.stl', [
    ('top',      frameZ(-D, D, -D, D, -h, h, -h, h, z1, 1)),
    ('bottom',   frameZ(-D, D, -D, D, -h, h, -h, h, z0, -1)),
    ('side',     sides(-D, D, -D, D, z0, z1, 1)),
    ('holewall', sides(-h, h, -h, h, z0, z1, -1)),
])

#regenerator solid: wider than the hole, caps 0.5 mm inside the hole ends
r = 0.006
writeStl('REGN-matrix.stl', [
    ('cold01',  rectZ(-r, r, -r, r, z1-0.0005, 1)),
    ('hot01',   rectZ(-r, r, -r, r, z0+0.0005, -1)),
    ('regside', sides(-r, r, -r, r, z0+0.0005, z1-0.0005, 1)),
])

#sample disc at the displacer mid height (select mode test)
writeStl('SMP-mid.stl', [('mid', rectZ(-0.0055, 0.0055, -0.0055, 0.0055, 0.0301, 1))])
#regenerator caps split per PID, for the select-mode outlet sample zones
writeStl('SMP-capCold.stl', [('cold01', rectZ(-r, r, -r, r, z1-0.0005, 1))])
writeStl('SMP-capHot.stl',  [('hot01',  rectZ(-r, r, -r, r, z0+0.0005, -1))])

#displacer side only, for the gap distance refinement
writeStl('GAPREF-dispSide.stl', [('side', sides(-D, D, -D, D, z0, z1, 1))])
print('geometry written')
