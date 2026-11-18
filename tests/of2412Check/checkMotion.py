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

#compares moved mesh points against the expected z-interpolation per band
import os
import re
import sys
import numpy as np

def readPoints(path):
    txt = open(path).read()
    body = txt[txt.index('(', txt.index('}')):]
    nums = re.findall(r'\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)', body)
    return np.array(nums, dtype=float)

p0 = readPoints('constant/polyMesh/points')
z0 = p0[:, 2]
zHot, zBot, zTop, zCold, zCrown = 0.0, 0.013, 0.047, 0.06, 0.08
f = 10.0
print('points: %d, z range %.5f .. %.5f' % (len(p0), z0.min(), z0.max()))
tol = 1e-6    #m
worst = 0.0

for tDir in ['0.025', '0.05', '0.075', '0.1']:
    pPath = '%s/polyMesh/points' % tDir
    if not os.path.isfile(pPath):
        print('%s: no moved points written' % tDir)
        continue
    p = readPoints(pPath)
    th = 2*np.pi*f*float(tDir)
    d, pis = 0.006*np.sin(th), 0.008*(1 - np.cos(th))
    expect = np.zeros_like(z0)
    b = z0 <= zBot
    expect[b] = (z0[b] - zHot)/(zBot - zHot)*d
    b = (z0 > zBot) & (z0 < zTop)
    expect[b] = d
    b = (z0 >= zTop) & (z0 <= zCold)
    expect[b] = d*(zCold - z0[b])/(zCold - zTop)
    b = z0 > zCold
    expect[b] = (z0[b] - zCold)/(zCrown - zCold)*pis
    dz = p[:, 2] - z0
    worst = max(worst, np.abs(dz - expect).max())
    dxy = np.abs(p[:, :2] - p0[:, :2]).max()
    print('t=%s disp=%+.4f piston=%+.4f | max |dz-expected| %.2e | max |dx,dy| %.2e'
          % (tDir, d, pis, np.abs(dz - expect).max(), dxy))
    for name, lo, hi in [('hot space', zHot, zBot), ('displacer band', zBot, zTop),
                         ('cold space', zTop, zCold), ('piston bore', zCold, zCrown)]:
        b = (z0 >= lo) & (z0 <= hi)
        print('    %-15s %7d pts  max err %.2e' % (name, b.sum(), np.abs(dz[b] - expect[b]).max()))

print('motion check %s (worst error %.2e m, tolerance %.0e m)' % ('OK' if worst < tol else 'FAILED', worst, tol))
sys.exit(0 if worst < tol else 1)
