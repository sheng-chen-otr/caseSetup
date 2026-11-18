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

#compares the mesh flux through the select-mode sample zones with the analytic
#displacer velocity x fluid area, and reports the gas mass drift
import sys
import math
import glob

f = 10.0
tol = 1e-3       #relative
zones = {        #fo name: (fluid area [m2], sign of the surface normal in z)
    'capColdMesh':  (1e-4, 1),
    'capHotMesh':   (1e-4, -1),
    'smpMidMesh':   (1e-4, 1),
    'planeMidMesh': (1e-4 + 0.1**2 - 0.096**2, 1),
}

def lastRow(fo):
    #run-time FOs write to <fo>/0, postProcess FOs to <fo>/<time>
    path = sorted(glob.glob('postProcessing/%s/*/%s.dat' % (fo, 'volFieldValue' if fo == 'gasMass' else 'surfaceFieldValue')))[-1]
    cols, rows = None, []
    for line in open(path):
        if line.startswith('#'):
            cols = line[1:].split()
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    return dict(zip(cols, rows[-1])), dict(zip(cols, rows[0]))

ok = True
for fo, (area, sign) in zones.items():
    last, first = lastRow(fo)
    t = last['Time']
    v = 0.006*2*math.pi*f*math.cos(2*math.pi*f*t)
    expect = sign*v*area
    got = last['sum(meshPhi)']
    err = abs(got - expect)/abs(expect)
    passed = err < tol
    ok = ok and passed
    print('%-13s t=%g  sum(meshPhi) %.5e  expected %.5e  rel err %.1e  %s'
          % (fo, t, got, expect, err, 'OK' if passed else 'FAIL'))

last, first = lastRow('gasMass')
drift = (last['volIntegrate(rho)'] - first['volIntegrate(rho)'])/first['volIntegrate(rho)']
print('gas mass drift over the run: %.2e' % drift)
sys.exit(0 if ok else 1)
