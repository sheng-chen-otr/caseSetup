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

#writes constant/tables for one 10 Hz cycle at 1 deg steps
import numpy as np

f = 10.0
theta = np.radians(np.arange(0, 361))
t = theta/(2*np.pi*f)
disp = {'zeroTable': 0*theta,
        'displacerTable': 0.006*np.sin(theta),          #+-6 mm
        'pistonTable': 0.008*(1 - np.cos(theta))}       #0 to 16 mm up

for name, dz in disp.items():
    with open('constant/tables/%s' % name, 'w') as fh:
        fh.write('FoamFile { version 2.0; format ascii; class scalarVectorTable; object %s; }\n' % name)
        fh.write('%d\n(\n' % len(t))
        for ti, di in zip(t, dz):
            fh.write('    (%.10g (0 0 %.10g))\n' % (ti, di))
        fh.write(')\n')
print('tables written')
