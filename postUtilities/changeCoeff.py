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

import os
import sys
import glob


dirlist = glob.glob(os.path.join('postProcessing','*'))
for dir in dirlist:
	if os.path.exists(os.path.join(dir,'0','coefficient_0.dat')):
		print(dir)
		os.system('mv %s %s' % (os.path.join(dir,'0','coefficient_0.dat'),os.path.join(dir,'0','coefficient.dat')))
