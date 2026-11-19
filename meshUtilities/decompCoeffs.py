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

#sets numberOfSubdomains and the hierarchical/simple split n (nx ny nz) in decomposeParDict for the
#processor count of the run (no scotch needed). the split is chosen so the subdomains of the background
#mesh box (system/blockMeshDict vertices) are as close to cubes as possible
#usage (from the case directory): python3 decompCoeffs.py system/decomposeParDict <nProcs> [nx ny nz]
#with nx ny nz (DECOMP in GLOBAL_COMPUTE_SETUP) that split is used and must give nProcs subdomains

import re
import sys


def blockMeshLengths(path):
    #box lengths from the vertices of blockMeshDict
    with open(path) as f:
        text = f.read()
    block = re.search(r'vertices\s*\((.*?)\);', text, re.S)
    if block is None:
        sys.exit('ERROR! no vertices in %s' % (path))
    points = [[float(v) for v in p.split()] for p in re.findall(r'\(([^()]+)\)', block.group(1))]
    return [max(p[i] for p in points) - min(p[i] for p in points) for i in range(3)]


def chooseSplit(nProcs, lengths):
    #n (nx ny nz) with nx*ny*nz = nProcs, subdomain aspect ratio max(L/n)/min(L/n) smallest
    best = None
    for nx in range(1, nProcs + 1):
        if nProcs % nx:
            continue
        for ny in range(1, nProcs//nx + 1):
            if (nProcs//nx) % ny:
                continue
            n = (nx, ny, nProcs//(nx*ny))
            sizes = [l/k for l, k in zip(lengths, n)]
            aspect = max(sizes)/max(min(sizes), 1e-300)
            if best is None or aspect < best[0] - 1e-12:
                best = (aspect, n)
    return best[1]


def setSplit(dictPath, nProcs, n):
    with open(dictPath) as f:
        text = f.read()
    text = re.sub(r'numberOfSubdomains\s+[^;]*;', 'numberOfSubdomains %d;' % (nProcs), text)
    text = re.sub(r'\n?coeffs\s*\{[^}]*\}\n?', '\n', text)
    entry = 'coeffs\n{\n    n               (%d %d %d);\n}\n' % n
    end = text.rfind('// ****')
    text = text[:end] + entry + '\n' + text[end:] if end != -1 else text + '\n' + entry
    with open(dictPath, 'w') as f:
        f.write(text)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) not in (2, 5):
        sys.exit('usage: decompCoeffs.py system/decomposeParDict <nProcs> [nx ny nz]')
    dictPath, nProcs = argv[0], int(argv[1])
    if len(argv) == 5:
        n = tuple(int(v) for v in argv[2:])
        if n[0]*n[1]*n[2] != nProcs:
            sys.exit('ERROR! DECOMP (%d %d %d) gives %d subdomains but the run has %d tasks: change DECOMP '
                     'or the task count' % (n + (n[0]*n[1]*n[2], nProcs)))
    else:
        n = chooseSplit(nProcs, blockMeshLengths('system/blockMeshDict'))
    setSplit(dictPath, nProcs, n)
    print('decomposition: %d subdomains, n (%d %d %d)' % ((nProcs,) + n))


if __name__ == '__main__':
    main()
