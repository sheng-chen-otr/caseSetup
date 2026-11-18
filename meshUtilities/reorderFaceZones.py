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

#puts the given faceZones first, in the given order, in constant/polyMesh/faceZones and in every
#processor*/constant/polyMesh/faceZones. works around the displacementInterpolation zone-index
#mismatch in ESI v2412 (zone ranges recorded by interpolationTables position, displacements looked
#up by faceZone index): the motion zones must be the first faceZones, in table order.
#usage: python3 reorderFaceZones.py zoneA zoneB ... [--case <dir>]
#stdlib only so it runs on cluster nodes

import os
import re
import sys
import glob


#byte sizes of binary list items (label/scalar sizes come from the header arch entry)
binarySizes = {'bool': 1}


def binarySkips(s, binary):
    #byte ranges of binary list data "List<type> N(...)", so braces inside them are not counted
    skips = {}
    if not binary:
        return skips
    for m in re.finditer(rb'List<(\w+)>\s*(\d+)\s*\(', s):
        size = binarySizes.get(m.group(1).decode())
        if size is None:
            raise ValueError('unknown binary list type %s' % (m.group(1).decode()))
        skips[m.end() - 1] = m.end() + int(m.group(2))*size
    return skips


def matchBrace(s, b, skips=None):
    #index of the brace closing the one at s[b], jumping over binary list data
    skips = skips or {}
    depth, j = 0, b
    while j < len(s):
        if j in skips:
            j = skips[j]
            continue
        c = s[j:j+1]
        if c == b'{':
            depth += 1
        elif c == b'}':
            depth -= 1
            if depth == 0:
                return j
        j += 1
    raise ValueError('unbalanced braces')


def reorderFile(path, first):
    #works on bytes so ascii and binary files are handled the same way
    with open(path, 'rb') as f:
        txt = f.read()
    fStart = txt.index(b'{', txt.index(b'FoamFile'))
    hEnd = matchBrace(txt, fStart) + 1
    header = txt[:hEnd].decode()
    binary = re.search(r'format\s+binary', header) is not None
    arch = re.search(r'label=(\d+);scalar=(\d+)', header)
    if arch:
        binarySizes['label'] = int(arch.group(1))//8
        binarySizes['scalar'] = int(arch.group(2))//8
    elif binary:
        binarySizes['label'], binarySizes['scalar'] = 4, 8
    m = re.search(rb'\s*(\d+)\s*\(', txt[hEnd:])
    if m is None:
        raise ValueError('%s: no zone list found' % (path))
    bodyStart = hEnd + m.end()
    bodyEnd = txt.rindex(b')')
    body = txt[bodyStart:bodyEnd]
    skips = binarySkips(body, binary)

    #split into "name { ... }" entries
    entries, i = [], 0
    while True:
        m = re.search(rb'\S', body[i:])
        if not m:
            break
        s = i + m.start()
        b = body.index(b'{', s)
        j = matchBrace(body, b, skips)
        entries.append((body[s:b].strip().decode(), body[s:j+1]))
        i = j + 1

    names = [n for n, e in entries]
    missing = [n for n in first if n not in names]
    if missing:
        raise ValueError('%s: faceZones %s not found (have %s)' % (path, missing, names))
    order = [names.index(n) for n in first] + [k for k, n in enumerate(names) if n not in first]
    newNames = [names[k] for k in order]

    #keep the header meta names list in the same order as the entries
    header = re.sub(r'names\s+\d+\s*\([^)]*\)', 'names %d ( %s )' % (len(newNames), ' '.join(newNames)), header)
    out = (header + '\n\n%d\n(\n' % len(entries)).encode() + b'\n'.join(entries[k][1] for k in order) + b'\n)\n'
    with open(path, 'wb') as f:
        f.write(out)
    return newNames


def main(argv):
    caseDir = '.'
    if '--case' in argv:
        k = argv.index('--case')
        caseDir = argv[k + 1]
        argv = argv[:k] + argv[k + 2:]
    first = argv
    if not first:
        sys.exit('usage: reorderFaceZones.py zoneA zoneB ... [--case <dir>]')
    paths = sorted(glob.glob(os.path.join(caseDir, 'processor*', 'constant', 'polyMesh', 'faceZones')))
    serial = os.path.join(caseDir, 'constant', 'polyMesh', 'faceZones')
    if os.path.isfile(serial):
        paths.append(serial)
    if not paths:
        sys.exit('ERROR! no faceZones file found in %s' % (caseDir))
    for path in paths:
        newNames = reorderFile(path, first)
    print('faceZones reordered in %d file(s): %s' % (len(paths), ' '.join(newNames)))


if __name__ == '__main__':
    main(sys.argv[1:])
