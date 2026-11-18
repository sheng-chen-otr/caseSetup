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

#puts the motion faceZones first, in interpolationTables order (works around the
#displacementInterpolation zone-index mismatch in v2412)
import re

first = ['hotAnchor', 'coldAnchor', 'displacerTop', 'displacerBottom', 'crown']
path = 'constant/polyMesh/faceZones'
txt = open(path).read()

def matchBrace(s, b):
    #index of the brace closing the one at s[b]
    depth = 0
    for j in range(b, len(s)):
        if s[j] == '{':
            depth += 1
        elif s[j] == '}':
            depth -= 1
            if depth == 0:
                return j
    raise ValueError('unbalanced braces')

#repair a header broken by the first version of this script (missing FoamFile close)
fStart = txt.index('{', txt.index('FoamFile'))
try:
    hEnd = matchBrace(txt, fStart) + 1
    m = re.search(r'\s*(\d+)\s*\(', txt[hEnd:])
    if m is None:
        raise ValueError
except ValueError:
    metaEnd = matchBrace(txt, txt.index('{', txt.index('meta')))
    txt = txt[:metaEnd+1] + '\n}\n' + txt[metaEnd+1:]
    hEnd = matchBrace(txt, fStart) + 1
    m = re.search(r'\s*(\d+)\s*\(', txt[hEnd:])

bodyStart = hEnd + m.end()
bodyEnd = txt.rindex(')')
body = txt[bodyStart:bodyEnd]

#split into "name { ... }" entries
entries, i = [], 0
while True:
    m = re.search(r'\S', body[i:])
    if not m:
        break
    s = i + m.start()
    b = body.index('{', s)
    j = matchBrace(body, b)
    entries.append((body[s:b].strip(), body[s:j+1]))
    i = j + 1

names = [n for n, e in entries]
order = [names.index(n) for n in first] + [k for k, n in enumerate(names) if n not in first]
newNames = [names[k] for k in order]

#keep the meta names list in the same order as the entries
header = re.sub(r'names\s+\d+\s*\([^)]*\)', 'names %d ( %s )' % (len(newNames), ' '.join(newNames)), txt[:hEnd])
out = header + '\n\n%d\n(\n' % len(entries) + '\n'.join(entries[k][1] for k in order) + '\n)\n'
open(path, 'w').write(out)
print('faceZones order: %s' % newNames)
