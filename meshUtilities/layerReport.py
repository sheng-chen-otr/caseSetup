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

#layer coverage of the MESH_GROUPS patches from the final layer table of a snappyHexMesh log.
#snappy drops layers where they do not fit (narrow gaps, sharp edges, quality checks) without stopping, so a
#thermal layer can be missing on part of a heat-exchanging face. With remeshing every segment mesh is checked:
#a face covered in one segment and not in the next changes the wall heat flux at the map.
#usage: python3 layerReport.py <log.snappyHexMesh> <system/meshGroups> [--name <mesh>] [--min <fraction>]
#never fails the run (exit 0), prints WARNING lines. stdlib only so it runs on cluster nodes

import sys

MIN_FRACTION = 0.5


def readMeshGroups(path):
    #{patch: (group, target layers)} from system/meshGroups
    groups = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3 and not parts[0].startswith('#'):
                groups[parts[0]] = (parts[1], int(parts[2]))
    return groups


def readLayerTable(path):
    #{patch: (faces, target, mesh layers, thickness m, thickness %)} from the last final layer table
    with open(path, errors='replace') as f:
        lines = f.read().splitlines()
    starts = [i for i, line in enumerate(lines) if 'overall thickness' in line and i + 1 < len(lines)
              and 'target' in lines[i + 1]]
    if not starts:
        return None
    table = {}
    for line in lines[starts[-1] + 3:]:
        parts = line.split()
        if len(parts) != 6:
            break
        try:
            table[parts[0]] = (int(parts[1]), int(parts[2]), float(parts[3]), float(parts[4]), float(parts[5]))
        except ValueError:
            break
    return table


def main(argv):
    args = [a for a in argv]
    name, minFraction = '', MIN_FRACTION
    if '--name' in args:
        k = args.index('--name')
        name = args[k + 1] + ': '
        del args[k:k + 2]
    if '--min' in args:
        k = args.index('--min')
        minFraction = float(args[k + 1])
        del args[k:k + 2]
    if len(args) != 2:
        print('usage: python3 layerReport.py <log.snappyHexMesh> <system/meshGroups> [--name <mesh>] [--min <fraction>]')
        return
    try:
        groups = readMeshGroups(args[1])
    except IOError:
        return
    if not groups:
        return
    try:
        table = readLayerTable(args[0])
    except IOError:
        print('\t%sWARNING: no %s, layer coverage not checked' % (name, args[0]))
        return
    if table is None:
        print('\t%sWARNING: no layer table in %s (addLayers off or snappy stopped early)' % (name, args[0]))
        return
    print('\t%slayers of the mesh groups (target, mean layers in the mesh, thickness of the target %%):' % (name))
    for patch, (group, target) in groups.items():
        if patch not in table:
            print('\t\tWARNING: %s (%s) not in the layer table: patch missing from the mesh?' % (patch, group))
            continue
        faces, _, mesh, thickness, percent = table[patch]
        warn = mesh < minFraction*target or percent < 100.0*minFraction
        print('\t\t%s%-32s %-14s %3d %6.2f %6.1f %%' % ('WARNING: ' if warn else '', patch, group, target, mesh, percent))


if __name__ == '__main__':
    main(sys.argv[1:])
