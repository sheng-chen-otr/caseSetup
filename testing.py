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

"""Scratch utility: print the bounding box of a surface file.

The OBJ/STL parsing now lives in the shared, hardened ``geometryIO`` module.
"""
import sys

import geometryIO


def getBoundingBoxOBJ(geomFile):
    objPath = 'constant/triSurface/%s' % (geomFile)
    return geometryIO.surfaceBoundingBox(objPath)


def getBoundingBoxSTL(geomFile):
    stlPath = 'constant/triSurface/%s' % (geomFile)
    return geometryIO.surfaceBoundingBox(stlPath)


def main():
    geomFile = sys.argv[1] if len(sys.argv) > 1 else 'ROTA-fr-wh-lhs.stl'
    bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ = getBoundingBoxSTL(geomFile)
    print(bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ)


if __name__ == '__main__':
    main()