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