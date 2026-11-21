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

#remeshing steps for the runScript solve loop (internalSetup MESH_MOTION REMESH auto)
#the case mesh is seg01, the other meshes are full cases in meshes/<seg> made at their start angle
#(system/remeshSegments). a remeshed time directory holds the new mesh (polyMesh incl. points0 and a
#remeshSegment marker), the solver restarts from it and stops at the end angle of that mesh.
#usage, from the case directory:
#  python3 remeshRun.py next               prints "solve <seg> <time>", "map <seg> <time>" or "done"
#  python3 remeshRun.py prepare <seg> <time> sets meshes/<seg>/system/controlDict to start at <time> (mapFieldsPar target)
#  python3 remeshRun.py nearest <seg> <time> renames the fields of the nearest-cell map to <field>_nearest, so the
#                                          volume-weighted map can be written next to them (remeshMassFix -combine)
#  python3 remeshRun.py graft <seg> <time> moves the fields mapped onto meshes/<seg> (mapFieldsPar) and the mesh
#                                          into processor*/<time>, ready for the restart
#stdlib only so it runs on cluster nodes

import os
import re
import sys
import glob
import gzip
import shutil


stateFile = '.remeshState'
markerName = 'remeshSegment'


def readSegments(caseDir='.'):
    #[(name, start angle, end angle)], angles in deg mod 360
    path = os.path.join(caseDir, 'system', 'remeshSegments')
    if not os.path.isfile(path):
        sys.exit('ERROR! %s not found: the case was not set up for remeshing' % (path))
    segments = []
    with open(path) as f:
        for line in f:
            words = line.split()
            if not words or words[0].startswith('#'):
                continue
            segments.append((words[0], float(words[1]), float(words[2])))
    if not segments:
        sys.exit('ERROR! %s lists no meshes' % (path))
    return segments


def readEntry(path, key):
    #value of a top-level "key value;" entry of an OpenFOAM dictionary, None if missing
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        text = f.read()
    m = re.search(r'^\s*%s\s+([^;]+);' % (re.escape(key)), text, re.MULTILINE)
    return m.group(1).strip() if m else None


def timeDirs(procDir):
    #numeric time directories, sorted: [(value, name)]
    times = []
    for name in os.listdir(procDir):
        if not os.path.isdir(os.path.join(procDir, name)):
            continue
        try:
            times.append((float(name), name))
        except ValueError:
            continue
    return sorted(times)


def findTime(procDir, value):
    #name of the time directory matching value (a name written with fewer digits still matches)
    for t, name in timeDirs(procDir):
        if abs(t - value) <= 1e-12 + 1e-6*abs(value):
            return name
    return None


def prepare(seg, timeName):
    #mapFieldsPar writes into the target's start time: the exact time name, with enough digits for the directory name
    #(foamDictionary -set rounds to 6 digits)
    path = os.path.join('meshes', seg, 'system', 'controlDict')
    with open(path) as f:
        text = f.read()
    for key, value in (('startFrom', 'startTime'), ('startTime', timeName), ('timePrecision', '15')):
        entry = r'^(\s*%s\s+)[^;]+;' % (key)
        if re.search(entry, text, flags=re.MULTILINE):
            text = re.sub(entry, r'\g<1>%s;' % (value), text, flags=re.MULTILINE)
        else:
            text = text.rstrip() + '\n%s %s;\n' % (key, value)
    with open(path, 'w') as f:
        f.write(text)


def segmentOf(theta, segments, tol=1e-6):
    #mesh whose crank angle range holds theta; a crank angle a hair before a segment start counts as that segment
    for name, a0, a1 in segments:
        span = (a1 - a0) % 360.0 or 360.0
        if (theta - a0 + tol) % 360.0 < span:
            return name
    sys.exit('ERROR! crank angle %.6g deg is in no remeshing segment' % (theta))


def meshSegment(procDir, timeValue):
    #mesh of a time: the latest remeshed time directory up to it (polyMesh/faces + marker), else the case mesh
    for t, name in reversed(timeDirs(procDir)):
        if t > timeValue + 1e-12 + 1e-9*abs(timeValue):
            continue
        meshDir = os.path.join(procDir, name, 'polyMesh')
        if os.path.isfile(os.path.join(meshDir, 'faces')) or os.path.isfile(os.path.join(meshDir, 'faces.gz')):
            marker = os.path.join(meshDir, markerName)
            if not os.path.isfile(marker):
                sys.exit('ERROR! %s has a mesh but no %s marker (not written by remeshRun.py)' % (meshDir, markerName))
            with open(marker) as f:
                return f.read().split()[0]
    return 'seg01'


def crankAngle(procDir, timeName):
    #accumulated crank angle (deg) written by the crank motion solver, None before the first write
    value = readEntry(os.path.join(procDir, timeName, 'uniform', 'crankState'), 'theta')
    return float(value) if value is not None else None


def useDynamicMeshDict(seg):
    shutil.copyfile(os.path.join('constant', 'dynamicMeshDict.%s' % (seg)), os.path.join('constant', 'dynamicMeshDict'))


def nextStep():
    segments = readSegments()
    procDirs = sorted(glob.glob('processor[0-9]*'))
    if not procDirs:
        sys.exit('ERROR! no processor directories: the case is not decomposed')
    procDir = procDirs[0]
    controlDict = os.path.join('system', 'controlDict')
    endTime = float(readEntry(controlDict, 'endTime'))
    tValue, tName = timeDirs(procDir)[-1]
    #the solver writes its last step and stops within one step of endTime (crank solver, remeshing)
    lastStep = float(readEntry(controlDict, 'maxDeltaT') or readEntry(controlDict, 'deltaT') or 0.0)
    if tValue >= endTime - lastStep - 1e-12:
        return 'done'
    theta = crankAngle(procDir, tName)
    if theta is None:
        #no crank state yet (start of the run): the start angle of the case mesh
        theta = segments[0][1]
    seg = segmentOf(theta, segments)
    mesh = meshSegment(procDir, tValue)
    if seg != mesh:
        return 'map %s %s' % (seg, tName)
    #the solver has to get past this time, otherwise the loop would run it again forever
    last = open(stateFile).read().split() if os.path.isfile(stateFile) else []
    if last[:2] == ['solve', tName]:
        sys.exit('ERROR! the solver stopped at t = %s again without reaching the next mesh, see the solver log' % (tName))
    with open(stateFile, 'w') as f:
        f.write('solve %s %s\n' % (tName, seg))
    useDynamicMeshDict(seg)
    return 'solve %s %s' % (seg, tName)


def renameField(path, newName):
    #field file under a new name: the header object entry has to match the file name
    with open(path, 'rb') as f:
        data = f.read()
    head = re.sub(rb'(object\s+)[^;\s]+;', b'\g<1>' + newName.encode() + b';', data[:2000], count=1)
    with open(os.path.join(os.path.dirname(path), newName), 'wb') as f:
        f.write(head + data[2000:])
    os.remove(path)


def stashNearest(seg, timeName):
    #fields of the nearest-cell map -> <field>_nearest
    segDir = os.path.join('meshes', seg)
    count = 0
    for procDir in sorted(glob.glob(os.path.join(segDir, 'processor[0-9]*'))):
        name = findTime(procDir, float(timeName))
        if name is None:
            sys.exit('ERROR! %s has no mapped time %s (mapFieldsPar failed?)' % (procDir, timeName))
        timeDir = os.path.join(procDir, name)
        for entry in os.listdir(timeDir):
            path = os.path.join(timeDir, entry)
            if os.path.isfile(path) and not entry.endswith('_nearest'):
                renameField(path, entry + '_nearest')
                count += 1
    print('%d nearest-map fields renamed in %s' % (count, segDir))


def writePoints0(meshDir):
    #points0 of the new mesh = its points (the motion solver reads points0 from the latest time that has one)
    for name, opener in (('points', open), ('points.gz', gzip.open)):
        src = os.path.join(meshDir, name)
        if not os.path.isfile(src):
            continue
        with opener(src, 'rb') as f:
            data = f.read()
        head = re.sub(rb'(object\s+)points;', rb'\1points0;', data[:2000], count=1)
        with opener(os.path.join(meshDir, name.replace('points', 'points0')), 'wb') as f:
            f.write(head + data[2000:])
        return
    sys.exit('ERROR! no points file in %s' % (meshDir))


def setDeltaT(timeDir, deltaT):
    #restart time step (uniform/time), the adjustable time step ramps it up again
    path = os.path.join(timeDir, 'uniform', 'time')
    if not os.path.isfile(path):
        return
    with open(path) as f:
        text = f.read()
    text = re.sub(r'^(\s*deltaT0?\s+)[^;]+;', r'\g<1>%s;' % (deltaT), text, flags=re.MULTILINE)
    with open(path, 'w') as f:
        f.write(text)


def graft(seg, timeName):
    segDir = os.path.join('meshes', seg)
    procDirs = sorted(glob.glob('processor[0-9]*'))
    #restart with a small time step after the map (INIT_DT, the controlDict deltaT) when the time step is adjustable
    adjust = (readEntry(os.path.join('system', 'controlDict'), 'adjustTimeStep') or 'no').lower() in ('yes', 'on', 'true')
    initDt = readEntry(os.path.join('system', 'controlDict'), 'deltaT') if adjust else None
    value = float(timeName)
    for procDir in procDirs:
        dst = os.path.join(procDir, timeName)
        srcName = findTime(os.path.join(segDir, procDir), value)
        if srcName is None:
            sys.exit('ERROR! %s has no mapped time %s (mapFieldsPar failed?)' % (os.path.join(segDir, procDir), timeName))
        src = os.path.join(segDir, procDir, srcName)
        #fields of the old mesh go (surface fluxes cannot be mapped), uniform/ (time, crank state) stays
        for entry in os.listdir(dst):
            if entry == 'uniform':
                continue
            path = os.path.join(dst, entry)
            if os.path.isdir(path) and not os.path.islink(path):
                shutil.rmtree(path)
            else:
                os.remove(path)
        for entry in os.listdir(src):
            if entry != 'uniform':
                shutil.move(os.path.join(src, entry), os.path.join(dst, entry))
        shutil.rmtree(src)
        meshDir = os.path.join(dst, 'polyMesh')
        shutil.copytree(os.path.join(segDir, procDir, 'constant', 'polyMesh'), meshDir)
        writePoints0(meshDir)
        with open(os.path.join(meshDir, markerName), 'w') as f:
            f.write('%s\n' % (seg))
        shutil.copyfile(os.path.join('system', 'pointDisplacement'), os.path.join(dst, 'pointDisplacement'))
        if initDt:
            setDeltaT(dst, initDt)
    useDynamicMeshDict(seg)
    print('%s mesh and mapped fields in processor*/%s%s' % (seg, timeName,
                                                         ', restart time step %s' % (initDt) if initDt else ''))


def main(argv):
    if len(argv) == 1 and argv[0] == 'next':
        print(nextStep())
    elif len(argv) == 3 and argv[0] == 'prepare':
        prepare(argv[1], argv[2])
    elif len(argv) == 3 and argv[0] == 'nearest':
        stashNearest(argv[1], argv[2])
    elif len(argv) == 3 and argv[0] == 'graft':
        graft(argv[1], argv[2])
    else:
        sys.exit('usage: remeshRun.py next | prepare <seg> <time> | nearest <seg> <time> | graft <seg> <time>')


if __name__ == '__main__':
    main(sys.argv[1:])
