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

#writers for internalSetup.py (internal flows / Stirling engines)
#all paths are relative to the case directory (the tool runs from inside the case)

import os
import sys
import json
import math

import numpy as np

import geometryIO
import motionUtilities
from setupCore import stripExt, __version__

#gas presets: molWeight, Cp, mu at Tref, Pr, Sutherland mu0/T0/S (As is worked out from them)
#approximate property values, override any of them in GLOBAL_MATERIAL
gasPresets = {
    'air':      {'MOL_WEIGHT': 28.96, 'CP': 1005.0,  'MU': 1.81e-5, 'PR': 0.71, 'MU0': 1.716e-5, 'T0': 273.15, 'S': 110.4},
    'nitrogen': {'MOL_WEIGHT': 28.01, 'CP': 1040.0,  'MU': 1.76e-5, 'PR': 0.71, 'MU0': 1.663e-5, 'T0': 273.15, 'S': 107.0},
    'helium':   {'MOL_WEIGHT': 4.003, 'CP': 5193.0,  'MU': 1.99e-5, 'PR': 0.67, 'MU0': 1.870e-5, 'T0': 273.15, 'S': 79.4},
    'hydrogen': {'MOL_WEIGHT': 2.016, 'CP': 14300.0, 'MU': 8.90e-6, 'PR': 0.69, 'MU0': 8.411e-6, 'T0': 273.15, 'S': 72.0},
}

#T options for BC groups: name -> number of arguments
tModes = {'adiabatic': 0, 'fixed': 1, 'convective': 2, 'flux': 1}
uModes = ['noSlip', 'slip']


def getVal(fullCaseSetupDict, section, key):
    return fullCaseSetupDict[section][key][0]


def copyTemplate(templateLoc, templateName, dest, replaceDict):
    #reads defaultDicts/<templateName>, replaces <KEY> placeholders, writes dest
    with open('%s/defaultDicts/%s' % (templateLoc, templateName)) as f:
        text = f.read()
    for key, val in replaceDict.items():
        text = text.replace('<%s>' % key, str(val))
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, 'w') as f:
        f.write(text)


#----- material

def getGasProperties(fullCaseSetupDict):
    #preset values for GAS, any key not 'default' overrides the preset
    material = fullCaseSetupDict['GLOBAL_MATERIAL']
    gas = material['GAS'][0].lower()
    if gas not in gasPresets and gas != 'custom':
        sys.exit('ERROR! GLOBAL_MATERIAL GAS must be one of %s or custom!' % (list(gasPresets.keys())))
    preset = gasPresets.get(gas, {})
    props = {}
    for key in ['MOL_WEIGHT', 'CP', 'MU', 'PR']:
        val = material[key][0]
        if val.lower() == 'default':
            if key not in preset:
                sys.exit('ERROR! GAS = custom needs GLOBAL_MATERIAL %s!' % (key))
            props[key] = preset[key]
        else:
            props[key] = float(val)
    #sutherland As from mu0*(T0+S)/T0^1.5 unless given
    aS, tS = material['SUTHERLAND_AS'][0], material['SUTHERLAND_TS'][0]
    if tS.lower() == 'default':
        props['TS'] = preset.get('S', 110.4)
    else:
        props['TS'] = float(tS)
    if aS.lower() == 'default':
        if gas == 'custom':
            props['AS'] = props['MU']*(300.0 + props['TS'])/300.0**1.5
        else:
            props['AS'] = preset['MU0']*(preset['T0'] + preset['S'])/preset['T0']**1.5
    else:
        props['AS'] = float(aS)
    transport = material['TRANSPORT'][0].lower()
    if transport not in ('const', 'sutherland'):
        sys.exit('ERROR! GLOBAL_MATERIAL TRANSPORT must be const or sutherland!')
    props['TRANSPORT'] = transport
    return props


#solver -> (flow, steady, buoyant)
solverTypes = {'rhoPimpleFoam': ('compressible', False, False), 'pimpleFoam': ('incompressible', False, False),
               'rhoSimpleFoam': ('compressible', True, False), 'buoyantSimpleFoam': ('compressible', True, True),
               'simpleFoam': ('incompressible', True, False)}


def getSolverType(fullCaseSetupDict):
    solver = getVal(fullCaseSetupDict, 'GLOBAL_SIM_CONTROL', 'SOLVER')
    if solver not in solverTypes:
        sys.exit('ERROR! GLOBAL_SIM_CONTROL SOLVER must be one of %s, got %s!' % (list(solverTypes.keys()), solver))
    return solverTypes[solver]


def getFlow(fullCaseSetupDict):
    #compressible (rhoPimpleFoam, rhoSimpleFoam, buoyantSimpleFoam) or incompressible (pimpleFoam, simpleFoam:
    #kinematic p, no energy equation)
    return getSolverType(fullCaseSetupDict)[0]


def isSteady(fullCaseSetupDict):
    return getSolverType(fullCaseSetupDict)[1]


def isBuoyant(fullCaseSetupDict):
    return getSolverType(fullCaseSetupDict)[2]


def getGasDensity(fullCaseSetupDict):
    #ideal gas density at CHARGE_PRESSURE and INIT_T
    props = getGasProperties(fullCaseSetupDict)
    material = fullCaseSetupDict['GLOBAL_MATERIAL']
    return float(material['CHARGE_PRESSURE'][0])*props['MOL_WEIGHT']/(8314.46*float(material['INIT_T'][0]))


def getNu(fullCaseSetupDict):
    #kinematic viscosity for incompressible runs: NU if given, else mu/rho of the gas at the charge state
    nu = fullCaseSetupDict['GLOBAL_MATERIAL'].get('NU', ['default'])[0]
    if nu.lower() != 'default':
        return float(nu)
    return getGasProperties(fullCaseSetupDict)['MU']/getGasDensity(fullCaseSetupDict)


#----- geometry PIDs and BC groups

def readGeomPIDs(geomDict, mshDir):
    #pidDict[geomName] = [pid, ...] read from the geometry files in mshDir (no OpenFOAM needed)
    pidDict = {}
    for geom in geomDict:
        geomPath = os.path.join(mshDir, geom)
        if not os.path.isfile(geomPath):
            sys.exit('ERROR! Geometry %s not found in %s!' % (geom, mshDir))
        pids = geometryIO.readSurfacePIDs(geomPath)
        if len(pids) == 0 or any(pid == '' for pid in pids):
            sys.exit('ERROR! %s has unnamed PIDs (STL "solid" without a name / OBJ without groups)!' % (geom))
        pidDict[stripExt(geom)] = pids
    return pidDict


def parseThermal(tokens, where):
    #T tokens: adiabatic | fixed <T> | convective <h> <Ta> | flux <q>
    if len(tokens) == 0:
        sys.exit('ERROR! %s: T needs a mode (%s)!' % (where, ' | '.join(tModes)))
    mode = tokens[0]
    if mode not in tModes:
        sys.exit('ERROR! %s: T mode %s is invalid, use one of %s!' % (where, mode, list(tModes.keys())))
    args = tokens[1:]
    if len(args) != tModes[mode]:
        sys.exit('ERROR! %s: T %s needs %d value(s), got %s!' % (where, mode, tModes[mode], args))
    try:
        args = [float(a) for a in args]
    except ValueError:
        sys.exit('ERROR! %s: T %s values must be numbers, got %s!' % (where, mode, tokens[1:]))
    return [mode] + args


def parseNumbers(tokens, where, what):
    try:
        return [float(t) for t in tokens]
    except ValueError:
        sys.exit('ERROR! %s: %s needs numbers, got %s!' % (where, what, ' '.join(tokens)))


def parseFlowSpec(uTokens, pTokens, where):
    #inlet/outlet from the U and P tokens, None for a wall
    #U massFlow <kg/s> | U velocity <m/s> (normal, into the domain) | U velocity <ux> <uy> <uz>
    #P total <Pa> (inlet) | P fixed <Pa> (outlet, static pressure)
    if uTokens and uTokens[0] in uModes:
        if len(uTokens) != 1 or pTokens:
            sys.exit('ERROR! %s: a wall (U %s) takes no other U values and no P!' % (where, uTokens[0]))
        return None
    if uTokens and pTokens:
        sys.exit('ERROR! %s: give either U (mass flow / velocity inlet) or P (pressure inlet / outlet), not both!' % (where))
    if uTokens:
        mode, values = uTokens[0], parseNumbers(uTokens[1:], where, 'U %s' % (uTokens[0]))
        if mode == 'massFlow' and len(values) == 1:
            if values[0] <= 0:
                sys.exit('ERROR! %s: U massFlow must be positive (into the domain)!' % (where))
            return {'kind': 'inlet', 'type': 'massFlow', 'value': values[0]}
        if mode == 'velocity' and len(values) == 1:
            return {'kind': 'inlet', 'type': 'velocity', 'value': values[0]}
        if mode == 'velocity' and len(values) == 3:
            return {'kind': 'inlet', 'type': 'velocityVector', 'value': values}
        sys.exit('ERROR! %s: U must be %s, massFlow <kg/s>, velocity <m/s> or velocity <ux> <uy> <uz>!'
                 % (where, ' | '.join(uModes)))
    if pTokens:
        values = parseNumbers(pTokens[1:], where, 'P %s' % (pTokens[0]))
        if len(values) == 1 and pTokens[0] == 'total':
            return {'kind': 'inlet', 'type': 'totalPressure', 'value': values[0]}
        if len(values) == 1 and pTokens[0] in ('fixed', 'static'):
            return {'kind': 'outlet', 'type': 'staticPressure', 'value': values[0]}
        sys.exit('ERROR! %s: P must be total <Pa> (inlet) or fixed <Pa> (outlet)!' % (where))
    return None


def parseGroupLine(tokens, where):
    #<group> = PIDS kw1 kw2 .. [U noSlip|slip] [T mode args]         wall
    #<group> = PIDS kw1 .. U massFlow|velocity .. | P total|fixed .. [T fixed <T>]   inlet / outlet
    tokens = [t for t in tokens if t != '']
    parts = {'PIDS': [], 'U': [], 'T': [], 'P': []}
    current = None
    for token in tokens:
        if token in parts:
            current = token
        elif current is None:
            sys.exit('ERROR! %s must start with PIDS, got: %s' % (where, ' '.join(tokens)))
        else:
            parts[current].append(token)
    if len(parts['PIDS']) == 0:
        sys.exit('ERROR! %s has no PIDS keywords!' % (where))
    group = {'keywords': parts['PIDS'], 'U': None, 'T': None, 'flow': parseFlowSpec(parts['U'], parts['P'], where)}
    if parts['U'] and group['flow'] is None:
        group['U'] = parts['U'][0]
    if parts['T']:
        group['T'] = parseThermal(parts['T'], where)
        if group['flow'] and group['T'][0] != 'fixed':
            sys.exit('ERROR! %s: inlets and outlets take T fixed <T> (inflow / backflow temperature)!' % (where))
    return group


def resolveBCGroups(wallName, wallSetup, pids):
    #matches the BC_GROUPS keywords (case-insensitive substring) to the PIDs of one WALL geometry
    #returns an ordered list of groups {name, pids, U, T}; unmatched PIDs go to <wall>_default
    uDefault = wallSetup['U_DEFAULT'][0]
    if uDefault not in uModes:
        sys.exit('ERROR! [%s] U_DEFAULT must be one of %s!' % (wallName, uModes))
    tDefault = parseThermal([t for t in wallSetup['T_DEFAULT'] if t != ''], '[%s] T_DEFAULT' % (wallName))
    groupNames = [g for g in wallSetup.get('BC_GROUPS', ['']) if g != '']
    groups = []
    owner = {}
    for groupName in groupNames:
        if groupName not in wallSetup:
            sys.exit('ERROR! [%s] BC_GROUPS lists %s but there is no "%s = PIDS ..." line!' % (wallName, groupName, groupName))
        group = parseGroupLine(wallSetup[groupName], '[%s] %s' % (wallName, groupName))
        matched = []
        for kw in group['keywords']:
            hits = [pid for pid in pids if kw.lower() in pid.lower()]
            if len(hits) == 0:
                sys.exit('ERROR! [%s] %s: keyword "%s" matches no PID. PIDs are: %s' % (wallName, groupName, kw, pids))
            for pid in hits:
                if pid in owner and owner[pid] != groupName:
                    sys.exit('ERROR! [%s] PID "%s" matches groups %s and %s!' % (wallName, pid, owner[pid], groupName))
                owner[pid] = groupName
                if pid not in matched:
                    matched.append(pid)
        if group['flow']:
            #inlet/outlet: T is the inflow (backflow) temperature, None = INIT_T
            groups.append({'name': groupName, 'pids': matched, 'U': None, 'T': group['T'], 'flow': group['flow']})
        else:
            groups.append({'name': groupName, 'pids': matched, 'flow': None,
                           'U': group['U'] or uDefault, 'T': group['T'] or tDefault})
    rest = [pid for pid in pids if pid not in owner]
    if len(rest) > 0:
        groups.append({'name': '%s_default' % (wallName), 'pids': rest, 'U': uDefault, 'T': tDefault, 'flow': None})
    return groups


def flowLabel(flow):
    #short text of an inlet/outlet for the BC table
    if flow['type'] == 'velocityVector':
        return 'velocity (%s)' % (' '.join('%g' % v for v in flow['value']))
    return '%s %g' % (flow['type'], flow['value'])


def getWallGroups(geomDict, fullCaseSetupDict, pidDict, movingPids=None):
    #wallGroups[geomName] = groups for WALL and MOVE geometry, prints the PID -> group -> U, T table
    #movingPids[geomName] = PIDs that move with a body (U becomes movingWallVelocity)
    movingPids = movingPids or {}
    wallGroups = {}
    print('\n\tBoundary condition groups:')
    print('\t\t%-32s %-20s %-14s %s' % ('patch', 'group', 'U', 'T'))
    for geom in geomDict:
        geomName = stripExt(geom)
        if not geomName.startswith(('WALL', 'MOVE')):
            continue
        groups = resolveBCGroups(geomName, fullCaseSetupDict[geomName], pidDict[geomName])
        moving = movingPids.get(geomName, [])
        for group in groups:
            group['moving'] = [pid for pid in group['pids'] if pid in moving]
            if group['flow'] and group['moving']:
                sys.exit('ERROR! [%s] %s is an %s but %s move with a body!' % (geomName, group['name'],
                                                                             group['flow']['kind'], group['moving']))
            for pid in group['pids']:
                if group['flow']:
                    uLabel = '%s: %s' % (group['flow']['kind'], flowLabel(group['flow']))
                    tLabel = 'T %g' % (group['T'][1]) if group['T'] else 'T INIT_T'
                else:
                    uLabel = 'movingWall' if pid in moving else group['U']
                    tLabel = ' '.join(str(t) for t in group['T'])
                print('\t\t%-32s %-20s %-14s %s' % ('%s_%s' % (geomName, pid), group['name'], uLabel, tLabel))
        wallGroups[geomName] = groups
    return wallGroups


#----- system dicts

def getBackgroundMesh(bbox, baseCell, nMargin=2):
    #bbox (xmin ymin zmin xmax ymax zmax) grown by nMargin cells and rounded to whole cells
    lo = np.floor(np.array(bbox[:3])/baseCell - nMargin)*baseCell
    hi = np.ceil(np.array(bbox[3:])/baseCell + nMargin)*baseCell
    n = np.rint((hi - lo)/baseCell).astype(int)
    return lo, hi, n


def caseGeomPath(geom):
    #geometry file as snappy reads it: shifted copies are written without .gz
    plain = os.path.join('constant', 'triSurface', geom.replace('.gz', ''))
    return plain if os.path.isfile(plain) else os.path.join('constant', 'triSurface', geom)


def writeBlockMeshInternal(templateLoc, geomDict, fullCaseSetupDict):
    print('\n\tWriting blockMeshDict...')
    baseCell = float(getVal(fullCaseSetupDict, 'GLOBAL_REFINEMENT', 'BASE_CELL'))
    allBox = []
    for geom in geomDict:
        scale = float(geomDict[geom]['scale'])
        allBox.append(np.array(geometryIO.surfaceBoundingBox(caseGeomPath(geom)))*scale)
    allBox = np.array(allBox)
    bbox = list(allBox[:, :3].min(axis=0)) + list(allBox[:, 3:].max(axis=0))
    lo, hi, n = getBackgroundMesh(bbox, baseCell)
    print('\t\tbackground box %s to %s, cells %s (base cell %g m)' % (lo, hi, n, baseCell))
    replaceDict = {'XMIN': '%.9g' % lo[0], 'YMIN': '%.9g' % lo[1], 'ZMIN': '%.9g' % lo[2],
                   'XMAX': '%.9g' % hi[0], 'YMAX': '%.9g' % hi[1], 'ZMAX': '%.9g' % hi[2],
                   'NX': n[0], 'NY': n[1], 'NZ': n[2]}
    copyTemplate(templateLoc, 'system/blockMeshDict', 'system/blockMeshDict', replaceDict)
    return bbox


def writeSnappyInternal(templateLoc, geomDict, fullCaseSetupDict, gapRefs=None):
    print('\n\tWriting snappyHexMeshDict...')
    refinement = fullCaseSetupDict['GLOBAL_REFINEMENT']
    locInMesh = [v for v in refinement['LOC_IN_MESH'] if v != '']
    if len(locInMesh) != 3:
        sys.exit('ERROR! GLOBAL_REFINEMENT LOC_IN_MESH must be 3 coordinates of a point inside the gas!')
    try:
        [float(v) for v in locInMesh]
    except ValueError:
        sys.exit('ERROR! GLOBAL_REFINEMENT LOC_IN_MESH must be 3 numbers, got %s!' % (locInMesh))
    #feature edges: surfaceFeatureExtract writes constant/triSurface/<geom>.eMesh (run before snappy),
    #refined at each geometry's surface level (FEATURE_LEVEL surface) or a fixed level
    featureLevel = refinement.get('FEATURE_LEVEL', ['surface'])[0]
    featureSnap = refinement.get('FEATURE_SNAP', ['explicit'])[0].lower()
    if featureSnap not in ('explicit', 'implicit'):
        sys.exit('ERROR! GLOBAL_REFINEMENT FEATURE_SNAP must be explicit or implicit!')
    if featureLevel != 'surface' and not featureLevel.isdigit():
        sys.exit('ERROR! GLOBAL_REFINEMENT FEATURE_LEVEL must be surface or a refinement level, got %s!' % (featureLevel))
    featAngle = refinement.get('FEAT_EDGE_ANGLE', ['150'])[0]
    geometryArray, surfaceArray, layerArray, featureArray, extractArray = [], [], [], [], []
    for geom in geomDict:
        geomName = stripExt(geom)
        geomFile = geom.replace('.gz', '')
        level = geomDict[geom]['refinement']
        geometryArray.append('    %s { type triSurfaceMesh; name %s; file "%s"; scale %s; }'
                             % (geomFile, geomName, geomFile, geomDict[geom]['scale']))
        featureArray.append('        { file "%s.eMesh"; scale %s; level %s; }'
                            % (geomName, geomDict[geom]['scale'], level if featureLevel == 'surface' else featureLevel))
        extractArray.append('%s { extractionMethod extractFromSurface; includedAngle %s; '
                            'subsetFeatures { nonManifoldEdges no; openEdges yes; } writeObj no; }' % (geomFile, featAngle))
        #fan disk: internal faceZone, made into a cyclic fan baffle by createBaffles after the zero directory
        if geomName.startswith('FAN'):
            surfaceArray.append('        %s { level (%s %s); faceZone %s; faceType internal; }'
                                % (geomName, level, level, geomName))
            continue
        #porous zone: the cells inside the closed POR surface become cellZone POR-<name> (no patch)
        if geomName.startswith('POR'):
            surfaceArray.append('        %s { level (%s %s); faceZone %s; cellZone %s; cellZoneInside inside; }'
                                % (geomName, level, level, geomName, geomName))
            continue
        #per-region levels by PID name (regex) from system/snappyRefinementDict
        surfaceArray.append('        %s { level (%s %s); patchInfo { type wall; } regions { #include "snappyRefinementDict" } }'
                            % (geomName, level, level))
        layerArray.append('        "%s_.*" { nSurfaceLayers %s; expansionRatio %s;}'
                          % (geomName, geomDict[geom]['layers'], geomDict[geom]['expansion']))
    #gap faces: no layers, written after the file-wide entry (last matching regex wins)
    for geom in geomDict:
        geomName = stripExt(geom)
        gapPids = [p for p in fullCaseSetupDict[geomName].get('GAP_PIDS', ['']) if p != '']
        for kw in gapPids:
            layerArray.append('        "%s_.*%s.*" { nSurfaceLayers 0; }' % (geomName, kw))
    #gap refinement: distance band around a per-body STL of the gap faces
    regionArray = []
    for gapRef in (gapRefs or []):
        geometryArray.append('    %s { type triSurfaceMesh; name %s; file "%s"; }'
                             % (gapRef['file'], gapRef['name'], gapRef['file']))
        regionArray.append('        %s { mode distance; levels ((%.6g %d)); }'
                           % (gapRef['name'], gapRef['distance'], gapRef['level']))
    replaceDict = {'GEOMETRY': '\n'.join(geometryArray),
                   'REFINEMENT_SURFACES': '\n'.join(surfaceArray),
                   'REFINEMENT_REGIONS': '\n'.join(regionArray),
                   'LAYERS': '\n'.join(layerArray),
                   'LOC_IN_MESH': ' '.join(locInMesh),
                   'RESOLVE_FEAT_ANGLE': refinement['RESOLVE_FEAT_ANGLE'][0],
                   'N_CELLS_BETWEEN_LEVELS': refinement['N_CELLS_BETWEEN_LEVELS'][0],
                   'FINAL_LAYER_THICKNESS': refinement['FINAL_LAYER_THICKNESS'][0],
                   'MIN_LAYER_THICKNESS': refinement['MIN_LAYER_THICKNESS'][0],
                   'FEATURES': '\n'.join(featureArray),
                   'EXPLICIT_SNAP': 'true' if featureSnap == 'explicit' else 'false',
                   'IMPLICIT_SNAP': 'false' if featureSnap == 'explicit' else 'true'}
    print('\t\tfeature edges: includedAngle %s, level %s, %s snapping' % (featAngle, featureLevel, featureSnap))
    copyTemplate(templateLoc, 'system/snappyHexMeshDict', 'system/snappyHexMeshDict', replaceDict)
    copyTemplate(templateLoc, 'system/surfaceFeatureExtractDict', 'system/surfaceFeatureExtractDict',
                 {'SURFACE_FEAT_EXTRACT': '\n'.join(extractArray)})
    copyTemplate(templateLoc, 'system/snappyRefinementDict', 'system/snappyRefinementDict', {})


def getTurbulence(fullCaseSetupDict):
    #laminar (default) or a RAS model name
    model = getVal(fullCaseSetupDict, 'GLOBAL_SIM_CONTROL', 'TURB_MODEL')
    aliases = {'laminar': 'laminar', 'kosst': 'kOmegaSST', 'komegasst': 'kOmegaSST',
               'kepsilon': 'kEpsilon', 'komega': 'kOmega'}
    if model.lower() not in aliases:
        sys.exit('ERROR! TURB_MODEL %s is not supported for internal cases, use one of %s!'
                 % (model, sorted(set(aliases.values()))))
    return aliases[model.lower()]


def pidArea(geom, geomDict, pids):
    #area (m2) of some PIDs of a geometry, from the case copy in constant/triSurface
    scale = float(geomDict[geom]['scale'])
    tris = geometryIO.readSurfaceTriangles(caseGeomPath(geom))
    area = 0.0
    for pid in pids:
        t = np.asarray(tris[pid])*scale
        area += 0.5*float(np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1).sum())
    return area


def inletTurbulence(fullCaseSetupDict, speed, area):
    #k, omega, epsilon at an inlet from TURB_INTENSITY and TURB_LENGTH (auto = 0.07 x hydraulic diameter of a
    #round inlet of the same area)
    simControl = fullCaseSetupDict['GLOBAL_SIM_CONTROL']
    intensity = float(simControl.get('TURB_INTENSITY', ['0.05'])[0])
    length = simControl.get('TURB_LENGTH', ['auto'])[0]
    length = 0.07*np.sqrt(4.0*area/np.pi) if length == 'auto' else float(length)
    k = 1.5*(intensity*speed)**2
    omega = np.sqrt(k)/(0.09**0.25*length)
    epsilon = 0.09**0.75*k**1.5/length
    return k, omega, epsilon, length


def flowEntry(geom, geomDict, group, fullCaseSetupDict, turbModel):
    #caseProperties entry of an inlet or outlet group
    #pressures in Pa: absolute for compressible runs, minus CHARGE_PRESSURE and over the gas density
    #(kinematic) for incompressible runs
    geomName = stripExt(geom)
    flow = group['flow']
    flowType = getFlow(fullCaseSetupDict)
    material = fullCaseSetupDict['GLOBAL_MATERIAL']
    pCharge = float(material['CHARGE_PRESSURE'][0])
    tFlow = group['T'][1] if group['T'] else float(material['INIT_T'][0])
    rhoCharge = getGasDensity(fullCaseSetupDict)
    if flowType == 'compressible':
        rhoFlow = pCharge*getGasProperties(fullCaseSetupDict)['MOL_WEIGHT']/(8314.46*tFlow)
    else:
        rhoFlow = rhoCharge

    def pressure(p):
        return p if flowType == 'compressible' else (p - pCharge)/rhoCharge

    patches = ' '.join('"%s_%s"' % (geomName, pid) for pid in group['pids'])
    area = pidArea(geom, geomDict, group['pids'])
    values = ['$:initialConditions;']
    speed = None
    if flow['type'] == 'massFlow':
        speed = flow['value']/(rhoFlow*area)
        if flowType == 'compressible':
            option, values = 'flowRate', values + ['massFlowRate constant %.9g;' % (flow['value'])]
        else:
            option, values = 'volumeFlowRate', values + ['volumetricFlowRate constant %.9g;' % (flow['value']/rhoFlow)]
    elif flow['type'] == 'velocity':
        speed = abs(flow['value'])
        option, values = 'normalVelocity', values + ['UnRef uniform %.9g;' % (-flow['value'])]
    elif flow['type'] == 'velocityVector':
        speed = float(np.linalg.norm(flow['value']))
        option, values = 'fixedVelocity', values + ['U uniform (%s);' % (' '.join('%.9g' % v for v in flow['value']))]
    else:
        option = 'fixedTotalPressure' if flow['type'] == 'totalPressure' else 'staticPressure'
        values = values + ['p uniform %.9g;' % (pressure(flow['value']))]
    if flowType == 'compressible':
        values.append('T uniform %.9g;' % (tFlow))
    text = '%s, area %.4g m2' % (flowLabel(flow), area)
    if speed is not None:
        text += ', mean velocity %.4g m/s' % (speed)
    if turbModel != 'laminar' and speed:
        k, omega, epsilon, length = inletTurbulence(fullCaseSetupDict, speed, area)
        values += ['k uniform %.6g;' % (k), 'omega uniform %.6g;' % (omega), 'epsilon uniform %.6g;' % (epsilon)]
        text += ', k %.3g omega %.3g (length %.3g m)' % (k, omega, length)
    print('\t\t%s %s: %s' % (flow['kind'], group['name'], text))
    key = 'flowSpecification' if flow['kind'] == 'inlet' else 'returnFlow'
    return ('    %s { category %s; type subSonic; patches (%s); options { %s %s; } values { %s } }'
            % (group['name'], flow['kind'], patches, key, option, ' '.join(values)))


def writeCasePropertiesInternal(templateLoc, geomDict, fullCaseSetupDict, wallGroups):
    print('\n\tWriting caseProperties...')
    material = fullCaseSetupDict['GLOBAL_MATERIAL']
    pInit = float(material['CHARGE_PRESSURE'][0])
    tInit = float(material['INIT_T'][0])
    #incompressible: kinematic p relative to pRefValue 0, no T (wall thermal settings are ignored)
    incompressible = getFlow(fullCaseSetupDict) == 'incompressible'
    if incompressible:
        initArray = ['U uniform (0 0 0);', 'p uniform 0;']
        thermal = sorted(g['name'] for name in wallGroups for g in wallGroups[name] if g['T'] and g['T'][0] != 'adiabatic')
        if thermal:
            print('\t\tNOTE: incompressible run, the T settings of %s are ignored' % (' '.join(thermal)))
    else:
        initArray = ['U uniform (0 0 0);', 'p uniform %.9g;' % (pInit), 'T uniform %.9g;' % (tInit)]
        if isBuoyant(fullCaseSetupDict):
            initArray.append('p_rgh uniform %.9g;' % (pInit))
    turbModel = getTurbulence(fullCaseSetupDict)
    wallFunctionNames = {'low': 'lowReynolds', 'high': 'highReynolds'}
    if turbModel != 'laminar':
        simControl = fullCaseSetupDict['GLOBAL_SIM_CONTROL']
        initArray += ['k uniform %s;' % (simControl['INIT_K'][0]), 'omega uniform %s;' % (simControl['INIT_OMEGA'][0]),
                      'epsilon uniform 1e-4;', 'nut uniform 0;', 'nuTilda uniform 0;', 'alphat uniform 0;']
    bcArray = []
    for geom in geomDict:
        geomName = stripExt(geom)
        if geomName not in wallGroups:
            continue
        wallModel = geomDict[geom]['wallmodel'].lower()
        if wallModel not in wallFunctionNames:
            sys.exit('ERROR! %s wall model column must be low or high!' % (geom))
        for group in wallGroups[geomName]:
            if group['flow']:
                bcArray.append(flowEntry(geom, geomDict, group, fullCaseSetupDict, turbModel))
                continue
            moving = group.get('moving', [])
            static = [pid for pid in group['pids'] if pid not in moving]
            for pids, motion, suffix in ((static, 'stationary', ''), (moving, 'movingMesh', '_moving')):
                if len(pids) == 0:
                    continue
                patches = ' '.join('"%s_%s"' % (geomName, pid) for pid in pids)
                mode = group['T'][0]
                heat = {'adiabatic': 'adiabatic', 'fixed': 'fixedTemperature',
                        'convective': 'convective', 'flux': 'heatFlux'}[mode]
                values = '$:initialConditions;'
                if mode == 'fixed':
                    values += ' T uniform %.9g;' % (group['T'][1])
                elif mode == 'convective':
                    values += ' h %.9g; Ta %.9g;' % (group['T'][1], group['T'][2])
                elif mode == 'flux':
                    values += ' q %.9g;' % (group['T'][1])
                uType = 'noSlip' if motion == 'movingMesh' else group['U']
                if uType == 'slip':
                    options = '' if incompressible else 'heatTransfer %s;' % (heat)
                else:
                    options = 'wallFunction %s; motion %s;' % (wallFunctionNames[wallModel], motion)
                    if not incompressible:
                        options += ' heatTransfer %s;' % (heat)
                if incompressible:
                    values = '$:initialConditions;'
                bcArray.append('    %s%s { category wall; type %s; patches (%s); options { %s } values { %s } }'
                               % (group['name'], suffix, uType, patches, options, values))
    with open('%s/defaultDicts/system/caseProperties' % (templateLoc)) as f:
        text = f.read()
    text = text.replace('<INITIAL_CONDITIONS>', '\n'.join('    ' + s for s in initArray))
    text = text.replace('<BOUNDARY_CONDITIONS>', '\n'.join(bcArray))
    with open('system/caseProperties', 'w') as f:
        f.write(text)


def writeThermoInternal(templateLoc, fullCaseSetupDict):
    #thermophysicalProperties (compressible) or transportProperties (incompressible)
    props = getGasProperties(fullCaseSetupDict)
    if getFlow(fullCaseSetupDict) == 'incompressible':
        print('\n\tWriting transportProperties...')
        nu = getNu(fullCaseSetupDict)
        print('\t\tgas %s, nu %.4g m2/s, reference density %.4g kg/m3 (forces only, p is kinematic)'
              % (getVal(fullCaseSetupDict, 'GLOBAL_MATERIAL', 'GAS'), nu, getGasDensity(fullCaseSetupDict)))
        copyTemplate(templateLoc, 'constant/transportProperties', 'constant/transportProperties', {'NU': '%.6g' % nu})
        return
    print('\n\tWriting thermophysicalProperties...')
    if props['TRANSPORT'] == 'const':
        coeffs = 'mu %.6g; Pr %.6g;' % (props['MU'], props['PR'])
    else:
        coeffs = 'As %.6g; Ts %.6g;' % (props['AS'], props['TS'])
    print('\t\tgas %s, molWeight %g, Cp %g, transport %s (%s)'
          % (getVal(fullCaseSetupDict, 'GLOBAL_MATERIAL', 'GAS'), props['MOL_WEIGHT'], props['CP'],
             props['TRANSPORT'], coeffs))
    replaceDict = {'TRANSPORT': props['TRANSPORT'], 'MOL_WEIGHT': '%.6g' % props['MOL_WEIGHT'],
                   'CP': '%.6g' % props['CP'], 'TRANSPORT_COEFFS': coeffs}
    copyTemplate(templateLoc, 'constant/thermophysicalProperties', 'constant/thermophysicalProperties', replaceDict)
    #buoyant solvers need a rho-based thermo and the gravity vector
    gPath = os.path.join('constant', 'g')
    if isBuoyant(fullCaseSetupDict):
        with open('constant/thermophysicalProperties') as f:
            text = f.read()
        with open('constant/thermophysicalProperties', 'w') as f:
            f.write(text.replace('hePsiThermo', 'heRhoThermo'))
        gravity = [v for v in fullCaseSetupDict['GLOBAL_SIM_CONTROL'].get('GRAVITY', ['0', '0', '-9.81']) if v != '']
        if len(gravity) != 3:
            sys.exit('ERROR! GLOBAL_SIM_CONTROL GRAVITY must be three numbers gx gy gz!')
        print('\t\tbuoyant: heRhoThermo, gravity (%s) m/s2' % (' '.join(gravity)))
        with open(gPath, 'w') as f:
            f.write('FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; object g; }\n\n'
                    'dimensions      [0 1 -2 0 0 0 0];\nvalue           (%s);\n' % (' '.join(gravity)))
    elif os.path.isfile(gPath):
        os.remove(gPath)


def writeFanBaffles(geomDict, fullCaseSetupDict):
    #system/createBafflesDict: each FAN-<name> faceZone becomes a cyclic pair <name>_master / <name>_slave with the
    #fan pressure jump FAN_PRESSURE_RISE (Pa; kinematic for incompressible runs) on p (p_rgh for buoyant runs)
    #returns the fan names (none: no dict, the run script skips createBaffles)
    fans = [stripExt(geom) for geom in geomDict if stripExt(geom).startswith('FAN')]
    path = os.path.join('system', 'createBafflesDict')
    if not fans:
        if os.path.isfile(path):
            os.remove(path)
        return []
    print('\n\tWriting createBafflesDict (fans)...')
    incompressible = getFlow(fullCaseSetupDict) == 'incompressible'
    pName = 'p_rgh' if isBuoyant(fullCaseSetupDict) else 'p'
    pInit = 0.0 if incompressible else float(getVal(fullCaseSetupDict, 'GLOBAL_MATERIAL', 'CHARGE_PRESSURE'))
    rho = getGasDensity(fullCaseSetupDict)
    baffles = []
    for fan in fans:
        try:
            rise = float(fullCaseSetupDict[fan]['FAN_PRESSURE_RISE'][0])
        except (KeyError, ValueError):
            sys.exit('ERROR! [%s] FAN_PRESSURE_RISE must be a number (Pa)!' % (fan))
        jump = rise/rho if incompressible else rise
        print('\t\t%s: pressure rise %g Pa%s' % (fan, rise, ' (%g m2/s2 kinematic)' % (jump) if incompressible else ''))
        fields = ('%s { type fan; patchType cyclic; jump uniform 0; value uniform %.9g; jumpTable constant %.9g; }'
                  % (pName, pInit, jump))
        lines = ['    %s' % (fan), '    {', '        type faceZone;', '        zoneName %s;' % (fan), '        patches',
                 '        {',
                 '            master { name %s_master; type cyclic; neighbourPatch %s_slave; patchFields { %s } }'
                 % (fan, fan, fields),
                 '            slave { name %s_slave; type cyclic; neighbourPatch %s_master; patchFields { %s } }'
                 % (fan, fan, fields),
                 '        }', '    }']
        baffles.append('\n'.join(lines))
    with open(path, 'w') as f:
        f.write('FoamFile { version 2.0; format ascii; class dictionary; object createBafflesDict; }\n\n')
        f.write('//fan baffles from the FAN faceZones (written by internalSetup)\ninternalFacesOnly true;\n\n')
        f.write('baffles\n{\n%s\n}\n' % ('\n'.join(baffles)))
    return fans


#----- porous media (POR-<name>: Darcy-Forchheimer pressure drop, heat exchange with the matrix)

def wireMeshCoeffs(porosity, wireDia, refVelocity, gas):
    #Gedeon-Wood woven screens: f = 129/Re + 2.91 Re^-0.103 on the hydraulic diameter dh = eps*dw/(1-eps) and the
    #pore velocity u = U/eps, Nu = (1 + 0.99 (Re Pr)^0.66) eps^1.79, wetted area per volume a = 4 (1-eps)/dw
    #returns Darcy d (1/m2), Forchheimer f (1/m) for the superficial velocity, and h*a (W/m3/K)
    #the Re^-0.103 term and Nu are taken at refVelocity (superficial)
    dh = porosity*wireDia/(1.0 - porosity)
    re = gas['RHO']*(refVelocity/porosity)*dh/gas['MU']
    d = 129.0/(2.0*porosity*dh**2)
    f = 2.91*re**-0.103/(porosity**2*dh)
    nu = (1.0 + 0.99*(re*gas['PR'])**0.66)*porosity**1.79
    kGas = gas['MU']*gas['CP']/gas['PR']
    ha = nu*kGas/dh*4.0*(1.0 - porosity)/wireDia
    return d, f, ha, re


def fitDarcyForchheimer(velocities, pressureDrops, length, gas):
    #least squares of dp/L = mu*d*U + 0.5*rho*f*U^2 through the origin
    U = np.asarray(velocities, dtype=float)
    dpl = np.asarray(pressureDrops, dtype=float)/length
    A = np.column_stack([gas['MU']*U, 0.5*gas['RHO']*U*np.abs(U)])
    d, f = np.linalg.lstsq(A, dpl, rcond=None)[0]
    return float(d), float(f), A.dot([d, f]) - dpl


def getPorous(geomDict, fullCaseSetupDict):
    #one entry per POR geometry: name, d, f (3 values along e1 e2 e3), e1, e2, heat mode, ha, Cs (J/m3/K), Ts
    porous = []
    zones = [stripExt(geom) for geom in geomDict if stripExt(geom).startswith('POR')]
    if not zones:
        return porous
    print('\n\tSetting up porous zones...')
    props = getGasProperties(fullCaseSetupDict)
    gas = {'RHO': getGasDensity(fullCaseSetupDict), 'MU': props['MU'], 'CP': props['CP'], 'PR': props['PR']}
    flowType = getFlow(fullCaseSetupDict)
    for name in zones:
        setup = fullCaseSetupDict[name]

        def number(key, count=1):
            values = [v for v in setup.get(key, ['']) if v != '']
            try:
                values = [float(v) for v in values]
            except ValueError:
                values = []
            if len(values) != count:
                sys.exit('ERROR! [%s] %s needs %d number(s)!' % (name, key, count))
            return values if count > 1 else values[0]

        zone = {'name': name, 'e1': number('E1', 3), 'e2': number('E2', 3)}
        resistance = setup['RESISTANCE'][0]
        if resistance == 'coeffs':
            zone['d'], zone['f'] = number('D', 3), number('F', 3)
            print('\t\t%s: d (%s) 1/m2, f (%s) 1/m' % (name, ' '.join('%g' % v for v in zone['d']),
                                                    ' '.join('%g' % v for v in zone['f'])))
        elif resistance == 'wireMesh':
            porosity, wireDia, refVelocity = number('POROSITY'), number('WIRE_DIA'), number('REF_VELOCITY')
            if not 0 < porosity < 1 or wireDia <= 0 or refVelocity <= 0:
                sys.exit('ERROR! [%s] wireMesh needs 0 < POROSITY < 1, WIRE_DIA > 0 and REF_VELOCITY > 0!' % (name))
            d, f, ha, re = wireMeshCoeffs(porosity, wireDia, refVelocity, gas)
            zone['d'], zone['f'], zone['wireHA'] = [d]*3, [f]*3, ha
            print('\t\t%s: wire mesh porosity %g, wire %g m, Re %.3g at %g m/s: d %.4g 1/m2, f %.4g 1/m, h*a %.4g W/m3/K'
                  % (name, porosity, wireDia, re, refVelocity, d, f, ha))
        elif resistance == 'fit':
            data = number('DP_DATA', len([v for v in setup.get('DP_DATA', ['']) if v != '']))
            data = data if isinstance(data, list) else [data]
            if len(data) < 4 or len(data) % 2:
                sys.exit('ERROR! [%s] DP_DATA needs at least two velocity / pressure drop pairs (m/s Pa)!' % (name))
            length = number('DP_LENGTH')
            d, f, misfit = fitDarcyForchheimer(data[0::2], data[1::2], length, gas)
            if d < 0 or f < 0:
                sys.exit('ERROR! [%s] DP_DATA fit gives a negative coefficient (d %g, f %g): check the data' % (name, d, f))
            zone['d'], zone['f'] = [d]*3, [f]*3
            print('\t\t%s: fit to %d points over %g m: d %.4g 1/m2, f %.4g 1/m, worst misfit %.3g Pa/m'
                  % (name, len(data)//2, length, d, f, float(np.abs(misfit).max())))
        else:
            sys.exit('ERROR! [%s] RESISTANCE must be coeffs, wireMesh or fit!' % (name))

        #heat exchange with the matrix (compressible only)
        heat = setup.get('HEAT', ['none'])[0]
        if heat not in ('none', 'fixed', 'solved'):
            sys.exit('ERROR! [%s] HEAT must be none, fixed or solved!' % (name))
        if heat != 'none' and flowType == 'incompressible':
            print('\t\tNOTE: %s HEAT %s ignored, incompressible run (no energy equation)' % (name, heat))
            heat = 'none'
        if heat == 'solved' and isSteady(fullCaseSetupDict):
            sys.exit('ERROR! [%s] HEAT solved needs a transient solver (rhoPimpleFoam): the matrix stores heat over time'
                     % (name))
        zone['heat'] = heat
        if heat != 'none':
            haSetting = setup.get('HA', ['auto'])[0]
            if haSetting == 'auto':
                if 'wireHA' not in zone:
                    sys.exit('ERROR! [%s] HA auto needs RESISTANCE wireMesh, else give HA (W/m3/K)' % (name))
                zone['ha'] = zone['wireHA']
            else:
                zone['ha'] = number('HA')
            zone['Ts'] = number('SOLID_T')
            zone['Cs'] = 0.0
            if heat == 'solved':
                zone['Cs'] = (1.0 - number('POROSITY'))*number('SOLID_RHO')*number('SOLID_CP')
            print('\t\t%s: heat %s, h*a %.4g W/m3/K, matrix T %g K%s'
                  % (name, heat, zone['ha'], zone['Ts'],
                     ', (1-eps)*rhos*cs %.4g J/m3/K' % (zone['Cs']) if heat == 'solved' else ' (fixed)'))
        porous.append(zone)
    return porous


def porousHeatSource(zone, cv):
    #scalarCodedSource on e: gas-matrix heat exchange h*a*(Ts - T), linearised in e (implicit part h*a/cv)
    #Ts_<zone> (written) and Qs_<zone> (W/m3 into the gas) are registered fields, Ts is updated once per time step
    #with the exact solution of Cs*dTs/dt = h*a*(T - Ts) (Cs 0: fixed matrix temperature)
    name = zone['name'].replace('-', '_')
    code = '''        codeAddSupRho
        #{
            const fvMesh& m = mesh_;
            const volScalarField& T = m.lookupObject<volScalarField>("T");
            const scalar ha = %.9g;
            const scalar Cv = %.9g;
            const scalar Cs = %.9g;
            const scalar Tinit = %.9g;
            const word TsName("Ts_%s");
            const word QsName("Qs_%s");
            if (!m.foundObject<volScalarField>(TsName))
            {
                volScalarField* TsPtr = new volScalarField
                (
                    IOobject(TsName, m.time().timeName(), m, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                    m,
                    dimensionedScalar(dimTemperature, Tinit)
                );
                TsPtr->store();
                volScalarField* QsPtr = new volScalarField
                (
                    IOobject(QsName, m.time().timeName(), m, IOobject::NO_READ, IOobject::NO_WRITE),
                    m,
                    dimensionedScalar(dimPower/dimVolume, Zero)
                );
                QsPtr->store();
            }
            volScalarField& Ts = const_cast<volScalarField&>(m.lookupObject<volScalarField>(TsName));
            volScalarField& Qs = const_cast<volScalarField&>(m.lookupObject<volScalarField>(QsName));
            static label lastIndex = -1;
            const bool update = (Cs > 0 && m.time().timeIndex() != lastIndex);
            lastIndex = m.time().timeIndex();
            const scalar dt = m.time().deltaTValue();
            const labelList& zoneCells = this->cells();
            const scalarField& V = m.V();
            const scalarField& e = eqn.psi();
            scalarField& diag = eqn.diag();
            scalarField& source = eqn.source();
            const scalar k = ha/Cv;
            forAll(zoneCells, i)
            {
                const label c = zoneCells[i];
                if (update)
                {
                    Ts[c] = T[c] + (Ts[c] - T[c])*Foam::exp(-ha*dt/Cs);
                }
                Qs[c] = ha*(Ts[c] - T[c]);
                diag[c] -= k*V[c];
                source[c] -= (Qs[c] + k*e[c])*V[c];
            }
        #};''' % (zone['ha'], cv, zone['Cs'], zone['Ts'], zone['name'], zone['name'])
    lines = ['%sHeat' % (name), '{', '    type            scalarCodedSource;', '    name            %sHeat;' % (name),
             '    scalarCodedSourceCoeffs', '    {', '        selectionMode   cellZone;',
             '        cellZone        %s;' % (zone['name']), '        fields          (e);', '',
             '        codeInclude', '        #{', '            #include "fvCFD.H"', '        #};', '',
             '        codeCorrect #{ #};', '        codeAddSup #{ #};', '        codeConstrain #{ #};', '', code,
             '    }', '}']
    return '\n'.join(lines)


def writeFvOptionsInternal(porous, fullCaseSetupDict):
    #constant/fvOptions: explicitPorositySource per POR zone, plus the heat source when HEAT is fixed or solved
    path = os.path.join('constant', 'fvOptions')
    if not porous:
        if os.path.isfile(path):
            os.remove(path)
        return
    print('\n\tWriting fvOptions (porous zones)...')
    props = getGasProperties(fullCaseSetupDict)
    cv = props['CP'] - 8314.46/props['MOL_WEIGHT']
    entries = []
    for zone in porous:
        name = zone['name'].replace('-', '_')
        entries.append('\n'.join([
            '%sPorosity' % (name), '{', '    type            explicitPorositySource;',
            '    explicitPorositySourceCoeffs', '    {', '        selectionMode   cellZone;',
            '        cellZone        %s;' % (zone['name']), '        type            DarcyForchheimer;',
            '        DarcyForchheimerCoeffs', '        {',
            '            d   (%s);' % (' '.join('%.9g' % v for v in zone['d'])),
            '            f   (%s);' % (' '.join('%.9g' % v for v in zone['f'])),
            '            coordinateSystem',
            '            {', '                origin  (0 0 0);',
            '                rotation { type axes; e1 (%s); e2 (%s); }' % (' '.join('%g' % v for v in zone['e1']),
                                                                          ' '.join('%g' % v for v in zone['e2'])),
            '            }', '        }', '    }', '}']))
        if zone['heat'] != 'none':
            entries.append(porousHeatSource(zone, cv))
    with open(path, 'w') as f:
        f.write('FoamFile { version 2.0; format ascii; class dictionary; object fvOptions; }\n\n')
        f.write('//porous zones (written by internalSetup)\n\n%s\n' % ('\n\n'.join(entries)))
    if any(zone['heat'] != 'none' for zone in porous):
        print('\t\tthe heat source is compiled at run time (needs a C++ compiler where the solver runs)')


def writeTurbulenceInternal(templateLoc, fullCaseSetupDict):
    turbModel = getTurbulence(fullCaseSetupDict)
    print('\n\tWriting turbulenceProperties (%s)...' % (turbModel))
    if turbModel == 'laminar':
        replaceDict = {'SIMULATION_TYPE': 'laminar', 'RAS_BLOCK': ''}
    else:
        replaceDict = {'SIMULATION_TYPE': 'RAS',
                       'RAS_BLOCK': 'RAS { RASModel %s; turbulence on; printCoeffs on; }' % (turbModel)}
    copyTemplate(templateLoc, 'constant/turbulenceProperties', 'constant/turbulenceProperties', replaceDict)


def getRunTimes(fullCaseSetupDict):
    #endTime and writeInterval from FREQ, N_CYCLES, WRITE_PER_CYCLE
    control = fullCaseSetupDict['GLOBAL_CONTROL']
    freq = float(control['FREQ'][0])
    nCycles = float(control['N_CYCLES'][0])
    writePerCycle = float(control['WRITE_PER_CYCLE'][0])
    if freq <= 0 or nCycles <= 0 or writePerCycle <= 0:
        sys.exit('ERROR! GLOBAL_CONTROL FREQ, N_CYCLES and WRITE_PER_CYCLE must be positive!')
    return nCycles/freq, 1.0/(freq*writePerCycle)


def getTimeStep(fullCaseSetupDict):
    #MAX_CO 0 or blank uses a fixed time step of MAX_DT, otherwise adjustable
    control = fullCaseSetupDict['GLOBAL_CONTROL']
    maxCo = ''.join(control.get('MAX_CO', [''])).strip()
    maxDt = ''.join(control.get('MAX_DT', [''])).strip()
    try:
        if float(maxDt) <= 0:
            raise ValueError
    except ValueError:
        sys.exit('ERROR! GLOBAL_CONTROL MAX_DT must be a positive number (got "%s")!' % (maxDt))
    try:
        fixed = maxCo == '' or float(maxCo) == 0
    except ValueError:
        sys.exit('ERROR! GLOBAL_CONTROL MAX_CO must be a number, 0 or blank (got "%s")!' % (maxCo))
    if fixed:
        return maxDt, '1e6', maxDt, 'no'
    if float(maxCo) < 0:
        sys.exit('ERROR! GLOBAL_CONTROL MAX_CO must not be negative!')
    return control['INIT_DT'][0], maxCo, maxDt, 'yes'


def getSteadyIterations(fullCaseSetupDict):
    #steady solvers: N_ITER iterations, written every WRITE_ITER
    control = fullCaseSetupDict['GLOBAL_CONTROL']
    nIter = int(control.get('N_ITER', ['2000'])[0])
    writeIter = int(control.get('WRITE_ITER', ['500'])[0])
    if nIter <= 0 or writeIter <= 0:
        sys.exit('ERROR! GLOBAL_CONTROL N_ITER and WRITE_ITER must be positive!')
    return nIter, writeIter


def writeControlDictInternal(templateLoc, fullCaseSetupDict):
    print('\n\tWriting controlDict...')
    control = fullCaseSetupDict['GLOBAL_CONTROL']
    if isSteady(fullCaseSetupDict):
        nIter, writeIter = getSteadyIterations(fullCaseSetupDict)
        print('\t\tsteady: %d iterations, write every %d' % (nIter, writeIter))
        replaceDict = {'SOLVER': getVal(fullCaseSetupDict, 'GLOBAL_SIM_CONTROL', 'SOLVER'),
                       'STARTFROM': control['STARTFROM'][0], 'ENDTIME': nIter, 'INIT_DT': 1,
                       'WRITE_INTERVAL': writeIter, 'MAX_CO': '1e6', 'MAX_DT': 1, 'ADJUST_DT': 'no'}
        copyTemplate(templateLoc, 'system/controlDict', 'system/controlDict', replaceDict)
        return
    endTime, writeInterval = getRunTimes(fullCaseSetupDict)
    print('\t\t%s cycles at %s Hz: endTime %.6g s, write every %.6g s'
          % (control['N_CYCLES'][0], control['FREQ'][0], endTime, writeInterval))
    initDt, maxCo, maxDt, adjustDt = getTimeStep(fullCaseSetupDict)
    if adjustDt == 'yes':
        print('\t\tadjustable time step: maxCo %s, maxDeltaT %s' % (maxCo, maxDt))
    else:
        print('\t\tfixed time step %s (MAX_CO 0 or blank), %.6g steps per cycle'
              % (maxDt, 1.0/(float(control['FREQ'][0])*float(maxDt))))
        #without time step adjustment the writes only land on the interval if it is a whole number of steps
        steps = writeInterval/float(maxDt)
        if abs(steps - round(steps)) > 1e-6*steps:
            print('\t\tWARNING! write interval %.6g s is not a whole number of time steps (%.4g): write times '
                  'will be the first step past each interval' % (writeInterval, steps))
    replaceDict = {'SOLVER': getVal(fullCaseSetupDict, 'GLOBAL_SIM_CONTROL', 'SOLVER'),
                   'STARTFROM': control['STARTFROM'][0], 'ENDTIME': '%.9g' % endTime,
                   'INIT_DT': initDt, 'WRITE_INTERVAL': '%.9g' % writeInterval,
                   'MAX_CO': maxCo, 'MAX_DT': maxDt, 'ADJUST_DT': adjustDt}
    copyTemplate(templateLoc, 'system/controlDict', 'system/controlDict', replaceDict)


def getRelaxation(fullCaseSetupDict):
    #under-relaxation: always for the steady solvers, transient only with outer correctors (N_OUTER_CORR > 1)
    #the last outer iteration (Final) is not relaxed
    simControl = fullCaseSetupDict['GLOBAL_SIM_CONTROL']
    nOuter = int(simControl['N_OUTER_CORR'][0])
    if nOuter < 1:
        sys.exit('ERROR! GLOBAL_SIM_CONTROL N_OUTER_CORR must be 1 or more!')
    steady = isSteady(fullCaseSetupDict)
    if nOuter == 1 and not steady:
        return ''
    relax = {}
    for key, default in (('RELAX_P', '0.3'), ('RELAX_U', '0.7'), ('RELAX_E', '0.7')):
        try:
            relax[key] = float(simControl.get(key, [default])[0])
        except ValueError:
            relax[key] = -1.0
        if not 0 < relax[key] <= 1:
            sys.exit('ERROR! GLOBAL_SIM_CONTROL %s must be > 0 and <= 1!' % (key))
    compressible = getFlow(fullCaseSetupDict) == 'compressible'
    pName = 'p_rgh' if isBuoyant(fullCaseSetupDict) else 'p'
    equations = ['        "(U|k|omega|epsilon)" %g;' % (relax['RELAX_U'])]
    if compressible:
        equations.append('        e               %g;' % (relax['RELAX_E']))
    fields = ['        %-15s %g;' % (pName, relax['RELAX_P']), '        %-15s 1;' % ('%sFinal' % (pName))]
    if compressible:
        fields.append('        rho             1;')
    if steady:
        print('\t\tsteady under-relaxation %s %g, U %g%s' % (pName, relax['RELAX_P'], relax['RELAX_U'],
                                                            ', e %g' % (relax['RELAX_E']) if compressible else ''))
    else:
        print('\t\t%d outer correctors, under-relaxation p %g, U %g%s (none on the last outer iteration)'
              % (nOuter, relax['RELAX_P'], relax['RELAX_U'], ', e %g' % (relax['RELAX_E']) if compressible else ''))
    return ('relaxationFactors\n{\n    fields\n    {\n%s\n    }\n    equations\n    {\n%s\n        ".*Final"       1;\n    }\n}\n'
            % ('\n'.join(fields), '\n'.join(equations)))


def getAlgorithm(fullCaseSetupDict, pRef, pressureBC):
    #SIMPLE (steady) or PIMPLE (transient) block of fvSolution
    simControl = fullCaseSetupDict['GLOBAL_SIM_CONTROL']
    if isSteady(fullCaseSetupDict):
        lines = ['SIMPLE', '{', '    nNonOrthogonalCorrectors %s;' % (simControl['N_NONORTH_CORR'][0])]
        if getFlow(fullCaseSetupDict) == 'compressible' and pressureBC:
            #limits p while the steady run settles (rhoSimpleFoam pressureControl)
            lines += ['    pMinFactor      0.1;', '    pMaxFactor      2;']
        if pRef:
            lines.append(pRef)
        lines += ['    residualControl', '    {', '        "(p|p_rgh)"     1e-4;', '        U               1e-5;',
                  '        e               1e-5;', '        "(k|omega|epsilon)" 1e-5;', '    }', '}']
        return '\n'.join(lines)
    lines = ['PIMPLE', '{', '    momentumPredictor yes;', '    nOuterCorrectors %s;' % (simControl['N_OUTER_CORR'][0]),
             '    nCorrectors %s;' % (simControl['N_CORR'][0]),
             '    nNonOrthogonalCorrectors %s;' % (simControl['N_NONORTH_CORR'][0]), '    correctPhi no;']
    if pRef:
        lines.append(pRef)
    return '\n'.join(lines + ['}'])


def writeNumericsInternal(templateLoc, fullCaseSetupDict, pressureBC=False):
    #pressureBC: an inlet or outlet fixes the pressure level, so no pressure reference is needed
    print('\n\tWriting fvSchemes, fvSolution, decomposeParDict...')
    compute = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP']
    steady = isSteady(fullCaseSetupDict)
    #divergence schemes per flow type, closed incompressible volume needs a pressure reference
    #linearUpwind keeps it stable at large Co (010 diverged at Co 2.4 with limitedLinearV)
    #steady: bounded convection terms (continuity is only met at convergence)
    bounded = 'bounded ' if steady else ''
    divTurb = ['div(phi,k)      %sGauss linearUpwind grad(k);' % (bounded),
               'div(phi,omega)  %sGauss linearUpwind grad(omega);' % (bounded),
               'div(phi,epsilon) %sGauss linearUpwind grad(epsilon);' % (bounded)]
    pRef = ''
    if getFlow(fullCaseSetupDict) == 'incompressible':
        divSchemes = ['div(phi,U)      %sGauss linearUpwind grad(U);' % (bounded)] + divTurb \
                     + ['div((nuEff*dev2(T(grad(U))))) Gauss linear;']
        if not pressureBC:
            locInMesh = ' '.join(v for v in fullCaseSetupDict['GLOBAL_REFINEMENT']['LOC_IN_MESH'] if v != '')
            pRef = '    pRefPoint ( %s );\n    pRefValue 0;' % (locInMesh)
    else:
        divSchemes = ['div(phi,U)      %sGauss linearUpwind grad(U);' % (bounded),
                      'div(phi,e)      %sGauss linearUpwind grad(e);' % (bounded),
                      'div(phi,K)      %sGauss linear;' % (bounded),
                      'div(phi,Ekp)    %sGauss linear;           // kinetic energy + p/rho in the steady e equation' % (bounded),
                      'div(phiv,p)     Gauss linear;           // pressure work in the e equation',
                      'div(meshPhi,p)  Gauss linear;           // pressure work from the moving mesh'] \
                     + divTurb + ['div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;']
    copyTemplate(templateLoc, 'system/fvSchemes', 'system/fvSchemes',
                 {'DIV_SCHEMES': '\n'.join('    ' + s for s in divSchemes),
                  'DDT_SCHEME': 'steadyState' if steady else 'Euler'})
    relaxation = getRelaxation(fullCaseSetupDict)
    #older custom template sets have no <RELAXATION> / <ALGORITHM> slots
    with open('%s/defaultDicts/system/fvSolution' % (templateLoc)) as f:
        template = f.read()
    for slot, needed in (('<RELAXATION>', bool(relaxation)), ('<ALGORITHM>', True)):
        if needed and slot not in template:
            print('\t\tWARNING: %s/defaultDicts/system/fvSolution has no %s line: update it from setupTemplates/internal'
                  % (templateLoc, slot))
    copyTemplate(templateLoc, 'system/fvSolution', 'system/fvSolution',
                 {'ALGORITHM': getAlgorithm(fullCaseSetupDict, pRef, pressureBC), 'RELAXATION': relaxation,
                  'P_RELTOL': '0.05' if steady else '0.01', 'U_RELTOL': '0.1' if steady else '0'})
    method, decomp = getDecomp(fullCaseSetupDict)
    nProcs, coeffs = compute['NUM_PROCS'][0], ''
    if decomp:
        nProcs = str(decomp[0]*decomp[1]*decomp[2])
        coeffs = 'coeffs\n{\n    n               (%d %d %d);\n}' % decomp
        print('\t\tNumber of processors: %s' % (nProcs))
        print('\t\tDecomposition: %d %d %d' % decomp)
        if nProcs != compute['LOCAL_PROCS'][0]:
            print('\t\tNOTE: local runs need LOCAL_PROCS = %s for this DECOMP' % (nProcs))
    elif method in ('hierarchical', 'simple'):
        print('\t\t%s decomposition: n set by the run scripts for the task count of each run' % (method))
    copyTemplate(templateLoc, 'system/decomposeParDict', 'system/decomposeParDict',
                 {'NUM_PROCS': nProcs, 'DECOMP_METHOD': method, 'DECOMP_COEFFS': coeffs})


def getDecomp(fullCaseSetupDict):
    #getting the decomposition method and the hierarchical/simple split (DECOMP = nx ny nz or auto)
    compute = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP']
    method = compute['DECOMP_METHOD'][0]
    decomp = [v for v in compute.get('DECOMP', ['auto']) if v != '']
    if method not in ('hierarchical', 'simple') or decomp in ([], ['auto']):
        return method, None
    try:
        n = tuple(int(v) for v in decomp)
    except ValueError:
        n = ()
    if len(n) != 3 or min(n) < 1:
        sys.exit('ERROR! GLOBAL_COMPUTE_SETUP DECOMP must be auto or three whole numbers nx ny nz, got %s' % (' '.join(decomp)))
    total = n[0]*n[1]*n[2]
    if total != int(compute['NUM_PROCS'][0]):
        sys.exit('ERROR! DECOMP %d %d %d gives %d subdomains but NUM_PROCS is %s: set NUM_PROCS = %d (job tasks)'
                 % (n + (total, compute['NUM_PROCS'][0], total)))
    return method, n


def functionObject(name, entries, fieldWriter=False):
    #one function object block, written every time step
    #fieldWriter: execute every step but only write the field at write times (else a time dir every step)
    lines = ['%s' % (name), '{'] + ['    %-15s %s' % (key, val) for key, val in entries]
    if fieldWriter:
        lines += ['    executeControl  timeStep;', '    executeInterval 1;', '    writeControl    writeTime;', '}']
    else:
        lines += ['    writeControl    timeStep;', '    writeInterval   1;', '}']
    return '\n'.join(lines)


def writeFunctionObjectsInternal(wallGroups, fullCaseSetupDict, motionInfo=None, porous=None):
    #gas averages, mass (must stay constant in a closed volume) and wall heat flux for compressible runs,
    #pressure force per moving body and area-averaged p on each motion face (both flow types)
    #gas volume for post: V = volIntegrate(rho)/volAverage(rho), no extra function object needed
    print('\n\tWriting internalFunctionObjects...')
    incompressible = getFlow(fullCaseSetupDict) == 'incompressible'
    allPatches, flowPatches = [], []
    for geomName in wallGroups:
        for group in wallGroups[geomName]:
            patches = ['%s_%s' % (geomName, pid) for pid in group['pids']]
            if group.get('flow'):
                flowPatches += patches
            else:
                allPatches += ['"%s"' % (patch) for patch in patches]
    volume = [('type', 'volFieldValue;'), ('libs', '(fieldFunctionObjects);'), ('regionType', 'all;')]
    fos = [functionObject('gasAvg', volume + [('operation', 'volAverage;'),
                                              ('fields', '(p);' if incompressible else '(p T rho);'),
                                              ('writeFields', 'false;')])]
    #inlets and outlets: flux (kg/s compressible, m3/s incompressible, + out of the domain) and area-averaged p, T
    for patch in flowPatches:
        surface = [('type', 'surfaceFieldValue;'), ('libs', '(fieldFunctionObjects);'), ('regionType', 'patch;'),
                   ('name', '%s;' % (patch)), ('writeFields', 'false;')]
        fos.append(functionObject('flow_%s' % (patch), surface + [('operation', 'sum;'), ('fields', '(phi);')]))
        fos.append(functionObject('avg_%s' % (patch), surface + [('operation', 'areaAverage;'),
                                                                 ('fields', '(p);' if incompressible else '(p T);')]))
    #porous zones: zone averages (incl. the matrix temperature) and the heat from the matrix into the gas (W)
    for zone in (porous or []):
        cellZone = [('type', 'volFieldValue;'), ('libs', '(fieldFunctionObjects);'), ('regionType', 'cellZone;'),
                    ('name', '%s;' % (zone['name'])), ('writeFields', 'false;')]
        fields = ['p'] if incompressible else ['p', 'T']
        if zone.get('heat', 'none') != 'none':
            fields.append('Ts_%s' % (zone['name']))
            fos.append(functionObject('porHeat_%s' % (zone['name']), cellZone + [('operation', 'volIntegrate;'),
                                                                               ('fields', '(Qs_%s);' % (zone['name']))]))
        fos.append(functionObject('por_%s' % (zone['name']), cellZone + [('operation', 'volAverage;'),
                                                                       ('fields', '(%s);' % (' '.join(fields)))]))
    #steady: initial residuals per iteration for the convergence plot
    if isSteady(fullCaseSetupDict):
        fields = ['U', 'p_rgh' if isBuoyant(fullCaseSetupDict) else 'p'] + ([] if incompressible else ['e'])
        if getTurbulence(fullCaseSetupDict) != 'laminar':
            fields += ['k', 'omega']
        fos.append(functionObject('residuals', [('type', 'solverInfo;'), ('libs', '(utilityFunctionObjects);'),
                                                ('fields', '(%s);' % (' '.join(fields)))]))
    if not incompressible:
        fos.append(functionObject('gasMass', volume + [('operation', 'volIntegrate;'), ('fields', '(rho);'),
                                                       ('writeFields', 'false;')]))
        fos.append(functionObject('wallHeatFlux', [('type', 'wallHeatFlux;'), ('libs', '(fieldFunctionObjects);'),
                                                   ('patches', '(%s);' % (' '.join(allPatches)))],
                                   fieldWriter=True))
    if motionInfo:
        #kinematic p for incompressible runs: forces scaled by the gas density at the charge state
        if incompressible:
            rho = [('rho', 'rhoInf;'), ('rhoInf', '%.6g;' % (getGasDensity(fullCaseSetupDict)))]
        else:
            rho = [('rho', 'rho;')]
        for body in motionInfo['bodies']:
            patches = ' '.join('"%s_%s"' % (body.source, pid) for pid in body.movedPids)
            fos.append(functionObject('force_%s' % (body.name), [('type', 'forces;'), ('libs', '(forces);'),
                                                                 ('patches', '(%s);' % (patches))] + rho
                                      + [('CofR', '(0 0 0);')]))
        for zone in motionInfo['zones']:
            if zone.body is None:
                continue
            fos.append(functionObject('p_%s' % (zone.patch), [('type', 'surfaceFieldValue;'),
                                                              ('libs', '(fieldFunctionObjects);'),
                                                              ('regionType', 'patch;'), ('name', '%s;' % (zone.patch)),
                                                              ('operation', 'areaAverage;'), ('fields', '(p);'),
                                                              ('writeFields', 'false;')]))
    with open('system/internalFunctionObjects', 'w') as f:
        f.write('\n'.join(fos) + '\n')


#----- mesh motion (MOVE bodies, displacementInterpolation)

def matchPids(keywords, pids, where):
    #PIDs containing any of the keywords (case-insensitive), every keyword must match
    #a keyword that is a whole PID name only takes that PID (MOVE-piston does not also take MOVE-piston-skirt)
    matched = []
    for kw in keywords:
        hits = [pid for pid in pids if kw.lower() == pid.lower()] or [pid for pid in pids if kw.lower() in pid.lower()]
        if len(hits) == 0:
            sys.exit('ERROR! %s: keyword "%s" matches no PID. PIDs are: %s' % (where, kw, pids))
        matched += [pid for pid in hits if pid not in matched]
    return matched


def getMotionFilePath(setup, section, caseDir):
    #MOTION_FILE is a path (no default): absolute, ~ and $VARS expanded, relative paths from the case directory
    value = ' '.join(v for v in setup['MOTION_FILE'] if v != '')
    if value in ('', 'required'):
        sys.exit('ERROR! [%s] MOTION_FILE is required: path to the motion csv (absolute, or relative to the case directory)'
                 % (section))
    path = os.path.expandvars(os.path.expanduser(value))
    if not os.path.isabs(path):
        path = os.path.join(caseDir, path)
    path = os.path.normpath(path)
    if not os.path.isfile(path):
        sys.exit('ERROR! [%s] MOTION_FILE %s not found (looked for %s)' % (section, value, path))
    return path


def getCrank(fullCaseSetupDict):
    #crank drive settings, None for prescribed motion
    meshMotion = fullCaseSetupDict['MESH_MOTION']
    drive = meshMotion.get('DRIVE', ['prescribed'])[0]
    if drive == 'prescribed':
        return None
    if drive != 'crank':
        sys.exit('ERROR! MESH_MOTION DRIVE must be prescribed or crank!')
    crank = {'LOAD': meshMotion['LOAD'][0]}
    if crank['LOAD'] not in ('speed', 'torque', 'damper', 'generator'):
        sys.exit('ERROR! MESH_MOTION LOAD must be speed, torque, damper or generator!')
    for key in ['INERTIA', 'LOAD_TORQUE', 'LOAD_DAMPING', 'GEN_K', 'GEN_R_LOAD', 'GEN_R_INT',
                'FRICTION_TORQUE', 'FRICTION_COEFF']:
        crank[key] = float(meshMotion[key][0])
    if crank['LOAD'] != 'speed' and crank['INERTIA'] <= 0:
        sys.exit('ERROR! MESH_MOTION INERTIA must be > 0 (flywheel, crank and the parts turning with it, kg m2)')
    if crank['LOAD'] == 'generator' and (crank['GEN_K'] <= 0 or crank['GEN_R_LOAD'] + crank['GEN_R_INT'] <= 0):
        sys.exit('ERROR! MESH_MOTION generator needs GEN_K > 0 and GEN_R_LOAD + GEN_R_INT > 0')
    crank['OMEGA0'] = 2.0*math.pi*float(getVal(fullCaseSetupDict, 'GLOBAL_CONTROL', 'FREQ'))

    #pressure behind the bodies (bounce space is not modelled)
    if getFlow(fullCaseSetupDict) == 'compressible':
        pBuffer = meshMotion['P_BUFFER'][0]
        crank['P_REF'] = float(getVal(fullCaseSetupDict, 'GLOBAL_MATERIAL', 'CHARGE_PRESSURE')) if pBuffer == 'charge' else float(pBuffer)
        crank['RHO_INF'] = 1.0
    else:
        crank['P_REF'] = 0.0
        crank['RHO_INF'] = getGasDensity(fullCaseSetupDict)
        if crank['LOAD'] != 'speed':
            print('\t\tWARNING: DRIVE crank with an incompressible solver: the gas is not compressed, so the piston '
                  'torque is not physical (use LOAD speed, or rhoPimpleFoam)')
    return crank


def crankEntries(crank, startAngle):
    #crank sub-dictionary of dynamicMeshDict
    lines = ['theta0          %.9g;   //crank angle at the start (meshing) angle, deg' % (startAngle),
             'omega0          %.9g;   //start speed, rad/s (FREQ)' % (crank['OMEGA0']),
             'inertia         %.9g;   //kg m2' % (crank['INERTIA']),
             'load            %s;' % (crank['LOAD'])]
    if crank['LOAD'] == 'torque':
        lines.append('torque          %.9g;   //N m' % (crank['LOAD_TORQUE']))
    elif crank['LOAD'] == 'damper':
        lines.append('damping         %.9g;   //N m s' % (crank['LOAD_DAMPING']))
    elif crank['LOAD'] == 'generator':
        lines.append('generator { k %.9g; RLoad %.9g; RInternal %.9g; }   //V s/rad, ohm, ohm'
                     % (crank['GEN_K'], crank['GEN_R_LOAD'], crank['GEN_R_INT']))
    lines += ['frictionTorque  %.9g;   //N m, against the rotation' % (crank['FRICTION_TORQUE']),
              'frictionCoeff   %.9g;   //N m s' % (crank['FRICTION_COEFF']),
              'pRef            %.9g;   //pressure behind the bodies (units of p)' % (crank['P_REF']),
              'rhoInf          %.9g;   //density for kinematic p' % (crank['RHO_INF'])]
    return '\n'.join('        %s' % (line) for line in lines)


def getMoveSections(fullCaseSetupDict):
    return [section for section in fullCaseSetupDict if section.startswith('MOVE')]


def buildMotion(geomDict, fullCaseSetupDict, pidDict, mshDir, caseDir):
    #zones, bodies, start angle, geometry shifts and gap refinement for the MOVE bodies
    #returns None for a static case
    moveSections = getMoveSections(fullCaseSetupDict)
    if len(moveSections) == 0:
        return None
    print('\n\tSetting up mesh motion...')
    meshMotion = fullCaseSetupDict['MESH_MOTION']
    if meshMotion['METHOD'][0] != 'displacementInterpolation' or meshMotion['AXIS'][0] != 'z':
        sys.exit('ERROR! MESH_MOTION: only METHOD displacementInterpolation along AXIS z is available')
    step = float(meshMotion['STEP_DEG'][0])
    if step <= 0 or step > 1:
        sys.exit('ERROR! MESH_MOTION STEP_DEG must be > 0 and <= 1 deg')
    maxDeformation = float(meshMotion['MAX_DEFORMATION'][0])
    flatTol = float(meshMotion['ZONE_FLAT_TOL'][0])
    control = fullCaseSetupDict['GLOBAL_CONTROL']
    freq, nCycles = float(control['FREQ'][0]), float(control['N_CYCLES'][0])
    baseCell = float(getVal(fullCaseSetupDict, 'GLOBAL_REFINEMENT', 'BASE_CELL'))
    getCrank(fullCaseSetupDict)

    #triangles per PID of every geometry, scaled to metres
    geomFiles = {stripExt(geom): geom for geom in geomDict}
    tris = {}
    for geomName, geom in geomFiles.items():
        scale = float(geomDict[geom]['scale'])
        tris[geomName] = {pid: t*scale for pid, t in
                          geometryIO.readSurfaceTriangles(os.path.join(mshDir, geom)).items()}

    #anchors: flat fixed faces of the WALL geometry
    zones = []
    for geomName in geomFiles:
        if not geomName.startswith('WALL'):
            continue
        keywords = [k for k in fullCaseSetupDict[geomName].get('ANCHOR_PIDS', ['']) if k != '']
        for pid in matchPids(keywords, pidDict[geomName], '[%s] ANCHOR_PIDS' % (geomName)):
            t = tris[geomName][pid]
            zones.append(motionUtilities.Zone('motion_%s_%s' % (geomName, pid), '%s_%s' % (geomName, pid),
                                              motionUtilities.flatZoneZ('%s_%s' % (geomName, pid), t, flatTol),
                                              motionUtilities.projectedArea(t)))
    if len(zones) == 0:
        sys.exit('ERROR! Mesh motion needs fixed anchor faces: set ANCHOR_PIDS in a WALL section (e.g. hot cold)')

    #bodies and their zones
    bodies, movingPids, profiles = [], {}, {}
    for section in moveSections:
        setup = fullCaseSetupDict[section]
        name = section[len('MOVE-'):]
        source = section if setup['FROM'][0] == 'self' else setup['FROM'][0]
        if source not in geomFiles:
            sys.exit('ERROR! [%s] FROM %s is not a geometry in GEOM' % (section, source))
        motionFile = getMotionFilePath(setup, section, caseDir)
        if motionFile not in profiles:
            print('\t\t%s motion from %s' % (name, motionFile))
            profiles[motionFile] = motionUtilities.readMotionFile(motionFile)
        theta, columns = profiles[motionFile]
        column = setup['MOTION_COLUMN'][0]
        if column not in columns:
            sys.exit('ERROR! [%s] MOTION_COLUMN %s not in %s (columns: %s)'
                     % (section, column, setup['MOTION_FILE'][0], list(columns.keys())))
        thetaC, pos = motionUtilities.periodicProfile(theta, columns[column], column)
        body = motionUtilities.Body(name, thetaC, pos, float(setup['REF_ANGLE'][0]), float(setup['TDC_CLEARANCE'][0]))
        body.source = source
        keywords = [k for k in setup['ZONE_PIDS'] if k != '']
        if len(keywords) == 0:
            sys.exit('ERROR! [%s] ZONE_PIDS is empty: list the flat faces that drive the motion (e.g. top bottom)' % (section))
        zonePids = matchPids(keywords, pidDict[source], '[%s] ZONE_PIDS' % (section))
        for pid in zonePids:
            t = tris[source][pid]
            zones.append(motionUtilities.Zone('motion_%s_%s' % (source, pid), '%s_%s' % (source, pid),
                                              motionUtilities.flatZoneZ('%s_%s' % (source, pid), t, flatTol),
                                              motionUtilities.projectedArea(t), body))
        #own file: the whole body moves, FROM: only the zone faces (e.g. a piston crown)
        body.movedPids = list(pidDict[source]) if source == section else zonePids
        movingPids.setdefault(source, [])
        movingPids[source] += [pid for pid in body.movedPids if pid not in movingPids[source]]
        bodies.append(body)

    anchorPatches = [zone.patch for zone in zones if zone.body is None]
    for body in bodies:
        clash = ['%s_%s' % (body.source, pid) for pid in body.movedPids if '%s_%s' % (body.source, pid) in anchorPatches]
        if clash:
            sys.exit('ERROR! %s are both anchors and part of the moving body %s' % (clash, body.name))

    #clearance, order, start angle, deformation
    offsets = motionUtilities.applyClearance(zones, bodies, step)
    for bodyName, (offset, addedVolume) in offsets.items():
        print('\t\t%s moved %+.4g m to keep its TDC clearance: dead volume +%.4g m3 (%.4g cm3)'
              % (bodyName, offset, addedVolume, addedVolume*1e6))
    motionUtilities.checkZoneOrder(zones, step)
    startSetting = meshMotion['START_ANGLE'][0]
    startAngle = motionUtilities.autoStartAngle(zones, step) if startSetting == 'auto' else float(startSetting)
    ratio, pair = motionUtilities.deformationRatio(zones, startAngle, step)
    print('\t\tstart (meshing) angle %.1f deg%s, worst squash/stretch %.2fx (%s - %s)'
          % (startAngle, ' (auto)' if startSetting == 'auto' else '', ratio, pair[0], pair[1]))
    if ratio > maxDeformation:
        sys.exit('ERROR! Squash/stretch %.2fx exceeds MAX_DEFORMATION %.2f: increase TDC_CLEARANCE or the limit'
                 % (ratio, maxDeformation))
    print('\t\tzones in z order at the start angle:')
    order = motionUtilities.sortedZones(zones, step)
    for zone in order:
        print('\t\t\t%-40s z %.6g m  %s' % (zone.name, float(motionUtilities.zoneZ(zone, startAngle)),
                                            'fixed' if zone.body is None else 'moves with %s' % (zone.body.name)))
    for body in bodies:
        kin = motionUtilities.bodyKinematics(body, freq, step)
        print('\t\t%s: stroke %.4g m, peak velocity %.4g m/s, peak acceleration %.4g m/s2'
              % (body.name, kin['stroke'], kin['vMax'], kin['aMax']))
        if kin['vJump'] > 0.2*kin['vMax']:
            print('\t\tWARNING: %s velocity jumps by %.4g m/s between samples: a sharp change in the profile '
                  'gives pressure spikes' % (body.name, kin['vJump']))

    #static wall faces dragged by the mesh (horizontal faces inside a stretching band)
    zonePatches = [zone.patch for zone in zones]
    staticTris = {}
    for geomName in geomFiles:
        if not geomName.startswith('WALL'):
            continue
        for pid, t in tris[geomName].items():
            if '%s_%s' % (geomName, pid) not in zonePatches and pid not in movingPids.get(geomName, []):
                staticTris['%s_%s' % (geomName, pid)] = t
    zoneZs = [float(motionUtilities.zoneZ(zone, startAngle)) for zone in zones]
    for label, zLo, zHi in motionUtilities.unanchoredHorizontalFaces(staticTris, zoneZs, zones, startAngle, step, flatTol):
        print('\t\tWARNING: %s has horizontal faces at z %.6g..%.6g m inside a stretching band: they will move '
              'with the mesh. Add them to ANCHOR_PIDS if they are fixed steps or port ends' % (label, zLo, zHi))

    #geometry shifts to the start angle (file units)
    shifts = {}
    for body in bodies:
        dz = float(body.shift(startAngle))
        scale = float(geomDict[geomFiles[body.source]]['scale'])
        for pid in body.movedPids:
            shifts.setdefault(body.source, {})[pid] = (0.0, 0.0, dz/scale)
    #porous zones inside a moving body (regenerator in the displacer) are moved with it (MOVES_WITH = <body>)
    bodyNames = dict((body.name, body) for body in bodies)
    for geomName in geomFiles:
        if not geomName.startswith('POR'):
            continue
        owner = fullCaseSetupDict[geomName].get('MOVES_WITH', ['none'])[0]
        if owner == 'none':
            continue
        if owner not in bodyNames:
            sys.exit('ERROR! [%s] MOVES_WITH %s is not a MOVE body (bodies: %s)' % (geomName, owner, list(bodyNames)))
        dz = float(bodyNames[owner].shift(startAngle))
        scale = float(geomDict[geomFiles[geomName]]['scale'])
        shifts[geomName] = dict((pid, (0.0, 0.0, dz/scale)) for pid in pidDict[geomName])
        print('\t\t%s moves with %s: shifted %+.4g m to the start angle' % (geomName, owner, dz))

    #gap refinement around the GAP_PIDS of each body
    wallGap = []
    for geomName in geomFiles:
        if geomName.startswith('WALL'):
            keywords = [k for k in fullCaseSetupDict[geomName].get('GAP_PIDS', ['']) if k != '']
            for pid in matchPids(keywords, pidDict[geomName], '[%s] GAP_PIDS' % (geomName)):
                wallGap.append(tris[geomName][pid])
    gapRefs = []
    for body, section in zip(bodies, moveSections):
        setup = fullCaseSetupDict[section]
        keywords = [k for k in setup.get('GAP_PIDS', ['']) if k != '']
        if len(keywords) == 0:
            continue
        if len(wallGap) == 0:
            sys.exit('ERROR! [%s] GAP_PIDS set but no WALL section has GAP_PIDS (the wall side of the gap)' % (section))
        dz = float(body.shift(startAngle))
        bodyGap = {pid: tris[body.source][pid] + np.array([0.0, 0.0, dz])
                   for pid in matchPids(keywords, pidDict[body.source], '[%s] GAP_PIDS' % (section))}
        gap = motionUtilities.gapWidth(np.vstack(list(bodyGap.values())), np.vstack(wallGap))
        gapCells = float(setup['GAP_CELLS'][0])
        level, distance, cell = motionUtilities.gapRefinement(gap, baseCell, gapCells)
        print('\t\t%s gap %.4g m: refinement level %d (cell %.4g m, %.1f cells across) within %.4g m'
              % (body.name, gap, level, cell, gap/cell, distance))
        gapRefs.append({'name': 'GAPREF-%s' % (body.name), 'file': 'GAPREF-%s.stl' % (body.name),
                        'distance': distance, 'level': level, 'tris': bodyGap})

    return {'zones': order, 'bodies': bodies, 'startAngle': startAngle, 'step': step, 'freq': freq,
            'nCycles': nCycles, 'shifts': shifts, 'movingPids': movingPids, 'gapRefs': gapRefs,
            'geomFiles': geomFiles}


def writeMotionGeometry(motionInfo, mshDir):
    #moved bodies: shifted copies in constant/triSurface (never write through a symlink), gap STLs
    os.makedirs(os.path.join('constant', 'triSurface'), exist_ok=True)
    for geomName, pidShift in motionInfo['shifts'].items():
        geom = motionInfo['geomFiles'][geomName]
        dest = os.path.join('constant', 'triSurface', geom.replace('.gz', ''))
        for path in (dest, os.path.join('constant', 'triSurface', geom)):
            if os.path.islink(path) or os.path.isfile(path):
                os.remove(path)
        geometryIO.translateSurfacePIDs(os.path.join(mshDir, geom), dest, pidShift)
        print('\t\t%s shifted to the start angle: %s' % (geom, ', '.join('%s %+.4g' % (pid, s[2])
                                                                          for pid, s in pidShift.items())))
    for gapRef in motionInfo['gapRefs']:
        motionUtilities.writeStlTriangles(os.path.join('constant', 'triSurface', gapRef['file']), gapRef['tris'])


def writeMotionDicts(templateLoc, motionInfo, fullCaseSetupDict):
    #tables, dynamicMeshDict, topoSetDict (motion zones), pointDisplacement, controlDictMotion
    print('\n\tWriting mesh motion dicts...')
    zones = motionInfo['zones']
    crank = getCrank(fullCaseSetupDict)
    #removing tables of the other drive from an earlier setup
    for body in motionInfo['bodies']:
        old = os.path.join('constant', 'tables', '%s%s' % (body.name, 'Table' if crank else 'Crank'))
        if os.path.isfile(old):
            os.remove(old)
    if crank:
        #tables vs crank angle, the solver works out the crank angle from the gas forces
        motionUtilities.writeCrankTables(motionInfo['bodies'], motionInfo['startAngle'], motionInfo['step'],
                                         'constant/tables')
        bodyLines = []
        for body in motionInfo['bodies']:
            patches = ' '.join('"%s_%s"' % (body.source, pid) for pid in body.movedPids)
            bodyLines.append('        %s\n        {\n            patches (%s);\n            position { type tableFile; '
                             'file "<constant>/tables/%sCrank"; outOfBounds repeat; }\n        }'
                             % (body.name, patches, body.name))
        zoneLines = ['        (%s %s)' % (zone.name, 'none' if zone.body is None else zone.body.name) for zone in zones]
        copyTemplate(templateLoc, 'constant/dynamicMeshDictCrank', 'constant/dynamicMeshDict',
                     {'BODIES': '\n'.join(bodyLines), 'ZONES': '\n'.join(zoneLines),
                      'CRANK': crankEntries(crank, motionInfo['startAngle'])})
        omega = crank['OMEGA0']
        print('\t\tcrank drive: start %.4g Hz, inertia %.4g kg m2, load %s' % (motionInfo['freq'], crank['INERTIA'],
                                                                           crank['LOAD']))
        if crank['LOAD'] == 'generator':
            rTot = crank['GEN_R_LOAD'] + crank['GEN_R_INT']
            current = crank['GEN_K']*omega/rTot
            print('\t\tgenerator at the start speed: torque %.4g N m, shaft power %.4g W, electrical power %.4g W'
                  % (crank['GEN_K']*current, crank['GEN_K']*current*omega, current**2*crank['GEN_R_LOAD']))
        elif crank['LOAD'] == 'damper':
            print('\t\tdamper at the start speed: torque %.4g N m, power %.4g W'
                  % (crank['LOAD_DAMPING']*omega, crank['LOAD_DAMPING']*omega**2))
        elif crank['LOAD'] == 'torque':
            print('\t\tload torque %.4g N m: power %.4g W at the start speed' % (crank['LOAD_TORQUE'], crank['LOAD_TORQUE']*omega))
    else:
        motionUtilities.writeTables(motionInfo['bodies'], motionInfo['startAngle'], motionInfo['freq'],
                                    motionInfo['nCycles'], motionInfo['step'], 'constant/tables')
        tableLines = ['        (%s %s)' % (zone.name, 'zeroTable' if zone.body is None else '%sTable' % (zone.body.name))
                      for zone in zones]
        copyTemplate(templateLoc, 'constant/dynamicMeshDict', 'constant/dynamicMeshDict',
                     {'INTERPOLATION_TABLES': '\n'.join(tableLines)})
    actions = []
    for zone in zones:
        actions.append('    { name %sSet; type faceSet; action new; source patchToFace; patch "%s"; }'
                       % (zone.name, zone.patch))
        actions.append('    { name %s; type faceZoneSet; action new; source setToFaceZone; faceSet %sSet; }'
                       % (zone.name, zone.name))
    with open('system/topoSetDict', 'w') as f:
        f.write('FoamFile { version 2.0; format ascii; class dictionary; object topoSetDict; }\n\n')
        f.write('//motion zones from the zone patches (written by internalSetup)\nactions\n(\n%s\n);\n' % ('\n'.join(actions)))
    copyTemplate(templateLoc, 'system/pointDisplacement', 'system/pointDisplacement', {})
    #one cycle of mesh motion only for the dry run: written MOTION_FRAMES times (mesh images), steps of 5 deg or less
    freq = motionInfo['freq']
    frames, steps = getMotionFrames(fullCaseSetupDict)
    with open('%s/defaultDicts/system/controlDict' % (templateLoc)) as f:
        text = f.read()
    replaceDict = {'<SOLVER>': 'moveDynamicMesh', '<STARTFROM>': 'startTime', '<ENDTIME>': '%.9g' % (1.0/freq),
                   '<INIT_DT>': '%.9g' % (1.0/(steps*freq)), '<WRITE_INTERVAL>': '%.9g' % (1.0/(frames*freq)),
                   '<MAX_CO>': '1e6', '<MAX_DT>': '1', '<ADJUST_DT>': 'no',
                   '#include "internalFunctionObjects"': ''}
    for key, val in replaceDict.items():
        text = text.replace(key, val)
    with open('system/controlDictMotion', 'w') as f:
        f.write(text)
    if crank:
        print('\t\t%d motion zones, crank angle tables (one cycle) from %.1f deg' % (len(zones), motionInfo['startAngle']))
    else:
        print('\t\t%d motion zones, tables for %g cycles from %.1f deg' % (len(zones), motionInfo['nCycles'],
                                                                          motionInfo['startAngle']))
    return [zone.name for zone in zones]


def getMotionFrames(fullCaseSetupDict):
    #moved meshes written per cycle in the dry run (MOTION_FRAMES, a multiple of 8 so checkMesh can run every 45 deg)
    #and the time steps per cycle (a multiple of the frames, at least 72 = 5 deg)
    try:
        frames = int(fullCaseSetupDict['MESH_MOTION'].get('MOTION_FRAMES', ['24'])[0])
    except ValueError:
        frames = 0
    if frames <= 0 or frames % 8:
        sys.exit('ERROR! MESH_MOTION MOTION_FRAMES must be a positive multiple of 8 (24 = every 15 deg, 72 = every 5 deg)')
    return frames, frames*int(np.ceil(72.0/frames))


def motionCheckTimes(fullCaseSetupDict):
    #dry run times every 45 deg (checkMesh -time)
    freq = float(getVal(fullCaseSetupDict, 'GLOBAL_CONTROL', 'FREQ'))
    return ','.join('%.9g' % (k/(8.0*freq)) for k in range(1, 9))


#----- scripts

def makeScriptsInternal(templateLoc, fullCaseSetupDict, execDir, motionZones=None, fanNames=None):
    #writing runScript, postScript, motionCheckScript and submitScript from clusterDict
    print('\n\tWriting scripts...')
    with open('%s/defaultCluster/slurm/clusterDict' % (templateLoc)) as f:
        clusterDict = json.load(f)
    solver = getVal(fullCaseSetupDict, 'GLOBAL_SIM_CONTROL', 'SOLVER')
    motionSteps = ['topoSet', 'reorderFaceZones', 'pointDisplacement', 'motionCheck', 'meshImages']
    skipSteps = list(motionSteps) if not motionZones else []
    #the crank drive library is compiled in the job (wmake skips it when up to date)
    if not motionZones or getCrank(fullCaseSetupDict) is None:
        skipSteps.append('compileLib')
    if not fanNames:
        skipSteps.append('createBaffles')
    export = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP'].get('EXPORT', ['ensight'])[0].lower()
    if export not in ('ensight', 'none'):
        sys.exit('ERROR! GLOBAL_COMPUTE_SETUP EXPORT must be ensight or none!')
    if export == 'none':
        skipSteps.append('export')
    #post job: summary, gsheet, images and videos (images need the EnSight export)
    gsheet = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP'].get('GSHEET', ['True'])[0].lower()
    if gsheet not in ('true', 'false'):
        sys.exit('ERROR! GLOBAL_COMPUTE_SETUP GSHEET must be True or False!')
    if gsheet == 'false':
        skipSteps.append('gsheet')
    if export == 'none':
        skipSteps += ['images', 'movies']
    scripts = [('run', 'runScript'), ('post', 'postScript'), ('submit', 'submitScript')]
    if motionZones:
        scripts.append(('motionCheck', 'motionCheckScript'))
    #scripts of earlier versions (separate meshing and solve jobs)
    for oldScript in ['meshingScript', 'solveScript'] + ([] if motionZones else ['motionCheckScript']):
        if os.path.isfile(oldScript):
            os.remove(oldScript)
    #setting the processor count in decomposeParDict, hierarchical/simple also get the split
    decompMethod, decomp = getDecomp(fullCaseSetupDict)
    if decompMethod in ('hierarchical', 'simple'):
        setDecomp = ("python3 '<EXEC_DIR>/meshUtilities/decompCoeffs.py' system/decomposeParDict $NP%s; "
                     "if [ $? != 0 ]; then exit 1; fi;" % (' %d %d %d' % decomp if decomp else ''))
    else:
        setDecomp = 'foamDictionary system/decomposeParDict -entry numberOfSubdomains -set $NP > /dev/null;'
    replaceDict = {'<SET_DECOMP>': setDecomp, '<CHECK_TIMES>': motionCheckTimes(fullCaseSetupDict) if motionZones else '',
                   '<ZERO_TEMPLATES>': os.path.join(execDir, 'zeroTemplates'), '<EXEC_DIR>': execDir,
                   '<SOLVER>': solver, '<MOTION_ZONES>': ' '.join(motionZones or [])}
    #slurm header from GLOBAL_COMPUTE_SETUP, job name defaults to the case directory name
    compute = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP']
    for key in compute:
        replaceDict['<%s>' % key] = ' '.join(compute[key])
    if replaceDict.get('<JOB_NAME>', 'auto') == 'auto':
        replaceDict['<JOB_NAME>'] = os.path.basename(os.getcwd())
    for section, scriptName in scripts:
        lines = [clusterDict[section][key] for key in clusterDict[section] if key not in skipSteps]
        script = '\n'.join(lines)
        for key, val in replaceDict.items():
            script = script.replace(key, val)
        with open(scriptName, 'w', newline='\n') as f:
            f.write(script + '\n')
        os.chmod(scriptName, 0o755)
        print('\t\t%s' % (scriptName))


#----- full dict extras for postRun

def hasPressureBC(wallGroups):
    #an inlet or outlet that fixes the pressure level (no pressure reference needed)
    return any(group.get('flow') and group['flow']['type'] in ('totalPressure', 'staticPressure')
               for groups in wallGroups.values() for group in groups)


def addSetupInfo(fullCaseSetupDict, wallGroups, motionInfo=None):
    #marker for postRun, resolved BC groups and the start crank angle for post
    fullCaseSetupDict['SETUP_INFO'] = {'TOOL': ['internal'], 'VERSION': [__version__],
                                        'FLOW': [getFlow(fullCaseSetupDict)],
                                        'STEADY': ['yes' if isSteady(fullCaseSetupDict) else 'no'],
                                        'RHO_REF': ['%.6g' % (getGasDensity(fullCaseSetupDict))]}
    if motionInfo:
        fullCaseSetupDict['SETUP_INFO']['START_ANGLE'] = ['%.6g' % (motionInfo['startAngle'])]
    groupDict, kindDict = {}, {}
    for geomName in wallGroups:
        for group in wallGroups[geomName]:
            groupDict[group['name']] = ['%s_%s' % (geomName, pid) for pid in group['pids']]
            #wall | inlet <type> <value> | outlet <type> <value>
            if group.get('flow'):
                kindDict[group['name']] = [group['flow']['kind']] + flowLabel(group['flow']).split()
            else:
                kindDict[group['name']] = ['wall']
    fullCaseSetupDict['BC_GROUPS_RESOLVED'] = groupDict
    fullCaseSetupDict['BC_GROUP_KINDS'] = kindDict
    return fullCaseSetupDict
