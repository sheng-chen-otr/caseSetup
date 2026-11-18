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


def parseGroupLine(tokens, where):
    #<group> = PIDS kw1 kw2 .. [U noSlip|slip] [T mode args]
    tokens = [t for t in tokens if t != '']
    parts = {'PIDS': [], 'U': [], 'T': []}
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
    group = {'keywords': parts['PIDS'], 'U': None, 'T': None}
    if parts['U']:
        if len(parts['U']) != 1 or parts['U'][0] not in uModes:
            sys.exit('ERROR! %s: U must be one of %s!' % (where, uModes))
        group['U'] = parts['U'][0]
    if parts['T']:
        group['T'] = parseThermal(parts['T'], where)
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
        groups.append({'name': groupName, 'pids': matched,
                       'U': group['U'] or uDefault, 'T': group['T'] or tDefault})
    rest = [pid for pid in pids if pid not in owner]
    if len(rest) > 0:
        groups.append({'name': '%s_default' % (wallName), 'pids': rest, 'U': uDefault, 'T': tDefault})
    return groups


def getWallGroups(geomDict, fullCaseSetupDict, pidDict):
    #wallGroups[geomName] = groups, prints the PID -> group -> U, T table
    wallGroups = {}
    print('\n\tBoundary condition groups:')
    print('\t\t%-28s %-20s %-8s %s' % ('patch', 'group', 'U', 'T'))
    for geom in geomDict:
        geomName = stripExt(geom)
        if not geomName.startswith('WALL'):
            continue
        groups = resolveBCGroups(geomName, fullCaseSetupDict[geomName], pidDict[geomName])
        for group in groups:
            for pid in group['pids']:
                print('\t\t%-28s %-20s %-8s %s' % ('%s_%s' % (geomName, pid), group['name'], group['U'],
                                                   ' '.join(str(t) for t in group['T'])))
        wallGroups[geomName] = groups
    return wallGroups


#----- system dicts

def getBackgroundMesh(bbox, baseCell, nMargin=2):
    #bbox (xmin ymin zmin xmax ymax zmax) grown by nMargin cells and rounded to whole cells
    lo = np.floor(np.array(bbox[:3])/baseCell - nMargin)*baseCell
    hi = np.ceil(np.array(bbox[3:])/baseCell + nMargin)*baseCell
    n = np.rint((hi - lo)/baseCell).astype(int)
    return lo, hi, n


def writeBlockMeshInternal(templateLoc, geomDict, fullCaseSetupDict, mshDir):
    print('\n\tWriting blockMeshDict...')
    baseCell = float(getVal(fullCaseSetupDict, 'GLOBAL_REFINEMENT', 'BASE_CELL'))
    allBox = []
    for geom in geomDict:
        scale = float(geomDict[geom]['scale'])
        allBox.append(np.array(geometryIO.surfaceBoundingBox(os.path.join(mshDir, geom)))*scale)
    allBox = np.array(allBox)
    bbox = list(allBox[:, :3].min(axis=0)) + list(allBox[:, 3:].max(axis=0))
    lo, hi, n = getBackgroundMesh(bbox, baseCell)
    print('\t\tbackground box %s to %s, cells %s (base cell %g m)' % (lo, hi, n, baseCell))
    replaceDict = {'XMIN': '%.9g' % lo[0], 'YMIN': '%.9g' % lo[1], 'ZMIN': '%.9g' % lo[2],
                   'XMAX': '%.9g' % hi[0], 'YMAX': '%.9g' % hi[1], 'ZMAX': '%.9g' % hi[2],
                   'NX': n[0], 'NY': n[1], 'NZ': n[2]}
    copyTemplate(templateLoc, 'system/blockMeshDict', 'system/blockMeshDict', replaceDict)
    return bbox


def writeSnappyInternal(templateLoc, geomDict, fullCaseSetupDict):
    print('\n\tWriting snappyHexMeshDict...')
    refinement = fullCaseSetupDict['GLOBAL_REFINEMENT']
    locInMesh = [v for v in refinement['LOC_IN_MESH'] if v != '']
    if len(locInMesh) != 3:
        sys.exit('ERROR! GLOBAL_REFINEMENT LOC_IN_MESH must be 3 coordinates of a point inside the gas!')
    try:
        [float(v) for v in locInMesh]
    except ValueError:
        sys.exit('ERROR! GLOBAL_REFINEMENT LOC_IN_MESH must be 3 numbers, got %s!' % (locInMesh))
    geometryArray, surfaceArray, layerArray = [], [], []
    for geom in geomDict:
        geomName = stripExt(geom)
        geomFile = geom.replace('.gz', '')
        level = geomDict[geom]['refinement']
        geometryArray.append('    %s { type triSurfaceMesh; name %s; file "%s"; scale %s; }'
                             % (geomFile, geomName, geomFile, geomDict[geom]['scale']))
        surfaceArray.append('        %s { level (%s %s); patchInfo { type wall; } }' % (geomName, level, level))
        layerArray.append('        "%s_.*" { nSurfaceLayers %s; expansionRatio %s; }'
                          % (geomName, geomDict[geom]['layers'], geomDict[geom]['expansion']))
    #gap faces: no layers, written after the file-wide entry (last matching regex wins)
    for geom in geomDict:
        geomName = stripExt(geom)
        gapPids = [p for p in fullCaseSetupDict[geomName].get('GAP_PIDS', ['']) if p != '']
        for kw in gapPids:
            layerArray.append('        "%s_.*%s.*" { nSurfaceLayers 0; }' % (geomName, kw))
    replaceDict = {'GEOMETRY': '\n'.join(geometryArray),
                   'REFINEMENT_SURFACES': '\n'.join(surfaceArray),
                   'REFINEMENT_REGIONS': '',
                   'LAYERS': '\n'.join(layerArray),
                   'LOC_IN_MESH': ' '.join(locInMesh),
                   'RESOLVE_FEAT_ANGLE': refinement['RESOLVE_FEAT_ANGLE'][0],
                   'N_CELLS_BETWEEN_LEVELS': refinement['N_CELLS_BETWEEN_LEVELS'][0],
                   'FINAL_LAYER_THICKNESS': refinement['FINAL_LAYER_THICKNESS'][0],
                   'MIN_LAYER_THICKNESS': refinement['MIN_LAYER_THICKNESS'][0]}
    copyTemplate(templateLoc, 'system/snappyHexMeshDict', 'system/snappyHexMeshDict', replaceDict)


def getTurbulence(fullCaseSetupDict):
    #laminar (default) or a RAS model name
    model = getVal(fullCaseSetupDict, 'GLOBAL_SIM_CONTROL', 'TURB_MODEL')
    aliases = {'laminar': 'laminar', 'kosst': 'kOmegaSST', 'komegasst': 'kOmegaSST',
               'kepsilon': 'kEpsilon', 'komega': 'kOmega'}
    if model.lower() not in aliases:
        sys.exit('ERROR! TURB_MODEL %s is not supported for internal cases, use one of %s!'
                 % (model, sorted(set(aliases.values()))))
    return aliases[model.lower()]


def writeCasePropertiesInternal(templateLoc, geomDict, fullCaseSetupDict, wallGroups):
    print('\n\tWriting caseProperties...')
    material = fullCaseSetupDict['GLOBAL_MATERIAL']
    pInit = float(material['CHARGE_PRESSURE'][0])
    tInit = float(material['INIT_T'][0])
    initArray = ['U uniform (0 0 0);', 'p uniform %.9g;' % (pInit), 'T uniform %.9g;' % (tInit)]
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
            patches = ' '.join('"%s_%s"' % (geomName, pid) for pid in group['pids'])
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
            if group['U'] == 'slip':
                options = 'heatTransfer %s;' % (heat)
            else:
                options = 'wallFunction %s; motion stationary; heatTransfer %s;' % (wallFunctionNames[wallModel], heat)
            bcArray.append('    %s { category wall; type %s; patches (%s); options { %s } values { %s } }'
                           % (group['name'], group['U'], patches, options, values))
    with open('%s/defaultDicts/system/caseProperties' % (templateLoc)) as f:
        text = f.read()
    text = text.replace('<INITIAL_CONDITIONS>', '\n'.join('    ' + s for s in initArray))
    text = text.replace('<BOUNDARY_CONDITIONS>', '\n'.join(bcArray))
    with open('system/caseProperties', 'w') as f:
        f.write(text)


def writeThermoInternal(templateLoc, fullCaseSetupDict):
    print('\n\tWriting thermophysicalProperties...')
    props = getGasProperties(fullCaseSetupDict)
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


def writeControlDictInternal(templateLoc, fullCaseSetupDict):
    print('\n\tWriting controlDict...')
    control = fullCaseSetupDict['GLOBAL_CONTROL']
    endTime, writeInterval = getRunTimes(fullCaseSetupDict)
    print('\t\t%s cycles at %s Hz: endTime %.6g s, write every %.6g s'
          % (control['N_CYCLES'][0], control['FREQ'][0], endTime, writeInterval))
    replaceDict = {'SOLVER': getVal(fullCaseSetupDict, 'GLOBAL_SIM_CONTROL', 'SOLVER'),
                   'STARTFROM': control['STARTFROM'][0], 'ENDTIME': '%.9g' % endTime,
                   'INIT_DT': control['INIT_DT'][0], 'WRITE_INTERVAL': '%.9g' % writeInterval,
                   'MAX_CO': control['MAX_CO'][0], 'MAX_DT': control['MAX_DT'][0]}
    copyTemplate(templateLoc, 'system/controlDict', 'system/controlDict', replaceDict)


def writeNumericsInternal(templateLoc, fullCaseSetupDict):
    print('\n\tWriting fvSchemes, fvSolution, decomposeParDict...')
    simControl = fullCaseSetupDict['GLOBAL_SIM_CONTROL']
    compute = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP']
    copyTemplate(templateLoc, 'system/fvSchemes', 'system/fvSchemes', {})
    copyTemplate(templateLoc, 'system/fvSolution', 'system/fvSolution',
                 {'N_OUTER_CORR': simControl['N_OUTER_CORR'][0], 'N_CORR': simControl['N_CORR'][0],
                  'N_NONORTH_CORR': simControl['N_NONORTH_CORR'][0]})
    copyTemplate(templateLoc, 'system/decomposeParDict', 'system/decomposeParDict',
                 {'NUM_PROCS': compute['NUM_PROCS'][0], 'DECOMP_METHOD': compute['DECOMP_METHOD'][0]})


def writeFunctionObjectsInternal(wallGroups):
    #gas averages and mass (mass must stay constant in a closed volume), wall heat flux per patch
    #gas volume for post: V = volIntegrate(rho)/volAverage(rho), no extra function object needed
    print('\n\tWriting internalFunctionObjects...')
    allPatches = []
    for geomName in wallGroups:
        for group in wallGroups[geomName]:
            allPatches += ['"%s_%s"' % (geomName, pid) for pid in group['pids']]
    text = '''gasAvg
{
    type            volFieldValue;
    libs            (fieldFunctionObjects);
    regionType      all;
    operation       volAverage;
    fields          (p T rho);
    writeFields     false;
    writeControl    timeStep;
    writeInterval   1;
}
gasMass
{
    type            volFieldValue;
    libs            (fieldFunctionObjects);
    regionType      all;
    operation       volIntegrate;
    fields          (rho);
    writeFields     false;
    writeControl    timeStep;
    writeInterval   1;
}
wallHeatFlux
{
    type            wallHeatFlux;
    libs            (fieldFunctionObjects);
    patches         (%s);
    writeFields     false;
    writeControl    timeStep;
    writeInterval   1;
}
''' % (' '.join(allPatches))
    with open('system/internalFunctionObjects', 'w') as f:
        f.write(text)


#----- scripts

def makeScriptsInternal(templateLoc, fullCaseSetupDict, zeroTemplatesDir):
    print('\n\tWriting meshingScript, solveScript...')
    with open('%s/defaultCluster/slurm/clusterDict' % (templateLoc)) as f:
        clusterDict = json.load(f)
    solver = getVal(fullCaseSetupDict, 'GLOBAL_SIM_CONTROL', 'SOLVER')
    for section, scriptName in [('meshing', 'meshingScript'), ('solve', 'solveScript')]:
        script = '\n'.join(clusterDict[section][line] for line in clusterDict[section])
        script = script.replace('<ZERO_TEMPLATES>', zeroTemplatesDir).replace('<SOLVER>', solver)
        with open(scriptName, 'w') as f:
            f.write(script + '\n')
        os.chmod(scriptName, 0o755)
        print('\t\t%s' % (scriptName))


#----- full dict extras for postRun

def addSetupInfo(fullCaseSetupDict, wallGroups):
    #marker read by postRun (internal vs aero) and the resolved BC groups for post
    fullCaseSetupDict['SETUP_INFO'] = {'TOOL': ['internal'], 'VERSION': [__version__]}
    groupDict = {}
    for geomName in wallGroups:
        for group in wallGroups[geomName]:
            groupDict[group['name']] = ['%s_%s' % (geomName, pid) for pid in group['pids']]
    fullCaseSetupDict['BC_GROUPS_RESOLVED'] = groupDict
    return fullCaseSetupDict
