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
import glob
import shutil
import hashlib

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


def resolveBCGroups(wallName, wallSetup, pids, fromBodies=()):
    #matches the BC_GROUPS keywords (case-insensitive substring) to the PIDs of one WALL geometry
    #returns an ordered list of groups {name, pids, U, T}; unmatched PIDs go to <wall>_default
    #fromBodies: (section, ZONE_PIDS keywords, U_DEFAULT, T_DEFAULT) of the FROM bodies in this file, each its own
    #group over its faces (blank U / T = this file's defaults)
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
    for section, keywords, uOwn, tOwn in fromBodies:
        matched = [pid for pid in pids if any(kw.lower() in pid.lower() for kw in keywords)]
        for pid in matched:
            if pid in owner:
                sys.exit('ERROR! PID "%s" of [%s] (FROM %s) is also in [%s] group %s: set its U and T in one place, '
                         'the T_DEFAULT / U_DEFAULT of [%s] or that group' % (pid, section, wallName, wallName, owner[pid], section))
            owner[pid] = section
        uOwn = [u for u in uOwn if u != '']
        if uOwn and uOwn[0] not in uModes:
            sys.exit('ERROR! [%s] U_DEFAULT must be one of %s or blank!' % (section, uModes))
        tOwn = [t for t in tOwn if t != '']
        if matched:
            groups.append({'name': section, 'pids': matched, 'flow': None, 'U': uOwn[0] if uOwn else uDefault,
                           'T': parseThermal(tOwn, '[%s] T_DEFAULT' % (section)) if tOwn else tDefault})
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
        #FROM bodies in this file (a piston crown in the chamber): their faces get the body's own U and T
        fromBodies = [(section, [k for k in body.get('ZONE_PIDS', ['']) if k != ''], body.get('U_DEFAULT', ['']),
                       body.get('T_DEFAULT', ['']))
                      for section, body in fullCaseSetupDict.items()
                      if section.startswith('MOVE') and ''.join(body.get('FROM', [''])).strip() == geomName]
        groups = resolveBCGroups(geomName, fullCaseSetupDict[geomName], pidDict[geomName], fromBodies)
        #since 2026-11-24 a FROM body's own T_DEFAULT is used (before, its faces took the file's default): say so when
        #they differ, so a case set up before and after does not change its crown temperature unnoticed
        fileT = parseThermal([t for t in fullCaseSetupDict[geomName]['T_DEFAULT'] if t != ''],
                             '[%s] T_DEFAULT' % (geomName))
        for group in groups:
            if group['name'] in [section for section, _, _, _ in fromBodies] and group['T'] != fileT:
                print('\t\tNOTE: [%s] faces use its own T_DEFAULT (%s), not the %s default (%s)'
                      % (group['name'], ' '.join(str(t) for t in group['T']), geomName, ' '.join(str(t) for t in fileT)))
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
    #BC_GROUPS of a FROM body are not used: its faces are one group, set by its own U_DEFAULT / T_DEFAULT
    for section, body in fullCaseSetupDict.items():
        if not section.startswith('MOVE') or not ''.join(body.get('FROM', [''])).strip():
            continue
        ownGroups = [g for g in body.get('BC_GROUPS', ['']) if g != '']
        if ownGroups:
            print('\t\tWARNING: [%s] BC_GROUPS %s is not used (FROM body): its faces take its T_DEFAULT / U_DEFAULT'
                  % (section, ' '.join(ownGroups)))
    return wallGroups


#----- mesh groups: surface level and wall layers per PID

#LAYERS auto: first layer = thermal penetration depth / AUTO_FIRST_DIV, layers growing at most AUTO_EXPANSION up to the
#final layer (FINAL x surface cell)
AUTO_FIRST_DIV = 10.0
AUTO_EXPANSION = 1.3
#and the whole stack at most AUTO_TOTAL_DEPTHS penetration depths thick (absolute sizes then)
AUTO_TOTAL_DEPTHS = 2.0

#layer sizes of a mesh group line and their snappy keys; any two pick the thicknessModel
layerSizeKeys = ('FIRST', 'FINAL', 'TOTAL', 'EXPANSION')
snappySizeKeys = {'FIRST': 'firstLayerThickness', 'FINAL': 'finalLayerThickness', 'TOTAL': 'thickness',
                  'EXPANSION': 'expansionRatio'}
thicknessModels = {('FIRST', 'TOTAL'): 'firstAndOverall', ('EXPANSION', 'FIRST'): 'firstAndExpansion',
                   ('FINAL', 'TOTAL'): 'finalAndOverall', ('EXPANSION', 'FINAL'): 'finalAndExpansion',
                   ('EXPANSION', 'TOTAL'): 'overallAndExpansion'}
#SIZES: relative = fractions of the local cell, absolute = m, mixed = FIRST in m and FINAL relative
sizeModes = ('relative', 'absolute', 'mixed')


def parseMeshGroupLine(tokens, where):
    #<group> = PIDS kw1 kw2 .. [LEVEL a [b]] [EDGE_LEVEL n] [LAYERS n|auto|0|-1] [SIZES relative|absolute|mixed]
    #          [FIRST x] [FINAL x] [TOTAL x] [EXPANSION r] [MIN x]
    tokens = [t for t in tokens if t != '']
    parts = dict((key, []) for key in ('PIDS', 'LEVEL', 'EDGE_LEVEL', 'LAYERS', 'SIZES', 'MIN') + layerSizeKeys)
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
    group = {'keywords': parts['PIDS'], 'level': None, 'edgeLevel': None, 'layers': None, 'mode': None, 'sizes': {},
             'min': None}
    if parts['LEVEL']:
        if len(parts['LEVEL']) > 2 or not all(v.isdigit() for v in parts['LEVEL']):
            sys.exit('ERROR! %s: LEVEL takes one or two refinement levels (min max), got %s!' % (where, parts['LEVEL']))
        level = [int(v) for v in parts['LEVEL']]
        group['level'] = (level[0], level[-1])
        if level[-1] < level[0]:
            sys.exit('ERROR! %s: LEVEL max below min!' % (where))
    if parts['EDGE_LEVEL']:
        if len(parts['EDGE_LEVEL']) != 1 or not parts['EDGE_LEVEL'][0].isdigit():
            sys.exit('ERROR! %s: EDGE_LEVEL takes one refinement level, got %s!' % (where, parts['EDGE_LEVEL']))
        group['edgeLevel'] = int(parts['EDGE_LEVEL'][0])
    if parts['LAYERS']:
        value = parts['LAYERS'][0]
        if len(parts['LAYERS']) != 1 or not (value == 'auto' or value.isdigit() or value == '-1'):
            sys.exit('ERROR! %s: LAYERS takes a number of layers, auto, 0 (none, wall points fixed) or -1 (none, wall '
                     'points slide), got %s!' % (where, parts['LAYERS']))
        group['layers'] = value if value == 'auto' else int(value)
    for key in layerSizeKeys + ('MIN',):
        if parts[key]:
            values = parseNumbers(parts[key], where, key)
            if len(values) != 1 or values[0] <= 0:
                sys.exit('ERROR! %s: %s takes one positive number!' % (where, key))
            if key == 'MIN':
                group['min'] = values[0]
            else:
                group['sizes'][key] = values[0]
    if parts['SIZES']:
        if len(parts['SIZES']) != 1 or parts['SIZES'][0] not in sizeModes:
            sys.exit('ERROR! %s: SIZES must be one of %s!' % (where, ' | '.join(sizeModes)))
        group['mode'] = parts['SIZES'][0]
    given = group['sizes'] or group['mode'] or group['min']
    if group['level'] is None and group['edgeLevel'] is None and group['layers'] is None:
        sys.exit('ERROR! %s sets nothing: give LEVEL, EDGE_LEVEL and/or LAYERS!' % (where))
    if given and group['layers'] is None:
        sys.exit('ERROR! %s: SIZES / FIRST / FINAL / TOTAL / EXPANSION / MIN need LAYERS!' % (where))
    if given and group['layers'] in (0, -1):
        sys.exit('ERROR! %s: LAYERS %d adds no layers, no sizes!' % (where, group['layers']))
    if group['layers'] in (None, 0, -1):
        return group
    #mixed by default with FIRST or auto (absolute first layer blending into the cells), else relative as the GEOM line
    if group['mode'] is None:
        group['mode'] = 'mixed' if group['layers'] == 'auto' or 'FIRST' in group['sizes'] else 'relative'
    sizes = sorted(group['sizes'])
    if len(sizes) > 2:
        sys.exit('ERROR! %s: give at most two of FIRST, FINAL, TOTAL, EXPANSION (got %s)!' % (where, ' '.join(sizes)))
    if group['layers'] == 'auto':
        if group['mode'] != 'mixed' or any(key != 'FINAL' for key in sizes):
            sys.exit('ERROR! %s: LAYERS auto sets the first layer itself (SIZES mixed, FINAL only)!' % (where))
    elif group['mode'] == 'mixed':
        if 'FIRST' not in sizes or any(key not in ('FIRST', 'FINAL') for key in sizes):
            sys.exit('ERROR! %s: SIZES mixed takes FIRST (m) and optionally FINAL (x the cell); for TOTAL / EXPANSION '
                     'use SIZES absolute or relative!' % (where))
    elif group['mode'] == 'absolute' and not any(key in sizes for key in ('FIRST', 'FINAL', 'TOTAL')):
        sys.exit('ERROR! %s: SIZES absolute needs FIRST, FINAL or TOTAL in m!' % (where))
    if 'EXPANSION' in sizes and group['sizes']['EXPANSION'] < 1.0:
        print('\t\tWARNING: %s: EXPANSION %g < 1, the layers get thinner away from the wall'
              % (where, group['sizes']['EXPANSION']))
    return group


def thermalDepth(fullCaseSetupDict, wallT):
    #thermal penetration depth sqrt(2 alpha/omega) of the gas at the charge pressure and wallT for a cyclic run at FREQ
    props = getGasProperties(fullCaseSetupDict)
    material = fullCaseSetupDict['GLOBAL_MATERIAL']
    mu = props['MU'] if props['TRANSPORT'] == 'const' else props['AS']*wallT**1.5/(wallT + props['TS'])
    rho = float(material['CHARGE_PRESSURE'][0])*props['MOL_WEIGHT']/(8314.46*wallT)
    alpha = mu/props['PR']/rho
    omega = 2.0*np.pi*float(getVal(fullCaseSetupDict, 'GLOBAL_CONTROL', 'FREQ'))
    return np.sqrt(2.0*alpha/omega)


def stackGeometry(sizes, n):
    #first, final, total thickness and growth ratio of n layers from two of FIRST / FINAL / TOTAL / EXPANSION (one length unit)
    span = lambda r: float(n) if abs(r - 1.0) < 1e-12 else (r**n - 1.0)/(r - 1.0)
    if 'EXPANSION' in sizes:
        ratio = sizes['EXPANSION']
    elif 'FIRST' in sizes and 'FINAL' in sizes:
        ratio = (sizes['FINAL']/sizes['FIRST'])**(1.0/(n - 1)) if n > 1 else 1.0
    else:
        #FIRST or FINAL with TOTAL: growth ratio by bisection (total grows with the ratio for a fixed first layer, and
        #falls for a fixed final layer)
        lo, hi = 1e-3, 1e3
        for _ in range(200):
            ratio = (lo*hi)**0.5
            if 'FIRST' in sizes:
                big = sizes['FIRST']*span(ratio) > sizes['TOTAL']
            else:
                big = sizes['FINAL']*span(ratio)/ratio**(n - 1) < sizes['TOTAL']
            lo, hi = (lo, ratio) if big else (ratio, hi)
    if 'FIRST' in sizes:
        first = sizes['FIRST']
    elif 'FINAL' in sizes:
        first = sizes['FINAL']/ratio**(n - 1)
    else:
        first = sizes['TOTAL']/span(ratio)
    return first, first*ratio**(n - 1), first*span(ratio), ratio


def layerEntry(group, cell, expansionDefault, finalDefault, minDefault):
    #snappy layers entry of a mesh group with layers > 0, and the stack in m (first, final, total, ratio)
    #every entry names its thicknessModel: without it a patch entry only reads the keys of the global model
    #(v2412: a per-patch firstLayerThickness under the global finalAndExpansion was ignored without a message)
    n = group['layers']
    sizes = dict(group['sizes'])
    if group['mode'] == 'mixed':
        final = sizes.get('FINAL', finalDefault)
        first, finalM, total, ratio = stackGeometry({'FIRST': sizes['FIRST'], 'FINAL': final*cell}, n)
        entry = ('nSurfaceLayers %d; thicknessModel firstAndRelativeFinal; relativeSizes true; firstLayerThickness %.6g; '
                 'finalLayerThickness %.6g; minThickness %.6g;' % (n, sizes['FIRST'], final, group['min'] or minDefault))
        return entry, (first, finalM, total, ratio)
    #relative / absolute: fill up to two sizes (FINAL default only when no thickness is given, relative)
    if not any(key in sizes for key in ('FIRST', 'FINAL', 'TOTAL')):
        sizes['FINAL'] = finalDefault
    if len(sizes) < 2:
        sizes['EXPANSION'] = expansionDefault
    scale = cell if group['mode'] == 'relative' else 1.0
    first, finalM, total, ratio = stackGeometry(dict((key, val*(1.0 if key == 'EXPANSION' else scale))
                                                     for key, val in sizes.items()), n)
    #FIRST and FINAL have no snappy model of their own: first layer and the growth ratio between them
    if 'FIRST' in sizes and 'FINAL' in sizes:
        sizes = {'FIRST': sizes['FIRST'], 'EXPANSION': ratio}
    model = thicknessModels[tuple(sorted(sizes))]
    #minThickness in the group's units: relative MIN_LAYER_THICKNESS, absolute a tenth of the stack
    minThickness = group['min'] or (minDefault if group['mode'] == 'relative' else 0.1*total)
    entry = ('nSurfaceLayers %d; thicknessModel %s; relativeSizes %s; %s minThickness %.6g;'
             % (n, model, 'true' if group['mode'] == 'relative' else 'false',
                ' '.join('%s %.6g;' % (snappySizeKeys[key], sizes[key]) for key in sorted(sizes)), minThickness))
    return entry, (first, finalM, total, ratio)


def pidWallT(wallGroups, geomName, pid, initT):
    #fixed wall T of the PID's BC group, else INIT_T
    for group in wallGroups.get(geomName, []):
        if pid in group['pids'] and group['T'] and group['T'][0] == 'fixed':
            return group['T'][1]
    return initT


def getMeshGroups(geomDict, fullCaseSetupDict, pidDict, wallGroups, motionInfo=None):
    #meshGroups[geomName] = [{name, pids, level, layers, mode, entry, stack, ...}], prints the PID -> mesh group table
    #PIDs match like BC_GROUPS (case-insensitive substring); a PID in two mesh groups is an error
    refinement = fullCaseSetupDict['GLOBAL_REFINEMENT']
    baseCell = float(getVal(fullCaseSetupDict, 'GLOBAL_REFINEMENT', 'BASE_CELL'))
    finalDefault = float(refinement['FINAL_LAYER_THICKNESS'][0])
    minDefault = float(refinement['MIN_LAYER_THICKNESS'][0])
    initT = float(getVal(fullCaseSetupDict, 'GLOBAL_MATERIAL', 'INIT_T'))
    cyclic = not isSteady(fullCaseSetupDict) and getEndTime(fullCaseSetupDict) is None
    meshGroups = {}
    rows = []
    for geom in geomDict:
        geomName = stripExt(geom)
        if not geomName.startswith(('WALL', 'MOVE')):
            continue
        setup = fullCaseSetupDict[geomName]
        groupNames = [g for g in setup.get('MESH_GROUPS', ['']) if g != '']
        owner = {}
        groups = []
        for groupName in groupNames:
            where = '[%s] %s' % (geomName, groupName)
            if groupName not in setup:
                sys.exit('ERROR! [%s] MESH_GROUPS lists %s but there is no "%s = PIDS ..." line!'
                         % (geomName, groupName, groupName))
            if groupName in [g for g in setup.get('BC_GROUPS', ['']) if g != '']:
                sys.exit('ERROR! [%s] %s is in BC_GROUPS and MESH_GROUPS: give the mesh group its own name'
                         % (geomName, groupName))
            group = parseMeshGroupLine(setup[groupName], where)
            group['name'] = groupName
            group['pids'] = []
            for kw in group['keywords']:
                hits = [pid for pid in pidDict[geomName] if kw.lower() in pid.lower()]
                if len(hits) == 0:
                    sys.exit('ERROR! %s: keyword "%s" matches no PID. PIDs are: %s' % (where, kw, pidDict[geomName]))
                for pid in hits:
                    if pid in owner and owner[pid] != groupName:
                        sys.exit('ERROR! [%s] PID "%s" matches mesh groups %s and %s!' % (geomName, pid, owner[pid], groupName))
                    owner[pid] = groupName
                    if pid not in group['pids']:
                        group['pids'].append(pid)
            #the mesh group entry is written after the GAP_PIDS one, so it wins there
            gapKeywords = [k for k in setup.get('GAP_PIDS', ['']) if k != '']
            gapHits = [pid for pid in group['pids'] if any(k.lower() in pid.lower() for k in gapKeywords)]
            if gapHits and group['layers'] not in (None, 0, -1):
                print('\t\tWARNING: %s puts layers on the GAP_PIDS faces %s (layers collapse in narrow gaps)'
                      % (where, gapHits))
            #edges of faces opening into the refined side gap: the rim cells are only cube-like (no highly skewed
            #faces) one level above the gap refinement (v2412 demo segment: displacer edges at the gap level 3 left
            #752 faces with skewness > 4, level 4 none)
            gapLevels = [gapRef['level'] for gapRef in (motionInfo['gapRefs'] if motionInfo else [])]
            if group['edgeLevel'] is not None and gapLevels and group['edgeLevel'] <= max(gapLevels):
                print('\t\tWARNING: %s: EDGE_LEVEL %d is not above the gap refinement level %d: edges opening into '
                      'the gap or an end-gap band keep flat, skewed cells, use EDGE_LEVEL %d'
                      % (where, group['edgeLevel'], max(gapLevels), max(gapLevels) + 1))
            #surface cell size the layers grow into: the group's min level, else the file's
            level = group['level'][0] if group['level'] else int(geomDict[geom]['refinement'])
            cell = baseCell/2**level
            group['cell'] = cell
            group['depth'] = None
            group['entry'], group['stack'] = None, None
            if group['layers'] == 'auto':
                if not cyclic:
                    sys.exit('ERROR! %s: LAYERS auto sizes the layers from the thermal penetration depth at FREQ '
                             '(cyclic runs only): give LAYERS n with FIRST' % (where))
                #coldest wall of the group = thinnest thermal layer (alpha grows with T)
                group['depth'] = min(thermalDepth(fullCaseSetupDict, pidWallT(wallGroups, geomName, pid, initT))
                                     for pid in group['pids'])
                group['sizes']['FIRST'] = group['depth']/AUTO_FIRST_DIV
                final = group['sizes'].get('FINAL', finalDefault)*cell
                if final <= group['sizes']['FIRST']:
                    #no layers needed: written like the walls without layers (ZERO_LAYERS)
                    group['layers'] = None
                    group['autoNone'] = True
                    print('\t\tNOTE: %s: surface cells (%.3g mm) already fine for the thermal layer (%.3g mm), no layers'
                          % (where, 1e3*cell, 1e3*group['depth']))
                else:
                    group['layers'] = int(np.ceil(np.log(final/group['sizes']['FIRST'])/np.log(AUTO_EXPANSION))) + 1
                    #stack capped at AUTO_TOTAL_DEPTHS thermal depths: a stack growing up to the cell size (3.3 mm for a
                    #0.84 mm layer on 2 mm cells) did not fit the end gap, snappy shrank the mesh to make room and
                    #inverted cells there (v2412 demo segment: 60 negative volumes); capped at 1.7 mm: none, all layers
                    first = group['sizes']['FIRST']
                    total = stackGeometry({'FIRST': first, 'FINAL': final}, group['layers'])[2]
                    cap = AUTO_TOTAL_DEPTHS*group['depth']
                    if total > cap:
                        group['mode'] = 'absolute'
                        group['sizes'] = {'FIRST': first, 'TOTAL': cap}
                        group['layers'] = int(np.ceil(np.log(1.0 + cap/first*(AUTO_EXPANSION - 1.0))/np.log(AUTO_EXPANSION)))
            if group['layers'] in (0, -1):
                group['entry'] = 'nSurfaceLayers %d;' % (group['layers'])
            elif group['layers'] is not None:
                group['entry'], group['stack'] = layerEntry(group, cell, float(geomDict[geom]['expansion']),
                                                            finalDefault, minDefault)
                first, final, total, ratio = group['stack']
                if ratio < 1.0:
                    print('\t\tWARNING: %s: the layers get thinner away from the wall (%.3g mm to %.3g mm)'
                          % (where, 1e3*first, 1e3*final))
                elif ratio > 1.5:
                    print('\t\tWARNING: %s: layers grow %.2fx per layer from %.3g mm to %.3g mm: use more LAYERS'
                          % (where, ratio, 1e3*first, 1e3*final))
                if final > cell:
                    print('\t\tWARNING: %s: final layer %.3g mm is thicker than the %.3g mm surface cell'
                          % (where, 1e3*final, 1e3*cell))
            groups.append(group)
        for group in groups:
            for pid in group['pids']:
                rows.append((geomName, pid, group))
        meshGroups[geomName] = groups
    #MESH_GROUPS of a FROM body are not used: its faces are PIDs of another file, set them in that file's section
    for section, body in fullCaseSetupDict.items():
        if section.startswith('MOVE') and ''.join(body.get('FROM', [''])).strip() not in ('', 'self'):
            ownGroups = [g for g in body.get('MESH_GROUPS', ['']) if g != '']
            if ownGroups:
                print('\t\tWARNING: [%s] MESH_GROUPS %s is not used (FROM body): set its faces in the MESH_GROUPS of [%s]'
                      % (section, ' '.join(ownGroups), ''.join(body['FROM']).strip()))
    if rows:
        print('\n\tMesh groups (layer sizes in mm at the surface cell size):')
        print('\t\t%-32s %-14s %-6s %-5s %-10s %-9s %-8s %-8s %-6s %-8s %s'
              % ('patch', 'group', 'level', 'edge', 'layers', 'sizes', 'first', 'final', 'ratio', 'total', 'thermal layer'))
        for geomName, pid, group in rows:
            level = ('%d %d' % group['level'] if group['level'] else '-').ljust(6) + ' ' + \
                    ('%d' % group['edgeLevel'] if group['edgeLevel'] is not None else '-').ljust(5)
            if group['stack']:
                first, final, total, ratio = group['stack']
                layer = ('%-10d %-9s %-8.3g %-8.3g %-6.2f %-8.3g'
                         % (group['layers'], group['mode'], 1e3*first, 1e3*final, ratio, 1e3*total))
            elif group['layers'] in (0, -1):
                layer = '%-10s' % ('0 (fixed)' if group['layers'] == 0 else '-1 (slip)')
            else:
                layer = '%-10s' % ('none' if group.get('autoNone') else 'GEOM')
            print(('\t\t%-32s %-14s %s %s %s' % ('%s_%s' % (geomName, pid), group['name'], level, layer,
                                                   '%.3g' % (1e3*group['depth']) if group['depth'] else '')).rstrip())
    return meshGroups


def getEdgeRefs(meshGroups, geomDict, mshDir):
    #EDGE_LEVEL groups: the triangles of their PIDs (file units) for EDGE-<geometry>-<group>.stl, whose feature edges
    #(the outline of the faces and their sharp edges) snappy refines at EDGE_LEVEL; sets group['edgeFile']
    edgeRefs = []
    for geom in geomDict:
        geomName = stripExt(geom)
        groups = [group for group in meshGroups.get(geomName, []) if group['edgeLevel'] is not None]
        if not groups:
            continue
        tris = geometryIO.readSurfaceTriangles(os.path.join(mshDir, geom))
        for group in groups:
            group['edgeFile'] = 'EDGE-%s-%s.stl' % (geomName, group['name'])
            edgeRefs.append({'file': group['edgeFile'], 'geomName': geomName,
                             'baseTris': dict((pid, tris[pid]) for pid in group['pids'])})
    return edgeRefs


def writeEdgeGeometry(edgeRefs, triDir=os.path.join('constant', 'triSurface'), shifts=None):
    #EDGE-*.stl of the EDGE_LEVEL groups in triDir, PIDs of moving bodies shifted like their geometry ({geometry: {pid:
    #(dx, dy, dz)}} in file units); old EDGE files of an earlier setup are removed. returns the files written
    shifts = shifts or {}
    os.makedirs(triDir, exist_ok=True)
    current = [edgeRef['file'] for edgeRef in edgeRefs]
    for name in os.listdir(triDir):
        if name.startswith('EDGE-') and name not in current and stripExt(name) + '.stl' not in current:
            os.remove(os.path.join(triDir, name))
    for edgeRef in edgeRefs:
        pidShift = shifts.get(edgeRef['geomName'], {})
        path = os.path.join(triDir, edgeRef['file'])
        if os.path.islink(path):
            os.remove(path)
        motionUtilities.writeStlTriangles(path, dict((pid, t + np.array(pidShift.get(pid, (0.0, 0.0, 0.0))))
                                                     for pid, t in edgeRef['baseTris'].items()))
    return current


#----- system dicts

def getBackgroundMesh(bbox, baseCell, nMargin=2):
    #bbox (xmin ymin zmin xmax ymax zmax) grown by nMargin cells and rounded to whole cells
    lo = np.floor(np.array(bbox[:3])/baseCell - nMargin)*baseCell
    hi = np.ceil(np.array(bbox[3:])/baseCell + nMargin)*baseCell
    n = np.rint((hi - lo)/baseCell).astype(int)
    return lo, hi, n


def gridAxes(point, baseCell, level, zExtra=0):
    #axes on which the point lies on a cell face of a background mesh of whole base cells refined up to level
    #(zExtra more levels in z: the end-gap bands split cells in z only, up to 4 times)
    axes = []
    for axis, v in enumerate(point):
        fine = baseCell/2**(level + (zExtra if axis == 2 else 0))
        if abs(v/fine - round(v/fine))*fine < 1e-9*baseCell:
            axes.append(axis)
    return axes


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


def writeSnappyInternal(templateLoc, geomDict, fullCaseSetupDict, motionInfo=None, meshGroups=None):
    print('\n\tWriting snappyHexMeshDict...')
    refinement = fullCaseSetupDict['GLOBAL_REFINEMENT']
    locInMesh = [v for v in refinement['LOC_IN_MESH'] if v != '']
    if len(locInMesh) != 3:
        sys.exit('ERROR! GLOBAL_REFINEMENT LOC_IN_MESH must be 3 coordinates of a point inside the gas!')
    try:
        [float(v) for v in locInMesh]
    except ValueError:
        sys.exit('ERROR! GLOBAL_REFINEMENT LOC_IN_MESH must be 3 numbers, got %s!' % (locInMesh))
    #LOC_IN_MESH on a cell face stops snappy ("not inside the mesh or on a face or edge"): the background box corners
    #are whole base cells, so faces sit at multiples of BASE_CELL/2^level from 0 (x = 0 on an axis through the origin).
    #On the base grid (faces everywhere) it is an error; on a finer level only where the point's cell is refined that
    #far (near a wall), a warning. The moved segment points sit in z-refined end gaps: shifted a little off any face
    baseCell = float(getVal(fullCaseSetupDict, 'GLOBAL_REFINEMENT', 'BASE_CELL'))
    meshGroups = meshGroups or {}
    maxLevel = max([int(geomDict[geom]['refinement']) for geom in geomDict]
                   + [int(gapRef['level']) for gapRef in (motionInfo['gapRefs'] if motionInfo else [])]
                   + [group['level'][1] for groups in meshGroups.values() for group in groups if group['level']]
                   + [group['edgeLevel'] for groups in meshGroups.values() for group in groups
                      if group['edgeLevel'] is not None])
    point = [float(v) for v in locInMesh]
    onGrid, onFine = gridAxes(point, baseCell, 0), gridAxes(point, baseCell, maxLevel)
    if onGrid or onFine:
        axes = onGrid or onFine
        suggest = [v + (0.3*baseCell/2**maxLevel if axis in axes else 0.0) for axis, v in enumerate(point)]
        message = ('GLOBAL_REFINEMENT LOC_IN_MESH %s lies on a cell face %s in %s (faces at multiples of BASE_CELL/2^level '
                   'from 0): move it off, e.g. LOC_IN_MESH = %s (check it is still in the gas)'
                   % (' '.join(locInMesh), 'of the background mesh' if onGrid else 'of refinement level %d' % (maxLevel),
                      ' '.join('xyz'[axis] for axis in axes), ' '.join('%.6g' % v for v in suggest)))
        if onGrid:
            sys.exit('ERROR! ' + message)
        print('\t\tWARNING: %s, if its cell is refined that far' % (message))
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
    #MESH_GROUPS EDGE_LEVEL: feature edges of the group's own STL (outline of its faces + sharp edges) at that level
    for geom in geomDict:
        geomName = stripExt(geom)
        if not geomName.startswith('REF-'):
            for group in meshGroups.get(geomName, []):
                if group.get('edgeFile'):
                    featureArray.append('        { file "%s.eMesh"; scale %s; level %d; }'
                                        % (stripExt(group['edgeFile']), geomDict[geom]['scale'], group['edgeLevel']))
                    extractArray.append('%s { extractionMethod extractFromSurface; includedAngle %s; '
                                        'subsetFeatures { nonManifoldEdges no; openEdges yes; } writeObj no; }'
                                        % (group['edgeFile'], featAngle))
    #walls without layers next to walls with layers: nSurfaceLayers -1 lets their points slide along the wall to make
    #room for the neighbour's layers, 0 holds them (v2412: zeroFixedValue). With 0 side walls a layer stack on the
    #end faces of a box was cut back from the edges (absolute sizes: removed completely), with -1 it was complete.
    #Only in cases with mesh-group layers, so the meshes of other cases do not change
    zeroMode = refinement.get('ZERO_LAYERS', ['slip'])[0]
    if zeroMode not in ('slip', 'fixed'):
        sys.exit('ERROR! GLOBAL_REFINEMENT ZERO_LAYERS must be slip or fixed!')
    groupLayers = any(group['stack'] for groups in meshGroups.values() for group in groups)
    zeroLayers = -1 if zeroMode == 'slip' and groupLayers else 0
    for geom in geomDict:
        geomName = stripExt(geom)
        geomFile = geom.replace('.gz', '')
        level = geomDict[geom]['refinement']
        geometryArray.append('    %s { type triSurfaceMesh; name %s; file "%s"; scale %s; }'
                             % (geomFile, geomName, geomFile, geomDict[geom]['scale']))
        if not geomName.startswith('REF-'):
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
        #per-region levels: MESH_GROUPS LEVEL by exact PID (an exact key wins over the regex keys), then by PID name
        #keyword (regex) from system/snappyRefinementDict
        #   regionLevels=''.join('%s { mode inside; levels ((1E15 %s));}')  #GEOM_NAME {mode inside; levels ((1E15 REF_LEVEL));}\n'
        if not geomName.startswith('REF-'):
            regionLevels = ''.join('%s { level (%d %d); } ' % ((pid,) + group['level'])
                                for group in meshGroups.get(geomName, []) if group['level'] for pid in group['pids'])
            surfaceArray.append('        %s { level (%s %s); patchInfo { type wall; } regions { %s#include "snappyRefinementDict" } }'
                                % (geomName, level, level, regionLevels))
            nLayers = int(geomDict[geom]['layers'])
            layerArray.append('        "%s_.*" { nSurfaceLayers %s; expansionRatio %s;}'
                            % (geomName, nLayers if nLayers > 0 else zeroLayers, geomDict[geom]['expansion']))
    #gap faces: no layers, written after the file-wide entry (last matching regex wins)
    for geom in geomDict:
        geomName = stripExt(geom)
        if not geomName.startswith('REF-'):
            gapPids = [p for p in fullCaseSetupDict[geomName].get('GAP_PIDS', ['']) if p != '']
            for kw in gapPids:
                layerArray.append('        "%s_.*%s.*" { nSurfaceLayers %d; }' % (geomName, kw, zeroLayers))
    #MESH_GROUPS layers by exact patch, last (after the file-wide and gap entries)
    reportLines = []
    for geom in geomDict:
        geomName = stripExt(geom)
        if not geomName.startswith('REF-'):
            for group in meshGroups.get(geomName, []):
                entry = 'nSurfaceLayers %d;' % (zeroLayers) if group.get('autoNone') else group['entry']
                if entry is None:
                    continue
                for pid in group['pids']:
                    patch = '%s_%s' % (geomName, pid)
                    layerArray.append('        "%s" { %s }' % (patch, entry))
                    if group['stack']:
                        reportLines.append('%s %s %d %.6g' % (patch, group['name'], group['layers'], group['stack'][0]))
    #patches with mesh-group layers for meshUtilities/layerReport.py (coverage after snappy, every remesh segment)
    os.makedirs('system', exist_ok=True)
    with open(os.path.join('system', 'meshGroups'), 'w') as f:
        f.write('#patch group layers firstLayerThickness (written by internalSetup, read by layerReport.py)\n')
        f.write(''.join(line + '\n' for line in reportLines))
    #gap refinement: distance band around a per-body STL of the gap faces
    regionArray = []
    for gapRef in (motionInfo['gapRefs'] if motionInfo else []):
        geometryArray.append('    %s { type triSurfaceMesh; name %s; file "%s"; }'
                             % (gapRef['file'], gapRef['name'], gapRef['file']))
        regionArray.append('        %s { mode distance; levels ((%.6g %d)); }'
                           % (gapRef['name'], gapRef['distance'], gapRef['level']))
    # add refinement region defined by geometry
    for geom in geomDict:
        if geom.startswith('REF-'):
            refLevel = geomDict[geom]['refinement']
            regionArray.append('        %s { mode inside; levels ((1E15 %s)); }' % (stripExt(geom),refLevel))
    
    #one dict per meshing angle: the case mesh, then meshes/segNN for remeshing (own LOC_IN_MESH and end-gap bands)
    meshSetups = motionInfo.get('meshSetups') if motionInfo else None
    if not meshSetups:
        meshSetups = [{'loc': locInMesh, 'bands': []}]
    names = segmentNames(motionInfo['segments']) if motionInfo and motionInfo.get('segments') else ['seg01']
    for name, meshSetup in zip(names, meshSetups):
        if name != 'seg01':
            point = [float(v) for v in meshSetup['loc']]
            for axis in gridAxes(point, baseCell, maxLevel, 4):
                point[axis] += 0.3*baseCell/2**(maxLevel + 4)
            meshSetup['loc'] = point
        geometryBands, regionBands = [], []
        for k, (zLo, zHi, splits) in enumerate(meshSetup['bands']):
            #z-only (directional) refinement of every cell in the band, after the other refinement
            geometryBands.append('    GAPZ-%d { type searchableBox; min (-1e3 -1e3 %.9g); max (1e3 1e3 %.9g); }' % (k, zLo, zHi))
            regionBands.append('        GAPZ-%d { mode inside; levels ((1e15 0)); levelIncrement (0 20 (0 0 %d)); }' % (k, splits))
        replaceDict = {'GEOMETRY': '\n'.join(geometryArray + geometryBands),
                       'REFINEMENT_SURFACES': '\n'.join(surfaceArray),
                       'REFINEMENT_REGIONS': '\n'.join(regionArray + regionBands),
                       'LAYERS': '\n'.join(layerArray),
                       'LOC_IN_MESH': ' '.join(locInMesh) if name == 'seg01' else ' '.join('%.9g' % v for v in meshSetup['loc']),
                       'RESOLVE_FEAT_ANGLE': refinement['RESOLVE_FEAT_ANGLE'][0],
                       'N_CELLS_BETWEEN_LEVELS': refinement['N_CELLS_BETWEEN_LEVELS'][0],
                       'FINAL_LAYER_THICKNESS': refinement['FINAL_LAYER_THICKNESS'][0],
                       'MIN_LAYER_THICKNESS': refinement['MIN_LAYER_THICKNESS'][0],
                       'FEATURES': '\n'.join(featureArray),
                       'EXPLICIT_SNAP': 'true' if featureSnap == 'explicit' else 'false',
                       'IMPLICIT_SNAP': 'false' if featureSnap == 'explicit' else 'true'}
        #seg01 is the case mesh, the other meshes get their own dict (runScript copies it over the case's)
        dest = 'system/snappyHexMeshDict' if name == 'seg01' else os.path.join('meshes', name, 'snappyHexMeshDict')
        copyTemplate(templateLoc, 'system/snappyHexMeshDict', dest, replaceDict)
        if name != 'seg01':
            print('\t\t%s: locationInMesh (%s)%s' % (name, replaceDict['LOC_IN_MESH'],
                                                  ', %d end-gap z bands' % (len(regionBands)) if regionBands else ''))
    print('\t\tfeature edges: includedAngle %s, level %s, %s snapping' % (featAngle, featureLevel, featureSnap))
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
    #low: resolved walls; lowReynoldsInternal sets omega by omegaWallFunction (the shared lowReynolds option fixes omega
    #at 1e-8, near zero at the wall instead of 6 nu/(beta1 y^2): kOmegaSST blew up in one omega solve)
    wallFunctionNames = {'low': 'lowReynoldsInternal', 'high': 'highReynolds'}
    if turbModel != 'laminar':
        simControl = fullCaseSetupDict['GLOBAL_SIM_CONTROL']
        initArray += ['k uniform %s;' % (simControl['INIT_K'][0]), 'omega uniform %s;' % (simControl['INIT_OMEGA'][0]),
                      'epsilon uniform 1e-4;', 'nut uniform 0;', 'nuTilda uniform 0;', 'alphat uniform 0;']
    #moving mesh: the points of static walls slide along them, slidingWallVelocity keeps the gas from leaking through
    #(noSlip leaks, movingWallVelocity drags the gas with the points)
    staticMotion = 'slidingMesh' if getMoveSections(fullCaseSetupDict) else 'stationary'
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
            #slip keeps U.n = 0 on static walls whose points slide with the mesh, but those faces sweep a little volume
            #each step: gas crosses them (noSlip walls get slidingWallVelocity instead, which follows the face)
            if static and group['U'] == 'slip' and staticMotion == 'slidingMesh':
                print('\t\tWARNING: %s group %s is slip on a moving mesh: its faces leak gas where their points slide '
                      '(about 0.5-1 %% of the mass per cycle in tests). Use noSlip for a sealed wall'
                      % (geomName, group['name']))
            for pids, motion, suffix in ((static, staticMotion, ''), (moving, 'movingMesh', '_moving')):
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


#random fibres: felt correlation up to the highest tested porosity, single fibres from fiberEpsCyl, log-linear between
fiberEpsFelt = 0.84
fiberEpsCyl = 0.97


def fiberBlend(porosity):
    return min(max((porosity - fiberEpsFelt)/(fiberEpsCyl - fiberEpsFelt), 0.0), 1.0)


def fiberCoeffs(porosity, wireDia, velocity, gas):
    #random fibres (steel wool, felt), two models blended by porosity (the C++ source does the same per cell):
    #  Gedeon-Wood metal felts (NASA CR-198442, eqs 2.5, 2.10, porosity up to 0.84 tested): f = 192/Re + 4.53 Re^-0.067,
    #  Nu = (1 + 1.16 Pe^0.66) eps^2.61 on dh = eps*dw/(1-eps), pore velocity u = U/eps;
    #  single fibres in cross flow (the high-porosity limit, where the felt forms fall off as (1-eps)^2): fibre length
    #  per volume 4(1-eps)/(pi dw^2), Cd = 1 + 10 Re_d^-2/3 (Lamb 8 pi/(Re_d (2.002 - ln Re_d)) below Re_d 1),
    #  Nu_d Churchill-Bernstein; felt up to eps 0.84, single fibres from 0.97, log-linear between (steel wool is
    #  ~0.97-0.99: the fibre model governs, it gives 4-10x the h*a of the felt forms extrapolated there)
    #returns drag coefficient c (force/volume = c*U, N s/m4) and h*a (W/m3/K): blended, felt, single fibres
    eps = porosity
    phi = 1.0 - eps
    dh = eps*wireDia/phi
    up = velocity/eps
    kGas = gas['MU']*gas['CP']/gas['PR']
    re = max(gas['RHO']*up*dh/gas['MU'], 1e-6)
    cFelt = gas['MU']*96.0/(eps*dh**2) + 0.5*gas['RHO']*velocity*4.53*re**-0.067/(eps**2*dh)
    haFelt = (1.0 + 1.16*(re*gas['PR'])**0.66)*eps**2.61*kGas/dh*4.0*phi/wireDia
    reD = max(gas['RHO']*up*wireDia/gas['MU'], 1e-6)
    cd = 1.0 + 10.0*reD**(-2.0/3.0)
    if reD < 1.0:
        cd = max(cd, 8.0*math.pi/(reD*(2.002 - math.log(reD))))
    cCyl = 0.5*gas['RHO']*up*wireDia*cd*4.0*phi/(math.pi*wireDia**2)/eps
    nuD = 0.3 + 0.62*reD**0.5*gas['PR']**(1.0/3.0)/(1.0 + (0.4/gas['PR'])**(2.0/3.0))**0.25
    haCyl = nuD*kGas/wireDia*4.0*phi/wireDia
    w = fiberBlend(porosity)
    return cFelt**(1.0 - w)*cCyl**w, haFelt**(1.0 - w)*haCyl**w, cFelt, cCyl, haFelt, haCyl


def fiberConductivity(solidFrac, wireDia, solidK, T, radiation=True):
    #effective conductivity of the fibre bed (solid side): conduction along random fibres k_s*phi/3, radiation
    #(Rosseland) 16 sigma T^3/(3 beta), extinction beta = phi/dw (mean projected fibre area per volume = surface/4)
    kc = solidK*solidFrac/3.0
    kr = 16.0*5.670374e-8*T**3*wireDia/(3.0*solidFrac) if radiation and solidFrac > 0 else 0.0
    return kc, kr


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
            #COEFFS local (default): the Re-dependent Forchheimer term and h*a per cell and step from the velocity
            #relative to the matrix (oscillating flow: h*a varies ~7x over a cycle), ref: constants at REF_VELOCITY
            coeffs = setup.get('COEFFS', ['local'])[0].lower()
            if coeffs not in ('local', 'ref'):
                sys.exit('ERROR! [%s] COEFFS must be local or ref!' % (name))
            zone['wire'] = {'eps': porosity, 'dw': wireDia, 'reRef': re, 'local': coeffs == 'local'}
            print('\t\t%s: wire mesh porosity %g, wire %g m, Re %.3g at %g m/s: d %.4g 1/m2, f %.4g 1/m, h*a %.4g W/m3/K'
                  % (name, porosity, wireDia, re, refVelocity, d, f, ha))
            if coeffs == 'local':
                print('\t\t\tf and h*a per cell from the local velocity relative to the matrix (COEFFS local):')
                for u in (0.1, 1.0, 5.0, 20.0):
                    dU, fU, haU, reU = wireMeshCoeffs(porosity, wireDia, u, gas)
                    print('\t\t\t\t%5g m/s: Re %8.3g  f %.4g 1/m  h*a %.4g W/m3/K  dp/L %.4g Pa/m'
                          % (u, reU, fU, haU, gas['MU']*dU*u + 0.5*gas['RHO']*fU*u*u))
            else:
                print('\t\t\tf and h*a fixed at REF_VELOCITY (COEFFS ref)')
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
        elif resistance == 'randomFiber':
            #random fibres (steel wool, felt): drag and heat transfer per cell from its porosity and the gas velocity
            #relative to the fibres, at run time (porosity changes when the zone deforms with the mesh)
            wireDia, solidRho = number('WIRE_DIA'), number('SOLID_RHO')
            deforming = setup.get('MOVES_WITH', ['none'])[0] == 'mesh'
            mass = [v for v in setup.get('SOLID_MASS', ['']) if v != '']
            if deforming and not mass:
                sys.exit('ERROR! [%s] MOVES_WITH mesh needs SOLID_MASS (kg of fibre in the zone)' % (name))
            if wireDia <= 0 or solidRho <= 0:
                sys.exit('ERROR! [%s] randomFiber needs WIRE_DIA > 0 (fibre diameter) and SOLID_RHO > 0' % (name))
            zone['fiber'] = {'dw': wireDia, 'rhos': solidRho, 'mass': number('SOLID_MASS') if mass else None,
                             'porosity': None if mass else number('POROSITY')}
            #heat through the fibres (conduction k_s*phi/3, radiation 16 sigma T^3 dw/(3 phi)) and from the walls it
            #touches (WALL_CONTACT_PIDS keywords in the patch names, WALL_CONTACT factor 0-1 on that conductance)
            zone['fiber']['ks'] = float(setup.get('SOLID_K', ['50'])[0] or 0)
            zone['fiber']['radiation'] = setup.get('RADIATION', ['yes'])[0].lower() == 'yes'
            zone['fiber']['contactPids'] = [k for k in setup.get('WALL_CONTACT_PIDS', ['']) if k != '']
            zone['fiber']['contact'] = float(setup.get('WALL_CONTACT', ['1'])[0] or 0)
            if zone['fiber']['ks'] < 0 or not 0 <= zone['fiber']['contact'] <= 1:
                sys.exit('ERROR! [%s] SOLID_K must be >= 0 and WALL_CONTACT between 0 and 1' % (name))
            if zone['fiber']['porosity'] is not None and not 0 < zone['fiber']['porosity'] < 1:
                sys.exit('ERROR! [%s] randomFiber needs 0 < POROSITY < 1 (or SOLID_MASS)' % (name))
            zone['d'], zone['f'] = [0.0]*3, [0.0]*3
            print('\t\t%s: random fibres %g m, %s' % (name, wireDia, 'solid mass %g kg (porosity from the zone volume '
                                                     'at run time, follows its deformation)' % (zone['fiber']['mass'])
                                                     if mass else 'porosity %g' % (zone['fiber']['porosity'])))
            refVelocity = float(setup.get('REF_VELOCITY', ['0.5'])[0] or 0.5)
            print('\t\t\tat %g m/s (charge state): porosity  drag c (N s/m4)  h*a (W/m3/K)   [felt / single fibres]'
                  % (refVelocity))
            for eps in ([zone['fiber']['porosity']] if not mass else [0.9, 0.95, 0.98, 0.99]):
                c, ha, cF, cC, hF, hC = fiberCoeffs(eps, wireDia, refVelocity, gas)
                kc, kr = fiberConductivity(1.0 - eps, wireDia, zone['fiber']['ks'], 600.0, zone['fiber']['radiation'])
                print('\t\t\t%30.3f  %10.3g  %12.3g   [%.3g / %.3g, %.3g / %.3g]   k %.3g W/m/K at 600 K'
                      ' (fibres %.3g, radiation %.3g)' % (eps, c, ha, cF, cC, hF, hC, kc + kr, kc, kr))
            fiber = zone['fiber']
            print('\t\t\tfibre conductivity %g W/m/K, radiation %s, wall contact %s'
                  % (fiber['ks'], 'on' if fiber['radiation'] else 'off',
                     '%s (factor %g)' % (' '.join(fiber['contactPids']), fiber['contact']) if fiber['contactPids'] else 'none'))
        else:
            sys.exit('ERROR! [%s] RESISTANCE must be coeffs, wireMesh, fit or randomFiber!' % (name))
        if setup.get('MOVES_WITH', ['none'])[0] == 'mesh' and 'fiber' not in zone:
            sys.exit('ERROR! [%s] MOVES_WITH mesh needs RESISTANCE randomFiber (porosity follows the deformation)' % (name))

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
            if 'fiber' in zone:
                if haSetting != 'auto':
                    sys.exit('ERROR! [%s] randomFiber works out h*a per cell, HA must be auto' % (name))
                zone['ha'] = 0.0
            elif haSetting == 'auto':
                if 'wireHA' not in zone:
                    sys.exit('ERROR! [%s] HA auto needs RESISTANCE wireMesh, else give HA (W/m3/K)' % (name))
                zone['ha'] = zone['wireHA']
                #per cell at run time with COEFFS local (zone['ha'] only for the printout)
                zone['haLocal'] = zone['wire']['local']
            else:
                zone['ha'] = number('HA')
            #SOLID_T: one value = uniform, two = linear from the -z to the +z end of the zone (start values only)
            nSolidT = len([v for v in setup.get('SOLID_T', ['']) if v != ''])
            if 'fiber' in zone and nSolidT == 2:
                sys.exit('ERROR! [%s] randomFiber takes one SOLID_T (uniform start)' % (name))
            solidT = number('SOLID_T', 2 if nSolidT == 2 else 1)
            zone['Ts'] = solidT if nSolidT == 2 else [solidT, solidT]
            zone['Cs'] = 0.0
            if heat == 'solved' and 'fiber' in zone:
                #per cell at run time: solid fraction * rhos * cs
                zone['fiber']['cs'] = number('SOLID_CP')
            elif heat == 'solved':
                zone['Cs'] = (1.0 - number('POROSITY'))*number('SOLID_RHO')*number('SOLID_CP')
            matrixT = ('%g K' % (zone['Ts'][0]) if zone['Ts'][0] == zone['Ts'][1]
                       else '%g K (-z end) to %g K (+z end), linear' % tuple(zone['Ts']))
            if 'fiber' in zone:
                print('\t\t%s: heat %s, h*a per cell (Gedeon-Wood random fibres), matrix T %s%s'
                      % (name, heat, matrixT, ', cs %g J/kg/K' % (zone['fiber']['cs']) if heat == 'solved' else ' (fixed)'))
            else:
                print('\t\t%s: heat %s, h*a %s, matrix T %s%s'
                      % (name, heat, 'per cell (Gedeon-Wood screens, local velocity)' if zone.get('haLocal')
                         else '%.4g W/m3/K' % (zone['ha']), matrixT,
                         ', (1-eps)*rhos*cs %.4g J/m3/K' % (zone['Cs']) if heat == 'solved' else ' (fixed)'))
        porous.append(zone)
    return porous


def porousHeatSource(zone, cv, props=None):
    #scalarCodedSource on e: gas-matrix heat exchange, gas and matrix implicit together (backward Euler), so the heat
    #the gas gets each step is exactly the heat the matrix gives (Cs 0: fixed matrix temperature):
    #  Cs*(Ts1 - Ts0)/dt = h*a*(T1 - Ts1)  ->  gas source h*a*(Ts1 - T1) = hEff*(Ts0 - T1), hEff = h*a*Cs/(Cs + h*a*dt)
    #addSup puts hEff*(Ts0 - T1) into the e equation (implicit in e, cv constant); correct (after the e solve) sets
    #Qs = hEff*(Ts0 - T1) (W/m3 into the gas) and Ts1 from the solved e; Ts0 moves on once per time step
    #Ts_<zone> is written (restarts), TsOld/Qs/eLin/TLin_<zone> are working fields
    #wire mesh with COEFFS local and HA auto: h*a per cell and step from the velocity relative to the matrix
    #(Gedeon-Wood screens, Nu = (1 + 0.99 (Re Pr)^0.66) eps^1.79), hEff_<zone> keeps it for the correct hook
    name = zone['name'].replace('-', '_')
    local = bool(zone.get('haLocal'))
    base = ['            const fvMesh& m = mesh_;',
            '            const scalar Cv = %.9g;' % (cv),
            '            const scalar Cs = %.9g;' % (zone['Cs']),
            '            const scalar dt = m.time().deltaTValue();',
            '            const word zoneName("%s");' % (zone['name']),
            '            const labelList& zoneCells = this->cells();']
    fields = ['Ts', 'TsOld', 'Qs', 'eLin', 'TLin'] + (['hEff'] if local else [])
    lookups = '\n'.join(['            volScalarField& %s = const_cast<volScalarField&>(m.lookupObject<volScalarField>'
                         '("%s_" + zoneName));' % (f, f) for f in fields])
    if local:
        gasLines, gasCell = gasCellCode(props, False, 0.0)
        wireConsts, wireCell = wireCellCode(zone)
        constants = '\n'.join(base)
        #ReRef is only used by the drag scaling
        addConstants = '\n'.join(base + gasLines + [w for w in wireConsts if 'ReRef' not in w] + [
            '            const scalar dwW = %.9g;' % (zone['wire']['dw']),
            '            const scalar Pr = %.9g;' % (props['PR']),
            '            const scalar cpGas = %.9g;' % (props['CP'])])
        hEffInit = '''
                (new volScalarField(IOobject("hEff_" + zoneName, now, m), m,
                    dimensionedScalar(dimPower/dimVolume/dimTemperature, Zero)))->store();'''
        kLine = ''
        cellCode = '\n'.join([
            '                const vector Ur(Ugas[c] - (Umat ? (*Umat)[c] : vector(Zero)));'] + gasCell + [wireCell[0],
            '                const scalar kC = muC*cpGas/Pr;',
            '                const scalar haC = (1 + 0.99*pow(ReC*Pr, 0.66))*pow(epsW, 1.79)*kC/dhW*4*(1 - epsW)/dwW;',
            '                hEff[c] = (Cs > 0 ? haC*Cs/(Cs + haC*dt) : haC);',
            '                const scalar k = hEff[c]/Cv;']) + '\n'
        velocity = '''
            const volVectorField& Ugas = m.lookupObject<volVectorField>("U");
            const volVectorField* Umat = m.findObject<volVectorField>("Umat_" + zoneName);'''
        hE = 'hEff[c]'
    else:
        constants = '\n'.join(base + ['            const scalar ha = %.9g;' % (zone['ha']),
                                      '            const scalar hEff = (Cs > 0 ? ha*Cs/(Cs + ha*dt) : ha);'])
        addConstants = constants
        hEffInit = ''
        kLine = '\n            const scalar k = hEff/Cv;'
        cellCode = ''
        velocity = ''
        hE = 'hEff'
    addSup = '''        codeAddSupRho
        #{
%s
            const scalar TinitLow = %.9g;
            const scalar TinitHigh = %.9g;
            const volScalarField& T = m.lookupObject<volScalarField>("T");
            if (!m.foundObject<volScalarField>("Ts_" + zoneName))
            {
                //first call is already in the first new time step: read the matrix temperature from the start
                //(restart) time, not the current one, else every restart begins again from the initial value
                const word startName(m.time().timeName(m.time().startTime().value()));
                IOobject TsIO("Ts_" + zoneName, startName, m, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE);
                const bool TsFound = TsIO.typeHeaderOk<volScalarField>(true);
                volScalarField* TsPtr = new volScalarField(TsIO, m, dimensionedScalar(dimTemperature, TinitLow));
                TsPtr->store();
                //no saved matrix temperature: linear from the -z to the +z end of the zone (cell centres)
                if (!TsFound && TinitHigh != TinitLow)
                {
                    const volVectorField& C = m.C();
                    scalar zLow = GREAT, zHigh = -GREAT;
                    forAll(zoneCells, i)
                    {
                        zLow = min(zLow, C[zoneCells[i]].z());
                        zHigh = max(zHigh, C[zoneCells[i]].z());
                    }
                    reduce(zLow, minOp<scalar>());
                    reduce(zHigh, maxOp<scalar>());
                    forAll(zoneCells, i)
                    {
                        const scalar s = (C[zoneCells[i]].z() - zLow)/max(zHigh - zLow, SMALL);
                        (*TsPtr)[zoneCells[i]] = TinitLow + s*(TinitHigh - TinitLow);
                    }
                    Info<< "    Ts_" << zoneName << ": linear start " << TinitLow << " K at z " << zLow << " to "
                        << TinitHigh << " K at z " << zHigh << endl;
                }
                //working fields (not written)
                const word now(m.time().timeName());
                (new volScalarField(IOobject("TsOld_" + zoneName, now, m), *TsPtr))->store();
                (new volScalarField(IOobject("Qs_" + zoneName, now, m), m,
                    dimensionedScalar(dimPower/dimVolume, Zero)))->store();
                (new volScalarField(IOobject("eLin_" + zoneName, now, m), m,
                    dimensionedScalar(dimEnergy/dimMass, Zero)))->store();
                (new volScalarField(IOobject("TLin_" + zoneName, now, m), m,
                    dimensionedScalar(dimTemperature, Zero)))->store();%s
            }
%s%s
            //new time step: the matrix temperature at the end of the last step is the old level
            static label lastIndex = -1;
            if (m.time().timeIndex() != lastIndex)
            {
                lastIndex = m.time().timeIndex();
                forAll(zoneCells, i)
                {
                    TsOld[zoneCells[i]] = Ts[zoneCells[i]];
                }
            }
            const scalarField& V = m.V();
            const scalarField& e = eqn.psi();
            scalarField& diag = eqn.diag();
            scalarField& source = eqn.source();%s
            forAll(zoneCells, i)
            {
                const label c = zoneCells[i];
%s                eLin[c] = e[c];
                TLin[c] = T[c];
                diag[c] -= k*V[c];
                source[c] -= (%s*(TsOld[c] - T[c]) + k*e[c])*V[c];
            }
        #};''' % (addConstants, zone['Ts'][0], zone['Ts'][1], hEffInit, lookups, velocity, kLine, cellCode, hE)
    correct = '''        codeCorrect
        #{
%s
            if (!m.foundObject<volScalarField>("Ts_" + zoneName))
            {
                return;
            }
%s
            //gas temperature after the e solve (cv constant), the heat it got, and the matrix temperature that gave it
            forAll(zoneCells, i)
            {
                const label c = zoneCells[i];
                const scalar T1 = TLin[c] + (fld[c] - eLin[c])/Cv;
                Qs[c] = %s*(TsOld[c] - T1);
                if (Cs > 0)
                {
                    Ts[c] = TsOld[c] - Qs[c]*dt/Cs;
                }
            }
        #};''' % (constants, lookups, hE)
    lines = ['%sHeat' % (name), '{', '    type            scalarCodedSource;', '    name            %sHeat;' % (name),
             '    scalarCodedSourceCoeffs', '    {', '        selectionMode   cellZone;',
             '        cellZone        %s;' % (zone['name']), '        fields          (e);', '',
             '        codeInclude', '        #{', '            #include "fvCFD.H"', '        #};', '',
             correct, '', '        codeAddSup #{ #};', '        codeConstrain #{ #};', '', addSup,
             '    }', '}']
    return '\n'.join(lines)


def getDissipation(fullCaseSetupDict):
    #DISSIPATION yes: viscous dissipation (ESI viscousDissipation) and porous drag work added to the gas energy,
    #compressible only. Default no: the rhoPimpleFoam energy equation is in total-energy form (K terms, no tau.U and no
    #momentum sources), so the kinetic energy that viscosity and drag destroy already ends up in e (globally; the
    #viscous part is placed by -U.div(tau) instead of tau:grad(U)). Switching it on counts that heat twice
    setting = fullCaseSetupDict['GLOBAL_SIM_CONTROL'].get('DISSIPATION', ['no'])[0].lower()
    if setting not in ('yes', 'no'):
        sys.exit('ERROR! GLOBAL_SIM_CONTROL DISSIPATION must be yes or no!')
    return setting == 'yes' and getFlow(fullCaseSetupDict) == 'compressible'


def porousTensors(zone):
    #Darcy-Forchheimer d, f in global axes: sum of d_i e_i e_i over the zone axes e1, e2, e3 = e1 x e2
    e1 = np.array(zone['e1'], dtype=float)
    e1 /= np.linalg.norm(e1)
    e3 = np.cross(e1, np.array(zone['e2'], dtype=float))
    e3 /= np.linalg.norm(e3)
    e2 = np.cross(e3, e1)
    axes = [e1, e2, e3]
    tensors = []
    for coeffs in (zone['d'], zone['f']):
        t = sum(c*np.outer(e, e) for c, e in zip(coeffs, axes))
        tensors.append('tensor(%s)' % (', '.join('%.9g' % v for v in t.flatten())))
    return tensors


def gasCellCode(props, incompressible, nu):
    #C++ lines: gas density rhoC and viscosity muC of cell c (incompressible: 1 and the kinematic viscosity)
    if incompressible:
        return [], ['                const scalar rhoC = 1;', '                const scalar muC = %.9g;' % (nu)]
    lines = ['            const volScalarField& rhoF = m.lookupObject<volScalarField>("rho");']
    cell = ['                const scalar rhoC = rhoF[c];']
    if props['TRANSPORT'] == 'sutherland':
        lines.append('            const volScalarField& TF = m.lookupObject<volScalarField>("T");')
        cell.append('                const scalar muC = %.9g*sqrt(TF[c])/(1 + %.9g/TF[c]);' % (props['AS'], props['TS']))
    else:
        cell.append('                const scalar muC = %.9g;' % (props['MU']))
    return lines, cell


def wireCellCode(zone):
    #C++ lines for a wire-mesh zone with COEFFS local: Reynolds number of cell c from |Ur| (pore velocity |Ur|/eps,
    #hydraulic diameter dh), Gedeon-Wood Forchheimer scaling (Re/ReRef)^-0.103 of the REF_VELOCITY tensor Ft.
    #Needs rhoC, muC, Ur in scope
    wire = zone['wire']
    dh = wire['eps']*wire['dw']/(1.0 - wire['eps'])
    consts = ['            const scalar epsW = %.9g;' % (wire['eps']), '            const scalar dhW = %.9g;' % (dh),
              '            const scalar ReRef = %.9g;' % (wire['reRef'])]
    cell = ['                const scalar ReC = max(rhoC*mag(Ur)/epsW*dhW/muC, scalar(1e-6));',
            '                const scalar fScale = pow(ReC/ReRef, -0.103);']
    return consts, cell


def porousDragCoeff(zone, props, incompressible, nu):
    #C++ lines: tensors Dt, Ft and the per-cell drag tensor C (force on the gas = -C & (U - Umat)) from rho, mu, |Ur|
    #wire mesh with COEFFS local: Ft (taken at REF_VELOCITY) scaled per cell with the local Reynolds number
    dT, fT = porousTensors(zone)
    lines = ['            const tensor Dt(%s);' % (dT), '            const tensor Ft(%s);' % (fT)]
    gasLines, cell = gasCellCode(props, incompressible, nu)
    lines += gasLines
    if zone.get('wire', {}).get('local'):
        wireConsts, wireCell = wireCellCode(zone)
        lines += wireConsts
        cell += wireCell
        cell.append('                const tensor C(muC*Dt + (0.5*rhoC*mag(Ur)*fScale)*Ft);')
    else:
        cell.append('                const tensor C(muC*Dt + (0.5*rhoC*mag(Ur))*Ft);')
    return '\n'.join(lines), '\n'.join(cell)


def porousDragSource(zone, props, incompressible, nu, rhoInf, moving=True):
    #vectorCodedSource on U for a porous zone moving with a body (MOVES_WITH): Darcy-Forchheimer drag on the velocity
    #relative to the matrix (explicitPorositySource uses the absolute U), implicit part on the diagonal as ESI does.
    #Umat_<zone> = mean motion of the cell points over the step (the zone moves rigidly with the body),
    #porDrag_<zone> = force of the gas on the matrix (N/m3), added to the body force by the crank drive and logged
    #static zones with per-cell coefficients (wire mesh, COEFFS local) use it too, with Umat_<zone> kept at zero
    name = zone['name'].replace('-', '_')
    constants, cell = porousDragCoeff(zone, props, incompressible, nu)
    head = '\n'.join([
        '            const fvMesh& m = mesh_;',
        '            const word zoneName("%s");' % (zone['name']),
        '            const labelList& zoneCells = this->cells();',
        constants])
    lookups = '\n'.join(['            volVectorField& %s = const_cast<volVectorField&>(m.lookupObject<volVectorField>'
                         '("%s_" + zoneName));' % (f, f) for f in ('Umat', 'porDrag')])
    addSup = '''        #{
%s
            if (!m.foundObject<volVectorField>("Umat_" + zoneName))
            {
                const word now(m.time().timeName());
                (new volVectorField(IOobject("Umat_" + zoneName, now, m), m, dimensionedVector(dimVelocity, Zero)))->store();
                (new volVectorField(IOobject("porDrag_" + zoneName, now, m), m,
                    dimensionedVector(dimForce/dimVolume, Zero)))->store();
            }
%s
            //new time step: matrix velocity from the motion of the cell points (rigid with the body), static zone: zero
            static label lastIndex = -1;
            if (%s && m.time().timeIndex() != lastIndex)
            {
                lastIndex = m.time().timeIndex();
                const scalar dt = m.time().deltaTValue();
                const pointField& pts = m.points();
                const pointField& pts0 = m.moving() ? m.oldPoints() : m.points();
                const labelListList& cellPoints = m.cellPoints();
                forAll(zoneCells, i)
                {
                    const label c = zoneCells[i];
                    vector d(Zero);
                    for (const label pointi : cellPoints[c])
                    {
                        d += pts[pointi] - pts0[pointi];
                    }
                    Umat[c] = d/(cellPoints[c].size()*dt);
                }
            }
            //drag -C & (U - Umat): -iso*U implicit, the rest explicit
            const vectorField& U = eqn.psi();
            const scalarField& V = m.V();
            scalarField& diag = eqn.diag();
            vectorField& source = eqn.source();
            forAll(zoneCells, i)
            {
                const label c = zoneCells[i];
                const vector Ur(U[c] - Umat[c]);
%s
                const scalar iso = tr(C)/3;
                diag[c] -= iso*V[c];
                source[c] -= ((C & Umat[c]) - ((C - iso*I) & U[c]))*V[c];
            }
        #};''' % (head, lookups.split('\n')[0], 'true' if moving else 'false', cell)
    correct = '''        codeCorrect
        #{
%s
            if (!m.foundObject<volVectorField>("Umat_" + zoneName))
            {
                return;
            }
%s
            //force of the gas on the matrix with the solved U (N/m3)
            forAll(zoneCells, i)
            {
                const label c = zoneCells[i];
                const vector Ur(fld[c] - Umat[c]);
%s
                porDrag[c] = %s(C & Ur);
            }
        #};''' % (head, lookups, cell, '%.9g*' % (rhoInf) if incompressible else '')
    hooks = ['        codeAddSup%s' % ('' if incompressible else 'Rho'), addSup, '',
             '        codeAddSup%s #{ #};' % ('Rho' if incompressible else '')]
    lines = ['%sDrag' % (name), '{', '    type            vectorCodedSource;', '    name            %sDrag;' % (name),
             '    vectorCodedSourceCoeffs', '    {', '        selectionMode   cellZone;',
             '        cellZone        %s;' % (zone['name']), '        fields          (U);', '',
             '        codeInclude', '        #{', '            #include "fvCFD.H"', '        #};', '',
             correct, '', '        codeConstrain #{ #};', ''] + hooks + ['    }', '}']
    return '\n'.join(lines)


def fiberCellCode(zone, props, incompressible, nu):
    #C++ for random-fibre zones: per cell drag coefficient cDrag (force on the gas = -cDrag*(U - Umat)) and h*a (haC)
    #from the cell's solid fraction and its velocity relative to the fibres Ur, as fiberCoeffs (felt / single fibres,
    #blended by porosity). Needs solidFrac, c, Ur in scope
    consts = ['            const scalar dw = %.9g;' % (zone['fiber']['dw']),
              '            const scalar Pr = %.9g;' % (props['PR']),
              '            const scalar cpGas = %.9g;' % (props['CP']),
              '            const scalar pi = constant::mathematical::pi;']
    if incompressible:
        cell = ['                const scalar rhoC = 1;', '                const scalar muC = %.9g;' % (nu)]
    else:
        consts.append('            const volScalarField& rhoF = m.lookupObject<volScalarField>("rho");')
        cell = ['                const scalar rhoC = rhoF[c];']
        if props['TRANSPORT'] == 'sutherland':
            consts.append('            const volScalarField& TF = m.lookupObject<volScalarField>("T");')
            cell.append('                const scalar muC = %.9g*sqrt(TF[c])/(1 + %.9g/TF[c]);' % (props['AS'], props['TS']))
        else:
            cell.append('                const scalar muC = %.9g;' % (props['MU']))
    cell += ['                const scalar phiC = min(max(solidFrac[c], scalar(0)), scalar(0.5));',
             '                const scalar eps = 1 - phiC;',
             '                scalar cDrag = 0;',
             '                scalar haC = 0;',
             '                if (phiC > 1e-9)',
             '                {',
             '                    const scalar up = mag(Ur)/eps;',
             '                    const scalar dh = eps*dw/phiC;',
             '                    const scalar kC = muC*cpGas/Pr;',
             '                    //Gedeon-Wood metal felts: f = 192/Re + 4.53 Re^-0.067, Nu = (1 + 1.16 Pe^0.66) eps^2.61',
             '                    const scalar Re = max(rhoC*up*dh/muC, scalar(1e-6));',
             '                    const scalar cFelt = muC*96/(eps*sqr(dh)) + 0.5*rhoC*mag(Ur)*4.53*pow(Re, -0.067)/(sqr(eps)*dh);',
             '                    const scalar haFelt = (1 + 1.16*pow(Re*Pr, 0.66))*pow(eps, 2.61)*kC/dh*4*phiC/dw;',
             '                    //single fibres in cross flow: Cd White (Lamb below Re_d 1), Nu_d Churchill-Bernstein',
             '                    const scalar ReD = max(rhoC*up*dw/muC, scalar(1e-6));',
             '                    scalar Cd = 1 + 10*pow(ReD, -2.0/3.0);',
             '                    if (ReD < 1)',
             '                    {',
             '                        Cd = max(Cd, 8*pi/(ReD*(2.002 - Foam::log(ReD))));',
             '                    }',
             '                    const scalar cCyl = 0.5*rhoC*up*dw*Cd*4*phiC/(pi*sqr(dw))/eps;',
             '                    const scalar NuD = 0.3 + 0.62*sqrt(ReD)*cbrt(Pr)/pow(1 + pow(0.4/Pr, 2.0/3.0), 0.25);',
             '                    const scalar haCyl = NuD*kC/dw*4*phiC/dw;',
             '                    //felt up to eps %g (tested), single fibres from %g, log-linear between' % (fiberEpsFelt, fiberEpsCyl),
             '                    const scalar w = min(max((eps - %.9g)/%.9g, scalar(0)), scalar(1));' % (fiberEpsFelt, fiberEpsCyl - fiberEpsFelt),
             '                    cDrag = pow(cFelt, 1 - w)*pow(cCyl, w);',
             '                    haC = pow(haFelt, 1 - w)*pow(haCyl, w);',
             '                }']
    return '\n'.join(consts), '\n'.join(cell)


def porousFiberDragSource(zone, props, incompressible, nu, rhoInf):
    #vectorCodedSource on U for a random-fibre zone: drag -cDrag*(U - Umat) per cell, implicit in U.
    #solidFrac_<zone> (written, restarts and remeshing maps): fibre volume fraction per cell, follows the cell volume
    #(each cell keeps its fibres when the zone deforms with the mesh: solidFrac *= V0/V every step); with SOLID_MASS it is
    #scaled at every (re)start so the zone holds exactly that much fibre. Umat_<zone> from the motion of the cell points
    #(zero for a static zone), porDrag_<zone> = force of the gas on the fibres (N/m3)
    name = zone['name'].replace('-', '_')
    fiber = zone['fiber']
    consts, cell = fiberCellCode(zone, props, incompressible, nu)
    solidVol = fiber['mass']/fiber['rhos'] if fiber['mass'] is not None else -1.0
    phiInit = 1.0 if fiber['mass'] is not None else 1.0 - fiber['porosity']
    head = '\n'.join([
        '            const fvMesh& m = mesh_;',
        '            const word zoneName("%s");' % (zone['name']),
        '            const labelList& zoneCells = this->cells();',
        consts])
    lookups = '\n'.join(['            volVectorField& %s = const_cast<volVectorField&>(m.lookupObject<volVectorField>'
                         '("%s_" + zoneName));' % (f, f) for f in ('Umat', 'porDrag')] +
                        ['            volScalarField& solidFrac = const_cast<volScalarField&>(m.lookupObject<volScalarField>'
                         '("solidFrac_" + zoneName));'])
    addSup = '''        #{
%s
            const scalarField& V = m.V();
            static label lastIndex = -1;
            if (!m.foundObject<volScalarField>("solidFrac_" + zoneName))
            {
                const word now(m.time().timeName());
                (new volVectorField(IOobject("Umat_" + zoneName, now, m), m, dimensionedVector(dimVelocity, Zero)))->store();
                (new volVectorField(IOobject("porDrag_" + zoneName, now, m), m,
                    dimensionedVector(dimForce/dimVolume, Zero)))->store();
                //fibre fraction from the start (restart) time, else uniform
                const word startName(m.time().timeName(m.time().startTime().value()));
                IOobject fIO("solidFrac_" + zoneName, startName, m, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE);
                const bool found = fIO.typeHeaderOk<volScalarField>(true);
                volScalarField* f = new volScalarField(fIO, m, dimensionedScalar(dimless, Zero));
                f->store();
                if (!found)
                {
                    forAll(zoneCells, i)
                    {
                        (*f)[zoneCells[i]] = %.9g;
                    }
                }
                //with SOLID_MASS: exactly that much fibre in the zone (also after a remeshing map)
                const scalar solidVol = %.9g;
                scalar vZone = 0, vSolid = 0;
                forAll(zoneCells, i)
                {
                    vZone += V[zoneCells[i]];
                    vSolid += (*f)[zoneCells[i]]*V[zoneCells[i]];
                }
                reduce(vZone, sumOp<scalar>());
                reduce(vSolid, sumOp<scalar>());
                if (solidVol > 0 && vSolid > 0)
                {
                    forAll(zoneCells, i)
                    {
                        (*f)[zoneCells[i]] *= solidVol/vSolid;
                    }
                }
                Info<< "    solidFrac_" << zoneName << ": " << (found ? "read" : "uniform start") << ", zone volume " << vZone
                    << " m3, fibre volume " << vSolid << (solidVol > 0 ? " -> " : "") << (solidVol > 0 ? solidVol : vSolid)
                    << " m3, mean porosity " << 1 - (solidVol > 0 ? solidVol : vSolid)/max(vZone, VSMALL) << endl;
                lastIndex = m.time().timeIndex();
            }
%s
            //new time step: fibres stay in their cells (fraction follows the cell volume), matrix velocity from the
            //motion of the cell points
            if (m.time().timeIndex() != lastIndex)
            {
                lastIndex = m.time().timeIndex();
                if (m.moving())
                {
                    const scalarField& V0 = m.V0();
                    forAll(zoneCells, i)
                    {
                        solidFrac[zoneCells[i]] *= V0[zoneCells[i]]/V[zoneCells[i]];
                    }
                }
            }
            {
                const scalar dt = m.time().deltaTValue();
                const pointField& pts = m.points();
                const pointField& pts0 = m.moving() ? m.oldPoints() : m.points();
                const labelListList& cellPoints = m.cellPoints();
                forAll(zoneCells, i)
                {
                    const label c = zoneCells[i];
                    vector d(Zero);
                    for (const label pointi : cellPoints[c])
                    {
                        d += pts[pointi] - pts0[pointi];
                    }
                    Umat[c] = d/(cellPoints[c].size()*dt);
                }
            }
            //drag -cDrag*(U - Umat): implicit in U
            const vectorField& U = eqn.psi();
            scalarField& diag = eqn.diag();
            vectorField& source = eqn.source();
            forAll(zoneCells, i)
            {
                const label c = zoneCells[i];
                const vector Ur(U[c] - Umat[c]);
%s
                diag[c] -= cDrag*V[c];
                source[c] -= cDrag*Umat[c]*V[c];
            }
        #};''' % (head, phiInit, solidVol, lookups, cell)
    correct = '''        codeCorrect
        #{
%s
            if (!m.foundObject<volScalarField>("solidFrac_" + zoneName))
            {
                return;
            }
%s
            //force of the gas on the fibres with the solved U (N/m3)
            forAll(zoneCells, i)
            {
                const label c = zoneCells[i];
                const vector Ur(fld[c] - Umat[c]);
%s
                porDrag[c] = %scDrag*Ur;
            }
        #};''' % (head, lookups, cell, '%.9g*' % (rhoInf) if incompressible else '')
    hooks = ['        codeAddSup%s' % ('' if incompressible else 'Rho'), addSup, '',
             '        codeAddSup%s #{ #};' % ('Rho' if incompressible else '')]
    lines = ['%sDrag' % (name), '{', '    type            vectorCodedSource;', '    name            %sDrag;' % (name),
             '    vectorCodedSourceCoeffs', '    {', '        selectionMode   cellZone;',
             '        cellZone        %s;' % (zone['name']), '        fields          (U);', '',
             '        codeInclude', '        #{', '            #include "fvCFD.H"', '        #};', '',
             correct, '', '        codeConstrain #{ #};', ''] + hooks + ['    }', '}']
    return '\n'.join(lines)


def fiberConductionCode(zone):
    #C++ (codeCorrect of the fibre heat source, after the gas-fibre exchange): heat through the fibre bed, conduction
    #along the fibres k_s*phi/3 plus radiation 16 sigma Ts^3 dw/(3 phi), between neighbouring zone cells (processor faces
    #included, conservative) and from the walls whose patch names contain a WALL_CONTACT_PIDS keyword (wall T, contact
    #factor on the conductance over the cell-centre distance). Explicit, sub-cycled to dt*sum(G)/C <= 0.4.
    #Qwall_<zone> = heat from the walls (W/m3, an extra heat input for energy balances). Ts is reset from TsOld before
    #this on every call (outer correctors), so the step is not applied twice
    fiber = zone['fiber']
    if zone['heat'] != 'solved' or (fiber['ks'] <= 0 and not fiber['radiation'] and not fiber['contactPids']):
        return ''
    keys = '\n'.join('                contactKeys.append(word("%s"));' % (k.lower()) for k in fiber['contactPids'])
    return '''            //heat through the fibre bed and from the walls it touches
            {
                volScalarField& kEff = const_cast<volScalarField&>(m.lookupObject<volScalarField>("kEff_" + zoneName));
                volScalarField& Qwall = const_cast<volScalarField&>(m.lookupObject<volScalarField>("Qwall_" + zoneName));
                const volScalarField& Tgas = m.lookupObject<volScalarField>("T");
                const scalar dw = %.9g;
                const scalar ks = %.9g;
                const scalar contact = %.9g;
                const bool radiation = %s;
                const scalar sigma = 5.670374e-8;
                wordList contactKeys;
%s
                const scalarField& V = m.V();
                const labelUList& own = m.owner();
                const labelUList& nei = m.neighbour();
                const surfaceScalarField& magSf = m.magSf();
                const volVectorField& C = m.C();
                kEff.primitiveFieldRef() = 0;
                Qwall.primitiveFieldRef() = 0;
                scalarField Cs(m.nCells(), 0);
                forAll(zoneCells, i)
                {
                    const label c = zoneCells[i];
                    const scalar phiC = max(solidFrac[c], scalar(0));
                    if (phiC > 1e-9)
                    {
                        kEff[c] = ks*phiC/3 + (radiation ? 16*sigma*pow3(Ts[c])*dw/(3*phiC) : 0);
                        Cs[c] = phiC*rhoCs*V[c];
                    }
                }
                kEff.correctBoundaryConditions();
                //face conductances (W/K): between zone cells (harmonic k), processor faces, wall contact
                scalarField Gi(own.size(), 0);
                forAll(own, f)
                {
                    const scalar kP = kEff[own[f]];
                    const scalar kN = kEff[nei[f]];
                    if (kP > 0 && kN > 0)
                    {
                        Gi[f] = 2*kP*kN/(kP + kN)*magSf[f]/mag(C[nei[f]] - C[own[f]]);
                    }
                }
                const fvBoundaryMesh& patches = m.boundary();
                List<scalarField> Gb(patches.size());
                boolList isContact(patches.size(), false);
                forAll(patches, p)
                {
                    const fvPatch& pp = patches[p];
                    const labelUList& fc = pp.faceCells();
                    Gb[p].setSize(pp.size(), 0);
                    const scalarField kP(kEff.boundaryField()[p].patchInternalField());
                    const vectorField d(pp.delta());
                    if (pp.coupled())
                    {
                        const scalarField kN(kEff.boundaryField()[p].patchNeighbourField());
                        forAll(fc, i)
                        {
                            if (kP[i] > 0 && kN[i] > 0)
                            {
                                Gb[p][i] = 2*kP[i]*kN[i]/(kP[i] + kN[i])*pp.magSf()[i]/mag(d[i]);
                            }
                        }
                        continue;
                    }
                    std::string pn(pp.name());
                    for (char& ch : pn)
                    {
                        ch = (ch >= 'A' && ch <= 'Z') ? char(ch + 32) : ch;
                    }
                    for (const word& k : contactKeys)
                    {
                        if (pn.find(k) != std::string::npos)
                        {
                            isContact[p] = true;
                        }
                    }
                    if (isContact[p])
                    {
                        forAll(fc, i)
                        {
                            if (kP[i] > 0)
                            {
                                Gb[p][i] = contact*kP[i]*pp.magSf()[i]/mag(d[i]);
                            }
                        }
                    }
                }
                //stable explicit sub-steps
                scalarField sumG(m.nCells(), 0);
                forAll(own, f)
                {
                    sumG[own[f]] += Gi[f];
                    sumG[nei[f]] += Gi[f];
                }
                forAll(patches, p)
                {
                    const labelUList& fc = patches[p].faceCells();
                    forAll(fc, i)
                    {
                        sumG[fc[i]] += Gb[p][i];
                    }
                }
                scalar rMax = 0;
                forAll(zoneCells, i)
                {
                    const label c = zoneCells[i];
                    if (Cs[c] > 0)
                    {
                        rMax = max(rMax, sumG[c]/Cs[c]);
                    }
                }
                reduce(rMax, maxOp<scalar>());
                const label nSub = max(label(1), label(std::ceil(dt*rMax/0.4)));
                const scalar dtSub = dt/nSub;
                scalarField dE(m.nCells(), 0);
                for (label s = 0; s < nSub; ++s)
                {
                    Ts.correctBoundaryConditions();
                    dE = 0;
                    forAll(own, f)
                    {
                        if (Gi[f] > 0)
                        {
                            const scalar q = Gi[f]*(Ts[nei[f]] - Ts[own[f]]);
                            dE[own[f]] += q;
                            dE[nei[f]] -= q;
                        }
                    }
                    forAll(patches, p)
                    {
                        const labelUList& fc = patches[p].faceCells();
                        if (patches[p].coupled())
                        {
                            const scalarField TN(Ts.boundaryField()[p].patchNeighbourField());
                            forAll(fc, i)
                            {
                                if (Gb[p][i] > 0)
                                {
                                    dE[fc[i]] += Gb[p][i]*(TN[i] - Ts[fc[i]]);
                                }
                            }
                        }
                        else if (isContact[p])
                        {
                            const scalarField& Tw = Tgas.boundaryField()[p];
                            forAll(fc, i)
                            {
                                if (Gb[p][i] > 0)
                                {
                                    const scalar q = Gb[p][i]*(Tw[i] - Ts[fc[i]]);
                                    dE[fc[i]] += q;
                                    Qwall[fc[i]] += q/V[fc[i]]/nSub;
                                }
                            }
                        }
                    }
                    forAll(zoneCells, i)
                    {
                        const label c = zoneCells[i];
                        if (Cs[c] > 0)
                        {
                            Ts[c] += dtSub*dE[c]/Cs[c];
                        }
                    }
                }
                Ts.correctBoundaryConditions();
            }''' % (fiber['dw'], fiber['ks'], fiber['contact'], 'true' if fiber['radiation'] else 'false', keys)


def porousFiberHeatSource(zone, props, cv):
    #scalarCodedSource on e for a random-fibre zone: as porousHeatSource (gas and fibres implicit together, heat into the
    #gas = heat the fibres give), with h*a per cell from the fibre fraction and the velocity relative to the fibres
    #(solidFrac_/Umat_ of the drag source) and the fibre heat capacity per cell solidFrac*rhos*cs; hEff_<zone> keeps the
    #per-cell coefficient of the step for the correct hook
    name = zone['name'].replace('-', '_')
    fiber = zone['fiber']
    consts, cell = fiberCellCode(zone, props, False, 0.0)
    rhoCs = fiber['rhos']*fiber.get('cs', 0.0) if zone['heat'] == 'solved' else 0.0
    head = '\n'.join([
        '            const fvMesh& m = mesh_;',
        '            const scalar Cv = %.9g;' % (cv),
        '            const scalar rhoCs = %.9g;' % (rhoCs),
        '            const scalar dt = m.time().deltaTValue();',
        '            const word zoneName("%s");' % (zone['name']),
        '            const labelList& zoneCells = this->cells();'])
    lookups = '\n'.join(['            volScalarField& %s = const_cast<volScalarField&>(m.lookupObject<volScalarField>'
                         '("%s_" + zoneName));' % (f, f) for f in ('Ts', 'TsOld', 'Qs', 'eLin', 'TLin', 'hEff')])
    addSup = '''        codeAddSupRho
        #{
%s
%s
            if (!m.foundObject<volScalarField>("solidFrac_" + zoneName))
            {
                return;
            }
            const volScalarField& solidFrac = m.lookupObject<volScalarField>("solidFrac_" + zoneName);
            const volVectorField& Umat = m.lookupObject<volVectorField>("Umat_" + zoneName);
            const volVectorField& Ugas = m.lookupObject<volVectorField>("U");
            const volScalarField& T = m.lookupObject<volScalarField>("T");
            if (!m.foundObject<volScalarField>("Ts_" + zoneName))
            {
                const word startName(m.time().timeName(m.time().startTime().value()));
                IOobject TsIO("Ts_" + zoneName, startName, m, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE);
                volScalarField* TsPtr = new volScalarField(TsIO, m, dimensionedScalar(dimTemperature, %.9g));
                TsPtr->store();
                const word now(m.time().timeName());
                (new volScalarField(IOobject("TsOld_" + zoneName, now, m), *TsPtr))->store();
                (new volScalarField(IOobject("Qs_" + zoneName, now, m), m,
                    dimensionedScalar(dimPower/dimVolume, Zero)))->store();
                (new volScalarField(IOobject("eLin_" + zoneName, now, m), m,
                    dimensionedScalar(dimEnergy/dimMass, Zero)))->store();
                (new volScalarField(IOobject("TLin_" + zoneName, now, m), m,
                    dimensionedScalar(dimTemperature, Zero)))->store();
                (new volScalarField(IOobject("hEff_" + zoneName, now, m), m,
                    dimensionedScalar(dimPower/dimVolume/dimTemperature, Zero)))->store();
                (new volScalarField(IOobject("kEff_" + zoneName, now, m), m,
                    dimensionedScalar(dimPower/dimLength/dimTemperature, Zero)))->store();
                (new volScalarField(IOobject("Qwall_" + zoneName, now, m), m,
                    dimensionedScalar(dimPower/dimVolume, Zero)))->store();
            }
%s
            static label lastIndex = -1;
            if (m.time().timeIndex() != lastIndex)
            {
                lastIndex = m.time().timeIndex();
                forAll(zoneCells, i)
                {
                    TsOld[zoneCells[i]] = Ts[zoneCells[i]];
                }
            }
            const scalarField& V = m.V();
            const scalarField& e = eqn.psi();
            scalarField& diag = eqn.diag();
            scalarField& source = eqn.source();
            forAll(zoneCells, i)
            {
                const label c = zoneCells[i];
                const vector Ur(Ugas[c] - Umat[c]);
%s
                const scalar Cs = solidFrac[c]*rhoCs;
                hEff[c] = (Cs > 0 ? haC*Cs/(Cs + haC*dt) : haC);
                const scalar k = hEff[c]/Cv;
                eLin[c] = e[c];
                TLin[c] = T[c];
                diag[c] -= k*V[c];
                source[c] -= (hEff[c]*(TsOld[c] - T[c]) + k*e[c])*V[c];
            }
        #};''' % (head, consts, zone['Ts'][0], lookups, cell)
    correct = '''        codeCorrect
        #{
%s
            if (!m.foundObject<volScalarField>("Ts_" + zoneName))
            {
                return;
            }
            const volScalarField& solidFrac = m.lookupObject<volScalarField>("solidFrac_" + zoneName);
%s
            forAll(zoneCells, i)
            {
                const label c = zoneCells[i];
                const scalar T1 = TLin[c] + (fld[c] - eLin[c])/Cv;
                Qs[c] = hEff[c]*(TsOld[c] - T1);
                const scalar Cs = solidFrac[c]*rhoCs;
                if (Cs > 0)
                {
                    Ts[c] = TsOld[c] - Qs[c]*dt/Cs;
                }
            }
%s
        #};''' % (head, lookups, fiberConductionCode(zone))
    lines = ['%sHeat' % (name), '{', '    type            scalarCodedSource;', '    name            %sHeat;' % (name),
             '    scalarCodedSourceCoeffs', '    {', '        selectionMode   cellZone;',
             '        cellZone        %s;' % (zone['name']), '        fields          (e);', '',
             '        codeInclude', '        #{', '            #include "fvCFD.H"', '        #};', '',
             correct, '', '        codeAddSup #{ #};', '        codeConstrain #{ #};', '', addSup,
             '    }', '}']
    return '\n'.join(lines)


def porousDragHeatSource(zone, props, moving):
    #scalarCodedSource on e: drag work Ur & C & Ur (W/m3) into the gas, the porosity sources drop it
    #Ur = U - Umat for a zone moving with a body (Umat_<zone> from its drag source), else U
    name = zone['name'].replace('-', '_')
    constants, cell = porousDragCoeff(zone, props, False, 0.0)
    head = '\n'.join([
        '            const fvMesh& m = mesh_;',
        '            const word zoneName("%s");' % (zone['name']),
        '            const labelList& zoneCells = this->cells();',
        constants])
    umat = ('m.lookupObject<volVectorField>("Umat_" + zoneName)[c]' if moving else 'vector(Zero)')
    addSup = '''        codeAddSupRho
        #{
%s
            if (!m.foundObject<volScalarField>("porDragHeat_" + zoneName))
            {
                (new volScalarField(IOobject("porDragHeat_" + zoneName, m.time().timeName(), m), m,
                    dimensionedScalar(dimPower/dimVolume, Zero)))->store();
            }
            volScalarField& heat = const_cast<volScalarField&>(m.lookupObject<volScalarField>("porDragHeat_" + zoneName));
            if (%s)
            {
                return;
            }
            const volVectorField& U = m.lookupObject<volVectorField>("U");
            const scalarField& V = m.V();
            scalarField& source = eqn.source();
            forAll(zoneCells, i)
            {
                const label c = zoneCells[i];
                const vector Ur(U[c] - %s);
%s
                heat[c] = Ur & (C & Ur);
                source[c] -= heat[c]*V[c];
            }
        #};''' % (head, '!m.foundObject<volVectorField>("Umat_" + zoneName)' if moving else 'false', umat, cell)
    lines = ['%sDragHeat' % (name), '{', '    type            scalarCodedSource;', '    name            %sDragHeat;' % (name),
             '    scalarCodedSourceCoeffs', '    {', '        selectionMode   cellZone;',
             '        cellZone        %s;' % (zone['name']), '        fields          (e);', '',
             '        codeInclude', '        #{', '            #include "fvCFD.H"', '        #};', '',
             '        codeCorrect #{ #};', '        codeAddSup #{ #};', '        codeConstrain #{ #};', '', addSup,
             '    }', '}']
    return '\n'.join(lines)


def getWarmStart(fullCaseSetupDict, porous):
    #INIT_T_FIELD: uniform (INIT_T everywhere) | matrix (the only POR zone with a linear SOLID_T) | <POR zone>
    #warm start: gas T along z as the matrix start (SOLID_T over the zone's z range, constant beyond its ends), so the
    #hot and cold end spaces start at the end temperatures; p uniform, raised so the gas mass is that of the charge state
    #(CHARGE_PRESSURE at INIT_T). Done at run time by remeshMassFix -warmStart before the first solve
    setting = fullCaseSetupDict['GLOBAL_MATERIAL'].get('INIT_T_FIELD', ['uniform'])[0]
    if setting == 'uniform':
        return None
    if getFlow(fullCaseSetupDict) == 'incompressible':
        print('\t\tNOTE: incompressible run, INIT_T_FIELD %s ignored (no T)' % (setting))
        return None
    linear = [zone for zone in porous if zone.get('heat', 'none') != 'none' and zone['Ts'][0] != zone['Ts'][1]]
    if setting == 'matrix':
        if len(linear) != 1:
            sys.exit('ERROR! GLOBAL_MATERIAL INIT_T_FIELD matrix needs exactly one POR zone with HEAT and two SOLID_T '
                     'values (found %d): give the zone name instead' % (len(linear)))
        zone = linear[0]
    else:
        match = [zone for zone in porous if zone['name'] == setting]
        if not match or match[0].get('heat', 'none') == 'none':
            sys.exit('ERROR! GLOBAL_MATERIAL INIT_T_FIELD must be uniform, matrix or a POR zone with HEAT (zones: %s)'
                     % (' '.join(zone['name'] for zone in porous) or 'none'))
        zone = match[0]
    warm = {'zone': zone['name'], 'T': list(zone['Ts'])}
    print('\n\tWarm start: gas T %g K (-z) to %g K (+z) along %s, constant beyond its ends; p raised to keep the gas '
          'mass of %s Pa at %s K (remeshMassFix -warmStart before the first solve)'
          % (warm['T'][0], warm['T'][1], zone['name'], getVal(fullCaseSetupDict, 'GLOBAL_MATERIAL', 'CHARGE_PRESSURE'),
             getVal(fullCaseSetupDict, 'GLOBAL_MATERIAL', 'INIT_T')))
    return warm


def writeFvOptionsInternal(porous, fullCaseSetupDict, porOwners=None):
    #constant/fvOptions: per POR zone explicitPorositySource (static) or the coded relative-velocity drag (moving with a
    #body), the heat source when HEAT is fixed or solved, the drag work (DISSIPATION); viscous dissipation for the gas
    path = os.path.join('constant', 'fvOptions')
    dissipation = getDissipation(fullCaseSetupDict)
    if not porous and not dissipation:
        if os.path.isfile(path):
            os.remove(path)
        return
    print('\n\tWriting fvOptions...')
    props = getGasProperties(fullCaseSetupDict)
    cv = props['CP'] - 8314.46/props['MOL_WEIGHT']
    incompressible = getFlow(fullCaseSetupDict) == 'incompressible'
    nu = getNu(fullCaseSetupDict) if incompressible else 0.0
    rhoInf = getGasDensity(fullCaseSetupDict)
    entries = []
    coded = False
    for zone in porous:
        name = zone['name'].replace('-', '_')
        moving = zone['name'] in (porOwners or {})
        if 'fiber' in zone:
            entries.append(porousFiberDragSource(zone, props, incompressible, nu, rhoInf))
            print('\t\t%s: random-fibre drag and heat per cell (porosity from the fibre fraction, velocity relative to '
                  'the fibres)' % (zone['name']))
            if zone['heat'] != 'none':
                entries.append(porousFiberHeatSource(zone, props, cv))
            coded = True
            continue
        local = zone.get('wire', {}).get('local', False)
        if moving or local:
            entries.append(porousDragSource(zone, props, incompressible, nu, rhoInf, moving))
            if moving:
                print('\t\t%s: drag on the velocity relative to the matrix (moves with %s), force added to the body'
                      % (zone['name'], porOwners[zone['name']].name))
            if local:
                print('\t\t%s: wire-mesh Forchheimer term%s per cell from the local Reynolds number'
                      % (zone['name'], ' and h*a' if zone.get('haLocal') and zone['heat'] != 'none' else ''))
            coded = True
        else:
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
            entries.append(porousHeatSource(zone, cv, props))
            coded = True
        if dissipation:
            entries.append(porousDragHeatSource(zone, props, moving or local))
            coded = True
    if dissipation:
        #ESI fvOption: -devRhoReff && grad(U) into the energy equation (he from the thermo), effective viscosity
        entries.append('viscousDissipation\n{\n    type            viscousDissipation;\n    enabled         true;\n}')
        print('\t\tviscous dissipation%s into the gas energy (DISSIPATION yes)' % (' and porous drag work' if porous else ''))
    with open(path, 'w') as f:
        f.write('FoamFile { version 2.0; format ascii; class dictionary; object fvOptions; }\n\n')
        f.write('//written by internalSetup\n\n%s\n' % ('\n\n'.join(entries)))
    if coded:
        print('\t\tthe coded sources are compiled at run time (needs a C++ compiler where the solver runs)')


def writeTurbulenceInternal(templateLoc, fullCaseSetupDict):
    turbModel = getTurbulence(fullCaseSetupDict)
    print('\n\tWriting turbulenceProperties (%s)...' % (turbModel))
    if turbModel == 'laminar':
        replaceDict = {'SIMULATION_TYPE': 'laminar', 'RAS_BLOCK': ''}
    else:
        replaceDict = {'SIMULATION_TYPE': 'RAS',
                       'RAS_BLOCK': 'RAS { RASModel %s; turbulence on; printCoeffs on; }' % (turbModel)}
    copyTemplate(templateLoc, 'constant/turbulenceProperties', 'constant/turbulenceProperties', replaceDict)


def getEndTime(fullCaseSetupDict):
    #END_TIME in GLOBAL_CONTROL (s): a transient run that is not cyclic, None when blank (run length from the cycles)
    value = ''.join(fullCaseSetupDict['GLOBAL_CONTROL'].get('END_TIME', [''])).strip()
    if value == '':
        return None
    try:
        endTime = float(value)
    except ValueError:
        sys.exit('ERROR! GLOBAL_CONTROL END_TIME must be a time in s or blank (got "%s")!' % (value))
    if endTime <= 0:
        sys.exit('ERROR! GLOBAL_CONTROL END_TIME must be positive!')
    if getMoveSections(fullCaseSetupDict):
        sys.exit('ERROR! GLOBAL_CONTROL END_TIME is for cases without moving bodies: with MOVE bodies the run '
                 'length is N_CYCLES at FREQ (leave END_TIME blank)!')
    return endTime


def getRunTimes(fullCaseSetupDict):
    #endTime and writeInterval from FREQ, N_CYCLES, WRITE_PER_CYCLE, or END_TIME, WRITE_INTERVAL (not cyclic)
    control = fullCaseSetupDict['GLOBAL_CONTROL']
    endTime = getEndTime(fullCaseSetupDict)
    if endTime is not None:
        value = ''.join(control.get('WRITE_INTERVAL', [''])).strip()
        try:
            writeInterval = float(value) if value else endTime/20.0
        except ValueError:
            sys.exit('ERROR! GLOBAL_CONTROL WRITE_INTERVAL must be a time in s or blank (got "%s")!' % (value))
        if writeInterval <= 0:
            sys.exit('ERROR! GLOBAL_CONTROL WRITE_INTERVAL must be positive!')
        return endTime, writeInterval
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
    cyclic = getEndTime(fullCaseSetupDict) is None
    if cyclic:
        print('\t\t%s cycles at %s Hz: endTime %.6g s, write every %.6g s'
              % (control['N_CYCLES'][0], control['FREQ'][0], endTime, writeInterval))
    else:
        print('\t\tnot cyclic (END_TIME): endTime %.6g s, write every %.6g s' % (endTime, writeInterval))
    initDt, maxCo, maxDt, adjustDt = getTimeStep(fullCaseSetupDict)
    if adjustDt == 'yes':
        print('\t\tadjustable time step: maxCo %s, maxDeltaT %s' % (maxCo, maxDt))
    elif cyclic:
        print('\t\tfixed time step %s (MAX_CO 0 or blank), %.6g steps per cycle'
              % (maxDt, 1.0/(float(control['FREQ'][0])*float(maxDt))))
    else:
        print('\t\tfixed time step %s (MAX_CO 0 or blank), %.6g steps' % (maxDt, endTime/float(maxDt)))
        #without time step adjustment the writes only land on the interval if it is a whole number of steps
        steps = writeInterval/float(maxDt)
        if abs(steps - round(steps)) > 1e-6*steps:
            print('\t\tWARNING! write interval %.6g s is not a whole number of time steps (%.4g): write times '
                  'will be the first step past each interval' % (writeInterval, steps))
    #one outer corrector (PISO) is only accurate around Co 1 on the moving mesh; a larger limit (or none, fixed step)
    #lets the step run at MAX_DT whatever the Courant number does (Co 10 seen in end gaps just after a remesh map)
    nOuter = int(fullCaseSetupDict['GLOBAL_SIM_CONTROL'].get('N_OUTER_CORR', ['1'])[0])
    if nOuter == 1 and (adjustDt == 'no' or float(maxCo) > 1):
        print('\t\tWARNING: N_OUTER_CORR 1 (PISO) with %s: the Courant number is not held near 1. Check COURANT_MAX '
              'and COURANT_ABOVE_1 in the summary (or the log); if Co stays above ~1, lower MAX_DT or use '
              'N_OUTER_CORR 3 (no relaxation)' % ('a fixed time step' if adjustDt == 'no' else 'MAX_CO %s' % (maxCo)))
    replaceDict = {'SOLVER': getVal(fullCaseSetupDict, 'GLOBAL_SIM_CONTROL', 'SOLVER'),
                   'STARTFROM': control['STARTFROM'][0], 'ENDTIME': '%.9g' % endTime,
                   'INIT_DT': initDt, 'WRITE_INTERVAL': '%.9g' % writeInterval,
                   'MAX_CO': maxCo, 'MAX_DT': maxDt, 'ADJUST_DT': adjustDt}
    copyTemplate(templateLoc, 'system/controlDict', 'system/controlDict', replaceDict)
    #moving mesh: static walls use slidingWallVelocity (foamLibs/stirlingMotion, compiled by the run scripts)
    if getMoveSections(fullCaseSetupDict):
        with open('system/controlDict') as f:
            text = f.read()
        text = text.replace('\nfunctions\n', '\nlibs            (stirlingMotion);   //slidingWallVelocity on the static walls '
                                              '(foamLibs/stirlingMotion)\n\nfunctions\n', 1)
        with open('system/controlDict', 'w') as f:
            f.write(text)


def getPMaxIter(fullCaseSetupDict):
    #cap on GAMG iterations for p: a near-static gas (very slow piston) stalls at a round-off floor of the normalised
    #residual (~2e-4) far above the 1e-8 tolerance, and every solve then ran the 1000 default
    val = ''.join(fullCaseSetupDict['GLOBAL_SIM_CONTROL'].get('P_MAX_ITER', ['1000'])).strip() or '1000'
    try:
        if int(val) < 1:
            raise ValueError
    except ValueError:
        sys.exit('ERROR! GLOBAL_SIM_CONTROL P_MAX_ITER must be a whole number >= 1 (got "%s")' % (val))
    return val


def getRelaxation(fullCaseSetupDict):
    #under-relaxation: steady solvers always, transient only with outer correctors (N_OUTER_CORR > 1)
    #the last outer iteration (Final) is not relaxed. RELAX_* auto: 0.3/0.7/0.7 steady, 1 (none) transient: with p
    #relaxed, rho = psi p lags continuity in the relaxed iterations and the final e solve still uses that rho, so the gas
    #energy goes wrong by about the lag times rho e every step (half the compression heating lost in a closed-cylinder
    #test with 3 outer correctors and p 0.3; within 0.2 % without relaxation, same cost)
    simControl = fullCaseSetupDict['GLOBAL_SIM_CONTROL']
    nOuter = int(simControl['N_OUTER_CORR'][0])
    if nOuter < 1:
        sys.exit('ERROR! GLOBAL_SIM_CONTROL N_OUTER_CORR must be 1 or more!')
    steady = isSteady(fullCaseSetupDict)
    if nOuter == 1 and not steady:
        return ''
    relax = {}
    for key, default in (('RELAX_P', '0.3'), ('RELAX_U', '0.7'), ('RELAX_E', '0.7')):
        value = simControl.get(key, ['auto'])[0]
        if value == 'auto':
            value = default if steady else '1'
        try:
            relax[key] = float(value)
        except ValueError:
            relax[key] = -1.0
        if not 0 < relax[key] <= 1:
            sys.exit('ERROR! GLOBAL_SIM_CONTROL %s must be auto or > 0 and <= 1!' % (key))
    compressible = getFlow(fullCaseSetupDict) == 'compressible'
    if not steady:
        if all(relax[key] == 1 for key in relax):
            print('\t\t%d outer correctors, no under-relaxation' % (nOuter))
            return ''
        #p relaxed: the gas energy error per process is about the p lag left after the relaxed outer iterations
        if compressible and relax['RELAX_P'] < 1:
            print('\t\tWARNING: p relaxed in a transient compressible run: rho lags p in the relaxed outer iterations '
                  'and the gas energy change (compression, heating) is off by about (1 - RELAX_P)^(N_OUTER_CORR - 1) '
                  '= %.1f %%; use RELAX_P auto (= 1), or more outer correctors' % (100*(1 - relax['RELAX_P'])**(nOuter - 1)))
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
    lines += outerResidualControl(fullCaseSetupDict)
    return '\n'.join(lines + ['}'])


def outerResidualControl(fullCaseSetupDict):
    #PIMPLE residualControl: with several outer correctors a time step ends once the initial residuals of the latest
    #solves are below OUTER_TOL_P (pressure) and OUTER_TOL (U, e), then one final iteration (v2412 pimpleControl), so
    #N_OUTER_CORR is an upper limit for the hard steps, not the cost of every step. Blank = always all of them
    simControl = fullCaseSetupDict['GLOBAL_SIM_CONTROL']
    if int(simControl['N_OUTER_CORR'][0]) < 2:
        return []
    tolP = ''.join(simControl.get('OUTER_TOL_P', ['1e-4'])).strip()
    tol = ''.join(simControl.get('OUTER_TOL', ['1e-5'])).strip()
    entries = []
    if tolP:
        entries.append(('p_rgh' if isBuoyant(fullCaseSetupDict) else 'p', float(tolP)))
    if tol:
        entries.append(('U', float(tol)))
        if getFlow(fullCaseSetupDict) != 'incompressible':
            entries.append(('e', float(tol)))
    if not entries:
        return []
    return (['    residualControl', '    {'] +
            ['        %-6s { tolerance %g; relTol 0; }' % (name, value) for name, value in entries] + ['    }'])


def writeNumericsInternal(templateLoc, fullCaseSetupDict, pressureBC=False):
    #pressureBC: an inlet or outlet fixes the pressure level, so no pressure reference is needed
    print('\n\tWriting fvSchemes, fvSolution, decomposeParDict...')
    compute = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP']
    steady = isSteady(fullCaseSetupDict)
    #divergence schemes per flow type, closed incompressible volume needs a pressure reference
    #linearUpwind keeps it stable at large Co (010 diverged at Co 2.4 with limitedLinearV)
    #steady: bounded convection terms (continuity is only met at convergence)
    bounded = 'bounded ' if steady else ''
    #turbulence fields upwind: linearUpwind took omega below zero next to the walls in every step (kOmegaSST, static
    #demo: bounding each step, min -42; upwind: none, mean gas T within 1e-3 K), and negative omega next to the large
    #wall values of omegaWallFunction blew up the omega solve in an engine run
    divTurb = ['div(phi,k)      %sGauss upwind;' % (bounded),
               'div(phi,omega)  %sGauss upwind;' % (bounded),
               'div(phi,epsilon) %sGauss upwind;' % (bounded)]
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
                  'P_RELTOL': '0.05' if steady else '0.01', 'U_RELTOL': '0.1' if steady else '0',
                  'P_MAX_ITER': getPMaxIter(fullCaseSetupDict)})
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


def minMaxFields(fullCaseSetupDict):
    #fields solved for this flow type and turbulence model (fieldMinMax)
    fields = ['U', 'p']
    if isBuoyant(fullCaseSetupDict):
        fields.append('p_rgh')
    compressible = getFlow(fullCaseSetupDict) != 'incompressible'
    if compressible:
        fields += ['T', 'rho']
    turbModel = getTurbulence(fullCaseSetupDict)
    if turbModel != 'laminar':
        fields += ['k', 'epsilon' if turbModel == 'kEpsilon' else 'omega', 'nut'] + (['alphat'] if compressible else [])
    return fields


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
    #min/max and their cell locations of every solved field each step: shows where a run starts to blow up
    fos.append(functionObject('fieldMinMax', [('type', 'fieldMinMax;'), ('libs', '(fieldFunctionObjects);'),
                                              ('mode', 'magnitude;'), ('location', 'true;'),
                                              ('fields', '(%s);' % (' '.join(minMaxFields(fullCaseSetupDict))))]))
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
            fos.append(functionObject('porHeatToGas_%s' % (zone['name']), cellZone + [('operation', 'volIntegrate;'),
                                                                               ('fields', '(Qs_%s);' % (zone['name']))]))
        #drag work into the gas (W) and, for a matrix moving with a body, the gas force on the matrix (N)
        if getDissipation(fullCaseSetupDict):
            fos.append(functionObject('porDragHeat_%s' % (zone['name']), cellZone + [('operation', 'volIntegrate;'),
                                                                              ('fields', '(porDragHeat_%s);' % (zone['name']))]))
        if motionInfo and zone['name'] in motionInfo.get('porOwners', {}) or 'fiber' in zone:
            fos.append(functionObject('porDrag_%s' % (zone['name']), cellZone + [('operation', 'volIntegrate;'),
                                                                          ('fields', '(porDrag_%s);' % (zone['name']))]))
        #random fibres: fibre volume in the zone (constant: the fibres stay in their cells) and the mean fraction
        if 'fiber' in zone:
            fos.append(functionObject('fiber_%s' % (zone['name']), cellZone + [('operation', 'volIntegrate;'),
                                                                        ('fields', '(solidFrac_%s);' % (zone['name']))]))
            #heat from the walls the fibres touch (W, into the fibres; the gas gets it through porHeatToGas)
            if zone['heat'] == 'solved' and zone['fiber']['contactPids']:
                fos.append(functionObject('porWallToFiber_%s' % (zone['name']),
                                          cellZone + [('operation', 'volIntegrate;'),
                                                      ('fields', '(Qwall_%s);' % (zone['name']))]))
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
        #one regex per body: a PID that snappy left without faces (e.g. a nogap displacer side on the chamber wall)
        #has no patch, separate names would each warn at run time
        for body in motionInfo['bodies']:
            patches = '"%s_(%s)"' % (body.source, '|'.join(body.movedPids))
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


def getRemesh(fullCaseSetupDict):
    #remeshing (MESH_MOTION REMESH auto): squash/stretch limit per mesh, None = one mesh for the whole run
    meshMotion = fullCaseSetupDict.get('MESH_MOTION', {})
    remesh = meshMotion.get('REMESH', ['no'])[0].lower()
    if remesh in ('no', ''):
        return None
    if remesh != 'auto':
        sys.exit('ERROR! MESH_MOTION REMESH must be no or auto!')
    try:
        limit = float(meshMotion.get('REMESH_DEFORMATION', ['3'])[0])
    except ValueError:
        limit = 0.0
    if limit <= 1.0:
        sys.exit('ERROR! MESH_MOTION REMESH_DEFORMATION must be a number above 1 (squash/stretch allowed per mesh, e.g. 3)')
    return limit


def getCrank(fullCaseSetupDict):
    #crank drive settings, None for prescribed motion
    meshMotion = fullCaseSetupDict['MESH_MOTION']
    drive = meshMotion.get('DRIVE', ['prescribed'])[0]
    load = None
    if drive == 'prescribed':
        if getRemesh(fullCaseSetupDict) is None:
            return None
        #remeshing stops each run at a crank angle: prescribed motion runs as the crank drive at a fixed speed
        drive, load = 'crank', 'speed'
    if drive != 'crank':
        sys.exit('ERROR! MESH_MOTION DRIVE must be prescribed or crank!')
    crank = {'LOAD': load or meshMotion['LOAD'][0]}
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


def crankEntries(crank, startAngle, stopAngle=None):
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
    if stopAngle is not None:
        lines.append('stopAngle       %.9g;   //remeshing: write and stop at this crank angle (deg), the next mesh '
                     'is made there' % (np.mod(stopAngle, 360.0)))
    return '\n'.join('        %s' % (line) for line in lines)


def segmentNames(segments):
    #remeshing: one mesh per segment, seg01 is the case mesh (start angle)
    return ['seg%02d' % (k + 1) for k in range(len(segments))]


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
    #MAX_DEFORMATION blank or 0: no limit (the ratio is still reported)
    try:
        maxDeformation = float(''.join(meshMotion.get('MAX_DEFORMATION', ['']))) or np.inf
    except ValueError:
        if ''.join(meshMotion.get('MAX_DEFORMATION', [''])).strip() != '':
            sys.exit('ERROR! MESH_MOTION MAX_DEFORMATION must be a number, 0 or blank (no limit)!')
        maxDeformation = np.inf
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
        mode = setup.get('CLEARANCE_MODE', ['soft'])[0]
        if mode not in ('soft', 'scale', 'shift'):
            sys.exit('ERROR! [%s] CLEARANCE_MODE must be soft, scale or shift!' % (section))
        body = motionUtilities.Body(name, thetaC, pos, float(setup['REF_ANGLE'][0]), float(setup['TDC_CLEARANCE'][0]),
                                    mode)
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
    #soft / scale: the profile ends that come too close to an anchor are adjusted, shift: the whole path is offset
    adjusted = motionUtilities.adjustProfiles(zones, bodies, step)
    zoneArea = dict((body.name, max([zone.area for zone in zones if zone.body is body] or [0.0])) for body in bodies)
    for bodyName, info in adjusted.items():
        ends = ' and '.join(text for text, value in [('low end +%.4g m' % (info['raise']), info['raise']),
                                                     ('high end -%.4g m' % (info['lower']), info['lower'])] if value > 0)
        print('\t\t%s profile adjusted for TDC_CLEARANCE (%s): %s, stroke %.4g -> %.4g m, max change %.4g m at %.1f deg, '
              'more than 0.02 mm over %.0f deg of the cycle, dead volume +%.4g cm3'
              % (bodyName, info['mode'], ends, info['strokeBefore'], info['strokeAfter'], info['maxDev'],
                 info['devAngle'], info['changedDeg'], (info['raise'] + info['lower'])*zoneArea[bodyName]*1e6))
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
    #remeshing: a new mesh wherever a gap has been squashed or stretched REMESH_DEFORMATION times since the last one
    segments = None
    remesh = getRemesh(fullCaseSetupDict)
    if remesh is not None:
        segments = motionUtilities.remeshAngles(zones, startAngle, step, remesh)
        if len(segments) == 1:
            print('\t\tremeshing: squash/stretch stays below REMESH_DEFORMATION %g over the cycle, one mesh' % (remesh))
            segments = None
        else:
            if any(geomName.startswith('FAN') for geomName in geomFiles):
                sys.exit('ERROR! Remeshing with FAN baffles is not supported yet (the baffles are made on the case mesh only)')
            print('\t\tremeshing: %d meshes per cycle (REMESH_DEFORMATION %g), each made at its start angle:'
                  % (len(segments), remesh))
            for name, (a0, a1) in zip(segmentNames(segments), segments):
                worst = motionUtilities.segmentRatio(zones, a0, a1, step)
                print('\t\t\t%s  %5.1f to %5.1f deg  worst squash/stretch %.2fx (%s - %s)'
                      % (name, np.mod(a0, 360.0), np.mod(a1, 360.0), worst[0], worst[1][0] if worst[1] else '-',
                         worst[1][1] if worst[1] else '-'))
            #fixed time step: OpenFOAM restarts with the controlDict deltaT (uniform/time is only read when adjustable)
            if getTimeStep(fullCaseSetupDict)[3] == 'no':
                print('\t\tNOTE: fixed time step, every segment restarts at the full MAX_DT (REMESH_DT is not used). '
                      'If a run diverges just after a map, use an adjustable step (MAX_CO > 0) to ramp up from REMESH_DT')
    if ratio > maxDeformation and segments is None:
        sys.exit('ERROR! Squash/stretch %.2fx exceeds MAX_DEFORMATION %.2f: increase TDC_CLEARANCE or the limit'
                 % (ratio, maxDeformation))
    #very flat cells near TDC: a run with a 76x squashed piston bore diverged at TDC without p under-relaxation
    if ratio > 20 and segments is None:
        print('\t\tWARNING: cells are squashed/stretched %.0fx over the cycle: runs can diverge where the gap is '
              'smallest; use START_ANGLE auto (meshes at the geometric mean gap), REMESH auto or a larger TDC_CLEARANCE'
              % (ratio))
    meshAngles = [startAngle] + ([a0 for a0, a1 in segments[1:]] if segments else [])
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

    #a body made of PIDs in a wall file (piston crown, FROM) must seal against its liner at every meshing angle: inside
    #the liner along z and covering the bore (after meshing the liner boundary stretches with the mesh)
    for body in bodies:
        if body.source.startswith('MOVE'):
            continue
        walls = dict((pid, t) for pid, t in tris[body.source].items() if pid not in body.movedPids)
        for zone in [z for z in zones if z.body is body]:
            pid = zone.patch[len(body.source) + 1:]
            for angle in meshAngles:
                problems = motionUtilities.sealProblems(tris[body.source][pid], float(motionUtilities.zoneZ(zone, angle)),
                                                        walls, flatTol)
                if problems:
                    sys.exit('ERROR! %s does not seal at the meshing angle %.1f deg: %s'
                             % (zone.patch, np.mod(angle, 360.0), '; '.join(problems)))

    #per meshing angle: LOC_IN_MESH moved with the mesh (segments only, the case mesh keeps the given point) and z-only
    #refinement of end spaces that would get fewer than GAP_CELLS_Z cells across (a 3 mm end gap at the stroke end
    #got one 2 mm cell, the deformed mesh had about 8 flat ones)
    try:
        loc = np.array([float(v) for v in fullCaseSetupDict['GLOBAL_REFINEMENT']['LOC_IN_MESH'] if v != ''])
    except (KeyError, ValueError):
        loc = None
    if loc is None or len(loc) != 3:
        sys.exit('ERROR! GLOBAL_REFINEMENT LOC_IN_MESH must be a point x y z in the gas (got: %s)'
                 % (' '.join(fullCaseSetupDict.get('GLOBAL_REFINEMENT', {}).get('LOC_IN_MESH', ['missing']))))
    gapCellsZ = float(''.join(meshMotion.get('GAP_CELLS_Z', ['6'])) or 0)
    levels = dict((geomName, int(geomDict[geom]['refinement'])) for geomName, geom in geomFiles.items())
    def cellOf(zone):
        geomName = [g for g in geomFiles if zone.patch.startswith(g + '_')]
        return baseCell/2**(levels[max(geomName, key=len)] if geomName else 0)
    meshSetups = []
    for k, angle in enumerate(meshAngles):
        point = list(loc) if k == 0 else motionUtilities.locationAt(loc, zones, startAngle, angle, step)
        bands = motionUtilities.gapBands(zones, angle, step, cellOf, gapCellsZ)
        for lower, upper, zLo, zHi, splits, before, after in bands:
            print('\t\t%.1f deg: %s - %s gap %.4g m, %.1f cells across -> %d z splits, %.1f cells'
                  % (np.mod(angle, 360.0), lower.name, upper.name, zHi - zLo, before, splits, after))
            if after < gapCellsZ:
                print('\t\tWARNING: fewer than GAP_CELLS_Z %g cells across that gap (at most 4 z splits)' % (gapCellsZ))
        meshSetups.append({'angle': angle, 'loc': point,
                           'bands': [(zLo - 0.5*max(cellOf(lower), cellOf(upper)), zHi + 0.5*max(cellOf(lower), cellOf(upper)), splits)
                                     for lower, upper, zLo, zHi, splits, before, after in bands]})

    #LOC_IN_MESH must stay in the gas at every meshing angle: warn when it falls inside the box of a body of its own
    if loc is not None and len(loc) == 3:
        for body in bodies:
            if not body.source.startswith('MOVE'):
                continue
            allTris = np.vstack([t.reshape(-1, 3) for t in tris[body.source].values()])
            for angle, meshSetup in zip(meshAngles, meshSetups):
                loc = np.array(meshSetup['loc'])
                dz = float(body.shift(angle))
                lo, hi = allTris.min(axis=0) + [0.0, 0.0, dz], allTris.max(axis=0) + [0.0, 0.0, dz]
                if np.all(loc > lo) and np.all(loc < hi):
                    print('\t\tWARNING: LOC_IN_MESH %s is inside the bounding box of %s at the meshing angle %.1f deg: '
                          'make sure it is in the gas (not in the body) at every meshing angle'
                          % (' '.join('%g' % v for v in loc), body.name, np.mod(angle, 360.0)))
                    break

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

    #porous zones inside a moving body (regenerator in the displacer) are moved with it (MOVES_WITH = <body>);
    #MOVES_WITH mesh: the matrix deforms with the mesh (steel wool in an end space, squashed by the displacer), its
    #geometry is deformed in z like the mesh to every meshing angle
    bodyNames = dict((body.name, body) for body in bodies)
    porOwners = {}
    porDeform = {}
    for geomName in geomFiles:
        if not geomName.startswith('POR'):
            continue
        owner = fullCaseSetupDict[geomName].get('MOVES_WITH', ['none'])[0]
        if owner == 'none':
            continue
        if owner == 'mesh':
            porDeform[geomName] = tris[geomName]
            z = np.vstack(list(tris[geomName].values()))[:, :, 2]
            dz = motionUtilities.meshDisplacementZ(z, zones, startAngle)
            print('\t\t%s deforms with the mesh: z %.6g..%.6g m at the start angle (geometry %.6g..%.6g m)'
                  % (geomName, float((z + dz).min()), float((z + dz).max()), float(z.min()), float(z.max())))
            continue
        if owner not in bodyNames:
            sys.exit('ERROR! [%s] MOVES_WITH %s is not a MOVE body (bodies: %s)' % (geomName, owner, list(bodyNames)))
        porOwners[geomName] = bodyNames[owner]
        print('\t\t%s moves with %s: shifted %+.4g m to the start angle' % (geomName, owner,
                                                                           float(bodyNames[owner].shift(startAngle))))

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
        baseGap = {pid: tris[body.source][pid] for pid in matchPids(keywords, pidDict[body.source], '[%s] GAP_PIDS' % (section))}
        gap = motionUtilities.gapWidth(np.vstack(list(baseGap.values())) + np.array([0.0, 0.0, dz]), np.vstack(wallGap))
        gapCells = float(setup['GAP_CELLS'][0])
        level, distance, cell = motionUtilities.gapRefinement(gap, baseCell, gapCells)
        print('\t\t%s gap %.4g m: refinement level %d (cell %.4g m, %.1f cells across) within %.4g m'
              % (body.name, gap, level, cell, gap/cell, distance))
        gapRefs.append({'name': 'GAPREF-%s' % (body.name), 'file': 'GAPREF-%s.stl' % (body.name),
                        'distance': distance, 'level': level, 'baseTris': baseGap, 'body': body})

    scales = dict((geomName, float(geomDict[geom]['scale'])) for geomName, geom in geomFiles.items())
    return {'zones': order, 'bodies': bodies, 'startAngle': startAngle, 'step': step, 'freq': freq,
            'nCycles': nCycles, 'movingPids': movingPids, 'gapRefs': gapRefs, 'geomFiles': geomFiles,
            'scales': scales, 'porOwners': porOwners, 'porPids': dict((g, list(pidDict[g])) for g in porOwners),
            'porDeform': porDeform, 'segments': segments, 'meshSetups': meshSetups}


def motionShifts(motionInfo, angle):
    #z shift (file units) of every moved PID at a crank angle: {geometry: {pid: (0, 0, dz)}}
    shifts = {}
    for body in motionInfo['bodies']:
        dz = float(body.shift(angle))/motionInfo['scales'][body.source]
        for pid in body.movedPids:
            shifts.setdefault(body.source, {})[pid] = (0.0, 0.0, dz)
    for geomName, owner in motionInfo['porOwners'].items():
        dz = float(owner.shift(angle))/motionInfo['scales'][geomName]
        shifts[geomName] = dict((pid, (0.0, 0.0, dz)) for pid in motionInfo['porPids'][geomName])
    return shifts


def writeMotionGeometry(motionInfo, mshDir, angle=None, triDir=os.path.join('constant', 'triSurface'), verbose=True):
    #moved bodies at a crank angle (default the start angle): shifted copies in triDir (never write through a
    #symlink), gap STLs. returns the files written
    angle = motionInfo['startAngle'] if angle is None else angle
    os.makedirs(triDir, exist_ok=True)
    written = []
    for geomName, pidShift in motionShifts(motionInfo, angle).items():
        geom = motionInfo['geomFiles'][geomName]
        dest = os.path.join(triDir, geom.replace('.gz', ''))
        for path in (dest, os.path.join(triDir, geom)):
            if os.path.islink(path) or os.path.isfile(path):
                os.remove(path)
        geometryIO.translateSurfacePIDs(os.path.join(mshDir, geom), dest, pidShift)
        written.append(os.path.basename(dest))
        if verbose:
            print('\t\t%s shifted to the start angle: %s' % (geom, ', '.join('%s %+.4g' % (pid, s[2])
                                                                              for pid, s in pidShift.items())))
    #porous zones deforming with the mesh (MOVES_WITH mesh): z of every vertex moved as the mesh moves it
    for geomName, pidTris in motionInfo.get('porDeform', {}).items():
        geom = motionInfo['geomFiles'][geomName]
        dest = os.path.join(triDir, geom.replace('.gz', ''))
        for path in (dest, os.path.join(triDir, geom)):
            if os.path.islink(path) or os.path.isfile(path):
                os.remove(path)
        scale = motionInfo['scales'][geomName]
        deformed = {}
        for pid, t in pidTris.items():
            t = np.array(t, dtype=float)
            t[:, :, 2] += motionUtilities.meshDisplacementZ(t[:, :, 2], motionInfo['zones'], angle)
            deformed[pid] = t/scale
        motionUtilities.writeStlTriangles(dest, deformed)
        written.append(os.path.basename(dest))
        if verbose:
            print('\t\t%s deformed with the mesh to the start angle' % (geom))
    for gapRef in motionInfo['gapRefs']:
        dz = np.array([0.0, 0.0, float(gapRef['body'].shift(angle))])
        motionUtilities.writeStlTriangles(os.path.join(triDir, gapRef['file']),
                                          dict((pid, t + dz) for pid, t in gapRef['baseTris'].items()))
        written.append(gapRef['file'])
    #edge refinement STLs (MESH_GROUPS EDGE_LEVEL) with the moving PIDs at this angle
    written += writeEdgeGeometry(motionInfo.get('edgeRefs', []), triDir, motionShifts(motionInfo, angle))
    return written


def writeRemeshMeshes(motionInfo, mshDir):
    #remeshing: system/remeshSegments (mesh name, start and end angle) and the geometry of each mesh at its start
    #angle in meshes/<name>/constant/triSurface (the other surfaces are links to the case's), seg01 = the case mesh
    listFile = os.path.join('system', 'remeshSegments')
    segments = motionInfo.get('segments')
    if not segments:
        if os.path.isfile(listFile):
            os.remove(listFile)
        if os.path.isdir('meshes'):
            print('\t\tWARNING: meshes/ is left over from an earlier remeshing setup (not used now)')
        return
    print('\n\tWriting the remeshing meshes...')
    names = segmentNames(segments)
    os.makedirs('system', exist_ok=True)
    with open(listFile, 'w', newline='\n') as f:
        f.write('#remeshing (written by internalSetup): mesh, start and end crank angle (deg, mod 360), each mesh is '
                'made at its start angle\n')
        for name, (a0, a1) in zip(names, segments):
            f.write('%s %.9g %.9g\n' % (name, np.mod(a0, 360.0), np.mod(a1, 360.0)))
    #meshes of an earlier setup with more segments
    if os.path.isdir('meshes'):
        for old in sorted(os.listdir('meshes')):
            if old.startswith('seg') and old not in names:
                shutil.rmtree(os.path.join('meshes', old))
                print('\t\tremoved meshes/%s (not in the segments any more)' % (old))
    caseTri = os.path.join('constant', 'triSurface')
    for name, (a0, a1) in zip(names[1:], segments[1:]):
        triDir = os.path.join('meshes', name, 'constant', 'triSurface')
        if os.path.isdir(triDir):
            shutil.rmtree(triDir)
        written = writeMotionGeometry(motionInfo, mshDir, a0, triDir, verbose=False)
        for entry in sorted(os.listdir(caseTri)):
            path = os.path.join(caseTri, entry)
            if entry in written or entry.endswith('.eMesh') or not os.path.isfile(path):
                continue
            os.symlink(os.path.relpath(path, triDir), os.path.join(triDir, entry))
        print('\t\tmeshes/%s: geometry at %.1f deg' % (name, np.mod(a0, 360.0)))


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
        segments = motionInfo.get('segments')
        #remeshing: one dynamicMeshDict per mesh (tables zero at its start angle, stops at its end angle), the run
        #script copies the one of the current mesh to constant/dynamicMeshDict; seg01 is also the case's own
        for old in glob.glob(os.path.join('constant', 'dynamicMeshDict.seg*')) + glob.glob(os.path.join('constant', 'tables', 'seg*')):
            if os.path.isdir(old):
                shutil.rmtree(old)
            else:
                os.remove(old)
        if segments:
            for name, (a0, a1) in zip(segmentNames(segments), segments):
                motionUtilities.writeCrankTables(motionInfo['bodies'], a0, motionInfo['step'], 'constant/tables/%s' % (name))
                writeCrankDict(templateLoc, motionInfo, crank, a0, 'constant/dynamicMeshDict.%s' % (name),
                               '<constant>/tables/%s' % (name), a1)
            shutil.copyfile('constant/dynamicMeshDict.seg01', 'constant/dynamicMeshDict')
        else:
            writeCrankDict(templateLoc, motionInfo, crank, motionInfo['startAngle'], 'constant/dynamicMeshDict',
                           '<constant>/tables')
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
    #remeshing: the dry run covers the first mesh only (its dynamicMeshDict stops at the end of the segment)
    span = motionCheckSpan(motionInfo)
    replaceDict = {'<SOLVER>': 'moveDynamicMesh', '<STARTFROM>': 'startTime', '<ENDTIME>': '%.9g' % (span/(360.0*freq)),
                   '<INIT_DT>': '%.9g' % (1.0/(steps*freq)), '<WRITE_INTERVAL>': '%.9g' % (1.0/(frames*freq)),
                   '<MAX_CO>': '1e6', '<MAX_DT>': '1', '<ADJUST_DT>': 'no',
                   '#include "internalFunctionObjects"': ''}
    for key, val in replaceDict.items():
        text = text.replace(key, val)
    with open('system/controlDictMotion', 'w') as f:
        f.write(text)
    if motionInfo.get('segments'):
        print('\t\t%d motion zones, crank angle tables per mesh (constant/tables/seg*), motion dry run on seg01 only '
              '(%.1f to %.1f deg)' % (len(zones), motionInfo['startAngle'], motionInfo['startAngle'] + span))
    elif crank:
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


def motionCheckSpan(motionInfo):
    #crank angle covered by the motion dry run: one cycle, or the first mesh when remeshing
    segments = motionInfo.get('segments') if motionInfo else None
    return segments[0][1] - segments[0][0] if segments else 360.0


def motionCheckTimes(fullCaseSetupDict, span=360.0):
    #dry run times every 45 deg and at the end (checkMesh -time)
    freq = float(getVal(fullCaseSetupDict, 'GLOBAL_CONTROL', 'FREQ'))
    angles = [45.0*k for k in range(1, 9) if 45.0*k < span - 1e-6] + [span]
    return ','.join('%.9g' % (angle/(360.0*freq)) for angle in angles)


def writeCrankDict(templateLoc, motionInfo, crank, angle, dest, tableDir, stopAngle=None):
    #dynamicMeshDict for the crank drive, tables zero at angle (the meshing angle)
    bodyLines = []
    for body in motionInfo['bodies']:
        patches = ' '.join('"%s_%s"' % (body.source, pid) for pid in body.movedPids)
        #porous zones moving with the body: the gas drag on the matrix (porDrag_<zone>) is part of the body force
        porZones = sorted(name for name, owner in motionInfo.get('porOwners', {}).items() if owner is body)
        porLine = '            porousZones (%s);\n' % (' '.join(porZones)) if porZones else ''
        bodyLines.append('        %s\n        {\n            patches (%s);\n%s            position { type tableFile; '
                         'file "%s/%sCrank"; outOfBounds repeat; }\n        }'
                         % (body.name, patches, porLine, tableDir, body.name))
    zoneLines = ['        (%s %s)' % (zone.name, 'none' if zone.body is None else zone.body.name)
                 for zone in motionInfo['zones']]
    copyTemplate(templateLoc, 'constant/dynamicMeshDictCrank', dest,
                 {'BODIES': '\n'.join(bodyLines), 'ZONES': '\n'.join(zoneLines),
                  'CRANK': crankEntries(crank, angle, stopAngle)})


#----- scripts

def writeMeshSignature():
    #system/meshSignature: md5 of everything the mesh is made from (meshing dicts, topoSetDict motion zones, the
    #geometry as meshed); the run scripts keep an existing mesh only while it matches the copy in .meshDone
    files = [os.path.join('system', name) for name in ['blockMeshDict', 'snappyHexMeshDict', 'snappyRefinementDict',
                                                       'surfaceFeatureExtractDict', 'topoSetDict', 'remeshSegments']]
    files += sorted(glob.glob(os.path.join('meshes', 'seg*', 'snappyHexMeshDict')))
    surfDir = os.path.join('constant', 'triSurface')
    if os.path.isdir(surfDir):
        files += [os.path.join(surfDir, name) for name in sorted(os.listdir(surfDir))
                  if not name.endswith('.eMesh') and os.path.isfile(os.path.join(surfDir, name))]
    lines = []
    for path in files:
        if not os.path.isfile(path):
            continue
        digest = hashlib.md5()
        with open(path, 'rb') as f:
            for block in iter(lambda: f.read(1 << 20), b''):
                digest.update(block)
        lines.append('%s  %s' % (digest.hexdigest(), path.replace('\\', '/')))
    with open(os.path.join('system', 'meshSignature'), 'w', newline='\n') as f:
        f.write('\n'.join(lines) + '\n')


#settings of the cluster (environment, launcher, post job), kept out of GLOBAL_COMPUTE_SETUP (per run settings)
CLUSTER_KEYS = ['ACCOUNT', 'MODULES', 'MPI_RUN', 'PARTITION', 'POST_MODULES', 'POST_MEM', 'POST_ENV', 'POST_CONFIG',
                'PVBATCH', 'GSHEET', 'COPY_BACK']


def clusterSettings(fullCaseSetupDict):
    #GLOBAL_COMPUTE_SETUP and CLUSTER_ENVIRONMENT together; a cluster key still in GLOBAL_COMPUTE_SETUP (case files
    #from before the split) wins, so those cases keep their values
    settings = dict(fullCaseSetupDict.get('CLUSTER_ENVIRONMENT', {}))
    compute = fullCaseSetupDict.get('GLOBAL_COMPUTE_SETUP', {})
    old = [key for key in CLUSTER_KEYS if key in compute]
    if old and 'CLUSTER_ENVIRONMENT' in fullCaseSetupDict:
        print('\t\tNOTE: %s set in GLOBAL_COMPUTE_SETUP (used); they belong in CLUSTER_ENVIRONMENT now' % (' '.join(old)))
    settings.update(compute)
    return settings


def makeScriptsInternal(templateLoc, fullCaseSetupDict, execDir, motionZones=None, fanNames=None, motionInfo=None,
                        warmStart=None):
    #writing runScript, postScript, motionCheckScript and submitScript from clusterDict
    print('\n\tWriting scripts...')
    os.makedirs('system', exist_ok=True)
    writeMeshSignature()
    with open('%s/defaultCluster/slurm/clusterDict' % (templateLoc)) as f:
        clusterDict = json.load(f)
    solver = getVal(fullCaseSetupDict, 'GLOBAL_SIM_CONTROL', 'SOLVER')
    motionSteps = ['topoSet', 'reorderFaceZones', 'pointDisplacement', 'motionCheck', 'meshImages']
    skipSteps = list(motionSteps) if not motionZones else []
    #the library (crank drive, slidingWallVelocity) is compiled in the job for moving meshes (wmake skips it when up to date)
    if not motionZones:
        skipSteps.append('compileLib')
    if not fanNames:
        skipSteps.append('createBaffles')
    #remeshing: the segment meshes are made with the case mesh, the solve loops over them (remeshSolve)
    if motionInfo and motionInfo.get('segments'):
        skipSteps.append('solve')
    else:
        skipSteps += ['remeshMesh', 'remeshSolve']
    #remeshMassFix also does the warm start (gas T profile, p for the charge mass) before the first solve
    if not (motionInfo and motionInfo.get('segments')) and not warmStart:
        skipSteps.append('compileMassFix')
    if not warmStart:
        skipSteps.append('warmStart')
    export = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP'].get('EXPORT', ['ensight'])[0].lower()
    if export not in ('ensight', 'none'):
        sys.exit('ERROR! GLOBAL_COMPUTE_SETUP EXPORT must be ensight or none!')
    if export == 'none':
        skipSteps.append('export')
    #post job: summary, gsheet, images and videos (images need the EnSight export)
    cluster = clusterSettings(fullCaseSetupDict)
    gsheet = cluster.get('GSHEET', ['True'])[0].lower()
    if gsheet not in ('true', 'false'):
        sys.exit('ERROR! CLUSTER_ENVIRONMENT GSHEET must be True or False!')
    if gsheet == 'false':
        skipSteps.append('gsheet')
    #copy the run data back to Google Drive at the end of the post job (copyBackInternal.py)
    copyBack = cluster.get('COPY_BACK', ['True'])[0].lower()
    if copyBack not in ('true', 'false'):
        sys.exit('ERROR! CLUSTER_ENVIRONMENT COPY_BACK must be True or False!')
    if copyBack == 'false':
        skipSteps.append('copyBack')
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
    checkTimes = motionCheckTimes(fullCaseSetupDict, motionCheckSpan(motionInfo)) if motionZones else ''
    replaceDict = {'<SET_DECOMP>': setDecomp, '<CHECK_TIMES>': checkTimes,
                   '<ZERO_TEMPLATES>': os.path.join(execDir, 'zeroTemplates'), '<EXEC_DIR>': execDir,
                   '<SOLVER>': solver, '<MOTION_ZONES>': ' '.join(motionZones or []),
                   '<MASS_FIX>': 'yes' if getFlow(fullCaseSetupDict) == 'compressible' else 'no',
                   '<REMESH_DT>': ''.join(fullCaseSetupDict.get('MESH_MOTION', {}).get('REMESH_DT', [''])),
                   '<WARM_ZONE>': warmStart['zone'] if warmStart else '',
                   '<WARM_T>': '%.9g %.9g' % tuple(warmStart['T']) if warmStart else ''}
    #slurm header from GLOBAL_COMPUTE_SETUP and CLUSTER_ENVIRONMENT, job name defaults to the case directory name
    compute = cluster
    for key in compute:
        replaceDict['<%s>' % key] = ' '.join(compute[key])
    #parallel launcher inside a Slurm job (srun, or e.g. mpirun where Open MPI cannot start under srun)
    if replaceDict.get('<MPI_RUN>', '').strip() == '':
        replaceDict['<MPI_RUN>'] = 'srun'
    if replaceDict.get('<JOB_NAME>', 'auto') == 'auto':
        replaceDict['<JOB_NAME>'] = os.path.basename(os.getcwd())
    #clusters without an account, with a partition / fixed node, or without environment modules (OpenFOAM from the
    #submitting shell, e.g. MPI_RUN = mpirun .../foamExec): ACCOUNT, MODULES, POST_MODULES blank, PARTITION, NODELIST
    #POST_MEM blank: no memory request for the post job (nodes without RealMemory in slurm.conf reject any --mem)
    extraHeader = ['#SBATCH --%s=%s' % (option, ' '.join(compute.get(key, ['']))) for key, option in
                   (('PARTITION', 'partition'), ('NODELIST', 'nodelist')) if ''.join(compute.get(key, [''])).strip()]
    #EXCLUSIVE False: share the node (small meshes on few cores), --mem=0 (all node memory) replaced by MEM_PER_CPU
    #(blank: the cluster default per core)
    exclusive = ''.join(compute.get('EXCLUSIVE', ['True'])).strip().lower() or 'true'
    if exclusive not in ('true', 'false'):
        sys.exit('ERROR! GLOBAL_COMPUTE_SETUP EXCLUSIVE must be True or False!')
    memPerCpu = ''.join(compute.get('MEM_PER_CPU', [''])).strip()
    shareLines = {'#SBATCH --exclusive': [], '#SBATCH --mem=0': ['#SBATCH --mem-per-cpu=%s' % (memPerCpu)] if memPerCpu else []}
    if exclusive == 'false':
        print('\t\tshared nodes (EXCLUSIVE False), memory %s' % ('%s per core' % (memPerCpu) if memPerCpu else 'cluster default per core'))
    #CPUS_PER_TASK > 1: each MPI rank gets that many cores, so fewer ranks are spread evenly over the node (one per
    #L3 cache group instead of all on the first cores; OpenFOAM is memory-bandwidth bound). srun does not take
    #--cpus-per-task from the job (Slurm 22.05+): SRUN_CPUS_PER_TASK. OpenFOAM runs one thread per rank
    cpusPerTask = ''.join(compute.get('CPUS_PER_TASK', ['1'])).strip() or '1'
    if not cpusPerTask.isdigit() or int(cpusPerTask) < 1:
        sys.exit('ERROR! GLOBAL_COMPUTE_SETUP CPUS_PER_TASK must be a whole number >= 1!')
    if int(cpusPerTask) > 1:
        print('\t\t%s ranks x %s cores each = %d cores per job' % (' '.join(compute['NUM_PROCS']), cpusPerTask,
                                                                  int(compute['NUM_PROCS'][0])*int(cpusPerTask)))
    for section, scriptName in scripts:
        lines = [clusterDict[section][key] for key in clusterDict[section] if key not in skipSteps]
        script = '\n'.join(lines)
        for key, val in replaceDict.items():
            script = script.replace(key, val)
        script = script.replace('module --force purge; module load ;', '').replace('module load ;', '')
        #an if with nothing left in it (postScript: only the module line) is a bash syntax error
        script = script.replace('if [ -n "$SLURM_JOB_ID" ]; then  fi;', '').replace('if [ -n "$SLURM_JOB_ID" ]; then fi;', '')
        script = '\n'.join(line for line in script.split('\n') if line.strip() not in ('#SBATCH --account=', '#SBATCH --mem='))
        if exclusive == 'false':
            script = '\n'.join(new for line in script.split('\n') for new in shareLines.get(line.strip(), [line]))
        if int(cpusPerTask) > 1 and section in ('run', 'motionCheck'):
            script = '\n'.join(new for line in script.split('\n') for new in
                               ([line, '#SBATCH --cpus-per-task=%s' % (cpusPerTask)] if line.startswith('#SBATCH --ntasks=') else [line]))
            script = script.replace('export UCX_VFS_ENABLE=n;', 'export UCX_VFS_ENABLE=n; export SRUN_CPUS_PER_TASK=%s OMP_NUM_THREADS=1;'
                                    % (cpusPerTask))
        if extraHeader:
            out, added = [], False
            for line in script.split('\n'):
                out.append(line)
                if not added and line.startswith('#SBATCH --job-name='):
                    out += extraHeader
                    added = True
            script = '\n'.join(out)
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
    #transient without cycles (END_TIME): post averages over the end of the run instead of whole cycles
    if not isSteady(fullCaseSetupDict) and getEndTime(fullCaseSetupDict) is not None:
        fullCaseSetupDict['SETUP_INFO']['CYCLIC'] = ['no']
    if motionInfo:
        fullCaseSetupDict['SETUP_INFO']['START_ANGLE'] = ['%.6g' % (motionInfo['startAngle'])]
        #remeshing: start angle of each mesh (deg, mod 360)
        if motionInfo.get('segments'):
            fullCaseSetupDict['SETUP_INFO']['REMESH_ANGLES'] = ['%.6g' % (np.mod(a0, 360.0))
                                                               for a0, a1 in motionInfo['segments']]
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
