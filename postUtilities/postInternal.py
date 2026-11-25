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

#post-processing of an internal (internalSetup) case: summary.csv, pV.png and cycleConvergence.png in the case
#(key,value rows), the database row pushed to the database gsheet by gsheetInternal.py
#usage (in the case directory): python3 postInternal.py [--avgCycles N] [--gsheet]
#
#fixed keys: JOB NAME DESCRIPTION RUN_DATE COMPLETION_DATE CELL_COUNT MESH_TIME SOLVE_TIME POST_TIME
#AVG_POWER AVG_ENERGY, then setup values, run status and results (each becomes a sheet column)
#times in hours, dates YYYY-MM-DD HH:MM
#AVG_ENERGY: net work per cycle done by the gas on the moving bodies (sum over bodies of the integral of
#F_z v dt over a cycle, J), AVG_POWER = AVG_ENERGY*FREQ (W). positive = power out of the gas (engine output),
#a displacer alone gives a negative value (the work needed to drive it, flow loss)

import argparse
import configparser
import csv
import datetime
import glob
import os
import re
import sys
import time

import numpy as np

FIXED_KEYS = ['JOB', 'NAME', 'DESCRIPTION', 'COMMENTS', 'RUN_DATE', 'COMPLETION_DATE', 'CELL_COUNT',
              'MESH_TIME', 'SOLVE_TIME', 'POST_TIME', 'AVG_POWER', 'AVG_ENERGY']
MESH_LOGS = ['log.blockMesh', 'log.decomposePar', 'log.snappyHexMesh', 'log.topoSet', 'log.checkMesh']



#----- case setup

def readCaseSetup(caseDir):
    path = os.path.join(caseDir, 'fullCaseSetupDict')
    if not os.path.isfile(path):
        sys.exit('ERROR! %s not found: run internalSetup.py in the case first' % (path))
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(path)
    if config.get('SETUP_INFO', 'TOOL', fallback='') != 'internal':
        sys.exit('ERROR! %s is not an internalSetup case (no [SETUP_INFO] TOOL = internal)' % (caseDir))
    return config


def getSetting(config, section, key, default=''):
    return ' '.join(config.get(section, key, fallback=default).split())


def getBodies(config):
    #moving bodies: name, motion table and the patches of their force function object
    return [section[len('MOVE-'):] for section in config.sections() if section.startswith('MOVE-')]


def getGroupSettings(config):
    #thermal setting per BC group: {group: 'fixed 600'} from the WALL/MOVE group lines (PIDS .. T ..)
    groups = {}
    for section in config.sections():
        if not (section.startswith('WALL-') or section.startswith('MOVE-')):
            continue
        names = getSetting(config, section, 'BC_GROUPS').split()
        for name in names:
            tokens = getSetting(config, section, name).split()
            if 'T' in tokens:
                rest = tokens[tokens.index('T') + 1:]
                if 'U' in rest:
                    rest = rest[:rest.index('U')]
                groups[name] = ' '.join(rest)
        groups['%s_default' % (section)] = getSetting(config, section, 'T_DEFAULT')
    return groups


def gasConstants(caseDir):
    #specific gas constant R and Cp of the perfect gas from constant/thermophysicalProperties, (None, None) if missing
    path = os.path.join(caseDir, 'constant', 'thermophysicalProperties')
    if not os.path.isfile(path):
        return None, None
    with open(path) as f:
        text = f.read()
    molWeight = re.search(r'\bmolWeight\s+([-+0-9.eE]+)\s*;', text)
    cp = re.search(r'\bCp\s+([-+0-9.eE]+)\s*;', text)
    if not molWeight or not cp:
        return None, None
    return 8314.46/float(molWeight.group(1)), float(cp.group(1))


def isothermalT(config):
    #reference T of the isothermal limit: the fixed wall T when every fixed-T group has the same one, else INIT_T
    fixed = set()
    for setting in getGroupSettings(config).values():
        tokens = setting.split()
        if len(tokens) == 2 and tokens[0] == 'fixed':
            fixed.add(float(tokens[1]))
    if len(fixed) == 1:
        return fixed.pop(), 'fixed wall T'
    return float(getSetting(config, 'GLOBAL_MATERIAL', 'INIT_T', '300') or 300), 'INIT_T'


def idealGasSummary(caseDir, config, t, mass, volume, p, gasT, t0, t1):
    #ideal-gas limits over the averaging cycles: polytropic exponent n of p V^n = const fitted to the loop (1 =
    #isothermal, gamma = adiabatic), p over the isothermal p = m R T_ref / V, gas T swing, loop work over p V
    sel = window(t, t0, t1)
    out = {}
    fit = np.polyfit(np.log(volume[sel]), np.log(p[sel]), 1)
    out['POLY_N'] = float(-fit[0])
    if gasT is not None:
        out['T_SWING'] = float(gasT[sel].max() - gasT[sel].min())
    pV = timeMean(t, p*volume, t0, t1)
    if pV > 0:
        out['PV_WORK_REL'] = float(integrate(volume[sel], p[sel])/pV)
    R, cp = gasConstants(caseDir)
    if R is None:
        return out
    tRef, source = isothermalT(config)
    dev = p[sel]*volume[sel]/(mass[sel]*R*tRef) - 1.0
    out['GAMMA'] = float(cp/(cp - R))
    out['T_ISO_REF'] = tRef
    out['P_ISO_DEV_MAX'] = float(dev[np.argmax(np.abs(dev))])
    print('\tideal gas: n = %.4f (isothermal 1, adiabatic %.4f), p/p_iso - 1 up to %+.3g %% (T_ref %g K, %s), '
          'gas T swing %.3g K' % (out['POLY_N'], out['GAMMA'], 100*out['P_ISO_DEV_MAX'], tRef, source,
                                  out.get('T_SWING', float('nan'))))
    return out


#----- logs

def clockTime(path):
    #last ClockTime of an OpenFOAM log in s, None if there is none
    last = None
    with open(path, errors='replace') as f:
        for line in f:
            m = re.search(r'ClockTime\s*=\s*([0-9.eE+-]+)', line)
            if m:
                last = float(m.group(1))
    return last


def logStart(path):
    #Date / Time from the OpenFOAM log header
    date, clock = None, None
    with open(path, errors='replace') as f:
        for k, line in enumerate(f):
            m = re.match(r'^Date\s*:\s*(.+)$', line)
            if m:
                date = m.group(1).strip()
            m = re.match(r'^Time\s*:\s*(\d+:\d+:\d+)\s*$', line)
            if m:
                clock = m.group(1)
            if (date and clock) or k > 60:
                break
    if not (date and clock):
        return None
    try:
        return datetime.datetime.strptime('%s %s' % (' '.join(date.split()), clock), '%b %d %Y %H:%M:%S')
    except ValueError:
        return None


def logDuration(path):
    #run time of a log in s: last ClockTime, or (meshing utilities print none) header start to last write
    seconds = clockTime(path)
    if seconds is None:
        startStamp = logStart(path)
        if startStamp is not None:
            seconds = max(0.0, os.path.getmtime(path) - time.mktime(startStamp.timetuple()))
    return seconds


def logEnded(path):
    #the solver writes End as the last word of a finished run
    with open(path, errors='replace') as f:
        tail = f.read()[-2000:]
    return re.search(r'^End\s*$', tail, re.M) is not None


def solverLogs(caseDir, solver):
    #log.<solver> plus the logs of earlier restarts (log.<solver>.1, .2, ...), oldest first
    logs = glob.glob(os.path.join(caseDir, 'log.%s.*' % (solver)))
    logs = sorted((p for p in logs if p.rsplit('.', 1)[-1].isdigit()), key=lambda p: int(p.rsplit('.', 1)[-1]))
    current = os.path.join(caseDir, 'log.%s' % (solver))
    if os.path.isfile(current):
        logs.append(current)
    return logs


def solverQuality(logs, t0, t1, pCap):
    #per time step in t0 < t <= t1 from the solver logs: the Courant number (the line printed before "Time = ", for
    #that step) and the iterations of the last p solve of the step (the final one). returns (Co max, share of steps
    #with Co > 1, share of steps whose final p solve stopped at the iteration cap, number of steps)
    coPending = None
    steps = []
    for path in logs:
        with open(path, errors='replace') as f:
            for line in f:
                if line.startswith('Courant Number mean'):
                    m = re.match(r'^Courant Number mean:\s*(\S+)\s+max:\s*(\S+)', line)
                    if m:
                        coPending = float(m.group(2))
                elif line.startswith('Time = '):
                    try:
                        t = float(line.split('=', 1)[1])
                    except ValueError:
                        continue
                    steps.append([t, coPending, None])
                    coPending = None
                elif 'Solving for p' in line and steps:
                    m = re.search(r'Solving for (p|p_rgh),.*No Iterations (\d+)', line)
                    if m:
                        steps[-1][2] = int(m.group(2))
    steps = [s for s in steps if t0 < s[0] <= t1 + 1e-12]
    if not steps:
        return None, None, None, 0
    co = [s[1] for s in steps if s[1] is not None]
    pIter = [s[2] for s in steps if s[2] is not None]
    coMax = max(co) if co else None
    above = sum(1 for c in co if c > 1.0)/len(co) if co else None
    capped = sum(1 for n in pIter if n >= pCap)/len(pIter) if pIter else None
    return coMax, above, capped, len(steps)


def cellCount(caseDir):
    path = os.path.join(caseDir, 'log.checkMesh')
    if not os.path.isfile(path):
        return ''
    with open(path, errors='replace') as f:
        for line in f:
            m = re.match(r'^\s*cells:\s*(\d+)', line)
            if m:
                return int(m.group(1))
    return ''


def formatDate(stamp):
    return stamp.strftime('%Y-%m-%d %H:%M') if stamp else ''


def hours(seconds):
    return round(seconds/3600.0, 4) if seconds is not None else ''


#----- function object output

def toFloat(value):
    try:
        return float(value)
    except ValueError:
        return float('nan')


def readDat(foDir, fileName, textColumns=False):
    #all <foDir>/<startTime>/<fileName>*.dat (restarts start a new directory or file), merged in time order;
    #later files replace overlapping times. returns (column names, array)
    #textColumns: text entries (solverInfo solver names, true/false) become nan instead of dropping the row
    files = []
    for timeDir in glob.glob(os.path.join(foDir, '*')):
        try:
            start = float(os.path.basename(timeDir))
        except ValueError:
            continue
        for path in glob.glob(os.path.join(timeDir, '%s*.dat' % (fileName))):
            files.append((start, os.path.getmtime(path), path))
    columns, data = None, []
    for start, mtime, path in sorted(files):
        names, rows = None, []
        with open(path, errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    if re.match(r'^#\s*Time\b', line):
                        names = line.lstrip('#').split()
                    continue
                #vector entries (x y z) become three columns
                values = line.replace('(', ' ').replace(')', ' ').split()
                if values:
                    try:
                        rows.append([float(v) for v in values])
                    except ValueError:
                        if textColumns:
                            rows.append([toFloat(v) for v in values])
        if not rows or names is None:
            continue
        rows = np.array(rows)
        columns = names
        if data:
            data = [d[d[:, 0] < rows[0, 0]] for d in data]
        data.append(rows)
    if not data:
        return None, None
    return columns, np.vstack(data)


def readWallHeatFlux(foDir):
    #wallHeatFlux rows: time patch min max integral -> {patch: (times, integral)}
    out = {}
    for timeDir in sorted(glob.glob(os.path.join(foDir, '*')), key=lambda p: os.path.basename(p)):
        for path in glob.glob(os.path.join(timeDir, 'wallHeatFlux*.dat')):
            with open(path, errors='replace') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    values = line.split()
                    if len(values) >= 5:
                        try:
                            t, q = float(values[0]), float(values[4])
                        except ValueError:
                            continue
                        out.setdefault(values[1], []).append((t, q))
    return {patch: np.array(sorted(dict(rows).items())) for patch, rows in out.items()}


def readTable(path):
    #scalarVectorTable (t (dx dy dz)) -> times, z displacement
    with open(path) as f:
        text = f.read()
    rows = re.findall(r'\(\s*([-0-9.eE+]+)\s+\(\s*([-0-9.eE+]+)\s+([-0-9.eE+]+)\s+([-0-9.eE+]+)\s*\)\s*\)', text)
    rows = np.array(rows, dtype=float)
    return rows[:, 0], rows[:, 3]


def readCrankTable(path):
    #constant/tables/<body>Crank (theta dz), crank drive and remeshing
    with open(path) as f:
        text = f.read()
    rows = np.array(re.findall(r'\(\s*([-0-9.eE+]+)\s+([-0-9.eE+]+)\s*\)', text), dtype=float)
    return rows[:, 0], rows[:, 1]


def tableVelocity(times, dz, t):
    #the solver interpolates the table linearly in time, so the velocity is constant between table points
    slopes = np.diff(dz)/np.diff(times)
    k = np.clip(np.searchsorted(times, t, side='right') - 1, 0, len(slopes) - 1)
    return slopes[k]


def window(times, t0, t1):
    return (times >= t0 - 1e-12) & (times <= t1 + 1e-12)


def clipped(x, y, x0, x1):
    #samples of y(x) between x0 and x1, end points interpolated
    inside = (x > x0) & (x < x1)
    return np.r_[x0, x[inside], x1], np.r_[np.interp(x0, x, y), y[inside], np.interp(x1, x, y)]


def integrate(t, y):
    #trapezoidal rule (np.trapz is deprecated in newer numpy)
    t, y = np.asarray(t, float), np.asarray(y, float)
    return float(np.sum(0.5*(y[1:] + y[:-1])*np.diff(t))) if len(t) > 1 else 0.0


def readCrank(caseDir, post):
    #reading the crank drive log, theta (deg) accumulated from theta0 (the start angle)
    #returns (columns, array, theta0) or Nones
    names, crank = readDat(os.path.join(post, 'crank'), 'crank')
    if crank is None:
        return None, None, None
    theta0 = None
    #remeshing: constant/dynamicMeshDict is the current segment's (theta0 = its start angle), seg01 has the run's
    path = os.path.join(caseDir, 'constant', 'dynamicMeshDict.seg01')
    if not os.path.isfile(path):
        path = os.path.join(caseDir, 'constant', 'dynamicMeshDict')
    if os.path.isfile(path):
        with open(path) as f:
            match = re.search(r'^\s*theta0\s+([-0-9.eE+]+)\s*;', f.read(), re.M)
        if match:
            theta0 = float(match.group(1))
    if theta0 is None:
        theta0 = float(crank[0, 1])
    #the run starts at t = 0 at theta0 (the log starts after the first step)
    if crank[0, 1] > theta0:
        crank = np.vstack([np.r_[0.0, theta0, crank[0, 2:]], crank])
    #a single repeated angle is a segment end, not a stop: runs before the crank solver requested its stop 1.5 steps
    #ahead could take one extra step at the stop angle with the crank at rest (no work, no motion)
    d = np.diff(crank[:, 1])
    repeat = [i + 1 for i in range(len(d)) if d[i] == 0 and (i == 0 or d[i - 1] > 0) and (i + 1 == len(d) or d[i + 1] > 0)]
    if repeat:
        crank = np.delete(crank, repeat, axis=0)
    if np.any(np.diff(crank[:, 1]) <= 0):
        print('\tWARNING: the crank stopped or turned backwards: cycles counted up to the first stop')
        stop = int(np.argmax(np.diff(crank[:, 1]) <= 0))
        crank = crank[:stop + 1]
    return names, crank, theta0


def patchArea(foDir):
    #patch area from the '# Area' header line of a surfaceFieldValue file, None if missing
    for path in sorted(glob.glob(os.path.join(foDir, '*', '*.dat'))):
        with open(path, errors='replace') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                match = re.match(r'^#\s*Area\s*:\s*([-0-9.eE+]+)', line)
                if match:
                    return float(match.group(1))
    return None


def flowGroups(config):
    #inlet and outlet groups: {group: (kind, [patches])} from BC_GROUP_KINDS / BC_GROUPS_RESOLVED
    kinds = dict(config.items('BC_GROUP_KINDS')) if config.has_section('BC_GROUP_KINDS') else {}
    patches = dict(config.items('BC_GROUPS_RESOLVED')) if config.has_section('BC_GROUPS_RESOLVED') else {}
    return dict((group, (kind.split()[0], patches.get(group, '').split())) for group, kind in kinds.items()
                if kind.split()[0] in ('inlet', 'outlet'))


def flowUnits(config, compressible):
    #(flow key, flux scale, p scale, p offset): incompressible flux and kinematic p times RHO_REF, p plus
    #CHARGE_PRESSURE, so both flow types give kg/s and Pa (raw m3/s and kinematic p without RHO_REF)
    if compressible:
        return 'MDOT', 1.0, 1.0, 0.0
    rhoRef = getSetting(config, 'SETUP_INFO', 'RHO_REF', '')
    if not rhoRef:
        return 'VOLUME_FLOW', 1.0, 1.0, 0.0
    return 'MDOT', float(rhoRef), float(rhoRef), float(getSetting(config, 'GLOBAL_MATERIAL', 'CHARGE_PRESSURE', '0'))


def flowSummary(caseDir, config, t0, t1, compressible):
    #per inlet/outlet group over t0..t1: MDOT (kg/s, + into the domain), area-averaged P (Pa) and T; then the mass
    #balance and the pressure drop (inlets minus outlets). incompressible: flux and kinematic p times RHO_REF
    #(gas density at the charge state), p plus CHARGE_PRESSURE; without RHO_REF the raw VOLUME_FLOW and kinematic p
    post = os.path.join(caseDir, 'postProcessing')
    flowKey, fluxScale, pScale, pOffset = flowUnits(config, compressible)
    results, total, inflow, pIn, pOut = {}, 0.0, 0.0, [], []
    for group, (kind, patches) in flowGroups(config).items():
        flux, found, pSum, tSum, weight = 0.0, False, 0.0, 0.0, 0.0
        for patch in patches:
            names, rows = readDat(os.path.join(post, 'flow_%s' % (patch)), 'surfaceFieldValue')
            if rows is not None:
                flux -= fluxScale*timeMean(rows[:, 0], rows[:, 1], t0, t1)
                found = True
            aNames, avg = readDat(os.path.join(post, 'avg_%s' % (patch)), 'surfaceFieldValue')
            if avg is not None:
                area = patchArea(os.path.join(post, 'avg_%s' % (patch))) or 1.0
                pSum += area*(pOffset + pScale*timeMean(avg[:, 0], avg[:, aNames.index('areaAverage(p)')], t0, t1))
                if 'areaAverage(T)' in aNames:
                    tSum += area*timeMean(avg[:, 0], avg[:, aNames.index('areaAverage(T)')], t0, t1)
                weight += area
        if not found:
            continue
        results['%s_%s' % (flowKey, group)] = flux
        total += flux
        inflow += max(flux, 0.0)
        if weight > 0:
            results['P_%s' % (group)] = pSum/weight
            (pIn if kind == 'inlet' else pOut).append(pSum/weight)
            if compressible and tSum:
                results['T_%s' % (group)] = tSum/weight
    if results:
        results['%s_IMBALANCE' % (flowKey)] = total/inflow if inflow > 0 else ''
        if pIn and pOut:
            results['P_DROP'] = float(np.mean(pIn) - np.mean(pOut))
    return results


def heatSummary(config, post, t0, t1):
    #mean heat rate per BC group over t0..t1 (W, + into the gas), from wallHeatFlux and BC_GROUPS_RESOLVED
    heat = readWallHeatFlux(os.path.join(post, 'wallHeatFlux'))
    groupPatches = dict(config.items('BC_GROUPS_RESOLVED')) if config.has_section('BC_GROUPS_RESOLVED') else {}
    results = {}
    for group, patches in groupPatches.items():
        total, found = 0.0, False
        for patch in patches.split():
            if patch in heat:
                rows = heat[patch]
                sel = window(rows[:, 0], t0, t1)
                if sel.sum() > 1:
                    total += integrate(rows[sel, 0], rows[sel, 1])/(t1 - t0)
                    found = True
        if found:
            results['HEAT_%s' % (group)] = total
    return results


def porousDragZ(config, post, body, times):
    #gas drag on the porous zones moving with the body (MOVES_WITH), z component at the given times (N), None if none
    total = None
    for section in config.sections():
        if not section.startswith('POR') or getSetting(config, section, 'MOVES_WITH') != body:
            continue
        names, rows = readDat(os.path.join(post, 'porDrag_%s' % (section)), 'volFieldValue')
        if rows is None or rows.shape[1] < 4:
            continue
        fz = np.interp(times, rows[:, 0], rows[:, 3])
        total = fz if total is None else total + fz
    return total


def porousSummary(post, t0, t1):
    #porous zones over t0..t1: HEAT_<zone> (W from the matrix into the gas), TS_MEAN_<zone> (matrix T), T_MEAN_<zone>
    results = {}
    #porHeatToGas_<zone> (porHeat_<zone> in older cases)
    for foDir in sorted(glob.glob(os.path.join(post, 'porHeat*_*'))):
        names, rows = readDat(foDir, 'volFieldValue')
        if rows is not None:
            zoneName = os.path.basename(foDir).split('_', 1)[1]
            results['HEAT_%s' % (zoneName)] = timeMean(rows[:, 0], rows[:, 1], t0, t1)
    #porDragHeat_<zone>: drag work put into the gas (W, DISSIPATION yes)
    for foDir in sorted(glob.glob(os.path.join(post, 'porDragHeat_*'))):
        names, rows = readDat(foDir, 'volFieldValue')
        if rows is not None:
            zoneName = os.path.basename(foDir).split('_', 1)[1]
            results['DRAG_HEAT_%s' % (zoneName)] = timeMean(rows[:, 0], rows[:, 1], t0, t1)
    for foDir in sorted(glob.glob(os.path.join(post, 'por_*'))):
        zone = os.path.basename(foDir)[len('por_'):]
        names, rows = readDat(foDir, 'volFieldValue')
        if rows is None:
            continue
        for k, name in enumerate(names[1:], start=1):
            field = re.sub(r'^volAverage\((.*)\)$', r'\1', name)
            key = 'TS_MEAN_%s' % (zone) if field.startswith('Ts_') else '%s_MEAN_%s' % (field.upper(), zone)
            results[key] = timeMean(rows[:, 0], rows[:, k], t0, t1)
    return results


#----- summary

def buildSummary(caseDir, avgCycles=1, avgIters=50, avgTime=None):
    start = time.time()
    caseDir = os.path.abspath(caseDir)
    config = readCaseSetup(caseDir)
    summary = dict((key, '') for key in FIXED_KEYS)
    summary.pop('COMMENTS')

    #identity: case directory inside <job>/CASES/<case>
    jobDir = os.path.dirname(os.path.dirname(caseDir))
    summary['JOB'] = os.path.basename(jobDir) if os.path.basename(os.path.dirname(caseDir)) == 'CASES' else ''
    summary['NAME'] = os.path.basename(caseDir)
    summary['DESCRIPTION'] = getSetting(config, 'TITLES', 'CASEDESCRP')

    #setup values
    solver = getSetting(config, 'GLOBAL_SIM_CONTROL', 'SOLVER')
    freq = float(getSetting(config, 'GLOBAL_CONTROL', 'FREQ') or 0)
    compressible = getSetting(config, 'SETUP_INFO', 'FLOW', 'compressible') == 'compressible'
    for key, section, item in [('SOLVER', 'GLOBAL_SIM_CONTROL', 'SOLVER'), ('TURB_MODEL', 'GLOBAL_SIM_CONTROL', 'TURB_MODEL'),
                               ('GAS', 'GLOBAL_MATERIAL', 'GAS'), ('CHARGE_PRESSURE', 'GLOBAL_MATERIAL', 'CHARGE_PRESSURE'),
                               ('INIT_T', 'GLOBAL_MATERIAL', 'INIT_T'), ('FREQ', 'GLOBAL_CONTROL', 'FREQ'),
                               ('N_CYCLES', 'GLOBAL_CONTROL', 'N_CYCLES'), ('BASE_CELL', 'GLOBAL_REFINEMENT', 'BASE_CELL'),
                               ('NUM_PROCS', 'GLOBAL_COMPUTE_SETUP', 'NUM_PROCS'), ('VERSION', 'SETUP_INFO', 'VERSION')]:
        summary[key] = getSetting(config, section, item)
    summary['FLOW'] = 'compressible' if compressible else 'incompressible'
    for group, setting in getGroupSettings(config).items():
        summary['BC_%s' % (group)] = setting
    bodies = getBodies(config)
    for body in bodies:
        summary['TDC_CLEARANCE_%s' % (body)] = getSetting(config, 'MOVE-%s' % (body), 'TDC_CLEARANCE')

    #mesh and run times from the logs
    summary['CELL_COUNT'] = cellCount(caseDir)
    meshTimes = [logDuration(os.path.join(caseDir, log)) for log in MESH_LOGS if os.path.isfile(os.path.join(caseDir, log))]
    meshTimes = [t for t in meshTimes if t is not None]
    summary['MESH_TIME'] = hours(sum(meshTimes)) if meshTimes else ''
    logs = solverLogs(caseDir, solver)
    if logs:
        solveTimes = [t for t in (logDuration(p) for p in logs) if t is not None]
        summary['SOLVE_TIME'] = hours(sum(solveTimes)) if solveTimes else ''
        starts = [s for s in (logStart(p) for p in logs) if s is not None]
        summary['RUN_DATE'] = formatDate(min(starts)) if starts else ''
        if logEnded(logs[-1]):
            summary['COMPLETION_DATE'] = formatDate(datetime.datetime.fromtimestamp(os.path.getmtime(logs[-1])))
    exportTime = logDuration(os.path.join(caseDir, 'log.foamToEnsight')) if os.path.isfile(os.path.join(caseDir, 'log.foamToEnsight')) else None

    #run status: complete cycles from t = 0 (the motion tables start at 0), by crank angle for DRIVE crank
    post = os.path.join(caseDir, 'postProcessing')
    names, gasAvg = readDat(os.path.join(post, 'gasAvg'), 'volFieldValue')
    endTime = float(gasAvg[-1, 0]) if gasAvg is not None else 0.0

    #steady solvers: results averaged over the last avgIters iterations; transient runs without cycles (END_TIME):
    #over the last avgTime seconds (default the last 10 % of the run)
    steady = getSetting(config, 'SETUP_INFO', 'STEADY', 'no') == 'yes'
    if steady or getSetting(config, 'SETUP_INFO', 'CYCLIC', 'yes') == 'no':
        for key in ['FREQ', 'N_CYCLES']:
            summary.pop(key, None)
        if steady:
            summary['N_ITER'] = getSetting(config, 'GLOBAL_CONTROL', 'N_ITER')
            summary['END_ITER'] = int(round(endTime))
        else:
            summary['END_TIME_SET'] = getSetting(config, 'GLOBAL_CONTROL', 'END_TIME')
            summary['END_TIME'] = round(endTime, 9)
        summary['FINISHED'] = 'yes' if summary['COMPLETION_DATE'] else 'no'
        if gasAvg is None:
            print('\tWARNING: no gasAvg output yet, results left empty')
            summary['POST_TIME'] = hours(time.time() - start + (exportTime or 0.0))
            return summary
        if steady:
            t0, t1 = max(float(gasAvg[0, 0]), endTime - avgIters), endTime
            summary['AVG_ITERS'] = int(round(t1 - t0))
            print('\taveraging over iterations %d to %d' % (int(t0), int(t1)))
        else:
            span = avgTime if avgTime else 0.1*endTime
            t0, t1 = max(float(gasAvg[0, 0]), endTime - span), endTime
            summary['AVG_TIME'] = t1 - t0
            print('\taveraging over t = %.6g to %.6g s' % (t0, t1))
        for k, name in enumerate(names[1:], start=1):
            field = re.sub(r'^volAverage\((.*)\)$', r'\1', name).upper()
            summary['%s_MEAN' % (field)] = timeMean(gasAvg[:, 0], gasAvg[:, k], t0, t1)
        summary.update(flowSummary(caseDir, config, t0, t1, compressible))
        summary.update(porousSummary(post, t0, t1))
        if compressible:
            summary.update(heatSummary(config, post, t0, t1))
            mNames, mass = readDat(os.path.join(post, 'gasMass'), 'volFieldValue')
            if mass is not None and not steady and not flowGroups(config):
                summary['MASS_DRIFT'] = float(mass[-1, 1]/mass[0, 1] - 1.0)
        summary['POST_TIME'] = hours(time.time() - start + (exportTime or 0.0))
        return summary

    crankNames, crank, theta0 = readCrank(caseDir, post)
    summary['DRIVE'] = getSetting(config, 'MESH_MOTION', 'DRIVE', 'prescribed') if bodies else ''
    if crank is not None:
        summary['LOAD'] = getSetting(config, 'MESH_MOTION', 'LOAD')
        summary['INERTIA'] = getSetting(config, 'MESH_MOTION', 'INERTIA')
        cyclesDone = int(np.floor((crank[-1, 1] - theta0)/360.0 + 1e-9))
    else:
        cyclesDone = int(np.floor(endTime*freq + 1e-6))
    summary['END_TIME'] = round(endTime, 6)
    summary['CYCLES_DONE'] = cyclesDone
    summary['FINISHED'] = 'yes' if summary['COMPLETION_DATE'] else 'no'
    nAvg = min(avgCycles, cyclesDone)
    summary['AVG_CYCLES'] = nAvg
    if nAvg == 0:
        print('\tWARNING: no complete cycle yet (end time %.6g s at %g Hz), results left empty' % (endTime, freq))
        summary['POST_TIME'] = hours(time.time() - start + (exportTime or 0.0))
        return summary
    if crank is not None:
        thetaEdges = theta0 + 360.0*np.array([cyclesDone - nAvg, cyclesDone])
        t0, t1 = (float(t) for t in np.interp(thetaEdges, crank[:, 1], crank[:, 0]))
    else:
        t0, t1 = (cyclesDone - nAvg)/freq, cyclesDone/freq
    print('\taveraging over cycles %d to %d (t = %.6g to %.6g s)' % (cyclesDone - nAvg + 1, cyclesDone, t0, t1))

    #solver quality over the same cycles: Courant number (PISO is only accurate around 1) and p solves cut off at the
    #iteration cap (a stalled pressure solve looks like any other step otherwise)
    pCap = int(float(getSetting(config, 'GLOBAL_SIM_CONTROL', 'P_MAX_ITER', '1000') or 1000))
    coMax, coAbove, pCapped, nSteps = solverQuality(solverLogs(caseDir, solver), t0, t1, pCap)
    if nSteps:
        summary['COURANT_MAX'] = round(coMax, 4) if coMax is not None else ''
        summary['COURANT_ABOVE_1'] = round(coAbove, 4) if coAbove is not None else ''
        summary['P_CAP_FRACTION'] = round(pCapped, 4) if pCapped is not None else ''
        nOuter = int(float(getSetting(config, 'GLOBAL_SIM_CONTROL', 'N_OUTER_CORR', '1') or 1))
        if coAbove is not None and coAbove > 0.05 and nOuter == 1:
            print('\tWARNING: Courant number above 1 in %.0f %% of the steps (max %.3g) with N_OUTER_CORR 1: lower '
                  'MAX_DT or use 3 outer correctors' % (100*coAbove, coMax))
        if pCapped is not None and pCapped > 0.02:
            print('\tWARNING: the final p solve stopped at P_MAX_ITER %d in %.0f %% of the steps: check the p '
                  'residuals in the log (round-off floor or a stalled solve)' % (pCap, 100*pCapped))

    #gas averages over the window (time-weighted, the time step varies)
    if gasAvg is not None:
        sel = window(gasAvg[:, 0], t0, t1)
        for k, name in enumerate(names[1:], start=1):
            field = re.sub(r'^volAverage\((.*)\)$', r'\1', name).upper()
            summary['%s_MEAN' % (field)] = timeMean(gasAvg[:, 0], gasAvg[:, k], t0, t1)
            if field == 'P':
                summary['P_MAX'] = float(gasAvg[sel, k].max())
                summary['P_MIN'] = float(gasAvg[sel, k].min())

    #crank drive: work per body = integral of its crank torque over the crank angle, speed and load
    if crank is not None:
        col = dict((name, k) for k, name in enumerate(crankNames))

        def overCycles(values):
            #values over the averaging cycles, ends interpolated to the exact crank angles
            return clipped(crank[:, 1], values, thetaEdges[0], thetaEdges[1])[1]

        theta = np.radians(clipped(crank[:, 1], crank[:, 1], thetaEdges[0], thetaEdges[1])[0])
        t = overCycles(crank[:, 0])
        totalWork = 0.0
        for body in bodies:
            if 'T_%s' % (body) not in col:
                print('\tWARNING: no crank torque for %s in the crank log, its work is left out' % (body))
                continue
            work = integrate(theta, overCycles(crank[:, col['T_%s' % (body)]]))/nAvg
            summary['WORK_%s' % (body)] = work
            fz = overCycles(crank[:, col['F_%s' % (body)]])
            summary['FORCE_Z_AMP_%s' % (body)] = float(0.5*(fz.max() - fz.min()))
            totalWork += work
        omega = overCycles(crank[:, col['omega']])
        freqMean = nAvg/(t1 - t0)
        summary['AVG_ENERGY'] = totalWork
        summary['AVG_POWER'] = totalWork*freqMean
        summary['FREQ_MEAN'] = freqMean
        summary['FREQ_MIN'] = float(omega.min()/(2.0*np.pi))
        summary['FREQ_MAX'] = float(omega.max()/(2.0*np.pi))
        summary['TORQUE_GAS_MEAN'] = integrate(theta, overCycles(crank[:, col['torqueGas']]))/(theta[-1] - theta[0])
        summary['POWER_LOAD'] = integrate(t, overCycles(crank[:, col['torqueLoad']]*crank[:, col['omega']]))/(t1 - t0)
        summary['POWER_ELEC'] = integrate(t, overCycles(crank[:, col['powerElec']]))/(t1 - t0)
        #convergence: mean speed of the last cycle minus the one before
        if cyclesDone >= 2:
            tc = np.interp(theta0 + 360.0*np.array([cyclesDone - 2, cyclesDone - 1, cyclesDone]), crank[:, 1], crank[:, 0])
            summary['FREQ_CHANGE'] = float(1.0/(tc[2] - tc[1]) - 1.0/(tc[1] - tc[0]))

    #prescribed motion: work on the moving bodies = integral of F_z v dt per cycle
    else:
        totalWork = 0.0
        for body in bodies:
            fNames, force = readDat(os.path.join(post, 'force_%s' % (body)), 'force')
            table = os.path.join(caseDir, 'constant', 'tables', '%sTable' % (body))
            if force is None or not os.path.isfile(table):
                print('\tWARNING: no force output or motion table for %s, its work is left out' % (body))
                continue
            tTable, dz = readTable(table)
            sel = window(force[:, 0], t0, t1)
            t = force[sel, 0]
            fz = force[sel, fNames.index('total_z')]
            #porous matrix moving with the body: the gas drag on it is part of the force on the body
            dragZ = porousDragZ(config, post, body, t)
            if dragZ is not None:
                summary['WORK_DRAG_%s' % (body)] = integrate(t, dragZ*tableVelocity(tTable, dz, t))/nAvg
                fz = fz + dragZ
            work = integrate(t, fz*tableVelocity(tTable, dz, t))/nAvg
            summary['WORK_%s' % (body)] = work
            summary['FORCE_Z_AMP_%s' % (body)] = float(0.5*(fz.max() - fz.min()))
            totalWork += work
        summary['AVG_ENERGY'] = totalWork
        summary['AVG_POWER'] = totalWork*freq

    #pressure on the moving faces (p_<patch> function objects)
    for foDir in sorted(glob.glob(os.path.join(post, 'p_*'))):
        pNames, pFace = readDat(foDir, 'surfaceFieldValue')
        if pFace is None:
            continue
        sel = window(pFace[:, 0], t0, t1)
        patch = os.path.basename(foDir)[2:]
        summary['P_%s_MEAN' % (patch)] = timeMean(pFace[:, 0], pFace[:, 1], t0, t1)
        summary['P_%s_AMP' % (patch)] = float(0.5*(pFace[sel, 1].max() - pFace[sel, 1].min()))

    #compressible: mass drift, p-V work, heat per BC group (positive = into the gas)
    if compressible:
        mNames, mass = readDat(os.path.join(post, 'gasMass'), 'volFieldValue')
        if mass is not None:
            summary['MASS_DRIFT'] = float(mass[-1, 1]/mass[0, 1] - 1.0)
            if gasAvg is not None and 'volAverage(rho)' in names:
                pAvg = np.interp(mass[:, 0], gasAvg[:, 0], gasAvg[:, names.index('volAverage(p)')])
                rhoAvg = np.interp(mass[:, 0], gasAvg[:, 0], gasAvg[:, names.index('volAverage(rho)')])
                volume = mass[:, 1]/rhoAvg
                sel = window(mass[:, 0], t0, t1)
                summary['PV_WORK'] = integrate(volume[sel], pAvg[sel])/nAvg
                summary['GAS_VOLUME_MIN'] = float(volume[sel].min())
                summary['GAS_VOLUME_MAX'] = float(volume[sel].max())
                gasT = None
                if 'volAverage(T)' in names:
                    gasT = np.interp(mass[:, 0], gasAvg[:, 0], gasAvg[:, names.index('volAverage(T)')])
                summary.update(idealGasSummary(caseDir, config, mass[:, 0], mass[:, 1], volume, pAvg, gasT,
                                               t0, t1))
        summary.update(heatSummary(config, post, t0, t1))

    #inlets and outlets of a through-flow case, porous zones
    summary.update(flowSummary(caseDir, config, t0, t1, compressible))
    summary.update(porousSummary(post, t0, t1))
    summary['POST_TIME'] = hours(time.time() - start + (exportTime or 0.0))
    return summary


#----- plots

def cycleEdges(freq, cyclesDone, crank, theta0):
    #start and end time of every complete cycle (by crank angle for the crank drive)
    if crank is not None:
        return np.interp(theta0 + 360.0*np.arange(cyclesDone + 1), crank[:, 1], crank[:, 0])
    return np.arange(cyclesDone + 1)/freq


def timeMean(t, y, t0, t1):
    #time average of y(t) between t0 and t1
    tc, yc = clipped(t, y, t0, t1)
    return integrate(tc, yc)/(t1 - t0)


def cycleSeries(caseDir):
    #values per complete cycle for the convergence plot: work per body (J), mean gas T, mass drift,
    #heat per BC group (J, + into the gas) and mean speed (crank drive)
    caseDir = os.path.abspath(caseDir)
    config = readCaseSetup(caseDir)
    post = os.path.join(caseDir, 'postProcessing')
    freq = float(getSetting(config, 'GLOBAL_CONTROL', 'FREQ'))
    bodies = getBodies(config)
    names, gasAvg = readDat(os.path.join(post, 'gasAvg'), 'volFieldValue')
    crankNames, crank, theta0 = readCrank(caseDir, post)
    if crank is not None:
        cyclesDone = int(np.floor((crank[-1, 1] - theta0)/360.0 + 1e-9))
    elif gasAvg is not None:
        cyclesDone = int(np.floor(gasAvg[-1, 0]*freq + 1e-6))
    else:
        return None
    if cyclesDone == 0:
        return None
    edges = cycleEdges(freq, cyclesDone, crank, theta0)
    series = {'cycle': np.arange(1, cyclesDone + 1), 'work': {}, 'heat': {}}

    #work per body and cycle, as in the summary
    if crank is not None:
        col = dict((name, k) for k, name in enumerate(crankNames))
        for body in bodies:
            if 'T_%s' % (body) not in col:
                continue
            work = []
            for k in range(cyclesDone):
                theta, torque = clipped(crank[:, 1], crank[:, col['T_%s' % (body)]],
                                        theta0 + 360.0*k, theta0 + 360.0*(k + 1))
                work.append(integrate(np.radians(theta), torque))
            series['work'][body] = np.array(work)
        series['freq'] = 1.0/np.diff(edges)
    else:
        for body in bodies:
            fNames, force = readDat(os.path.join(post, 'force_%s' % (body)), 'force')
            table = os.path.join(caseDir, 'constant', 'tables', '%sTable' % (body))
            if force is None or not os.path.isfile(table):
                continue
            tTable, dz = readTable(table)
            fz = force[:, fNames.index('total_z')]
            dragZ = porousDragZ(config, post, body, force[:, 0])
            if dragZ is not None:
                fz = fz + dragZ
            power = fz*tableVelocity(tTable, dz, force[:, 0])
            work = []
            for k in range(cyclesDone):
                t, p = clipped(force[:, 0], power, edges[k], edges[k + 1])
                work.append(integrate(t, p))
            series['work'][body] = np.array(work)

    #gas state per cycle
    if gasAvg is not None and 'volAverage(T)' in names:
        series['T'] = np.array([timeMean(gasAvg[:, 0], gasAvg[:, names.index('volAverage(T)')], edges[k], edges[k + 1])
                                for k in range(cyclesDone)])
    mNames, mass = readDat(os.path.join(post, 'gasMass'), 'volFieldValue')
    if mass is not None:
        series['massDrift'] = np.interp(edges[1:], mass[:, 0], mass[:, 1])/mass[0, 1] - 1.0

    #heat per BC group and cycle
    heat = readWallHeatFlux(os.path.join(post, 'wallHeatFlux'))
    groupPatches = dict(config.items('BC_GROUPS_RESOLVED')) if config.has_section('BC_GROUPS_RESOLVED') else {}
    for group, patches in groupPatches.items():
        rows = [heat[patch] for patch in patches.split() if patch in heat]
        if not rows:
            continue
        series['heat'][group] = np.array([sum(integrate(*clipped(r[:, 0], r[:, 1], edges[k], edges[k + 1])) for r in rows)
                                          for k in range(cyclesDone)])
    return series


def pvLoops(caseDir, nLoops=3):
    #p and gas volume V = mass/rho of the last nLoops complete cycles (compressible runs)
    caseDir = os.path.abspath(caseDir)
    config = readCaseSetup(caseDir)
    post = os.path.join(caseDir, 'postProcessing')
    freq = float(getSetting(config, 'GLOBAL_CONTROL', 'FREQ'))
    names, gasAvg = readDat(os.path.join(post, 'gasAvg'), 'volFieldValue')
    mNames, mass = readDat(os.path.join(post, 'gasMass'), 'volFieldValue')
    if gasAvg is None or mass is None or 'volAverage(rho)' not in names:
        return []
    crankNames, crank, theta0 = readCrank(caseDir, post)
    if crank is not None:
        cyclesDone = int(np.floor((crank[-1, 1] - theta0)/360.0 + 1e-9))
    else:
        cyclesDone = int(np.floor(gasAvg[-1, 0]*freq + 1e-6))
    edges = cycleEdges(freq, cyclesDone, crank, theta0)
    p = np.interp(mass[:, 0], gasAvg[:, 0], gasAvg[:, names.index('volAverage(p)')])
    volume = mass[:, 1]/np.interp(mass[:, 0], gasAvg[:, 0], gasAvg[:, names.index('volAverage(rho)')])
    loops = []
    for k in range(max(cyclesDone - nLoops, 0), cyclesDone):
        sel = window(mass[:, 0], edges[k], edges[k + 1])
        loops.append((k + 1, volume[sel], p[sel]))
    return loops


def bodyDisplacement(caseDir, config, body, t, theta):
    #z displacement of a body at times t (crank angle theta in deg, accumulated), None without a motion table
    tables = os.path.join(caseDir, 'constant', 'tables')
    if os.path.isfile(os.path.join(tables, '%sCrank' % (body))):
        angle, dz = readCrankTable(os.path.join(tables, '%sCrank' % (body)))
        return np.interp(np.mod(theta, 360.0), angle, dz, period=360.0)
    if os.path.isfile(os.path.join(tables, '%sTable' % (body))):
        times, dz = readTable(os.path.join(tables, '%sTable' % (body)))
        return np.interp(t, times, dz)
    return None


def crownPressure(caseDir, config, post, t, theta, volume):
    #area-averaged p on the faces of the body that changes the gas volume (the piston crown): its p-V loop is the
    #piston work. The body is the one whose table displacement follows the gas volume, its faces the p_<source>_<pid>
    #function objects (ZONE_PIDS keywords). Returns (label, p) or (None, None)
    best, bestCorr = None, 0.9
    for body in getBodies(config):
        dz = bodyDisplacement(caseDir, config, body, t, theta)
        if dz is None or np.std(dz) == 0 or np.std(volume) == 0:
            continue
        corr = abs(float(np.corrcoef(dz, volume)[0, 1]))
        if corr > bestCorr:
            best, bestCorr = body, corr
    if best is None:
        return None, None
    section = 'MOVE-%s' % (best)
    source = getSetting(config, section, 'FROM') or section
    keywords = [k.lower() for k in getSetting(config, section, 'ZONE_PIDS').split()]
    pSum, areaSum, patches = np.zeros_like(t), 0.0, []
    for foDir in sorted(glob.glob(os.path.join(post, 'p_%s_*' % (source)))):
        pid = os.path.basename(foDir)[len('p_%s_' % (source)):]
        if keywords and not any(k in pid.lower() for k in keywords):
            continue
        names, rows = readDat(foDir, 'surfaceFieldValue')
        if rows is None:
            continue
        area = patchArea(foDir) or 1.0
        pSum += area*np.interp(t, rows[:, 0], rows[:, 1])
        areaSum += area
        patches.append(pid)
    if not patches:
        return None, None
    return 'pressure on the %s faces (%s)' % (best, ' '.join(patches)), pSum/areaSum


def pvData(caseDir, nLoops=3):
    #p-V data of the last nLoops cycles (compressible runs): per cycle the crank angle in the cycle, gas volume
    #V = mass/rho, volume-averaged p and the piston crown p; an unfinished last cycle is kept (closed False)
    caseDir = os.path.abspath(caseDir)
    config = readCaseSetup(caseDir)
    post = os.path.join(caseDir, 'postProcessing')
    names, gasAvg = readDat(os.path.join(post, 'gasAvg'), 'volFieldValue')
    mNames, mass = readDat(os.path.join(post, 'gasMass'), 'volFieldValue')
    if gasAvg is None or mass is None or 'volAverage(rho)' not in names:
        return None
    t = mass[:, 0]
    crankNames, crank, theta0 = readCrank(caseDir, post)
    if crank is not None:
        theta = np.interp(t, crank[:, 0], crank[:, 1])
    else:
        #prescribed motion: the crank turns at FREQ from the start angle
        freq = float(getSetting(config, 'GLOBAL_CONTROL', 'FREQ'))
        theta0 = float(getSetting(config, 'SETUP_INFO', 'START_ANGLE', '0') or 0)
        theta = theta0 + 360.0*freq*t
    volume = mass[:, 1]/np.interp(t, gasAvg[:, 0], gasAvg[:, names.index('volAverage(rho)')])
    pAvg = np.interp(t, gasAvg[:, 0], gasAvg[:, names.index('volAverage(p)')])
    crownLabel, pCrown = crownPressure(caseDir, config, post, t, theta, volume)
    phase = theta - theta0
    nCycles = int(np.ceil(phase[-1]/360.0 - 1e-6))
    cycles = []
    for k in range(max(nCycles - nLoops, 0), nCycles):
        sel = (phase >= 360.0*k - 1e-9) & (phase <= 360.0*(k + 1) + 1e-9)
        if sel.sum() < 10:
            continue
        #closed: a later cycle follows, or the data reach the cycle end within two samples
        closed = k < nCycles - 1 or phase[-1] >= 360.0*(k + 1) - 2.0*float(np.max(np.diff(phase[sel])))
        cycles.append({'cycle': k + 1, 'angle': phase[sel] - 360.0*k, 'volume': volume[sel], 'pAvg': pAvg[sel],
                       'pCrown': None if pCrown is None else pCrown[sel], 'mass': mass[sel, 1],
                       'closed': bool(closed)})
    return {'theta0': theta0, 'crown': crownLabel, 'cycles': cycles}


def flowSeries(caseDir):
    #steady convergence: initial residuals and, per inlet/outlet group, flow (+ into the domain), p and T vs iteration
    caseDir = os.path.abspath(caseDir)
    config = readCaseSetup(caseDir)
    post = os.path.join(caseDir, 'postProcessing')
    series = {'residuals': {}, 'flow': {}, 'p': {}, 'T': {}}
    compressible = getSetting(config, 'SETUP_INFO', 'FLOW', 'compressible') == 'compressible'
    flowKey, fluxScale, pScale, pOffset = flowUnits(config, compressible)
    series['flowKey'] = flowKey
    rNames, residuals = readDat(os.path.join(post, 'residuals'), 'solverInfo', textColumns=True)
    if residuals is not None:
        series['residualIter'] = residuals[:, 0]
        for k, name in enumerate(rNames):
            if name.endswith('_initial'):
                series['residuals'][name[:-len('_initial')]] = residuals[:, k]
    for group, (kind, patches) in flowGroups(config).items():
        iters, flux, pSum, tSum, weight = None, None, None, None, 0.0
        for patch in patches:
            names, rows = readDat(os.path.join(post, 'flow_%s' % (patch)), 'surfaceFieldValue')
            aNames, avg = readDat(os.path.join(post, 'avg_%s' % (patch)), 'surfaceFieldValue')
            if rows is None or avg is None:
                continue
            if iters is None:
                iters = rows[:, 0]
                flux, pSum, tSum = np.zeros_like(iters), np.zeros_like(iters), np.zeros_like(iters)
            area = patchArea(os.path.join(post, 'avg_%s' % (patch))) or 1.0
            flux -= fluxScale*np.interp(iters, rows[:, 0], rows[:, 1])
            pSum += area*(pOffset + pScale*np.interp(iters, avg[:, 0], avg[:, aNames.index('areaAverage(p)')]))
            if 'areaAverage(T)' in aNames:
                tSum += area*np.interp(iters, avg[:, 0], avg[:, aNames.index('areaAverage(T)')])
            weight += area
        if iters is None:
            continue
        series['flow'][group] = (iters, flux)
        series['p'][group] = (iters, pSum/weight)
        if tSum.any():
            series['T'][group] = (iters, tSum/weight)
    return series


def writeFlowPlot(caseDir, plt):
    #flowConvergence.png for steady runs and transient runs without cycles (against time)
    series = flowSeries(caseDir)
    panels = []
    if series['residuals']:
        panels.append(('residuals', 'initial residual'))
    if series['flow']:
        massFlow = series['flowKey'] == 'MDOT'
        panels.append(('flow', 'mass flow [kg/s], + in' if massFlow else 'volume flow [m$^3$/s], + in'))
        panels.append(('p', 'area-averaged p [Pa]' if massFlow else 'area-averaged p [m$^2$/s$^2$]'))
    if series['T']:
        panels.append(('T', 'area-averaged T [K]'))
    if not panels:
        return None
    fig, axes = plt.subplots(len(panels), 1, figsize=(8, 2.8*len(panels)), sharex=True, squeeze=False)
    for ax, (key, label) in zip(axes[:, 0], panels):
        if key == 'residuals':
            for name, values in series['residuals'].items():
                ax.semilogy(series['residualIter'], np.maximum(values, 1e-12), label=name)
        else:
            for group, (iters, values) in series[key].items():
                ax.plot(iters, values, label=group)
            ax.ticklabel_format(axis='y', useOffset=False)
        ax.set_ylabel(label, fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(True)
    steady = getSetting(readCaseSetup(os.path.abspath(caseDir)), 'SETUP_INFO', 'STEADY', 'no') == 'yes'
    axes[-1, 0].set_xlabel('iteration' if steady else 'time [s]')
    axes[0, 0].set_title('%s convergence' % (os.path.basename(os.path.abspath(caseDir))))
    path = os.path.join(os.path.abspath(caseDir), 'flowConvergence.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return path


#plot colours: categorical order (cycles), text and grid inks
PLOT_COLOURS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
PLOT_INK, PLOT_INK2, PLOT_GRID = '#0b0b0b', '#52514e', '#e4e3df'


def writePvPlot(caseDir, plt, nLoops=3):
    #pV.png: p-V loops of the last cycles (area = piston work, + out of the gas), p and gas volume vs crank angle
    #and the loop shape (deviation from a polytropic fit of each cycle); crown p when the piston faces are known
    data = pvData(caseDir, nLoops)
    if not data or not data['cycles']:
        return None
    caseName = os.path.basename(os.path.abspath(caseDir))
    key, pLabel = ('pCrown', data['crown']) if data['crown'] else ('pAvg', 'volume-averaged p')
    style = {'axes.edgecolor': PLOT_INK2, 'axes.labelcolor': PLOT_INK, 'xtick.color': PLOT_INK2,
             'ytick.color': PLOT_INK2, 'grid.color': PLOT_GRID, 'axes.spines.top': False,
             'axes.spines.right': False, 'legend.frameon': False}
    with plt.rc_context(style):
        fig, ax = plt.subplots(2, 2, figsize=(12, 9))
        for i, c in enumerate(data['cycles']):
            colour = PLOT_COLOURS[i % len(PLOT_COLOURS)]
            v, p = c['volume']*1e6, c[key]/1e5
            if c['closed']:
                label = 'cycle %d: %+.2f J' % (c['cycle'], integrate(c['volume'], c[key]))
            else:
                label = 'cycle %d (%.0f deg, not closed)' % (c['cycle'], c['angle'][-1] - c['angle'][0])
            ax[0, 0].plot(v, p, color=colour, lw=1.6, label=label)
            ax[0, 1].plot(c['angle'], p, color=colour, lw=1.6, label='cycle %d' % (c['cycle']))
            #loop shape: p over the polytropic fit p = C V^-n of this cycle
            fit = np.linalg.lstsq(np.column_stack([np.ones(len(v)), np.log(v)]), np.log(p), rcond=None)[0]
            ax[1, 1].plot(v, 100.0*(p/np.exp(fit[0] + fit[1]*np.log(v)) - 1.0), color=colour, lw=1.6,
                          label='cycle %d (n = %.3f)' % (c['cycle'], -fit[1]))
        last = data['cycles'][-1]
        #ideal-gas limits through the last cycle: isothermal p = m R T_ref / V, adiabatic p V^gamma through its largest volume
        R, cp = gasConstants(os.path.abspath(caseDir))
        if R is not None:
            tRef, source = isothermalT(readCaseSetup(os.path.abspath(caseDir)))
            gamma = cp/(cp - R)
            v = np.linspace(last['volume'].min(), last['volume'].max(), 200)
            k = int(np.argmax(last['volume']))
            ax[0, 0].plot(v*1e6, float(np.mean(last['mass']))*R*tRef/v/1e5, '--', color=PLOT_INK2, lw=1.0,
                          label='isothermal (%g K)' % (tRef))
            ax[0, 0].plot(v*1e6, last[key][k]*(last['volume'][k]/v)**gamma/1e5, ':', color=PLOT_INK2, lw=1.2,
                          label='adiabatic (gamma = %.3f)' % (gamma))
        ax[1, 0].plot(last['angle'], last['volume']*1e6, color=PLOT_INK2, lw=1.6)
        ax[0, 0].set(xlabel='gas volume [cm$^3$]', ylabel='p [bar]', title='p-V loops (area = work per cycle, + out of the gas)')
        ax[0, 1].set(xlabel='crank angle from the start angle [deg]', ylabel='p [bar]', title='p vs crank angle')
        ax[1, 0].set(xlabel='crank angle from the start angle [deg]', ylabel='gas volume [cm$^3$]',
                     title='gas volume vs crank angle (cycle %d)' % (last['cycle']))
        ax[1, 1].axhline(0.0, color=PLOT_INK2, lw=0.8)
        ax[1, 1].set(xlabel='gas volume [cm$^3$]', ylabel='p / polytropic fit - 1 [%]',
                     title='loop shape: deviation from p V$^n$ = const')
        for a in ax.flat:
            a.grid(True, lw=0.6)
        for a in (ax[0, 0], ax[0, 1], ax[1, 1]):
            a.legend(fontsize=9)
        for a in (ax[0, 1], ax[1, 0]):
            a.set_xlim(0, 360)
            a.set_xticks(range(0, 361, 45))
        fig.suptitle('%s: %s, start angle %g deg' % (caseName, pLabel, data['theta0']), color=PLOT_INK)
        fig.tight_layout()
        path = os.path.join(os.path.abspath(caseDir), 'pV.png')
        fig.savefig(path, dpi=150)
        plt.close(fig)
    return path


def writePlots(caseDir, nLoops=3):
    #pV.png (last cycles overlaid) and cycleConvergence.png in the case directory
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('\tWARNING: no matplotlib, no plots')
        return []
    caseDir = os.path.abspath(caseDir)
    caseName = os.path.basename(caseDir)
    written = []
    #steady runs and transient runs without cycles: flowConvergence.png only
    config = readCaseSetup(caseDir)
    if (getSetting(config, 'SETUP_INFO', 'STEADY', 'no') == 'yes'
            or getSetting(config, 'SETUP_INFO', 'CYCLIC', 'yes') == 'no'):
        path = writeFlowPlot(caseDir, plt)
        return [path] if path else []

    path = writePvPlot(caseDir, plt, nLoops)
    if path:
        written.append(path)

    series = cycleSeries(caseDir)
    if series is not None:
        panels = [('work', 'work per cycle [J]')]
        if 'T' in series:
            panels.append(('T', 'mean gas T [K]'))
        if 'massDrift' in series:
            panels.append(('massDrift', 'gas mass drift [-]'))
        if series['heat']:
            panels.append(('heat', 'heat per cycle [J], + into gas'))
        if 'freq' in series:
            panels.append(('freq', 'mean speed [Hz]'))
        fig, axes = plt.subplots(len(panels), 1, figsize=(8, 2.8*len(panels)), sharex=True, squeeze=False)
        cycle = series['cycle']
        for ax, (key, label) in zip(axes[:, 0], panels):
            if key in ('work', 'heat'):
                for name, values in series[key].items():
                    ax.plot(cycle, values, marker='o', label=name)
                if key == 'work' and len(series['work']) > 1:
                    ax.plot(cycle, sum(series['work'].values()), marker='o', color='k', label='total')
                if series[key]:
                    ax.legend(fontsize=8)
            else:
                ax.plot(cycle, series[key], marker='o')
            ax.set_ylabel(label, fontsize=9)
            ax.ticklabel_format(axis='y', useOffset=False)
            ax.grid(True)
        axes[-1, 0].set_xlabel('cycle')
        axes[-1, 0].xaxis.get_major_locator().set_params(integer=True)
        axes[0, 0].set_title('%s cycle convergence' % (caseName))
        path = os.path.join(caseDir, 'cycleConvergence.png')
        fig.savefig(path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        written.append(path)
    return written


def writeSummary(summary, caseDir):
    outDir = os.path.join(caseDir)
    os.makedirs(outDir, exist_ok=True)
    path = os.path.join(outDir, 'summary.csv')
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        for key, value in summary.items():
            writer.writerow([key, value])
    return path


def readSummary(caseDir):
    #summary.csv back as {key: value}, numbers as float/int (the identity columns stay text: case 003 is not 3)
    path = os.path.join(caseDir, 'summary.csv')
    if not os.path.isfile(path):
        sys.exit('ERROR! %s not found: run postInternal.py in the case first' % (path))
    summary = {}
    with open(path, newline='') as f:
        for row in csv.reader(f):
            if not row:
                continue
            value = row[1] if len(row) > 1 else ''
            if row[0] in ('JOB', 'NAME', 'DESCRIPTION', 'COMMENTS'):
                summary[row[0]] = value
                continue
            try:
                value = int(value) if re.match(r'^-?\d+$', value) else float(value)
            except ValueError:
                pass
            summary[row[0]] = value
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description='Summary and plots of an internalSetup case (summary.csv, pV.png, cycleConvergence.png)')
    parser.add_argument('cases', nargs='*', default=['.'], help='case directories (default: current directory)')
    parser.add_argument('--avgCycles', type=int, default=1, help='average results over the last N complete cycles (default 1)')
    parser.add_argument('--avgTime', type=float, default=None,
                        help='transient runs without cycles (END_TIME): average over the last AVGTIME s (default: last 10 %%)')
    parser.add_argument('--pvCycles', type=int, default=3, help='number of last cycles overlaid in pV.png (default 3)')
    parser.add_argument('--noPlots', action='store_true', help='summary only')
    parser.add_argument('--gsheet', action='store_true', help='push the summary to the database gsheet (gsheetInternal.py)')
    args = parser.parse_args(argv)
    for case in args.cases:
        print('\n\tPost-processing %s...' % (os.path.abspath(case)))
        summary = buildSummary(case, args.avgCycles, avgTime=args.avgTime)
        path = writeSummary(summary, case)
        print('\twritten %s' % (path))
        for key in ['CYCLES_DONE', 'AVG_ENERGY', 'AVG_POWER', 'MASS_DRIFT']:
            if summary.get(key, '') != '':
                print('\t\t%-12s %s' % (key, summary[key]))
        if not args.noPlots:
            for path in writePlots(case, args.pvCycles):
                print('\twritten %s' % (path))
    if args.gsheet:
        sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
        import gsheetInternal
        gsheetInternal.main(args.cases)


if __name__ == '__main__':
    main()
