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

#post-processing of an internal (internalSetup) case: one summary per run in postInternal/summary.csv
#(key,value rows), the database row pushed to the Ekstera gsheet by gsheetEkstera.py
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
SUMMARY_DIR = 'postInternal'


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

def readDat(foDir, fileName):
    #all <foDir>/<startTime>/<fileName>*.dat (restarts start a new directory or file), merged in time order;
    #later files replace overlapping times. returns (column names, array)
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
                values = line.split()
                if values:
                    try:
                        rows.append([float(v) for v in values])
                    except ValueError:
                        continue
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


def tableVelocity(times, dz, t):
    #the solver interpolates the table linearly in time, so the velocity is constant between table points
    slopes = np.diff(dz)/np.diff(times)
    k = np.clip(np.searchsorted(times, t, side='right') - 1, 0, len(slopes) - 1)
    return slopes[k]


def window(times, t0, t1):
    return (times >= t0 - 1e-12) & (times <= t1 + 1e-12)


def integrate(t, y):
    #trapezoidal rule (np.trapz is deprecated in newer numpy)
    t, y = np.asarray(t, float), np.asarray(y, float)
    return float(np.sum(0.5*(y[1:] + y[:-1])*np.diff(t))) if len(t) > 1 else 0.0


#----- summary

def buildSummary(caseDir, avgCycles=1):
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
    freq = float(getSetting(config, 'GLOBAL_CONTROL', 'FREQ'))
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

    #run status: complete cycles from t = 0 (the motion tables start at 0)
    post = os.path.join(caseDir, 'postProcessing')
    names, gasAvg = readDat(os.path.join(post, 'gasAvg'), 'volFieldValue')
    endTime = float(gasAvg[-1, 0]) if gasAvg is not None else 0.0
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
    t0, t1 = (cyclesDone - nAvg)/freq, cyclesDone/freq
    print('\taveraging over cycles %d to %d (t = %.6g to %.6g s)' % (cyclesDone - nAvg + 1, cyclesDone, t0, t1))

    #gas averages over the window
    if gasAvg is not None:
        sel = window(gasAvg[:, 0], t0, t1)
        for k, name in enumerate(names[1:], start=1):
            field = re.sub(r'^volAverage\((.*)\)$', r'\1', name).upper()
            summary['%s_MEAN' % (field)] = float(gasAvg[sel, k].mean())
            if field == 'P':
                summary['P_MAX'] = float(gasAvg[sel, k].max())
                summary['P_MIN'] = float(gasAvg[sel, k].min())

    #work on the moving bodies: integral of F_z v dt per cycle
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
        summary['P_%s_MEAN' % (patch)] = float(pFace[sel, 1].mean())
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
        heat = readWallHeatFlux(os.path.join(post, 'wallHeatFlux'))
        groupPatches = dict(config.items('BC_GROUPS_RESOLVED')) if config.has_section('BC_GROUPS_RESOLVED') else {}
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
                summary['HEAT_%s' % (group)] = total

    summary['POST_TIME'] = hours(time.time() - start + (exportTime or 0.0))
    return summary


def writeSummary(summary, caseDir):
    outDir = os.path.join(caseDir, SUMMARY_DIR)
    os.makedirs(outDir, exist_ok=True)
    path = os.path.join(outDir, 'summary.csv')
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        for key, value in summary.items():
            writer.writerow([key, value])
    return path


def readSummary(caseDir):
    #summary.csv back as {key: value}, numbers as float/int
    path = os.path.join(caseDir, SUMMARY_DIR, 'summary.csv')
    if not os.path.isfile(path):
        sys.exit('ERROR! %s not found: run postInternal.py in the case first' % (path))
    summary = {}
    with open(path, newline='') as f:
        for row in csv.reader(f):
            if not row:
                continue
            value = row[1] if len(row) > 1 else ''
            try:
                value = int(value) if re.match(r'^-?\d+$', value) else float(value)
            except ValueError:
                pass
            summary[row[0]] = value
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description='Summary of an internalSetup case (postInternal/summary.csv)')
    parser.add_argument('cases', nargs='*', default=['.'], help='case directories (default: current directory)')
    parser.add_argument('--avgCycles', type=int, default=1, help='average results over the last N complete cycles (default 1)')
    parser.add_argument('--gsheet', action='store_true', help='push the summary to the Ekstera gsheet (gsheetEkstera.py)')
    args = parser.parse_args(argv)
    for case in args.cases:
        print('\n\tPost-processing %s...' % (os.path.abspath(case)))
        summary = buildSummary(case, args.avgCycles)
        path = writeSummary(summary, case)
        print('\twritten %s' % (path))
        for key in ['CYCLES_DONE', 'AVG_ENERGY', 'AVG_POWER', 'MASS_DRIFT']:
            if summary.get(key, '') != '':
                print('\t\t%-12s %s' % (key, summary[key]))
    if args.gsheet:
        sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
        import gsheetEkstera
        gsheetEkstera.main(args.cases)


if __name__ == '__main__':
    main()
