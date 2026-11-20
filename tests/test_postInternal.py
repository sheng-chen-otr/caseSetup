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

#tests for the internal post (summary) and the database gsheet push logic (no OpenFOAM, no Google API)
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'postUtilities'))
import postInternal
import gsheetInternal

FREQ = 10.0
AMP = 0.01          #piston amplitude (m)
DAMP = 50.0         #F_z = -DAMP*v (N s/m): the gas takes work from the piston
P_AMP, V_AMP = 2000.0, 1e-4


SETUP = """[TITLES]
CASENAME = run01
CASEDESCRP = synthetic test

[GLOBAL_CONTROL]
FREQ = 10
N_CYCLES = 2

[GLOBAL_SIM_CONTROL]
SOLVER = rhoPimpleFoam
TURB_MODEL = laminar

[GLOBAL_MATERIAL]
GAS = air
CHARGE_PRESSURE = 1e5
INIT_T = 300

[GLOBAL_REFINEMENT]
BASE_CELL = 0.004

[GLOBAL_COMPUTE_SETUP]
NUM_PROCS = 4

[WALL-chamber]
T_DEFAULT = adiabatic
BC_GROUPS = hotEnd
hotEnd = PIDS hot T fixed 600

[MOVE-piston]
FROM = WALL-chamber
TDC_CLEARANCE = 0.002

[SETUP_INFO]
TOOL = internal
VERSION = 5.0.0-dev
FLOW = compressible

[BC_GROUPS_RESOLVED]
hotEnd = WALL-chamber_hot
WALL-chamber_default = WALL-chamber_side
"""


def _log(path, date, clock, clockTimes, end=True):
    lines = ['Date   : %s' % (date), 'Time   : %s' % (clock), '']
    lines += ['ExecutionTime = %g s  ClockTime = %g s' % (c, c) for c in clockTimes]
    lines += ['', 'End', ''] if end else []
    path.write_text('\n'.join(lines))


def _dat(path, header, rows, comment='# Region : all'):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = [comment, '# %s' % ('\t'.join(header))]
    text += ['\t'.join('%.12g' % v for v in row) for row in rows]
    path.write_text('\n'.join(text) + '\n')


def makeCase(root, cycles=2.0):
    case = root / 'JOB01' / 'CASES' / 'run01'
    case.mkdir(parents=True)
    (case / 'fullCaseSetupDict').write_text(SETUP)
    (case / 'log.checkMesh').write_text('    cells:            123456\n')
    _log(case / 'log.snappyHexMesh', 'Nov 18 2026', '10:00:00', [1800])
    _log(case / 'log.blockMesh', 'Nov 18 2026', '09:59:00', [60])
    #two solver segments (restart), 1 h each
    _log(case / 'log.rhoPimpleFoam.1', 'Nov 18 2026', '11:00:00', [100, 3600], end=False)
    _log(case / 'log.rhoPimpleFoam', 'Nov 18 2026', '13:00:00', [3600])

    w = 2*math.pi*FREQ
    t = np.linspace(0.0, cycles/FREQ, 4001)[1:]
    #motion table: z = AMP sin(wt), fine enough that the piecewise-linear velocity is close to the analytic one
    tt = np.linspace(0.0, cycles/FREQ, 20001)
    tables = case / 'constant' / 'tables'
    tables.mkdir(parents=True)
    body = '\n'.join('(%.12g (0 0 %.12g))' % (a, AMP*math.sin(w*a)) for a in tt)
    (tables / 'pistonTable').write_text('FoamFile { class scalarVectorTable; }\n%d\n(\n%s\n)\n' % (len(tt), body))

    v = AMP*w*np.cos(w*t)
    fz = -DAMP*v
    forceRows = np.column_stack([t, 0*t, 0*t, fz, 0*t, 0*t, fz, 0*t, 0*t, 0*t])
    header = ['Time', 'total_x', 'total_y', 'total_z', 'pressure_x', 'pressure_y', 'pressure_z',
              'viscous_x', 'viscous_y', 'viscous_z']
    #restart at the half-way time: second directory starts there and overlaps one sample
    half = len(t)//2
    post = case / 'postProcessing'
    _dat(post / 'force_piston' / '0' / 'force.dat', header, forceRows[:half + 5], '# Force')
    _dat(post / 'force_piston' / ('%.12g' % t[half]) / 'force.dat', header, forceRows[half:], '# Force')

    #gas: V = V0 + V_AMP sin, p = p0 + P_AMP cos -> closed p-V loop of area pi*P_AMP*V_AMP
    volume = 1e-3 + V_AMP*np.sin(w*t)
    p = 1e5 + P_AMP*np.cos(w*t)
    mass = np.full_like(t, 1.2e-3)
    rho = mass/volume
    _dat(post / 'gasAvg' / '0' / 'volFieldValue.dat', ['Time', 'volAverage(p)', 'volAverage(T)', 'volAverage(rho)'],
         np.column_stack([t, p, 0*t + 300, rho]))
    _dat(post / 'gasMass' / '0' / 'volFieldValue.dat', ['Time', 'volIntegrate(rho)'], np.column_stack([t, mass]))
    _dat(post / 'p_WALL-chamber_crown' / '0' / 'surfaceFieldValue.dat', ['Time', 'areaAverage(p)'],
         np.column_stack([t, p]))
    heat = post / 'wallHeatFlux' / '0' / 'wallHeatFlux.dat'
    heat.parent.mkdir(parents=True)
    lines = ['# Wall heat-flux', '# Time  patch  min  max  integral']
    for a in t:
        lines.append('%.12g\tWALL-chamber_hot\t1\t2\t150' % (a))
        lines.append('%.12g\tWALL-chamber_side\t1\t2\t-20' % (a))
    heat.write_text('\n'.join(lines) + '\n')
    return case


#----- summary

def test_summary_values(tmp_path):
    case = makeCase(tmp_path)
    s = postInternal.buildSummary(str(case), avgCycles=1)
    assert s['JOB'] == 'JOB01' and s['NAME'] == 'run01' and s['DESCRIPTION'] == 'synthetic test'
    assert s['CELL_COUNT'] == 123456
    assert s['MESH_TIME'] == pytest.approx(1860/3600.0, abs=1e-4)
    #both solver segments count, the start date is the earliest one
    assert s['SOLVE_TIME'] == pytest.approx(2.0, abs=1e-4)
    assert s['RUN_DATE'] == '2026-11-18 11:00'
    assert s['COMPLETION_DATE'] != '' and s['FINISHED'] == 'yes'
    assert s['CYCLES_DONE'] == 2 and s['AVG_CYCLES'] == 1
    #damper: work by the gas on the piston per cycle = -DAMP*AMP^2*w^2*T/2
    w = 2*math.pi*FREQ
    work = -DAMP*AMP**2*w**2/(2*FREQ)
    assert s['WORK_piston'] == pytest.approx(work, rel=2e-3)
    assert s['AVG_ENERGY'] == pytest.approx(work, rel=2e-3)
    assert s['AVG_POWER'] == pytest.approx(work*FREQ, rel=2e-3)
    assert s['PV_WORK'] == pytest.approx(math.pi*P_AMP*V_AMP, rel=2e-3)
    assert s['MASS_DRIFT'] == pytest.approx(0.0, abs=1e-12)
    assert s['P_MEAN'] == pytest.approx(1e5, rel=1e-4)
    assert s['P_WALL-chamber_crown_AMP'] == pytest.approx(P_AMP, rel=1e-3)
    assert s['HEAT_hotEnd'] == pytest.approx(150.0)
    assert s['HEAT_WALL-chamber_default'] == pytest.approx(-20.0)
    assert s['BC_hotEnd'] == 'fixed 600' and s['TDC_CLEARANCE_piston'] == '0.002'


def test_summary_no_complete_cycle(tmp_path):
    case = makeCase(tmp_path, cycles=0.5)
    s = postInternal.buildSummary(str(case))
    assert s['CYCLES_DONE'] == 0 and s['AVG_ENERGY'] == '' and s['AVG_POWER'] == ''


def test_summary_csv_round_trip(tmp_path):
    case = makeCase(tmp_path)
    s = postInternal.buildSummary(str(case))
    postInternal.writeSummary(s, str(case))
    back = postInternal.readSummary(str(case))
    assert back['CELL_COUNT'] == 123456 and isinstance(back['AVG_POWER'], float)
    assert back['NAME'] == 'run01' and list(back)[:3] == ['JOB', 'NAME', 'DESCRIPTION']


#crank drive: speeding up from 10 Hz, piston torque B + A sin(theta), generator c*omega
CRANK = {'theta0': 30.0, 'w0': 2*math.pi*FREQ, 'alpha': 2*math.pi*5.0, 'A': 3.0, 'B': 0.5, 'c': 0.01}


def makeCrankCase(tmp_path):
    case = tmp_path / 'JOB01' / 'CASES' / 'crank01'
    (case / 'constant').mkdir(parents=True)
    (case / 'fullCaseSetupDict').write_text(SETUP + '\n[MESH_MOTION]\nDRIVE = crank\nLOAD = generator\nINERTIA = 0.02\n')
    (case / 'constant' / 'dynamicMeshDict').write_text('crank\n{\n        theta0          30;   //deg\n}\n')
    theta0, w0, alpha, A, B, c = (CRANK[k] for k in ['theta0', 'w0', 'alpha', 'A', 'B', 'c'])
    t = np.linspace(0.0, 0.33, 20001)[1:]
    theta = theta0 + np.degrees(w0*t + 0.5*alpha*t**2)
    omega = w0 + alpha*t
    torque = B + A*np.sin(np.radians(theta))
    rows = np.column_stack([t, theta, omega, torque, c*omega, 0.9*c*omega**2, 10*torque, torque])
    header = ['Time', 'theta', 'omega', 'torqueGas', 'torqueLoad', 'powerElec', 'F_piston', 'T_piston']
    half = len(t)//2
    post = case / 'postProcessing' / 'crank'
    _dat(post / '0' / 'crank.dat', header, rows[:half + 3], '# crank')
    _dat(post / ('%.12g' % t[half]) / 'crank.dat', header, rows[half:], '# crank')
    return case


def test_summary_crank_drive(tmp_path):
    case = makeCrankCase(tmp_path)
    w0, alpha, A, B, c = CRANK['w0'], CRANK['alpha'], CRANK['A'], CRANK['B'], CRANK['c']
    s = postInternal.buildSummary(str(case), avgCycles=1)
    #revolutions in 0.33 s: (w0 t + alpha t^2/2)/(2 pi) = 3.57
    assert s['DRIVE'] == 'crank' and s['LOAD'] == 'generator' and s['CYCLES_DONE'] == 3
    assert s['WORK_piston'] == pytest.approx(2*math.pi*B, rel=1e-4)
    assert s['AVG_ENERGY'] == pytest.approx(2*math.pi*B, rel=1e-4)
    assert s['FORCE_Z_AMP_piston'] == pytest.approx(10*A, rel=1e-3)
    #third cycle: phase 4 pi .. 6 pi
    tc = [(-w0 + math.sqrt(w0**2 + 2*alpha*a))/alpha for a in (4*math.pi, 6*math.pi)]
    assert s['FREQ_MEAN'] == pytest.approx(1/(tc[1] - tc[0]), rel=1e-4)
    assert s['AVG_POWER'] == pytest.approx(2*math.pi*B/(tc[1] - tc[0]), rel=1e-4)
    assert s['FREQ_MIN'] == pytest.approx((w0 + alpha*tc[0])/(2*math.pi), rel=1e-3)
    assert s['FREQ_MAX'] == pytest.approx((w0 + alpha*tc[1])/(2*math.pi), rel=1e-3)
    assert s['FREQ_CHANGE'] > 0
    assert s['TORQUE_GAS_MEAN'] == pytest.approx(B, rel=1e-3)
    wMeanSq = (((w0 + alpha*tc[1])**3 - (w0 + alpha*tc[0])**3)/(3*alpha))/(tc[1] - tc[0])
    assert s['POWER_ELEC'] == pytest.approx(0.9*c*wMeanSq, rel=1e-4)
    assert s['POWER_LOAD'] == pytest.approx(c*wMeanSq, rel=1e-4)


#----- plots

def test_cycle_series_and_plots(tmp_path):
    case = makeCase(tmp_path)
    series = postInternal.cycleSeries(str(case))
    #every cycle of the damper does the same work, as WORK_piston in the summary
    w = 2*math.pi*FREQ
    work = -DAMP*AMP**2*w**2/(2*FREQ)
    assert list(series['cycle']) == [1, 2]
    assert series['work']['piston'] == pytest.approx([work, work], rel=2e-3)
    assert series['T'] == pytest.approx([300.0, 300.0])
    assert series['massDrift'] == pytest.approx([0.0, 0.0], abs=1e-12)
    #heat per cycle = mean heat rate / FREQ
    assert series['heat']['hotEnd'] == pytest.approx([15.0, 15.0], rel=1e-3)
    assert series['heat']['WALL-chamber_default'] == pytest.approx([-2.0, -2.0], rel=1e-3)
    assert 'freq' not in series
    #p-V: one closed loop per cycle, area pi*P_AMP*V_AMP
    loops = postInternal.pvLoops(str(case), nLoops=3)
    assert [loop[0] for loop in loops] == [1, 2]
    assert postInternal.integrate(loops[-1][1], loops[-1][2]) == pytest.approx(math.pi*P_AMP*V_AMP, rel=2e-3)
    written = postInternal.writePlots(str(case))
    assert sorted(os.path.basename(p) for p in written) == ['cycleConvergence.png', 'pV.png']
    assert (case / 'pV.png').stat().st_size > 1000 and (case / 'cycleConvergence.png').stat().st_size > 1000


def test_cycle_series_crank_drive(tmp_path):
    case = makeCrankCase(tmp_path)
    series = postInternal.cycleSeries(str(case))
    w0, alpha = CRANK['w0'], CRANK['alpha']
    #work per cycle = integral of the torque over one turn = 2 pi B, speed per cycle from the cycle times
    assert series['work']['piston'] == pytest.approx([2*math.pi*CRANK['B']]*3, rel=1e-4)
    tc = [(-w0 + math.sqrt(w0**2 + 2*alpha*a))/alpha for a in (0.0, 2*math.pi, 4*math.pi, 6*math.pi)]
    assert series['freq'] == pytest.approx(1/np.diff(tc), rel=1e-4)
    #no gas output in this log: only the convergence plot, with work and speed
    assert postInternal.pvLoops(str(case)) == []
    written = postInternal.writePlots(str(case))
    assert [os.path.basename(p) for p in written] == ['cycleConvergence.png']


def test_not_internal_case(tmp_path):
    (tmp_path / 'fullCaseSetupDict').write_text('[TITLES]\nCASENAME = x\n')
    with pytest.raises(SystemExit):
        postInternal.readCaseSetup(str(tmp_path))


#----- gsheet update logic

def test_plan_empty_sheet_writes_fixed_header():
    headers, row, updates, added = gsheetInternal.planUpdate([], [], {'JOB': 'J', 'NAME': 'a', 'AVG_POWER': 5.0, 'X_NEW': 1})
    assert headers[:12] == postInternal.FIXED_KEYS and headers[-1] == 'X_NEW'
    assert row is None and updates[headers.index('AVG_POWER')] == 5.0


def test_plan_existing_row_new_column_and_blanks():
    headers = ['JOB', 'NAME', 'DESCRIPTION', 'COMMENTS', 'SOLVE_TIME ', 'AVG_POWER', '']
    rows = [['J', 'a', 'old', 'my note', '1', '2'], ['J', 'b', '', '', '', '']]
    summary = {'JOB': 'J', 'NAME': 'b', 'DESCRIPTION': '', 'SOLVE_TIME': 3.5, 'AVG_POWER': float('nan'),
               'P_crown_AMP': 1200.0, 'COMMENTS': 'ignored'}
    new, row, updates, added = gsheetInternal.planUpdate(headers, rows, summary)
    assert row == 3 and added == ['P_crown_AMP']
    #trailing-space header still matches, empty and nan values are not written, COMMENTS untouched
    assert updates == {1: 'b', 0: 'J', 4: 3.5, 6: 1200.0}
    new, row, updates, added = gsheetInternal.planUpdate(headers, rows, summary, comment='rerun')
    assert updates[3] == 'rerun'


class FakeWorksheet:
    #enough of a gspread worksheet for pushSummary
    def __init__(self, values):
        self.values = [list(r) for r in values]
        self.title = 'runs'
        self.col_count, self.row_count = 26, 100

    def get_all_values(self):
        return [list(r) for r in self.values]

    def add_cols(self, n):
        self.col_count += n

    def add_rows(self, n):
        self.row_count += n

    def batch_update(self, batch, value_input_option='RAW'):
        import re
        for item in batch:
            m = re.match(r'([A-Z]+)(\d+)', item['range'])
            col = 0
            for ch in m.group(1):
                col = col*26 + ord(ch) - 64
            row = int(m.group(2))
            for k, value in enumerate(item['values'][0]):
                while len(self.values) < row:
                    self.values.append([])
                line = self.values[row - 1]
                while len(line) < col + k:
                    line.append('')
                line[col - 1 + k] = value


def test_push_appends_then_updates():
    ws = FakeWorksheet([['JOB', 'NAME', 'AVG_POWER'], ['J', 'a', 1.0]])
    gsheetInternal.pushSummary(ws, {'JOB': 'J', 'NAME': 'b', 'AVG_POWER': 2.0, 'CYCLES_DONE': 3})
    assert ws.values[0] == ['JOB', 'NAME', 'AVG_POWER', 'CYCLES_DONE']
    assert ws.values[2] == ['J', 'b', 2.0, 3]
    gsheetInternal.pushSummary(ws, {'JOB': 'J', 'NAME': 'a', 'AVG_POWER': 9.0})
    assert ws.values[1][:3] == ['J', 'a', 9.0] and len(ws.values) == 3
