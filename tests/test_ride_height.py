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

"""Tests for the ride height tool: hardpoint cfg parsing, part matching, wheel corner
assignment, INIT_RH, heave, child setup command, cornering yaw, and the kinematic solver
(pushrod mount, camber, failure reporting, pushrod transform consistency)."""
import math
import sys

import numpy as np
import pandas as pd
import pytest

import utilities
from rideHeightUtils.doubleWishbonePushrodKinematics import (
    DoubleWishbonePushrodSolver,
    RESIDUAL_TOL,
    resolve_pushrod_mount,
)


# one corner of real-looking hardpoints (mm), left side (y < 0). pushrod on the UCA (pullrod)
FL_POINTS = """\
WHEEL_CENTER_STATIC = 0 -619.1 201.78
UCA_F_INNER = -172.4 -246.77 307.4
UCA_R_INNER = 131.3 -245.8 307.5
LCA_F_INNER = -124.4 -209.1 137.3
LCA_R_INNER = 142.3 -209.3 153.7
TIE_INNER = 80 -213.51 194.4
UCA_OUTER_STATIC = 17.8 -549.5 283.6
LCA_OUTER_STATIC = 4.32 -556 118.4
TIE_OUTER_STATIC = 78.8155 -525.2 163.6
PUSHROD_OUTER_STATIC = 18.95 -484.44 339.61
ROCKER_PIVOT = 18.855 -233.7 507.5
ROCKER_AXIS = 1.0 0.0 0.0
ROCKER_PUSHROD_JOINT_REF = 19.49 -273.77 544.15
ROCKER_DAMPER_JOINT_REF = 18.8 -240.6 562.9
ROCKER_DAMPER_CHASSIS = 19.635 -46.9 565.3
WHEEL_AXIS_LOCAL = 0.0 1.0 0.0
WHEEL_FORWARD_LOCAL = 1.0 0.0 0.0
"""


def _mirror(points_text):
    # mirror a corner block across y = 0 (FL -> FR)
    out = []
    for line in points_text.splitlines():
        key, vals = line.split('=')
        x, y, z = [float(v) for v in vals.split()]
        out.append('%s= %s %s %s' % (key, x, -y, z))
    return '\n'.join(out) + '\n'


def _write_cfg(tmp_path, fl_pids, fr_pids='', extra=''):
    text = '[SUSP_FL]\n' + FL_POINTS + fl_pids + '\n[SUSP_FR]\n' + _mirror(FL_POINTS) + fr_pids + '\n' + extra
    p = tmp_path / 'hp.cfg'
    p.write_text(text)
    return str(p)


def _hp(points_text=FL_POINTS, **overrides):
    # hardpoint dict in metres for the solver
    hp = {}
    rocker = {}
    rockerMap = {'ROCKER_PIVOT': 'pivot', 'ROCKER_AXIS': 'axis', 'ROCKER_PUSHROD_JOINT_REF': 'pushrod_joint_ref',
                 'ROCKER_DAMPER_JOINT_REF': 'damper_joint_ref', 'ROCKER_DAMPER_CHASSIS': 'damper_chassis'}
    for line in points_text.splitlines():
        key, vals = [s.strip() for s in line.split('=')]
        vec = [float(v) for v in vals.split()]
        if key in ('ROCKER_AXIS', 'WHEEL_AXIS_LOCAL', 'WHEEL_FORWARD_LOCAL'):
            scaled = vec
        else:
            scaled = [v * 0.001 for v in vec]
        if key in rockerMap:
            rocker[rockerMap[key]] = scaled
        else:
            hp[key.lower()] = scaled
    hp['rocker'] = rocker
    hp.update(overrides)
    return hp


# ---------------------------------------------------------------------------
# step 1: hardpoint cfg parsing
# ---------------------------------------------------------------------------
class TestHardpointCfgParsing:
    def test_key_decides_component(self, tmp_path):
        path = _write_cfg(tmp_path, 'FL_PUSHROD_PID = fl_link_a\nFL_UPRIGHT_PID = flupright\n')
        _corners, kw = utilities._load_suspension_hardpoints_cfg(path)
        assert kw['fl']['PUSHROD'] == ['fl_link_a']
        assert kw['fl']['UPRIGHT'] == ['flupright']  # contains "pr" but is not a pushrod

    def test_invalid_component_key_errors(self, tmp_path):
        path = _write_cfg(tmp_path, 'FL_PUSHRD_PID = flprod\n')
        with pytest.raises(ValueError, match='not a valid part key'):
            utilities._load_suspension_hardpoints_cfg(path)

    def test_wrong_corner_in_key_errors(self, tmp_path):
        path = _write_cfg(tmp_path, 'FR_UCA_PID = fruca\n')
        with pytest.raises(ValueError, match='belongs to corner FL'):
            utilities._load_suspension_hardpoints_cfg(path)

    def test_duplicate_keyword_errors(self, tmp_path):
        path = _write_cfg(tmp_path, 'FL_UCA_PID = arm\nFL_LCA_PID = arm\n')
        with pytest.raises(ValueError, match='used by both'):
            utilities._load_suspension_hardpoints_cfg(path)

    def test_extra_section_uses_its_corner(self, tmp_path):
        extra = '[FL_BRAKE]\nCORNER = fl\nUPRIGHT_PID = flbrakeduct\n'
        path = _write_cfg(tmp_path, 'FL_UCA_PID = fluca\n', extra=extra)
        _corners, kw = utilities._load_suspension_hardpoints_cfg(path)
        assert kw['fl']['UPRIGHT'] == ['flbrakeduct']

    def test_extra_section_without_component_errors(self, tmp_path):
        extra = '[FL_BRAKE]\nCORNER = fl\nPID = flbrakeduct\n'
        path = _write_cfg(tmp_path, 'FL_UCA_PID = fluca\n', extra=extra)
        with pytest.raises(ValueError, match='not a valid part key'):
            utilities._load_suspension_hardpoints_cfg(path)


# ---------------------------------------------------------------------------
# step 2: part matching
# ---------------------------------------------------------------------------
class TestPidMatching:
    def test_clash_between_components_errors(self):
        uca = {'translation': [0, 0, 0]}
        lca = {'translation': [0, 0, 1]}
        patterns = {'fluca': uca, 'flarm': lca}
        info = {'fluca': {'corner': 'fl', 'component': 'UCA', 'keyword': 'fluca'},
                'flarm': {'corner': 'fl', 'component': 'LCA', 'keyword': 'flarm'}}
        with pytest.raises(ValueError, match='FL_UCA_PID = fluca'):
            utilities._build_pid_regex_mapping(['fluca_flarm_part'], patterns, info)

    def test_same_component_overlap_is_fine(self):
        uca = {'translation': [0, 0, 0]}
        patterns = {'fluca': uca, 'fluca_bracket': uca}
        mapping, unmatched, pidToPatterns = utilities._build_pid_regex_mapping(['fluca_bracket_1'], patterns)
        assert mapping['fluca_bracket_1'] is uca
        assert set(pidToPatterns['fluca_bracket_1']) == {'fluca', 'fluca_bracket'}
        assert unmatched == []

    def test_matching_is_case_sensitive(self):
        uca = {'translation': [0, 0, 0]}
        mapping, unmatched, _ = utilities._build_pid_regex_mapping(['FLUCA_main'], {'fluca': uca})
        assert mapping == {}
        assert unmatched == ['FLUCA_main']


# ---------------------------------------------------------------------------
# step 3: simple mode wheel corner
# ---------------------------------------------------------------------------
class TestWheelCorner:
    REFCOR = ['1.5', '0', '0']

    @pytest.mark.parametrize('cx,cy,expected', [
        (0.0, -0.6, 'fl'), (0.0, 0.6, 'fr'), (3.0, -0.6, 'rl'), (3.0, 0.6, 'rr')])
    def test_corner_from_position(self, cx, cy, expected):
        verts = [[cx - 0.3, cy - 0.15, 0.0], [cx + 0.3, cy + 0.15, 0.6]]
        assert utilities._assign_wheel_corner(verts, self.REFCOR) == expected

    def test_file_with_two_wheels_errors(self):
        verts = [[-0.3, -0.75, 0.0], [0.3, 0.75, 0.6]]  # both front wheels in one file
        with pytest.raises(ValueError, match='single wheel'):
            utilities._assign_wheel_corner(verts, self.REFCOR)


# ---------------------------------------------------------------------------
# step 4: INIT_RH
# ---------------------------------------------------------------------------
def _rh_table():
    return pd.DataFrame({'point': [1, 2], 'fl': [30.0, 25.0], 'fr': [30.0, 25.0],
                         'rl': [60.0, 70.0], 'rr': [60.0, 70.0]})


class TestInitialRideHeight:
    def test_blank_means_changes(self):
        df, init = utilities._apply_initial_ride_height(_rh_table(), {'RH_UNIT': ['mm'], 'INIT_RH': ['']})
        assert init == [0.0, 0.0, 0.0, 0.0]
        assert df['dz_fl'].tolist() == pytest.approx([0.030, 0.025])

    def test_values_mean_targets(self):
        setup = {'RH_UNIT': ['mm'], 'INIT_RH': ['30', '30', '60', '60']}
        df, _ = utilities._apply_initial_ride_height(_rh_table(), setup)
        assert df['dz_fl'].tolist() == pytest.approx([0.0, -0.005])
        assert df['dz_rr'].tolist() == pytest.approx([0.0, 0.010])
        assert df['fl'].tolist() == pytest.approx([0.030, 0.025])  # targets kept, in m

    def test_metres(self):
        t = _rh_table()
        t[['fl', 'fr', 'rl', 'rr']] = t[['fl', 'fr', 'rl', 'rr']] / 1000.0
        df, _ = utilities._apply_initial_ride_height(t, {'RH_UNIT': ['m'], 'INIT_RH': ['0.03', '0.03', '0.06', '0.06']})
        assert df['dz_rl'].tolist() == pytest.approx([0.0, 0.010])

    def test_wrong_count_errors(self):
        with pytest.raises(SystemExit):
            utilities._apply_initial_ride_height(_rh_table(), {'RH_UNIT': ['mm'], 'INIT_RH': ['30', '30']})

    def test_bad_unit_errors(self):
        with pytest.raises(SystemExit):
            utilities._apply_initial_ride_height(_rh_table(), {'RH_UNIT': ['inch'], 'INIT_RH': ['']})


# ---------------------------------------------------------------------------
# step 5: heave
# ---------------------------------------------------------------------------
class TestHeave:
    def test_zero_change_gives_zero_heave(self):
        pitch, heave = utilities.calculatePitchAngles(0.0, 0.0, 3.6)
        assert pitch == 0.0 and heave == 0.0

    def test_heave_is_mean_change(self):
        _pitch, heave = utilities.calculatePitchAngles(0.01, 0.03, 3.6)
        assert heave == pytest.approx(0.02)


# ---------------------------------------------------------------------------
# step 6: child case command
# ---------------------------------------------------------------------------
def test_child_command_uses_parent_setup_and_python():
    cmd = utilities._child_case_setup_command('f1-2026')
    assert cmd[0] == sys.executable
    assert cmd[-3:] == ['-s', 'f1-2026', '--rideHeightMode']


# ---------------------------------------------------------------------------
# step 6b: executeAll.sh for the child cases
# ---------------------------------------------------------------------------
def _make_children(tmp_path, monkeypatch, points, header):
    parent = tmp_path / '001'
    parent.mkdir()
    for p in points:
        child = parent / ('001_%s' % p)
        child.mkdir()
        for script in utilities.EXECUTE_ALL_SCRIPTS:
            (child / script).write_text(header + 'echo run\n')
    monkeypatch.chdir(parent)
    return parent


def _rh_dict(useDepend):
    return {'RIDE_HEIGHT_SETUP': {'USE_DEPEND': [useDepend]}}


class TestExecuteAll:
    @pytest.mark.parametrize('useDepend, flag', [('True', 'USE_DEPEND=1'), ('False', 'USE_DEPEND=0')])
    def test_slurm_script_follows_use_depend(self, tmp_path, monkeypatch, useDepend, flag):
        parent = _make_children(tmp_path, monkeypatch, [1, 3], '#!/bin/bash\n#SBATCH --nodes=1\n')
        rh = pd.DataFrame({'point': [1, 3]})
        path = utilities.writeExecuteAllScript(rh, _rh_dict(useDepend))
        text = open(path).read()
        assert path == str(parent / utilities.EXECUTE_ALL_NAME)
        assert 'CHILDREN=(001_1 001_3)' in text
        assert 'SCRIPTS=(meshingScript solveScript exportScript)' in text
        assert flag in text
        assert 'sbatch' in text and 'afterany' in text and 'afterok' in text
        assert '\r' not in open(path, 'rb').read().decode()

    def test_plain_bash_scripts_run_in_order(self, tmp_path, monkeypatch):
        _make_children(tmp_path, monkeypatch, [2], '#!/bin/bash\n')
        path = utilities.writeExecuteAllScript(pd.DataFrame({'point': [2]}), _rh_dict('True'))
        text = open(path).read()
        assert 'sbatch' not in text
        assert 'bash "$script"' in text
        assert 'CHILDREN=(001_2)' in text

    def test_skips_children_without_scripts(self, tmp_path, monkeypatch):
        parent = _make_children(tmp_path, monkeypatch, [1], '#!/bin/bash\n#SBATCH -N1\n')
        (parent / '001_2').mkdir()
        path = utilities.writeExecuteAllScript(pd.DataFrame({'point': [1, 2]}), _rh_dict('True'))
        assert 'CHILDREN=(001_1)' in open(path).read()

    def test_no_children_writes_nothing(self, tmp_path, monkeypatch):
        parent = _make_children(tmp_path, monkeypatch, [], '')
        assert utilities.writeExecuteAllScript(pd.DataFrame({'point': [1]}), _rh_dict('True')) is None
        assert not (parent / utilities.EXECUTE_ALL_NAME).exists()


# ---------------------------------------------------------------------------
# step 7: cornering yaw
# ---------------------------------------------------------------------------
def _corner_dict(direction, yaw, radius=30.0):
    return {
        'BC_SETUP': {'REFCOR': ['1.0', '0.0', '0.0'], 'DOMAIN_PITCH': ['0'], 'DOMAIN_ROLL': ['0'],
                     'INLET_MAG': ['30'], 'YAW': [str(yaw)]},
        'CORNERING_SETUP': {'CORNER_RADIUS': [str(radius)], 'CORNER_DIR': [direction],
                            'CORNER_CENTER': ['default']},
    }


def _local_flow_at_refcor(d):
    # relative air velocity at REFCOR in the rotating frame: -(omega x r)
    omega, axis, centre = utilities.corneringFrame(d)
    r = np.array([1.0, 0.0, 0.0]) - centre
    return -np.cross(omega * axis, r)


class TestCorneringYaw:
    @pytest.mark.parametrize('direction', ['left', 'right'])
    def test_zero_yaw_gives_plus_x_flow(self, direction):
        flow = _local_flow_at_refcor(_corner_dict(direction, 0.0))
        assert flow[0] > 0 and flow[1] == pytest.approx(0.0, abs=1e-9)

    @pytest.mark.parametrize('direction', ['left', 'right'])
    def test_yaw_rotates_flow_like_straight_line_yaw(self, direction):
        flow = _local_flow_at_refcor(_corner_dict(direction, 5.0))
        assert math.degrees(math.atan2(flow[1], flow[0])) == pytest.approx(5.0)

    def test_yaw_keeps_radius(self):
        _, _, centre = utilities.corneringFrame(_corner_dict('left', 5.0))
        assert np.linalg.norm(centre - np.array([1.0, 0.0, 0.0])) == pytest.approx(30.0)


# ---------------------------------------------------------------------------
# step 8: solver
# ---------------------------------------------------------------------------
class TestPushrodMount:
    def test_pullrod_on_uca(self):
        assert resolve_pushrod_mount(_hp()) == 'uca'

    def test_on_lca_joint(self):
        hp = _hp()
        hp['pushrod_outer_static'] = list(hp['lca_outer_static'])
        assert resolve_pushrod_mount(hp) == 'lca'

    def test_ambiguous_point_errors(self):
        hp = _hp()
        mid = 0.5 * (np.array(hp['uca_outer_static']) + np.array(hp['lca_outer_static']))
        hp['pushrod_outer_static'] = (mid + np.array([0.0, -0.05, 0.0])).tolist()  # outboard, between arms
        with pytest.raises(ValueError, match="can't tell which wishbone"):
            resolve_pushrod_mount(hp)


class TestCamber:
    def test_axis_direction_does_not_change_camber(self):
        a = DoubleWishbonePushrodSolver(_hp()).solve_for_wheel_travel(0.02)
        b = DoubleWishbonePushrodSolver(_hp(wheel_axis_local=[0.0, -1.0, 0.0])).solve_for_wheel_travel(0.02)
        assert a['camber_deg'] == pytest.approx(b['camber_deg'])

    def test_left_and_right_match(self):
        left = DoubleWishbonePushrodSolver(_hp()).solve_for_wheel_travel(0.02)
        right = DoubleWishbonePushrodSolver(_hp(_mirror(FL_POINTS))).solve_for_wheel_travel(0.02)
        assert left['camber_deg'] == pytest.approx(right['camber_deg'], abs=1e-9)

    def test_negative_when_top_leans_inboard(self):
        # bump (wheel up relative to body) on this geometry gains negative camber
        bump = DoubleWishbonePushrodSolver(_hp()).solve_for_wheel_travel(0.02)
        assert bump['camber_deg'] < 0.0


class TestSolverFailure:
    def test_reachable_travel_converges(self):
        r = DoubleWishbonePushrodSolver(_hp()).solve_for_wheel_travel(-0.02)
        assert r['residual_norm'] <= RESIDUAL_TOL

    def test_unreachable_travel_reports_cause(self):
        s = DoubleWishbonePushrodSolver(_hp())
        r = s.solve_for_wheel_travel(0.6)  # 600 mm of bump is beyond the linkage
        assert r['residual_norm'] > RESIDUAL_TOL
        msg = s.describe_failure(r, 0.6, 0.0, 'none')
        assert '600.000 mm' in msg and 'off by' in msg


# ---------------------------------------------------------------------------
# step 9: geometry pushrod matches solver pushrod
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('mount', ['uca', 'lca'])
def test_pushrod_transform_matches_solver(mount):
    hp = _hp()
    if mount == 'lca':
        hp['pushrod_outer_static'] = list(hp['lca_outer_static'])
    hp['pushrod_mount'] = resolve_pushrod_mount(hp)
    solver = DoubleWishbonePushrodSolver(hp)
    r = solver.solve_for_wheel_travel(0.015)

    row = pd.Series({
        'fl_wc_x': r['wheel_center'][0], 'fl_wc_y': r['wheel_center'][1], 'fl_wc_z': r['wheel_center'][2],
        'fl_rvec_x': r['state'].rvec[0], 'fl_rvec_y': r['state'].rvec[1], 'fl_rvec_z': r['state'].rvec[2],
        'fl_rocker_deg': r['rocker_theta_deg'], 'fl_damper_delta': r['damper_delta'],
        'fl_rack_travel': r['rack_travel'],
    })
    rideHeights = pd.DataFrame([row])
    kinSetup = {'corners': {'fl': hp}}
    t = utilities._compute_component_transforms(kinSetup, 'fl', row, rideHeights)['PUSHROD']

    from rideHeightUtils.doubleWishbonePushrodKinematics import _rodrigues
    outerStatic = np.array(hp['pushrod_outer_static'])
    pivot = np.array(t['pivot'])
    moved = _rodrigues(np.array(t['rotation_rvec'])) @ (outerStatic - pivot) + pivot + np.array(t['translation'])
    assert np.allclose(moved, r['upright_points']['pushrod'], atol=1e-6)
    assert 'UPRIGHT' in utilities._compute_component_transforms(kinSetup, 'fl', row, rideHeights)


def _obj_text(parts):
    """OBJ text with one triangle per (name, offset) part."""
    lines, n = [], 0
    for name, off in parts:
        lines.append('g %s' % name)
        for v in ([0, 0, 0], [0.1, 0, 0], [0, 0.1, 0]):
            lines.append('v %f %f %f' % (v[0] + off[0], v[1] + off[1], v[2] + off[2]))
        lines.append('f %d %d %d' % (n + 1, n + 2, n + 3))
        n += 3
    return '\n'.join(lines) + '\n'


def _run_case(tmp_path, monkeypatch, fl_pids, files, n_points=1, kin=True, pitch=0.0, roll=0.0, heave=0.0):
    """Build a tiny case from {file name: OBJ text} and run transformGeom. Returns the case folder."""
    import matplotlib
    matplotlib.use('Agg')  # the linkage plot must not open a window
    base = tmp_path / 'base'
    tri = base / 'constant' / 'triSurface'
    tri.mkdir(parents=True)
    for name, text in files.items():
        (tri / name).write_bytes(text.encode())
    (base / 'hp.cfg').write_text('[SUSP_FL]\n' + FL_POINTS + fl_pids + '\n[SUSP_FR]\n' + _mirror(FL_POINTS) + '\n')
    monkeypatch.chdir(base)
    d = {'BC_SETUP': {'REFCOR': ['1.7', '0', '0.1'], 'DOMAIN_PITCH': ['0'], 'DOMAIN_ROLL': ['0']},
         'GLOBAL_REFINEMENT': {'TEMPLATE_TYPE': ['snappy']},
         'RIDE_HEIGHT_SETUP': {'USE_KINEMATIC_SOLVER': ['True' if kin else 'False'], 'HARDPOINT_FILE': ['hp.cfg'],
                               'HARDPOINT_SCALE': ['0.001'], 'SUSP_REQUIRE_ALL_PIDS_MATCHED': ['False']}}
    pts = list(range(1, n_points + 1))
    df = pd.DataFrame({'point': pts, 'wheel_fl': [0.02] * n_points, 'wheel_fr': [0.02] * n_points,
                       'wheel_rl': [0.0] * n_points, 'wheel_rr': [0.0] * n_points,
                       'tunnel_pitch': [pitch] * n_points, 'tunnel_roll': [roll] * n_points,
                       'tunnel_heave': [heave] * n_points})
    rh = utilities._append_suspension_kinematics(df, d) if kin else df
    utilities.transformGeom(d, rh, {name: {} for name in files})
    return base


def _run_transform_geom(tmp_path, monkeypatch, fl_pids, parts):
    """One OBJ, one point, kinematic mode. Returns the child car.obj text."""
    base = _run_case(tmp_path, monkeypatch, fl_pids, {'car.obj': _obj_text(parts)})
    return (base / 'base_1' / 'constant' / 'triSurface' / 'car.obj').read_text()


class TestUnmatchedKeywords:
    def test_unmatched_keyword_warns_and_continues(self, tmp_path, monkeypatch, capsys):
        out = _run_transform_geom(tmp_path, monkeypatch, 'FL_UCA_PID = fluca\nFL_TIE_PID = fltie\n',
                                  [('fluca_arm', (0.0, -0.5, 0.28)), ('body', (1, 0, 0.3))])
        log = capsys.readouterr().out
        assert 'WARNING: These hardpoint-file keywords matched no part' in log
        assert 'FL TIE keyword "fltie" (key FL_TIE_PID)' in log
        assert 'NO SUSPENSION PART WILL BE MOVED' not in log
        # the matched UCA moved, the body did not
        assert 'v 0.000000 -0.500000 0.280000' not in out
        assert 'v 1.000000 0.000000 0.300000' in out

    def test_no_keyword_matches_gives_loud_warning(self, tmp_path, monkeypatch, capsys):
        out = _run_transform_geom(tmp_path, monkeypatch, 'FL_UCA_PID = fluca\n',
                                  [('body', (1, 0, 0.3))])
        log = capsys.readouterr().out
        assert 'NO SUSPENSION PART WILL BE MOVED' in log
        assert 'v 1.000000 0.000000 0.300000' in out


# ---------------------------------------------------------------------------
# faster geometry translation: same results, unmoved data kept as is
# ---------------------------------------------------------------------------
def _obj_vertices(text):
    return np.array([[float(t) for t in l.split()[1:4]] for l in text.splitlines() if l.startswith('v ')])


class TestFastGeometry:
    def test_body_file_copied_and_static_lines_kept(self, tmp_path, monkeypatch):
        car = _obj_text([('fluca_arm', (0.0, -0.5, 0.28)), ('bracket', (0.0, -0.3, 0.3))])
        car = car.replace('v 0.000000 -0.300000 0.300000', 'v 0.0 -0.3 0.3')  # unrounded text must survive
        body = _obj_text([('body', (1, 0, 0.3))])
        base = _run_case(tmp_path, monkeypatch, 'FL_UCA_PID = fluca\n', {'car.obj': car, 'body.obj': body},
                         n_points=2)
        for p in (1, 2):
            tri = base / ('base_%d' % p) / 'constant' / 'triSurface'
            assert (tri / 'body.obj').read_bytes() == body.encode()
            out = (tri / 'car.obj').read_text().splitlines()
            src = car.splitlines()
            assert len(out) == len(src)
            changed = [i for i, (x, y) in enumerate(zip(src, out)) if x != y]
            assert changed == [1, 2, 3]  # only the UCA vertex lines
            assert 'v 0.0 -0.3 0.3' in out

    def test_unmatched_warning_printed_once(self, tmp_path, monkeypatch, capsys):
        _run_case(tmp_path, monkeypatch, 'FL_UCA_PID = fluca\nFL_TIE_PID = fltie\n',
                  {'car.obj': _obj_text([('fluca_arm', (0.0, -0.5, 0.28))])}, n_points=3)
        assert capsys.readouterr().out.count('matched no part') == 1

    def test_vertex_shared_by_two_components_errors(self, tmp_path, monkeypatch):
        text = ('g fluca_arm\nv 0 -0.5 0.28\nv 0.1 -0.5 0.28\nv 0 -0.4 0.28\nf 1 2 3\n'
                'g fllca_arm\nv 0 -0.5 0.12\nf 1 2 4\n')
        with pytest.raises(SystemExit, match='shared by parts of different suspension components'):
            _run_case(tmp_path, monkeypatch, 'FL_UCA_PID = fluca\nFL_LCA_PID = fllca\n', {'car.obj': text})

    def test_ground_single_pass_matches_three_steps(self, tmp_path, monkeypatch):
        grnd = _obj_text([('grnd', (1.0, 0.5, 0.0)), ('grnd2', (3.0, -0.4, 0.02))])
        pitch, roll, heave = 0.4, -0.3, 0.007
        base = _run_case(tmp_path, monkeypatch, '', {'GRND_zone.obj': grnd}, kin=False,
                         pitch=pitch, roll=roll, heave=heave)
        out = _obj_vertices((base / 'base_1' / 'constant' / 'triSurface' / 'GRND_zone.obj').read_text())
        # old order: heave, then roll, then pitch about the raised REFCOR
        c = np.array([1.7, 0.0, 0.1 + heave])
        v = _obj_vertices(grnd) + [0, 0, heave]
        v = (v - c) @ utilities._rodrigues_matrix_from_rvec([math.radians(roll), 0, 0]).T + c
        v = (v - c) @ utilities._rodrigues_matrix_from_rvec([0, math.radians(pitch), 0]).T + c
        assert np.allclose(out, v, atol=1e-6)

    @pytest.mark.parametrize('morph', [
        {'mode': 'displacement', 'vector': [0.001, 0.002, 0.003]},
        {'mode': 'axis_scale', 'axis': 'z', 'factor': 1.1, 'origin': [0, 0, 0.2], 'limits': {'min': 0.1, 'max': 0.35}},
        {'mode': 'vector_stroke', 'vector': [0, 1, 0], 'distance': 0.05, 'origin': [0, -0.7, 0.2],
         'target_length': 0.3, 'limits': {'min': 0.0, 'max': 0.4}},
        {'mode': 'vector_scale', 'vector': [0.1, -1, 0], 'factor': 0.9, 'origin': [0, -0.1, 0.5],
         'limits': {'min': 0.0, 'max': 0.3}},
    ])
    def test_block_morph_matches_per_vertex(self, morph):
        rng = np.random.default_rng(3)
        block = rng.random((200, 3)) - [0.5, 0.8, 0.0]
        cfg = {'translation': [0.01, 0.0, -0.02], 'rotation_rvec': [0.02, -0.01, 0.03], 'pivot': [0.0, -0.2, 0.3],
               'morph': morph}
        fn = utilities._build_vertex_transformer(cfg, None)
        expected = np.array([fn(v, i) for i, v in enumerate(block)])
        got = utilities._transform_vertex_block(cfg, None, block, np.arange(len(block)))
        assert np.allclose(got, expected, atol=1e-12)

    @pytest.mark.parametrize('name', ['g.obj', 'g.obj.gz', 'g.stl', 'g.stl.gz', 'bin.stl'])
    def test_preserve_pid_transform_all_formats(self, tmp_path, name):
        import gzip
        import struct
        rng = np.random.default_rng(4)
        tris = rng.random((10, 3, 3))
        path = tmp_path / name
        if name == 'bin.stl':
            buf = bytearray(b'solid x'.ljust(80, b' ')) + struct.pack('<I', len(tris))
            for t in tris:
                buf += struct.pack('<12f', 0, 0, 0, *t.ravel()) + b'\x07\x00'
            path.write_bytes(bytes(buf))
        else:
            if '.obj' in name:
                text = 'g part\n' + ''.join('v %.9f %.9f %.9f\n' % tuple(p) for p in tris.reshape(-1, 3))
                text += ''.join('f %d %d %d\n' % (3 * k + 1, 3 * k + 2, 3 * k + 3) for k in range(len(tris)))
            else:
                text = 'solid part\n' + ''.join(
                    '  facet normal 0 0 0\n    outer loop\n'
                    + ''.join('      vertex %.9e %.9e %.9e\n' % tuple(p) for p in t)
                    + '    endloop\n  endfacet\n' for t in tris) + 'endsolid part\n'
            data = text.encode()
            path.write_bytes(gzip.compress(data) if name.endswith('.gz') else data)
        out = tmp_path / ('out_' + name)
        utilities.transformGeometryPreservePID(str(path), str(out), rotation={'x': 1.0, 'y': 2.0, 'z': -3.0},
                                               translation=[0.1, 0.2, -0.3], scale=1.01,
                                               rotation_rvec=[0.0, 0.05, 0.0], pivot=[1.7, 0.0, 0.0])

        def ref(v):
            v = v * 1.01
            v = utilities._rodrigues_matrix_from_rvec([0, 0.05, 0]) @ (v - [1.7, 0, 0]) + [1.7, 0, 0]
            v = utilities.x_rotation(v, math.radians(1.0))
            v = utilities.y_rotation(v, math.radians(2.0))
            v = utilities.z_rotation(v, math.radians(-3.0))
            return v + [0.1, 0.2, -0.3]
        expected = np.array([ref(p) for p in tris.reshape(-1, 3)]).reshape(-1, 3, 3)
        if name == 'bin.stl':
            raw = out.read_bytes()
            rec = np.frombuffer(raw, dtype=utilities._STL_BINARY_DTYPE, count=len(tris), offset=84)
            assert raw[:84] == bytes(buf[:84]) and (rec['a'] == 7).all()
            got, normals = rec['v'].astype(float), rec['n'].astype(float)
            tol = 1e-6
        else:
            data = out.read_bytes()
            text = (gzip.decompress(data) if name.endswith('.gz') else data).decode()
            key = 'v ' if '.obj' in name else 'vertex'
            got = np.array([[float(t) for t in l.split()[1:4]] for l in text.splitlines()
                            if l.strip().startswith(key)]).reshape(-1, 3, 3)
            normals = np.array([[float(t) for t in l.split()[2:5]] for l in text.splitlines()
                                if l.strip().startswith('facet normal')])
            tol = 2e-6
        assert np.allclose(got, expected, atol=tol)
        if '.stl' in name:
            n = np.cross(expected[:, 1] - expected[:, 0], expected[:, 2] - expected[:, 0])
            assert np.allclose(normals, n / np.linalg.norm(n, axis=1)[:, None], atol=1e-5)
