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

#unit tests for remeshing: segment angles (motionUtilities) and the run loop steps (meshUtilities/remeshRun.py)
import os
import importlib.util

import numpy as np
import pytest

import motionUtilities as mu
import writeInternal


def _engine():
    #hot 0 < displacer 0.013..0.047 < cold 0.06 < crown, piston gap 2..18 mm
    theta = np.arange(0.0, 360.0, 1.0)
    disp = mu.Body('displacer', *mu.periodicProfile(theta, 0.006*np.sin(np.radians(theta)), 'd'), 0.0, 0.0)
    piston = mu.Body('piston', *mu.periodicProfile(theta, 0.008*(1 - np.cos(np.radians(theta))), 'p'), 0.0, 0.0)
    zones = [mu.Zone('hot', 'w_hot', 0.0, 0.01), mu.Zone('cold', 'w_cold', 0.06, 0.01),
             mu.Zone('dBot', 'd_bottom', 0.013, 0.009, disp), mu.Zone('dTop', 'd_top', 0.047, 0.009, disp),
             mu.Zone('crown', 'w_crown', 0.062, 4e-4, piston)]
    return zones


def test_segment_ratio_matches_cycle_ratio():
    zones = _engine()
    full = mu.deformationRatio(zones, 0.0, 1.0)[0]
    #squash/stretch only grows away from the meshing angle, so the whole cycle gives the cycle value
    assert mu.segmentRatio(zones, 0.0, 360.0, 1.0)[0] == pytest.approx(full)
    assert mu.segmentRatio(zones, 0.0, 0.0, 1.0)[0] == 1.0


def test_remesh_angles_cover_the_cycle_within_the_limit():
    zones = _engine()
    segments = mu.remeshAngles(zones, 30.0, 1.0, 2.0)
    assert segments[0][0] == 30.0 and segments[-1][1] == pytest.approx(390.0)
    assert all(a1 == b0 for (a0, a1), (b0, b1) in zip(segments[:-1], segments[1:]))
    for a0, a1 in segments:
        assert mu.segmentRatio(zones, a0, a1, 1.0)[0] <= 2.0 + 1e-12
    #each segment ends where one more step would break the limit (greedy), except the last one
    for a0, a1 in segments[:-1]:
        assert mu.segmentRatio(zones, a0, a1 + 1.0, 1.0)[0] > 2.0
    #a looser limit needs fewer meshes, a limit above the cycle ratio one
    assert len(mu.remeshAngles(zones, 30.0, 1.0, 3.0)) < len(segments)
    assert len(mu.remeshAngles(zones, 30.0, 1.0, 100.0)) == 1


def test_remesh_limit_below_one_step_is_an_error():
    zones = _engine()
    with pytest.raises(SystemExit):
        mu.remeshAngles(zones, 0.0, 1.0, 1.0001)


def test_getRemesh_keys():
    assert writeInternal.getRemesh({'MESH_MOTION': {'REMESH': ['no']}}) is None
    assert writeInternal.getRemesh({'MESH_MOTION': {'REMESH': ['auto'], 'REMESH_DEFORMATION': ['2.5']}}) == 2.5
    with pytest.raises(SystemExit):
        writeInternal.getRemesh({'MESH_MOTION': {'REMESH': ['auto'], 'REMESH_DEFORMATION': ['1']}})
    with pytest.raises(SystemExit):
        writeInternal.getRemesh({'MESH_MOTION': {'REMESH': ['yes']}})


def test_crank_entries_stop_angle():
    crank = {'OMEGA0': 62.8, 'INERTIA': 0.05, 'LOAD': 'speed', 'FRICTION_TORQUE': 0.0, 'FRICTION_COEFF': 0.0,
             'P_REF': 1e5, 'RHO_INF': 1.0}
    assert 'stopAngle' not in writeInternal.crankEntries(crank, 88.0)
    assert 'stopAngle       28;' in writeInternal.crankEntries(crank, 300.0, 388.0)


#----- run loop helper

def _remeshRun(repo_root):
    spec = importlib.util.spec_from_file_location(
        'remeshRun', os.path.join(repo_root, 'meshUtilities', 'remeshRun.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_segment_of_wraps_and_tolerates_round_off(repo_root):
    rr = _remeshRun(repo_root)
    segments = [('seg01', 88.0, 99.0), ('seg02', 99.0, 318.0), ('seg03', 318.0, 88.0)]
    assert rr.segmentOf(88.0, segments) == 'seg01'
    assert rr.segmentOf(98.9, segments) == 'seg01'
    assert rr.segmentOf(99.0, segments) == 'seg02'
    assert rr.segmentOf(98.9999999999, segments) == 'seg02'
    assert rr.segmentOf(10.0, segments) == 'seg03'
    assert rr.segmentOf(360.0 + 50.0, segments) == 'seg03'
    assert rr.segmentOf(360.0 + 88.0, segments) == 'seg01'


def _fakeCase(path, nProcs=2):
    os.makedirs(path / 'system')
    os.makedirs(path / 'constant')
    (path / 'system' / 'remeshSegments').write_text('#test\nseg01 88 99\nseg02 99 88\n')
    (path / 'system' / 'controlDict').write_text('FoamFile {}\nendTime 0.1;\ndeltaT 1e-06;\nadjustTimeStep yes;\n')
    (path / 'system' / 'pointDisplacement').write_text('pointDisplacement\n')
    for seg in ['seg01', 'seg02']:
        (path / 'constant' / ('dynamicMeshDict.%s' % seg)).write_text('%s\n' % seg)
    for k in range(nProcs):
        os.makedirs(path / ('processor%d' % k) / '0')
        meshDir = path / 'meshes' / 'seg02' / ('processor%d' % k) / 'constant' / 'polyMesh'
        os.makedirs(meshDir)
        (meshDir / 'faces').write_text('faces\n')
        (meshDir / 'points').write_bytes(b'FoamFile\n{\n    class       vectorField;\n    object      points;\n}\n'
                                         b'\x00\x01binary')


def _writeTime(path, name, theta, nProcs=2):
    for k in range(nProcs):
        t = path / ('processor%d' % k) / name
        os.makedirs(t / 'uniform' / 'functionObjects')
        os.makedirs(t / 'polyMesh')
        (t / 'polyMesh' / 'points').write_text('moved points\n')
        (t / 'uniform' / 'crankState').write_text('theta %s;\nomega 62.8;\n' % theta)
        (t / 'uniform' / 'time').write_text('value %s;\ndeltaT 0.0001;\ndeltaT0 0.0001;\n' % name)
        for field in ['p', 'U', 'phi', 'meshPhi']:
            (t / field).write_text('old %s\n' % field)


def test_run_loop_steps(tmp_path, repo_root, monkeypatch):
    rr = _remeshRun(repo_root)
    _fakeCase(tmp_path)
    monkeypatch.chdir(tmp_path)
    #start: case mesh, crank angle from the first segment
    assert rr.nextStep() == 'solve seg01 0'
    assert (tmp_path / 'constant' / 'dynamicMeshDict').read_text() == 'seg01\n'
    #the solver did not get past the start: stop instead of looping
    with pytest.raises(SystemExit):
        rr.nextStep()
    #stopped at the end of seg01 (99 deg, a hair short): map onto seg02
    _writeTime(tmp_path, '0.003', '98.99999999999')
    assert rr.nextStep() == 'map seg02 0.003'
    #remeshCover written on the old mesh before mapping (job died before the graft): still the whole map
    for k in range(2):
        (tmp_path / ('processor%d' % k) / '0.003' / 'remeshCover').write_text('old cover\n')
    assert rr.nextStep() == 'map seg02 0.003'
    #mapFieldsPar output in meshes/seg02/processor*/0.003 (remeshCover mapped too), then graft
    for k in range(2):
        mapped = tmp_path / 'meshes' / 'seg02' / ('processor%d' % k) / '0.003'
        os.makedirs(mapped)
        for field in ['p', 'U', 'rho', 'remeshCover']:
            (mapped / field).write_text('mapped %s\n' % field)
    rr.graft('seg02', '0.003', '5e-07')
    t = tmp_path / 'processor1' / '0.003'
    assert sorted(os.listdir(t)) == ['U', 'p', 'pointDisplacement', 'polyMesh', 'remeshCover', 'rho', 'uniform']
    #the fill / mass fix has not run (remeshCover still there): finish the map before solving
    assert rr.nextStep() == 'combine seg02 0.003'
    #remeshMassFix -combine removes remeshCover as its last step
    for k in range(2):
        (tmp_path / ('processor%d' % k) / '0.003' / 'remeshCover').unlink()
    assert sorted(os.listdir(t)) == ['U', 'p', 'pointDisplacement', 'polyMesh', 'rho', 'uniform']
    assert (t / 'p').read_text() == 'mapped p\n'
    assert sorted(os.listdir(t / 'polyMesh')) == ['faces', 'points', 'points0', 'remeshSegment']
    assert b'object      points0;' in (t / 'polyMesh' / 'points0').read_bytes()
    assert (t / 'polyMesh' / 'points0').read_bytes().endswith(b'\x00\x01binary')
    #mesh files hard-linked from the segment mesh (kept for the next cycle), points0 and the marker only in the time
    segMesh = tmp_path / 'meshes' / 'seg02' / 'processor1' / 'constant' / 'polyMesh'
    assert os.path.samefile(t / 'polyMesh' / 'faces', segMesh / 'faces')
    assert sorted(os.listdir(segMesh)) == ['faces', 'points']
    assert 'deltaT 5e-07;' in (t / 'uniform' / 'time').read_text()
    assert (t / 'uniform' / 'crankState').is_file()
    assert not os.path.exists(tmp_path / 'meshes' / 'seg02' / 'processor1' / '0.003')
    assert (tmp_path / 'constant' / 'dynamicMeshDict').read_text() == 'seg02\n'
    #now the time holds seg02: solve on it
    assert rr.nextStep() == 'solve seg02 0.003'
    #a later time still uses the remeshed mesh of 0.003 (moved points only)
    _writeTime(tmp_path, '0.05', '268.0')
    assert rr.meshSegment(str(tmp_path / 'processor0'), 0.05) == 'seg02'
    assert rr.nextStep() == 'solve seg02 0.05'
    #end time reached
    _writeTime(tmp_path, '0.1', '448.0')
    assert rr.nextStep() == 'done'


def test_prepare_keeps_all_digits_of_the_map_time(tmp_path, repo_root, monkeypatch):
    rr = _remeshRun(repo_root)
    os.makedirs(tmp_path / 'meshes' / 'seg02' / 'system')
    (tmp_path / 'meshes' / 'seg02' / 'system' / 'controlDict').write_text(
        'FoamFile {}\nstartFrom       latestTime;\nstartTime       0;\nendTime 0.1;\n')
    monkeypatch.chdir(tmp_path)
    rr.prepare('seg02', '0.003053740817')
    text = (tmp_path / 'meshes' / 'seg02' / 'system' / 'controlDict').read_text()
    assert 'startFrom       startTime;' in text and 'startTime       0.003053740817;' in text
    assert 'timePrecision 15;' in text and 'endTime 0.1;' in text


def test_done_within_one_step_of_end_time(tmp_path, repo_root, monkeypatch):
    #the crank solver writes its last step within one step of endTime: that counts as done, earlier does not
    rr = _remeshRun(repo_root)
    _fakeCase(tmp_path)
    with open(tmp_path / 'system' / 'controlDict', 'a') as f:
        f.write('maxDeltaT 0.0001;\n')
    monkeypatch.chdir(tmp_path)
    _writeTime(tmp_path, '0.0998', '98.0')
    assert rr.nextStep() == 'solve seg01 0.0998'
    _writeTime(tmp_path, '0.09995', '98.5')
    assert rr.nextStep() == 'done'


def test_unfinished_write_step_set_aside(tmp_path, repo_root, monkeypatch, capsys):
    #a solver that crashed in its last (write) step leaves a time with the crank state at the stop angle and no fields:
    #renamed to <time>.incomplete in every processor, the loop goes on from the last finished time
    rr = _remeshRun(repo_root)
    _fakeCase(tmp_path)
    monkeypatch.chdir(tmp_path)
    _writeTime(tmp_path, '0.002', '95.0')
    _writeTime(tmp_path, '0.003', '98.99999999999')
    for k in range(2):
        for field in ['p', 'U', 'phi', 'meshPhi']:
            (tmp_path / ('processor%d' % k) / '0.003' / field).unlink()
    assert rr.nextStep() == 'solve seg01 0.002'
    for k in range(2):
        assert sorted(os.listdir(tmp_path / ('processor%d' % k))) == ['0', '0.002', '0.003.incomplete']
    assert 'time 0.003 had no fields' in capsys.readouterr().err
