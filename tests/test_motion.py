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

#unit tests for the mesh motion helpers (no OpenFOAM needed)
import os
import sys
import importlib.util

import numpy as np
import pytest

import geometryIO
import motionUtilities as mu


def _body(name, amp, clearance=0.0, phase=0.0, kind='sin'):
    theta = np.arange(0.0, 360.0, 1.0)
    if kind == 'sin':
        pos = amp*np.sin(np.radians(theta + phase))
    else:
        pos = amp*(1 - np.cos(np.radians(theta + phase)))
    t, p = mu.periodicProfile(theta, pos, name)
    return mu.Body(name, t, p, 0.0, clearance)


def _engine(pistonClearance=0.002):
    #hot 0 < displacer bottom 0.013 < top 0.047 < cold 0.06 < crown (flush at TDC before clearance)
    disp = _body('displacer', 0.006)
    piston = _body('piston', 0.008, clearance=pistonClearance, kind='cos')
    zones = [mu.Zone('hot', 'w_hot', 0.0, 0.01),
             mu.Zone('cold', 'w_cold', 0.06, 0.01),
             mu.Zone('dBot', 'd_bottom', 0.013, 0.009, disp),
             mu.Zone('dTop', 'd_top', 0.047, 0.009, disp),
             mu.Zone('crown', 'w_crown', 0.06, 4e-4, piston)]
    return zones, [disp, piston]


#----- motion file

def test_readMotionFile_and_periodic(tmp_path):
    path = tmp_path / 'm.csv'
    rows = ['theta,displacer,piston'] + ['%g,%g,%g' % (a, np.sin(np.radians(a)), 0.0) for a in range(0, 361, 10)]
    path.write_text('\n'.join(rows))
    theta, cols = mu.readMotionFile(str(path))
    assert theta[0] == 0 and theta[-1] == 360 and set(cols) == {'displacer', 'piston'}
    t, p = mu.periodicProfile(theta, cols['displacer'], 'displacer')
    assert len(t) == len(theta)
    #wrap-around interpolation
    assert mu.positionAt(370.0, t, p) == pytest.approx(mu.positionAt(10.0, t, p))


def test_periodic_error_and_closing():
    theta = np.array([0.0, 180.0, 360.0])
    with pytest.raises(SystemExit):
        mu.periodicProfile(theta, np.array([0.0, 1.0, 0.5]), 'x')
    t, p = mu.periodicProfile(np.array([0.0, 180.0]), np.array([0.2, 1.0]), 'x')
    assert t[-1] == 360.0 and p[-1] == 0.2


def test_readMotionFile_needs_theta(tmp_path):
    path = tmp_path / 'm.csv'
    path.write_text('angle,a\n0,0\n')
    with pytest.raises(SystemExit):
        mu.readMotionFile(str(path))


#----- zones

def test_flat_zone_and_area():
    tri = np.array([[[0, 0, 0.06], [1, 0, 0.06], [0, 1, 0.06]]], dtype=float)
    assert mu.flatZoneZ('z', tri, 1e-6) == pytest.approx(0.06)
    assert mu.projectedArea(tri) == pytest.approx(0.5)
    tri[0, 2, 2] = 0.061
    with pytest.raises(SystemExit):
        mu.flatZoneZ('z', tri, 1e-4)


def test_clearance_moves_flush_piston_up():
    zones, bodies = _engine(0.002)
    report = mu.applyClearance(zones, bodies, 1.0)
    assert report['piston'][0] == pytest.approx(0.002)
    assert report['piston'][1] == pytest.approx(0.002*4e-4)
    gaps = {(lo.name, up.name): g for lo, up, g in mu.zoneGaps(zones, 1.0)}
    assert gaps[('cold', 'crown')].min() == pytest.approx(0.002)
    mu.checkZoneOrder(zones, 1.0)


def test_zone_crossing_is_an_error():
    disp = _body('displacer', 0.02)
    zones = [mu.Zone('hot', 'h', 0.0, 1.0), mu.Zone('dBot', 'b', 0.013, 1.0, disp)]
    with pytest.raises(SystemExit):
        mu.checkZoneOrder(zones, 1.0)


def test_auto_start_angle_minimises_deformation():
    zones, bodies = _engine(0.002)
    mu.applyClearance(zones, bodies, 1.0)
    best = mu.autoStartAngle(zones, 1.0)
    rBest, pair = mu.deformationRatio(zones, best, 1.0)
    for angle in (0.0, 90.0, 180.0, 270.0):
        assert rBest <= mu.deformationRatio(zones, angle, 1.0)[0] + 1e-12
    #piston gap 2..18 mm: best possible ratio is sqrt(18/2) = 3
    assert rBest == pytest.approx(3.0, rel=0.02)
    assert pair == ('cold', 'crown')


def test_kinematics_sine():
    body = _body('d', 0.006)
    kin = mu.bodyKinematics(body, 10.0, 1.0)
    w = 2*np.pi*10.0
    assert kin['stroke'] == pytest.approx(0.012, rel=1e-3)
    assert kin['vMax'] == pytest.approx(0.006*w, rel=1e-3)
    assert kin['aMax'] == pytest.approx(0.006*w**2, rel=1e-2)


def test_unanchored_horizontal_face_found():
    zones, bodies = _engine(0.002)
    mu.applyClearance(zones, bodies, 1.0)
    #at 180 deg (BDC) the crown is at 0.078: a step at 0.07 is inside the cold face - crown band
    step = np.array([[[0, 0, 0.07], [1, 0, 0.07], [0, 1, 0.07]]], dtype=float)
    wall = np.array([[[0, 0, 0.0], [0, 1, 0.0], [0, 0, 0.06]]], dtype=float)        #vertical
    zoneZs = [float(mu.zoneZ(z, 180.0)) for z in zones]
    found = mu.unanchoredHorizontalFaces({'w_step': step, 'w_side': wall}, zoneZs, zones, 180.0, 1.0, 1e-6)
    assert [f[0] for f in found] == ['w_step']
    #at TDC the crown (0.062) is below the step: not in a band
    zoneZs = [float(mu.zoneZ(z, 0.0)) for z in zones]
    assert mu.unanchoredHorizontalFaces({'w_step': step}, zoneZs, zones, 0.0, 1.0, 1e-6) == []


#----- gap

def test_point_triangle_distance_plane_and_edge():
    tri = np.array([[[0, 0, 0], [1, 0, 0], [0, 1, 0]]], dtype=float)
    pts = np.array([[0.2, 0.2, 0.5], [2.0, 0.0, 0.0], [-1.0, -1.0, 0.0]])
    d = mu.pointTriangleDistance(pts, tri)
    assert d == pytest.approx([0.5, 1.0, np.sqrt(2.0)])


def test_gap_width_and_refinement():
    #two parallel walls 2 mm apart
    wall = np.array([[[0.05, -1, 0], [0.05, 1, 0], [0.05, 0, 1]]], dtype=float)
    body = np.array([[[0.048, -0.5, 0.1], [0.048, 0.5, 0.1], [0.048, 0, 0.5]]], dtype=float)
    gap = mu.gapWidth(body, wall)
    assert gap == pytest.approx(0.002)
    level, distance, cell = mu.gapRefinement(gap, 0.004, 6)
    assert level == 4 and cell == pytest.approx(0.00025) and gap/cell >= 6
    assert distance == pytest.approx(0.0025)


#----- tables

def test_tables_start_at_zero_and_share_times(tmp_path):
    zones, bodies = _engine(0.002)
    names = mu.writeTables(bodies, 90.0, 10.0, 2, 1.0, str(tmp_path))
    assert names == ['zeroTable', 'displacerTable', 'pistonTable']
    rows = {}
    for name in names:
        lines = open(tmp_path / name).read().splitlines()
        assert 'scalarVectorTable' in lines[0]
        assert int(lines[1]) == 721
        rows[name] = [l for l in lines if l.strip().startswith('((') or l.strip().startswith('(0') or l.strip().startswith('(')][1:]
    first = open(tmp_path / 'displacerTable').read().splitlines()[3]
    assert first.strip() == '(0 (0 0 0))'
    last = open(tmp_path / 'displacerTable').read().splitlines()[-2]
    assert float(last.strip()[1:].split()[0]) == pytest.approx(0.2)


def test_crank_tables_one_cycle_zero_at_start(tmp_path):
    zones, bodies = _engine(0.002)
    #0.7 deg does not divide 360: the grid must still end at exactly 360 (the table period)
    names = mu.writeCrankTables(bodies, 90.0, 0.7, str(tmp_path))
    assert names == ['displacerCrank', 'pistonCrank']
    for body, name in zip(bodies, names):
        lines = open(tmp_path / name).read().splitlines()
        rows = np.array([[float(v) for v in l.strip()[1:-1].split()] for l in lines[2:-1]])
        assert int(lines[0]) == len(rows)
        assert rows[0, 0] == 0.0 and rows[-1, 0] == 360.0 and np.all(np.diff(rows[:, 0]) > 0)
        assert rows[0, 1] == pytest.approx(rows[-1, 1])
        assert np.interp(90.0, rows[:, 0], rows[:, 1]) == pytest.approx(0.0, abs=1e-9)
        #same displacement as the time tables: shift(theta) - shift(start)
        assert rows[100, 1] == pytest.approx(float(body.shift(rows[100, 0]) - body.shift(90.0)))


#----- geometry translation

STL2 = """solid crown
 facet normal 0 0 1
  outer loop
   vertex 0 0 0.06
   vertex 1 0 0.06
   vertex 0 1 0.06
  endloop
 endfacet
endsolid crown
solid liner
 facet normal 1 0 0
  outer loop
   vertex 0 0 0.06
   vertex 0 1 0.06
   vertex 0 0 0.1
  endloop
 endfacet
endsolid liner
"""


def test_triangles_and_translate_stl(tmp_path):
    src, out = tmp_path / 'a.stl', tmp_path / 'b.stl'
    src.write_text(STL2)
    tris = geometryIO.readSurfaceTriangles(str(src))
    assert list(tris) == ['crown', 'liner'] and tris['crown'].shape == (1, 3, 3)
    geometryIO.translateSurfacePIDs(str(src), str(out), {'crown': (0, 0, 0.002)})
    moved = geometryIO.readSurfaceTriangles(str(out))
    assert np.allclose(moved['crown'][:, :, 2], 0.062)
    assert np.allclose(moved['liner'], tris['liner'])


def test_translate_obj_shared_vertex_error(tmp_path):
    src, out = tmp_path / 'a.obj', tmp_path / 'b.obj'
    src.write_text('v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\ng crown\nf 1 2 3\ng liner\nf 1 3 4\n')
    with pytest.raises(ValueError):
        geometryIO.translateSurfacePIDs(str(src), str(out), {'crown': (0, 0, 1)})
    #separate shells move fine
    src.write_text('v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\nv 0 1 1\nv 0 0 2\ng crown\nf 1 2 3\ng liner\nf 4 5 6\n')
    geometryIO.translateSurfacePIDs(str(src), str(out), {'crown': (0, 0, 1)})
    moved = geometryIO.readSurfaceTriangles(str(out))
    assert np.allclose(moved['crown'][:, :, 2], 1.0) and np.allclose(moved['liner'][0, :, 2], [1, 1, 2])


#----- faceZone reorder tool

FACEZONES = """FoamFile
{
    version     2.0;
    format      ascii;
    class       regIOobject;
    object      faceZones;
    meta
    {
        names           3 ( REGN coldAnchor hotAnchor );
    }
}

3
(
REGN
{
    type faceZone;
    faceLabels List<label> 2(1 2);
    flipMap List<bool> 2{0};
}
coldAnchor
{
    type faceZone;
    faceLabels List<label> 1(5);
    flipMap List<bool> 1{0};
}
hotAnchor
{
    type faceZone;
    faceLabels List<label> 1(7);
    flipMap List<bool> 1{0};
}
)
"""


def _reorderTool(repo_root):
    spec = importlib.util.spec_from_file_location(
        'reorderFaceZones', os.path.join(repo_root, 'meshUtilities', 'reorderFaceZones.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_reorderFaceZones_serial_and_processors(tmp_path, repo_root):
    tool = _reorderTool(repo_root)
    for sub in ['constant/polyMesh', 'processor0/constant/polyMesh', 'processor1/constant/polyMesh']:
        os.makedirs(tmp_path / sub)
        (tmp_path / sub / 'faceZones').write_text(FACEZONES)
    tool.main(['hotAnchor', 'coldAnchor', '--case', str(tmp_path)])
    for sub in ['constant/polyMesh', 'processor0/constant/polyMesh', 'processor1/constant/polyMesh']:
        txt = (tmp_path / sub / 'faceZones').read_text()
        assert txt.index('hotAnchor\n{') < txt.index('coldAnchor\n{') < txt.index('REGN\n{')
        assert 'names 3 ( hotAnchor coldAnchor REGN )' in txt
        assert 'faceLabels List<label> 1(7);' in txt


def test_reorderFaceZones_missing_zone(tmp_path, repo_root):
    tool = _reorderTool(repo_root)
    os.makedirs(tmp_path / 'constant/polyMesh')
    (tmp_path / 'constant/polyMesh/faceZones').write_text(FACEZONES)
    with pytest.raises(ValueError):
        tool.main(['crown', '--case', str(tmp_path)])


def _binaryZone(name, labels):
    #binary faceZone entry as snappy/topoSet write it (labels 32 bit, bools 1 byte)
    import struct
    if not labels:
        return name.encode() + b'\n{\n    type faceZone;\n    faceLabels List<label> 0;\n    flipMap List<bool> 0;\n}'
    n = len(labels)
    return (name.encode() + b'\n{\n    type faceZone;\n    faceLabels List<label> \n%d\n(' % n
            + struct.pack('<%di' % n, *labels) + b');\n    flipMap List<bool> \n%d\n(' % n
            + bytes([0]*n) + b');\n}')


def test_reorderFaceZones_binary(tmp_path, repo_root):
    #label bytes 123 '{', 125 '}' and 41 ')' inside the binary data must not break the brace matching
    tool = _reorderTool(repo_root)
    header = (b'FoamFile\n{\n    version 2.0;\n    format binary;\n    arch "LSB;label=32;scalar=64";\n'
              b'    class regIOobject;\n    object faceZones;\n    meta\n    {\n'
              b'        names 3 ( REGN coldAnchor hotAnchor );\n    }\n}\n\n3\n(')
    zones = {'REGN': [123, 125, 41], 'coldAnchor': [125], 'hotAnchor': []}
    body = b'\n'.join(_binaryZone(n, zones[n]) for n in ['REGN', 'coldAnchor', 'hotAnchor'])
    os.makedirs(tmp_path / 'constant/polyMesh')
    path = tmp_path / 'constant/polyMesh/faceZones'
    path.write_bytes(header + body + b'\n)\n\n// ***** //\n')
    tool.main(['hotAnchor', 'coldAnchor', '--case', str(tmp_path)])
    out = path.read_bytes()
    assert b'names 3 ( hotAnchor coldAnchor REGN )' in out
    assert out.index(b'hotAnchor\n{') < out.index(b'coldAnchor\n{') < out.index(b'REGN\n{')
    for name, labels in zones.items():
        assert _binaryZone(name, labels) in out
    #reordering back gives the zones in the original order again
    tool.main(['REGN', 'coldAnchor', 'hotAnchor', '--case', str(tmp_path)])
    assert path.read_bytes().endswith(body + b'\n)\n')


#----- clearance by profile and sealing

def _gap(zones, lower, upper):
    return dict(((lo.name, up.name), g) for lo, up, g in mu.zoneGaps(zones, 1.0))[(lower, upper)]


def test_soft_clearance_keeps_the_far_end():
    zones, bodies = _engine(0.002)
    piston = bodies[1]
    piston.mode = 'soft'
    before = piston.pos.copy()
    report = mu.adjustProfiles(zones, bodies, 1.0)['piston']
    #crown stops 2 mm from the cold face, the other end of the stroke is where the profile puts it
    assert _gap(zones, 'cold', 'crown').min() == pytest.approx(0.002, abs=1e-8)
    assert piston.pos.max() == pytest.approx(before.max(), abs=1e-12)
    assert report['raise'] == pytest.approx(0.002, abs=1e-8) and report['lower'] == 0.0
    assert report['strokeAfter'] == pytest.approx(report['strokeBefore'] - 0.002, abs=1e-8)
    #smooth: the local floor at most about doubles the peak acceleration (no velocity jump)
    assert np.abs(np.diff(np.diff(piston.pos))).max() < 2.5*np.abs(np.diff(np.diff(before))).max()
    assert np.abs(np.diff(piston.pos)).max() <= np.abs(np.diff(before)).max() + 1e-12
    #nothing left for the shift step
    assert mu.applyClearance(zones, bodies, 1.0) == {}
    #shift mode is unchanged: the whole path moves
    zones, bodies = _engine(0.002)
    assert mu.adjustProfiles(zones, bodies, 1.0) == {}
    assert mu.applyClearance(zones, bodies, 1.0)['piston'][0] == pytest.approx(0.002)


def test_soft_clearance_both_ends_and_scale():
    for mode in ('soft', 'scale'):
        disp = _body('displacer', 0.0125)
        disp.mode = mode
        disp.clearance = 0.002
        zones = [mu.Zone('hot', 'h', 0.0, 1.0), mu.Zone('cold', 'c', 0.06, 1.0),
                 mu.Zone('dBot', 'b', 0.013, 1.0, disp), mu.Zone('dTop', 't', 0.047, 1.0, disp)]
        before = disp.pos.copy()
        report = mu.adjustProfiles(zones, [disp], 1.0)['displacer']
        assert _gap(zones, 'hot', 'dBot').min() == pytest.approx(0.002, abs=1e-8)
        assert _gap(zones, 'dTop', 'cold').min() == pytest.approx(0.002, abs=1e-8)
        assert report['strokeAfter'] == pytest.approx(report['strokeBefore'] - 0.003, abs=1e-8)
        if mode == 'scale':
            #shape kept exactly
            shape = (disp.pos - disp.pos.min())/(disp.pos.max() - disp.pos.min())
            assert shape == pytest.approx((before - before.min())/(before.max() - before.min()), abs=1e-12)
    #a displacer that leaves no room for the clearance on both sides
    disp = _body('displacer', 0.0005)
    disp.mode, disp.clearance = 'soft', 0.002
    zones = [mu.Zone('hot', 'h', 0.0, 1.0), mu.Zone('cold', 'c', 0.06, 1.0),
             mu.Zone('dBot', 'b', 0.001, 1.0, disp), mu.Zone('dTop', 't', 0.059, 1.0, disp)]
    with pytest.raises(SystemExit):
        mu.adjustProfiles(zones, [disp], 1.0)


def test_soft_clearance_on_a_dwell_profile():
    #displacer holds at the hot face (gap 0) for 120 deg, smooth transitions: the dwell stays flat at the clearance
    theta = np.arange(0.0, 360.0, 0.5)
    x = np.clip((theta - 60.0)/60.0, 0, 1) - np.clip((theta - 240.0)/60.0, 0, 1)
    smooth = x*x*(3 - 2*x)
    t, p = mu.periodicProfile(theta, 0.01*smooth, 'displacer')
    disp = mu.Body('displacer', t, p, 0.0, 0.002, 'soft')
    zones = [mu.Zone('hot', 'h', 0.0, 1.0), mu.Zone('cold', 'c', 0.06, 1.0),
             mu.Zone('dBot', 'b', 0.0, 1.0, disp), mu.Zone('dTop', 't', 0.034, 1.0, disp)]
    mu.adjustProfiles(zones, [disp], 1.0)
    dwell = (t <= 60.0) | (t >= 300.0)
    assert np.ptp(disp.pos[dwell]) < 1e-9 and disp.pos[dwell].mean() == pytest.approx(0.002, abs=1e-8)
    #top dwell (gap to the cold face 16 mm) untouched apart from the tiny far-end correction
    top = (t >= 120.0) & (t <= 240.0)
    assert disp.pos[top] == pytest.approx(0.01, abs=1e-9)


def test_seal_check_crown_in_liner():
    n = 24
    ang = np.linspace(0, 2*np.pi, n + 1)
    ring = np.column_stack([0.05*np.cos(ang), 0.05*np.sin(ang)])

    def disc(r, z):
        return np.array([[[0, 0, z], [r*ring[i, 0]/0.05, r*ring[i, 1]/0.05, z],
                          [r*ring[i + 1, 0]/0.05, r*ring[i + 1, 1]/0.05, z]] for i in range(n)])
    liner = np.array([[[ring[i, 0], ring[i, 1], 0.03], [ring[i + 1, 0], ring[i + 1, 1], 0.03], [ring[i, 0], ring[i, 1], 0.18]]
                      for i in range(n)])
    #the chamber side wall encloses the crown but is not its liner
    chamber = np.array([[[-0.2, -0.2, -0.03], [0.2, -0.2, -0.03], [0.2, -0.2, 0.03]],
                        [[-0.2, 0.2, -0.03], [0.2, 0.2, -0.03], [0.2, 0.2, 0.03]]])
    walls = {'liner': liner, 'chamber': chamber,
             'cold': np.array([[[-0.2, -0.2, 0.03], [0.2, -0.2, 0.03], [0.2, 0.2, 0.03]]])}
    assert mu.sideWallRange(disc(0.0502, 0.18), walls, 1e-6)[2] == ['liner']
    crown = disc(0.0502, 0.18)
    assert mu.sealProblems(crown, 0.1, walls, 1e-6) == []
    assert mu.sealProblems(crown, 0.18, walls, 1e-6) == []
    #2 mm past the liner end
    problems = mu.sealProblems(crown, 0.182, walls, 1e-6)
    assert len(problems) == 1 and 'outside its liner' in problems[0] and '0.002' in problems[0]
    #crown smaller than the bore
    problems = mu.sealProblems(disc(0.048, 0.18), 0.1, walls, 1e-6)
    assert len(problems) == 1 and 'does not cover the bore' in problems[0]
    #no liner around it: nothing to check
    assert mu.sealProblems(crown, 0.5, {'cold': walls['cold']}, 1e-6) == []


def test_soft_floor_is_local():
    #2 mm raise on a 152 mm stroke: exact at both ends, below 0.05 mm from 5 mm on, monotonic
    x = np.linspace(0.0, 0.152, 15201)
    y = mu.softFloor(x, 0.0, 0.002, 0.152)
    assert y[0] == pytest.approx(0.002, abs=1e-15) and y[-1] == pytest.approx(0.152, abs=1e-15)
    assert np.abs(y - x)[x >= 0.005].max() < 5e-5
    assert np.all(np.diff(y) >= -1e-15) and y.min() >= 0.002 - 1e-15
    z = mu.softCeiling(x, 0.152, 0.002, 0.0)
    assert z[-1] == pytest.approx(0.150, abs=1e-15) and z[0] == pytest.approx(0.0, abs=1e-15)
