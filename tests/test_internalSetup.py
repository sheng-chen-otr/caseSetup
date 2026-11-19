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

#unit tests for the internal flow tool (no OpenFOAM needed)
import gzip
import os
import re

import numpy as np
import pytest

import geometryIO
import setupCore
import writeInternal
import internalSetup


STL = """solid hot
 facet normal 0 0 -1
  outer loop
   vertex 0 0 0
   vertex 1 0 0
   vertex 1 1 0
  endloop
 endfacet
endsolid hot
solid cold
 facet normal 0 0 1
  outer loop
   vertex 0 0 1
   vertex 1 0 1
   vertex 1 1 1
  endloop
 endfacet
endsolid cold
solid sideWall
 facet normal 1 0 0
  outer loop
   vertex 1 0 0
   vertex 1 1 0
   vertex 1 1 1
  endloop
 endfacet
endsolid sideWall
"""


def wall(**kw):
    #WALL section as the config reader gives it (values are lists of strings)
    section = {'U_DEFAULT': ['noSlip'], 'T_DEFAULT': ['adiabatic'], 'BC_GROUPS': ['']}
    for key, val in kw.items():
        section[key] = val.split(' ')
    return section


#----- PID reader

def test_readSurfacePIDs_stl(tmp_path):
    path = tmp_path / 'WALL-a.stl'
    path.write_text(STL)
    assert geometryIO.readSurfacePIDs(str(path)) == ['hot', 'cold', 'sideWall']


def test_readSurfacePIDs_stl_gz(tmp_path):
    path = tmp_path / 'WALL-a.stl.gz'
    with gzip.open(str(path), 'wt') as f:
        f.write(STL)
    assert geometryIO.readSurfacePIDs(str(path)) == ['hot', 'cold', 'sideWall']


def test_readSurfacePIDs_obj_groups_and_objects(tmp_path):
    path = tmp_path / 'WALL-a.obj'
    path.write_text('v 0 0 0\nv 1 0 0\nv 1 1 0\ng crown\nf 1 2 3\ng liner\nf 1 2 3\ng crown\nf 1 2 3\n')
    assert geometryIO.readSurfacePIDs(str(path)) == ['crown', 'liner']
    path.write_text('v 0 0 0\nv 1 0 0\nv 1 1 0\no bodyA\nf 1 2 3\n')
    assert geometryIO.readSurfacePIDs(str(path)) == ['bodyA']


#----- BC groups

def test_parseGroupLine_full():
    group = writeInternal.parseGroupLine('PIDS hot ring U slip T fixed 600'.split(' '), 'g')
    assert group == {'keywords': ['hot', 'ring'], 'U': 'slip', 'T': ['fixed', 600.0]}


@pytest.mark.parametrize('line', ['hot T fixed 600', 'PIDS', 'PIDS hot U moving',
                                  'PIDS hot T fixed', 'PIDS hot T convective 15',
                                  'PIDS hot T fixed hot', 'PIDS hot T cold 3'])
def test_parseGroupLine_errors(line):
    with pytest.raises(SystemExit):
        writeInternal.parseGroupLine(line.split(' '), 'g')


def test_resolveBCGroups_groups_and_default():
    setup = wall(BC_GROUPS='hotEnd coldEnd', hotEnd='PIDS hot T fixed 600',
                 coldEnd='PIDS COLD T convective 15 300')
    groups = writeInternal.resolveBCGroups('WALL-a', setup, ['hot', 'cold', 'sideWall', 'liner'])
    assert [g['name'] for g in groups] == ['hotEnd', 'coldEnd', 'WALL-a_default']
    assert groups[0]['pids'] == ['hot'] and groups[0]['T'] == ['fixed', 600.0] and groups[0]['U'] == 'noSlip'
    assert groups[1]['pids'] == ['cold'] and groups[1]['T'] == ['convective', 15.0, 300.0]
    assert groups[2]['pids'] == ['sideWall', 'liner'] and groups[2]['T'] == ['adiabatic']


def test_resolveBCGroups_no_default_when_all_matched():
    setup = wall(BC_GROUPS='all', all='PIDS o i T fixed 300')
    groups = writeInternal.resolveBCGroups('WALL-a', setup, ['hot', 'cold', 'liner'])
    assert [g['name'] for g in groups] == ['all']


def test_resolveBCGroups_keyword_matches_nothing():
    setup = wall(BC_GROUPS='g', g='PIDS crown')
    with pytest.raises(SystemExit):
        writeInternal.resolveBCGroups('WALL-a', setup, ['hot', 'cold'])


def test_resolveBCGroups_pid_in_two_groups():
    setup = wall(BC_GROUPS='a b', a='PIDS hot', b='PIDS ho')
    with pytest.raises(SystemExit):
        writeInternal.resolveBCGroups('WALL-a', setup, ['hot', 'cold'])


def test_resolveBCGroups_missing_group_line():
    setup = wall(BC_GROUPS='a')
    with pytest.raises(SystemExit):
        writeInternal.resolveBCGroups('WALL-a', setup, ['hot'])


#----- materials, mesh, times

def _material(**kw):
    material = {'GAS': ['air'], 'TRANSPORT': ['const'], 'MOL_WEIGHT': ['default'], 'CP': ['default'],
                'MU': ['default'], 'PR': ['default'], 'SUTHERLAND_AS': ['default'], 'SUTHERLAND_TS': ['default']}
    for key, val in kw.items():
        material[key] = [val]
    return {'GLOBAL_MATERIAL': material}


def test_gas_presets_and_overrides():
    props = writeInternal.getGasProperties(_material())
    assert props['MOL_WEIGHT'] == 28.96 and props['TRANSPORT'] == 'const'
    #sutherland As for air is the textbook 1.458e-6
    assert props['AS'] == pytest.approx(1.458e-6, rel=5e-3)
    props = writeInternal.getGasProperties(_material(GAS='helium', CP='5200', TRANSPORT='sutherland'))
    assert props['CP'] == 5200.0 and props['MOL_WEIGHT'] == 4.003 and props['TRANSPORT'] == 'sutherland'


def test_gas_errors():
    with pytest.raises(SystemExit):
        writeInternal.getGasProperties(_material(GAS='argonish'))
    with pytest.raises(SystemExit):
        writeInternal.getGasProperties(_material(GAS='custom'))
    with pytest.raises(SystemExit):
        writeInternal.getGasProperties(_material(TRANSPORT='polynomial'))


def test_background_mesh_whole_cells_and_margin():
    lo, hi, n = writeInternal.getBackgroundMesh([-0.05, -0.05, 0.0, 0.05, 0.05, 0.06], 0.004)
    assert np.allclose(lo, [-0.06, -0.06, -0.008]) and np.allclose(hi, [0.06, 0.06, 0.068])
    assert list(n) == [30, 30, 19]
    assert np.allclose(lo + n*0.004, hi)


def test_run_times():
    full = {'GLOBAL_CONTROL': {'FREQ': ['10'], 'N_CYCLES': ['3'], 'WRITE_PER_CYCLE': ['20']}}
    endTime, writeInterval = writeInternal.getRunTimes(full)
    assert endTime == pytest.approx(0.3) and writeInterval == pytest.approx(0.005)


def test_turbulence_aliases():
    full = {'GLOBAL_SIM_CONTROL': {'TURB_MODEL': ['kosst']}}
    assert writeInternal.getTurbulence(full) == 'kOmegaSST'
    full['GLOBAL_SIM_CONTROL']['TURB_MODEL'] = ['spalart']
    with pytest.raises(SystemExit):
        writeInternal.getTurbulence(full)


#----- config plumbing

def test_fillAddonSection_keepExtra_keeps_group_lines():
    w = {'WALL-a': {'BC_GROUPS': ['hotEnd'], 'hotEnd': ['PIDS', 'hot']}}
    flag = setupCore.fillAddonSection(w, 'WALL-a', [('U_DEFAULT', 'noSlip'), ('BC_GROUPS', '')], keepExtra=True)
    assert flag is True
    assert list(w['WALL-a']) == ['U_DEFAULT', 'BC_GROUPS', 'hotEnd']
    assert w['WALL-a']['hotEnd'] == ['PIDS', 'hot']


def test_internal_templates_load(repo_root):
    defaults = setupCore.getGlobalDefaults(os.path.join(repo_root, 'setupTemplates', 'internal'))
    for section in ['DEFAULTS', 'TITLES', 'GEOMETRY', 'GLOBAL_CONTROL', 'GLOBAL_SIM_CONTROL',
                    'GLOBAL_MATERIAL', 'GLOBAL_REFINEMENT', 'GLOBAL_COMPUTE_SETUP']:
        assert section in defaults
    assert defaults['GLOBAL_SIM_CONTROL']['TURB_MODEL'] == ['laminar']
    assert defaults['GLOBAL_MATERIAL']['TRANSPORT'] == ['const']


def test_findFromSections_keeps_valid_from_drops_orphan():
    #MOVE-piston: FROM = WALL-chamber, still in GEOM -> kept
    #MOVE-old: leftover section from a removed part, FROM = self -> not kept (cleanUpCaseSetup deletes it)
    w = {'MOVE-piston': {'FROM': ['WALL-chamber']}, 'MOVE-old': {'FROM': ['self']}, 'WALL-chamber': {}}
    assert internalSetup.findFromSections(w, ['WALL-chamber']) == ['MOVE-piston']


def test_geom_prefix_check():
    internalSetup.checkGeomPrefixes({'WALL-chamber.stl': {}})
    with pytest.raises(SystemExit):
        internalSetup.checkGeomPrefixes({'chamber.stl': {}})
    internalSetup.checkGeomPrefixes({'MOVE-displacer.stl': {}})
    with pytest.raises(SystemExit):
        internalSetup.checkGeomPrefixes({'REGN-matrix.stl': {}})


#----- scripts

def _writeScripts(tmp_path, monkeypatch, repo_root, motionZones):
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    case = tmp_path / 'engineRun01'
    case.mkdir()
    monkeypatch.chdir(case)
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', motionZones)
    return case


def test_scripts_slurm_header_and_placeholders(tmp_path, monkeypatch, repo_root):
    case = _writeScripts(tmp_path, monkeypatch, repo_root, ['motion_a', 'motion_b'])
    for name in ['runScript', 'motionCheckScript', 'submitScript']:
        text = (case / name).read_bytes().decode()
        assert '\r' not in text and text.startswith('#!/bin/bash\n')
        assert re.search(r'<[A-Z_]+>', text) is None
    run = (case / 'runScript').read_text()
    assert '#SBATCH --account=def-brendan8' in run and '#SBATCH --ntasks=192' in run
    assert '#SBATCH --job-name=engineRun01\n' in run and '#SBATCH --time=0-12:00' in run
    assert 'module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 openfoam/v2412' in run
    assert 'NP=4' in run and '$RUN snappyHexMesh' in run
    assert "reorderFaceZones.py' motion_a motion_b" in run
    #mesh, solve, export in order, the mesh is skipped once done
    assert run.index('snappyHexMesh') < run.index('createZeroDirectory') < run.index('foamToEnsight')
    assert '.meshDone' in run and '--meshOnly' in run
    motion = (case / 'motionCheckScript').read_text()
    assert '#SBATCH --job-name=engineRun01_motion' in motion and 'SLURM_NTASKS != $NP' in motion
    submit = (case / 'submitScript').read_text()
    assert 'runScript --meshOnly' in submit and 'motionCheckScript' in submit


def test_scripts_no_motion_and_no_export(tmp_path, monkeypatch, repo_root):
    case = _writeScripts(tmp_path, monkeypatch, repo_root, None)
    assert not (case / 'motionCheckScript').exists()
    run = (case / 'runScript').read_text()
    assert 'topoSet' not in run and 'pointDisplacement' not in run and 'foamToEnsight' in run
    assert 'motionCheckScript' not in (case / 'submitScript').read_text()
    #EXPORT none drops the export step, old separate scripts are removed
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    full['GLOBAL_COMPUTE_SETUP']['EXPORT'] = ['none']
    (case / 'meshingScript').write_text('old')
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', None)
    assert 'foamToEnsight' not in (case / 'runScript').read_text()
    assert not (case / 'meshingScript').exists()


#----- incompressible (pimpleFoam) option

def _fullDefaults(repo_root, solver):
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    full['GLOBAL_SIM_CONTROL']['SOLVER'] = [solver]
    full['GLOBAL_REFINEMENT']['LOC_IN_MESH'] = ['0.001', '0.002', '0.003']
    return templateLoc, full


def test_flow_type_and_nu(repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'pimpleFoam')
    assert writeInternal.getFlow(full) == 'incompressible'
    #air at 1 bar, 300 K: rho 1.161, nu = mu/rho
    assert writeInternal.getGasDensity(full) == pytest.approx(1.161, rel=1e-3)
    assert writeInternal.getNu(full) == pytest.approx(writeInternal.getGasProperties(full)['MU']/1.161, rel=1e-3)
    full['GLOBAL_MATERIAL']['NU'] = ['1.5e-5']
    assert writeInternal.getNu(full) == 1.5e-5
    full['GLOBAL_SIM_CONTROL']['SOLVER'] = ['simpleFoam']
    with pytest.raises(SystemExit):
        writeInternal.getFlow(full)


def test_incompressible_dicts(tmp_path, monkeypatch, repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'pimpleFoam')
    monkeypatch.chdir(tmp_path)
    writeInternal.writeNumericsInternal(templateLoc, full)
    writeInternal.writeThermoInternal(templateLoc, full)
    writeInternal.writeFunctionObjectsInternal({}, full)
    schemes = (tmp_path / 'system/fvSchemes').read_text()
    assert 'div((nuEff*dev2(T(grad(U)))))' in schemes and 'meshPhi' not in schemes and 'div(phi,e)' not in schemes
    assert 'pRefPoint ( 0.001 0.002 0.003 );' in (tmp_path / 'system/fvSolution').read_text()
    assert 'nu ' in (tmp_path / 'constant/transportProperties').read_text()
    assert not (tmp_path / 'constant/thermophysicalProperties').exists()
    fos = (tmp_path / 'system/internalFunctionObjects').read_text()
    assert 'fields          (p);' in fos and 'gasMass' not in fos and 'wallHeatFlux' not in fos


def test_compressible_dicts_unchanged(tmp_path, monkeypatch, repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    monkeypatch.chdir(tmp_path)
    writeInternal.writeNumericsInternal(templateLoc, full)
    writeInternal.writeFunctionObjectsInternal({}, full)
    schemes = (tmp_path / 'system/fvSchemes').read_text()
    assert 'div(meshPhi,p)' in schemes and 'div(((rho*nuEff)*dev2(T(grad(U)))))' in schemes
    assert 'pRef' not in (tmp_path / 'system/fvSolution').read_text()
    fos = (tmp_path / 'system/internalFunctionObjects').read_text()
    assert 'gasMass' in fos and 'fields          (p T rho);' in fos
    #wallHeatFlux: .dat every step, the field only at the write times
    heat = fos[fos.index('wallHeatFlux\n{'):]
    heat = heat[:heat.index('}')]
    assert 'executeControl  timeStep;' in heat and 'writeControl    writeTime;' in heat


#----- motion file path

def test_motion_file_path(tmp_path, monkeypatch):
    case = tmp_path / 'job' / 'CASES' / 'run1'
    case.mkdir(parents=True)
    ref = tmp_path / 'job' / '02_reference'
    ref.mkdir()
    (ref / 'engine motion.csv').write_text('theta,d\n0,0\n')
    (case / 'local.csv').write_text('theta,d\n0,0\n')
    get = writeInternal.getMotionFilePath
    #relative to the case directory, path with a space (config values are split on spaces)
    assert get({'MOTION_FILE': ['../../02_reference/engine', 'motion.csv']}, 'MOVE-d', str(case)) \
        == os.path.normpath(str(ref / 'engine motion.csv'))
    assert get({'MOTION_FILE': ['local.csv']}, 'MOVE-d', str(case)) == os.path.normpath(str(case / 'local.csv'))
    #absolute and environment variables
    assert get({'MOTION_FILE': [str(case / 'local.csv')]}, 'MOVE-d', '/elsewhere') == os.path.normpath(str(case / 'local.csv'))
    monkeypatch.setenv('MOTIONDIR', str(ref))
    assert get({'MOTION_FILE': ['$MOTIONDIR/engine', 'motion.csv']}, 'MOVE-d', '/elsewhere') \
        == os.path.normpath(str(ref / 'engine motion.csv'))
    for value in (['required'], [''], ['missing.csv']):
        with pytest.raises(SystemExit):
            get({'MOTION_FILE': value}, 'MOVE-d', str(case))


#----- crank drive

def test_crank_option(repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    assert writeInternal.getCrank(full) is None
    full['MESH_MOTION']['DRIVE'] = ['crank']
    crank = writeInternal.getCrank(full)
    assert crank['LOAD'] == 'generator' and crank['OMEGA0'] == pytest.approx(2*np.pi*10)
    #compressible: the bounce space is at the charge pressure
    assert crank['P_REF'] == float(full['GLOBAL_MATERIAL']['CHARGE_PRESSURE'][0]) and crank['RHO_INF'] == 1.0
    text = writeInternal.crankEntries(crank, 86.0)
    assert 'theta0          86;' in text and 'load            generator;' in text and 'generator { k 1;' in text
    full['MESH_MOTION']['P_BUFFER'] = ['2e5']
    assert writeInternal.getCrank(full)['P_REF'] == 2e5
    #incompressible: kinematic p, forces with the gas density
    full['GLOBAL_SIM_CONTROL']['SOLVER'] = ['pimpleFoam']
    crank = writeInternal.getCrank(full)
    assert crank['P_REF'] == 0.0 and crank['RHO_INF'] == pytest.approx(1.161, rel=1e-3)
    for key, value in (('DRIVE', 'free'), ('LOAD', 'brake'), ('INERTIA', '0'), ('GEN_K', '0')):
        bad = {s: dict(v) for s, v in full.items()}
        bad['MESH_MOTION'][key] = [value]
        with pytest.raises(SystemExit):
            writeInternal.getCrank(bad)
    #fixed speed needs no inertia
    full['MESH_MOTION']['LOAD'] = ['speed']
    full['MESH_MOTION']['INERTIA'] = ['0']
    assert writeInternal.getCrank(full)['LOAD'] == 'speed'


def test_crank_scripts_compile_the_library(tmp_path, monkeypatch, repo_root):
    case = _writeScripts(tmp_path, monkeypatch, repo_root, ['motion_a'])
    assert 'wmake' not in (case / 'runScript').read_text()
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    full['MESH_MOTION']['DRIVE'] = ['crank']
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    run = (case / 'runScript').read_text()
    assert "cd '/opt/caseSetup/foamLibs/stirlingMotion' && wmake libso" in run
    assert run.index('wmake libso') < run.index('createZeroDirectory')
    motion = (case / 'motionCheckScript').read_text()
    assert motion.index('wmake libso') < motion.index('moveDynamicMesh -parallel')
