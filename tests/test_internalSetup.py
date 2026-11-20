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
    assert group == {'keywords': ['hot', 'ring'], 'U': 'slip', 'T': ['fixed', 600.0], 'flow': None}


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


def test_move_pids_in_wall_files_become_sections(repo_root):
    pidDict = {'WALL-chamber-full-system': ['MOVE-piston', 'WALL-piston-wall', 'WALL-chamber-hot', 'move-rod'],
               'MOVE-displacer': ['MOVE-displacer-hot', 'side']}
    #MOVE- PIDs of non-MOVE files only, one body per PID (prefix case-insensitive)
    bodies = internalSetup.findPidBodies(pidDict)
    assert bodies == {'MOVE-piston': ('WALL-chamber-full-system', 'MOVE-piston'),
                      'MOVE-rod': ('WALL-chamber-full-system', 'move-rod')}
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    moveItems = setupCore.readAddonTemplate(templateLoc, 'defaultMove', 'MOVE_SETUP')
    #an existing section is kept, a new one is filled in template order and flags the update
    w = {'MOVE-rod': {'FROM': ['WALL-chamber-full-system'], 'ZONE_PIDS': ['move-rod'], 'MOTION_COLUMN': ['myRod']}}
    assert internalSetup.addPidBodies(w, bodies, ['WALL-chamber-full-system', 'MOVE-displacer'], moveItems)
    assert w['MOVE-rod']['MOTION_COLUMN'] == ['myRod']
    piston = w['MOVE-piston']
    assert list(piston) == [key for key, val in moveItems]
    assert piston['FROM'] == ['WALL-chamber-full-system'] and piston['ZONE_PIDS'] == ['MOVE-piston']
    assert piston['MOTION_COLUMN'] == ['piston'] and piston['TDC_CLEARANCE'] == ['0.002']
    assert not internalSetup.addPidBodies(w, bodies, ['WALL-chamber-full-system'], moveItems)
    #two files with the same MOVE- PID
    with pytest.raises(SystemExit):
        internalSetup.findPidBodies({'WALL-a': ['MOVE-piston'], 'WALL-b': ['MOVE-piston']})
    #the PID is gone from the file: the section is dropped, a manual FROM section with a matching PID stays
    w['MOVE-crown'] = {'FROM': ['WALL-chamber-full-system'], 'ZONE_PIDS': ['hot']}
    kept = internalSetup.findFromSections(w, ['WALL-chamber-full-system'],
                                          {'WALL-chamber-full-system': ['WALL-chamber-hot', 'move-rod']})
    assert kept == ['MOVE-rod', 'MOVE-crown']


def test_match_pids_whole_name_first():
    pids = ['MOVE-piston', 'MOVE-piston-skirt', 'hot']
    assert writeInternal.matchPids(['MOVE-piston'], pids, 'x') == ['MOVE-piston']
    assert writeInternal.matchPids(['piston'], pids, 'x') == ['MOVE-piston', 'MOVE-piston-skirt']
    with pytest.raises(SystemExit):
        writeInternal.matchPids(['cold'], pids, 'x')


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
    for name in ['runScript', 'motionCheckScript', 'postScript', 'submitScript']:
        text = (case / name).read_bytes().decode()
        assert '\r' not in text and text.startswith('#!/bin/bash\n')
        assert re.search(r'<[A-Z_]+>', text) is None
    run = (case / 'runScript').read_text()
    assert '#SBATCH --account=def-account' in run and '#SBATCH --ntasks=192' in run
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


def test_post_script_and_chaining(tmp_path, monkeypatch, repo_root):
    case = _writeScripts(tmp_path, monkeypatch, repo_root, ['motion_a'])
    post = (case / 'postScript').read_text()
    assert '#SBATCH --job-name=engineRun01_post' in post and '#SBATCH --ntasks=1' in post
    assert '#SBATCH --cpus-per-task=4' in post and '#SBATCH --mem=32G' in post and '#SBATCH --time=0-03:00' in post
    assert 'module load StdEnv/2023 python scipy-stack paraview/5.13.1 ffmpeg' in post
    #summary, gsheet, images, videos in order, from the repo
    assert (post.index("/opt/caseSetup/postUtilities/postInternal.py'") < post.index("/opt/caseSetup/postUtilities/gsheetInternal.py'")
            < post.index("pvbatch --force-offscreen-rendering '/opt/caseSetup/postUtilities/pvPostInternal.py' > log.pvPostInternal")
            < post.index("/opt/caseSetup/postUtilities/imagesToMovies.py' > log.imagesToMovies"))
    #motion cases: postScript --mesh makes the mesh images and videos of the dry run, submitScript --check chains it
    assert "pvPostInternal.py' --mesh > log.pvPostMesh" in post and '--imageDir postProcessing/meshImages' in post
    assert post.index('"--mesh"') < post.index('log.postInternal')
    submit = (case / 'submitScript').read_text()
    assert submit.index('motionCheckScript)') < submit.index('postScript --mesh)') < submit.index('"--post"')
    motion = (case / 'motionCheckScript').read_text()
    assert "-name EnSightMotion -fields '(none)'" in motion and motion.index('EnSightMotion') < motion.index('motion check done')
    assert "checkMesh -parallel -time '0.0125,0.025,0.0375,0.05,0.0625,0.075,0.0875,0.1'" in motion
    submit = (case / 'submitScript').read_text()
    #runScript then postScript after it, --post alone, --noPost without
    assert 'jid=$(sbatch --parsable --dependency=afterok:$jid postScript)' in submit
    assert '"--post"' in submit and '"--noPost"' in submit
    assert submit.index('sbatch --parsable runScript)') < submit.index('afterok:$jid postScript)')
    #no gsheet push, no EnSight export: no gsheet / image / video steps
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    full['GLOBAL_COMPUTE_SETUP']['GSHEET'] = ['False']
    full['GLOBAL_COMPUTE_SETUP']['EXPORT'] = ['none']
    full['GLOBAL_COMPUTE_SETUP']['PVBATCH'] = ['/opt/pv/bin/pvbatch']
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', None)
    post = (case / 'postScript').read_text()
    assert "postUtilities/postInternal.py'" in post and "postUtilities/gsheetInternal.py'" not in post
    assert '"--mesh"' not in post and 'postScript --mesh' not in (case / 'submitScript').read_text()
    assert '--force-offscreen-rendering' not in post and "postUtilities/imagesToMovies.py' > log.imagesToMovies" not in post
    full['GLOBAL_COMPUTE_SETUP']['GSHEET'] = ['maybe']
    with pytest.raises(SystemExit):
        writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', None)


def test_post_setup_copied_once(tmp_path, repo_root):
    import internalSetup
    internalSetup.copyPostSetup(repo_root, str(tmp_path))
    dest = tmp_path / 'pvPostInternalSetup'
    assert '[SLICE]' in dest.read_text()
    dest.write_text('[SLICE]\nNAMES = mine\n')
    internalSetup.copyPostSetup(repo_root, str(tmp_path))
    assert dest.read_text() == '[SLICE]\nNAMES = mine\n'


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
    full['GLOBAL_SIM_CONTROL']['SOLVER'] = ['interFoam']
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


def test_relaxation_only_with_outer_correctors(tmp_path, monkeypatch, repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    monkeypatch.chdir(tmp_path)
    writeInternal.writeNumericsInternal(templateLoc, full)
    text = (tmp_path / 'system/fvSolution').read_text()
    assert 'nOuterCorrectors 1;' in text and 'relaxationFactors' not in text and '<RELAXATION>' not in text
    full['GLOBAL_SIM_CONTROL']['N_OUTER_CORR'] = ['3']
    writeInternal.writeNumericsInternal(templateLoc, full)
    text = (tmp_path / 'system/fvSolution').read_text()
    assert 'nOuterCorrectors 3;' in text and 'relaxationFactors' in text
    assert 'p               0.3;' in text and 'pFinal          1;' in text
    assert '"(U|k|omega|epsilon)" 0.7;' in text and 'e               0.7;' in text and '".*Final"       1;' in text
    #incompressible: no energy equation
    full['GLOBAL_SIM_CONTROL']['SOLVER'] = ['pimpleFoam']
    full['GLOBAL_SIM_CONTROL']['RELAX_U'] = ['0.5']
    writeInternal.writeNumericsInternal(templateLoc, full)
    text = (tmp_path / 'system/fvSolution').read_text()
    assert '"(U|k|omega|epsilon)" 0.5;' in text and ' e  ' not in text
    full['GLOBAL_SIM_CONTROL']['RELAX_P'] = ['1.5']
    with pytest.raises(SystemExit):
        writeInternal.writeNumericsInternal(templateLoc, full)


def test_fixed_time_step_when_max_co_blank_or_zero(tmp_path, monkeypatch, repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    monkeypatch.chdir(tmp_path)
    full['GLOBAL_CONTROL']['MAX_DT'] = ['2e-5']
    writeInternal.writeControlDictInternal(templateLoc, full)
    text = (tmp_path / 'system/controlDict').read_text()
    assert 'adjustTimeStep  yes;' in text and 'maxCo           1;' in text
    for maxCo in (['0'], [''], ['0.0']):
        full['GLOBAL_CONTROL']['MAX_CO'] = maxCo
        writeInternal.writeControlDictInternal(templateLoc, full)
        text = (tmp_path / 'system/controlDict').read_text()
        assert 'adjustTimeStep  no;' in text and 'deltaT          2e-5;' in text
    full['GLOBAL_CONTROL']['MAX_CO'] = ['x']
    with pytest.raises(SystemExit):
        writeInternal.writeControlDictInternal(templateLoc, full)
    full['GLOBAL_CONTROL']['MAX_CO'] = ['']
    full['GLOBAL_CONTROL']['MAX_DT'] = ['']
    with pytest.raises(SystemExit):
        writeInternal.writeControlDictInternal(templateLoc, full)


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


#----- decomposition without scotch

def _decompTool(repo_root):
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'decompCoeffs', os.path.join(repo_root, 'meshUtilities', 'decompCoeffs.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_decomp_split_and_dict(tmp_path, monkeypatch, repo_root):
    tool = _decompTool(repo_root)
    #flat box (Stirling chamber): no splits across the thin direction while the others can take them
    assert tool.chooseSplit(4, [0.4, 0.4, 0.07]) == (2, 2, 1)
    n = tool.chooseSplit(192, [0.4, 0.4, 0.07])
    assert n[0]*n[1]*n[2] == 192 and n[2] < n[0] and n[2] < n[1]
    assert tool.chooseSplit(7, [1.0, 1.0, 1.0]) in [(7, 1, 1), (1, 7, 1), (1, 1, 7)]
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'system').mkdir()
    (tmp_path / 'system' / 'blockMeshDict').write_text(
        'vertices\n(\n    (-0.2 -0.2 0) (0.2 -0.2 0) (0.2 0.2 0) (-0.2 0.2 0)\n'
        '    (-0.2 -0.2 0.07) (0.2 -0.2 0.07) (0.2 0.2 0.07) (-0.2 0.2 0.07)\n);\n')
    path = tmp_path / 'system' / 'decomposeParDict'
    path.write_text('FoamFile {}\nnumberOfSubdomains 192;\nmethod          hierarchical;\n\n// ***** //\n')
    tool.main([str(path), '4'])
    #running again (other task count) replaces the entries, no duplicates
    tool.main([str(path), '6'])
    text = path.read_text()
    assert 'numberOfSubdomains 6;' in text and text.count('coeffs') == 1
    assert re.search(r'n\s+\((\d+) (\d+) 1\);', text) and text.rstrip().endswith('// ***** //')


def test_scripts_hierarchical_decomposition(tmp_path, monkeypatch, repo_root):
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    case = tmp_path / 'run'
    case.mkdir()
    monkeypatch.chdir(case)
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    run = (case / 'runScript').read_text()
    assert run.count('-entry numberOfSubdomains -set $NP') == 2 and 'decompCoeffs' not in run
    full['GLOBAL_COMPUTE_SETUP']['DECOMP_METHOD'] = ['hierarchical']
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    run = (case / 'runScript').read_text()
    motion = (case / 'motionCheckScript').read_text()
    #before the meshing (decomposePar) and before the solve, and in the motion check
    call = "python3 '/opt/caseSetup/meshUtilities/decompCoeffs.py' system/decomposeParDict $NP"
    assert run.count(call) == 2 and run.index(call) < run.index('decomposePar -force')
    assert call in motion and 'numberOfSubdomains -set' not in run + motion


def test_decomp_override(tmp_path, monkeypatch, repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    compute = full['GLOBAL_COMPUTE_SETUP']
    assert writeInternal.getDecomp(full) == ('scotch', None)
    compute['DECOMP_METHOD'] = ['hierarchical']
    assert writeInternal.getDecomp(full) == ('hierarchical', None)
    compute['DECOMP'] = ['8', '6', '4']
    assert writeInternal.getDecomp(full) == ('hierarchical', (8, 6, 4))
    #dict as caseSetup: numberOfSubdomains = nx*ny*nz, n in the coeffs
    monkeypatch.chdir(tmp_path)
    writeInternal.writeNumericsInternal(templateLoc, full)
    text = (tmp_path / 'system/decomposeParDict').read_text()
    assert 'numberOfSubdomains 192;' in text and 'method          hierarchical;' in text
    assert re.search(r'coeffs\s*\{\s*n\s+\(8 6 4\);\s*\}', text)
    #scripts pass the split on, the helper checks it against the task count
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    assert "decompCoeffs.py' system/decomposeParDict $NP 8 6 4;" in (tmp_path / 'runScript').read_text()
    for bad in (['8', '6'], ['8', '6', 'x'], ['0', '6', '4'], ['8', '6', '2']):
        compute['DECOMP'] = bad
        with pytest.raises(SystemExit):
            writeInternal.getDecomp(full)
    #other methods ignore DECOMP
    compute['DECOMP'], compute['DECOMP_METHOD'] = ['8', '6', '2'], ['scotch']
    assert writeInternal.getDecomp(full) == ('scotch', None)


def test_decomp_helper_given_split(tmp_path, monkeypatch, repo_root):
    tool = _decompTool(repo_root)
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'system').mkdir()
    path = tmp_path / 'system' / 'decomposeParDict'
    path.write_text('numberOfSubdomains 1;\nmethod          hierarchical;\ncoeffs\n{\n    n (1 1 1);\n}\n\n// ***** //\n')
    #no blockMeshDict needed when the split is given
    tool.main([str(path), '4', '2', '2', '1'])
    text = path.read_text()
    assert 'numberOfSubdomains 4;' in text and '(2 2 1);' in text and text.count('coeffs') == 1
    with pytest.raises(SystemExit):
        tool.main([str(path), '6', '2', '2', '1'])


#----- through-flow (inlets, outlets, steady solvers)

def test_flow_group_lines():
    group = writeInternal.parseGroupLine('PIDS inlet U massFlow 0.5 T fixed 600'.split(), 'g')
    assert group['flow'] == {'kind': 'inlet', 'type': 'massFlow', 'value': 0.5} and group['T'] == ['fixed', 600.0]
    assert writeInternal.parseGroupLine('PIDS in U velocity 5'.split(), 'g')['flow']['type'] == 'velocity'
    vector = writeInternal.parseGroupLine('PIDS in U velocity 0 0 3'.split(), 'g')['flow']
    assert vector == {'kind': 'inlet', 'type': 'velocityVector', 'value': [0.0, 0.0, 3.0]}
    assert writeInternal.parseGroupLine('PIDS in P total 1.01e5'.split(), 'g')['flow']['type'] == 'totalPressure'
    assert writeInternal.parseGroupLine('PIDS out P fixed 1e5'.split(), 'g')['flow']['kind'] == 'outlet'
    assert writeInternal.parseGroupLine('PIDS wall U slip'.split(), 'g')['flow'] is None
    for bad in ['PIDS in U massFlow 0.5 P fixed 1e5', 'PIDS in U massFlow -1', 'PIDS in U massFlow 0.5 T convective 10 300',
                'PIDS out P total', 'PIDS in U speed 5', 'PIDS w U slip P fixed 1e5']:
        with pytest.raises(SystemExit):
            writeInternal.parseGroupLine(bad.split(), bad)
    #inlet/outlet T defaults to None (INIT_T), walls keep their default
    groups = writeInternal.resolveBCGroups('WALL-pipe', {'U_DEFAULT': ['noSlip'], 'T_DEFAULT': ['adiabatic'],
                                                         'BC_GROUPS': ['pIn', 'pOut'],
                                                         'pIn': 'PIDS inlet U massFlow 0.02'.split(),
                                                         'pOut': 'PIDS outlet P fixed 1e5'.split()},
                                           ['inlet', 'outlet', 'wall'])
    assert [g['name'] for g in groups] == ['pIn', 'pOut', 'WALL-pipe_default']
    assert groups[0]['T'] is None and groups[2]['T'] == ['adiabatic'] and groups[2]['flow'] is None


def _pipeCase(tmp_path, monkeypatch):
    #constant/triSurface/WALL-pipe.stl with a 0.1 x 0.1 inlet square (area 0.01), and its wall groups
    monkeypatch.chdir(tmp_path)
    surf = tmp_path / 'constant' / 'triSurface'
    surf.mkdir(parents=True)
    (tmp_path / 'system').mkdir()
    square = lambda z: ['facet normal 0 0 1\n outer loop\n  vertex 0 0 %g\n  vertex 0.1 0 %g\n  vertex 0.1 0.1 %g\n endloop\nendfacet\n' % (z, z, z),
                        'facet normal 0 0 1\n outer loop\n  vertex 0 0 %g\n  vertex 0.1 0.1 %g\n  vertex 0 0.1 %g\n endloop\nendfacet\n' % (z, z, z)]
    text = 'solid inlet\n%sendsolid inlet\nsolid outlet\n%sendsolid outlet\nsolid wall\n%sendsolid wall\n' \
           % (''.join(square(0)), ''.join(square(1)), ''.join(square(0.5)))
    (surf / 'WALL-pipe.stl').write_text(text)
    geomDict = {'WALL-pipe.stl': {'scale': '1', 'refinement': '1', 'layers': '0', 'expansion': '1.2', 'wallmodel': 'high'}}
    groups = writeInternal.resolveBCGroups('WALL-pipe', {'U_DEFAULT': ['noSlip'], 'T_DEFAULT': ['flux', '-500'],
                                                         'BC_GROUPS': ['pIn', 'pOut'],
                                                         'pIn': 'PIDS inlet U massFlow 0.02 T fixed 600'.split(),
                                                         'pOut': 'PIDS outlet P fixed 1e5'.split()},
                                           ['inlet', 'outlet', 'wall'])
    for group in groups:
        group['moving'] = []
    return geomDict, {'WALL-pipe': groups}


def test_flow_entries_compressible_and_incompressible(tmp_path, monkeypatch, repo_root):
    geomDict, wallGroups = _pipeCase(tmp_path, monkeypatch)
    templateLoc, full = _fullDefaults(repo_root, 'rhoSimpleFoam')
    full['GLOBAL_SIM_CONTROL']['TURB_MODEL'] = ['kOmegaSST']
    writeInternal.writeCasePropertiesInternal(templateLoc, geomDict, full, wallGroups)
    text = (tmp_path / 'system/caseProperties').read_text()
    inlet = [l for l in text.splitlines() if l.strip().startswith('pIn ')][0]
    assert 'category inlet; type subSonic;' in inlet and 'flowSpecification flowRate;' in inlet
    assert 'massFlowRate constant 0.02;' in inlet and 'T uniform 600;' in inlet
    #k from 5 % intensity of the mean inlet velocity mdot/(rho A), rho of air at 1 bar 600 K
    rho = 1e5*28.96/(8314.46*600)
    speed = 0.02/(rho*0.01)
    k = float(re.search(r'k uniform ([0-9.e+-]+);', inlet).group(1))
    assert k == pytest.approx(1.5*(0.05*speed)**2, rel=1e-4)
    outlet = [l for l in text.splitlines() if l.strip().startswith('pOut ')][0]
    assert 'category outlet;' in outlet and 'returnFlow staticPressure;' in outlet and 'p uniform 100000;' in outlet
    assert 'T uniform 300;' in outlet
    wall = [l for l in text.splitlines() if 'WALL-pipe_default' in l][0]
    assert 'category wall;' in wall and 'heatTransfer heatFlux;' in wall
    #incompressible: volume flow from the charge density, kinematic gauge p
    templateLoc, full = _fullDefaults(repo_root, 'simpleFoam')
    writeInternal.writeCasePropertiesInternal(templateLoc, geomDict, full, wallGroups)
    text = (tmp_path / 'system/caseProperties').read_text()
    rhoCharge = writeInternal.getGasDensity(full)
    inlet = [l for l in text.splitlines() if l.strip().startswith('pIn ')][0]
    flow = float(re.search(r'volumetricFlowRate constant ([0-9.e+-]+);', inlet).group(1))
    assert 'flowSpecification volumeFlowRate;' in inlet and flow == pytest.approx(0.02/rhoCharge, rel=1e-6)
    assert 'T uniform' not in inlet and 'p uniform 0;' in [l for l in text.splitlines() if l.strip().startswith('pOut ')][0]
    #buoyant: p_rgh starts at the charge pressure
    templateLoc, full = _fullDefaults(repo_root, 'buoyantSimpleFoam')
    writeInternal.writeCasePropertiesInternal(templateLoc, geomDict, full, wallGroups)
    assert 'p_rgh uniform 100000;' in (tmp_path / 'system/caseProperties').read_text()


def test_steady_numerics_and_control(tmp_path, monkeypatch, repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'rhoSimpleFoam')
    monkeypatch.chdir(tmp_path)
    writeInternal.writeNumericsInternal(templateLoc, full, pressureBC=True)
    schemes = (tmp_path / 'system/fvSchemes').read_text()
    assert 'default steadyState;' in schemes and 'div(phi,U)      bounded Gauss linearUpwind grad(U);' in schemes
    assert 'div(phi,Ekp)' in schemes
    solution = (tmp_path / 'system/fvSolution').read_text()
    assert 'SIMPLE' in solution and 'PIMPLE' not in solution and 'residualControl' in solution
    assert 'pMinFactor' in solution and 'pRefPoint' not in solution and '<' not in solution
    assert 'p               0.3;' in solution and 'rho             1;' in solution
    #closed incompressible steady volume keeps the pressure reference, buoyant relaxes p_rgh
    templateLoc, full = _fullDefaults(repo_root, 'simpleFoam')
    writeInternal.writeNumericsInternal(templateLoc, full)
    solution = (tmp_path / 'system/fvSolution').read_text()
    assert 'pRefPoint ( 0.001 0.002 0.003 );' in solution and 'pMinFactor' not in solution
    templateLoc, full = _fullDefaults(repo_root, 'buoyantSimpleFoam')
    writeInternal.writeNumericsInternal(templateLoc, full, pressureBC=True)
    assert 'p_rgh           0.3;' in (tmp_path / 'system/fvSolution').read_text()
    #transient runs keep PIMPLE and Euler
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    writeInternal.writeNumericsInternal(templateLoc, full)
    assert 'PIMPLE' in (tmp_path / 'system/fvSolution').read_text()
    assert 'default Euler;' in (tmp_path / 'system/fvSchemes').read_text()
    #steady controlDict: N_ITER iterations of deltaT 1
    templateLoc, full = _fullDefaults(repo_root, 'rhoSimpleFoam')
    full['GLOBAL_CONTROL']['N_ITER'] = ['400']
    full['GLOBAL_CONTROL']['WRITE_ITER'] = ['100']
    writeInternal.writeControlDictInternal(templateLoc, full)
    control = (tmp_path / 'system/controlDict').read_text()
    assert 'endTime         400;' in control and 'deltaT          1;' in control
    assert 'writeInterval   100;' in control and 'adjustTimeStep  no;' in control


def test_buoyant_thermo_and_gravity(tmp_path, monkeypatch, repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'buoyantSimpleFoam')
    monkeypatch.chdir(tmp_path)
    full['GLOBAL_SIM_CONTROL']['GRAVITY'] = ['0', '-9.81', '0']
    writeInternal.writeThermoInternal(templateLoc, full)
    assert 'heRhoThermo' in (tmp_path / 'constant/thermophysicalProperties').read_text()
    assert 'value           (0 -9.81 0);' in (tmp_path / 'constant/g').read_text()
    templateLoc, full = _fullDefaults(repo_root, 'rhoSimpleFoam')
    writeInternal.writeThermoInternal(templateLoc, full)
    assert 'hePsiThermo' in (tmp_path / 'constant/thermophysicalProperties').read_text()
    assert not (tmp_path / 'constant/g').exists()


def test_flow_function_objects(tmp_path, monkeypatch, repo_root):
    geomDict, wallGroups = _pipeCase(tmp_path, monkeypatch)
    templateLoc, full = _fullDefaults(repo_root, 'rhoSimpleFoam')
    writeInternal.writeFunctionObjectsInternal(wallGroups, full)
    fos = (tmp_path / 'system/internalFunctionObjects').read_text()
    assert 'flow_WALL-pipe_inlet' in fos and 'avg_WALL-pipe_outlet' in fos and 'residuals' in fos
    #wall heat flux on the walls only
    heat = fos[fos.index('wallHeatFlux\n{'):]
    heat = heat[:heat.index('}')]
    assert 'WALL-pipe_wall' in heat and 'WALL-pipe_inlet' not in heat
    #bookkeeping for post: group kinds and the reference density
    info = writeInternal.addSetupInfo(full, wallGroups)
    assert info['BC_GROUP_KINDS']['pIn'] == ['inlet', 'massFlow', '0.02']
    assert info['BC_GROUP_KINDS']['pOut'][0] == 'outlet' and info['SETUP_INFO']['STEADY'] == ['yes']
    assert writeInternal.hasPressureBC(wallGroups)


#----- fans and porous zones

def _porDefaults(repo_root, solver, **keys):
    templateLoc, full = _fullDefaults(repo_root, solver)
    items = dict(setupCore.readAddonTemplate(templateLoc, 'defaultPorous', 'POR_SETUP'))
    section = dict((key, str(value).split()) for key, value in items.items())
    for key, value in keys.items():
        section[key] = value.split()
    full['POR-matrix'] = section
    return templateLoc, full, {'POR-matrix.stl': {}}


def test_porous_resistance_inputs(repo_root):
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoSimpleFoam', D='1e6 2e6 3e6', F='10 20 30')
    zone = writeInternal.getPorous(geomDict, full)[0]
    assert zone['d'] == [1e6, 2e6, 3e6] and zone['f'] == [10, 20, 30] and zone['heat'] == 'none'
    #wire mesh: Gedeon-Wood, Darcy part 129/(2 eps dh^2)
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoSimpleFoam', RESISTANCE='wireMesh', POROSITY='0.8',
                                               WIRE_DIA='1e-4', REF_VELOCITY='2', HEAT='fixed', SOLID_T='550')
    zone = writeInternal.getPorous(geomDict, full)[0]
    dh = 0.8*1e-4/0.2
    assert zone['d'][0] == pytest.approx(129.0/(2*0.8*dh**2)) and zone['f'][0] > 0
    assert zone['heat'] == 'fixed' and zone['ha'] > 0 and zone['Ts'] == 550 and zone['Cs'] == 0
    #fit: data made from d 5e6, f 400 is fitted back
    gas = {'RHO': writeInternal.getGasDensity(full), 'MU': writeInternal.getGasProperties(full)['MU']}
    pairs = []
    for u in (0.5, 1.0, 2.0, 4.0):
        pairs += [u, 0.05*(gas['MU']*5e6*u + 0.5*gas['RHO']*400*u*u)]
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoSimpleFoam', RESISTANCE='fit', DP_LENGTH='0.05',
                                               DP_DATA=' '.join('%.12g' % v for v in pairs))
    zone = writeInternal.getPorous(geomDict, full)[0]
    assert zone['d'][0] == pytest.approx(5e6, rel=1e-6) and zone['f'][0] == pytest.approx(400, rel=1e-6)
    #errors: HA auto without wire mesh, solved matrix in a steady run, bad RESISTANCE
    for solver, keys in [('rhoSimpleFoam', {'HEAT': 'fixed'}), ('rhoSimpleFoam', {'HEAT': 'solved', 'HA': '1e4'}),
                         ('rhoSimpleFoam', {'RESISTANCE': 'magic'}), ('rhoSimpleFoam', {'D': '1 2'})]:
        templateLoc, full, geomDict = _porDefaults(repo_root, solver, **keys)
        with pytest.raises(SystemExit):
            writeInternal.getPorous(geomDict, full)
    #incompressible: heat ignored
    templateLoc, full, geomDict = _porDefaults(repo_root, 'simpleFoam', HEAT='fixed', HA='1e4')
    assert writeInternal.getPorous(geomDict, full)[0]['heat'] == 'none'


def test_porous_fvoptions_and_function_objects(tmp_path, monkeypatch, repo_root):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'constant').mkdir()
    (tmp_path / 'system').mkdir()
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoPimpleFoam', HEAT='solved', HA='2e4', SOLID_T='600')
    porous = writeInternal.getPorous(geomDict, full)
    assert porous[0]['Cs'] == pytest.approx(0.1*7900*500)
    writeInternal.writeFvOptionsInternal(porous, full)
    text = (tmp_path / 'constant/fvOptions').read_text()
    assert 'type            explicitPorositySource;' in text and 'cellZone        POR-matrix;' in text
    assert 'd   (10000000 10000000 10000000);' in text and 'rotation { type axes; e1 (1 0 0); e2 (0 1 0); }' in text
    assert 'type            scalarCodedSource;' in text and 'fields          (e);' in text
    assert 'const scalar ha = 20000;' in text and 'const scalar Cs = 395000;' in text
    assert 'Ts_POR-matrix' in text and text.count('#{') == text.count('#}')
    #function objects: zone averages incl. the matrix T, heat into the gas
    writeInternal.writeFunctionObjectsInternal({}, full, None, porous)
    fos = (tmp_path / 'system/internalFunctionObjects').read_text()
    assert 'porHeat_POR-matrix' in fos and '(Qs_POR-matrix);' in fos and '(p T Ts_POR-matrix);' in fos
    #no porous zone: fvOptions removed
    writeInternal.writeFvOptionsInternal([], full)
    assert not (tmp_path / 'constant/fvOptions').exists()


def test_fan_baffles(tmp_path, monkeypatch, repo_root):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'system').mkdir()
    templateLoc, full = _fullDefaults(repo_root, 'simpleFoam')
    full['FAN-disk'] = {'FAN_PRESSURE_RISE': ['50']}
    fans = writeInternal.writeFanBaffles({'WALL-pipe.stl': {}, 'FAN-disk.stl': {}}, full)
    text = (tmp_path / 'system/createBafflesDict').read_text()
    jump = 50.0/writeInternal.getGasDensity(full)
    assert fans == ['FAN-disk'] and 'zoneName FAN-disk;' in text and 'FAN-disk_master' in text
    assert 'jumpTable constant %.9g;' % (jump) in text and 'p { type fan; patchType cyclic;' in text
    #buoyant: jump on p_rgh in Pa
    templateLoc, full = _fullDefaults(repo_root, 'buoyantSimpleFoam')
    full['FAN-disk'] = {'FAN_PRESSURE_RISE': ['50']}
    writeInternal.writeFanBaffles({'FAN-disk.stl': {}}, full)
    assert 'p_rgh { type fan; patchType cyclic; jump uniform 0; value uniform 100000; jumpTable constant 50; }' \
        in (tmp_path / 'system/createBafflesDict').read_text()
    assert writeInternal.writeFanBaffles({'WALL-pipe.stl': {}}, full) == []
    assert not (tmp_path / 'system/createBafflesDict').exists()


def test_motion_frames(repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    full['MESH_MOTION']['MOTION_FRAMES'] = ['72']
    assert writeInternal.getMotionFrames(full) == (72, 72)
    full['MESH_MOTION']['MOTION_FRAMES'] = ['24']
    assert writeInternal.getMotionFrames(full) == (24, 72)
    full['MESH_MOTION']['MOTION_FRAMES'] = ['144']
    assert writeInternal.getMotionFrames(full) == (144, 144)
    for bad in ['30', '0', 'many']:
        full['MESH_MOTION']['MOTION_FRAMES'] = [bad]
        with pytest.raises(SystemExit):
            writeInternal.getMotionFrames(full)


def test_mpi_launcher_key(tmp_path, monkeypatch, repo_root):
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    monkeypatch.chdir(tmp_path)
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    for name in ['runScript', 'motionCheckScript']:
        assert 'RUN="srun"' in (tmp_path / name).read_text()
    #a cluster whose Open MPI cannot start under srun: every script uses the given launcher
    full['GLOBAL_COMPUTE_SETUP']['MPI_RUN'] = '/usr/bin/mpirun /opt/of/bin/tools/foamExec'.split()
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    for name in ['runScript', 'motionCheckScript']:
        text = (tmp_path / name).read_text()
        assert 'RUN="/usr/bin/mpirun /opt/of/bin/tools/foamExec"' in text and 'srun' not in text
    #missing key (older template sets): srun
    del full['GLOBAL_COMPUTE_SETUP']['MPI_RUN']
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    assert '<MPI_RUN>' not in (tmp_path / 'motionCheckScript').read_text()


def test_mesh_signature_and_stale_check(tmp_path, monkeypatch, repo_root):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'system').mkdir()
    surf = tmp_path / 'constant' / 'triSurface'
    surf.mkdir(parents=True)
    (tmp_path / 'system/snappyHexMeshDict').write_text('a')
    (tmp_path / 'system/topoSetDict').write_text('zones 4')
    (surf / 'WALL-x.stl').write_text('solid a\nendsolid a\n')
    (surf / 'WALL-x.eMesh').write_text('made by surfaceFeatureExtract')
    writeInternal.writeMeshSignature()
    first = (tmp_path / 'system/meshSignature').read_text()
    assert 'system/topoSetDict' in first and 'constant/triSurface/WALL-x.stl' in first and 'eMesh' not in first
    #generated feature edges do not change it, a new motion zone does
    (surf / 'WALL-x.eMesh').write_text('changed')
    writeInternal.writeMeshSignature()
    assert (tmp_path / 'system/meshSignature').read_text() == first
    (tmp_path / 'system/topoSetDict').write_text('zones 5')
    writeInternal.writeMeshSignature()
    assert (tmp_path / 'system/meshSignature').read_text() != first
    #the scripts compare it with the copy kept in .meshDone
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    run = (tmp_path / 'runScript').read_text()
    assert '! cmp -s system/meshSignature .meshDone' in run and 'cp system/meshSignature .meshDone;' in run
    assert run.index('cmp -s system/meshSignature') < run.index('blockMesh')
    assert '! cmp -s system/meshSignature .meshDone' in (tmp_path / 'motionCheckScript').read_text()


def test_copy_back_step(tmp_path, monkeypatch, repo_root):
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    monkeypatch.chdir(tmp_path)
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', None)
    post = (tmp_path / 'postScript').read_text()
    assert post.index("imagesToMovies.py' > log.imagesToMovies") < post.index("copyBackInternal.py' > log.copyBackInternal")
    full['GLOBAL_COMPUTE_SETUP']['COPY_BACK'] = ['False']
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', None)
    assert "copyBackInternal.py' > log.copyBackInternal" not in (tmp_path / 'postScript').read_text()
    full['GLOBAL_COMPUTE_SETUP']['COPY_BACK'] = ['sometimes']
    with pytest.raises(SystemExit):
        writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', None)
