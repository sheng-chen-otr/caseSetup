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
import subprocess
import types

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
    assert list(piston) == [key for key, val in moveItems if key not in internalSetup.FROM_UNUSED]
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


def test_scripts_without_account_or_modules(tmp_path, monkeypatch, repo_root):
    #in-house cluster: partition and node fixed, no account, no environment modules, mpirun through foamExec
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    #cluster settings in their own module, the node in the per-run compute settings
    env = full['CLUSTER_ENVIRONMENT']
    env['ACCOUNT'], env['MODULES'], env['POST_MODULES'] = [''], [''], ['']
    env['PARTITION'], env['MPI_RUN'], env['POST_MEM'] = ['solve2'], ['/usr/bin/mpirun', '/of/bin/tools/foamExec'], ['']
    full['GLOBAL_COMPUTE_SETUP']['NODELIST'] = ['cfd07']
    assert 'ACCOUNT' not in full['GLOBAL_COMPUTE_SETUP'] and 'NODELIST' not in env
    case = tmp_path / '008'
    case.mkdir()
    monkeypatch.chdir(case)
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    for name in ['runScript', 'motionCheckScript', 'postScript']:
        text = (case / name).read_text()
        assert '--account' not in text and 'module' not in text
        assert '#SBATCH --job-name=008' in text and '#SBATCH --partition=solve2\n#SBATCH --nodelist=cfd07\n' in text
    assert 'RUN="/usr/bin/mpirun /of/bin/tools/foamExec"' in (case / 'runScript').read_text()
    post = (case / 'postScript').read_text()
    assert '--mem' not in post and '#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=8\n' in post
    #no if left empty by the removed module line (bash syntax error, the post job would stop at once)
    for name in ['runScript', 'motionCheckScript', 'postScript', 'submitScript']:
        if (case / name).exists():
            assert not re.search(r'then\s*(fi|else)', (case / name).read_text())
            assert subprocess.run(['bash', '-n', str(case / name)]).returncode == 0


def test_scripts_shared_node(tmp_path, monkeypatch, repo_root):
    #EXCLUSIVE False: no --exclusive, --mem=0 replaced by MEM_PER_CPU (blank: no memory line), post job unchanged
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    assert full['GLOBAL_COMPUTE_SETUP']['EXCLUSIVE'] == ['True']
    case = tmp_path / '002'
    case.mkdir()
    monkeypatch.chdir(case)
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    assert '#SBATCH --mem=0\n' in (case / 'runScript').read_text() and '#SBATCH --exclusive' in (case / 'runScript').read_text()
    full['GLOBAL_COMPUTE_SETUP']['EXCLUSIVE'], full['GLOBAL_COMPUTE_SETUP']['MEM_PER_CPU'] = ['False'], ['4G']
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    for name in ['runScript', 'motionCheckScript']:
        text = (case / name).read_text()
        assert '--exclusive' not in text and '--mem=0' not in text and '#SBATCH --mem-per-cpu=4G\n' in text
        assert subprocess.run(['bash', '-n', str(case / name)]).returncode == 0
    assert '#SBATCH --mem=32G' in (case / 'postScript').read_text()
    full['GLOBAL_COMPUTE_SETUP']['MEM_PER_CPU'] = ['']
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    assert '--mem' not in (case / 'runScript').read_text()
    full['GLOBAL_COMPUTE_SETUP']['EXCLUSIVE'] = ['maybe']
    with pytest.raises(SystemExit):
        writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])


def test_scripts_cpus_per_task(tmp_path, monkeypatch, repo_root):
    #CPUS_PER_TASK 4: ranks spread over the node, srun told explicitly, post job unchanged
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    assert full['GLOBAL_COMPUTE_SETUP']['CPUS_PER_TASK'] == ['1']
    case = tmp_path / 'engine'
    case.mkdir()
    monkeypatch.chdir(case)
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    assert 'SRUN_CPUS_PER_TASK' not in (case / 'runScript').read_text()
    full['GLOBAL_COMPUTE_SETUP']['NUM_PROCS'], full['GLOBAL_COMPUTE_SETUP']['CPUS_PER_TASK'] = ['48'], ['4']
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    for name in ['runScript', 'motionCheckScript']:
        text = (case / name).read_text()
        assert '#SBATCH --ntasks=48\n#SBATCH --cpus-per-task=4\n' in text
        assert 'export SRUN_CPUS_PER_TASK=4 OMP_NUM_THREADS=1;' in text
        assert subprocess.run(['bash', '-n', str(case / name)]).returncode == 0
    assert '#SBATCH --cpus-per-task=8' in (case / 'postScript').read_text()
    full['GLOBAL_COMPUTE_SETUP']['CPUS_PER_TASK'] = ['0']
    with pytest.raises(SystemExit):
        writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])


def test_cluster_keys_in_old_compute_section(tmp_path, monkeypatch, repo_root, capsys):
    #case files from before CLUSTER_ENVIRONMENT: cluster keys left in GLOBAL_COMPUTE_SETUP are still used
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    full['GLOBAL_COMPUTE_SETUP']['MODULES'] = ['openfoam/old']
    full['GLOBAL_COMPUTE_SETUP']['GSHEET'] = ['False']
    case = tmp_path / 'old'
    case.mkdir()
    monkeypatch.chdir(case)
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    run, post = (case / 'runScript').read_text(), (case / 'postScript').read_text()
    assert 'module load openfoam/old' in run and 'gsheetInternal.py' not in post
    assert 'belong in CLUSTER_ENVIRONMENT now' in capsys.readouterr().out


def test_post_script_and_chaining(tmp_path, monkeypatch, repo_root):
    case = _writeScripts(tmp_path, monkeypatch, repo_root, ['motion_a'])
    post = (case / 'postScript').read_text()
    assert '#SBATCH --job-name=engineRun01_post' in post and '#SBATCH --ntasks=1' in post
    assert '#SBATCH --cpus-per-task=8' in post and '#SBATCH --mem=32G' in post and '#SBATCH --time=0-03:00' in post
    assert 'module load StdEnv/2023 python/3.11 scipy-stack paraview/5.13.1 ffmpeg' in post
    #summary, gsheet, images, videos in order, from the repo
    assert (post.index("/opt/caseSetup/postUtilities/postInternal.py'") < post.index("/opt/caseSetup/postUtilities/gsheetInternal.py'")
            < post.index("pvbatch --force-offscreen-rendering '/opt/caseSetup/postUtilities/pvPostInternal.py' > log.pvPostInternal")
            < post.index("/opt/caseSetup/postUtilities/imagesToMovies.py' > log.imagesToMovies"))
    #motion cases: postScript --mesh makes the mesh images and videos of the dry run, submitScript --check chains it
    assert "pvPostInternal.py' --mesh > log.pvPostMesh" in post and '--imageDir postProcessing/meshImages' in post
    assert post.index('"--mesh"') < post.index('log.postInternal')
    #no post environment set: nothing activated, settings from the default folder
    assert post.count('(if [ -n "" ]; then source /bin/activate; fi; export CASESETUP_CONFIG=;') == 2
    assert subprocess.run(['bash', '-n', str(case / 'postScript')]).returncode == 0
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
    full['GLOBAL_REFINEMENT']['LOC_IN_MESH'] = ['0.0013', '0.0021', '0.0033']
    return templateLoc, full


def test_segment_snappy_dicts(tmp_path, monkeypatch, repo_root):
    #remeshing: the case dict keeps LOC_IN_MESH as given, each other mesh gets its own dict with its location and
    #z-only end-gap bands; the run script copies it over the case's
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    full['WALL-chamber'] = {}
    geomDict = {'WALL-chamber.stl': {'refinement': '1', 'scale': '1', 'layers': '0', 'expansion': '1.2'}}
    motionInfo = {'gapRefs': [], 'segments': [(0.0, 90.0), (90.0, 360.0)],
                  'meshSetups': [{'loc': [0.0013, 0.0021, 0.0033], 'bands': []},
                                 {'loc': [0.0013, 0.0021, 0.0287936], 'bands': [(0.026, 0.0315, 2)]}]}
    monkeypatch.chdir(tmp_path)
    writeInternal.writeSnappyInternal(templateLoc, geomDict, full, motionInfo)
    case = (tmp_path / 'system/snappyHexMeshDict').read_text()
    seg = (tmp_path / 'meshes/seg02/snappyHexMeshDict').read_text()
    assert 'locationInMesh (0.0013 0.0021 0.0033);' in case and 'GAPZ' not in case
    assert 'locationInMesh (0.0013 0.0021 0.0287936);' in seg
    assert 'GAPZ-0 { type searchableBox; min (-1e3 -1e3 0.026); max (1e3 1e3 0.0315); }' in seg
    assert 'GAPZ-0 { mode inside; levels ((1e15 0)); levelIncrement (0 20 (0 0 2)); }' in seg
    assert not (tmp_path / 'meshes/seg01').exists()
    #a moved segment point on a cell face is shifted off it (z faces every 4 mm / 2^(1+4) here), the given point is an error
    motionInfo['meshSetups'][1]['loc'] = [0.0013, 0.0021, 0.03]
    writeInternal.writeSnappyInternal(templateLoc, geomDict, full, motionInfo)
    seg = (tmp_path / 'meshes/seg02/snappyHexMeshDict').read_text()
    assert 'locationInMesh (0.0013 0.0021 0.0300375);' in seg
    full['GLOBAL_REFINEMENT']['LOC_IN_MESH'] = ['0', '0', '0.0033']
    with pytest.raises(SystemExit, match='lies on a cell face.*in x y'):
        writeInternal.writeSnappyInternal(templateLoc, geomDict, full, motionInfo)


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
    assert 'pRefPoint ( 0.0013 0.0021 0.0033 );' in (tmp_path / 'system/fvSolution').read_text()
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
    #turbulence fields upwind (linearUpwind took omega below zero next to the walls every step)
    assert 'div(phi,k)      Gauss upwind;' in schemes and 'div(phi,omega)  Gauss upwind;' in schemes
    assert 'pRef' not in (tmp_path / 'system/fvSolution').read_text()
    fos = (tmp_path / 'system/internalFunctionObjects').read_text()
    assert 'gasMass' in fos and 'fields          (p T rho);' in fos
    #wallHeatFlux: .dat every step, the field only at the write times
    heat = fos[fos.index('wallHeatFlux\n{'):]
    heat = heat[:heat.index('}')]
    assert 'executeControl  timeStep;' in heat and 'writeControl    writeTime;' in heat


def test_relaxation_auto_steady(tmp_path, monkeypatch, repo_root):
    #steady solvers relax by default (auto = p 0.3, U and e 0.7)
    templateLoc, full = _fullDefaults(repo_root, 'rhoSimpleFoam')
    monkeypatch.chdir(tmp_path)
    assert full['GLOBAL_SIM_CONTROL']['RELAX_P'] == ['auto']
    text = writeInternal.getRelaxation(full)
    assert 'p               0.3;' in text and 'e               0.7;' in text


def test_relaxation_only_with_outer_correctors(tmp_path, monkeypatch, repo_root):
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    monkeypatch.chdir(tmp_path)
    writeInternal.writeNumericsInternal(templateLoc, full)
    text = (tmp_path / 'system/fvSolution').read_text()
    assert 'nOuterCorrectors 1;' in text and 'relaxationFactors' not in text and '<RELAXATION>' not in text
    full['GLOBAL_SIM_CONTROL']['N_OUTER_CORR'] = ['3']
    #transient, RELAX_* auto: no under-relaxation (it breaks the gas energy balance)
    writeInternal.writeNumericsInternal(templateLoc, full)
    text = (tmp_path / 'system/fvSolution').read_text()
    assert 'nOuterCorrectors 3;' in text and 'relaxationFactors' not in text and '<RELAXATION>' not in text
    for key, value in (('RELAX_P', '0.3'), ('RELAX_U', '0.7'), ('RELAX_E', '0.7')):
        full['GLOBAL_SIM_CONTROL'][key] = [value]
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
    #moving meshes always compile it (slidingWallVelocity on the static walls), static cases do not
    os.makedirs(tmp_path / 'static')
    case = _writeScripts(tmp_path / 'static', monkeypatch, repo_root, None)
    assert 'wmake' not in (case / 'runScript').read_text()
    case = _writeScripts(tmp_path, monkeypatch, repo_root, ['motion_a'])
    assert 'wmake libso' in (case / 'runScript').read_text()
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    full['MESH_MOTION']['DRIVE'] = ['crank']
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'])
    run = (case / 'runScript').read_text()
    assert "cd '/opt/caseSetup/foamLibs/stirlingMotion' && wmake libso" in run
    assert run.index('wmake libso') < run.index('createZeroDirectory')
    motion = (case / 'motionCheckScript').read_text()
    assert motion.index('wmake libso') < motion.index('moveDynamicMesh -parallel')


def test_crank_dict_porous_zones(tmp_path, monkeypatch, repo_root):
    #a porous zone moving with a body is listed under it: its gas drag is part of the body force
    monkeypatch.chdir(tmp_path)
    os.makedirs('constant')
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    full['MESH_MOTION']['DRIVE'] = ['crank']
    crank = writeInternal.getCrank(full)
    disp = types.SimpleNamespace(name='displacer', source='MOVE-displacer', movedPids=['top', 'bottom'])
    piston = types.SimpleNamespace(name='piston', source='WALL-chamber', movedPids=['crown'])
    zone = types.SimpleNamespace(name='displacerTop', body=disp)
    motionInfo = {'bodies': [disp, piston], 'zones': [zone], 'porOwners': {'POR-matrix': disp}}
    writeInternal.writeCrankDict(templateLoc, motionInfo, crank, 0.0, 'constant/dynamicMeshDict', '<constant>/tables')
    text = open('constant/dynamicMeshDict').read()
    displacer = text[text.index('displacer\n'):text.index('piston\n')]
    assert 'porousZones (POR-matrix);' in displacer and text.count('porousZones') == 1


def test_warm_start_script(tmp_path, monkeypatch, repo_root):
    #warm start: remeshMassFix compiled and run before the solve, only while the case has nothing but the 0 time
    case = _writeScripts(tmp_path, monkeypatch, repo_root, ['motion_a'])
    assert 'warmStart' not in (case / 'runScript').read_text() and 'remeshMassFix' not in (case / 'runScript').read_text()
    templateLoc = os.path.join(repo_root, 'setupTemplates', 'internal')
    full = setupCore.getGlobalDefaults(templateLoc)
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', ['motion_a'], None, None,
                                      {'zone': 'POR-matrix', 'T': [600, 300]})
    run = (case / 'runScript').read_text()
    assert "remeshMassFix -parallel -time 0 -warmStart POR-matrix -warmStartT '(600 300)'" in run
    assert "cd '/opt/caseSetup/foamLibs/remeshMassFix' && wmake" in run
    assert run.index('createZeroDirectory') < run.index('-warmStart') < run.index('running rhoPimpleFoam')
    assert 'ls -d processor0/[0-9]* | wc -l) = 1' in run


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
    assert 'pRefPoint ( 0.0013 0.0021 0.0033 );' in solution and 'pMinFactor' not in solution
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


def test_end_time_not_cyclic(tmp_path, monkeypatch, repo_root):
    #END_TIME: transient run of a set length, write every WRITE_INTERVAL (blank = END_TIME/20), not with MOVE bodies
    monkeypatch.chdir(tmp_path)
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    assert writeInternal.getEndTime(full) is None
    assert writeInternal.getRunTimes(full) == pytest.approx((3/10.0, 1/200.0))
    full['GLOBAL_CONTROL']['END_TIME'] = ['2.5']
    assert writeInternal.getRunTimes(full) == pytest.approx((2.5, 0.125))
    full['GLOBAL_CONTROL']['WRITE_INTERVAL'] = ['0.5']
    writeInternal.writeControlDictInternal(templateLoc, full)
    control = (tmp_path / 'system/controlDict').read_text()
    assert 'endTime         2.5;' in control and 'writeInterval   0.5;' in control
    assert writeInternal.addSetupInfo(full, {})['SETUP_INFO']['CYCLIC'] == ['no']
    full['MOVE-piston'] = {'MOTION_FILE': ['m.csv']}
    with pytest.raises(SystemExit):
        writeInternal.getEndTime(full)


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
    assert zone['heat'] == 'fixed' and zone['ha'] > 0 and zone['Ts'] == [550, 550] and zone['Cs'] == 0
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


def test_wire_mesh_local_coefficients(tmp_path, monkeypatch, repo_root):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'constant').mkdir()
    #static wire mesh, COEFFS local (default): coded drag with the Forchheimer term scaled by the local Re (no
    #explicitPorositySource), h*a per cell in the heat source, Umat kept at zero
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoPimpleFoam', RESISTANCE='wireMesh', POROSITY='0.9',
                                               WIRE_DIA='1e-4', REF_VELOCITY='1', HEAT='solved', SOLID_T='600 300')
    zone = writeInternal.getPorous(geomDict, full)[0]
    assert zone['wire']['local'] and zone['haLocal'] and zone['wire']['reRef'] > 0
    writeInternal.writeFvOptionsInternal([zone], full)
    text = (tmp_path / 'constant/fvOptions').read_text()
    assert 'explicitPorositySource' not in text and 'type            vectorCodedSource;' in text
    assert 'const scalar fScale = pow(ReC/ReRef, -0.103);' in text
    assert 'const tensor C(muC*Dt + (0.5*rhoC*mag(Ur)*fScale)*Ft);' in text
    assert 'if (false && m.time().timeIndex() != lastIndex)' in text
    assert 'pow(epsW, 1.79)' in text and 'hEff[c] = (Cs > 0 ? haC*Cs/(Cs + haC*dt) : haC);' in text
    assert 'Qs[c] = hEff[c]*(TsOld[c] - T1);' in text and 'const scalar ha = ' not in text
    assert text.count('#{') == text.count('#}')
    #the per-cell h*a at the reference velocity is the setup value (same correlation)
    props = writeInternal.getGasProperties(full)
    wire = zone['wire']
    dh = wire['eps']*wire['dw']/(1 - wire['eps'])
    k = props['MU']*props['CP']/props['PR']
    ha = (1 + 0.99*(wire['reRef']*props['PR'])**0.66)*wire['eps']**1.79*k/dh*4*(1 - wire['eps'])/wire['dw']
    assert ha == pytest.approx(zone['ha'], rel=1e-9)
    #COEFFS ref: the constants at REF_VELOCITY as before
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoPimpleFoam', RESISTANCE='wireMesh', POROSITY='0.9',
                                               WIRE_DIA='1e-4', REF_VELOCITY='1', HEAT='solved', COEFFS='ref')
    zone = writeInternal.getPorous(geomDict, full)[0]
    assert not zone['wire']['local'] and not zone['haLocal']
    writeInternal.writeFvOptionsInternal([zone], full)
    text = (tmp_path / 'constant/fvOptions').read_text()
    assert 'type            explicitPorositySource;' in text and 'const scalar ha = ' in text and 'fScale' not in text
    #given HA: constant h*a, drag still local
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoPimpleFoam', RESISTANCE='wireMesh', POROSITY='0.9',
                                               WIRE_DIA='1e-4', HEAT='solved', HA='5e5')
    zone = writeInternal.getPorous(geomDict, full)[0]
    assert zone['wire']['local'] and not zone.get('haLocal') and zone['ha'] == 5e5
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoPimpleFoam', RESISTANCE='wireMesh', COEFFS='fast')
    with pytest.raises(SystemExit):
        writeInternal.getPorous(geomDict, full)


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
    assert 'const word zoneName("POR-matrix");' in text and text.count('#{') == text.count('#}')
    #gas and matrix implicit together: the same hEff*(TsOld - T1) heats the gas and cools the matrix
    assert 'hEff = (Cs > 0 ? ha*Cs/(Cs + ha*dt) : ha)' in text and 'Qs[c] = hEff*(TsOld[c] - T1);' in text
    assert 'Ts[c] = TsOld[c] - Qs[c]*dt/Cs;' in text and 'codeCorrect\n        #{' in text
    assert 'const scalar TinitLow = 600;' in text and 'const scalar TinitHigh = 600;' in text
    #two SOLID_T values: linear start from the -z to the +z end
    templateLoc, full2, geomDict2 = _porDefaults(repo_root, 'rhoPimpleFoam', HEAT='solved', HA='2e4', SOLID_T='600 300')
    porous2 = writeInternal.getPorous(geomDict2, full2)
    assert porous2[0]['Ts'] == [600, 300]
    text2 = writeInternal.porousHeatSource(porous2[0], 718.0)
    assert 'const scalar TinitLow = 600;' in text2 and 'const scalar TinitHigh = 300;' in text2
    assert 'typeHeaderOk<volScalarField>(true)' in text2
    #function objects: zone averages incl. the matrix T, heat into the gas
    writeInternal.writeFunctionObjectsInternal({}, full, None, porous)
    fos = (tmp_path / 'system/internalFunctionObjects').read_text()
    assert 'porHeatToGas_POR-matrix' in fos and '(Qs_POR-matrix);' in fos and '(p T Ts_POR-matrix);' in fos
    #DISSIPATION no (default): no extra energy sources
    assert 'DragHeat' not in text and 'viscousDissipation' not in text
    #DISSIPATION yes: drag work of the static zone and viscous dissipation into e, logged
    full['GLOBAL_SIM_CONTROL']['DISSIPATION'] = ['yes']
    writeInternal.writeFvOptionsInternal(porous, full)
    text = (tmp_path / 'constant/fvOptions').read_text()
    writeInternal.writeFunctionObjectsInternal({}, full, None, porous)
    fos = (tmp_path / 'system/internalFunctionObjects').read_text()
    assert 'POR_matrixDragHeat' in text and 'heat[c] = Ur & (C & Ur);' in text and 'vector(Zero)' in text
    assert 'type            viscousDissipation;' in text
    assert 'const tensor Dt(tensor(10000000, 0, 0, 0, 10000000, 0, 0, 0, 10000000));' in text
    assert '(porDragHeat_POR-matrix);' in fos and 'porDrag_POR-matrix' not in fos.replace('porDragHeat', '')
    #no porous zone: viscous dissipation only, DISSIPATION no: fvOptions removed
    writeInternal.writeFvOptionsInternal([], full)
    text = (tmp_path / 'constant/fvOptions').read_text()
    assert 'viscousDissipation' in text and 'Porosity' not in text
    full['GLOBAL_SIM_CONTROL']['DISSIPATION'] = ['no']
    writeInternal.writeFvOptionsInternal([], full)
    assert not (tmp_path / 'constant/fvOptions').exists()
    writeInternal.writeFvOptionsInternal(porous, full)
    text = (tmp_path / 'constant/fvOptions').read_text()
    assert 'viscousDissipation' not in text and 'DragHeat' not in text


def test_porous_moving_drag(tmp_path, monkeypatch, repo_root):
    #a zone moving with a body: coded drag on U - Umat instead of explicitPorositySource, force logged and in the crank
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'constant').mkdir()
    (tmp_path / 'system').mkdir()
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoPimpleFoam', HEAT='solved', HA='2e4', SOLID_T='600 300',
                                              E1='0 0 1', E2='1 0 0', D='1 2 3', F='4 5 6')
    full['GLOBAL_SIM_CONTROL']['DISSIPATION'] = ['yes']
    porous = writeInternal.getPorous(geomDict, full)
    body = types.SimpleNamespace(name='displacer')
    writeInternal.writeFvOptionsInternal(porous, full, {'POR-matrix': body})
    text = (tmp_path / 'constant/fvOptions').read_text()
    assert 'explicitPorositySource' not in text and 'type            vectorCodedSource;' in text
    assert 'fields          (U);' in text and 'codeAddSupRho\n        #{' in text and 'codeAddSup #{ #};' in text
    assert 'Umat[c] = d/(cellPoints[c].size()*dt);' in text and 'porDrag[c] = (C & Ur);' in text
    assert 'diag[c] -= iso*V[c];' in text and 'source[c] -= ((C & Umat[c]) - ((C - iso*I) & U[c]))*V[c];' in text
    #d along e1 = z (1), e2 = x (2), e3 = e1 x e2 = y (3)
    assert 'const tensor Dt(tensor(2, 0, 0, 0, 3, 0, 0, 0, 1));' in text
    assert 'const tensor Ft(tensor(5, 0, 0, 0, 6, 0, 0, 0, 4));' in text
    assert 'm.lookupObject<volVectorField>("Umat_" + zoneName)[c]' in text and text.count('#{') == text.count('#}')
    motionInfo = {'bodies': [], 'zones': [], 'porOwners': {'POR-matrix': body}}
    writeInternal.writeFunctionObjectsInternal({}, full, motionInfo, porous)
    fos = (tmp_path / 'system/internalFunctionObjects').read_text()
    assert '(porDrag_POR-matrix);' in fos and '(porDragHeat_POR-matrix);' in fos
    #incompressible: kinematic drag (nu), force scaled by the gas density, no energy sources
    full['GLOBAL_SIM_CONTROL']['SOLVER'] = ['pimpleFoam']
    porous = writeInternal.getPorous(geomDict, full)
    writeInternal.writeFvOptionsInternal(porous, full, {'POR-matrix': body})
    text = (tmp_path / 'constant/fvOptions').read_text()
    rho = writeInternal.getGasDensity(full)
    assert 'codeAddSup\n        #{' in text and 'codeAddSupRho #{ #};' in text and 'const scalar rhoC = 1;' in text
    assert 'porDrag[c] = %.9g*(C & Ur);' % (rho) in text
    assert 'viscousDissipation' not in text and 'scalarCodedSource' not in text


def test_fiber_coeffs_blend():
    #felt correlation where tested, single fibres at steel-wool porosity, log-linear between; Darcy limit of the felt
    #form = 192/(2 eps dh^2)
    gas = {'RHO': 1.0, 'MU': 2e-5, 'CP': 1000.0, 'PR': 0.7}
    c, ha, cF, cC, hF, hC = writeInternal.fiberCoeffs(0.8, 5e-5, 1e-9, gas)
    dh = 0.8*5e-5/0.2
    assert c == pytest.approx(cF) and ha == pytest.approx(hF) and cF == pytest.approx(2e-5*96/(0.8*dh**2), rel=1e-6)
    c, ha, cF, cC, hF, hC = writeInternal.fiberCoeffs(0.99, 5e-5, 0.5, gas)
    assert c == pytest.approx(cC) and ha == pytest.approx(hC) and hC > 3*hF
    c, ha, cF, cC, hF, hC = writeInternal.fiberCoeffs(0.905, 5e-5, 0.5, gas)
    assert c == pytest.approx(cF**0.5*cC**0.5) and ha == pytest.approx(hF**0.5*hC**0.5)


def test_porous_fiber_deforming(tmp_path, monkeypatch, repo_root):
    #steel wool: random fibres with a fixed mass, deforming with the mesh: fibre fraction field, per-cell drag and h*a
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'constant').mkdir()
    (tmp_path / 'system').mkdir()
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoPimpleFoam', RESISTANCE='randomFiber', WIRE_DIA='5e-5',
                                              SOLID_RHO='7900', SOLID_CP='500', SOLID_MASS='0.02', HEAT='solved',
                                              SOLID_T='550', MOVES_WITH='mesh')
    zone = writeInternal.getPorous(geomDict, full)[0]
    assert zone['fiber']['mass'] == 0.02 and zone['fiber']['cs'] == 500 and zone['Cs'] == 0
    writeInternal.writeFvOptionsInternal([zone], full, {})
    text = (tmp_path / 'constant/fvOptions').read_text()
    assert 'explicitPorositySource' not in text and text.count('#{') == text.count('#}')
    assert 'const scalar solidVol = %.9g;' % (0.02/7900) in text
    assert 'solidFrac[zoneCells[i]] *= V0[zoneCells[i]]/V[zoneCells[i]];' in text
    assert 'diag[c] -= cDrag*V[c];' in text and 'source[c] -= cDrag*Umat[c]*V[c];' in text
    assert 'const scalar Cs = solidFrac[c]*rhoCs;' in text and 'const scalar rhoCs = %.9g;' % (7900*500) in text
    assert 'pow(eps, 2.61)' in text and 'muC*96/(eps*sqr(dh))' in text
    writeInternal.writeFunctionObjectsInternal({}, full, {'bodies': [], 'zones': [], 'porOwners': {}}, [zone])
    fos = (tmp_path / 'system/internalFunctionObjects').read_text()
    assert '(solidFrac_POR-matrix);' in fos and '(porDrag_POR-matrix);' in fos and '(Qs_POR-matrix);' in fos
    #heat through the fibre bed (conduction + radiation) and from the hot wall it touches
    assert 'kEff[c] = ks*phiC/3 + (radiation ? 16*sigma*pow3(Ts[c])*dw/(3*phiC) : 0);' in text
    assert 'contactKeys.append' not in text and 'const label nSub' in text
    full['POR-matrix']['WALL_CONTACT_PIDS'] = ['Hot']
    full['POR-matrix']['WALL_CONTACT'] = ['0.5']
    zone = writeInternal.getPorous(geomDict, full)[0]
    writeInternal.writeFvOptionsInternal([zone], full, {})
    text = (tmp_path / 'constant/fvOptions').read_text()
    assert 'contactKeys.append(word("hot"));' in text and 'const scalar contact = 0.5;' in text
    assert 'Qwall[fc[i]] += q/V[fc[i]]/nSub;' in text and text.count('#{') == text.count('#}')
    writeInternal.writeFunctionObjectsInternal({}, full, {'bodies': [], 'zones': [], 'porOwners': {}}, [zone])
    assert '(Qwall_POR-matrix);' in (tmp_path / 'system/internalFunctionObjects').read_text()
    kc, kr = writeInternal.fiberConductivity(0.01, 5e-5, 50.0, 600.0)
    assert kc == pytest.approx(50*0.01/3) and kr == pytest.approx(16*5.670374e-8*600**3*5e-5/(3*0.01))
    #deforming needs a fibre mass, and only random fibres can deform
    full['POR-matrix']['WALL_CONTACT_PIDS'] = ['']
    full['POR-matrix']['SOLID_MASS'] = ['']
    with pytest.raises(SystemExit):
        writeInternal.getPorous(geomDict, full)
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoPimpleFoam', RESISTANCE='wireMesh', POROSITY='0.8',
                                              WIRE_DIA='1e-4', REF_VELOCITY='2', MOVES_WITH='mesh')
    with pytest.raises(SystemExit):
        writeInternal.getPorous(geomDict, full)


def test_warm_start(tmp_path, monkeypatch, repo_root):
    monkeypatch.chdir(tmp_path)
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoPimpleFoam', HEAT='solved', HA='2e4', SOLID_T='600 300')
    porous = writeInternal.getPorous(geomDict, full)
    assert writeInternal.getWarmStart(full, porous) is None
    full['GLOBAL_MATERIAL']['INIT_T_FIELD'] = ['matrix']
    assert writeInternal.getWarmStart(full, porous) == {'zone': 'POR-matrix', 'T': [600, 300]}
    full['GLOBAL_MATERIAL']['INIT_T_FIELD'] = ['POR-matrix']
    assert writeInternal.getWarmStart(full, porous)['zone'] == 'POR-matrix'
    full['GLOBAL_MATERIAL']['INIT_T_FIELD'] = ['POR-other']
    with pytest.raises(SystemExit):
        writeInternal.getWarmStart(full, porous)
    #uniform matrix: 'matrix' has no linear zone to follow
    templateLoc, full, geomDict = _porDefaults(repo_root, 'rhoPimpleFoam', HEAT='solved', HA='2e4', SOLID_T='600')
    full['GLOBAL_MATERIAL']['INIT_T_FIELD'] = ['matrix']
    with pytest.raises(SystemExit):
        writeInternal.getWarmStart(full, writeInternal.getPorous(geomDict, full))


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
    assert post.index("imagesToMovies.py' > log.imagesToMovies") < post.index("copyBackInternal.py') > log.copyBackInternal")
    full['GLOBAL_COMPUTE_SETUP']['COPY_BACK'] = ['False']
    writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', None)
    assert "copyBackInternal.py') > log.copyBackInternal" not in (tmp_path / 'postScript').read_text()
    full['GLOBAL_COMPUTE_SETUP']['COPY_BACK'] = ['sometimes']
    with pytest.raises(SystemExit):
        writeInternal.makeScriptsInternal(templateLoc, full, '/opt/caseSetup', None)


def test_static_walls_slide_in_a_moving_mesh(tmp_path, monkeypatch, repo_root):
    #static noSlip walls of a moving-mesh case get slidingWallVelocity (no leak through walls whose points slide),
    #the controlDict loads the library; without MOVE bodies they stay stationary
    geomDict, wallGroups = _pipeCase(tmp_path, monkeypatch)
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    writeInternal.writeCasePropertiesInternal(templateLoc, geomDict, full, wallGroups)
    wall = [l for l in (tmp_path / 'system/caseProperties').read_text().splitlines() if 'category wall;' in l][0]
    assert 'motion stationary;' in wall
    writeInternal.writeControlDictInternal(templateLoc, full)
    assert 'stirlingMotion' not in (tmp_path / 'system/controlDict').read_text()
    full['MOVE-displacer'] = {'FROM': ['self']}
    writeInternal.writeCasePropertiesInternal(templateLoc, geomDict, full, wallGroups)
    wall = [l for l in (tmp_path / 'system/caseProperties').read_text().splitlines() if 'category wall;' in l][0]
    assert 'motion slidingMesh;' in wall
    writeInternal.writeControlDictInternal(templateLoc, full)
    control = (tmp_path / 'system/controlDict').read_text()
    assert 'libs            (stirlingMotion);' in control and control.index('libs') < control.index('functions')
    options = open(os.path.join(repo_root, 'zeroTemplates', 'boundaryConditions', 'fluid', 'wallOptions')).read()
    assert 'slidingMesh' in options and 'slidingWallVelocity' in options
    #resolved (low) walls: omega from omegaWallFunction, not the shared lowReynolds option's fixed 1e-8
    geomDict['WALL-pipe.stl']['wallmodel'] = 'low'
    writeInternal.writeCasePropertiesInternal(templateLoc, geomDict, full, wallGroups)
    wall = [l for l in (tmp_path / 'system/caseProperties').read_text().splitlines() if 'category wall;' in l][0]
    assert 'wallFunction lowReynoldsInternal;' in wall
    block = options[options.index('lowReynoldsInternal'):]
    block = block[:block.index('nuTilda')]
    assert 'omegaWallFunction' in block and '1e-8' not in block


def test_slip_wall_on_moving_mesh_warns(tmp_path, monkeypatch, repo_root, capsys):
    #slip has U.n = 0 on faces that sweep volume as their points slide: gas leaks, so the setup says so
    geomDict, wallGroups = _pipeCase(tmp_path, monkeypatch)
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    walls = [g for groups in wallGroups.values() for g in groups if not g['flow']]
    for group in walls:
        group['U'] = 'slip'
    capsys.readouterr()
    writeInternal.writeCasePropertiesInternal(templateLoc, geomDict, full, wallGroups)
    assert 'leak' not in capsys.readouterr().out
    full['MOVE-displacer'] = {'FROM': ['self']}
    writeInternal.writeCasePropertiesInternal(templateLoc, geomDict, full, wallGroups)
    assert 'is slip on a moving mesh: its faces leak gas' in capsys.readouterr().out


def test_from_body_own_bc(capsys):
    #a FROM body's faces are one group with the body's own U_DEFAULT / T_DEFAULT (blank = the file's defaults)
    geomDict = {'WALL-cyl.stl': {}}
    pidDict = {'WALL-cyl': ['MOVE-crown', 'WALL-liner', 'WALL-bottom']}
    full = {'WALL-cyl': {'U_DEFAULT': ['noSlip'], 'T_DEFAULT': ['fixed', '295'], 'BC_GROUPS': ['']},
            'MOVE-crown': {'FROM': ['WALL-cyl'], 'ZONE_PIDS': ['MOVE-crown'], 'U_DEFAULT': [''], 'T_DEFAULT': ['adiabatic']}}
    groups = dict((g['name'], g) for g in writeInternal.getWallGroups(geomDict, full, pidDict, {'WALL-cyl': ['MOVE-crown']})['WALL-cyl'])
    assert groups['MOVE-crown']['pids'] == ['MOVE-crown'] and groups['MOVE-crown']['T'] == ['adiabatic']
    assert groups['MOVE-crown']['U'] == 'noSlip' and groups['MOVE-crown']['moving'] == ['MOVE-crown']
    assert groups['WALL-cyl_default']['pids'] == ['WALL-liner', 'WALL-bottom'] and groups['WALL-cyl_default']['T'] == ['fixed', 295.0]
    out = capsys.readouterr().out
    assert 'movingWall' in out
    #a T different from the file's default is pointed out (its meaning changed on 2026-11-24)
    assert 'NOTE: [MOVE-crown] faces use its own T_DEFAULT (adiabatic), not the WALL-cyl default (fixed 295.0)' in out
    #blank T: the file's default
    full['MOVE-crown']['T_DEFAULT'] = ['']
    groups = dict((g['name'], g) for g in writeInternal.getWallGroups(geomDict, full, pidDict, {'WALL-cyl': ['MOVE-crown']})['WALL-cyl'])
    assert groups['MOVE-crown']['T'] == ['fixed', 295.0]
    assert 'NOTE' not in capsys.readouterr().out
    #the same faces in a group of the wall file as well: error (one place for U and T)
    full['WALL-cyl']['BC_GROUPS'], full['WALL-cyl']['crown'] = ['crown'], ['PIDS', 'crown', 'T', 'adiabatic']
    with pytest.raises(SystemExit, match='set its U and T in one place'):
        writeInternal.getWallGroups(geomDict, full, pidDict, {'WALL-cyl': ['MOVE-crown']})
    #BC_GROUPS in the FROM body's own section are not used
    full['WALL-cyl']['BC_GROUPS'] = ['']
    full['MOVE-crown']['BC_GROUPS'] = ['face']
    capsys.readouterr()
    writeInternal.getWallGroups(geomDict, full, pidDict, {'WALL-cyl': ['MOVE-crown']})
    assert 'WARNING: [MOVE-crown] BC_GROUPS face is not used (FROM body)' in capsys.readouterr().out


#----- mesh groups (surface level and wall layers per PID)

def test_parseMeshGroupLine():
    group = writeInternal.parseMeshGroupLine(['PIDS', 'hot', 'cold', 'LEVEL', '2', '3', 'LAYERS', '6', 'FIRST', '1e-4',
                                              'FINAL', '0.5'], 'g')
    assert group['keywords'] == ['hot', 'cold'] and group['level'] == (2, 3) and group['layers'] == 6
    #FIRST without SIZES: mixed (absolute first layer, final relative)
    assert group['sizes'] == {'FIRST': 1e-4, 'FINAL': 0.5} and group['mode'] == 'mixed'
    assert writeInternal.parseMeshGroupLine(['PIDS', 'hot', 'LEVEL', '3'], 'g')['level'] == (3, 3)
    auto = writeInternal.parseMeshGroupLine(['PIDS', 'hot', 'LAYERS', 'auto'], 'g')
    assert auto['layers'] == 'auto' and auto['mode'] == 'mixed'
    assert writeInternal.parseMeshGroupLine(['PIDS', 'hot', 'LAYERS', '3'], 'g')['mode'] == 'relative'
    absolute = writeInternal.parseMeshGroupLine('PIDS hot LAYERS 5 SIZES absolute FIRST 1e-4 EXPANSION 1.2 MIN 5e-5'.split(), 'g')
    assert absolute['mode'] == 'absolute' and absolute['min'] == 5e-5
    assert writeInternal.parseMeshGroupLine(['PIDS', 'side', 'LAYERS', '-1'], 'g')['layers'] == -1


@pytest.mark.parametrize('line, message', [
    ('hot LAYERS 3', 'must start with PIDS'),
    ('PIDS hot', 'sets nothing'),
    ('PIDS hot LEVEL 3 2', 'max below min'),
    ('PIDS hot LEVEL a', 'refinement levels'),
    ('PIDS hot LAYERS -2', 'number of layers, auto, 0'),
    ('PIDS hot LEVEL 2 FIRST 1e-4', 'need LAYERS'),
    ('PIDS hot LAYERS -1 FIRST 1e-4', 'adds no layers, no sizes'),
    ('PIDS hot LAYERS 4 SIZES other', 'SIZES must be one of'),
    ('PIDS hot LAYERS 4 FIRST 1e-4 EXPANSION 1.2', 'SIZES mixed takes FIRST'),
    ('PIDS hot LAYERS 4 SIZES absolute FIRST 1e-4 FINAL 1e-3 TOTAL 2e-3', 'at most two of'),
    ('PIDS hot LAYERS 4 SIZES absolute EXPANSION 1.2', 'SIZES absolute needs FIRST, FINAL or TOTAL'),
    ('PIDS hot LAYERS auto FIRST 1e-4', 'LAYERS auto sets the first layer'),
    ('PIDS hot LAYERS auto SIZES relative', 'LAYERS auto sets the first layer'),
    ('PIDS hot LAYERS 4 FIRST 0', 'one positive number')])
def test_parseMeshGroupLine_errors(line, message):
    with pytest.raises(SystemExit, match=message):
        writeInternal.parseMeshGroupLine(line.split(' '), 'g')


def _meshGroupCase(repo_root, **groups):
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    full['GLOBAL_MATERIAL']['GAS'] = ['helium']
    full['GLOBAL_MATERIAL']['CHARGE_PRESSURE'] = ['1.261e6']
    full['GLOBAL_CONTROL']['FREQ'] = ['5']
    full['WALL-chamber'] = wall(BC_GROUPS='hotEnd coldEnd', hotEnd='PIDS hot T fixed 600', coldEnd='PIDS cold T fixed 300',
                                MESH_GROUPS=' '.join(groups.keys()), GAP_PIDS='')
    for name, line in groups.items():
        full['WALL-chamber'][name] = line.split(' ')
    geomDict = {'WALL-chamber.stl': {'refinement': '1', 'scale': '1', 'layers': '0', 'expansion': '1.2'}}
    pidDict = {'WALL-chamber': ['hot', 'cold', 'sideWall']}
    wallGroups = writeInternal.getWallGroups(geomDict, full, pidDict)
    return templateLoc, full, geomDict, pidDict, wallGroups


def test_thermal_depth():
    #sqrt(2 alpha/omega), helium at 12.61 bar and 300 K, 5 Hz: alpha = mu cp/Pr/(rho cp)
    full = _material(GAS='helium', CHARGE_PRESSURE='1.261e6')
    full['GLOBAL_CONTROL'] = {'FREQ': ['5']}
    rho = 1.261e6*4.003/(8314.46*300.0)
    alpha = 1.99e-5/0.67/rho
    assert writeInternal.thermalDepth(full, 300.0) == pytest.approx(np.sqrt(2.0*alpha/(2.0*np.pi*5.0)))
    #sutherland: thicker at the hot wall (alpha ~ T^1.7)
    full['GLOBAL_MATERIAL']['TRANSPORT'] = ['sutherland']
    assert writeInternal.thermalDepth(full, 600.0) > 1.5*writeInternal.thermalDepth(full, 300.0)


@pytest.mark.parametrize('sizes', [{'FIRST': 1e-4, 'EXPANSION': 1.2}, {'FIRST': 1e-4, 'TOTAL': 7.44e-4},
                                   {'FINAL': 2.0736e-4, 'TOTAL': 7.44e-4}, {'FINAL': 2.0736e-4, 'EXPANSION': 1.2},
                                   {'TOTAL': 7.44e-4, 'EXPANSION': 1.2}, {'FIRST': 1e-4, 'FINAL': 2.0736e-4}])
def test_stack_geometry(sizes):
    #5 layers from 0.1 mm growing 1.2x: final 0.207 mm, total 0.744 mm, from any two of the four sizes
    first, final, total, ratio = writeInternal.stackGeometry(sizes, 5)
    assert first == pytest.approx(1e-4, rel=1e-3) and final == pytest.approx(2.0736e-4, rel=1e-3)
    assert total == pytest.approx(7.44e-4, rel=1e-3) and ratio == pytest.approx(1.2, rel=1e-3)


def test_layer_entries(repo_root):
    #every entry names its thicknessModel (a patch entry without it only reads the global model's keys)
    entry = lambda line, cell=0.002: writeInternal.layerEntry(writeInternal.parseMeshGroupLine(line.split(), 'g'), cell,
                                                              1.2, 0.4, 0.1)
    text, stack = entry('PIDS hot LAYERS 5 SIZES absolute FIRST 1e-4 EXPANSION 1.2')
    assert text == ('nSurfaceLayers 5; thicknessModel firstAndExpansion; relativeSizes false; expansionRatio 1.2; '
                    'firstLayerThickness 0.0001; minThickness 7.4416e-05;')
    assert stack[2] == pytest.approx(7.4416e-4)
    #absolute FIRST + FINAL: no snappy model of its own, the growth ratio between them
    text, stack = entry('PIDS hot LAYERS 5 SIZES absolute FIRST 1e-4 FINAL 2.0736e-4 MIN 5e-5')
    assert 'thicknessModel firstAndExpansion;' in text and 'expansionRatio 1.2;' in text and 'minThickness 5e-05;' in text
    #relative: GEOM-like, sizes x the cell; only EXPANSION given -> FINAL_LAYER_THICKNESS
    text, stack = entry('PIDS hot LAYERS 3 EXPANSION 1.3')
    assert text == ('nSurfaceLayers 3; thicknessModel finalAndExpansion; relativeSizes true; expansionRatio 1.3; '
                    'finalLayerThickness 0.4; minThickness 0.1;')
    assert stack[1] == pytest.approx(0.4*0.002)
    text, stack = entry('PIDS hot LAYERS 3 SIZES relative TOTAL 1.5')
    assert 'thicknessModel overallAndExpansion;' in text and 'thickness 1.5;' in text and 'expansionRatio 1.2;' in text
    assert stack[2] == pytest.approx(0.003)
    #mixed: absolute first, final relative
    text, stack = entry('PIDS hot LAYERS 6 FIRST 1e-4')
    assert text == ('nSurfaceLayers 6; thicknessModel firstAndRelativeFinal; relativeSizes true; '
                    'firstLayerThickness 0.0001; finalLayerThickness 0.4; minThickness 0.1;')
    assert stack[0] == pytest.approx(1e-4) and stack[1] == pytest.approx(8e-4)


def test_mesh_groups_auto_layers(repo_root, capsys):
    templateLoc, full, geomDict, pidDict, wallGroups = _meshGroupCase(
        repo_root, heat='PIDS hot cold LAYERS auto', sides='PIDS side LEVEL 2 3')
    groups = dict((g['name'], g) for g in writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)['WALL-chamber'])
    heat = groups['heat']
    assert heat['pids'] == ['hot', 'cold']
    #the colder wall (300 K) sets the depth, the first layer is a tenth of it; growing to 0.4 x the 2 mm cell would make
    #a 3.2 mm stack, so it is capped at 2 depths (absolute sizes), growing at most 1.3x
    depth = writeInternal.thermalDepth(full, 300.0)
    first, final, total, ratio = heat['stack']
    assert heat['depth'] == pytest.approx(depth) and first == pytest.approx(depth/10.0)
    assert heat['cell'] == pytest.approx(0.002) and heat['mode'] == 'absolute'
    assert total == pytest.approx(2.0*depth) and ratio <= 1.3 + 1e-9
    assert writeInternal.stackGeometry({'FIRST': first, 'TOTAL': total}, heat['layers'] - 1)[3] > 1.3
    assert 'thicknessModel firstAndOverall; relativeSizes false;' in heat['entry']
    #finer surface cells: the stack up to 0.4 x the cell stays below the cap, mixed sizes as before
    full['WALL-chamber']['heat'] = 'PIDS hot cold LEVEL 3 LAYERS auto'.split()
    fine = dict((g['name'], g) for g in writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)['WALL-chamber'])['heat']
    assert fine['mode'] == 'mixed' and fine['stack'][1] == pytest.approx(0.4*0.0005) and fine['stack'][2] < 2.0*depth
    assert 'thicknessModel firstAndRelativeFinal;' in fine['entry']
    full['WALL-chamber']['heat'] = 'PIDS hot cold LAYERS auto'.split()
    assert groups['sides']['level'] == (2, 3) and groups['sides']['layers'] is None
    out = capsys.readouterr().out
    assert 'Mesh groups' in out and 'WALL-chamber_hot' in out and 'WALL-chamber_sideWall' in out
    #cells already finer than the first layer: no layers, written like the walls without layers
    full['WALL-chamber']['heat'] = 'PIDS hot cold LEVEL 6 LAYERS auto'.split()
    groups = dict((g['name'], g) for g in writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)['WALL-chamber'])
    assert groups['heat']['autoNone'] and groups['heat']['layers'] is None
    assert 'already fine for the thermal layer' in capsys.readouterr().out


def test_mesh_groups_errors(repo_root, capsys):
    templateLoc, full, geomDict, pidDict, wallGroups = _meshGroupCase(
        repo_root, heat='PIDS hot LAYERS 4', both='PIDS ho LEVEL 2')
    with pytest.raises(SystemExit, match='matches mesh groups heat and both'):
        writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    templateLoc, full, geomDict, pidDict, wallGroups = _meshGroupCase(repo_root, heat='PIDS top LAYERS 4')
    with pytest.raises(SystemExit, match='keyword "top" matches no PID'):
        writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    templateLoc, full, geomDict, pidDict, wallGroups = _meshGroupCase(repo_root, heat='PIDS hot LAYERS 4')
    full['WALL-chamber']['MESH_GROUPS'] = ['hotEnd']
    with pytest.raises(SystemExit, match='in BC_GROUPS and MESH_GROUPS'):
        writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    #auto needs a cyclic run
    templateLoc, full, geomDict, pidDict, wallGroups = _meshGroupCase(repo_root, heat='PIDS hot LAYERS auto')
    full['GLOBAL_CONTROL']['END_TIME'] = ['2']
    with pytest.raises(SystemExit, match='cyclic runs only'):
        writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    #layers on a gap face, layers growing too fast, final layer above the cell: warnings
    templateLoc, full, geomDict, pidDict, wallGroups = _meshGroupCase(
        repo_root, heat='PIDS side LAYERS 2 FIRST 1e-5', thick='PIDS hot LAYERS 2 SIZES absolute FINAL 0.003')
    full['WALL-chamber']['GAP_PIDS'] = ['side']
    writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    out = capsys.readouterr().out
    assert "puts layers on the GAP_PIDS faces ['sideWall']" in out and 'use more LAYERS' in out
    assert 'thicker than the 2 mm surface cell' in out


def test_mesh_groups_snappy_dict(tmp_path, monkeypatch, repo_root):
    templateLoc, full, geomDict, pidDict, wallGroups = _meshGroupCase(
        repo_root, heat='PIDS hot cold LAYERS 6 FIRST 1e-4', sides='PIDS side LEVEL 2 3')
    meshGroups = writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    full['WALL-chamber']['GAP_PIDS'] = ['side']
    monkeypatch.chdir(tmp_path)
    writeInternal.writeSnappyInternal(templateLoc, geomDict, full, None, meshGroups)
    snappy = (tmp_path / 'system/snappyHexMeshDict').read_text()
    #exact region levels before the keyword include (an exact key wins over the regex keys)
    assert 'regions { sideWall { level (2 3); } #include "snappyRefinementDict" }' in snappy
    #per-patch layers after the file-wide and the gap entries (last matching entry wins), absolute first layer with
    #the final layer relative to the cell (firstAndRelativeFinal)
    hot = ('"WALL-chamber_hot" { nSurfaceLayers 6; thicknessModel firstAndRelativeFinal; relativeSizes true; '
           'firstLayerThickness 0.0001; finalLayerThickness 0.4; minThickness 0.1; }')
    assert hot in snappy
    #walls without layers slide (-1) next to the layered ones: the file-wide GEOM 0 and the gap faces
    assert '"WALL-chamber_.*" { nSurfaceLayers -1; expansionRatio 1.2;}' in snappy
    assert '"WALL-chamber_.*side.*" { nSurfaceLayers -1; }' in snappy
    assert snappy.index('"WALL-chamber_.*" {') < snappy.index('"WALL-chamber_.*side.*"') < snappy.index(hot)
    #patches with layers for the coverage report
    report = (tmp_path / 'system/meshGroups').read_text().splitlines()
    assert report[1:] == ['WALL-chamber_hot heat 6 0.0001', 'WALL-chamber_cold heat 6 0.0001']
    #ZERO_LAYERS fixed: 0 as before; an explicit LAYERS 0 / -1 in a group is kept
    full['GLOBAL_REFINEMENT']['ZERO_LAYERS'] = ['fixed']
    full['WALL-chamber']['sides'] = 'PIDS side LAYERS -1'.split()
    meshGroups = writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    writeInternal.writeSnappyInternal(templateLoc, geomDict, full, None, meshGroups)
    snappy = (tmp_path / 'system/snappyHexMeshDict').read_text()
    assert '"WALL-chamber_.*" { nSurfaceLayers 0; expansionRatio 1.2;}' in snappy
    assert '"WALL-chamber_sideWall" { nSurfaceLayers -1; }' in snappy
    #no mesh-group layers: GEOM 0 stays 0 (the dicts of existing cases do not change)
    full['GLOBAL_REFINEMENT']['ZERO_LAYERS'] = ['slip']
    full['WALL-chamber']['heat'] = 'PIDS hot cold LEVEL 2'.split()
    meshGroups = writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    writeInternal.writeSnappyInternal(templateLoc, geomDict, full, None, meshGroups)
    snappy = (tmp_path / 'system/snappyHexMeshDict').read_text()
    assert '"WALL-chamber_.*" { nSurfaceLayers 0; expansionRatio 1.2;}' in snappy
    full['GLOBAL_REFINEMENT']['ZERO_LAYERS'] = ['maybe']
    with pytest.raises(SystemExit, match='ZERO_LAYERS must be slip or fixed'):
        writeInternal.writeSnappyInternal(templateLoc, geomDict, full, None, meshGroups)


def test_layer_report(tmp_path, repo_root, capsys):
    import importlib.util
    spec = importlib.util.spec_from_file_location('layerReport', os.path.join(repo_root, 'meshUtilities', 'layerReport.py'))
    layerReport = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(layerReport)
    (tmp_path / 'meshGroups').write_text('#patch group layers first\nWALL-chamber_hot heat 8 0.0001\n'
                                         'WALL-chamber_cold heat 8 0.0001\nWALL-chamber_crown heat 8 0.0001\n')
    (tmp_path / 'log.snappyHexMesh').write_text(
        'patch        faces        layers        overall thickness\n'
        '                      target   mesh     [m]       [%]\n'
        '-----        -----    -----    ----     ---       ---\n'
        'WALL-chamber_hot  1024     8        5.05     0.00217   66.9\n'
        'WALL-chamber_cold 1024     8        1.2      0.0004    12.0\n'
        'WALL-chamber_side 4096     2        1.65     0.00118   64.7\n'
        '\nLayer mesh : cells:29924\n')
    layerReport.main([str(tmp_path / 'log.snappyHexMesh'), str(tmp_path / 'meshGroups'), '--name', 'seg02'])
    out = capsys.readouterr().out
    lines = dict((patch, [line for line in out.splitlines() if patch + ' ' in line][0])
                 for patch in ('WALL-chamber_hot', 'WALL-chamber_cold'))
    assert 'seg02: layers of the mesh groups' in out
    assert 'WARNING' not in lines['WALL-chamber_hot'] and '5.05' in lines['WALL-chamber_hot']
    assert 'WARNING' in lines['WALL-chamber_cold']
    assert 'WALL-chamber_crown (heat) not in the layer table' in out
    #no mesh groups with layers: nothing printed
    (tmp_path / 'meshGroups').write_text('#patch group layers first\n')
    layerReport.main([str(tmp_path / 'log.snappyHexMesh'), str(tmp_path / 'meshGroups')])
    assert capsys.readouterr().out == ''


def test_scripts_run_layer_report(tmp_path, monkeypatch, repo_root):
    case = _writeScripts(tmp_path, monkeypatch, repo_root, ['motion_a'])
    run = (case / 'runScript').read_text()
    assert run.index('$RUN snappyHexMesh') < run.index("layerReport.py' log.snappyHexMesh system/meshGroups;")
    #negative cell volumes stop the meshing (case mesh and every remeshing segment), other checkMesh failures do not
    assert "grep -q 'Zero or negative cell volume' log.checkMesh" in run
    import json
    for setup in ('internal', 'internal-relaxed'):
        with open(os.path.join(repo_root, 'setupTemplates', setup, 'defaultCluster', 'slurm', 'clusterDict')) as f:
            remesh = json.load(f)['run']['remeshMesh']
        assert "grep -q 'Zero or negative cell volume' meshes/$s/log.checkMesh" in remesh
        assert 'ERROR: negative cell volumes in meshes/$s' in remesh


#----- mesh groups: EDGE_LEVEL (feature edges of some PIDs, e.g. the hole rims of the displacer faces)

def test_edge_level_parse():
    group = writeInternal.parseMeshGroupLine('PIDS bottom top EDGE_LEVEL 4'.split(), 'g')
    assert group['edgeLevel'] == 4 and group['layers'] is None and group['level'] is None
    with pytest.raises(SystemExit, match='EDGE_LEVEL takes one refinement level'):
        writeInternal.parseMeshGroupLine('PIDS bottom EDGE_LEVEL 4 5'.split(), 'g')


def test_edge_level_gap_warning(repo_root, capsys):
    #rims opening into a gap refined to level 3: only level 4 removes the skewed faces
    templateLoc, full, geomDict, pidDict, wallGroups = _meshGroupCase(repo_root, rims='PIDS hot EDGE_LEVEL 3')
    motionInfo = {'gapRefs': [{'level': 3}]}
    writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups, motionInfo)
    assert 'EDGE_LEVEL 3 is not above the gap refinement level 3' in capsys.readouterr().out
    full['WALL-chamber']['rims'] = 'PIDS hot EDGE_LEVEL 4'.split()
    writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups, motionInfo)
    assert 'not above the gap refinement' not in capsys.readouterr().out


def test_edge_level_geometry_and_snappy(tmp_path, monkeypatch, repo_root):
    templateLoc, full, geomDict, pidDict, wallGroups = _meshGroupCase(repo_root, rims='PIDS hot cold EDGE_LEVEL 4')
    msh = tmp_path / 'MSH'
    msh.mkdir()
    (msh / 'WALL-chamber.stl').write_text(STL)
    meshGroups = writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    edgeRefs = writeInternal.getEdgeRefs(meshGroups, geomDict, str(msh))
    assert [e['file'] for e in edgeRefs] == ['EDGE-WALL-chamber-rims.stl']
    case = tmp_path / 'case'
    case.mkdir()
    monkeypatch.chdir(case)
    #an EDGE file of an earlier setup is removed
    (case / 'constant/triSurface').mkdir(parents=True)
    (case / 'constant/triSurface/EDGE-WALL-chamber-old.stl').write_text('solid a\nendsolid a\n')
    #moving PIDs are shifted like their geometry (cold moved by +0.5 here)
    writeInternal.writeEdgeGeometry(edgeRefs, shifts={'WALL-chamber': {'cold': (0.0, 0.0, 0.5)}})
    assert sorted(os.listdir(case / 'constant/triSurface')) == ['EDGE-WALL-chamber-rims.stl']
    tris = geometryIO.readSurfaceTriangles(str(case / 'constant/triSurface/EDGE-WALL-chamber-rims.stl'))
    assert list(tris) == ['hot', 'cold']
    assert np.allclose(tris['hot'][:, :, 2], 0.0) and np.allclose(tris['cold'][:, :, 2], 1.5)
    writeInternal.writeSnappyInternal(templateLoc, geomDict, full, None, meshGroups)
    snappy = (case / 'system/snappyHexMeshDict').read_text()
    assert '{ file "EDGE-WALL-chamber-rims.eMesh"; scale 1; level 4; }' in snappy
    extract = (case / 'system/surfaceFeatureExtractDict').read_text()
    assert 'EDGE-WALL-chamber-rims.stl { extractionMethod extractFromSurface;' in extract
    assert 'openEdges yes;' in extract


def test_mesh_groups_from_self_no_warning(repo_root, capsys):
    #MESH_GROUPS of a body with its own file (FROM = self) are used; only a body made of another file's PIDs warns
    templateLoc, full, geomDict, pidDict, wallGroups = _meshGroupCase(repo_root, heat='PIDS hot LEVEL 2')
    full['MOVE-displacer'] = {'FROM': ['self'], 'MESH_GROUPS': ['rims']}
    writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    assert 'is not used (FROM body)' not in capsys.readouterr().out
    full['MOVE-displacer']['FROM'] = ['WALL-chamber']
    writeInternal.getMeshGroups(geomDict, full, pidDict, wallGroups)
    assert 'MESH_GROUPS rims is not used (FROM body)' in capsys.readouterr().out


def test_field_min_max_function_object(tmp_path, monkeypatch, repo_root):
    #fieldMinMax of every solved field with its cell location, each step: where a run starts to blow up
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    assert writeInternal.minMaxFields(full) == ['U', 'p', 'T', 'rho']
    full['GLOBAL_SIM_CONTROL']['TURB_MODEL'] = ['kOmegaSST']
    assert writeInternal.minMaxFields(full) == ['U', 'p', 'T', 'rho', 'k', 'omega', 'nut', 'alphat']
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'system').mkdir()
    writeInternal.writeFunctionObjectsInternal({}, full)
    fos = (tmp_path / 'system/internalFunctionObjects').read_text()
    block = fos[fos.index('fieldMinMax\n{'):]
    block = block[:block.index('}')]
    assert 'type            fieldMinMax;' in block and 'location        true;' in block
    assert 'fields          (U p T rho k omega nut alphat);' in block and 'writeControl    timeStep;' in block
    full['GLOBAL_SIM_CONTROL']['SOLVER'] = ['pimpleFoam']
    full['GLOBAL_SIM_CONTROL']['TURB_MODEL'] = ['kEpsilon']
    assert writeInternal.minMaxFields(full) == ['U', 'p', 'k', 'epsilon', 'nut']
    full['GLOBAL_SIM_CONTROL']['SOLVER'] = ['buoyantSimpleFoam']
    full['GLOBAL_SIM_CONTROL']['TURB_MODEL'] = ['laminar']
    assert writeInternal.minMaxFields(full) == ['U', 'p', 'p_rgh', 'T', 'rho']


def test_pimple_residual_control(tmp_path, monkeypatch, repo_root):
    #several outer correctors: the step ends once U, e and p are below OUTER_TOL / OUTER_TOL_P, so a large
    #N_OUTER_CORR only costs where a step needs it; one outer corrector (PISO): no residualControl
    templateLoc, full = _fullDefaults(repo_root, 'rhoPimpleFoam')
    monkeypatch.chdir(tmp_path)
    writeInternal.writeNumericsInternal(templateLoc, full)
    assert 'residualControl' not in (tmp_path / 'system/fvSolution').read_text()
    full['GLOBAL_SIM_CONTROL']['N_OUTER_CORR'] = ['8']
    writeInternal.writeNumericsInternal(templateLoc, full)
    text = (tmp_path / 'system/fvSolution').read_text()
    pimple = text[text.index('PIMPLE'):]
    assert 'nOuterCorrectors 8;' in pimple
    for line in ['p      { tolerance 0.0001; relTol 0; }', 'U      { tolerance 1e-05; relTol 0; }',
                 'e      { tolerance 1e-05; relTol 0; }']:
        assert line in pimple
    #incompressible: no energy; blank tolerances: no residualControl (every outer corrector runs)
    full['GLOBAL_SIM_CONTROL']['SOLVER'] = ['pimpleFoam']
    assert [l.split()[0] for l in writeInternal.outerResidualControl(full)[2:-1]] == ['p', 'U']
    full['GLOBAL_SIM_CONTROL']['OUTER_TOL_P'] = ['']
    full['GLOBAL_SIM_CONTROL']['OUTER_TOL'] = ['']
    assert writeInternal.outerResidualControl(full) == []
