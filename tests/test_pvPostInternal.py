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

#pvPostInternal helpers (no ParaView needed) and imagesToMovies (ffmpeg mocked)

import os
import sys
import configparser

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'postUtilities'))
import pvPostInternal
import imagesToMovies


def _setup(text):
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read_string(text)
    return config


def test_default_setup_and_case_override(tmp_path):
    setup = pvPostInternal.getPVSetup(str(tmp_path))
    assert setup['SLICE']['NAMES'].split() == ['yMid', 'xMid', 'zMid'] and setup['MOVIES']['FPS'] == '10'
    (tmp_path / 'pvPostInternalSetup').write_text('[SLICE]\nVARIABLES = Umag\n\n[MOVIES]\nFPS = 25\n')
    setup = pvPostInternal.getPVSetup(str(tmp_path))
    #case values win, keys the case file leaves out come from the default
    assert setup['SLICE']['VARIABLES'] == 'Umag' and setup['SLICE']['NORMALS'] == 'Y X Z'
    assert setup['MOVIES']['FPS'] == '25' and setup['PV_POST_MAIN']['RESOLUTION'] == '1920 1080'


def test_variable_dict_pref():
    setup = _setup('[PV_POST_MAIN]\nVARIABLE_DICT = default\n')
    varDict = pvPostInternal.getVariableDict(setup, 1.5e6)
    assert varDict['sliceVariables']['pGauge']['equation'] == 'p-1500000'
    assert varDict['surfaceVariables']['wallHeatFlux']['range'] == 'symmetric'


def test_slices_parse_and_errors():
    slices = pvPostInternal.getSlices({'NAMES': 'a b c', 'NORMALS': 'y Z x', 'ORIGINS': 'center loc 0,0,0.01',
                                       'FOLLOW': 'none displacer none'})
    assert [s['normal'] for s in slices] == ['Y', 'Z', 'X']
    assert slices[2]['origin'] == [0.0, 0.0, 0.01] and slices[1]['follow'] == 'displacer' and slices[0]['follow'] is None
    bounds = [-1, 1, -2, 2, 0, 0.5]
    assert pvPostInternal.sliceOrigin('center', bounds, [9, 9, 9]) == [0, 0, 0.25]
    assert pvPostInternal.sliceOrigin('loc', bounds, [9, 8, 7]) == [9, 8, 7]
    for bad in [{'NAMES': 'a b', 'NORMALS': 'Y', 'ORIGINS': 'center'},
                {'NAMES': 'a', 'NORMALS': 'W', 'ORIGINS': 'center'},
                {'NAMES': 'a', 'NORMALS': 'Y', 'ORIGINS': '1,2'}]:
        with pytest.raises(SystemExit):
            pvPostInternal.getSlices(bad)


def test_select_times_keeps_the_last():
    times = [0.0] + [0.005*i for i in range(1, 61)]
    assert pvPostInternal.selectTimes(times, 10.0)[0] == 0.005
    last = pvPostInternal.selectTimes(times, 10.0, lastCycles='1')
    assert len(last) == 21 and last[0] == pytest.approx(0.2)
    every = pvPostInternal.selectTimes(times, 10.0, every=5, lastCycles='1')
    assert every[-1] == pytest.approx(0.3) and len(every) == 5


def test_needed_fields():
    fields = ['U', 'p', 'T', 'rho']
    assert pvPostInternal.neededFields('mag(U)', fields) == ['U']
    assert pvPostInternal.neededFields('U_Z', fields) == ['U']
    assert pvPostInternal.neededFields('p-100000', fields) == ['p']
    assert pvPostInternal.neededFields('wallHeatFlux', fields) == []


def test_slice_camera_frames_the_cut():
    bounds = [-0.2, 0.2, -0.2, 0.2, -0.03, 0.03]
    cam = pvPostInternal.sliceCamera('Y', [0, 0.05, 0], bounds, [1920, 1080])
    #looks along y at the plane, z up, wide box -> width sets the scale
    assert cam['up'] == [0, 0, 1] and cam['focal'][1] == 0.05 and cam['position'][1] < 0.05
    assert cam['scale'] == pytest.approx(0.5*0.4*1080/1920*1.05)
    cam = pvPostInternal.sliceCamera('Z', [0, 0, 0.01], bounds, [1920, 1080])
    assert cam['up'] == [0, 1, 0] and cam['focal'][2] == 0.01
    with pytest.raises(SystemExit):
        pvPostInternal.surfaceCamera('Diagonal', bounds, [1920, 1080])


def test_motion_follow_prescribed_and_crank(tmp_path):
    case = tmp_path / 'case'
    tables = case / 'constant' / 'tables'
    tables.mkdir(parents=True)
    (tables / 'dispTable').write_text('FoamFile { }\n3\n(\n    (0 (0 0 0))\n    (0.05 (0 0 0.01))\n    (0.1 (0 0 0))\n)\n')
    config = _setup('[GLOBAL_CONTROL]\nFREQ = 10\n[MESH_MOTION]\nDRIVE = prescribed\n[MOVE-disp]\nFROM = self\n'
                    '[SETUP_INFO]\nSTART_ANGLE = 90\n')
    motion = pvPostInternal.Motion(str(case), config)
    assert motion.dz('disp', 0.025) == pytest.approx(0.005)
    assert motion.angle(0.025) == pytest.approx(180.0)
    with pytest.raises(SystemExit):
        motion.dz('other', 0.0)
    #crank drive: theta(t) from the crank log, dz from the crank angle table
    (tables / 'dispCrank').write_text('3\n(\n    (0 0)\n    (180 0.02)\n    (360 0)\n)\n')
    post = case / 'postProcessing' / 'crank' / '0'
    post.mkdir(parents=True)
    (post / 'crank.dat').write_text('# Time theta omega torqueGas torqueLoad powerElec\n'
                                    '0.01 36 62.8 0 0 0\n0.05 180 62.8 0 0 0\n0.1 360 62.8 0 0 0\n')
    config['MESH_MOTION']['DRIVE'] = 'crank'
    motion = pvPostInternal.Motion(str(case), config)
    assert motion.angle(0.05) == pytest.approx(180.0)
    assert motion.dz('disp', 0.05) == pytest.approx(0.02)


def test_movies_frames_and_command(tmp_path, monkeypatch):
    images = tmp_path / 'postProcessing' / 'images'
    for panel in ['yMid_Umag', 'zMid_p']:
        (images / panel).mkdir(parents=True)
        for i in [2, 0, 1]:
            (images / panel / ('case_%s_%05d.png' % (panel, i))).write_text('x')
    (images / 'empty').mkdir()
    frames = imagesToMovies.findFrames(str(images))
    assert list(frames.keys()) == ['yMid_Umag', 'zMid_p']
    assert [os.path.basename(f)[-9:] for f in frames['zMid_p']] == ['00000.png', '00001.png', '00002.png']
    command = imagesToMovies.ffmpegCommand('list.txt', 12, 'out.mp4')
    assert 'setpts=N/(12*TB),scale=trunc(iw/2)*2:trunc(ih/2)*2' in command and command[-1] == 'out.mp4'

    #ffmpeg mocked: one call per panel at the FPS from the case setup, frame lists removed afterwards
    calls = []

    class Result:
        returncode = 0
        stderr = ''

    def fakeRun(command, **kwargs):
        listFile = command[command.index('-i') + 1]
        calls.append((command, open(listFile).read()))
        return Result()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(imagesToMovies.shutil, 'which', lambda name: '/usr/bin/ffmpeg')
    monkeypatch.setattr(imagesToMovies.subprocess, 'run', fakeRun)
    (tmp_path / 'pvPostInternalSetup').write_text('[MOVIES]\nENABLE = True\nFPS = 24\n')
    assert imagesToMovies.main([]) == 0
    assert len(calls) == 2 and '-r' in calls[0][0] and calls[0][0][calls[0][0].index('-r') + 1] == '24'
    assert calls[0][1].count("file '") == 3
    assert not list((tmp_path / 'postProcessing' / 'movies').glob('*.txt'))
    #ENABLE False: no videos
    calls.clear()
    (tmp_path / 'pvPostInternalSetup').write_text('[MOVIES]\nENABLE = False\n')
    assert imagesToMovies.main([]) == 0 and calls == []
