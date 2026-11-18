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

"""Unit tests for setupCore (shared config/geometry plumbing)."""
import os

import pytest

import setupCore


ADDONS = ['POR', 'ROTA']


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(text)


@pytest.fixture
def templateLoc(tmp_path):
    loc = tmp_path / 'tpl'
    _write(str(loc / 'defaultSetup' / 'defaultOrder'),
           '[DEFAULT_ORDER]\nORDER = defaultDefaults defaultGeometry defaultGlobalControl\n')
    _write(str(loc / 'defaultSetup' / 'defaultDefaults'),
           '[DEFAULTS]\nDEFAULT_MODULES = defaultGlobalControl\nADDON_MODULES =\n')
    _write(str(loc / 'defaultSetup' / 'defaultGeometry'),
           '[GEOMETRY]\nGEOM = a.obj,1,5,1,1.1,high\n')
    _write(str(loc / 'defaultSetup' / 'defaultGlobalControl'),
           '[GLOBAL_CONTROL]\nENDTIME = 500\nWRITE = 1 2\n')
    _write(str(loc / 'defaultBCTemplates' / 'defaultPorous'),
           '[POR_SETUP]\nDCOEFFS = 1000\nVEC1 = 1 0 0\nSAMPLE_POR = True\n')
    return str(loc)


def test_version_is_single_string():
    assert setupCore.__version__.startswith('5.')


def test_getGlobalDefaults_reads_order_and_splits_values(templateLoc):
    d = setupCore.getGlobalDefaults(templateLoc)
    assert list(d['MODULE_LINK']) == ['defaultDefaults', 'defaultGeometry', 'defaultGlobalControl']
    assert d['MODULE_LINK']['defaultGlobalControl'] == ['GLOBAL_CONTROL']
    assert d['GLOBAL_CONTROL']['WRITE'] == ['1', '2']


def test_geomToDict_columns_and_duplicates():
    cols = ['scale', 'refinement', 'layers', 'expansion', 'wallmodel']
    g = setupCore.geomToDict({}, ['a.obj,1,5,1,1.1,high\nPOR-b.obj.gz,2,6,0,1.2,low'], cols)
    assert g['POR-b.obj.gz'] == {'scale': '2', 'refinement': '6', 'layers': '0',
                                 'expansion': '1.2', 'wallmodel': 'low'}
    with pytest.raises(SystemExit):
        setupCore.geomToDict({}, ['a.obj,1,5,1,1.1,high\na.obj,1,5,1,1.1,high'], cols)
    with pytest.raises(SystemExit):
        setupCore.geomToDict({}, ['a.obj,1,5'], cols)


def test_fillAddonSection_new_section_sets_flag():
    w = {}
    items = [('A', '1'), ('B', '2 3')]
    assert setupCore.fillAddonSection(w, 'POR-x', items) is True
    assert w['POR-x'] == {'A': ['1'], 'B': ['2 3']}


def test_fillAddonSection_keeps_user_values_fills_missing_drops_extra():
    w = {'POR-x': {'EXTRA': ['9'], 'B': ['7']}}
    flag = setupCore.fillAddonSection(w, 'POR-x', [('A', '1'), ('B', '2')])
    assert flag is True
    assert list(w['POR-x']) == ['A', 'B']          # template order, EXTRA dropped
    assert w['POR-x'] == {'A': ['1'], 'B': ['7']}


def test_fillAddonSection_complete_section_no_flag():
    w = {'POR-x': {'A': ['5'], 'B': ['6']}}
    assert setupCore.fillAddonSection(w, 'POR-x', [('A', '1'), ('B', '2')]) is False
    assert w['POR-x'] == {'A': ['5'], 'B': ['6']}


def test_checkAddonDefaults_only_touches_matching_prefix(templateLoc):
    geom = {'body.obj': {}, 'POR-rad.obj.gz': {}}
    w, flag = setupCore.checkAddonDefaults(geom, {}, templateLoc, 'POR', 'defaultPorous', 'POR_SETUP')
    assert flag is True
    assert list(w) == ['POR-rad']
    assert w['POR-rad']['VEC1'] == ['1 0 0']


def test_checkAddonDefaults_bad_section_header_exits(templateLoc):
    with pytest.raises(SystemExit):
        setupCore.checkAddonDefaults({'POR-a.obj': {}}, {}, templateLoc, 'POR', 'defaultPorous', 'NOPE')


def test_getCaseSetupDefaultModules_hidden_exposed_and_addons(templateLoc):
    defaults = setupCore.getGlobalDefaults(templateLoc)
    case = {'DEFAULTS': {'DEFAULT_MODULES': ['defaultGlobalControl'], 'ADDON_MODULES': ['']},
            'GEOMETRY': {'GEOM': ['a.obj,1,5,1,1.1,high']},
            'POR-x': {'A': ['1']}}
    w, full, flag = setupCore.getCaseSetupDefaultModules(case, defaults, ADDONS)
    assert 'GLOBAL_CONTROL' not in w                 # hidden default module
    assert full['GLOBAL_CONTROL']['ENDTIME'] == ['500']
    assert w['POR-x'] == full['POR-x'] == {'A': ['1']}
    assert flag is False

    # exposing a module with a missing key fills it and flags an update
    case['DEFAULTS']['DEFAULT_MODULES'] = ['']
    case['GLOBAL_CONTROL'] = {'ENDTIME': ['100']}
    w, full, flag = setupCore.getCaseSetupDefaultModules(case, defaults, ADDONS)
    assert flag is True
    assert w['GLOBAL_CONTROL'] == {'ENDTIME': ['100'], 'WRITE': ['1', '2']}


def test_getCaseSetupDefaultModules_invalid_module_exits(templateLoc):
    defaults = setupCore.getGlobalDefaults(templateLoc)
    case = {'DEFAULTS': {'DEFAULT_MODULES': ['notAModule']}}
    with pytest.raises(SystemExit):
        setupCore.getCaseSetupDefaultModules(case, defaults, ADDONS)


def test_cleanUpCaseSetup_removes_orphan_addons():
    geom = {'POR-a.obj': {}}
    w = {'POR-a': {}, 'POR-old': {}, 'TITLES': {}}
    full = {'POR-a': {}, 'POR-old': {}, 'TITLES': {}}
    setupCore.cleanUpCaseSetup(geom, w, full, ADDONS)
    assert set(w) == set(full) == {'POR-a', 'TITLES'}


def test_writeToCaseSetup_roundtrip_and_update_exit(tmp_path, templateLoc):
    d = {'DEFAULTS': {'DEFAULT_MODULES': ['defaultGlobalControl'], 'ADDON_MODULES': ['']},
         'GEOMETRY': {'GEOM': ['a.obj,1,5,1,1.1,high']}}
    setupCore.writeToCaseSetup(d, str(tmp_path), 'internalSetup')
    defaults = setupCore.getGlobalDefaults(templateLoc)
    case, w, full, flag = setupCore.getCaseSetup(defaults, str(tmp_path), ADDONS, 'internalSetup')
    assert case == d
    assert flag is False

    with pytest.raises(SystemExit):
        setupCore.writeToCaseSetup(d, str(tmp_path), 'internalSetup', updateFlag=True)
    assert os.path.isfile(str(tmp_path / 'internalSetup'))   # still written before exit


def test_writeNewCaseSetup_writes_only_starter_sections(tmp_path, templateLoc):
    defaults = setupCore.getGlobalDefaults(templateLoc)
    setupCore.writeNewCaseSetup(defaults, str(tmp_path))
    text = (tmp_path / 'caseSetup').read_text()
    assert '[DEFAULTS]' in text and '[GEOMETRY]' in text
    assert '[GLOBAL_CONTROL]' not in text
