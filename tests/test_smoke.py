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

"""Smoke tests: the core modules must import and compile.

These give the CI a green baseline before the math-specific unit tests land.
They intentionally avoid importing the postUtilities modules, whose optional
dependencies (scipy, python-pptx, gspread, ...) may not be installed.
"""
import importlib
import py_compile

import pytest

CORE_MODULES = [
    "utilities",
    "writeConstant",
    "writeScripts",
    "writeSystem",
]


@pytest.mark.parametrize("module_name", CORE_MODULES)
def test_core_module_imports(module_name):
    """Each core module imports without error."""
    assert importlib.import_module(module_name) is not None


def test_caseSetup_compiles(repo_root):
    """caseSetup.py byte-compiles (it runs main() under a __main__ guard)."""
    import os

    py_compile.compile(
        os.path.join(repo_root, "caseSetup.py"), doraise=True
    )
