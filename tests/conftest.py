"""Shared pytest fixtures for the caseSetup test suite.

These fixtures provide small in-memory geometry samples and paths so unit
tests for the pure-Python math (kinematics, cornering, geometry I/O, radius
calculation) can run without an OpenFOAM installation.
"""
import os
import sys
import textwrap

import pytest

# Make the repository root importable so tests can `import utilities`, etc.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)


@pytest.fixture
def repo_root():
    """Absolute path to the repository root."""
    return REPO_ROOT


@pytest.fixture
def sample_obj_text():
    """A minimal ASCII OBJ with mixed v / vn / vt lines and multiple spaces.

    Deliberately includes the parsing traps that have bitten this codebase:
    - `vn` / `vt` lines that must NOT be read as vertices
    - multiple/irregular spaces between tokens
    """
    return textwrap.dedent(
        """\
        # sample object
        v 0.0 0.0 0.0
        v  1.0   0.0 0.0
        vn 0.0 0.0 1.0
        vt 0.5 0.5
        v 1.0 1.0 0.0
        v 0.0 1.0 2.0
        f 1//1 2//1 3//1
        """
    )


@pytest.fixture
def sample_obj_file(tmp_path, sample_obj_text):
    """Write the sample OBJ to a temp file and return its path."""
    p = tmp_path / "sample.obj"
    p.write_text(sample_obj_text)
    return str(p)


@pytest.fixture
def sample_obj_vertices():
    """Expected vertex coordinates for `sample_obj_text` (v lines only)."""
    return [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 2.0),
    ]


@pytest.fixture
def sample_stl_text():
    """A minimal ASCII STL with irregular leading whitespace on vertex lines."""
    return textwrap.dedent(
        """\
        solid sample
          facet normal 0 0 1
            outer loop
              vertex 0.0 0.0 0.0
              vertex   1.0   0.0 0.0
              vertex 0.0 1.0 2.0
            endloop
          endfacet
        endsolid sample
        """
    )


@pytest.fixture
def sample_stl_file(tmp_path, sample_stl_text):
    """Write the sample STL to a temp file and return its path."""
    p = tmp_path / "sample.stl"
    p.write_text(sample_stl_text)
    return str(p)


@pytest.fixture
def sample_stl_vertices():
    """Expected vertex coordinates for `sample_stl_text`."""
    return [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 2.0),
    ]

