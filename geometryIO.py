"""Centralized, hardened surface-geometry readers (OBJ / ASCII STL).

This module provides a single place to parse vertex data from Wavefront OBJ and
ASCII STL surface files, transparently handling gzip compression. It replaces
several divergent ad-hoc parsers that used fragile patterns such as
``line.startswith('v')`` (which also matches ``vn``/``vt``/``vp``) and
``line.split(' ')`` (which breaks on runs of multiple spaces).

Rules enforced here (matching the production ``load_obj`` reader in utilities):
  * OBJ vertices are lines beginning with ``v `` (a ``v`` followed by whitespace),
    which excludes vertex normals (``vn``), texture coords (``vt``) and parameter
    space vertices (``vp``).
  * Tokens are split on arbitrary whitespace via ``str.split()``.
  * ASCII STL vertices are lines whose first token is ``vertex``.
"""
import gzip

import numpy as np


def _open_text_maybe_gz(path, mode="rt"):
    """Open ``path`` for text I/O, transparently handling ``.gz`` files."""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _looks_like_stl(path):
    """Best-effort format sniff: True for ASCII STL, False otherwise (OBJ)."""
    lower = str(path).lower()
    # strip a trailing .gz to inspect the real extension
    if lower.endswith(".gz"):
        lower = lower[:-3]
    if lower.endswith(".stl"):
        return True
    if lower.endswith(".obj"):
        return False
    # fall back to content sniffing
    with _open_text_maybe_gz(path, "rt") as handle:
        for _ in range(5):
            line = handle.readline()
            if not line:
                break
            stripped = line.strip().lower()
            if not stripped:
                continue
            if stripped.startswith("solid"):
                return True
            if stripped.startswith(("v ", "vn", "vt", "vp", "f ", "o ", "g ", "#")):
                return False
    return False


def readObjVertices(path):
    """Return an ``(N, 3)`` float array of vertex coordinates from an OBJ file.

    Only ``v `` lines are read; ``vn``/``vt``/``vp`` are ignored. Tokens are split
    on arbitrary whitespace so irregular spacing is handled correctly.
    """
    vertices = []
    with _open_text_maybe_gz(path, "rt") as handle:
        for line in handle:
            if line.startswith("v "):
                vertices.append([float(token) for token in line.split()[1:4]])
    if not vertices:
        return np.empty((0, 3), dtype=np.float64)
    return np.asarray(vertices, dtype=np.float64)


def readStlVertices(path):
    """Return an ``(N, 3)`` float array of vertex coordinates from an ASCII STL.

    Only lines whose first token is ``vertex`` are read; tokens are split on
    arbitrary whitespace.
    """
    vertices = []
    with _open_text_maybe_gz(path, "rt") as handle:
        for line in handle:
            tokens = line.split()
            if tokens and tokens[0].lower() == "vertex":
                vertices.append([float(token) for token in tokens[1:4]])
    if not vertices:
        return np.empty((0, 3), dtype=np.float64)
    return np.asarray(vertices, dtype=np.float64)


def readSurfaceVertices(path):
    """Return an ``(N, 3)`` vertex array from an OBJ or ASCII STL file.

    The format is inferred from the file extension (``.obj`` / ``.stl``, with an
    optional ``.gz`` suffix), falling back to content sniffing when the extension
    is ambiguous.
    """
    if _looks_like_stl(path):
        return readStlVertices(path)
    return readObjVertices(path)


def surfaceBoundingBox(path):
    """Return the axis-aligned bounding box of a surface file.

    :returns: ``(minX, minY, minZ, maxX, maxY, maxZ)`` as floats.
    :raises ValueError: if the file contains no vertices.
    """
    vertices = readSurfaceVertices(path)
    if vertices.shape[0] == 0:
        raise ValueError("No vertices found in surface file: %s" % path)
    mins = vertices.min(axis=0)
    maxs = vertices.max(axis=0)
    return (
        float(mins[0]), float(mins[1]), float(mins[2]),
        float(maxs[0]), float(maxs[1]), float(maxs[2]),
    )
