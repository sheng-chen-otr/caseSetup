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


def readSurfacePIDs(path):
    """Return the PID (region) names of an OBJ or ASCII STL file, in file order.

    STL: the name after each ``solid`` keyword. OBJ: ``g`` group names, or ``o``
    object names when the file has no groups. These become the snappyHexMesh
    regions, i.e. patches ``<surfaceName>_<pid>``. Repeated names are listed once.
    """
    stl = _looks_like_stl(path)
    groups, objects = [], []
    with _open_text_maybe_gz(path, "rt") as handle:
        for line in handle:
            tokens = line.split()
            if not tokens:
                continue
            key = tokens[0].lower() if stl else tokens[0]
            name = " ".join(tokens[1:])
            if stl and key == "solid":
                groups.append(name)
            elif not stl and key == "g" and name:
                groups.append(name)
            elif not stl and key == "o" and name:
                objects.append(name)
    names = groups if groups or stl else objects
    return list(dict.fromkeys(names))


def _objIndex(token, nVertices):
    """0-based vertex index from an OBJ face token (``7``, ``7/1``, ``7//2``, ``-1``)."""
    idx = int(token.split("/")[0])
    return idx - 1 if idx > 0 else nVertices + idx


def readSurfaceTriangles(path):
    """Return ``{pid: (N, 3, 3) array}`` of triangles per PID, in file order.

    STL: facets grouped by ``solid`` name. OBJ: faces grouped by the current
    ``g`` (or ``o`` when there are no groups) name; polygons are fan
    triangulated. Faces before any group go to PID ``""``.
    """
    tris = {}
    if _looks_like_stl(path):
        current, facet = "", []
        with _open_text_maybe_gz(path, "rt") as handle:
            for line in handle:
                tokens = line.split()
                if not tokens:
                    continue
                key = tokens[0].lower()
                if key == "solid":
                    current = " ".join(tokens[1:])
                    tris.setdefault(current, [])
                elif key == "vertex":
                    facet.append([float(t) for t in tokens[1:4]])
                    if len(facet) == 3:
                        tris.setdefault(current, []).append(facet)
                        facet = []
    else:
        hasGroups = False
        with _open_text_maybe_gz(path, "rt") as handle:
            for line in handle:
                if line.startswith("g ") and line.split()[1:]:
                    hasGroups = True
                    break
        vertices, current = [], ""
        with _open_text_maybe_gz(path, "rt") as handle:
            for line in handle:
                tokens = line.split()
                if not tokens:
                    continue
                if tokens[0] == "v":
                    vertices.append([float(t) for t in tokens[1:4]])
                elif (tokens[0] == "g" and hasGroups) or (tokens[0] == "o" and not hasGroups):
                    current = " ".join(tokens[1:])
                    tris.setdefault(current, [])
                elif tokens[0] == "f":
                    idx = [_objIndex(t, len(vertices)) for t in tokens[1:]]
                    for k in range(1, len(idx) - 1):
                        tris.setdefault(current, []).append(
                            [vertices[idx[0]], vertices[idx[k]], vertices[idx[k + 1]]])
    return {pid: np.asarray(t, dtype=np.float64).reshape(-1, 3, 3) for pid, t in tris.items()}


def translateSurfacePIDs(path, outPath, shifts):
    """Write a copy of an OBJ / ASCII STL with the vertices of some PIDs translated.

    ``shifts`` maps PID name -> (dx, dy, dz) in file units. Output is plain text
    (not gzipped). For OBJ, a vertex shared by PIDs with different shifts is an
    error (the parts must be separate shells so they can move apart).
    """
    shifts = {pid: np.asarray(s, dtype=np.float64) for pid, s in shifts.items()}
    with _open_text_maybe_gz(path, "rt") as handle:
        lines = handle.readlines()
    out = []
    if _looks_like_stl(path):
        current = ""
        for line in lines:
            tokens = line.split()
            if tokens and tokens[0].lower() == "solid":
                current = " ".join(tokens[1:])
            if tokens and tokens[0].lower() == "vertex" and current in shifts:
                p = np.array([float(t) for t in tokens[1:4]]) + shifts[current]
                indent = line[:len(line) - len(line.lstrip())]
                line = "%svertex %.9g %.9g %.9g\n" % (indent, p[0], p[1], p[2])
            out.append(line)
    else:
        hasGroups = any(l.startswith("g ") and l.split()[1:] for l in lines)
        #which shift each vertex gets, from the faces that use it
        vertexShift, vertexCount, current = {}, 0, ""
        for line in lines:
            tokens = line.split()
            if not tokens:
                continue
            if tokens[0] == "v":
                vertexCount += 1
            elif (tokens[0] == "g" and hasGroups) or (tokens[0] == "o" and not hasGroups):
                current = " ".join(tokens[1:])
            elif tokens[0] == "f":
                for t in tokens[1:]:
                    vi = _objIndex(t, vertexCount)
                    prev = vertexShift.get(vi)
                    if prev is not None and prev != current and (
                            prev in shifts or current in shifts):
                        raise ValueError("Vertex %d is shared by PIDs %s and %s, which move "
                                         "differently: make them separate shells in %s"
                                         % (vi + 1, prev, current, path))
                    vertexShift[vi] = current
        vi = 0
        for line in lines:
            if line.startswith("v "):
                pid = vertexShift.get(vi)
                if pid in shifts:
                    p = np.array([float(t) for t in line.split()[1:4]]) + shifts[pid]
                    line = "v %.9g %.9g %.9g\n" % (p[0], p[1], p[2])
                vi += 1
            out.append(line)
    with open(outPath, "w") as handle:
        handle.writelines(out)


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
