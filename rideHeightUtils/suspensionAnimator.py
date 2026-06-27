#!/usr/bin/env python3
'''
suspensionAnimator.py

Standalone visualisation tool that animates the motion of the suspension
components across the ride-height map using the SAME kinematic pipeline as the
case generator (utilities.calculateRideHeights ->
utilities._compute_component_transforms). It loads the actual triSurface STL/OBJ
geometry, keeps only the suspension components (matched by the per-corner PID
keywords from the hardpoint CFG), and rigidly moves each component in-memory for
every ride-height point.

Output is one animation file per view (front/side/top), each written separately:
    - <base>_front.gif : looking in the +X direction
    - <base>_side.gif  : looking in the +Y direction
    - <base>_top.gif   : looking in the -Z direction

This is intended to be invoked by caseSetup.py only when the --animateSuspension
flag is supplied, but it can also be run standalone:

    python rideHeightUtils/suspensionAnimator.py -c <caseDir>

It reuses utilities for the kinematics so the animation matches exactly what the
ride-height case generator produces. Rendering uses PyVista (VTK) off-screen,
which is far faster than matplotlib's mplot3d. Install with:

    pip install pyvista imageio imageio-ffmpeg
'''

import os
import sys
import re
import gzip
import argparse
import time

import numpy as np

# PyVista (VTK) is the rendering backend: it updates mesh coordinates in C++ and
# renders off-screen, which is far faster than matplotlib's mplot3d (the latter
# re-sorts every polygon in Python on each frame).
try:
    import pyvista as pv
    pv.OFF_SCREEN = True
    _HAVE_PYVISTA = True
except Exception:  # pragma: no cover - import-time environment check
    pv = None
    _HAVE_PYVISTA = False


def _ensure_offscreen_display():
    '''
    Make sure VTK has a usable GL context on headless machines.

    On Linux without a DISPLAY, VTK's X backend segfaults. PyVista can spin up a
    virtual framebuffer (xvfb) to give it an off-screen target. Returns True if a
    display is available (or not needed), False if rendering would crash.
    '''
    if not _HAVE_PYVISTA:
        return False
    # macOS/Windows do not use X11; rendering works without a DISPLAY.
    if sys.platform != 'linux':
        return True
    if os.environ.get('DISPLAY'):
        return True
    # Headless Linux: try to start a virtual framebuffer.
    try:
        pv.start_xvfb()
        return bool(os.environ.get('DISPLAY'))
    except Exception as exc:
        print('\t\tWARNING! Could not start xvfb virtual framebuffer: %s' % exc)
        return False

# Make the repo root importable so we can reuse the kinematics in utilities.py.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
import utilities  # noqa: E402


# Colour per component type so the GIF is readable.
_COMPONENT_COLORS = {
    'WHEEL': '#222222',
    'UCA': '#d62728',
    'LCA': '#1f77b4',
    'ROCKER': '#2ca02c',
    'PUSHROD': '#9467bd',
    'DAMPER': '#8c564b',
    'TIE': '#ff7f0e',
}
_DEFAULT_COLOR = '#7f7f7f'

# Cap triangles per component so a heavy mesh does not make the render crawl.
# Faces are strided (decimated), not removed, so the silhouette is preserved.
_MAX_FACES_PER_COMPONENT = 40000


# --------------------------------------------------------------------------- #
# caseSetup loading (standalone use)
# --------------------------------------------------------------------------- #
def loadCaseSetup(caseDir):
    '''
    Read a case's ``caseSetup`` file into the same list-valued dict structure that
    caseSetup.getCaseSetup produces: ``{SECTION: {KEY: [tokens...]}}``.
    '''
    import configparser

    caseSetupPath = os.path.join(caseDir, 'caseSetup')
    if not os.path.isfile(caseSetupPath):
        raise FileNotFoundError('No caseSetup file found at: %s' % caseSetupPath)

    config = configparser.ConfigParser()
    config.optionxform = str
    with open(caseSetupPath) as fh:
        config.read_file(fh)

    fullCaseSetupDict = {}
    for section in config.sections():
        fullCaseSetupDict[section] = {}
        for key, value in config.items(section):
            fullCaseSetupDict[section][key] = value.split(' ')
    return fullCaseSetupDict


# --------------------------------------------------------------------------- #
# Geometry parsing (mirrors transformGeometryByPIDRegex conventions)
# --------------------------------------------------------------------------- #
def _open_maybe_gz(path, mode='rt'):
    if path.lower().endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def _parse_obj_groups(path):
    '''
    Parse an OBJ file into per-PID triangle groups.

    PID = current ``g``/``o`` group name (matching transformGeometryByPIDRegex).
    Returns dict: pid -> {'verts': (N,3) float array, 'faces': (M,3) int array}
    where faces index into the local per-PID vertex array.
    '''
    vertices = []
    currentPid = '__DEFAULT__'
    pidGlobalFaces = {currentPid: []}

    with _open_maybe_gz(path, 'rt') as fh:
        for line in fh:
            if line.startswith('v '):
                parts = line.split()
                vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith('f '):
                nVerts = len(vertices)
                refs = []
                for ref in line.split()[1:]:
                    if '/' in ref:
                        ref = ref.split('/', 1)[0]
                    idx = int(ref)
                    refs.append(idx - 1 if idx > 0 else nVerts + idx)
                # fan-triangulate any polygon into triangles
                for i in range(1, len(refs) - 1):
                    pidGlobalFaces[currentPid].append((refs[0], refs[i], refs[i + 1]))
            elif line.startswith('g ') or line.startswith('o '):
                toks = line.split(maxsplit=1)
                currentPid = toks[1].strip() if len(toks) > 1 else '__EMPTY__'
                if currentPid not in pidGlobalFaces:
                    pidGlobalFaces[currentPid] = []

    vertices = np.asarray(vertices, dtype=np.float64) if vertices else np.zeros((0, 3))
    return _localise_groups(vertices, pidGlobalFaces)


def _localise_groups(vertices, pidGlobalFaces):
    groups = {}
    for pid, faces in pidGlobalFaces.items():
        if len(faces) == 0:
            continue
        facesArr = np.asarray(faces, dtype=np.int64)
        used = np.unique(facesArr)
        remap = {g: l for l, g in enumerate(used)}
        localVerts = vertices[used]
        localFaces = np.vectorize(remap.get)(facesArr)
        groups[pid] = {'verts': localVerts, 'faces': localFaces.reshape(-1, 3)}
    return groups


def _parse_ascii_stl_groups(path):
    '''
    Parse an ASCII STL into per-PID triangle groups. PID = ``solid <name>``.
    Returns dict: pid -> {'verts': (N,3) array, 'faces': (M,3) int array}.
    '''
    groups = {}
    currentPid = None
    verts = []
    faces = []
    pending = []

    def _flush(pid):
        if pid is None or len(faces) == 0:
            return
        groups[pid] = {
            'verts': np.asarray(verts, dtype=np.float64),
            'faces': np.asarray(faces, dtype=np.int64).reshape(-1, 3),
        }

    with _open_maybe_gz(path, 'rt') as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith('solid'):
                currentPid = line[5:].strip() or '__EMPTY__'
                verts = []
                faces = []
                pending = []
            elif line.startswith('endsolid'):
                _flush(currentPid)
                currentPid = None
            elif line.startswith('vertex'):
                parts = line.split()
                pending.append((float(parts[1]), float(parts[2]), float(parts[3])))
                if len(pending) == 3:
                    base = len(verts)
                    verts.extend(pending)
                    faces.append((base, base + 1, base + 2))
                    pending = []
    # handle files with a single unterminated solid
    if currentPid is not None:
        _flush(currentPid)
    return groups


def _load_pid_groups(path):
    lower = path.lower()
    if '.obj' in lower:
        return _parse_obj_groups(path)
    if '.stl' in lower:
        if not utilities._is_ascii_stl_path(path):
            raise ValueError('Binary STL not supported for animation: %s' % path)
        return _parse_ascii_stl_groups(path)
    raise ValueError('Unsupported geometry format (need OBJ/STL): %s' % path)


# --------------------------------------------------------------------------- #
# PID -> component matching and rigid transform application
# --------------------------------------------------------------------------- #
def _match_pid_to_component(pidName, categorizedKeywords):
    '''
    Find the (corner, compType) that owns ``pidName`` by case-sensitive regex
    search of the escaped keywords (same matching basis as the generator).
    Returns (corner, compType) or None.
    '''
    for corner, compMap in categorizedKeywords.items():
        for compType, keywords in compMap.items():
            for kw in keywords:
                if re.search(re.escape(kw), pidName):
                    return corner, compType
    return None


def _apply_rigid_transform(verts, transform):
    '''
    Apply ``v' = R(rvec) . (v - pivot) + pivot + translation`` to an (N,3) array,
    matching utilities.transformGeometryByPIDRegex's rigid convention.
    '''
    rvec = np.asarray(transform.get('rotation_rvec', [0.0, 0.0, 0.0]), dtype=np.float64)
    translation = np.asarray(transform.get('translation', [0.0, 0.0, 0.0]), dtype=np.float64)
    pivot = transform.get('pivot', None)

    out = verts
    if float(np.linalg.norm(rvec)) > 1e-15:
        R = utilities._rodrigues_matrix_from_rvec(rvec)
        if pivot is not None:
            piv = np.asarray(pivot, dtype=np.float64)
            out = (verts - piv) @ R.T + piv
        else:
            out = verts @ R.T
    return out + translation


def _interpolate_row(rowA, rowB, alpha):
    '''
    Linearly blend the numeric columns of two ride-height rows.

    Every kinematic-solver column read by ``_compute_component_transforms``
    (``*_wc_x/y/z``, ``*_rvec_x/y/z``, ``*_rocker_deg``, ``*_damper_delta``,
    ``*_rack_travel`` ...) is numeric, so a per-column lerp produces a
    kinematically-consistent in-between pose. Non-numeric columns (caseName,
    steer mode, ...) are carried through from ``rowA`` unchanged.
    '''
    out = rowA.copy()
    for col in rowA.index:
        try:
            fa = float(rowA[col])
            fb = float(rowB[col])
        except (TypeError, ValueError):
            continue
        out[col] = fa + (fb - fa) * float(alpha)
    return out


# --------------------------------------------------------------------------- #
# Main animation builder
# --------------------------------------------------------------------------- #
def generateSuspensionAnimation(fullCaseSetupDict, caseDir=None, outputPath=None,
                                intervalMs=80, transitionFrames=12, holdFrames=6,
                                resolution=1600):
    '''
    Build the suspension-motion GIF for ``caseDir`` (defaults to cwd).

    The animation holds ``holdFrames`` frames at each ride-height point and inserts
    ``transitionFrames`` interpolated poses between consecutive points so the motion is
    smooth. ``intervalMs`` is the per-frame delay (GIF fps = 1000 / intervalMs).

    Returns the output GIF path, or None if it could not be generated.
    '''
    if caseDir is None:
        caseDir = os.getcwd()
    caseDir = os.path.abspath(caseDir)

    kinSetup = utilities._load_suspension_kinematics_setup(fullCaseSetupDict)
    if kinSetup is None:
        print('\tSuspension animation skipped: kinematic solver is not enabled '
              '(set [RIDE_HEIGHT_SETUP] USE_KINEMATIC_SOLVER true and HARDPOINT_FILE).')
        return None

    # Per-corner keyword lists -> per-corner/per-component keyword map.
    cornerPidKeywords = {'fl': [], 'fr': [], 'rl': [], 'rr': []}
    for corner in cornerPidKeywords:
        cornerPidKeywords[corner].extend(list(kinSetup['corner_pid_keywords'].get(corner, [])))
        # de-duplicate while preserving order
        seen = set()
        dedup = []
        for kw in cornerPidKeywords[corner]:
            if kw not in seen:
                dedup.append(kw)
                seen.add(kw)
        cornerPidKeywords[corner] = dedup
    categorizedKeywords = utilities._categorize_pid_keywords_by_component(cornerPidKeywords)

    # Ride-height table with kinematic solver columns (same call as the generator).
    print('\tSolving suspension kinematics across the ride-height map...')
    rideHeights = utilities.calculateRideHeights(fullCaseSetupDict)

    # Filter to RUN_RH_POINTS, preserving the table order.
    runPoints = fullCaseSetupDict['RIDE_HEIGHT_SETUP'].get('RUN_RH_POINTS', [''])
    if len(runPoints) == 1 and runPoints[0] == '':
        selected = rideHeights
    else:
        mask = rideHeights['point'].apply(lambda p: str(int(p)) in runPoints)
        selected = rideHeights[mask]
    selected = selected.reset_index(drop=True)
    if len(selected) == 0:
        print('\tSuspension animation skipped: no ride-height points selected.')
        return None

    # Load suspension geometry from the case triSurface dir, keep only matched PIDs.
    triSurfaceDir = os.path.join(caseDir, 'constant', 'triSurface')
    if not os.path.isdir(triSurfaceDir):
        print('\tSuspension animation skipped: no geometry directory at %s' % triSurfaceDir)
        return None

    components = []  # list of dicts: {pid, corner, comp, baseVerts, faces}
    geomFiles = [f for f in sorted(os.listdir(triSurfaceDir))
                 if ('.obj' in f.lower() or '.stl' in f.lower())]
    print('\tParsing %d geometry file(s) from %s ...' % (len(geomFiles), triSurfaceDir))
    for fIdx, fname in enumerate(geomFiles):
        print('\t\t[%d/%d] %s' % (fIdx + 1, len(geomFiles), fname))
        fpath = os.path.join(triSurfaceDir, fname)
        try:
            groups = _load_pid_groups(fpath)
        except Exception as exc:
            print('\t\tSkipping %s for animation: %s' % (fname, exc))
            continue
        for pid, group in groups.items():
            match = _match_pid_to_component(pid, categorizedKeywords)
            if match is None:
                continue
            corner, comp = match
            faces = group['faces']
            if len(faces) > _MAX_FACES_PER_COMPONENT:
                stride = int(np.ceil(len(faces) / _MAX_FACES_PER_COMPONENT))
                faces = faces[::stride]
            components.append({
                'pid': pid,
                'corner': corner,
                'comp': comp,
                'baseVerts': group['verts'],
                'faces': faces,
            })

    if len(components) == 0:
        print('\tSuspension animation skipped: no geometry PIDs matched the '
              'suspension keywords from the hardpoint CFG.')
        return None

    print('\tBuilding suspension animation over %d ride-height point(s), %d component(s)...'
          % (len(selected), len(components)))

    # Build the playback schedule: hold a few frames at each ride-height point so the
    # values are readable, then linearly interpolate the kinematic-solver columns into
    # ``transitionFrames`` in-between poses so the motion to the next point is smooth.
    # Each schedule entry is (row, srcPoint, dstPoint, alpha).
    transitionFrames = max(1, int(transitionFrames))
    holdFrames = max(1, int(holdFrames))
    schedule = []
    nPts = len(selected)
    for i in range(nPts):
        rowI = selected.iloc[i]
        ptI = int(rowI['point'])
        for _ in range(holdFrames):
            schedule.append((rowI, ptI, ptI, 0.0))
        if i < nPts - 1:
            rowNext = selected.iloc[i + 1]
            ptNext = int(rowNext['point'])
            for s in range(1, transitionFrames):
                alpha = s / float(transitionFrames)
                schedule.append((_interpolate_row(rowI, rowNext, alpha), ptI, ptNext, alpha))

    # Static bounding box (padded) for a fixed, equal-aspect camera across all frames.
    allVerts = np.vstack([c['baseVerts'] for c in components])
    mins = allVerts.min(axis=0)
    maxs = allVerts.max(axis=0)
    center = (mins + maxs) / 2.0
    span = float((maxs - mins).max()) * 1.25
    if span <= 0:
        span = 1.0
    half = span / 2.0

    rhUnit = fullCaseSetupDict['RIDE_HEIGHT_SETUP'].get('RH_UNIT', ['mm'])[0].lower()

    def _ride_height_label(row, srcPoint, dstPoint, alpha):
        scale = 1000.0 if rhUnit == 'm' else 1.0  # table fl/fr/rl/rr are in metres
        unit = 'mm'
        parts = []
        for corner in ('fl', 'fr', 'rl', 'rr'):
            if corner in row.index:
                parts.append('%s=%.1f' % (corner.upper(), float(row[corner]) * scale))
        rh = '  '.join(parts)
        extra = []
        for corner in ('fl', 'fr', 'rl', 'rr'):
            camCol = '%s_camber_deg' % corner
            if camCol in row.index:
                extra.append('%s cam=%.2f' % (corner.upper(), float(row[camCol])))
        if srcPoint == dstPoint or alpha <= 0.0:
            header = 'Point %d' % srcPoint
        else:
            header = 'Point %d -> %d  (%.0f%%)' % (srcPoint, dstPoint, alpha * 100.0)
        line = '%s   RH(%s): %s' % (header, unit, rh)
        if extra:
            line += '\n' + '  '.join(extra)
        return line

    if not _HAVE_PYVISTA:
        print('\tSuspension animation skipped: PyVista is not installed. Install it with\n'
              '\t\tpip install pyvista imageio imageio-ffmpeg')
        return None

    if not _ensure_offscreen_display():
        print('\tSuspension animation skipped: no X display and xvfb is unavailable.\n'
              '\t\tInstall xvfb (e.g. `sudo apt-get install xvfb`) or run on a node with a\n'
              '\t\tdisplay so the off-screen GL context can be created.')
        return None

    # Pre-build one VTK PolyData per component. VTK face format is a flat array of
    # [n, i0, i1, i2, n, ...]; here every face is a triangle (n == 3).
    meshes = []
    for comp in components:
        faces = np.asarray(comp['faces'], dtype=np.int64)
        nFaces = len(faces)
        faceArr = np.empty((nFaces, 4), dtype=np.int64)
        faceArr[:, 0] = 3
        faceArr[:, 1:] = faces
        poly = pv.PolyData(np.asarray(comp['baseVerts'], dtype=np.float64),
                           faceArr.ravel())
        meshes.append(poly)

    # Three orthographic views, each written to its OWN file:
    #   front : camera on -X looking +X   (Y-Z plane), Z up
    #   side  : camera on -Y looking +Y   (X-Z plane), Z up
    #   top   : camera on +Z looking -Z   (X-Y plane), Y up
    dist = span * 2.0
    views = [
        ('front', 'Front view (looking +X)',
         (center + np.array([-dist, 0.0, 0.0])), (0.0, 0.0, 1.0)),
        ('side', 'Side view (looking +Y)',
         (center + np.array([0.0, -dist, 0.0])), (0.0, 0.0, 1.0)),
        ('top', 'Top view (looking -Z)',
         (center + np.array([0.0, 0.0, dist])), (0.0, 1.0, 0.0)),
    ]

    # Resolve per-view output paths from the (optional) base output path.
    if outputPath is None:
        base, ext = os.path.join(caseDir, 'suspensionAnimation'), '.mp4'
    else:
        base, ext = os.path.splitext(outputPath)
        if ext.lower() not in ('.gif', '.mp4'):
            ext = '.mp4'
    fps = max(1.0, 1000.0 / float(intervalMs))
    parallelScale = half * 1.1
    # ffmpeg/h264 needs even (ideally /16) dimensions; round up to avoid a quality-
    # degrading auto-resize on mp4 export.
    resolution = int(np.ceil(resolution / 16.0) * 16)

    totalFrames = len(schedule)
    print('\tRendering %d frame(s) x %d view(s) (%d ride-height point(s), %d component(s))...'
          % (totalFrames, len(views), len(selected), len(components)))

    outputs = []
    for viewKey, viewTitle, camPos, viewUp in views:
        outPath = '%s_%s%s' % (base, viewKey, ext)
        plotter = pv.Plotter(off_screen=True, window_size=(resolution, resolution))
        plotter.set_background('white')
        plotter.enable_parallel_projection()
        # Smooth anti-aliased edges for higher visual quality.
        try:
            plotter.enable_anti_aliasing('ssaa')
        except Exception:
            pass
        for comp, poly in zip(components, meshes):
            color = _COMPONENT_COLORS.get(comp['comp'], _DEFAULT_COLOR)
            # Fully opaque solid surfaces (no transparency).
            plotter.add_mesh(poly, color=color, smooth_shading=True,
                             opacity=1.0, ambient=0.25, diffuse=0.7,
                             specular=0.3, specular_power=15,
                             reset_camera=False)
        plotter.camera.position = tuple(camPos)
        plotter.camera.focal_point = tuple(center)
        plotter.camera.up = viewUp
        plotter.camera.parallel_scale = parallelScale
        plotter.add_text(viewTitle, position='upper_edge', font_size=12,
                         color='black', name='title')

        if ext.lower() == '.mp4':
            # quality 1-10 (10 = best); high bitrate for crisp output.
            plotter.open_movie(outPath, framerate=int(round(fps)), quality=9)
        else:
            plotter.open_gif(outPath, fps=fps)

        renderStart = time.time()
        for frameIdx in range(totalFrames):
            row, srcPoint, dstPoint, alpha = schedule[frameIdx]
            cornerTransforms = {}
            for corner in ('fr', 'fl', 'rr', 'rl'):
                if corner in kinSetup['corners']:
                    cornerTransforms[corner] = utilities._compute_component_transforms(
                        kinSetup, corner, row, selected)
            for cIdx, comp in enumerate(components):
                transform = cornerTransforms.get(comp['corner'], {}).get(comp['comp'])
                if transform is None:
                    meshes[cIdx].points = np.asarray(comp['baseVerts'], dtype=np.float64)
                else:
                    meshes[cIdx].points = _apply_rigid_transform(comp['baseVerts'], transform)
            plotter.add_text(_ride_height_label(row, srcPoint, dstPoint, alpha),
                             position='lower_edge', font_size=11, color='black',
                             name='info')
            plotter.write_frame()

            done = frameIdx + 1
            elapsed = time.time() - renderStart
            pct = 100.0 * done / totalFrames
            eta = elapsed * (100.0 - pct) / pct if pct > 0 else 0.0
            sys.stdout.write('\r\t\t[%s] frame %d/%d (%5.1f%%) elapsed %4.0fs ETA %4.0fs'
                             % (viewKey, done, totalFrames, pct, elapsed, eta))
            sys.stdout.flush()
        sys.stdout.write('\n')
        plotter.close()
        print('\t\tSaved %s view (%.0fs): %s' % (viewKey, time.time() - renderStart, outPath))
        outputs.append(outPath)

    return outputs


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main():
    parser = argparse.ArgumentParser(
        prog='suspensionAnimator',
        description='Animate suspension component motion across the ride-height map '
                    'into a front/side/top GIF.')
    parser.add_argument('-c', '--case', default='.',
                        help='Case directory containing caseSetup and constant/triSurface.')
    parser.add_argument('-o', '--output', default=None,
                        help='Output GIF path (default: <case>/suspensionAnimation.gif).')
    parser.add_argument('--interval', type=int, default=80,
                        help='Per-frame interval in milliseconds (default: 80, ~12.5 fps).')
    parser.add_argument('--transition-frames', type=int, default=12, dest='transitionFrames',
                        help='Interpolated in-between frames per ride-height transition '
                             '(default: 12). Higher = smoother, longer GIF.')
    parser.add_argument('--hold-frames', type=int, default=6, dest='holdFrames',
                        help='Frames to pause at each ride-height point (default: 6).')
    parser.add_argument('--resolution', type=int, default=1600,
                        help='Square pixel resolution per view (default: 1600 -> 1600x1600).')
    args = parser.parse_args()

    caseDir = os.path.abspath(args.case)
    fullCaseSetupDict = loadCaseSetup(caseDir)
    out = generateSuspensionAnimation(fullCaseSetupDict, caseDir=caseDir,
                                      outputPath=args.output, intervalMs=args.interval,
                                      transitionFrames=args.transitionFrames,
                                      holdFrames=args.holdFrames, resolution=args.resolution)
    if out is None:
        sys.exit(1)


if __name__ == '__main__':
    main()
