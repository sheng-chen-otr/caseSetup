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

#makes images of an internal case from the foamToEnsight output, one image per written time per panel
#panels are slices x variables, surface views x variables, volume renderings (views x variables) and clipped
#half volumes (clips x variables), set up like pvPost with pvPostInternalSetup
#colour ranges are kept fixed over the run, slices can follow a moving body with FOLLOW
#images go to postProcessing/images/<panel>, imagesToMovies.py makes the videos
#usage (in the case directory): pvbatch --force-offscreen-rendering pvPostInternal.py [--every N] [--lastCycles N]
#  [--panels slice,surface,volume,clip]

import os
import sys
import re
import glob
import json
import time
import argparse
import configparser

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import postInternal

templateLoc = os.path.dirname(os.path.realpath(__file__))
normalVectors = {'X': [1, 0, 0], 'Y': [0, 1, 0], 'Z': [0, 0, 1]}


def parseArgs(argv=None):
    parser = argparse.ArgumentParser(prog='pvPostInternal', description='Images of an internalSetup case (pvbatch)')
    parser.add_argument('--every', type=int, default=None, help='use every N-th written time (overrides EVERY)')
    parser.add_argument('--lastCycles', default=None, help='only the last N cycles, or all (overrides LAST_CYCLES)')
    parser.add_argument('--mesh', action='store_true',
                        help='mesh images of the motion dry run (EnSightMotion -> postProcessing/meshImages)')
    parser.add_argument('--panels', default=None,
                        help='panel types to make, e.g. volume,clip (overrides SLICE_IMG SURFACE_IMG VOLUME_IMG CLIP_IMG)')
    return parser.parse_args(argv)


panelTypes = ('slice', 'surface', 'volume', 'clip')


def getPanelTypes(mainSetup, panels=None):
    #panel types switched on in PV_POST_MAIN (<TYPE>_IMG), or the --panels list
    if panels:
        chosen = [p.strip().lower() for p in panels.split(',') if p.strip()]
        for p in chosen:
            if p not in panelTypes:
                sys.exit('ERROR! --panels %s unknown, use %s' % (p, ','.join(panelTypes)))
        return chosen
    return [p for p in panelTypes if mainSetup.get('%s_IMG' % (p.upper()), 'False').lower() == 'true']


#----- settings

def readConfig(path):
    config = configparser.ConfigParser()
    config.optionxform = str
    with open(path) as f:
        config.read_file(f)
    return config


def getPVSetup(caseDir):
    #getting pvPostInternalSetup from the case, missing keys from the default file
    defaultPath = os.path.join(templateLoc, 'default', 'pvPostInternalSetup')
    setup = readConfig(defaultPath)
    casePath = os.path.join(caseDir, 'pvPostInternalSetup')
    if os.path.isfile(casePath):
        print('\t\tFound pvPostInternalSetup in the case directory')
        caseSetup = readConfig(casePath)
        for section in caseSetup.sections():
            if not setup.has_section(section):
                setup.add_section(section)
            for key, value in caseSetup[section].items():
                setup[section][key] = value
    else:
        print('\t\tNo pvPostInternalSetup in the case directory, using the default settings')
    return setup


def getVariableDict(setup, pRef):
    #reading the variable dict, PREF is replaced by the charge pressure
    path = setup['PV_POST_MAIN']['VARIABLE_DICT']
    if path.lower() == 'default':
        path = os.path.join(templateLoc, 'default', 'defaultInternalVariables')
    with open(path) as f:
        text = f.read().replace('PREF', '%.9g' % (pRef))
    return json.loads(text)


def getSlices(sliceSetup):
    #reading the slice lists (NAMES NORMALS ORIGINS FOLLOW)
    names = sliceSetup['NAMES'].split()
    normals = [n.upper() for n in sliceSetup['NORMALS'].split()]
    origins = sliceSetup['ORIGINS'].split()
    follow = sliceSetup.get('FOLLOW', ' '.join(['none']*len(names))).split()
    if not len(names) == len(normals) == len(origins) == len(follow):
        sys.exit('ERROR! pvPostInternalSetup [SLICE]: NAMES, NORMALS, ORIGINS and FOLLOW need the same number of entries!')
    slices = []
    for name, normal, origin, body in zip(names, normals, origins, follow):
        if normal not in normalVectors:
            sys.exit('ERROR! pvPostInternalSetup [SLICE] %s: NORMAL must be X, Y or Z!' % (name))
        if origin.lower() not in ('center', 'loc'):
            try:
                origin = [float(v) for v in origin.split(',')]
            except ValueError:
                origin = []
            if len(origin) != 3:
                sys.exit('ERROR! pvPostInternalSetup [SLICE] %s: ORIGIN must be center, loc or x,y,z!' % (name))
        slices.append({'name': name, 'normal': normal, 'origin': origin,
                       'follow': None if body.lower() == 'none' else body})
    return slices


def getClips(clipSetup):
    #reading the clip lists (NAMES NORMALS ORIGINS), the normal points to the removed half (the camera side)
    names = clipSetup['NAMES'].split()
    normals = [n.upper() for n in clipSetup['NORMALS'].split()]
    origins = clipSetup['ORIGINS'].split()
    if not len(names) == len(normals) == len(origins):
        sys.exit('ERROR! pvPostInternalSetup [CLIP]: NAMES, NORMALS and ORIGINS need the same number of entries!')
    clips = []
    for name, normal, origin in zip(names, normals, origins):
        sign = -1.0 if normal.startswith('-') else 1.0
        axis = normal.lstrip('+-')
        if axis not in normalVectors:
            sys.exit('ERROR! pvPostInternalSetup [CLIP] %s: NORMAL must be X, Y, Z, -X, -Y or -Z!' % (name))
        if origin.lower() not in ('center', 'loc'):
            try:
                origin = [float(v) for v in origin.split(',')]
            except ValueError:
                origin = []
            if len(origin) != 3:
                sys.exit('ERROR! pvPostInternalSetup [CLIP] %s: ORIGIN must be center, loc or x,y,z!' % (name))
        clips.append({'name': name, 'normal': [sign*v for v in normalVectors[axis]], 'origin': origin})
    return clips


def sliceOrigin(origin, bounds, locInMesh):
    #slice origin from the mesh box center, LOC_IN_MESH or the given point
    if isinstance(origin, list):
        return list(origin)
    if origin.lower() == 'loc':
        return list(locInMesh)
    return [0.5*(bounds[0] + bounds[1]), 0.5*(bounds[2] + bounds[3]), 0.5*(bounds[4] + bounds[5])]


def selectTimes(times, freq, every=1, lastCycles='all'):
    #picking the written times to use
    times = [t for t in times if t > 0]
    if str(lastCycles).lower() != 'all' and times:
        tStart = times[-1] - float(lastCycles)/freq - 1e-9
        times = [t for t in times if t >= tStart]
    #counting back from the last time so the last time is always used
    return times[::-1][::max(int(every), 1)][::-1]


def neededFields(equation, fields):
    #fields of the case used in the equation
    names = set(re.findall(r'[A-Za-z_][A-Za-z0-9_]*', equation))
    names |= set(name.split('_')[0] for name in names)
    return [field for field in fields if field in names]


#----- motion

class Motion:
    #z displacement of each body and the crank angle vs time
    def __init__(self, caseDir, config):
        self.bodies = [s[len('MOVE-'):] for s in config.sections() if s.startswith('MOVE-')]
        self.freq = float(config['GLOBAL_CONTROL']['FREQ'])
        self.startAngle = None
        if config.has_section('SETUP_INFO') and config['SETUP_INFO'].get('START_ANGLE', '') != '':
            self.startAngle = float(config['SETUP_INFO']['START_ANGLE'])
        self.crank = None
        self.tables = {}
        drive = config['MESH_MOTION'].get('DRIVE', 'prescribed') if config.has_section('MESH_MOTION') else 'prescribed'
        #crank tables are against the crank angle: without a crank log (motion dry run) the crank turns at FREQ
        self.crankTables = bool(self.bodies) and drive == 'crank'
        tableDir = os.path.join(caseDir, 'constant', 'tables')
        if self.bodies and drive == 'crank':
            names, crank, theta0 = postInternal.readCrank(caseDir, os.path.join(caseDir, 'postProcessing'))
            if crank is not None:
                self.crank = (crank[:, 0], crank[:, 1])
            for body in self.bodies:
                path = os.path.join(tableDir, '%sCrank' % (body))
                if os.path.isfile(path):
                    self.tables[body] = readCrankTable(path)
        else:
            for body in self.bodies:
                path = os.path.join(tableDir, '%sTable' % (body))
                if os.path.isfile(path):
                    self.tables[body] = postInternal.readTable(path)

    def angle(self, t):
        #crank angle at time t (deg), None if unknown
        if self.crank is not None:
            return float(np.mod(np.interp(t, self.crank[0], self.crank[1]), 360.0))
        if self.startAngle is not None:
            return float(np.mod(self.startAngle + 360.0*self.freq*t, 360.0))
        return None

    def dz(self, body, t):
        #z displacement of a body from where it was meshed
        if body not in self.tables:
            sys.exit('ERROR! FOLLOW %s: no motion table for this body (bodies: %s)' % (body, self.bodies))
        x, z = self.tables[body]
        if self.crankTables:
            theta = self.angle(t)
            if theta is None:
                sys.exit('ERROR! FOLLOW %s: no crank log and no START_ANGLE, the crank angle is unknown' % (body))
            return float(np.interp(theta, x, z))
        return float(np.interp(t, x, z))


def readCrankTable(path):
    #reading constant/tables/<body>Crank (theta dz)
    with open(path) as f:
        text = f.read()
    rows = np.array(re.findall(r'\(\s*([-0-9.eE+]+)\s+([-0-9.eE+]+)\s*\)', text), dtype=float)
    return rows[:, 0], rows[:, 1]


#----- cameras

def sliceCamera(normal, origin, bounds, size, barSpace=0.15):
    #camera looking along the slice normal, fitted to the cut through the mesh box
    #barSpace leaves room for the colour bar at the bottom
    centre = [0.5*(bounds[0] + bounds[1]), 0.5*(bounds[2] + bounds[3]), 0.5*(bounds[4] + bounds[5])]
    lengths = [bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]]
    if normal == 'Z':
        width, height, up, axis = lengths[0], lengths[1], [0, 1, 0], 2
    elif normal == 'Y':
        width, height, up, axis = lengths[0], lengths[2], [0, 0, 1], 1
    else:
        width, height, up, axis = lengths[1], lengths[2], [0, 0, 1], 0
    focal = list(centre)
    focal[axis] = origin[axis]
    #shifting the picture up so the colour bar sits below the geometry
    upAxis = up.index(1)
    focal[upAxis] -= barSpace*height
    position = list(focal)
    position[axis] -= 10.0*max(lengths)
    scale = 0.5*max(height*(1.0 + 2.0*barSpace), width*size[1]/size[0])*1.05
    return {'position': position, 'focal': focal, 'up': up, 'scale': scale}


def surfaceCamera(view, bounds, size):
    #fixed views of the whole mesh box
    directions = {'Iso': [1, -1, 1], 'Top': [0, 0, 1], 'Bottom': [0, 0, -1], 'Front': [0, -1, 0],
                  'Back': [0, 1, 0], 'Left': [-1, 0, 0], 'Right': [1, 0, 0]}
    if view not in directions:
        sys.exit('ERROR! pvPostInternalSetup [SURFACE] VIEWS: %s unknown, use %s' % (view, list(directions.keys())))
    centre = np.array([0.5*(bounds[0] + bounds[1]), 0.5*(bounds[2] + bounds[3]), 0.5*(bounds[4] + bounds[5])])
    lengths = np.array([bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]])
    direction = np.array(directions[view], dtype=float)
    direction /= np.linalg.norm(direction)
    up = [0, 1, 0] if view in ('Top', 'Bottom') else [0, 0, 1]
    radius = 0.5*float(np.linalg.norm(lengths))
    return {'position': list(centre + 10.0*radius*direction), 'focal': list(centre), 'up': up,
            'scale': 1.1*radius*max(1.0, size[1]/size[0])}


def fittedCamera(direction, up, bounds, size, barSpace=0.12, margin=1.08):
    #parallel camera from a direction, fitted to the mesh box corners as seen on the screen
    #barSpace leaves room for the colour bar at the bottom
    direction = np.array(direction, dtype=float)
    direction /= np.linalg.norm(direction)
    upVec = np.array(up, dtype=float)
    upVec -= np.dot(upVec, direction)*direction
    upVec /= np.linalg.norm(upVec)
    right = np.cross(upVec, direction)
    corners = np.array([[x, y, z] for x in bounds[0:2] for y in bounds[2:4] for z in bounds[4:6]], dtype=float)
    u, v = corners.dot(right), corners.dot(upVec)
    halfWidth, halfHeight = 0.5*(u.max() - u.min()), 0.5*(v.max() - v.min())
    scale = max(halfHeight*(1.0 + 2.0*barSpace), halfWidth*size[1]/size[0])*margin
    centre = corners.mean(axis=0) - barSpace*halfHeight*2.0*upVec
    radius = float(np.linalg.norm(corners.max(axis=0) - corners.min(axis=0)))
    return {'position': list(centre + 10.0*radius*direction), 'focal': list(centre), 'up': list(upVec), 'scale': scale}


def clipCamera(normal, bounds, size):
    #looking at the cut face from the removed side, tilted so the depth of the half volume shows
    normal = np.array(normal, dtype=float)
    direction = normal + 0.6*(np.ones(3) - np.abs(normal))
    up = [0, 1, 0] if abs(normal[2]) > 0.5 else [0, 0, 1]
    return fittedCamera(direction, up, bounds, size)


def volumeCamera(view, bounds, size):
    #the surface views (Iso, Top, ...) fitted to the mesh box
    directions = {'Iso': [1, -1, 1], 'Top': [0, 0, 1], 'Bottom': [0, 0, -1], 'Front': [0, -1, 0],
                  'Back': [0, 1, 0], 'Left': [-1, 0, 0], 'Right': [1, 0, 0]}
    if view not in directions:
        sys.exit('ERROR! pvPostInternalSetup [VOLUME] VIEWS: %s unknown, use %s' % (view, list(directions.keys())))
    return fittedCamera(directions[view], [0, 1, 0] if view in ('Top', 'Bottom') else [0, 0, 1], bounds, size)


def resampleDims(bounds, points):
    #regular grid for the volume rendering, about the same spacing along every axis
    lengths = np.array([bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]], dtype=float)
    spacing = lengths.max()/max(points - 1, 1)
    return [max(int(round(l/spacing)) + 1, 10) for l in lengths]


def barLabels(lo, hi, n=5):
    #evenly spaced colour bar labels (the automatic ticks overlap on long ranges)
    return [float(v) for v in np.linspace(lo, hi, n)]


#----- paraview

def meshImages(pvs, caseDir, caseName, setup, motion, locInMesh, size, every=None):
    #mesh motion dry run: the moved meshes in EnSightMotion (exported by motionCheckScript) as grey
    #surface-with-edges slices ([SLICE] NAMES / NORMALS / ORIGINS / FOLLOW) and body surfaces ([SURFACE] PATCHES / VIEWS),
    #one image per written crank angle in postProcessing/meshImages/<panel>/<case>_<panel>_NNNNN.png
    begin = time.time()
    cases = glob.glob(os.path.join(caseDir, 'EnSightMotion', '*.case'))
    if not cases:
        sys.exit('ERROR! no EnSightMotion/*.case: run motionCheckScript first (it exports the moved meshes)')
    print('\n\tLoading %s...' % (cases[0]))
    reader = pvs.EnSightReader(CaseFileName=cases[0])
    times = list(reader.TimestepValues)
    if every:
        times = times[::-1][::max(int(every), 1)][::-1]
    if not times:
        sys.exit('ERROR! no moved meshes in EnSightMotion')
    print('\t\t%d meshes from %.6g to %.6g s' % (len(times), times[0], times[-1]))
    reader.UpdatePipeline(times[0])
    blocks = getBlocks(reader)
    volume = pvs.ExtractBlock(Input=reader)
    volume.Selectors = ['/Root/internalMesh']
    volume.UpdatePipeline(times[0])
    bounds = volume.GetDataInformation().GetBounds()

    #panels: (name, source, camera, follow, origin, label)
    panels = []
    for s in getSlices(setup['SLICE']):
        origin = sliceOrigin(s['origin'], bounds, locInMesh)
        cut = pvs.Slice(Input=volume)
        cut.SliceType = 'Plane'
        cut.SliceType.Origin = origin
        cut.SliceType.Normal = normalVectors[s['normal']]
        cut.Triangulatetheslice = 0
        if s['follow']:
            motion.dz(s['follow'], times[0])
        panels.append({'name': s['name'], 'source': cut, 'camera': sliceCamera(s['normal'], origin, bounds, size, 0.0),
                       'follow': s['follow'], 'origin': origin, 'label': s['name']})
    keywords = setup['SURFACE']['PATCHES'].split()
    selectors = ['/Root/%s' % (b) for b in blocks
                 if b != 'internalMesh' and any(k.replace('-', '').lower() in b.lower() for k in keywords)]
    if selectors:
        surface = pvs.ExtractBlock(Input=reader)
        surface.Selectors = selectors
        for view in setup['SURFACE']['VIEWS'].split():
            panels.append({'name': 'surface%s' % (view), 'source': surface, 'camera': surfaceCamera(view, bounds, size),
                           'follow': None, 'origin': None, 'label': 'surface %s' % (view)})

    #one render view per panel, grey cells with dark edges
    for panel in panels:
        view = pvs.CreateRenderView()
        view.ViewSize = size
        view.OrientationAxesVisibility = 0
        try:
            view.UseColorPaletteForBackground = 0
        except AttributeError:
            pass
        view.Background = [1, 1, 1]
        view.CameraParallelProjection = 1
        view.CameraPosition = panel['camera']['position']
        view.CameraFocalPoint = panel['camera']['focal']
        view.CameraViewUp = panel['camera']['up']
        view.CameraParallelScale = panel['camera']['scale']
        display = pvs.Show(panel['source'], view)
        display.ColorArrayName = ['POINTS', '']
        display.SetRepresentationType('Surface With Edges')
        display.AmbientColor = [0.85, 0.85, 0.85]
        display.DiffuseColor = [0.85, 0.85, 0.85]
        display.EdgeColor = [0.15, 0.15, 0.3]
        title = pvs.Text()
        titleDisplay = pvs.Show(title, view)
        titleDisplay.Color = [0, 0, 0]
        titleDisplay.FontSize = 24
        titleDisplay.WindowLocation = 'Upper Left Corner'
        panel['view'], panel['title'] = view, title
        os.makedirs(os.path.join('postProcessing', 'meshImages', panel['name']), exist_ok=True)

    #all panels at each crank angle
    print('\n\tWriting %d mesh images (%d panels x %d meshes)...' % (len(panels)*len(times), len(panels), len(times)))
    for i, t in enumerate(times):
        angle = motion.angle(t)
        stamp = 't = %.5f s' % (t) + ('   crank %.1f deg' % (angle) if angle is not None else '')
        for panel in panels:
            if panel['follow']:
                origin = list(panel['origin'])
                origin[2] += motion.dz(panel['follow'], t)
                panel['source'].SliceType.Origin = origin
            panel['view'].ViewTime = t
            panel['title'].Text = '%s   mesh   %s   %s' % (caseName, panel['label'], stamp)
            pvs.SaveScreenshot(os.path.join('postProcessing', 'meshImages', panel['name'],
                                            '%s_%s_%05d.png' % (caseName, panel['name'], i)),
                               panel['view'], ImageResolution=size)
        print('\t\t%d/%d  %s' % (i + 1, len(times), stamp))
    print('\n\tMesh images in postProcessing/meshImages (%d panels), total time %.1f s' % (len(panels), time.time() - begin))


def main(argv=None):
    #line by line into the log (redirected output is block-buffered: the log stayed empty for hours of rendering)
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except AttributeError:
        pass
    args = parseArgs(argv)
    begin = time.time()
    caseDir = os.getcwd()
    caseName = os.path.basename(caseDir)
    print('\n\t####\t\tpvPostInternal\t\t####\n')
    print('\tProcessing case: %s' % (caseName))

    #case settings
    print('\n\tReading fullCaseSetupDict...')
    config = readConfig(os.path.join(caseDir, 'fullCaseSetupDict'))
    if config.has_section('SETUP_INFO') and config['SETUP_INFO'].get('FLOW', ''):
        compressible = config['SETUP_INFO']['FLOW'] == 'compressible'
    else:
        compressible = config['GLOBAL_SIM_CONTROL']['SOLVER'] in ('rhoPimpleFoam', 'rhoSimpleFoam', 'buoyantSimpleFoam')
    pRef = float(config['GLOBAL_MATERIAL']['CHARGE_PRESSURE']) if compressible else 0.0
    locInMesh = [float(v) for v in config['GLOBAL_REFINEMENT']['LOC_IN_MESH'].split()]
    motion = Motion(caseDir, config)
    steady = config.has_section('SETUP_INFO') and config['SETUP_INFO'].get('STEADY', 'no') == 'yes'

    print('\tGetting pvPostInternalSetup...')
    setup = getPVSetup(caseDir)
    for section in setup.sections():
        print('\t\t%s' % (section))
        for key, value in setup[section].items():
            print('\t\t\t%s: %s' % (key, value))
    varDict = getVariableDict(setup, pRef)
    mainSetup = setup['PV_POST_MAIN']
    size = [int(v) for v in mainSetup['RESOLUTION'].split()]
    every = args.every if args.every is not None else int(mainSetup['EVERY'])
    lastCycles = args.lastCycles if args.lastCycles is not None else mainSetup['LAST_CYCLES']
    scalarBar = mainSetup['CREATE_SCALAR_BAR'].lower() == 'true'
    kinds = getPanelTypes(mainSetup, args.panels)

    import paraview.simple as pvs
    pvs._DisableFirstRenderCameraReset()
    if args.mesh:
        meshImages(pvs, caseDir, caseName, setup, motion, locInMesh, size, args.every)
        return

    #loading the data
    cases = glob.glob(os.path.join(caseDir, 'EnSight', '*.case'))
    if not cases:
        sys.exit('ERROR! no EnSight/*.case: has foamToEnsight run (EXPORT = ensight)?')
    print('\n\tLoading %s...' % (cases[0]))
    reader = pvs.EnSightReader(CaseFileName=cases[0])
    times = selectTimes(list(reader.TimestepValues), motion.freq, every, lastCycles)
    if not times:
        sys.exit('ERROR! no written times to post-process')
    print('\t\t%d times from %.6g to %.6g s (every %d, last cycles %s)' % (len(times), times[0], times[-1], every, lastCycles))
    fields = list(reader.CellArrays)
    print('\t\tfields: %s' % (' '.join(fields)))
    points = pvs.CellDatatoPointData(Input=reader)

    reader.UpdatePipeline(times[0])
    blocks = getBlocks(reader)
    print('\t\tblocks: %s' % (' '.join(blocks)))
    volume = pvs.ExtractBlock(Input=points)
    volume.Selectors = ['/Root/internalMesh']
    volume.UpdatePipeline(times[0])
    bounds = volume.GetDataInformation().GetBounds()
    print('\t\tmesh box: %s' % (' '.join('%.4g' % b for b in bounds)))

    #building the panels
    panels = []
    if 'slice' in kinds:
        slices = getSlices(setup['SLICE'])
        sliceVars = setup['SLICE']['VARIABLES'].split()
        for s in slices:
            origin = sliceOrigin(s['origin'], bounds, locInMesh)
            cut = pvs.Slice(Input=volume)
            cut.SliceType = 'Plane'
            cut.SliceType.Origin = origin
            cut.SliceType.Normal = normalVectors[s['normal']]
            #checking the body has a motion table before writing any images
            if s['follow']:
                motion.dz(s['follow'], times[0])
            camera = sliceCamera(s['normal'], origin, bounds, size)
            for variable in sliceVars:
                if variable not in varDict['sliceVariables']:
                    print('\t\t%s is not in the slice variable dict, skipping' % (variable))
                    continue
                if not neededFields(varDict['sliceVariables'][variable]['equation'], fields):
                    print('\t\t%s: field not in the case, skipping' % (variable))
                    continue
                panels.append({'name': '%s_%s' % (s['name'], variable), 'source': cut, 'variable': variable,
                               'spec': varDict['sliceVariables'][variable], 'camera': camera, 'kind': 'slice',
                               'follow': s['follow'], 'origin': origin, 'label': s['name']})
    if 'surface' in kinds:
        keywords = setup['SURFACE']['PATCHES'].split()
        selectors = ['/Root/%s' % (b) for b in blocks
                     if b != 'internalMesh' and any(k.replace('-', '').lower() in b.lower() for k in keywords)]
        if not selectors:
            print('\t\tno patches match SURFACE PATCHES %s, no surface images' % (keywords))
        else:
            print('\t\tsurface patches: %s' % (' '.join(selectors)))
            surface = pvs.ExtractBlock(Input=points)
            surface.Selectors = selectors
            for view in setup['SURFACE']['VIEWS'].split():
                camera = surfaceCamera(view, bounds, size)
                for variable in setup['SURFACE']['VARIABLES'].split():
                    if variable not in varDict['surfaceVariables']:
                        print('\t\t%s is not in the surface variable dict, skipping' % (variable))
                        continue
                    if not neededFields(varDict['surfaceVariables'][variable]['equation'], fields):
                        print('\t\t%s: field not in the case, skipping' % (variable))
                        continue
                    panels.append({'name': 'surface%s_%s' % (view, variable), 'source': surface,
                                   'variable': variable, 'spec': varDict['surfaceVariables'][variable], 'kind': 'surface',
                                   'camera': camera, 'follow': None, 'origin': None, 'label': 'surface %s' % (view)})
    #context for the volume renderings: the bodies (CONTEXT_PATCHES) in light grey, the mesh box outline with OUTLINE
    context = None
    if 'volume' in kinds:
        keywords = setup['VOLUME'].get('CONTEXT_PATCHES', '').split()
        selectors = ['/Root/%s' % (b) for b in blocks
                     if b != 'internalMesh' and any(k.replace('-', '').lower() in b.lower() for k in keywords)]
        if selectors:
            context = pvs.ExtractBlock(Input=reader)
            context.Selectors = selectors
        #volume rendering needs a regular grid (the unstructured mapper skips snappy's polyhedra): the gas is
        #resampled to RESAMPLE points along the longest box side, points outside the gas are 0 (clear in the ramps)
        dims = resampleDims(bounds, int(setup['VOLUME'].get('RESAMPLE', '250')))
        print('\t\tvolume rendering on a %d x %d x %d grid' % tuple(dims))
        grid = pvs.ResampleToImage(Input=volume)
        grid.UseInputBounds = 1
        grid.SamplingDimensions = dims
        for view in setup['VOLUME']['VIEWS'].split():
            camera = volumeCamera(view, bounds, size)
            for variable in setup['VOLUME']['VARIABLES'].split():
                if not volumeVariable(variable, varDict, fields):
                    continue
                panels.append({'name': 'volume%s_%s' % (view, variable), 'source': grid, 'variable': variable,
                               'spec': varDict['sliceVariables'][variable], 'camera': camera, 'kind': 'volume',
                               'follow': None, 'origin': None, 'label': 'volume %s' % (view)})
    if 'clip' in kinds:
        for c in getClips(setup['CLIP']):
            origin = sliceOrigin(c['origin'], bounds, locInMesh)
            half = pvs.Clip(Input=volume)
            half.ClipType = 'Plane'
            half.ClipType.Origin = origin
            half.ClipType.Normal = c['normal']
            #keeping the half behind the plane (away from the normal and the camera)
            half.Invert = 1
            half.Crinkleclip = 0
            camera = clipCamera(c['normal'], bounds, size)
            for variable in setup['CLIP']['VARIABLES'].split():
                if not volumeVariable(variable, varDict, fields):
                    continue
                panels.append({'name': '%s_%s' % (c['name'], variable), 'source': half, 'variable': variable,
                               'spec': varDict['sliceVariables'][variable], 'camera': camera, 'kind': 'clip',
                               'follow': None, 'origin': origin, 'label': c['name']})
    if not panels:
        sys.exit('ERROR! nothing to plot: check SLICE_IMG / SURFACE_IMG / VOLUME_IMG / CLIP_IMG and the variables')

    #colour ranges fixed over the selected times, slices use the whole volume so every slice has the same range
    print('\n\tColour ranges over the run...')
    rangeTime = time.time()
    rangeCalcs = {}
    for panel in panels:
        calc = pvs.Calculator(Input=panel['source'])
        calc.AttributeType = 'Point Data'
        calc.ResultArrayName = panel['variable']
        calc.Function = str(panel['spec']['equation'])
        panel['calc'] = calc
        panel['rangeKey'] = (panel['variable'], 'surface' if panel['kind'] == 'surface' else 'volume')
        if panel['rangeKey'] not in rangeCalcs:
            calc = pvs.Calculator(Input=panel['source'] if panel['rangeKey'][1] == 'surface' else volume)
            calc.AttributeType = 'Point Data'
            calc.ResultArrayName = panel['variable']
            calc.Function = str(panel['spec']['equation'])
            rangeCalcs[panel['rangeKey']] = calc
    ranges = {}
    for t in times:
        for key, calc in rangeCalcs.items():
            calc.UpdatePipeline(t)
            info = calc.GetPointDataInformation().GetArray(key[0])
            if info is None:
                continue
            lo, hi = info.GetRange(0)
            old = ranges.get(key, (lo, hi))
            ranges[key] = (min(old[0], lo), max(old[1], hi))
    for panel in panels:
        lo, hi = ranges.get(panel['rangeKey'], (0.0, 1.0))
        mode = panel['spec']['range']
        if isinstance(mode, list):
            lo, hi = float(mode[0]), float(mode[1])
        elif mode == 'symmetric':
            hi = max(abs(lo), abs(hi))
            lo = -hi
        if hi <= lo:
            hi = lo + max(abs(lo)*1e-6, 1e-12)
        panel['range'] = (lo, hi)
        print('\t\t%-28s %.6g to %.6g' % (panel['name'], lo, hi))
    for calc in rangeCalcs.values():
        pvs.Delete(calc)
    print('\t\tdone in %.1f s' % (time.time() - rangeTime))

    #one render view per panel
    for panel in panels:
        view = pvs.CreateRenderView()
        view.ViewSize = size
        view.OrientationAxesVisibility = 0
        try:
            view.UseColorPaletteForBackground = 0
        except AttributeError:
            pass
        view.Background = [1, 1, 1]
        view.CameraParallelProjection = 1
        view.CameraPosition = panel['camera']['position']
        view.CameraFocalPoint = panel['camera']['focal']
        view.CameraViewUp = panel['camera']['up']
        view.CameraParallelScale = panel['camera']['scale']
        display = pvs.Show(panel['calc'], view)
        display.Representation = 'Volume' if panel['kind'] == 'volume' else 'Surface'
        display.UseSeparateColorMap = 1
        pvs.ColorBy(display, ('POINTS', panel['variable']), separate=True)
        lut = display.LookupTable
        lut.ApplyPreset(panel['spec']['color'], True)
        lut.AutomaticRescaleRangeMode = 'Never'
        lut.RescaleTransferFunction(*panel['range'])
        panel['pwf'] = None
        if panel['kind'] == 'volume':
            panel['pwf'] = volumeOpacity(pvs, display, panel, setup['VOLUME'], bounds)
            if context is not None:
                contextDisplay = pvs.Show(context, view)
                contextDisplay.ColorArrayName = ['POINTS', '']
                contextDisplay.AmbientColor = [0.6, 0.6, 0.6]
                contextDisplay.DiffuseColor = [0.6, 0.6, 0.6]
                contextDisplay.Opacity = float(setup['VOLUME'].get('CONTEXT_OPACITY', '0.15'))
            #mesh box outline only on request (OUTLINE = True)
            if setup['VOLUME'].get('OUTLINE', 'False').lower() == 'true':
                outline = pvs.Show(volume, view)
                outline.SetRepresentationType('Outline')
                outline.ColorArrayName = ['POINTS', '']
                outline.AmbientColor = [0.3, 0.3, 0.3]
                outline.DiffuseColor = [0.3, 0.3, 0.3]
        display.SetScalarBarVisibility(view, scalarBar)
        if scalarBar:
            bar = pvs.GetScalarBar(lut, view)
            bar.Title = panel['spec']['label']
            bar.ComponentTitle = ''
            bar.TitleColor = [0, 0, 0]
            bar.LabelColor = [0, 0, 0]
            bar.TitleFontSize = 22
            bar.LabelFontSize = 20
            bar.AutomaticLabelFormat = 0
            bar.UseCustomLabels = 1
            bar.CustomLabels = barLabels(*panel['range'])
            bar.AddRangeLabels = 0
            bar.LabelFormat = '%-#6.4g'
            bar.RangeLabelFormat = '%-#6.4g'
            bar.Orientation = 'Horizontal'
            bar.WindowLocation = 'Lower Center'
            bar.ScalarBarLength = 0.4
        title = pvs.Text()
        titleDisplay = pvs.Show(title, view)
        titleDisplay.Color = [0, 0, 0]
        titleDisplay.FontSize = 24
        titleDisplay.WindowLocation = 'Upper Left Corner'
        panel['view'], panel['lut'], panel['title'] = view, lut, title
        os.makedirs(os.path.join('postProcessing', 'images', panel['name']), exist_ok=True)

    #writing the images, all panels per time so the data is only read once per time
    print('\n\tWriting %d images (%d panels x %d times)...' % (len(panels)*len(times), len(panels), len(times)))
    for i, t in enumerate(times):
        angle = motion.angle(t)
        if steady:
            stamp = 'iteration %d' % (int(round(t)))
        else:
            stamp = 't = %.5f s' % (t) + ('   crank %.1f deg' % (angle) if angle is not None else '')
        for panel in panels:
            if panel['follow']:
                origin = list(panel['origin'])
                origin[2] += motion.dz(panel['follow'], t)
                panel['source'].SliceType.Origin = origin
            panel['lut'].RescaleTransferFunction(*panel['range'])
            if panel['pwf'] is not None:
                panel['pwf'].Points = panel['pwfPoints']
            panel['view'].ViewTime = t
            panel['title'].Text = '%s   %s   %s   %s' % (caseName, panel['label'], panel['variable'], stamp)
            fileName = os.path.join('postProcessing', 'images', panel['name'],
                                    '%s_%s_%05d.png' % (caseName, panel['name'], i))
            pvs.SaveScreenshot(fileName, panel['view'], ImageResolution=size)
        print('\t\t%d/%d  %s' % (i + 1, len(times), stamp))
    print('\n\tImages in postProcessing/images (%d panels), total time %.1f s' % (len(panels), time.time() - begin))


def volumeVariable(variable, varDict, fields):
    #volume and clip panels use the slice variables (volume fields)
    if variable not in varDict['sliceVariables']:
        print('\t\t%s is not in the slice variable dict, skipping' % (variable))
        return False
    if not neededFields(varDict['sliceVariables'][variable]['equation'], fields):
        print('\t\t%s: field not in the case, skipping' % (variable))
        return False
    return True


def opacityPoints(lo, hi, ramp, opacity, symmetric=False):
    #opacity transfer function points (value, opacity, midpoint, sharpness)
    #rising: clear at the low end, opaque at the high end (hot gas, fast flow), falling: the other way,
    #ends: clear in the middle, opaque at both ends (signed variables), auto: ends for symmetric ranges, else rising
    if ramp == 'auto':
        ramp = 'ends' if symmetric else 'rising'
    mid = 0.5*(lo + hi)
    if ramp == 'rising':
        return [lo, 0.0, 0.5, 0.0, hi, opacity, 0.5, 0.0]
    if ramp == 'falling':
        return [lo, opacity, 0.5, 0.0, hi, 0.0, 0.5, 0.0]
    if ramp == 'ends':
        return [lo, opacity, 0.5, 0.0, mid, 0.0, 0.5, 0.0, hi, opacity, 0.5, 0.0]
    sys.exit('ERROR! pvPostInternalSetup [VOLUME] OPACITY_RAMP must be auto, rising, falling or ends!')


def volumeOpacity(pvs, display, panel, volumeSetup, bounds):
    #opacity of the volume rendering: ramp over the fixed colour range, unit distance from the mesh box
    #(auto = 0.05 x the shortest box side, so thin gas spaces still show; a larger UNIT_DISTANCE is more transparent)
    ramp = volumeSetup.get('OPACITY_RAMP', 'auto').lower()
    opacity = float(volumeSetup.get('OPACITY', '0.8'))
    unit = volumeSetup.get('UNIT_DISTANCE', 'auto')
    shortest = min(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
    display.ScalarOpacityUnitDistance = 0.05*shortest if unit.lower() == 'auto' else float(unit)
    try:
        pwf = display.ScalarOpacityFunction
    except AttributeError:
        pwf = None
    if pwf is None:
        pwf = pvs.GetOpacityTransferFunction(panel['variable'], display, separate=True)
        display.ScalarOpacityFunction = pwf
    lo, hi = panel['range']
    panel['pwfPoints'] = opacityPoints(lo, hi, ramp, opacity, panel['spec']['range'] == 'symmetric')
    pwf.Points = panel['pwfPoints']
    pwf.AllowDuplicateScalars = 1
    return pwf


def getBlocks(source):
    #block names as the selectors use them (data hierarchy in ParaView 5.10+, old block names otherwise)
    info = source.GetDataInformation().DataInformation
    names = []
    try:
        hierarchy = info.GetHierarchy()
        for node in hierarchy.GetChildNodes(hierarchy.GetRootNode(), False):
            names.append(hierarchy.GetNodeName(node))
    except AttributeError:
        for i in range(1, info.GetNumberOfDataSets() + 1):
            name = info.GetBlockName(i)
            if name:
                names.append(name.replace('-', ''))
    return names


if __name__ == '__main__':
    main()
