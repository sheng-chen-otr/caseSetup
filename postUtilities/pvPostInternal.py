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
#panels are slices x variables and surface views x variables, set up like pvPost with pvPostInternalSetup
#colour ranges are kept fixed over the run, slices can follow a moving body with FOLLOW
#images go to postProcessing/images/<panel>, imagesToMovies.py makes the videos
#usage (in the case directory): pvbatch --force-offscreen-rendering pvPostInternal.py [--every N] [--lastCycles N]

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
    return parser.parse_args(argv)


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
    if mainSetup['SLICE_IMG'].lower() == 'true':
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
                               'spec': varDict['sliceVariables'][variable], 'camera': camera,
                               'follow': s['follow'], 'origin': origin, 'label': s['name']})
    if mainSetup['SURFACE_IMG'].lower() == 'true':
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
                                   'variable': variable, 'spec': varDict['surfaceVariables'][variable],
                                   'camera': camera, 'follow': None, 'origin': None, 'label': 'surface %s' % (view)})
    if not panels:
        sys.exit('ERROR! nothing to plot: check SLICE_IMG / SURFACE_IMG and the variables')

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
        panel['rangeKey'] = (panel['variable'], 'volume' if panel['origin'] is not None else 'surface')
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
        display.Representation = 'Surface'
        display.UseSeparateColorMap = 1
        pvs.ColorBy(display, ('POINTS', panel['variable']), separate=True)
        lut = display.LookupTable
        lut.ApplyPreset(panel['spec']['color'], True)
        lut.AutomaticRescaleRangeMode = 'Never'
        lut.RescaleTransferFunction(*panel['range'])
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
            panel['view'].ViewTime = t
            panel['title'].Text = '%s   %s   %s   %s' % (caseName, panel['label'], panel['variable'], stamp)
            fileName = os.path.join('postProcessing', 'images', panel['name'],
                                    '%s_%s_%05d.png' % (caseName, panel['name'], i))
            pvs.SaveScreenshot(fileName, panel['view'], ImageResolution=size)
        print('\t\t%d/%d  %s' % (i + 1, len(times), stamp))
    print('\n\tImages in postProcessing/images (%d panels), total time %.1f s' % (len(panels), time.time() - begin))


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
