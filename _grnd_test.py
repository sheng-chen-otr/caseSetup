import os, sys, shutil, tempfile
import numpy as np

# import the modules under test
import utilities as U
import writeSystem as W

REPO = os.path.dirname(os.path.abspath(__file__))
fails = []

def check(name, cond):
    print(('PASS' if cond else 'FAIL') + ' - ' + name)
    if not cond:
        fails.append(name)

# ---------------------------------------------------------------------------
# helpers to write/read a tiny ascii STL cube with a named solid (PID)
# ---------------------------------------------------------------------------
def write_cube_stl(path, name, lo=(0,0,0), hi=(1,1,1)):
    x0,y0,z0 = lo; x1,y1,z1 = hi
    v = [(x0,y0,z0),(x1,y0,z0),(x1,y1,z0),(x0,y1,z0),
         (x0,y0,z1),(x1,y0,z1),(x1,y1,z1),(x0,y1,z1)]
    tris = [(0,1,2),(0,2,3),(4,6,5),(4,7,6),(0,4,5),(0,5,1),
            (1,5,6),(1,6,2),(2,6,7),(2,7,3),(3,7,4),(3,4,0)]
    with open(path,'w') as f:
        f.write('solid %s\n' % name)
        for a,b,c in tris:
            f.write(' facet normal 0 0 0\n  outer loop\n')
            for idx in (a,b,c):
                f.write('   vertex %f %f %f\n' % v[idx])
            f.write('  endloop\n endfacet\n')
        f.write('endsolid %s\n' % name)

def read_stl_verts(path):
    verts = []; solids = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith('vertex'):
                p = s.split()
                verts.append([float(p[1]),float(p[2]),float(p[3])])
            elif s.startswith('solid'):
                solids.append(s.split(None,1)[1] if len(s.split(None,1))>1 else '')
    return np.array(verts), solids

# ---------------------------------------------------------------------------
# TEST 1: ground STL transform matches the domain transform for a ride-height
#         child point: translate by heave, then roll(x), then pitch(y) about
#         the heaved REFCOR - independent numpy reference
# ---------------------------------------------------------------------------
def test_stl_tilt():
    tmp = tempfile.mkdtemp()
    src = os.path.join(tmp,'GRND-belt.stl')
    write_cube_stl(src,'GRND-belt',lo=(-2,-2,-0.5),hi=(2,2,0.5))

    baseREFCOR = [1.0, 0.0, 0.3]
    tunnel_heave = 0.05
    roll = np.radians(2.0)   # about x
    pitch = np.radians(-3.0) # about y
    pivot = [baseREFCOR[0], baseREFCOR[1], baseREFCOR[2] + tunnel_heave]

    # replicate the transformGeom GRND branch exactly
    out = os.path.join(tmp,'child.stl')
    tmp2 = out + '.grndtmp'
    U.transformGeometryPreservePID(src, outputFile=out, translation=[0,0,tunnel_heave])
    U.transformGeometryPreservePID(out, outputFile=tmp2, rotation_rvec=[roll,0,0], pivot=pivot)
    U.transformGeometryPreservePID(tmp2, outputFile=out, rotation_rvec=[0,pitch,0], pivot=pivot)

    got, solids = read_stl_verts(out)
    src_v, _ = read_stl_verts(src)

    # reference: same order the blockMesh points get (Rx then Ry about pivot)
    Rx = np.array([[1,0,0],[0,np.cos(roll),-np.sin(roll)],[0,np.sin(roll),np.cos(roll)]])
    Ry = np.array([[np.cos(pitch),0,np.sin(pitch)],[0,1,0],[-np.sin(pitch),0,np.cos(pitch)]])
    p = np.array(pivot)
    ref = []
    for v in src_v:
        v1 = v + np.array([0,0,tunnel_heave])
        v2 = Rx @ (v1 - p) + p
        v3 = Ry @ (v2 - p) + p
        ref.append(v3)
    ref = np.array(ref)

    check('stl tilt matches Rx->Ry about heaved REFCOR', np.allclose(got, ref, atol=1e-4))
    check('stl PID/solid name preserved', solids and solids[0] == 'GRND-belt')
    shutil.rmtree(tmp)

# ---------------------------------------------------------------------------
# TEST 2: createGroundPatch writes topoSetDictGround + createPatchDict
# ---------------------------------------------------------------------------
def test_ground_patch():
    tmp = tempfile.mkdtemp()
    cwd = os.getcwd(); os.chdir(tmp)
    try:
        os.makedirs('constant/triSurface', exist_ok=True)
        geomDict = {'GRND-belt.stl':{}, 'car.stl':{}, 'GRND-plate.stl':{}}
        fcsd = {'GLOBAL_REFINEMENT':{'LOC_IN_MESH':['-1','0','0.5'],
                                     'TEMPLATE_TYPE':['snappyHexMesh']},
                'GRND-belt':{'GRND_TYPE':['moving'],'GRND_VEL':['default'],'GRND_WALL_MODEL':['high']},
                'GRND-plate':{'GRND_TYPE':['slip'],'GRND_VEL':['default'],'GRND_WALL_MODEL':['high']}}
        check('hasGroundZones true', W.hasGroundZones(fcsd) is True)
        check('hasGroundZones false', W.hasGroundZones({'BC_SETUP':{}}) is False)

        ret = W.createGroundPatch('x', geomDict, fcsd)
        check('createGroundPatch returned True', ret is True)
        tsd = open('system/topoSetDictGround').read()
        cpd = open('system/createPatchDict').read()

        check('topoSet has belt faceSet', 'groundZone-GRND-belt' in tsd)
        check('topoSet has plate faceSet', 'groundZone-GRND-plate' in tsd)
        check('topoSet uses surfaceToCell', 'surfaceToCell' in tsd)
        check('topoSet subsets z-min patch', 'patchToFace' in tsd and 'z-min' in tsd)
        check('topoSet uses LOC_IN_MESH outsidePoints', '(-1 0 0.5)' in tsd)
        check('topoSet points at triSurface file', 'constant/triSurface/GRND-belt.stl' in tsd)
        check('topoSet ignores non-GRND geom', 'car.stl' not in tsd)
        check('createPatch builds belt patch', 'groundZone-GRND-belt' in cpd and 'type wall' in cpd)
        check('createPatch builds plate patch', 'groundZone-GRND-plate' in cpd)
    finally:
        os.chdir(cwd); shutil.rmtree(tmp)

# ---------------------------------------------------------------------------
# TEST 3: writeBoundaries emits per-zone ground caseProperties entries using
#         the real default caseProperties template
# ---------------------------------------------------------------------------
def test_write_boundaries():
    tmp = tempfile.mkdtemp()
    cwd = os.getcwd(); os.chdir(tmp)
    try:
        os.makedirs('system', exist_ok=True)
        gsc = {'SIM_SYM':['full'],'INIT_P':['default'],'INIT_K':['default'],
               'INIT_OMEGA':['default'],'INIT_NUT':['default'],'INIT_NUTILDA':['default']}
        fcsd = {'BC_SETUP':{'INLET_MAG':['30'],'YAW':['0'],'DOMAIN_PITCH':['0'],
                            'DOMAIN_ROLL':['0'],'REFCOR':['0','0','0'],'GROUND':['stationary']},
                'GLOBAL_SIM_CONTROL':gsc,
                'GRND-belt':{'GRND_TYPE':['moving'],'GRND_VEL':['default'],'GRND_WALL_MODEL':['high']},
                'GRND-plate':{'GRND_TYPE':['slip'],'GRND_VEL':['default'],'GRND_WALL_MODEL':['high']},
                'GRND-fixed':{'GRND_TYPE':['noslip'],'GRND_VEL':['default'],'GRND_WALL_MODEL':['high']}}
        geomDict = {}  # no geometry patches, just domain + ground zones
        templateLoc = os.path.join(REPO,'setupTemplates','default','defaultDicts')
        W.writeBoundaries(templateLoc, geomDict, fcsd)
        cp = open('system/caseProperties').read()

        check('belt moving-wall entry', 'groundZone-GRND-belt' in cp and 'motion moving' in cp)
        check('plate slip entry', 'groundZone-GRND-plate' in cp and 'type slip' in cp)
        check('fixed noSlip stationary entry',
              'groundZone-GRND-fixed' in cp and 'motion stationary' in cp)
        check('belt patch pattern', '("groundZone-GRND-belt.*")' in cp)
    finally:
        os.chdir(cwd); shutil.rmtree(tmp)

test_stl_tilt()
test_ground_patch()
test_write_boundaries()

print('\n' + ('ALL TESTS PASSED' if not fails else 'FAILURES: %s' % fails))
sys.exit(1 if fails else 0)
