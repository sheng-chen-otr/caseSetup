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

import sys
import fileinput
import os
import getopt
import re
import glob
import argparse as argparse
import configparser
import numpy as np
import pandas as pd
import math
from glob import glob
import subprocess as sp
from writeSystem import *
from utilities import *
from writeConstant import *
from writeScripts import *
import shutil
import gzip
import re
import multiprocessing
import setupCore
from setupCore import __version__

#these lines get the path which the program is run from
path = os.path.split(os.getcwd())[0] #path of the case
case = os.path.split(os.getcwd())[1] #case name
jobPath = os.path.abspath(os.path.join(path,os.pardir)) #path of the job i.e. 100001

# Get the directory where the executable is running from
if getattr(sys, 'frozen', False):
    # If running as a PyInstaller bundle (i.e., an executable)
    execDir = os.path.dirname(sys.executable)
else:
    # If running as a script
    execDir = os.path.dirname(os.path.realpath(__file__))

#location of the template file, should be in the same directory as _internal when compiled
templateBaseLoc = "%s/setupTemplates" % (execDir) 
print('EXECUTING FROM: ' + execDir)

#checks the location of the templates to make sure they are good
if not os.path.isdir(templateBaseLoc):
    sys.exit("ERROR: Template path incorrect!")

parser = argparse.ArgumentParser(prog='caseSetup-v%s' % __version__,description='Set us the case based on settings written out in the caseSetup')
                    
parser.add_argument("-s","--setup", default='default', 
                   
                    help='Identifies which setup templates to use. Available setups are: default, default-local, otr, otr2606')
# parser.add_argument("-m","--mesh", action="store_true",
                    # help='Sets up meshing dicts and links geometry files. Will copy template case files if the trial folder is empty.')
parser.add_argument('-d',"--controlDict", action="store_true", 
                    help='Writes only controlDict for solver.')
parser.add_argument("--new", action="store_true", 
                    help='Writes a new caseSetup, will overwrite what is currently there if it exists!')
parser.add_argument("--modules", action="store_true", 
                    help='Shows all possible modules.')
parser.add_argument("--postProDict", action="store_true", 
                    help='Copies post-processing config file into case folder.')
parser.add_argument("--rideHeightMode", action="store_true", 
                    help='Run in ride height mode, does not link geometry files.')
parser.add_argument("--animateSuspension", action="store_true",
                    help='Generate a suspension-motion GIF (front/side/top views) across '
                         'the ride-height map using the kinematic solver. Requires '
                         '[RIDE_HEIGHT_SETUP] USE_KINEMATIC_SOLVER true.')
parser.add_argument("--gsheetPrep", action="store_true",
                    help='After setup completes, reserve/prepare gsheet row(s) for this case (and its '
                         'ride-height children, if any, in parent-then-child order) with the setup-time '
                         'values already known (velocity, yaw, density, symmetry, wheelbase), without '
                         'touching result columns. Intended to run before the case is submitted.')
args = parser.parse_args()
CONTROLDICT = args.controlDict
NEWCS = args.new
SETUP = args.setup
PR_MODULES = args.modules
updateCaseSetupFlag = False

#available addon keywords
addonKeyWords = ['POR','FAN','REFX','WAKE','GEOMX','ROTA','MOVG','IDOM','MRFG','GRND']

#(prefix, defaultBCTemplates file, section header), in the order they are checked.
#FAN is a placeholder for checkFan (extra validation); REFX is handled by checkRefinements
addonTemplates = [
    ('GEOMX', 'defaultGeomx',          'GEOMX_SETUP'),
    ('POR',   'defaultPorous',         'POR_SETUP'),
    ('FAN',   'defaultFan',            'FAN_SETUP'),
    ('MRFG',  'defaultMRFG',           'MRF_SETUP'),
    ('ROTA',  'defaultWheel',          'WH_SETUP'),
    ('IDOM',  'defaultInternalDomain', 'INTERNAL_DOMAIN'),
    ('GRND',  'defaultGrnd',           'GRND_SETUP'),
]


#getting default values from template
def main():
    global updateCaseSetupFlag
    titleText = '''\t##############################\n\t\tcaseSetup-v%s\t\n\t##############################''' % __version__
    print(titleText)
    getTemplateType(SETUP)
    caseDir = os.path.join(path, case)



    defaultDict = setupCore.getGlobalDefaults(templateLoc, PR_MODULES)
    if NEWCS == True:
        setupCore.writeNewCaseSetup(defaultDict, caseDir)
    else:
        #copy over the pvPostSetup
        if args.postProDict:
            copyPvPostSetup()
            return
        caseSetupDict,writeCaseSetupDict,fullCaseSetupDict,flag = setupCore.getCaseSetup(
            defaultDict, caseDir, addonKeyWords)
        updateCaseSetupFlag = updateCaseSetupFlag or flag

        #srf cornering: force a full (non-symmetric) model in-memory, leaves caseSetup SIM_SYM alone
        if 'CORNERING_SETUP' in fullCaseSetupDict and \
           fullCaseSetupDict['CORNERING_SETUP']['RUN_CORNERING'][0].lower() == 'true':
            print('\tCornering (SRF) enabled: validating domain and simulation symmetry...')
            checkCorneringDomain(fullCaseSetupDict)
            if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0].lower() == 'half':
                print('\t\tWARNING! SIM_SYM=half is invalid for cornering (the flow field is not '
                      'laterally symmetric). Forcing SIM_SYM=full for this run.')
                fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'] = ['full']

        #animateSuspension: link geometry, make the animation, then skip the rest
        if args.animateSuspension:
            from rideHeightUtils.suspensionAnimator import generateSuspensionAnimation
            print('\n\tanimateSuspension enabled: linking geometry then generating the '
                  'suspension motion animation only (the rest of caseSetup is skipped)...')
            writeCaseSetupDict,geomDict,fullCaseSetupDict = getGeometry(fullCaseSetupDict,writeCaseSetupDict)
            setupCore.linkGeomFiles(geomDict, caseDir, jobPath)
            generateSuspensionAnimation(fullCaseSetupDict, caseDir=os.path.join(path, case))
            print('\nSuspension animation completed.')
            return

        writeCaseSetupDict,geomDict,fullCaseSetupDict = getGeometry(fullCaseSetupDict,writeCaseSetupDict) 
        
        setupCore.cleanUpCaseSetup(geomDict,writeCaseSetupDict,fullCaseSetupDict,addonKeyWords)
        setupCore.writeToCaseSetup(writeCaseSetupDict, caseDir, 'caseSetup', updateCaseSetupFlag)
        if args.rideHeightMode == False:
            setupCore.linkGeomFiles(geomDict, caseDir, jobPath)
                
        geomDict ,fullCaseSetupDict = writeSnappy(geomDict,fullCaseSetupDict)
        
        writeSurfaceFeatureExtract(templateLoc,geomDict,fullCaseSetupDict)
        writeBlockMesh(templateLoc,fullCaseSetupDict)
        transformBlockMesh(fullCaseSetupDict)
        writeDecomposeParDict(templateLoc, fullCaseSetupDict)
        writeControlDict(templateLoc, fullCaseSetupDict)
        copyFunctionObjects(templateLoc,foList,geomDict,fullCaseSetupDict)
        writeSurfaceFieldAverage(geomDict,fullCaseSetupDict)
        writeBoundaries(templateLoc,geomDict,fullCaseSetupDict)
        
        writePostProSurfaceList(fullCaseSetupDict)
        writeSchemes(templateLoc,fullCaseSetupDict)
        writeSolution(templateLoc,fullCaseSetupDict)
        
        writeForceCoeff(templateLoc,geomDict,fullCaseSetupDict)
        writeOptions(templateLoc,geomDict,fullCaseSetupDict)
        writeMRFG(templateLoc,geomDict,fullCaseSetupDict)
        createTopoSet(templateLoc,geomDict,fullCaseSetupDict)
        createGroundPatch(templateLoc,geomDict,fullCaseSetupDict)
        writeSRFProperties(templateLoc,fullCaseSetupDict)
        writeSurfaces(templateLoc, geomDict,fullCaseSetupDict)
        writeCorneringRelativeFields(fullCaseSetupDict)
        writeTransportProperties(templateLoc, fullCaseSetupDict)
        writeTurbulenceProperties(templateLoc, fullCaseSetupDict)
        setupCore.writeToCaseSetup(fullCaseSetupDict, caseDir, 'fullCaseSetupDict', updateCaseSetupFlag)
        
        if args.rideHeightMode:
            makeScripts(templateLoc,fullCaseSetupDict)
            copyScripts(templateLoc, fullCaseSetupDict,case)
        else:
            if fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RUN_RIDE_HEIGHT'][0].lower() == 'true':
                rideHeights = calculateRideHeights(fullCaseSetupDict)
                rideHeights = createRideHeightCases(rideHeights, fullCaseSetupDict)
                rideHeights = transformGeom(fullCaseSetupDict,rideHeights,geomDict)
                runRHCaseSetup(rideHeights,fullCaseSetupDict,SETUP)
                writeExecuteAllScript(rideHeights,fullCaseSetupDict)
            
            else:
                makeScripts(templateLoc,fullCaseSetupDict)
                copyScripts(templateLoc, fullCaseSetupDict,case)
            #tilt the base ground zones (child points are tilted+heaved in transformGeom)
            transformGroundGeom(geomDict,fullCaseSetupDict)

            if args.gsheetPrep:
                print('\n\tPreparing gsheet row(s)...')
                try:
                    gsheetExportScript = os.path.join(execDir, 'postUtilities', 'gsheetExport.py')
                    sp.run([sys.executable, gsheetExportScript, '--prep'], cwd=os.path.join(path, case), check=True)
                except Exception as e:
                    print('\tWARNING! Unable to prepare gsheet row(s): %s' % (e))
        
        
        
        #getClusterType(templateLoc,fullCaseSetupDict)
        print('\nCase setup completed successfully!')
  

       
   
def getTemplateType(SETUP):
    global templateLoc
    templateLoc = "%s/%s" % (templateBaseLoc,SETUP)
    
    #checking if setup template exists
    if os.path.isdir(templateLoc):
        print('\n\tUsing %s setup template in path: %s' % (SETUP,templateLoc))
    else:
        print('ERROR! %s is not a valid setup template location!' % (templateLoc))
        exit()
   
def resolveFanGeometry(geom, geomDict, fullCaseSetupDict):
    """Resolve FAN defaults after geometry has been linked into triSurface."""
    fanName = stripExt(geom)
    setup = fullCaseSetupDict[fanName]
    try:
        vertices, faces = readGeomFile(geom)
        surface = calculate_planar_surface_geometry(vertices, faces)
    except Exception as error:
        sys.exit('ERROR! Unable to resolve planar FAN geometry %s: %s' % (fanName, error))

    centerSetting = setup['FAN_CENTER']
    if len(centerSetting) == 1 and centerSetting[0].lower() == 'default':
        center = surface['center']
    else:
        center = np.asarray([float(value) for value in centerSetting], dtype=float)

    directionSetting = setup['FAN_DIRECTION']
    if len(directionSetting) == 1 and directionSetting[0].lower() == 'default':
        targetName = setup['FAN_TARGET_GEOMETRY'][0]
        targetGeom = next(
            (candidate for candidate in geomDict
             if stripExt(candidate).lower() == targetName.lower()),
            None
        )
        if targetGeom is None:
            sys.exit('ERROR! %s -> FAN_TARGET_GEOMETRY %s was not found in GEOMETRY.' %
                     (fanName, targetName))
        try:
            targetVertices, targetFaces = readGeomFile(targetGeom)
            targetCenter = calculate_surface_centroid(targetVertices, targetFaces)
            direction = orient_surface_normal(surface['normal'], center, targetCenter)
        except Exception as error:
            sys.exit('ERROR! Unable to orient %s FAN_DIRECTION toward %s: %s' %
                     (fanName, targetName, error))
    else:
        direction = np.asarray([float(value) for value in directionSetting], dtype=float)
        norm = np.linalg.norm(direction)
        if norm <= np.finfo(float).eps:
            sys.exit('ERROR! %s -> FAN_DIRECTION cannot be a zero vector.' % fanName)
        direction = direction / norm

    diameterSetting = setup['FAN_DIAMETER']
    if len(diameterSetting) == 1 and diameterSetting[0].lower() == 'default':
        diameter = surface['diameter']
    else:
        diameter = float(diameterSetting[0])
    if diameter <= 0.0:
        sys.exit('ERROR! %s -> FAN_DIAMETER must be positive.' % fanName)

    setup['FAN_CENTER'] = ['%.12g' % value for value in center]
    setup['FAN_DIRECTION'] = ['%.12g' % value for value in direction]
    setup['FAN_DIAMETER'] = ['%.12g' % diameter]
    setup['FAN_AREA'] = ['%.12g' % (np.pi * (diameter / 2.0) ** 2)]
    setup['FAN_PLANARITY_ERROR'] = ['%.12g' % surface['planarity_error']]
    print('\t\tResolved %s center=%s direction=%s diameter=%s' %
          (fanName, setup['FAN_CENTER'], setup['FAN_DIRECTION'], setup['FAN_DIAMETER'][0]))


#writes to snappyHexMeshDict
def writeSnappy(geomDict,fullCaseSetupDict):
    print('\n\tPreparing geometry for writing to snappyHexMeshDict...')
    snappyTemplate = fullCaseSetupDict['GLOBAL_REFINEMENT']['TEMPLATE_TYPE'][0]
    
    snappyTemplatePath = '%s/defaultDicts/system/%s' % (templateLoc,snappyTemplate)
    snappyRefinementPath = '%s/defaultDicts/system/snappyRefinementDict' % (templateLoc)
    snappyQualityPath = '%s/defaultDicts/system/meshQualityDict' % (templateLoc)
    #check if case already has system folder
    if not os.path.exists('%s/%s/system' % (path,case)):
        print('\t\tsystem directory not found, creating!')
        os.system('mkdir %s/%s/system' % (path,case))
        print('\t\tUsing snappyHexMesh template: %s' % (snappyTemplate))
        if 'ansa' in snappyTemplate.lower():
            print('\t\tUsing ansaMesh!')
        elif 'fidelity' in snappyTemplate.lower():
            print('\t\tUsing fidelityMesh!')
        else:
            print('\t\tCopying snappyHexMesh template from: %s' % (snappyTemplatePath))
            os.system('cp %s %s/%s/system/snappyHexMeshDict' % (snappyTemplatePath,path,case))
    else:
        print('\t\tUsing snappyHexMesh template: %s' % (snappyTemplate))
        if 'ansa' in snappyTemplate.lower():
            print('\t\tUsing ansaMesh!')
        elif 'fidelity' in snappyTemplate.lower():
            print('\t\tUsing fidelityMesh!')
        else:
            print('\t\tCopying snappyHexMesh template from: %s' % (snappyTemplatePath))
            os.system('cp %s %s/%s/system/snappyHexMeshDict' % (snappyTemplatePath,path,case))
    
    os.system('cp %s %s/%s/system/snappyRefinementDict' % (snappyRefinementPath,path,case))
    os.system('cp %s %s/%s/system/meshQualityDict' % (snappyQualityPath,path,case))
    ############
    #if using ansa, doesn't write out snappy things!
    ############
    if 'ansa' in snappyTemplate.lower():
        return geomDict,fullCaseSetupDict
    elif 'fidelity' in snappyTemplate.lower(): 
        return geomDict,fullCaseSetupDict
    if not os.path.exists(snappyTemplatePath):
        print('ERROR! TEMPLATE_TYPE in path %s is invalid!' % (snappyTemplatePath))
        exit()   
    snappyDictSections = ['GEOMETRY','FEATURE_EDGE','REFINEMENT_SURFACES','REFINEMENT_REGIONS','FACE_ZONE_CONTROLS','LOC_IN_MESH','LAYERS','LOC_IN_MESH','REF_ANGLE','DEF_EX_RATIO']
    snappyDict = {}
    
    #making sections in snappyDict
    for sec in snappyDictSections:
        snappyDict[sec] = []
        
    #setting up the geometries
    geometryStrings = {'GEOM':'\tGEOM_NAME {type distributedTriSurfaceMesh; scale GEOM_SCALE; file "GEOM_FILE"; regions{".*";}}\n',
                       'FAN':'\tGEOM_NAME {type distributedTriSurfaceMesh; scale GEOM_SCALE; file "GEOM_FILE"; regions{".*";}}\n',
                       'POR':'\tGEOM_NAME {type distributedTriSurfaceMesh; scale GEOM_SCALE; file "GEOM_FILE"; regions{".*";}}\n',
                       'REF':'\tGEOM_NAME {type triSurfaceMesh; scale GEOM_SCALE; file "GEOM_FILE"; regions{".*";}}\n',
                       'REFX':'\tGEOM_NAME {type triSurfaceMesh; scale GEOM_SCALE; file "GEOM_FILE"; regions{".*";}}\n',
                       'MRFG':'\tGEOM_NAME {type triSurfaceMesh; scale GEOM_SCALE; file "GEOM_FILE"; regions{".*";}}\n',
                       }
    refinementRegionStrings = {'GEOM':'\t"GEOM_NAME.*" {mode distance; levels (REF_LEVEL);}\n',
                               'MRFG':'\t"GEOM_NAME.*" {mode distance; levels (REF_LEVEL);}\n',
                               'REFX':{'inside':'\tGEOM_NAME {mode inside; levels ((1E15 REF_LEVEL));}\n',
                                       'outside':'\tGEOM_NAME {mode outside; levels ((1E15 REF_LEVEL));}\n',
                                       'distance':'\tGEOM_NAME {mode distance; levels (REF_LEVELS);}\n'
                                       },
                               'REF':'\tGEOM_NAME {mode inside; levels ((1E15 REF_LEVEL));}\n',
                               }
    refinementSurfaceStrings = {'GEOM': '\tGEOM_NAME {level (GEOM_LEVEL GEOM_LEVEL); regions{#include"snappyRefinementDict"}}\n',
                                'GEOMX': '\tGEOM_NAME {level (GEOMX_LEVEL_MIN GEOMX_LEVEL_MAX); regions{#include"snappyRefinementDict"}}\n',
                                'POR': '\tGEOM_NAME {level (GEOM_LEVEL GEOM_LEVEL); faceZone GEOM_NAME; cellZone GEOM_NAME_INTERNAL; cellZoneInside insidePoint; insidePoint (POR_POINT);}\n',
                                'FAN': '\tGEOM_NAME{level (GEOM_LEVEL GEOM_LEVEL); faceZone GEOM_NAME; faceType baffle;}\n',
                                'MRFG': '\tGEOM_NAME {level (GEOM_LEVEL GEOM_LEVEL); faceZone toint-GEOM_NAME; cellZone fluid-GEOM_NAME; cellZoneInside insidePoint; insidePoint (MRF_POINT);}\n',
                                
                                }
    edgeRefinementStrings = {'GEOM': '\t{file "GEOM_NAME.eMesh"; scale GEOM_SCALE; level EDGE_LEVEL;}\n',
                             'MRFG': '\t{file "GEOM_NAME.eMesh"; scale GEOM_SCALE; level EDGE_LEVEL;}\n',
                            'POR': '',
                            'REF':"",
                            'REFX':""
                            }
    layerStrings = {'GEOM': '''\t"GEOM_NAME.*" {nSurfaceLayers NLAYERS; expansionRatio EXP_RATIO;}\n''',
                    'POR':'''\t"GEOM_NAME.*" {nSurfaceLayers -1; expansionRatio 1;}\n''',
                    'MRFG':'''\t"GEOM_NAME.*" {nSurfaceLayers -1; expansionRatio 1;}\n''',
                    'REF':"",
                    'REFX':""
                    }
    noSnappyGeom = ['SMP','GRND']
    
    for geom in geomDict.keys():
        if geom.startswith(tuple(noSnappyGeom)):
            continue
        geomName = stripExt(geom)
        geomFile = geom.replace('.gz','')
        geomLayers = geomDict[geom]['layers']
        geomLayerExp = geomDict[geom]['expansion']
        try:
            geomScale = float(geomDict[geom]['scale'])
        except:
            print('ERROR! Geometry scale value is invalid for %s, must be float or int!' % (geomName))
            exit()
        
        #writing out geometry to geometry section
        if geomName.split('-')[0] in geometryStrings.keys():
            geomString = geometryStrings[geomName.split('-')[0]]
            geomString = geomString.replace('GEOM_NAME',geomName).replace('GEOM_FILE',geomFile).replace('GEOM_SCALE',str(geomScale))
        else:
            geomString = geometryStrings['GEOM']
            geomString = geomString.replace('GEOM_NAME',geomName).replace('GEOM_FILE',geomFile).replace('GEOM_SCALE',str(geomScale))
        
        snappyDict['GEOMETRY'].append(geomString)
        
        
        #writing out layer to layers section
        if geomName.split('-')[0] in layerStrings.keys():
            layerString = layerStrings[geomName.split('-')[0]]
            layerString = layerString.replace('GEOM_NAME',geomName)\
                                     .replace('GEOM_FILE',geomFile)\
                                     .replace('NLAYERS',str(geomLayers))\
                                     .replace('EXP_RATIO',str(geomLayerExp))
        elif geomName.startswith('GEOMX'):
            layerString = layerStrings['GEOM']
            layerString = layerString.replace('GEOM_NAME',geomName)\
                                     .replace('GEOM_FILE',geomFile)\
                                     .replace('NLAYERS',str(fullCaseSetupDict[geomName]['GEOMX_NLAYERS'][0]))\
                                     .replace('EXP_RATIO',str(fullCaseSetupDict[geomName]['GEOMX_LAYER_RATIO'][0]))
        else:
            layerString = layerStrings['GEOM']
            layerString = layerString.replace('GEOM_NAME',geomName)\
                                     .replace('GEOM_FILE',geomFile)\
                                     .replace('NLAYERS',str(geomLayers))\
                                     .replace('EXP_RATIO',str(geomLayerExp))

        snappyDict['LAYERS'].append(layerString)
        
        #writing out refinement surfaces to refinement surface section
        geomLevel = geomDict[geom]['refinement'] 
        refType = ''
        refLevel = ''
        if geomName.split('-')[0] == 'FAN':
            resolveFanGeometry(geom, geomDict, fullCaseSetupDict)
            fanString = refinementSurfaceStrings['FAN']
            fanString = fanString.replace('GEOM_NAME', geomName).replace('GEOM_LEVEL', geomLevel)
            snappyDict['REFINEMENT_SURFACES'].append(fanString)
            snappyDict['FACE_ZONE_CONTROLS'].append(
                '\t%s { type faceZone; faceZone %s; faceType baffle; }\n' %
                (geomName, geomName)
            )
            continue
        elif geomName.split('-')[0] == 'POR':
            porString = refinementSurfaceStrings[geomName.split('-')[0]]
            #test if it's porous media
            try:
                porPoint = fullCaseSetupDict[geomName]['POINT'][0]
                if porPoint.lower() != 'default':
                    porPoint = (fullCaseSetupDict[geomName]['POINT'])
                else:
                    bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ = getBoundingBox(geom.replace('.gz',''))
                    xcenter,ycenter,zcenter,radius = getRotaCoordinates(bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ)
                    porPoint = (str(xcenter),str(ycenter),str(zcenter))
                for coord in porPoint:
                    try:
                        coord = float(coord)
                    except:
                        print('\nERROR! Inside point for %s is not valid!' % (geomName))
                        return
                        sys.exit() 
                porPoint = " ".join(porPoint)
                print('\t\t\tSetting inside point for %s: %s' % (geomName,porPoint))
                porPIDs = getGeomPID("constant/triSurface/%s" % (geom.replace('.gz','')))       
            except Exception as E:
                print(E)
                porPoint = ''
                sys.exit('ERROR! Inside point for %s not valid!' % (geom.replace('.gz','')))

            # #test if refinement geometry
            porString = porString.replace('GEOM_NAME',geomName)\
                                            .replace('GEOM_LEVEL',geomLevel)\
                                            .replace('POR_POINT',str(porPoint))\
                                            .replace('REF_TYPE',str(refType))\
                                            .replace('REF_LEVEL',str(refLevel))   
            snappyDict['REFINEMENT_SURFACES'].append(porString)
            continue
        elif geomName.split('-')[0] == 'MRFG':
            mrfString = refinementSurfaceStrings[geomName.split('-')[0]]
            #test if it's mrf region
            try:
                mrfPoint = fullCaseSetupDict[geomName]['POINT'][0]
                if mrfPoint.lower() != 'default':
                    mrfPoint = (fullCaseSetupDict[geomName]['POINT'])
                else:
                    bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ = getBoundingBox(geom.replace('.gz',''))
                    xcenter,ycenter,zcenter,radius = getRotaCoordinates(bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ)
                    mrfPoint = (str(xcenter),str(ycenter),str(zcenter))
                for coord in mrfPoint:
                    try:
                        coord = float(coord)
                    except:
                        print('\nERROR! Inside point for %s is not valid!' % (geomName))
                        sys.exit() 
                mrfPoint = " ".join(mrfPoint)
                print('\t\t\tSetting inside point for %s: %s' % (geomName,mrfPoint))
                mrfPIDs = getGeomPID("constant/triSurface/%s" % (geom.replace('.gz','')))       
            except Exception as E:
                print(E)
                mrfPoint = ''
                sys.exit('ERROR! Inside point for %s not valid!' % (geom.replace('.gz','')))

           
            mrfString = mrfString.replace('GEOM_NAME',geomName)\
                                            .replace('GEOM_LEVEL',geomLevel)\
                                            .replace('MRF_POINT',str(mrfPoint))\
                                            .replace('REF_TYPE',str(refType))\
                                            .replace('REF_LEVEL',str(refLevel))
            
            
            distanceStringArray = []
            if fullCaseSetupDict[geomName]['MRF_REF_LEVELS'][0].lower() == 'default' or fullCaseSetupDict[geomName]['MRF_REF_DIST'][0].lower() == 'default':
                mrfLevels = fullCaseSetupDict['GLOBAL_REFINEMENT']['GEOM_REF_LEVELS']
                mrfDists = fullCaseSetupDict['GLOBAL_REFINEMENT']['GEOM_REF_DIST']
            else:
                mrfLevels = fullCaseSetupDict[geomName]['GEOM_REF_LEVELS']
                mrfDists = fullCaseSetupDict[geomName]['GEOM_REF_DIST']
            for dist,level in zip(mrfDists,mrfLevels):
                distanceString = '(%s %s)' % (dist, level)
                distanceStringArray.append(distanceString)
            MRFRefLevels = ' '.join(distanceStringArray)
            mrfRefString = refinementRegionStrings['MRFG']
            mrfRefString = refinementRegionStrings['MRFG'].replace('GEOM_NAME',geomName)\
                                                            .replace('REF_LEVEL',MRFRefLevels)
            snappyDict['REFINEMENT_REGIONS'].append(mrfRefString)
            snappyDict['REFINEMENT_SURFACES'].append(mrfString)
            
        elif geomName.split('-')[0] in ['REFX','REF']:
            regRefString = refinementRegionStrings[geomName.split('-')[0]]
            
            #test what refinement
            refType = ''
            refLevel = ''
            refLevels = ''
            
            if geomName.startswith('REFX'):
                refType = fullCaseSetupDict[geomName]['REF_TYPE'][0]
                if refType not in refinementRegionStrings['REFX'].keys():
                    print('ERROR! Invalid refinement type for %s' % (geomName))
                    print('Possible refinement types are: %s' % (list(refinementRegionStrings['REFX'].keys())))
                    sys.exit()
                
                if refType == 'distance':
                    distanceStringArray = []
                    refLevels = fullCaseSetupDict[geomName]['REF_LEVELS']
                    refDists = fullCaseSetupDict[geomName]['REF_DIST']
                    regRefString = refinementRegionStrings['REFX']['distance']
                    for dist,level in zip(refDists,refLevels):
                        distanceString = '(%s %s)' % (dist, level)
                        distanceStringArray.append(distanceString)
                    refLevels = ' '.join(distanceStringArray)
                    
                    regRefString = regRefString.replace('GEOM_NAME',geomName)\
                                               .replace('GEOM_LEVEL',geomLevel)\
                                               .replace('REF_TYPE',str(refType))\
                                               .replace('REF_LEVELS',str(refLevels))\
                                               .replace('REF_LEVELS',str(refLevels))
                    snappyDict['REFINEMENT_REGIONS'].append(regRefString)
                    continue
                    
            elif geomName.startswith('REF'):
                refLevel = fullCaseSetupDict[geomName]['REF_LEVEL'][0]
                try:
                    refLevel = int(refLevel)
                    #continue
                except:
                    print('ERROR! Invalid refinement values for %s' % (geomName))
                    
                    #return
                    sys.exit()


            #test if refinement geometry
            regRefString = regRefString.replace('GEOM_NAME',geomName)\
                                            .replace('GEOM_LEVEL',geomLevel)\
                                            .replace('REF_TYPE',str(refType))\
                                            .replace('REF_LEVEL',str(refLevel))                     
        else:

            if geomName.startswith('GEOMX'):
                
                geomRefString = refinementSurfaceStrings['GEOMX']
                geomLevelMin = fullCaseSetupDict[geomName]['GEOMX_MIN_MAX_LEVEL'][0]
                geomLevelMax = fullCaseSetupDict[geomName]['GEOMX_MIN_MAX_LEVEL'][1]
                geomRefString = geomRefString.replace('GEOM_NAME',geomName)\
                                             .replace('GEOMX_LEVEL_MIN',geomLevelMin)\
                                             .replace('GEOMX_LEVEL_MAX',geomLevelMax)
                distanceStringArray = []
                refLevels = fullCaseSetupDict[geomName]['GEOMX_REF_DIST']
                refDists = fullCaseSetupDict[geomName]['GEOMX_REF_LEVELS']
                for dist,level in zip(refDists,refLevels):
                    distanceString = '(%s %s)' % (dist, level)
                    distanceStringArray.append(distanceString)
                refLevels = ' '.join(distanceStringArray)
                regRefString = refinementRegionStrings['GEOM']
                regRefString = refinementRegionStrings['GEOM'].replace('GEOM_NAME',geomName)\
                                                              .replace('REF_LEVEL',refLevels)
            
            else:
                geomRefString = refinementSurfaceStrings['GEOM']
                geomRefString = geomRefString.replace('GEOM_NAME',geomName)\
                                            .replace('GEOM_LEVEL',geomLevel)
                distanceStringArray = []
                refLevels = fullCaseSetupDict['GLOBAL_REFINEMENT']['GEOM_REF_LEVELS']
                refDists = fullCaseSetupDict['GLOBAL_REFINEMENT']['GEOM_REF_DIST']
                for dist,level in zip(refDists,refLevels):
                    distanceString = '(%s %s)' % (dist, level)
                    distanceStringArray.append(distanceString)
                refLevels = ' '.join(distanceStringArray)
                regRefString = refinementRegionStrings['GEOM']
                regRefString = refinementRegionStrings['GEOM'].replace('GEOM_NAME',geomName)\
                                                              .replace('REF_LEVEL',refLevels)
        
        if geomRefString != None:
            snappyDict['REFINEMENT_SURFACES'].append(geomRefString)
        if regRefString != None:
            snappyDict['REFINEMENT_REGIONS'].append(regRefString)

        geomRefString = None
        regRefString = None
        
        
        #writing out feature edges to feature edges section
        edgeLevel = int(geomDict[geom]['refinement'][0]) + int(fullCaseSetupDict['GLOBAL_REFINEMENT']['FEAT_EDGE_LEVEL_INC'][0])
        if geomName.split('-')[0] in edgeRefinementStrings.keys():
            edgeString = edgeRefinementStrings[geomName.split('-')[0]]
            edgeString = edgeString.replace('GEOM_NAME',geomName)\
                                            .replace('GEOM_FILE',geomFile)\
                                            .replace('GEOM_SCALE',str(geomScale))\
                                            .replace('EDGE_LEVEL',str(edgeLevel))
                                            
        elif geomName.split('-')[0] == 'GEOMX':
            minRef = fullCaseSetupDict[geomName]['GEOMX_MIN_MAX_LEVEL'][0]
            refInc = fullCaseSetupDict[geomName]['GEOMX_FEAT_EDGE_LEVEL_INC'][0]
            edgeString = edgeRefinementStrings['GEOM']
            edgeString = edgeString.replace('GEOM_NAME',geomName)\
                                            .replace('GEOM_FILE',geomFile)\
                                            .replace('GEOM_SCALE',str(geomScale))\
                                            .replace('EDGE_LEVEL',str(int(int(minRef) + int(refInc))))
        else:
            edgeString = edgeRefinementStrings['GEOM']
            edgeString = edgeString.replace('GEOM_NAME',geomName)\
                                            .replace('GEOM_FILE',geomFile)\
                                            .replace('GEOM_SCALE',str(geomScale))\
                                            .replace('EDGE_LEVEL',str(edgeLevel))
        
        snappyDict['FEATURE_EDGE'].append(edgeString)
        
        #get location in mesh
        snappyDict['LOC_IN_MESH'] = fullCaseSetupDict['GLOBAL_REFINEMENT']['LOC_IN_MESH']
        snappyDict['LOC_IN_MESH'] = '(%s %s %s)' % (snappyDict['LOC_IN_MESH'][0], snappyDict['LOC_IN_MESH'][1], snappyDict['LOC_IN_MESH'][2])
        snappyDict['REF_ANGLE'] = fullCaseSetupDict['GLOBAL_REFINEMENT']['REF_ANGLE']
        snappyDict['DEF_EX_RATIO'] = fullCaseSetupDict['GLOBAL_REFINEMENT']['DEF_EX_RATIO']
        
    #check if default wake is selected, if true, get default
    try:
        wakeRefType = fullCaseSetupDict['GLOBAL_REFINEMENT']['DEFAULT_WAKE_REF'][0]
    except:
        print('ERROR! Entry in GLOBAL_REFINEMENT -> DEFAULT_WAKE_REF is not a valid string!')
        sys.exit('Please correct entry and re-run!')
    if wakeRefType.lower() == 'default' or wakeRefType == 'true':
        print('\t\tUsing default wake boxes!')
        wakeRefConfigPath = "%s/defaultRefinements/defaultWake" % (templateLoc)
    else:
        print('\t\tUsing %s wake boxes!' % (wakeRefType))
        wakeRefConfigPath = "%s/defaultRefinements/%s" % (templateLoc,wakeRefType)
    #check if wake config path exists
    if os.path.exists(wakeRefConfigPath):
        wakeRefConfig = configparser.ConfigParser()
        wakeRefConfig.optionxform = str
        try:
            wakeRefConfig.read_file(open(wakeRefConfigPath))
        except:
            print('ERROR! wake configuration is invalid!\nPath: %s' % (wakeRefConfigPath))
            sys.exit()
        wakeRefConfigSections = wakeRefConfig.sections()
        print('\n\t\tCreating wake boxes: %s' % (", ".join(wakeRefConfigSections)))
        
        #setting up wake values
        boxTypes = {'box':'searchableBox'}
        for section in wakeRefConfigSections:
            
            geomString = "\tWAKE_NAME{type BOX_TYPE; min (BOX_MIN); max (BOX_MAX);}\n"
            refString = "\tWAKE_NAME{mode inside; levels ((1E15 REF_LEVEL));}\n"
            wakeName = section
            boxType = boxTypes[wakeRefConfig[section]['TYPE']]
            boxMin = wakeRefConfig[section]['MIN']
            boxMax = wakeRefConfig[section]['MAX']
            refLevel = wakeRefConfig[section]['REF_LEVEL']
            
            geomString = geomString \
                        .replace('WAKE_NAME',wakeName)\
                        .replace('BOX_TYPE',boxType)\
                        .replace('BOX_MIN',boxMin)\
                        .replace('BOX_MAX',boxMax)
            
            refString = refString \
                        .replace('WAKE_NAME',wakeName)\
                        .replace('REF_LEVEL',refLevel)
                        
            snappyDict['GEOMETRY'].append(geomString)
            snappyDict['REFINEMENT_REGIONS'].append(refString)
    else:
        sys.exit('ERROR! %s is an invalid wake refinement configuration!\n\tMake sure %s exists!' % (wakeRefType,wakeRefConfigPath))
        
    #writing to snappyHexMeshDict
    print('\t\tWriting to snappyHexMeshDict:')
    search_and_replace("%s/%s/system/snappyHexMeshDict" % (path,case),
                       '<FACE_ZONE_CONTROLS>',
                       ''.join(snappyDict['FACE_ZONE_CONTROLS']))
    for snappySec in snappyDict.keys():
        print('\t\t\t%s' % (snappySec))
        search_and_replace("%s/%s/system/snappyHexMeshDict" % (path,case), '<%s>' % (snappySec),''.join(snappyDict[snappySec]))
    return geomDict,fullCaseSetupDict

def checkRefinements(geomDict,writeCaseSetupDict,fullCaseSetupDict):
    global updateCaseSetupFlag
    
    #load defaultRefinements
    
    refConfigRead = configparser.ConfigParser()
    refConfigRead.optionxform = str
    try:
        refConfigRead.read_file(open("%s/defaultBCTemplates/defaultRefinementRegions" % (templateLoc)))
    except:
        print('ERROR! defaultRefinementRegions template is invalid!')
        exit()
    
    refSections = refConfigRead.sections()
    
    for geom in geomDict:
        if geom.startswith('REF'):
            if geom.startswith('REFX'):
                refName = stripExt(geom)
                if refName in writeCaseSetupDict.keys():
                    refTypes = {'inside':'REF_SETUP','outside':'REF_SETUP','distance':'REF_SETUP_DISTANCE'}
                    #checking what type of refinement then making sure the required variables are in
                    try:
                        refType = writeCaseSetupDict[refName]['REF_TYPE'][0]                        
                        if refType in list(refTypes.keys()):                 
                            refTemplate = refTypes[refType]
                        else:
                            print('ERROR! Available REF_TYPE for [%s] are %s.' % (refName, list(refTypes.keys())))
                            sys.exit()
                    except: 
                        
                        print('WARNING: Invalid input for [%s]: REF_TYPE' % (refName))
                        refTemplate = 'REF_SETUP'
                        
                    if setupCore.fillAddonSection(writeCaseSetupDict, refName,
                                                  refConfigRead.items(refTemplate)):
                        updateCaseSetupFlag = True

                else:   
                    try:
                        writeCaseSetupDict[refName] = {}
                        for var in list(refConfigRead.items('REF_SETUP')):
                            var = list(var)
                            writeCaseSetupDict[refName][var[0]] = [var[1]]
                        updateCaseSetupFlag = True
                    except:
                        print('ERROR! defaultRefinementRegions section header is not [REF_SETUP]!')
                        exit()
            else:
                refName = stripExt(geom)
                refLevel = str(geomDict[geom]['refinement'])
                fullCaseSetupDict[refName] = {}
                for var in list(refConfigRead.items('REF_SETUP')):
                    var = list(var)
                    if var[0] == 'REF_LEVEL':
                        fullCaseSetupDict[refName][var[0]] = [refLevel]
                    else:
                        fullCaseSetupDict[refName][var[0]] = [var[1]]
            
            
    return writeCaseSetupDict,fullCaseSetupDict
                
def writeSurfaceFieldAverage(geomDict,fullCaseSetupDict):
    print('\n\tWriting surfaceFieldAverage...')
    
    sampleArray = []
    fieldAverageArray = []
    for geom in geomDict.keys():
        if geom.startswith('POR'):
            geomName = geom.split('.')[0]
            geomFile = geom.replace('.gz','')
            geomPath = "constant/triSurface/%s" % (geomFile)
            if fullCaseSetupDict[geomName]['SAMPLE_POR'][0].lower() == 'true':
                print('\t\tGetting PIDs for %s' % (geomName))
                pidArray = getGeomPID(geomPath)
                if len(pidArray) > 0:
                    for pid in pidArray:
                        
                        pidName = pid.split('/')[-1]
                        print('\t\t\tAdding %s to surfaceFieldAverage' % (pidName))
                        sampleArray.append(pidName)
                       

                else:
                    print('\t\tUnable to get PIDs for %s, skipping writing to it to surfaceFieldAverage!' % (geomFile))
        elif geom.startswith('SMP'):

            geomName = geom.split('.')[0]
            geomFile = geom.replace('.gz','')
            sampleArray.append(geomFile)
            
            print('\n\t\tAdding %s to surfaceFieldAverage' % (geomName))

    for surface in sampleArray:
        fieldAverageString = """%s{type surfaceFieldValue; surfaceFormat vtk; libs (fieldFunctionObjects);fields (U UMean p pMean);operation  areaNormalAverage;regionType sampledSurface;name sampledSurface; sampledSurfaceDict{type sampledTriSurfaceMesh; surface %s; source cells; interpolate true;} writeFields true;writeToFile true;writeControl timeStep;}\n""" % (surface.split('.')[0],surface)
        fieldAverageArray.append(fieldAverageString)
    
    fieldAverageArray = ''.join(fieldAverageArray)
    search_and_replace('system/surfaceFieldAverage','<POROUS_MEDIA>',fieldAverageArray)
        
def getGeometry(fullCaseSetupDict,writeCaseSetupDict):
    global updateCaseSetupFlag
    geomDict = {}
    geomColumnNames = ['scale','refinement','layers','expansion','wallmodel']
    #getting geometry list from geometry section
    try:
        geometryList = fullCaseSetupDict['GEOMETRY']['GEOM']
    except:
        print('ERROR! [GEOMETRY] section invalid!')
        exit()
    #breaking down geometry list into geometry dict
    geomDict = setupCore.geomToDict(geomDict,geometryList,geomColumnNames)
    
    for geom in geomDict.keys():
        print('\n\t\t\t%s' % (geom))
        for col in geomDict[geom].keys():
            val = geomDict[geom][col]
            print('\t\t\t%s: %s' % (col,val))
    #simple add-ons: fill missing keys from defaultBCTemplates/<template> [<section>]
    #FAN and REFX need extra validation and keep their own check functions
    for prefix, template, section in addonTemplates:
        if prefix == 'FAN':
            writeCaseSetupDict = checkFan(geomDict,writeCaseSetupDict)
            continue
        writeCaseSetupDict, flag = setupCore.checkAddonDefaults(
            geomDict, writeCaseSetupDict, templateLoc, prefix, template, section)
        updateCaseSetupFlag = updateCaseSetupFlag or flag
    writeCaseSetupDict,fullCaseSetupDict = checkRefinements(geomDict,writeCaseSetupDict,fullCaseSetupDict)

    return writeCaseSetupDict,geomDict,fullCaseSetupDict

def _fan_vector_setting(section, key, values, allow_default=True):
    """Validate a FAN three-component setting and return normalized tokens."""
    if len(values) == 1 and values[0].lower() == 'default':
        if allow_default:
            return ['default']
        sys.exit('ERROR! %s -> %s cannot be default.' % (section, key))
    if len(values) != 3:
        sys.exit('ERROR! %s -> %s must be default or exactly three numeric values.' %
                 (section, key))
    try:
        vector = np.asarray([float(value) for value in values], dtype=float)
    except (TypeError, ValueError):
        sys.exit('ERROR! %s -> %s must contain three numeric values.' % (section, key))
    if key == 'FAN_DIRECTION' and np.linalg.norm(vector) <= np.finfo(float).eps:
        sys.exit('ERROR! %s -> FAN_DIRECTION cannot be a zero vector.' % section)
    return [str(value) for value in vector]


def checkFan(geomDict,writeCaseSetupDict):
    """Populate and validate configuration for geometry-driven FAN components.

    FAN is intentionally independent of POR: a POR geometry retains its
    existing porous-media behavior, while a FAN geometry requests a planar
    baffle and the ESI fan pressure-jump model.
    """
    global updateCaseSetupFlag
    fanConfigRead = configparser.ConfigParser()
    fanConfigRead.optionxform = str
    try:
        fanConfigRead.read_file(open('%s/defaultBCTemplates/defaultFan' % templateLoc))
    except Exception:
        print('ERROR! defaultFan template is invalid!')
        sys.exit()
    if 'FAN_SETUP' not in fanConfigRead.sections():
        sys.exit('ERROR! defaultFan section header must be [FAN_SETUP]!')

    defaults = dict(fanConfigRead.items('FAN_SETUP'))
    for geom in geomDict:
        if not geom.startswith('FAN'):
            continue
        fanName = stripExt(geom)
        if fanName not in writeCaseSetupDict:
            writeCaseSetupDict[fanName] = {
                key: [value] for key, value in defaults.items()
            }
            updateCaseSetupFlag = True
        else:
            for key, value in defaults.items():
                if key not in writeCaseSetupDict[fanName]:
                    writeCaseSetupDict[fanName][key] = [value]
                    updateCaseSetupFlag = True

        setup = writeCaseSetupDict[fanName]
        model = setup['FAN_MODEL'][0].lower()
        if model not in ('fan', 'pressurejump'):
            sys.exit('ERROR! %s -> FAN_MODEL must be fan or pressureJump for ESI OpenFOAM.' % fanName)
        _fan_vector_setting(fanName, 'FAN_CENTER', setup['FAN_CENTER'])
        direction = _fan_vector_setting(fanName, 'FAN_DIRECTION', setup['FAN_DIRECTION'])
        if direction[0].lower() == 'default':
            target = setup.get('FAN_TARGET_GEOMETRY', ['default'])
            if len(target) != 1 or target[0].lower() == 'default':
                sys.exit('ERROR! %s -> FAN_TARGET_GEOMETRY is required when FAN_DIRECTION=default.' % fanName)
        elif 'FAN_TARGET_GEOMETRY' in setup:
            print('\t	%s: explicit FAN_DIRECTION overrides FAN_TARGET_GEOMETRY.' % fanName)
        diameter = setup['FAN_DIAMETER']
        if len(diameter) == 1 and diameter[0].lower() == 'default':
            pass
        elif len(diameter) == 1:
            try:
                if float(diameter[0]) <= 0.0:
                    raise ValueError
            except (TypeError, ValueError):
                sys.exit('ERROR! %s -> FAN_DIAMETER must be default or one positive scalar.' % fanName)
        else:
            sys.exit('ERROR! %s -> FAN_DIAMETER must be default or one positive scalar.' % fanName)
        try:
            pressureRise = setup['FAN_PRESSURE_RISE']
            if len(pressureRise) != 1 or float(pressureRise[0]) < 0.0:
                raise ValueError
        except (KeyError, TypeError, ValueError):
            sys.exit('ERROR! %s -> FAN_PRESSURE_RISE must be one non-negative scalar.' % fanName)

    return writeCaseSetupDict

def copyPvPostSetup():
    print('\tCopying pvPostSetup into case directory...')
    pvSetupPath = os.path.join(execDir,'postUtilities','default','pvPostSetup')
    destPath = os.path.join(path,case,'pvPostSetup')
    print('\t\t%s -> %s.' % (pvSetupPath,destPath))
    shutil.copyfile(pvSetupPath,destPath) 
    
    
  
           
        
if __name__ == '__main__':
    main()





