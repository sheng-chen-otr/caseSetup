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

#these lines get the path which the program is run from
path = os.path.split(os.getcwd())[0] #path of the case
case = os.path.split(os.getcwd())[1] #case name
jobPath = os.path.abspath(os.path.join(path,os.pardir)) #path of the job i.e. 100001

#location of the template file, should be in the same directory as _internal when compiled
templateBaseLoc = "%s/setupTemplates" % (os.path.dirname(os.path.realpath(__file__))) 

#checks the location of the templates to make sure they are good
if not os.path.isdir(templateBaseLoc):
    sys.exit("ERROR: Template path incorrect!")

parser = argparse.ArgumentParser(prog='caseSetup-v4.0',description='Set us the case based on settings written out in the caseSetup')
                    
parser.add_argument("-s","--setup", default='default', choices=['otr','otrwt','acewt'],
                    help='Identifies which setup templates to use.')
# parser.add_argument("-m","--mesh", action="store_true",
                    # help='Sets up meshing dicts and links geometry files. Will copy template case files if the trial folder is empty.')
parser.add_argument('-d',"--controlDict", action="store_true", 
                    help='Writes only controlDict for solver.')
parser.add_argument("--new", action="store_true", 
                    help='Writes a new caseSetup, will overwrite what is currently there if it exists!')
parser.add_argument("--modules", action="store_true", 
                    help='Shows all possible modules.')

args = parser.parse_args()
CONTROLDICT = args.controlDict
NEWCS = args.new
SETUP = args.setup
PR_MODULES = args.modules
updateCaseSetupFlag = False

#available addon keywords
addonKeyWords = ['POR','REFX','WAKE','GEOMX','ROTA','MOVG','IDOM']


#getting default values from template
def main():
    titleText = '''\t##############################\n\t######\tcaseSetup-v4.0\t######\n\t##############################'''
    print(titleText)
    getTemplateType(SETUP)
    
    
    
    defaultDict = getGlobalDefaults()
    if NEWCS == True:
        writeNewCaseSetup(defaultDict)
    else:
        caseSetupDict,writeCaseSetupDict,fullCaseSetupDict = getCaseSetup(defaultDict)
        writeCaseSetupDict,geomDict,fullCaseSetupDict = getGeometry(fullCaseSetupDict,writeCaseSetupDict) 
        
        cleanUpCaseSetup(geomDict,writeCaseSetupDict,fullCaseSetupDict,defaultDict)
        writeToCaseSetup(writeCaseSetupDict)
        linkGeomFiles(geomDict)
                
        geomDict ,fullCaseSetupDict = writeSnappy(geomDict,fullCaseSetupDict)
        writeSurfaceFeatureExtract(templateLoc,geomDict,fullCaseSetupDict)
        writeBlockMesh(templateLoc,fullCaseSetupDict)
        writeDecomposeParDict(templateLoc, fullCaseSetupDict)
        writeControlDict(templateLoc, fullCaseSetupDict)
        copyFunctionObjects(templateLoc,foList)
        writeSurfaceFieldAverage(geomDict,fullCaseSetupDict)
        writeBoundaries(templateLoc,geomDict,fullCaseSetupDict)
        
        writePostProSurfaceList(fullCaseSetupDict)
        writeSchemes(templateLoc,fullCaseSetupDict)
        writeSolution(templateLoc,fullCaseSetupDict)
        
        writeForceCoeff(geomDict,fullCaseSetupDict)
        writeOptions(templateLoc,geomDict,fullCaseSetupDict)
  

       
   
def getTemplateType(SETUP):
    global templateLoc
    templateLoc = "%s/%s" % (templateBaseLoc,SETUP)
    
    #checking if setup template exists
    if os.path.isdir(templateLoc):
        print('\n\tUsing %s setup template in path: %s' % (SETUP,templateLoc))
    else:
        print('ERROR! %s is not a valid setup template location!' % (templateLoc))
        exit()
   
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
        print('\t\tCopying snappyHexMesh template from: %s' % (snappyTemplatePath))
        os.system('cp %s %s/%s/system/snappyHexMeshDict' % (snappyTemplatePath,path,case))
    else:
        print('\t\tUsing snappyHexMesh template: %s' % (snappyTemplate))
        print('\t\tCopying snappyHexMesh template from: %s' % (snappyTemplatePath))
        os.system('cp %s %s/%s/system/snappyHexMeshDict' % (snappyTemplatePath,path,case))
    
    os.system('cp %s %s/%s/system/snappyRefinementDict' % (snappyRefinementPath,path,case))
    os.system('cp %s %s/%s/system/meshQualityDict' % (snappyQualityPath,path,case))
    if not os.path.exists(snappyTemplatePath):
        print('ERROR! TEMPLATE_TYPE in path %s is invalid!' % (snappyTemplatePath))
        exit()
    
    snappyDictSections = ['GEOMETRY','FEATURE_EDGE','REFINEMENT_SURFACES','REFINEMENT_REGIONS','LOC_IN_MESH','LAYERS','LOC_IN_MESH','REF_ANGLE','DEF_EX_RATIO']
    snappyDict = {}
    
    #making sections in snappyDict
    for sec in snappyDictSections:
        snappyDict[sec] = []
        
    #setting up the geometries
    geometryStrings = {'GEOM':'\tGEOM_NAME {type distributedTriSurfaceMesh; scale GEOM_SCALE; file "GEOM_FILE"; regions{".*";}}\n',
                       'POR':'\tGEOM_NAME {type distributedTriSurfaceMesh; scale GEOM_SCALE; file "GEOM_FILE"; regions{".*";}}\n',
                       'REF':'\tGEOM_NAME {type triSurfaceMesh; scale GEOM_SCALE; file "GEOM_FILE"; regions{".*";}}\n',
                       'REFX':'\tGEOM_NAME {type triSurfaceMesh; scale GEOM_SCALE; file "GEOM_FILE"; regions{".*";}}\n',
                       }
    refinementRegionStrings = {'GEOM':'\t"GEOM_NAME.*" {mode distance; levels (REF_LEVEL);}\n',
                               'REFX':{'inside':'\tGEOM_NAME {mode inside; levels ((1E15 REF_LEVEL));}\n',
                                       'outside':'\tGEOM_NAME {mode outside; levels ((1E15 REF_LEVEL));}\n',
                                       'distance':'\tGEOM_NAME {mode distance; levels (REF_LEVELS);}\n'
                                       },
                               'REF':'\tGEOM_NAME {mode inside; levels ((1E15 REF_LEVEL));}\n'}
    refinementSurfaceStrings = {'GEOM': '\tGEOM_NAME {level (GEOM_LEVEL GEOM_LEVEL); regions{#include"snappyRefinementDict"}}\n',
                                'GEOMX': '\tGEOM_NAME {level (GEOMX_LEVEL_MIN GEOMX_LEVEL_MAX); regions{#include"snappyRefinementDict"}}\n',
                                'POR': '\tGEOM_NAME {level (GEOM_LEVEL GEOM_LEVEL); faceZone GEOM_NAME; cellZone GEOM_NAME_INTERNAL; cellZoneInside insidePoint; insidePoint (POR_POINT);}\n'
                                
                                }
    edgeRefinementStrings = {'GEOM': '\t{file "GEOM_NAME.eMesh"; scale GEOM_SCALE; level EDGE_LEVEL;}\n',
                            'POR': '',
                            'REF':"",
                            'REFX':""
                            }
    layerStrings = {'GEOM': '''\t"GEOM_NAME.*" {nSurfaceLayers NLAYERS; expansionRatio EXP_RATIO;}\n''',
                    'POR':'''\t"GEOM_NAME.*" {nSurfaceLayers -1; expansionRatio 1;}\n''',
                    'REF':"",
                    'REFX':""
                    }
    noSnappyGeom = ['SMP']
    
    for geom in geomDict.keys():
        if geom.startswith(tuple(noSnappyGeom)):
            continue
        geomName = stripExt(geom)
        geomFile = geom.replace('.gz','')
        geomLayers = geomDict[geom]['layers'][0]
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
        if geomName.split('-')[0] == 'POR':
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
                                               .replace('POR_POINT',str(porPoint))\
                                               .replace('REF_TYPE',str(refType))\
                                               .replace('REF_LEVELS',str(refLevels))\
                                               .replace('REF_LEVELS',str(refLevels))
                    snappyDict['REFINEMENT_REGIONS'].append(regRefString)
                    continue
                    
            elif geomName.startswith('REF'):
                refLevel = fullCaseSetupDict[geomName]['REF_LEVEL'][0]
                #test if values are valid
                try:
                    refLevel = int(refLevel)
                except:
                    print('ERROR! Invalid refinement values for %s' % (geomName))
                    return
                    exit()

            #test if refinement geometry
            regRefString = regRefString.replace('GEOM_NAME',geomName)\
                                            .replace('GEOM_LEVEL',geomLevel)\
                                            .replace('POR_POINT',str(porPoint))\
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
        
        snappyDict['REFINEMENT_SURFACES'].append(geomRefString)
        snappyDict['REFINEMENT_REGIONS'].append(regRefString)
        
        
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
        
    #check if default wake is selected, if true, get defaults
    if fullCaseSetupDict['GLOBAL_REFINEMENT']['DEFAULT_WAKE_REF'][0].lower() == 'true':
        wakeRefConfig = configparser.ConfigParser()
        wakeRefConfig.optionxform = str
        try:
            wakeRefConfig.read_file(open("%s/defaultRefinements/defaultWake" % (templateLoc)))
        except:
            print('ERROR! defaultWake is invalid!')
            sys.exit()
        wakeRefConfigSections = wakeRefConfig.sections()
        print('\n\t\tCreating default wake boxes: %s' % (", ".join(wakeRefConfigSections)))
        
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
   
        
    
        
        
    #writing to snappyHexMeshDict
    print('\t\tWriting to snappyHexMeshDict:')
    for snappySec in snappyDict.keys():
        print('\t\t\t%s' % (snappySec))
        search_and_replace("%s/%s/system/snappyHexMeshDict" % (path,case), '<%s>' % (snappySec),''.join(snappyDict[snappySec]))
    return geomDict,fullCaseSetupDict

def linkGeomFiles(geomDict):
    #have to check if constant/triSurface folder exists
    if os.path.exists("%s/%s/constant" % (path,case)):
        if os.path.exists("%s/%s/constant/triSurface" % (path,case)):
            print('\n\t\ttriSurface directory OK!')
        else:
            print('\n\t\tWARNING: triSurface directory not found, creating!')
            try:
                os.system('mkdir %s/%s/constant/triSurface' % (path,case))
            except Exception as e:
                print('\n\t\tERROR: unable to make triSurface directory due to error:')
                print(e)
                
    else:
        print('\n\t\tWARNING: constant directory not found, creating constant and constant/triSurface!')
        try:
            os.system('mkdir %s/%s/constant' % (path,case))
        except Exception as e:
                print('\n\t\tERROR: unable to make constant directory due to error:')
                print(e)
        try:
            os.system('mkdir %s/%s/constant/triSurface' % (path,case))
        except Exception as e:
                print('\n\t\tERROR: unable to make triSurface due to error:')
                print(e)
        
    print('\n\t\tLinking geometry files:')
        
    for geom in geomDict.keys():
        if not os.path.isfile("%s/02_reference/MSH/%s" % (jobPath,geom)):
            sys.exit("ERROR! TriSurface: %s not found in MSH folder!" % (geom))
        elif os.path.isfile("constant/triSurface/%s" % (geom)):
            print("\t\t\t%s is already linked, skipping..." % (geom))
            continue
        else:
            print('\t\t\t%s' % (geom))
            cmd = "ln -s ../../../../02_reference/MSH/%s constant/triSurface/%s" % (geom,geom)
            os.system(cmd)
            
#removes any unnecessary sections
def cleanUpCaseSetup(geomDict,writeCaseSetupDict,fullCaseSetupDict,defaultDict):
    sectionsToDelete = []
    #porous and ref check
    geomNames = [stripExt(i) for i in geomDict.keys()]
    for key in writeCaseSetupDict.keys():
        if any(key.startswith(i) for i in addonKeyWords) and key in geomNames:
            continue
        elif any(key.startswith(i) for i in addonKeyWords) and key not in geomNames:
            if key not in sectionsToDelete:
                sectionsToDelete.append(key)
 
    for key in fullCaseSetupDict.keys():
        if any(key.startswith(i) for i in addonKeyWords) and key in geomNames:
            continue
        elif any(key.startswith(i) for i in addonKeyWords) and key not in geomNames:
            if key not in sectionsToDelete:
                sectionsToDelete.append(key)

    for sectionsDel in sectionsToDelete:
        del writeCaseSetupDict[sectionsDel]
        del fullCaseSetupDict[sectionsDel]

    return writeCaseSetupDict,fullCaseSetupDict

    

       
def writeNewCaseSetup(defaultDict):
    print('\tWriting new caseSetup file...')
    defaultDictKeys = defaultDict.keys()
    defaultModules = defaultDict['MODULE_LINK']
    defaultModuleNames = defaultDict['MODULE_LINK'].keys()
    writeCaseSetupDict = {}
    toWriteNew = ['DEFAULTS','TITLES','GEOMETRY']
    for module in defaultModules:
        
        sectionName = defaultModules[module][0]
        if sectionName in toWriteNew:
            writeCaseSetupDict[sectionName] = {}
            writeCaseSetupDict[sectionName] = defaultDict[sectionName]
    
    writeToCaseSetup(writeCaseSetupDict)
        
def writeToCaseSetup(writeCaseSetupDict):
    writeConfig = configparser.ConfigParser()
    writeConfig.optionxform = str
    for module in writeCaseSetupDict.keys():
        writeConfig.add_section(module)
        for key in writeCaseSetupDict[module].keys():
            #print(" ".join(list(writeCaseSetupDict[module][key])))
            writeConfig.set(module,key," ".join(list(writeCaseSetupDict[module][key])))
    
    try:
        #caseSetupConfig.read_file(open("%s/%s/caseSetup" % (path,case)))
        with open("%s/%s/caseSetup" % (path,case),'w') as caseSetupFile:
            writeConfig.write(caseSetupFile)
    except:
        print('\nERROR! Unable to write caseSetup!')
        
    if updateCaseSetupFlag == True:
        sys.exit('\nWARNING: caseSetup has been updated with default values, please check caseSetup and rerun.')
        
        
def getGlobalDefaults():
    defaultDict = {}
    print('\n\tGetting global templates...')
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read_file(open("%s/defaultSetup/defaultOrder" % (templateLoc)))
    configSections = config.sections()
    
    #get the order so that the when writing to the caseSetup file, orders will be right
    defaultOrder = np.array(config.items('DEFAULT_ORDER')[0])
    
    defaultOrder = defaultOrder[1].split(' ')
    if PR_MODULES == True:
        print('\tAVAILABLE MODULES:' % (defaultOrder))
        for module in defaultOrder:
            print('\t\t' + module)
        sys.exit()
    defaultDict['MODULE_LINK'] = {}
    for file in defaultOrder:
        print('\t\t%s' % (file))
        defaultConfigRead = configparser.ConfigParser()
        defaultConfigRead.optionxform = str
        try:
            defaultConfigRead.read_file(open("%s/defaultSetup/%s" % (templateLoc,file)))
        except:
            print('ERROR! Default %s cannot be opened!' % (file))
            exit()
        defaultSections = defaultConfigRead.sections()
        
        for section in defaultSections:
            defaultDict[section] = {}
            defaultDict['MODULE_LINK'][file] = [section]
            defaultLine = np.array(defaultConfigRead.items(section))
            defaultVars = defaultLine[:,0]
            for i in range(len(defaultVars)):
                defaultDict[section][defaultVars[i]] = defaultLine[i,1].split(' ')
    return defaultDict
    
#reading in caseSetup
def getCaseSetup(defaultDict):

    print('\n\tReading caseSetup from case directory...')
    caseSetupConfig = configparser.ConfigParser()
    caseSetupConfig.optionxform = str
    #tries to read the caseSetup file, if it can't find it, it throws an error
    try:
        caseSetupConfig.read_file(open("%s/%s/caseSetup" % (path,case)))
    except Exception as e:
        print('\n\n\nERROR! caseSetup invalid!')
        print(e)
        exit()
    caseSetupDict = {}
    caseSetupConfigSections = caseSetupConfig.sections()
    #getting the modules that are already there
    for section in caseSetupConfigSections:
            caseSetupLine = np.array(caseSetupConfig.items(section))
            if len(caseSetupLine[:]) < 1:
                updateCaseSetupFlag = True
            else:
                caseSetupDict[section] = {}
                caseSetupVars = caseSetupLine[:,0]
                for i in range(len(caseSetupVars)):
                    caseSetupDict[section][caseSetupVars[i]] = caseSetupLine[i,1].split(' ')
    
    writeCaseSetupDict,fullCaseSetupDict = getCaseSetupDefaultModules(caseSetupDict,defaultDict)
    
    return caseSetupDict,writeCaseSetupDict,fullCaseSetupDict

def getCaseSetupDefaultModules(caseSetupDict,defaultDict):  
    global updateCaseSetupFlag
    caseSetupKeys = caseSetupDict.keys()
    defaultDictKeys = defaultDict.keys()
    try:
        defaultSection = caseSetupDict['DEFAULTS']   
    except:
        print('\n\nERROR! caseSetup does not have [DEFAULTS] section! To generate a new caseSetup, use the --new flag.')
        exit()
    caseSetupModules = caseSetupDict['DEFAULTS']['DEFAULT_MODULES']
    defaultModules = defaultDict['MODULE_LINK']
    defaultModuleNames = defaultDict['MODULE_LINK'].keys()
    
    #check that modules are valid
    invalidModules = []
    for inputModule in caseSetupModules:
        if inputModule in defaultModuleNames:
            continue
        elif inputModule != '':
            invalidModules.append(inputModule)
    if len(invalidModules) > 0:
        
        print('\n\n\nERROR! The following modules are invalid: ')
        for invalidModule in invalidModules:
            print('\t' + invalidModule)
        exit()
        
    
    #init dict for casesetup to read
    fullCaseSetupDict = {}
    #init dict for casesetup to write back
    writeCaseSetupDict = {}
    
    #if loops for default modules if they are there or not
    for module in defaultModules:
        sectionName = defaultModules[module][0]
        if module in caseSetupModules and sectionName in caseSetupKeys:
            fullCaseSetupDict[sectionName] = {}
            fullCaseSetupDict[sectionName] = defaultDict[sectionName]    
        elif module in caseSetupModules and sectionName not in caseSetupKeys:
            fullCaseSetupDict[sectionName] = {}
            fullCaseSetupDict[sectionName] = defaultDict[sectionName]
        elif module not in caseSetupModules and sectionName in caseSetupKeys:
            fullCaseSetupDict[sectionName] = {}
            fullCaseSetupDict[sectionName] = caseSetupDict[sectionName]
            writeCaseSetupDict[sectionName] = {}
            writeCaseSetupDict[sectionName] = caseSetupDict[sectionName]
            
            caseSetupSectionVars = caseSetupDict[sectionName].keys()
            defaultSectionVars = defaultDict[sectionName].keys()
            for var in defaultSectionVars:
                if var in caseSetupSectionVars:
                    continue
                elif var not in caseSetupSectionVars:
                    writeCaseSetupDict[sectionName][var] = defaultDict[sectionName][var]
                    updateCaseSetupFlag = True
            
        elif module not in caseSetupModules and sectionName not in caseSetupKeys:
            fullCaseSetupDict[sectionName] = {}
            fullCaseSetupDict[sectionName] = defaultDict[sectionName]
            writeCaseSetupDict[sectionName] = {}
            writeCaseSetupDict[sectionName] = defaultDict[sectionName]
            updateCaseSetupFlag = True
        elif sectionName.startswith('POR'):
            writeCaseSetupDict[sectionName] = {}
            writeCaseSetupDict[sectionName] = caseSetupDict[sectionName]
        else:
            print('Some weird error happened!')
    
    
    #addon check sections
    
    for section in caseSetupKeys:
        if any(section.startswith(i) for i in addonKeyWords):
            writeCaseSetupDict[section] = {}
            writeCaseSetupDict[section] = caseSetupDict[section]
            fullCaseSetupDict[section] = {}
            fullCaseSetupDict[section] = caseSetupDict[section]
            
    return writeCaseSetupDict,fullCaseSetupDict
    
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
                    tempDict = {}
                    refTypes = {'inside':'REF_SETUP','outside':'REF_SETUP','distance':'REF_SETUP_DISTANCE'}
                    #checking what type of refinement then making sure the required variables are in
                    try:
                        refType = writeCaseSetupDict[refName]['REF_TYPE'][0]                        
                        if refType in list(refTypes.keys()):                 
                            refTemplate = refTypes[refType]
                        else:
                            print('ERROR! Available REF_TYPE for [%s] are %s.' % (refName, list(refTypes.keys())))
                            exit()
                    except: 
                        
                        print('WARNING: Invalid input for [%s]: REF_TYPE' % (refName))
                        refTemplate = 'REF_SETUP'
                        
                    for defaultRefVar in refConfigRead.items(refTemplate):
                        defaultVar = defaultRefVar[0]
                        defaultVal = defaultRefVar[1]
                        try:
                            tempDict[defaultVar] = writeCaseSetupDict[refName][defaultVar]      
                        except:
                            updateCaseSetupFlag = True
                            tempDict[defaultVar] = [str(defaultVal)]
                
                    writeCaseSetupDict[refName] = tempDict
                    
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
    geomDict = {}
    geomColumnNames = ['scale','refinement','layers','expansion','wallmodel']
    #getting geometry list from geometry section
    try:
        geometryList = fullCaseSetupDict['GEOMETRY']['GEOM']
    except:
        print('ERROR! [GEOMETRY] section invalid!')
        exit()
    #breaking down geometry list into geometry dict
    geomDict = geomToDict(geomDict,geometryList,geomColumnNames)
    
    for geom in geomDict.keys():
        print('\n\t\t\t%s' % (geom))
        for col in geomDict[geom].keys():
            val = geomDict[geom][col]
            print('\t\t\t%s: %s' % (col,val))
    writeCaseSetupDict = checkGeom(geomDict,writeCaseSetupDict)
    writeCaseSetupDict = checkPorous(geomDict,writeCaseSetupDict)
    writeCaseSetupDict = checkRotation(geomDict,writeCaseSetupDict)
    writeCaseSetupDict = checkInternalDomain(geomDict,writeCaseSetupDict)
    writeCaseSetupDict,fullCaseSetupDict = checkRefinements(geomDict,writeCaseSetupDict,fullCaseSetupDict)

    return writeCaseSetupDict,geomDict,fullCaseSetupDict




def checkPorous(geomDict,writeCaseSetupDict):
    global updateCaseSetupFlag
    #getting porous media default module
    porousConfigRead = configparser.ConfigParser()
    porousConfigRead.optionxform = str
    try:
        porousConfigRead.read_file(open("%s/defaultBCTemplates/defaultPorous" % (templateLoc)))
    except:
        print('ERROR! defaultPorous template is invalid!')
        exit()
    
    porousSections = porousConfigRead.sections()
    
    for geom in geomDict.keys():
        if geom.startswith('POR'):
            porousName = stripExt(geom)
            if porousName in writeCaseSetupDict.keys():               
                tempDict = {}
                for defaultPorousVar in porousConfigRead.items('POR_SETUP'):
                    defaultVar = defaultPorousVar[0]
                    defaultVal = defaultPorousVar[1]
                    try:
                        tempDict[defaultVar] = writeCaseSetupDict[porousName][defaultVar]
                    except:
                        updateCaseSetupFlag = True
                        tempDict[defaultVar] = [str(defaultVal)]
                
                writeCaseSetupDict[porousName] = tempDict
                    
            else:   
                try:
                    writeCaseSetupDict[porousName] = {}
                    for var in list(porousConfigRead.items('POR_SETUP')):
                        var = list(var)
                        writeCaseSetupDict[porousName][var[0]] = [var[1]]
                    updateCaseSetupFlag = True
                except:
                    print('ERROR! defaultPorous section header is not [POR_SETUP]!')
                    exit()
            
        
    return writeCaseSetupDict    

def checkGeom(geomDict,writeCaseSetupDict):
    global updateCaseSetupFlag
    #getting custom geom default template
    geomConfigRead = configparser.ConfigParser()
    geomConfigRead.optionxform = str
    try:
        geomConfigRead.read_file(open("%s/defaultBCTemplates/defaultGeomx" % (templateLoc)))
    except:
        print('ERROR! defaultGeomx template is invalid!')
        exit()
    
    porousSections = geomConfigRead.sections()
    
    for geom in geomDict.keys():
        if geom.startswith('GEOMX'):
            geomName = stripExt(geom)
            if geomName in writeCaseSetupDict.keys():               
                tempDict = {}
                for defaultGeomVar in geomConfigRead.items('GEOMX_SETUP'):
                    defaultVar = defaultGeomVar[0]
                    defaultVal = defaultGeomVar[1]
                    try:
                        tempDict[defaultVar] = writeCaseSetupDict[geomName][defaultVar]
                    except:
                        
                        updateCaseSetupFlag = True
                        tempDict[defaultVar] = [str(defaultVal)]
                
                writeCaseSetupDict[geomName] = tempDict
                    
            else:   
                try:
                    writeCaseSetupDict[geomName] = {}
                    for var in list(geomConfigRead.items('GEOMX_SETUP')):
                        var = list(var)
                        writeCaseSetupDict[geomName][var[0]] = [var[1]]
                    updateCaseSetupFlag = True
                except:
                    print('ERROR! defaultGeomx section header is not [GEOMX_SETUP]!')
                    exit()
            
        
    return writeCaseSetupDict        

def checkInternalDomain(geomDict,writeCaseSetupDict):
    global updateCaseSetupFlag
    #getting internal domain default module
    idomConfigRead = configparser.ConfigParser()
    idomConfigRead.optionxform = str
    try:
        idomConfigRead.read_file(open("%s/defaultBCTemplates/defaultInternalDomain" % (templateLoc)))
    except:
        print('ERROR! defaultInternalDomain template is invalid!')
        exit()
        
    idomconfig = idomConfigRead.sections()
    for geom in geomDict.keys():
        if geom.startswith('IDOM'):
            domName = stripExt(geom)
            
            if domName in writeCaseSetupDict.keys():               
                tempDict = {}
                for defaultDomVar in idomConfigRead.items('INTERNAL_DOMAIN'):
                    defaultVar = defaultDomVar[0]
                    defaultVal = defaultDomVar[1]
                    try:
                        tempDict[defaultVar] = writeCaseSetupDict[domName][defaultVar]
                    except:
                        updateCaseSetupFlag = True
                        tempDict[defaultVar] = [str(defaultVal)]
                
                writeCaseSetupDict[domName] = tempDict
                    
            else:   
                try:
                    writeCaseSetupDict[domName] = {}
                    for var in list(idomConfigRead.items('INTERNAL_DOMAIN')):
                        var = list(var)
                        writeCaseSetupDict[domName][var[0]] = [var[1]]
                    updateCaseSetupFlag = True
                except:
                    print('ERROR! defaultInternalDomain section header is not [INTERNAL_DOMAIN]!')
                    exit()
        
            
        
    return writeCaseSetupDict
def checkRotation(geomDict,writeCaseSetupDict):
    global updateCaseSetupFlag

    #getting wheel default module
    rotationConfigRead = configparser.ConfigParser()
    rotationConfigRead.optionxform = str
    try:
        rotationConfigRead.read_file(open("%s/defaultBCTemplates/defaultWheel" % (templateLoc)))
    except:
        print('ERROR! defaultWheel template is invalid!')
        exit()
    
    rotationConfig = rotationConfigRead.sections()
    for geom in geomDict.keys():
        if geom.startswith('ROTA'):
            rotName = stripExt(geom)
            
            if rotName in writeCaseSetupDict.keys():               
                tempDict = {}
                for defaultWheelVar in rotationConfigRead.items('WH_SETUP'):
                    defaultVar = defaultWheelVar[0]
                    defaultVal = defaultWheelVar[1]
                    try:
                        tempDict[defaultVar] = writeCaseSetupDict[rotName][defaultVar]
                    except:
                        updateCaseSetupFlag = True
                        tempDict[defaultVar] = [str(defaultVal)]
                
                writeCaseSetupDict[rotName] = tempDict
                    
            else:   
                try:
                    writeCaseSetupDict[rotName] = {}
                    for var in list(rotationConfigRead.items('WH_SETUP')):
                        var = list(var)
                        writeCaseSetupDict[rotName][var[0]] = [var[1]]
                    updateCaseSetupFlag = True
                except:
                    print('ERROR! defaultWheel section header is not [WH_SETUP]!')
                    exit()
            
        
    return writeCaseSetupDict    
    
def geomToDict(geomDict,geometryList,geomColumnNames):
    
    try:
        print('\t\tGetting geometry...')
        geometryList = geometryList[0].split('\n')
        nGeoms = len(geometryList)
        print('\t\t\tFound: %1.0f' % (nGeoms))
    except:
        print('ERROR! Geometry input invalid!')
        
    #splitting geoms into columns
    for geometry in geometryList:
        geometry = geometry.split(',')
        
        if len(geometry) != len(geomColumnNames) + 1:
            print('ERROR! Geometry input invalid! Insufficient options input for: %s' % (geometry[0]))
            exit()
        else:
            geometryName = geometry[0]
            if geometryName in geomDict.keys():
                print('ERROR! Duplicate geometry name found: %s' % (geometryName))
                exit()
            geomDict[geometryName] = {}
            n = 0
            while n < (len(geomColumnNames)):
                geomCol = geomColumnNames[n]
                geomDict[geometryName][geomCol] = geometry[n+1]
                n = n + 1
    
    
    return geomDict


  
           
        
main()   


