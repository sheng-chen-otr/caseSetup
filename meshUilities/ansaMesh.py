import sys
import fileinput
import os
import getopt
import re
import glob
import argparse as argparse
import configparser
import numpy as np
import math
from pathlib import Path


import ansa
from ansa import base
from ansa import constants
from ansa import utils
from ansa import mesh

path = os.path.split(os.getcwd())[0] #path of the case
case = os.path.split(os.getcwd())[1] #case name
jobPath = os.path.abspath(os.path.join(path,os.pardir)) #path of the job i.e. 100001
rootDir = Path(__file__).parents[1]
templateBaseLoc = "%s/setupTemplates" % (rootDir) 

addonKeyWords = ['POR','REFX','WAKE','GEOMX','ROTA','MOVG','IDOM','MRFG']

parser = argparse.ArgumentParser(prog='ansaMesh-v1.0.0',description='ANSA Meshing process for caseSetup')
                    
parser.add_argument("-s","--setup", default='default', 
                    help='Identifies which setup templates to use.')
parser.add_argument("-r","--refinement", default='default', 
                    help='Identified regional refinemetn boxes templates')

# parser.add_argument('-d',"--controlDict", action="store_true", 
#                     help='Writes only controlDict for solver.')
# parser.add_argument("--new", action="store_true", 
#                     help='Writes a new caseSetup, will overwrite what is currently there if it exists!')
# parser.add_argument("--modules", action="store_true", 
#                     help='Shows all possible modules.')
# parser.add_argument("--postProDict", action="store_true", 
#                     help='Copies post-processing config file into case folder.')
args = parser.parse_args()
SETUP = args.setup

def main():
    getTemplateType(SETUP)
    defaultDict = getGlobalDefaults()
   
    caseSetupDict,writeCaseSetupDict,fullCaseSetupDict = getCaseSetup(defaultDict)
    for section in fullCaseSetupDict.keys():
        if 'GEOMETRY' in section:
            continue
        print('\t' + section)
        for var in fullCaseSetupDict[section].keys():
            print('\t\t' + var + ' = ' + str(fullCaseSetupDict[section][var]).replace('[','')
                                                                             .replace(']','')
                                                                             .replace(',','')
                                                                             .replace('\'',''))
    geomDict = geomToDict(fullCaseSetupDict)

    print('\t\tFound %s geometries!' % (len(geomDict.keys())))
    for geom in geomDict.keys():
        print('\t\t' + geom)
        for var in geomDict[geom].keys():
            print('\t\t\t%s: %s' % (var,geomDict[geom][var]))
    
    importGeometry(geomDict)
    importDomain()
    createSizeField(fullCaseSetupDict)
    
    # getAllParts(geomDict)

def saveAnsa(caseName):
    print('\t\tSaving ANSA file...')
    base.SaveAs('%s.ansa.gz' % (caseName))

def createSizeField(fullCaseSetupDict):
    print('\t\tCreating offset size fields...')
    sf = base.CreateEntity(constants.NASTRAN, "SIZE FIELD")
    baseSize = float(fullCaseSetupDict['BC_SETUP']['BASE_CELL_SIZE'][0])
    refinementLevels = fullCaseSetupDict['GLOBAL_REFINEMENT']['GEOM_REF_LEVELS']
    refinementDist = fullCaseSetupDict['GLOBAL_REFINEMENT']['GEOM_REF_DIST']
    pidList = getAllPID()
    n = 0
    for level,dist in zip(refinementLevels,refinementDist):
        size = calculateCellSize(int(level),baseSize)
        dist = float(dist)*1000
        sfname = 'level_%s' % (level)
        print('\t\t\t%s' % (sfname))
        global_offset = mesh.GetNewSizeFieldSurfaceOffsetRule(size_field=sf,
                                                              name=sfname,
                                                              max_surf_len = size,
                                                              max_vol_len = size,
                                                              offset = dist)
        
        for pid in pidList:
            pidEntity = base.GetEntity(constants.NASTRAN,"PSHELL",pid)
            pidName = base.GetEntity(constants.NASTRAN,"PSHELL",pid)._name
            boundaryList = ['-min','-max']
            if not any(x in pidName for x in boundaryList):
                print('\t\t\t\t%s' % (pidName))
                mesh.AddContentsToSizeFieldRule( entities = pidEntity , rule = global_offset)

        
        # sf_rules = mesh.GetRulesFromSizeField(sf)
        # mesh.AddFilterToSizeFieldRule(field = "PIDs", expression = "does not contain", value = '-min', match='all', rule = sf_rules[n])
        # mesh.AddFilterToSizeFieldRule(field = "PIDs", expression = "does not contain", value = '-max', match='all', rule = sf_rules[n])
        
        
        # n = n + 1
    # mesh.ApplySizeFieldFilters(sf) 
    sf_rules = mesh.GetRulesFromSizeField(sf)
    print('\t\tAdded offsets:')
    for rule in sf_rules:
        print('\t\t\t' + rule._name)
        print('\t\t\t\tNumber of PIDs: ',str(len(mesh.GetContentsFromSizeFieldRule(rule))))
        mesh.SaveSizeFieldRuleParams( rule = rule, mpar_file = '%s.ansa_mpar' % (rule._name))

    print('\t\t\tBuilding size field!')
    mesh.BuildSizeField(sf)
    saveAnsa(case)

    

   
        
    
        
    
def calculateCellSize(level,baseSize):
    return (baseSize*1000)/(2**level)





def importDomain():
    print('\t\tCreating blockMesh to import to ANSA...')
    os.system('blockMesh;foamToSurface constant/triSurface/domain.stl')
    base.InputStereoLithography(filename='constant/triSurface/domain.stl',
                                unit_system=utils.UnitSystem(length='meter'))
    

def getAllPID():

    deck = constants.NASTRAN
    pidList = []

    # for geom in geomDict.keys():

    #     part = base.GetPartFromName(geom)
    #     print('\t\tPart: %s, %s' % (part._name,part._id))
        
    pids = base.CollectEntities(constants.NASTRAN, None, "PSHELL")
    

    for pid in pids:
    
        pidList.append(pid._id)

    return pidList
        
        
    



    
           

    


def importGeometry(geomDict):
    print('\n\n\t\tImporting geometry!')
    for geom in geomDict:
        geomPath = os.path.join('constant','triSurface',geom)
        if not os.path.exists(geomPath):
            sys.exit('ERROR! %s not in triSurface directory!' % (geom))
        print('\t\t\t' + geom)

        scale = geomDict[geom]['scale']
        
        if float(scale) == 1:
            units = utils.UnitSystem(length='meter')
        elif float(scale) == 0.001:
            units = utils.UnitSystem(length='millimeter')

        if '.obj' in geom:
            base.InputWaveFront(filename = geomPath,
                                unit_system=units)
        elif '.stl.' in geom:
            base.InputStereoLithography(filename = geomPath,
                                        unit_system=units)
    








        
        

def geomToDict(fullCaseSetupDict):
    geomDict = {}
    geomColumnNames = ['scale','refinement','layers','expansion','wallmodel']
    #getting geometry list from geometry section
    try:
        geometryList = fullCaseSetupDict['GEOMETRY']['GEOM']
    except:
        print('ERROR! [GEOMETRY] section invalid!')
        exit()
    
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


def getTemplateType(SETUP):
    global templateLoc
    templateLoc = "%s/%s" % (templateBaseLoc,SETUP)
    
    #checking if setup template exists
    if os.path.isdir(templateLoc):
        print('\n\tUsing %s setup template in path: %s' % (SETUP,templateLoc))
    else:
        print('ERROR! %s is not a valid setup template location!' % (templateLoc))
        exit()


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
        sys.exit()
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

main()