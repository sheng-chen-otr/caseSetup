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



def main():
    pass


def getInitPoints():
    pass


def importGeometry(geomDict):
    print('\n\n\t\tImporting geometry!')
    for geom in geomDict:
        geomPath = os.path.join('constant','triSurface',geom)
        if not os.path.exists(geomPath):
            sys.exit('ERROR! %s not in triSurface directory!' % (geom))
        print('\t\t\t' + geom)

        scale = geomDict[geom]['scale']
        
        if float(scale) == 1:
            print('\t\t\t\tScale from: meters')
            units = utils.UnitSystem(length='meter')
        elif float(scale) == 0.001:
            print('\t\t\t\tScale from: millimeters')
            units = utils.UnitSystem(length='millimeter')
        else:
            print('WARNING! Unkown scale unit for %s, using default (m)!' % (geom))

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


def getTemplateType():
    global templateLoc, ansaTemplateLoc
    fullCaseSetupDict = getFullCaseSetup()
    SETUP = fullCaseSetupDict['GLOBAL_REFINEMENT']['TEMPLATE_TYPE'][0].split('_')[0]
    templateLoc = "%s/%s" % (templateBaseLoc,SETUP)
    ansaTemplateLoc = os.path.join(templateLoc,
                                   'defaultANSA',
                                   fullCaseSetupDict['GLOBAL_REFINEMENT']['TEMPLATE_TYPE'][0].replace(SETUP+"_ansa_",''))
    
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

def getFullCaseSetup():
    '''
    same as getCaseSetup() but pulls in the existing fullCaseSetup output by caseSetup
    '''

    caseSetupConfig = configparser.ConfigParser()
    caseSetupConfig.optionxform = str
    #tries to read the caseSetup file, if it can't find it, it throws an error
    try:
        caseSetupConfig.read_file(open("%s/%s/fullCaseSetupDict" % (path,case)))
    except Exception as e:
        print('\n\n\nERROR! fullCaseSetup invalid, make sure to run caseSetup first!')
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
    return caseSetupDict
    
