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

#shared config/geometry functions for caseSetup.py and internalSetup.py
#no module globals: templateLoc, caseDir, config name and addon prefixes are passed in,
#functions that add default values return updateFlag instead of setting a global

import os
import sys
import configparser

import numpy as np

__version__ = "5.0.0-dev"

#config values are lists of strings: dict[section][key] = value.split(' ')


def stripExt(geom):
    #same as utilities.stripExt, kept here so setupCore has no heavy imports
    return geom.split('.')[0]


#reads defaultSetup/defaultOrder and every default module in it
#MODULE_LINK maps module file name -> [section name]
def getGlobalDefaults(templateLoc, printModules=False):
    defaultDict = {}
    print('\n\tGetting global templates...')
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read_file(open("%s/defaultSetup/defaultOrder" % (templateLoc)))

    #get the order so that the when writing to the caseSetup file, orders will be right
    defaultOrder = np.array(config.items('DEFAULT_ORDER')[0])

    defaultOrder = defaultOrder[1].split(' ')
    if printModules == True:
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


#reads the case config file and merges it with the defaults
#returns caseSetupDict, writeCaseSetupDict, fullCaseSetupDict, updateFlag
def getCaseSetup(defaultDict, caseDir, addonKeyWords, configName='caseSetup'):
    print('\n\tReading %s from case directory...' % (configName))
    caseSetupConfig = configparser.ConfigParser()
    caseSetupConfig.optionxform = str
    #tries to read the caseSetup file, if it can't find it, it throws an error
    try:
        caseSetupConfig.read_file(open(os.path.join(caseDir, configName)))
    except Exception as e:
        print('\n\n\nERROR! %s invalid!' % (configName))
        print(e)
        sys.exit()
    caseSetupDict = {}
    caseSetupConfigSections = caseSetupConfig.sections()
    #getting the modules that are already there
    for section in caseSetupConfigSections:
            caseSetupLine = np.array(caseSetupConfig.items(section))
            if len(caseSetupLine[:]) < 1:
                #empty section: dropped. (pre-v5 code set a local flag here that
                #had no effect; kept as a no-op so aero output is unchanged)
                pass
            else:
                caseSetupDict[section] = {}
                caseSetupVars = caseSetupLine[:,0]
                for i in range(len(caseSetupVars)):
                    caseSetupDict[section][caseSetupVars[i]] = caseSetupLine[i,1].split(' ')

    writeCaseSetupDict,fullCaseSetupDict,updateFlag = getCaseSetupDefaultModules(
        caseSetupDict,defaultDict,addonKeyWords)

    return caseSetupDict,writeCaseSetupDict,fullCaseSetupDict,updateFlag


#modules in DEFAULT_MODULES use the template values and are hidden from the written file,
#other modules are written out with missing keys filled from the defaults
#addon sections (addonKeyWords) are passed through as is
def getCaseSetupDefaultModules(caseSetupDict, defaultDict, addonKeyWords):
    updateFlag = False
    caseSetupKeys = caseSetupDict.keys()
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
            fullCaseSetupDict[sectionName] = defaultDict[sectionName]
        elif module in caseSetupModules and sectionName not in caseSetupKeys:
            fullCaseSetupDict[sectionName] = defaultDict[sectionName]
        elif module not in caseSetupModules and sectionName in caseSetupKeys:
            fullCaseSetupDict[sectionName] = caseSetupDict[sectionName]
            writeCaseSetupDict[sectionName] = caseSetupDict[sectionName]

            caseSetupSectionVars = caseSetupDict[sectionName].keys()
            defaultSectionVars = defaultDict[sectionName].keys()
            for var in defaultSectionVars:
                if var in caseSetupSectionVars:
                    continue
                elif var not in caseSetupSectionVars:
                    writeCaseSetupDict[sectionName][var] = defaultDict[sectionName][var]
                    updateFlag = True

        elif module not in caseSetupModules and sectionName not in caseSetupKeys:
            fullCaseSetupDict[sectionName] = defaultDict[sectionName]
            writeCaseSetupDict[sectionName] = defaultDict[sectionName]
            updateFlag = True
        else:
            print('Some weird error happened!')

    #addon check sections
    for section in caseSetupKeys:
        if any(section.startswith(i) for i in addonKeyWords):
            writeCaseSetupDict[section] = caseSetupDict[section]
            fullCaseSetupDict[section] = caseSetupDict[section]

    return writeCaseSetupDict,fullCaseSetupDict,updateFlag


#writes a config dict to caseDir/filename, exits after writing if updateFlag is set
def writeToCaseSetup(writeCaseSetupDict, caseDir, filename='caseSetup', updateFlag=False):
    writeConfig = configparser.ConfigParser()
    writeConfig.optionxform = str
    for module in writeCaseSetupDict.keys():
        writeConfig.add_section(module)
        for key in writeCaseSetupDict[module].keys():
            writeConfig.set(module,key," ".join(list(writeCaseSetupDict[module][key])))

    try:
        with open(os.path.join(caseDir, filename),'w') as caseSetupFile:
            writeConfig.write(caseSetupFile)
    except:
        print('\nERROR! Unable to write %s!' % (filename))

    if updateFlag == True:
        sys.exit('\nWARNING: %s has been updated with default values, please check %s and rerun.'
                 % (filename, filename))


#writes a new config with only the toWriteNew sections
def writeNewCaseSetup(defaultDict, caseDir, configName='caseSetup',
                      toWriteNew=('DEFAULTS','TITLES','GEOMETRY')):
    print('\tWriting new %s file...' % (configName))
    defaultModules = defaultDict['MODULE_LINK']
    writeCaseSetupDict = {}
    for module in defaultModules:
        sectionName = defaultModules[module][0]
        if sectionName in toWriteNew:
            writeCaseSetupDict[sectionName] = defaultDict[sectionName]

    writeToCaseSetup(writeCaseSetupDict, caseDir, configName)


#splits the GEOM line into geomDict[geometryName][column]
def geomToDict(geomDict, geometryList, geomColumnNames):
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


#reads defaultBCTemplates/templateName and returns the [sectionHeader] items in file order
def readAddonTemplate(templateLoc, templateName, sectionHeader):
    configRead = configparser.ConfigParser()
    configRead.optionxform = str
    try:
        configRead.read_file(open("%s/defaultBCTemplates/%s" % (templateLoc, templateName)))
    except:
        print('ERROR! %s template is invalid!' % (templateName))
        exit()
    if not configRead.has_section(sectionHeader):
        print('ERROR! %s section header is not [%s]!' % (templateName, sectionHeader))
        exit()
    return list(configRead.items(sectionHeader))


#fills one addon section from the template defaults
#existing section: keeps user values, adds missing keys, drops keys not in the template (template order)
#missing section: written from the defaults
#keepExtra: keep keys that are not in the template (after the template keys), e.g. BC group lines
#returns updateFlag
def fillAddonSection(writeCaseSetupDict, sectionName, defaultItems, keepExtra=False):
    updateFlag = False
    if sectionName in writeCaseSetupDict.keys():
        tempDict = {}
        for defaultVar, defaultVal in defaultItems:
            try:
                tempDict[defaultVar] = writeCaseSetupDict[sectionName][defaultVar]
            except:
                updateFlag = True
                tempDict[defaultVar] = [str(defaultVal)]
        if keepExtra:
            for var in writeCaseSetupDict[sectionName]:
                if var not in tempDict:
                    tempDict[var] = writeCaseSetupDict[sectionName][var]
        writeCaseSetupDict[sectionName] = tempDict
    else:
        writeCaseSetupDict[sectionName] = {}
        for defaultVar, defaultVal in defaultItems:
            writeCaseSetupDict[sectionName][defaultVar] = [defaultVal]
        updateFlag = True
    return updateFlag


#generic addon check: every geometry starting with prefix gets a [geometry name] section
#with all keys from defaultBCTemplates/templateName [sectionHeader]
def checkAddonDefaults(geomDict, writeCaseSetupDict, templateLoc, prefix, templateName, sectionHeader,
                       keepExtra=False):
    defaultItems = readAddonTemplate(templateLoc, templateName, sectionHeader)
    updateFlag = False
    for geom in geomDict.keys():
        if geom.startswith(prefix):
            if fillAddonSection(writeCaseSetupDict, stripExt(geom), defaultItems, keepExtra):
                updateFlag = True
    return writeCaseSetupDict, updateFlag


#removes addon sections that no longer have a geometry in GEOM
#keepSections: addon sections kept without a geometry (e.g. internal MOVE sections with FROM)
def cleanUpCaseSetup(geomDict, writeCaseSetupDict, fullCaseSetupDict, addonKeyWords, keepSections=()):
    sectionsToDelete = []
    geomNames = [stripExt(i) for i in geomDict.keys()]
    for key in list(writeCaseSetupDict.keys()) + list(fullCaseSetupDict.keys()):
        if any(key.startswith(i) for i in addonKeyWords) and key not in geomNames and key not in keepSections:
            if key not in sectionsToDelete:
                sectionsToDelete.append(key)

    for sectionsDel in sectionsToDelete:
        del writeCaseSetupDict[sectionsDel]
        del fullCaseSetupDict[sectionsDel]

    return writeCaseSetupDict,fullCaseSetupDict


#links each geometry from jobPath/02_reference/MSH into constant/triSurface
def linkGeomFiles(geomDict, caseDir, jobPath):
    #have to check if constant/triSurface folder exists
    if os.path.exists("%s/constant" % (caseDir)):
        if os.path.exists("%s/constant/triSurface" % (caseDir)):
            print('\n\t\ttriSurface directory OK!')
        else:
            print('\n\t\tWARNING: triSurface directory not found, creating!')
            try:
                os.system('mkdir %s/constant/triSurface' % (caseDir))
            except Exception as e:
                print('\n\t\tERROR: unable to make triSurface directory due to error:')
                print(e)

    else:
        print('\n\t\tWARNING: constant directory not found, creating constant and constant/triSurface!')
        try:
            os.system('mkdir %s/constant' % (caseDir))
        except Exception as e:
                print('\n\t\tERROR: unable to make constant directory due to error:')
                print(e)
        try:
            os.system('mkdir %s/constant/triSurface' % (caseDir))
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
