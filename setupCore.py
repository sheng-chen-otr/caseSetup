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

"""Shared config / geometry plumbing for the setup tools (caseSetup.py,
internalSetup.py).

Nothing in here reads module globals or parses the command line: the template
location, case directory, config file name and add-on prefixes are passed in,
and functions that can add default values return an ``updateFlag`` instead of
setting a global. The caller decides what to do with the flag (normally pass
it to ``writeToCaseSetup``, which exits asking the user to check the file).
"""

import os
import sys
import configparser

import numpy as np

__version__ = "5.0.0-dev"

#config values are lists of strings: dict[section][key] = value.split(' ')


def stripExt(geom):
    #same as utilities.stripExt, kept here so setupCore has no heavy imports
    return geom.split('.')[0]


def getGlobalDefaults(templateLoc, printModules=False):
    """Read setupTemplates/<setup>/defaultSetup/defaultOrder and every default
    module it lists. Returns defaultDict with a MODULE_LINK entry mapping
    module file name -> [section name]."""
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


def getCaseSetup(defaultDict, caseDir, addonKeyWords, configName='caseSetup'):
    """Read <caseDir>/<configName> and merge it with the defaults.

    Returns (caseSetupDict, writeCaseSetupDict, fullCaseSetupDict, updateFlag).
    """
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


def getCaseSetupDefaultModules(caseSetupDict, defaultDict, addonKeyWords):
    """Build the dict used by the writers (fullCaseSetupDict) and the dict
    written back to the case config (writeCaseSetupDict).

    Modules listed in DEFAULTS/DEFAULT_MODULES use the template values and are
    hidden from the written file; other modules are exposed, with any missing
    keys filled from the defaults. Sections starting with one of
    addonKeyWords are passed through unchanged.

    Returns (writeCaseSetupDict, fullCaseSetupDict, updateFlag).
    """
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


def writeToCaseSetup(writeCaseSetupDict, caseDir, filename='caseSetup', updateFlag=False):
    """Write a config dict to <caseDir>/<filename>. If updateFlag is set,
    exit afterwards asking the user to check the updated file."""
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


def writeNewCaseSetup(defaultDict, caseDir, configName='caseSetup',
                      toWriteNew=('DEFAULTS','TITLES','GEOMETRY')):
    """Write a fresh config containing only the toWriteNew sections."""
    print('\tWriting new %s file...' % (configName))
    defaultModules = defaultDict['MODULE_LINK']
    writeCaseSetupDict = {}
    for module in defaultModules:
        sectionName = defaultModules[module][0]
        if sectionName in toWriteNew:
            writeCaseSetupDict[sectionName] = defaultDict[sectionName]

    writeToCaseSetup(writeCaseSetupDict, caseDir, configName)


def geomToDict(geomDict, geometryList, geomColumnNames):
    """Split the GEOM line (one geometry per line, comma separated columns)
    into geomDict[geometryName][column] = value."""
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


def readAddonTemplate(templateLoc, templateName, sectionHeader):
    """Read defaultBCTemplates/<templateName> and return the items of
    [sectionHeader] as a list of (key, value) tuples, in file order."""
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


def fillAddonSection(writeCaseSetupDict, sectionName, defaultItems):
    """Fill one add-on section from template defaults.

    An existing section is rebuilt in template order: user values are kept,
    missing keys get the default, and keys not in the template are dropped
    from the written file. A missing section is created from the defaults.
    Returns updateFlag (True if anything was added).
    """
    updateFlag = False
    if sectionName in writeCaseSetupDict.keys():
        tempDict = {}
        for defaultVar, defaultVal in defaultItems:
            try:
                tempDict[defaultVar] = writeCaseSetupDict[sectionName][defaultVar]
            except:
                updateFlag = True
                tempDict[defaultVar] = [str(defaultVal)]
        writeCaseSetupDict[sectionName] = tempDict
    else:
        writeCaseSetupDict[sectionName] = {}
        for defaultVar, defaultVal in defaultItems:
            writeCaseSetupDict[sectionName][defaultVar] = [defaultVal]
        updateFlag = True
    return updateFlag


def checkAddonDefaults(geomDict, writeCaseSetupDict, templateLoc, prefix, templateName, sectionHeader):
    """Generic add-on check: for every geometry whose name starts with prefix,
    make sure its [<geometry name>] section has every key of
    defaultBCTemplates/<templateName> [sectionHeader].

    Returns (writeCaseSetupDict, updateFlag).
    """
    defaultItems = readAddonTemplate(templateLoc, templateName, sectionHeader)
    updateFlag = False
    for geom in geomDict.keys():
        if geom.startswith(prefix):
            if fillAddonSection(writeCaseSetupDict, stripExt(geom), defaultItems):
                updateFlag = True
    return writeCaseSetupDict, updateFlag


def cleanUpCaseSetup(geomDict, writeCaseSetupDict, fullCaseSetupDict, addonKeyWords):
    """Remove add-on sections whose geometry is no longer in the GEOM list."""
    sectionsToDelete = []
    geomNames = [stripExt(i) for i in geomDict.keys()]
    for key in list(writeCaseSetupDict.keys()) + list(fullCaseSetupDict.keys()):
        if any(key.startswith(i) for i in addonKeyWords) and key not in geomNames:
            if key not in sectionsToDelete:
                sectionsToDelete.append(key)

    for sectionsDel in sectionsToDelete:
        del writeCaseSetupDict[sectionsDel]
        del fullCaseSetupDict[sectionsDel]

    return writeCaseSetupDict,fullCaseSetupDict


def linkGeomFiles(geomDict, caseDir, jobPath):
    """Symlink each geometry from <jobPath>/02_reference/MSH into
    constant/triSurface (relative to the current directory, which is the case)."""
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
