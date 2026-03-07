'''
version 1.6

Edit log:
- added pid refinement based on key words
- changed self proximity refinement
- turned on orientation proximity, so it shouldn't refine if two faces facing opposite ways are too close
'''

import sys
import os
import argparse as argparse
import configparser
import numpy as np
import ast
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



def main():
    
    getTemplateType()
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
    geomDict = getParts(geomDict) 
    sf = createSizeField(fullCaseSetupDict,geomDict)
    createOctree2(fullCaseSetupDict,geomDict)
    exportAnsaMesh()
    
    
    # getAllParts(geomDict)

def layerCoverage(geomDict,fullCaseSetupDict):
        #calculate layer coverage on each PID

        print("\n\t\tCalculating layer coverage per PID...")

        #first set result of number of layers
        results = base.CollectEntities(ansa.constants.OPENFOAM, None, 'RESULT')
        for result in results:
            ret = base.GetEntityCardValues(ansa.constants.OPENFOAM, result, ['__id__', 'Name', 'Status'])
            if ret['Status'] == 'Built' and ret['Name'] == "Hextreme Octree Mesh - Number of layers":
                base.SetCurrentEntity(result)

        totalLayerCellsGrown = 0
        totalLayerCellsRequested = 0

        #collect mesh parts of Hextreme Octree

        octree = base.GetEntity(ansa.constants.NASTRAN, "HEXTREME OCTREE", 1)
        octreeAreas = mesh.GetAreasFromOctree(octree)
        for area in octreeAreas:
            parts = mesh.GetPartsFromOctreeArea(area = area)
            layersBin = {}
            for part in parts:
                partName = part._name.split('_')[0]
                pidName = part._name
                base.Or(part)
                #getting requested layers from part dict
                if 'GEOMX-' in partName:
                    requestedLayers = int(fullCaseSetupDict[partName]['GEOMX_NLAYERS'])
                else:
                    for geomPart in geomDict.keys():
                        if partName in geomPart:
                            requestedLayers = int(geomDict[geomPart]['layers'])
                            break

                
                #requestedLayers = int(geomDict[partName]['layers'])

                if requestedLayers == 0:
                    print('\t\t\tRequested layers is 0, skipping!')
                if requestedLayers != 0:
                    print("\t\t\tRequested prism layers for %s region: %s" % (pidName, str(requestedLayers)))
                    #create bins for reporting layer growth
                    for i in range(0, requestedLayers+1):
                        layersBin.update({i: 0})
                    #gather all visible shells
                    partShells = ansa.base.CollectEntities(constants.OPENFOAM, None, ['POLYGON', 'SHELL'], filter_visible=True)

                    #count of surface elements
                    faceElements = len(partShells)

                    #pull total number of layer results for each surface shell
                    partShellsResults = ansa.base.ShellResult(partShells)

                    layerCellsGrown = 0
                    #assume requested layers is equal to number of surface elements times requested number of layers
                    layerCellsRequested = requestedLayers*faceElements

                    for shellResult in partShellsResults:
                        #for each element, update bin for number of layer cells grown
                        layersBin.update({int(shellResult[1]): (layersBin[shellResult[1]]+1)})
                        #update total number of layer cells grown for the PID
                        layerCellsGrown+=shellResult[1]

                    #update total number of prism layer cells grown and requested
                    totalLayerCellsGrown+=layerCellsGrown
                    totalLayerCellsRequested+=layerCellsRequested

                    print("\n\t\t\t%s" % pidName)
                    print("\t\t\t\tLayer cells requested: %s" % str(layerCellsRequested))
                    print("\t\t\t\tLayer cells delivered: %s" % str(int(layerCellsGrown)))
                    if layerCellsRequested == 0:
                        print("\t\t\t\tZero prism layer cells grown.")
                    if layerCellsRequested != 0:
                        print("\t\t\t\tDelivery rate: %0.2f%%" % (100*layerCellsGrown/layerCellsRequested))
                        for line in layersBin:
                            print("\t\t\t\tFaces with %0.0f layer(s) grown: %0.0f       %0.2f%% of faces in PID" % (line, layersBin[line], 100*layersBin[line]/faceElements))
                


                
        return
        octreeParts = ansa.mesh.GetPartsFromOctreeArea(octree)
        for part in octreeParts:
            partName = ansa.base.GetEntityCardValues(ansa.constants.OPENFOAM, part, ['Name'])['Name']
            ansa.base.Or(part)

            #match partName with meshArea in SESSION to grab target number of prism layers
            for area in self.meshAreas:
                if partName == area.octreeFilterName:
                    targetLayers = int(area.meshParametersDict['nLayers'][0])
                    layersBin = {}

                    if targetLayers == 0:
                        logger.debug("\n%s" % partName)
                        logger.debug("     Zero requested prism layers. Skipping...")
                    if targetLayers != 0:
                        print("Requested prism layers for %s region: %s" % (partName, str(targetLayers)))
                        #show only the surface mesh that is part of the volume mesh
                        ansa.base.Not(area.ansaGroupEntity)

                        #create bins for reporting layer growth
                        for i in range(0, targetLayers+1):
                            layersBin.update({i: 0})

                        #gather all visible shells
                        partShells = ansa.base.CollectEntities(ansa.constants.OPENFOAM, None, ['POLYGON', 'SHELL'], filter_visible=True)

                        #count of surface elements
                        faceElements = len(partShells)

                        #pull total number of layer results for each surface shell
                        partShellsResults = ansa.base.ShellResult(partShells)

                        layerCellsGrown = 0
                        #assume requested layers is equal to number of surface elements times requested number of layers
                        layerCellsRequested = targetLayers*faceElements

                        for shellResult in partShellsResults:
                            #for each element, update bin for number of layer cells grown
                            layersBin.update({int(shellResult[1]): (layersBin[shellResult[1]]+1)})
                            #update total number of layer cells grown for the PID
                            layerCellsGrown+=shellResult[1]

                        #update total number of prism layer cells grown and requested
                        totalLayerCellsGrown+=layerCellsGrown
                        totalLayerCellsRequested+=layerCellsRequested

                        logger.debug("\n%s" % partName)
                        logger.debug("     Layer cells requested: %s" % str(layerCellsRequested))
                        logger.debug("     Layer cells delivered: %s" % str(int(layerCellsGrown)))
                        if layerCellsRequested == 0:
                            logger.debug("     Zero prism layer cells grown.")
                        if layerCellsRequested != 0:
                            logger.debug("     Delivery rate: %0.2f%%" % (100*layerCellsGrown/layerCellsRequested))
                            for line in layersBin:
                                logger.debug("     Faces with %0.0f layer(s) grown: %0.0f       %0.2f%% of faces in PID" % (line, layersBin[line], 100*layersBin[line]/faceElements))

        logger.debug("\nTotal layer cells requested: %s" % str(totalLayerCellsRequested))
        logger.debug("Total layer cells delivered: %s" % str(int(totalLayerCellsGrown)))
        logger.debug("Total layer cells requested: %0.2f%%\n" % (100*totalLayerCellsGrown/totalLayerCellsRequested))


def runQualityImprovement(method='fixQuality'):
    #runs either fixQuality or reconstruct
    #method = self.globalMeshDefaults['globalQualityCriteriaFix'][1]
    #show only volume mesh entity
    volumeEnt = ansa.base.CollectEntities(ansa.constants.OPENFOAM, None, 'VOLUME')[0]
    ansa.base.Or(volumeEnt)

    #run deck report pre quality fix and save to log directory
    #ansa.utils.DeckInfo(os.path.join(log_dir, "deckinfo_preQual.html"), 'VISIBLE', 'HTML')

    if method == 'fixQuality':
        print("\n\t\tImproving quality criteria by appying %s" % method)    
        ansa.mesh.FixQualSolids()
        ansa.utils.DeckInfo("deckinfo_fixQuality.html", 'VISIBLE', 'HTML')

    if method == 'reconstruct':
        print("\n\t\tImproving quality criteria by appying %s" % method)
        ansa.mesh.ReconstructSolids()
        ansa.utils.DeckInfo("deckinfo_fixQuality.html", 'VISIBLE', 'HTML')


def saveAnsa(caseName):
    print('\t\tSaving ANSA file...')
    base.SaveAs('%s.ansa.gz' % (caseName))

def createRefinementRegion(templateLoc,fullCaseSetupDict,geomDict):

    refinementName = fullCaseSetupDict['GLOBAL_REFINEMENT']['DEFAULT_WAKE_REF'][0]
    if refinementName.lower() == 'true':
        refinementName = 'defaultWake'
    refinementPath = os.path.join(templateLoc,'defaultRefinements',refinementName)
    #check if wake config path exists
    if os.path.exists(refinementPath):
        wakeRefConfig = configparser.ConfigParser()
        wakeRefConfig.optionxform = str
        try:
            wakeRefConfig.read_file(open(refinementPath))
        except:
            print('ERROR! wake configuration is invalid!\nPath: %s' % (refinementPath))
            sys.exit()
        wakeRefConfigSections = wakeRefConfig.sections()
        print('\n\t\tCreating wake boxes: %s' % (", ".join(wakeRefConfigSections)))
    else:
        print('ERROR! wake configuration is invalid!\nPath: %s' % (refinementPath))
        sys.exit()

    wakeDict = {}
    for wake in wakeRefConfig.sections():
        print('\t\t\t%s' %(wake))
        min = np.array(wakeRefConfig[wake]['MIN'].split(' '),dtype=float)
        max = np.array(wakeRefConfig[wake]['MAX'].split(' '),dtype=float)
        print('\t\t\t\tMin: %s'% (min))
        print('\t\t\t\tMax: %s'% (max))
        size = int(wakeRefConfig[wake]['REF_LEVEL'][0])
        baseSize = float(fullCaseSetupDict['BC_SETUP']['BASE_CELL_SIZE'][0])
        wakeDict[wake] = base.SizeBoxMinMax(None, min*1000, max*1000,
                                        calculateCellSize(size,baseSize=baseSize),
                                        calculateCellSize(size,baseSize=baseSize))
        
    
        
    return wakeDict

    


def createOctree(fullCaseSetupDict):
    vals = {'Name': 'octree'}
    octree = base.CreateEntity(constants.NASTRAN, "HEXTREME OCTREE", vals)
    pidList = getAllPID()
    usedList = []
    areaDict = {}
    baseCellSize = float(fullCaseSetupDict['BC_SETUP']['BASE_CELL_SIZE'][0])
    print('\t\tCreating octree...')
    print('\t\t\tdomain...')
    areaDict['domain'] = mesh.GetNewOctreeArea( octree, "domain",
                                               baseCellSize,
                                               baseCellSize,
                                               baseCellSize,) #create domain octree
    for pid in pidList:
        pidEntity = base.GetEntity(constants.NASTRAN,"PSHELL",pid)
        boundaryList = ['-min','-max']
        if any(x in pidEntity._name for x in boundaryList):
            mesh.AddPartsToOctreeArea(entities = pidEntity , area = areaDict['domain'])
            usedList.append(pid)

    print('\t\t\tbase...')
    areaDict['base'] = mesh.GetNewOctreeArea( octree, "base",
                                               calculateCellSize(8,baseCellSize),
                                               calculateCellSize(8,baseCellSize),
                                               calculateCellSize(8,baseCellSize),) #create base octree
    for pid in pidList:
        pidEntity = base.GetEntity(constants.NASTRAN,"PSHELL",pid)
        boundaryList = ['-lev8','-coarse']
        if any(x in pidEntity._name for x in boundaryList):
            mesh.AddPartsToOctreeArea(entities = pidEntity , area = areaDict['base'])
            usedList.append(pid)
    print('\t\t\tfine...')
    areaDict['fine'] = mesh.GetNewOctreeArea( octree, "fine",
                                               calculateCellSize(10,baseCellSize),
                                               calculateCellSize(10,baseCellSize),
                                               calculateCellSize(10,baseCellSize),) #create fine octree
    for pid in pidList:
        pidEntity = base.GetEntity(constants.NASTRAN,"PSHELL",pid)
        boundaryList = ['-fine','-lev10']
        if any(x in pidEntity._name for x in boundaryList):
            mesh.AddPartsToOctreeArea(entities = pidEntity , area = areaDict['fine'])
            usedList.append(pid)

    print('\t\t\textra-fine...')
    areaDict['extraFine'] = mesh.GetNewOctreeArea( octree, "extraFine",
                                               calculateCellSize(11,baseCellSize),
                                               calculateCellSize(11,baseCellSize),
                                               calculateCellSize(11,baseCellSize),) #create extrafine octree
    for pid in pidList:
        pidEntity = base.GetEntity(constants.NASTRAN,"PSHELL",pid)
        boundaryList = ['-extra-fine','-lev11']
        if any(x in pidEntity._name for x in boundaryList):
            mesh.AddPartsToOctreeArea(entities = pidEntity , area = areaDict['extraFine'])
            usedList.append(pid)
    
    print('\t\t\tregular...')
    areaDict['regular'] = mesh.GetNewOctreeArea( octree, "domain",
                                               calculateCellSize(9,baseCellSize),
                                               calculateCellSize(9,baseCellSize),
                                               calculateCellSize(9,baseCellSize),) #create regular octree
    for pid in pidList:
        pidEntity = base.GetEntity(constants.NASTRAN,"PSHELL",pid)
        boundaryList = ['-lev9']
        if any(x in pidEntity._name for x in boundaryList) or pid not in usedList:
            mesh.AddPartsToOctreeArea(entities = pidEntity , area = areaDict['regular'])
    print('\n\n\t\tSummary of octree PIDs:')
    for area in areaDict.keys():
        print('\t\t\t%s' % (area))
        ansa.mesh.SaveOctreeAreaParams( areaDict[area], case + '_' + area + '.ansa_mpar' )
        for part in mesh.GetPartsFromOctreeArea( area = areaDict[area] ):
            print('\t\t\t\t',part._name)
        
    saveAnsa(case)

def createOctree2(fullCaseSetupDict,geomDict):
    
    print('\t\tCreating octree...')
    #get mpar file
    mparfile = readMparFile(os.path.join(ansaTemplateLoc,'octree.ansa_mpar'))
    qualityFile = mesh.ReadQualityCriteria(os.path.join(ansaTemplateLoc,'octree_quality.ansa_qual'))
    baseSize = float(fullCaseSetupDict['BC_SETUP']['BASE_CELL_SIZE'][0])

    #starting the octree creation

    #start with global settings
    print('\t\tSetting up global octree settings, getting template from: %s' % (os.path.join(ansaTemplateLoc,'octree.ansa_mpar')))
    globalMparFile = mparfile.copy()
    globalMparFile['octree_parameters_name'] = 'global_octree'
    globalMparFile['curvature_minimum_length'] = baseSize*1000
    globalMparFile['maximum_surface_length'] = baseSize*1000
    globalMparFile['maximum_length'] = baseSize*1000
    globalMparFile['hextreme_layers_flag'] = 'false'
    globalMparFile['distortion_angle_check_box'] = 'false'
    globalMparFile['sharp_angle_check_box'] = 'false'
    globalMparFile['sharp_edges_check_box'] = 'false'
    globalMparFile['free_edges_check_box'] = 'false'
    globalMparFile['intersection_lines_check_box'] = 'false'
    globalMparFile['self_proximity'] = 'false'
    globalMparFile['detect_proximity_with_any_property'] = 'false'
    globalMparFile['pid_bounds'] = 'false'
    globalMparFile['feature_bounds'] = 'false'
    globalMparFile['proximity'] = 'false'
    globalMparFile['oriented_proximity'] = 'false'
    
    if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0].lower() == 'half':
        print('\t\t\tSymmetry mode selected as half!')
        globalMparFile['hextreme_connect_to_symmetry'] = 'true'
        globalMparFile['symmetry_direction'] = 'Y-'
        #set y-max to symmetryPlane type
        
    elif fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0].lower() == 'full':
        print('\t\t\tSymmetry mode selected as full!')
        globalMparFile['hextreme_connect_to_symmetry'] = 'false'
    else:
        sys.exit('ERROR! Wrong symmetry type in GLOBAL_SIM_CONTROL -> SIM_SYM! Valid options: (half/full)')

    writeMparfile(globalMparFile,'global_octree.ansa_mpar')

    #creating octree definition for global
    global_octree = base.CreateEntity(constants.NASTRAN, "HEXTREME OCTREE", {'Name':'global_octree'})
    mesh.ReadOctreeAreaParams( area = global_octree, mpar_file = 'global_octree.ansa_mpar')

    geomPartList = []
    for geom in geomDict.keys():
        geomPrefixList = ['REFX-','REF-','MRF-','MRFG-']
        if not any(x in geom for x in geomPrefixList):
            geomPartList.append(geomDict[geom]['part_ent'])
            
    base.Or(geomPartList)
    assignBCType(geomDict,fullCaseSetupDict)
    pids = base.CollectEntities(constants.NASTRAN,None,'__PROPERTIES__',filter_visible=True)
    mesh.AddPartsToOctreeArea(pids,global_octree)

    #setting z-min layers
    zmin_area = mesh.GetNewOctreeArea(global_octree,name="zmin_area")
    #zminMparFile = globalMparFile.copy()
    zminMparFile = readMparFile(os.path.join(ansaTemplateLoc,'octree_zmin.ansa_mpar'))
    zminMparFile['octree_parameters_name'] = 'zmin_octree'
    zminMparFile['curvature_minimum_length'] = baseSize*1000
    zminMparFile['maximum_surface_length'] = baseSize*1000
    zminMparFile['maximum_length'] = baseSize*1000
    if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0].lower() == 'half':
        print('\t\t\tSymmetry mode selected as half!')
        zminMparFile['hextreme_connect_to_symmetry'] = 'true'
        zminMparFile['symmetry_direction'] = 'Y-'
        #set y-max to symmetryPlane type
        
    elif fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0].lower() == 'full':
        print('\t\t\tSymmetry mode selected as full!')
        zminMparFile['hextreme_connect_to_symmetry'] = 'false'
    else:
        sys.exit('ERROR! Wrong symmetry type in GLOBAL_SIM_CONTROL -> SIM_SYM! Valid options: (half/full)')
    
    writeMparfile(zminMparFile,'zmin_octree.ansa_mpar')
    mesh.ReadOctreeAreaParams( area = zmin_area, mpar_file = 'zmin_octree.ansa_mpar')
    zminpids = []
    for pid in pids:
        if 'z-min' in pid._name:
            zminpids.append(pid)
    mesh.AddPartsToOctreeArea(zminpids,zmin_area) #add zmin pids to octree
    base.SetEntityCardValues(constants.NASTRAN,zmin_area,{'Name':'domain_zmin'})



    #preparing regular geometry mpar file
    geom_area = mesh.GetNewOctreeArea(global_octree,name="geom_area")
    minGeomLevel = int(fullCaseSetupDict['GLOBAL_REFINEMENT']['MIN_GEOM_LEVEL'][0])
    firstLayerHeight = float(fullCaseSetupDict['GLOBAL_REFINEMENT']['FIRST_LAYER_HEIGHT'][0])*1000
    expansionRatio = float(fullCaseSetupDict['GLOBAL_REFINEMENT']['DEF_EX_RATIO'][0])
    refAngle = float(fullCaseSetupDict['GLOBAL_REFINEMENT']['REF_FEAT_ANGLE'][0])
    layerType = fullCaseSetupDict['GLOBAL_REFINEMENT']['LAYER_TYPE'][0]
    geomMparfile = readMparFile(os.path.join(ansaTemplateLoc,'octree_geom.ansa_mpar'))
    geomMparfile['octree_parameters_name'] = 'geom_octree'
    geomMparfile['curvature_minimum_length'] = calculateCellSize(minGeomLevel+1,baseSize)
    geomMparfile['maximum_surface_length'] = calculateCellSize(minGeomLevel,baseSize)
    geomMparfile['maximum_length'] = calculateCellSize(minGeomLevel,baseSize)
    geomMparfile['distortion_angle'] = refAngle
    geomMparfile['hextreme_layers_flag'] = 'true'
    geomMparfile['distortion_angle_check_box'] = 'true'
    geomMparfile['sharp_angle_check_box'] = 'true'
    geomMparfile['sharp_edges_check_box'] = 'true'
    geomMparfile['free_edges_check_box'] = 'true'
    geomMparfile['intersection_lines_check_box'] = 'true'
    geomMparfile['self_proximity'] = 'true'
    geomMparfile['oriented_proximity'] = 'true'
    geomMparfile['hextreme_number_of_layers_value'] = 6
    geomMparfile['hextreme_layers_growth_rate'] = expansionRatio
    geomMparfile['hextreme_layers_first_layer_height'] = firstLayerHeight
    geomMparfile['hextreme_layers_size_mode'] = layerType.title()
    geomMparfile['detect_proximity_with_any_property'] = 'true'
    geomMparfile['pid_bounds'] = 'true'
    geomMparfile['feature_bounds'] = 'true'
    geomMparfile['proximity'] = 'true'
    #geomMparfile['proximity_minimum_length'] = calculateCellSize(minGeomLevel+1,baseSize)
    geomMparfile['proximity_minimum_length'] = calculateCellSize(minGeomLevel,baseSize)
    if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0].lower() == 'half':
        print('\t\t\tSymmetry mode selected as half!')
        geomMparfile['hextreme_connect_to_symmetry'] = 'true'
        geomMparfile['symmetry_direction'] = 'Y-'
        #set y-max to symmetryPlane type
        
    elif fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0].lower() == 'full':
        print('\t\t\tSymmetry mode selected as full!')
        geomMparfile['hextreme_connect_to_symmetry'] = 'false'
    else:
        sys.exit('ERROR! Wrong symmetry type in GLOBAL_SIM_CONTROL -> SIM_SYM! Valid options: (half/full)')

    geompids = []
    for pid in pids:
        if not 'domain_' in pid._name:
            geompids.append(pid)
    writeMparfile(geomMparfile,'geom_octree.ansa_mpar')
    mesh.ReadOctreeAreaParams( area = geom_area, mpar_file = 'geom_octree.ansa_mpar')
    mesh.AddPartsToOctreeArea(geompids,geom_area) #add geom pids to octree
    base.SetEntityCardValues(constants.NASTRAN,geom_area,{'Name':'global_geom_level_%s' % (minGeomLevel)})

    #getting pid refinement table
    refTable = getPIDRefTable(os.path.join(ansaTemplateLoc,'ref_tables'),baseSize)

    if os.path.exists('pidMpars'):
        print('pidMpars directory exists, not making!')
    else:
        print('pidMpars does not exists, making!')
        os.mkdir('pidMpars')
    #starting general geometry settings
   
    for geom in geomDict.keys():
        geomPrefixList = ['domain','REFX-','REF-','MRF-','MRFG-']
        if not any(x in geom for x in geomPrefixList):
            partMparfile = geomMparfile.copy()
            
            if 'GEOMX-' in geom:
                part_area = mesh.GetNewOctreeArea(global_octree,name=geom.split('.')[0])
                refLevel = fullCaseSetupDict[geom.split('.')[0]]['GEOMX_MIN_MAX_LEVEL']
                minLevel = refLevel[0]
                maxLevel = refLevel[1]
                edgeIncrement = fullCaseSetupDict[geom.split('.')[0]]['GEOMX_FEAT_EDGE_LEVEL_INC'][0]
                refAngle = fullCaseSetupDict[geom.split('.')[0]]['GEOMX_REFINEMENT_ANGLE'][0]
                nlayers = fullCaseSetupDict[geom.split('.')[0]]['GEOMX_NLAYERS'][0]
                firstLayerHeight = float(fullCaseSetupDict[geom.split('.')[0]]['FIRST_LAYER_HEIGHT'][0])*1000
                layerType = fullCaseSetupDict[geom.split('.')[0]]['LAYER_TYPE'][0]
                growthRate = fullCaseSetupDict[geom.split('.')[0]]['GEOMX_LAYER_RATIO'][0]

                print('\t\t',geom,'\sNLAYERS: ',nlayers)
                print('\t\t',geom,'\sFIRST LAYER HEIGHT: ',firstLayerHeight)

                
                partMparfile.update({'octree_parameters_name':geom.split('.')[0],
                                     'curvature_minimum_length':calculateCellSize(int(maxLevel),baseSize),
                                     'maximum_surface_length':calculateCellSize(int(minLevel),baseSize),
                                     'maximum_length':calculateCellSize(int(minLevel),baseSize),
                                     'hextreme_number_of_layers_value':int(nlayers),
                                     'hextreme_layers_growth_rate':float(growthRate),
                                     'hextreme_layers_first_layer_height':float(firstLayerHeight),
                                     'hextreme_layers_size_mode':layerType.title(),
                                     'distortion_angle':refAngle,
                                     'proximity_minimum_length':calculateCellSize(int(maxLevel)+int(edgeIncrement),baseSize),
                                     'sharp_edges_length':calculateCellSize(int(maxLevel)+int(edgeIncrement),baseSize),
                                     'free_edges_length':calculateCellSize(int(maxLevel)+int(edgeIncrement),baseSize),

                })
                writeMparfile(partMparfile,'%s_octree.ansa_mpar' % (geom.split('.')[0]))
                mesh.ReadOctreeAreaParams( area = part_area, mpar_file = '%s_octree.ansa_mpar' % (geom.split('.')[0]))
                geompids = getPIDsFromPart(geomDict[geom]['part_ent'])
                mesh.AddPartsToOctreeArea(geompids,part_area) #add geom pids to octree
                base.SetEntityCardValues(constants.NASTRAN,part_area,{'Name':geom.split('.')[0]})
                print('\n\t\tAdding PID refinements...')
                n = 0
                for geompid in geompids:
                    refPartMparfile = partMparfile.copy()
                    pidRefTable = None
                    for key in refTable.keys():
                        #going through each key word, if key word is in pid name, that ref table is assigned
                        if key in geompid._name:
                            #print(key)
                            pidRefTable = refTable[key]
                            #will keep looping until the last one is used (if exists) avoids duplicates
                    if pidRefTable != None:
                        print('\t\t\tAdding %s to %s refinement' % (geompid._name,pidRefTable['octree_parameters_name']))
                        refPartMparfile.update(pidRefTable)
                        ref_part_area = mesh.GetNewOctreeArea(global_octree,name=geompid._name)
                        writeMparfile(refPartMparfile,'pidMpars/%s_octree.ansa_mpar' % (geompid._name))
                        mesh.ReadOctreeAreaParams( area = ref_part_area, mpar_file = 'pidMpars/%s_octree.ansa_mpar' % (geompid._name))
                        mesh.AddPartsToOctreeArea(geompid,ref_part_area) #add geom pids to octree
                        base.SetEntityCardValues(constants.NASTRAN,ref_part_area,{'Name':geompid._name})
                        n = n + 1
                if n > 0:
                    print('\n\t\tAdded %s PID refinements!' % (n))

                
            else:
                part_area = mesh.GetNewOctreeArea(global_octree,name=geom.split('.')[0])
                refLevel = geomDict[geom]['refinement']
                nlayers = geomDict[geom]['layers']
                growthRate = geomDict[geom]['expansion']

                
                partMparfile.update({'octree_parameters_name':geom.split('.')[0],
                                     'curvature_minimum_length':calculateCellSize(int(refLevel)+1,baseSize),
                                     'sharp_edges_length':calculateCellSize(int(refLevel)+1,baseSize),
                                     'proximity_minimum_length':calculateCellSize(int(refLevel)+1,baseSize),
                                     'maximum_surface_length':calculateCellSize(int(refLevel),baseSize),
                                     'maximum_length':calculateCellSize(int(refLevel),baseSize),
                                     'hextreme_number_of_layers_value':int(nlayers),
                                     'hextreme_layers_growth_rate':float(growthRate),

                })
                writeMparfile(partMparfile,'%s_octree.ansa_mpar' % (geom.split('.')[0]))
                mesh.ReadOctreeAreaParams( area = part_area, mpar_file = '%s_octree.ansa_mpar' % (geom.split('.')[0]))
                geompids = getPIDsFromPart(geomDict[geom]['part_ent'])
                mesh.AddPartsToOctreeArea(geompids,part_area) #add geom pids to octree
                base.SetEntityCardValues(constants.NASTRAN,part_area,{'Name':'%s_level_%s' % (geom.split('.')[0],refLevel)})
                print('\n\t\tAdding PID refinements...')
                n = 0
                for geompid in geompids:
                    refPartMparfile = partMparfile.copy()
                    pidRefTable = None
                    for key in refTable.keys():
                        #going through each key word, if key word is in pid name, that ref table is assigned
                        if key in geompid._name:
                            #print(key)
                            pidRefTable = refTable[key]
                            #will keep looping until the last one is used (if exists) avoids duplicates
                    if pidRefTable != None:
                        print('\t\t\tAdding %s to %s refinement' % (geompid._name,pidRefTable['octree_parameters_name']))
                        refPartMparfile.update(pidRefTable)
                        ref_part_area = mesh.GetNewOctreeArea(global_octree,name=geompid._name)
                        writeMparfile(refPartMparfile,'pidMpars/%s_octree.ansa_mpar' % (geompid._name))
                        mesh.ReadOctreeAreaParams( area = ref_part_area, mpar_file = 'pidMpars/%s_octree.ansa_mpar' % (geompid._name))
                        mesh.AddPartsToOctreeArea(geompid,ref_part_area) #add geom pids to octree
                        base.SetEntityCardValues(constants.NASTRAN,ref_part_area,{'Name':geompid._name})
                        n = n + 1
                if n > 0:
                    print('\n\t\tAdded %s PID refinements!' % (n))

                    



    print('\t\tSetting seed point(s):')
    seedPoints = np.array(fullCaseSetupDict['GLOBAL_REFINEMENT']['LOC_IN_MESH'],dtype=float)*1000
    print('\t\t\t',seedPoints/1000)
    mesh.OctreeSeedPoints(global_octree, seedPoints[0], seedPoints[1], seedPoints[2], "global_octree_volume", baseSize*1000)
    print('\t\tOctree Summary:')
    for area in mesh.GetAreasFromOctree( octree = global_octree ):
        print('\n\t\t\t',area._name)
        parts = mesh.GetPartsFromOctreeArea( area = area )
        for part in parts:
            print('\t\t\t\t',part._name)
    
    saveAnsa(case)
    mesh.RunOctree(global_octree)
    runQualityImprovement()
    saveAnsa(case)

def getPIDRefTable(refTablePath,baseSize):
    #geomMparfile = readMparFile(os.path.join(ansaTemplateLoc,'octree_geom.ansa_mpar'))
    print('\t\tImporting PID refinement tables!')
    with open(refTablePath, "r") as data:
        refTable = ast.literal_eval(data.read())
    
    for key in refTable.keys():
        for var in refTable[key].keys():
            if 'name' in var:
                continue
            try:
                level = int(refTable[key][var])
            except:
                sys.exit('ERROR! Unable to convert level string to int() in ref_table, please check entry!')
            size = calculateCellSize(level,baseSize)
            refTable[key][var] = size

    for keyWord in refTable.keys():
        print('\t\t\t',keyWord)
        for var in refTable[keyWord].keys():
            print('\t\t\t\t',var,": ",refTable[keyWord][var])
    return refTable

def assignBCType(geomDict,fullCaseSetupDict):
    print('\t\tAssigning boundary conditions to PIDs')

    bcDict = {'half':{'x-min':'patch',
                      'x-max':'patch',
                      'y-min':'wall',
                      'y-max':'symmetry',
                      'z-min':'wall',
                      'z-max':'wall',},
              'full':{'x-min':'patch',
                      'x-max':'patch',
                      'y-min':'patch',
                      'y-max':'patch',
                      'z-min':'wall',
                      'z-max':'wall',}}
    caseSym = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0].lower()
    pids = base.CollectEntities(constants.NASTRAN,None,'__PROPERTIES__',filter_visible=True)
    for bc in bcDict[caseSym]:
        for pid in pids:
            if bc in pid._name:
                print('\t\t\tSetting: %s to %s' % (pid._name,bcDict[caseSym][bc]))
                base.SetEntityCardValues(constants.OPENFOAM, pid, {"TYPE":bcDict[caseSym][bc]})
                
    


    
    


def writeMparfile(mparfile,filePath):
    mparlines = []
    for var in mparfile.keys():
        val = mparfile[var]
        mparlines.append('%s=%s\n' % (var,val))


    with open(filePath,'w') as file:
        file.write("".join(mparlines))
    
    return 

def readMparFile(filePath):
    mparDict = {}
    with open(filePath, 'r') as file:
        for line in file:
            # Process each line here
            line = line.strip()
            if not line.startswith('#'):
                
                lineVar = line.split('=')[0].strip()
                if lineVar == '':
                    continue
                try:
                    lineVal = line.split('=')[1].strip()
                except:
                    lineVal = ''
                mparDict[lineVar] = lineVal
    return mparDict
    




def createSizeField(fullCaseSetupDict,geomDict):
    print('\t\tCreating size fields...')
    sf = base.CreateEntity(constants.NASTRAN, "SIZE FIELD")
    baseSize = float(fullCaseSetupDict['BC_SETUP']['BASE_CELL_SIZE'][0])
    refinementLevels = fullCaseSetupDict['GLOBAL_REFINEMENT']['GEOM_REF_LEVELS']
    refinementDist = fullCaseSetupDict['GLOBAL_REFINEMENT']['GEOM_REF_DIST']
    #create domain size field
    #print('\t\t\tDomain')
    #createRefSizeField(sf,geomDict['domain.stl']['part_ent'],0,baseSize)
    n = 0
    print('\t\t\tglobal size fields...')
    for level,dist in zip(refinementLevels,refinementDist):
        size = calculateCellSize(int(level),baseSize)
        dist = float(dist)*1000
        sfname = 'level_%s' % (level)
        #print('\t\t\t%s' % (sfname))
        global_offset = mesh.GetNewSizeFieldSurfaceOffsetRule(size_field=sf,
                                                              name=sfname,
                                                              max_surf_len = size,
                                                              max_vol_len = size,
                                                              offset = dist)
        for part in geomDict.keys():
            part_ent = geomDict[part]['part_ent']
            pids = getPIDsFromPart(part_ent)
            for pidEntity in pids:
                pidName = pidEntity._name
                boundaryList = ['-min','-max','REFX-','REF-','MRFG-'] # don't create offset size fields for these types
                if not any(x in pidName for x in boundaryList):
                    #print('\t\t\t\t%s' % (pidName))
                    mesh.AddContentsToSizeFieldRule( entities = pidEntity , rule = global_offset)

    for geom in geomDict:
        if geom.startswith('GEOMX-'):
            geomxRefDist = fullCaseSetupDict[geom.split('.')[0]]['GEOMX_REF_DIST']
            geomxRefLevel = fullCaseSetupDict[geom.split('.')[0]]['GEOMX_REF_LEVELS']
            for level,dist in zip(geomxRefLevel,geomxRefDist):
                size = calculateCellSize(int(level),baseSize)
                dist = float(dist)*1000
                sfname = '%s_level_%s' % (geom.split('.')[0],level)
                geomx_offset = mesh.GetNewSizeFieldSurfaceOffsetRule(size_field=sf,
                                                                    name=sfname,
                                                                    max_surf_len = size,
                                                                    max_vol_len = size,
                                                                    offset = dist)
                geomxPids = getPIDsFromPart(geomDict[geom]['part_ent'])
                mesh.AddContentsToSizeFieldRule( entities = geomxPids , rule = geomx_offset)


            

    sf_rules = mesh.GetRulesFromSizeField(sf)
    print('\t\tSize field summary:')
    for rule in sf_rules:
        print('\t\t\t' + rule._name)
        print('\t\t\t\tNumber of parts: ',str(len(mesh.GetContentsFromSizeFieldRule(rule))))
        for ent in mesh.GetContentsFromSizeFieldRule(rule):
            print('\t\t\t\t',ent._name)
        mesh.SaveSizeFieldRuleParams( rule = rule, mpar_file = '%s.ansa_mpar' % (rule._name))

    print('\t\tAdding wake boxes to size field...')
    wakeDict = createRefinementRegion(templateLoc,fullCaseSetupDict,geomDict)
    wakeList = []
    for wake in wakeDict.keys():
        wakeList.append(wakeDict[wake])
    mesh.ConvertSizeBoxesToSizeField( size_boxes = wakeList , size_field = sf)

    print('\n\t\tChecking for refinement regions...')
    refinementList = []
    refinementKeys = ['REF-','REFX-']
    for geom in geomDict.keys():
        if any(x in geom for x in refinementKeys):
            refinementList.append(geomDict[geom]['part_ent'])
    
    if len(refinementList) > 0:
        print('\t\t\tFound %s refinement regions:' % (len(refinementList)))
        for refRegion in refinementList:
            print('\t\t\t\tCreating size field for: %s' % (refRegion._name))
            if refRegion._name.startswith('REFX-'):
                refinement = fullCaseSetupDict[refRegion._name.split('.')[0]]['REF_LEVEL'][0]
            else:
                refinement = geomDict[refRegion._name]['refinement']
            
            createRefSizeField(sf,refRegion,refinement,baseSize)
           


    else:
        print('\t\t\tNo refinement regions found!')




    print('\t\t\tBuilding size field!')
    mesh.BuildSizeField(sf)
    saveAnsa(case)

    return sf

def createRefSizeField(size_field,part_ent,refinement,baseSize):
    try:
        refinement = int(refinement)
    except Exception as E:
        print('ERROR! Refinement region %s does not have a valid refinement level, must be int!' % (part_ent._name))
        print(E)

    
    closed_surf = mesh.GetNewSizeFieldClosedSurfaceRule( size_field = size_field,
                                                        name = part_ent._name.split('.')[0],
                                                        max_surf_len = calculateCellSize(int(refinement),baseSize),
                                                        max_vol_len =  calculateCellSize(int(refinement),baseSize))
    
    mesh.AddContentsToSizeFieldRule(entities = getPIDsFromPart(part_ent) , rule = closed_surf)
    
     
def getPIDsFromPart(part_ent):
    base.Or(part_ent)
    return base.CollectEntities(constants.NASTRAN,None,'__PROPERTIES__',filter_visible=True)
def calculateCellSize(level,baseSize):
    return (baseSize*1000)/(2**level)





def importDomain():
    print('\t\tCreating blockMesh to import to ANSA...')
    os.system('cp system/controlDictPotential system/controlDict;blockMesh;foamToSurface constant/triSurface/domain.stl')
    base.InputStereoLithography(filename='constant/triSurface/domain.stl',
                                unit_system=utils.UnitSystem(length='meter'))
    

def getAllPID():

    deck = constants.NASTRAN
    pidList = []

        
    pids = base.CollectEntities(constants.NASTRAN, None, "__PROPERTIES__")
    for pid in pids:
        part = base.GetEntityPart(pid)
        pidList.append(pid)

    return pidList
        

def getParts(geomDict):
    
    deck = constants.OPENFOAM
    parts = base.CollectEntities(deck,None,'ANSAGROUP')
        
    for part in parts:
        part_vals = base.GetEntityCardValues(
        deck,
        part,
        ["Name", "Module Id"]
        )
        #make only the part visible
        base.Or(part)
        part_name = part_vals["Name"]
        print('\t\t\t%s' % (part_name))
        #get the pids for only the part that is visible
        #pids = base.CollectEntities(deck,part,'__PROPERTIES__')
        pids = base.CollectEntities(deck,None,'__PROPERTIES__',filter_visible=True)
        for pid in pids:
            #opid = base.GetEntity(constants.NASTRAN,'SHELL',pid._id)
            base.SetEntityCardValues(
                deck,
                pid,
                {"Name": "%s_%s" % (part_name.replace('.obj','')
                                    .replace('.stl','')
                                    .replace('.gz',''),
                                    pid._name)} # add part name to the PID for OF use
            )

            print('\t\t\t\t%s' %  (pid._name))

        if 'domain.stl' in part_name:
            geomDict['domain.stl'] = {}
            geomDict['domain.stl']['scale'] = 1
            geomDict['domain.stl']['refinement'] = 0
            geomDict['domain.stl']['wallmodel'] = 'high'
            geomDict['domain.stl']['layers'] = 6
            geomDict['domain.stl']['expansion'] = 1.2

        geomDict[part_name]['part_ent'] = part #adds the part entity into the geomDict for future use
 
    return geomDict
    
    



def exportAnsaMesh():

    base.All()
    ansaParts = base.CollectEntities(constants.OPENFOAM,None,'ANSAGROUP')
    for ansaPart in ansaParts:
        if 'global_octree' in ansaPart._name:
            base.Or(ansaPart)

    print('Writing out OpenFOAM mesh!')
    base.OutputOpenFoam("%s/%s/" % (path,case), mode="visible",
                        version="ESI",
                        solver_info="Off",
                        initial_conditions_folder="Off",
                        unit_system = utils.UnitSystem(length='meter'),
                        binary_io_64bit_integers="Off",
                        )
    



    


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
    

main()