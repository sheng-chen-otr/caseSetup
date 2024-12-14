import os
import sys
import numpy as np
import pandas as pd
import time as time
import re
import configparser

#import paraview modules
from paraview.simple import *
from paraview.numpy_support import vtk_to_numpy


#disable first views
paraview.simple._DisableFirstRenderCameraReset()

#getting case information
casePath = os.getcwd() #path of directory this script is being run in
caseName = casePath.split('/')[-1] #directory name 


#debug options
USE_PRE_DEF_VARS = False #using only predefined variables in list PRE_DEF_VAR_LIST
PRE_DEF_VAR_LIST = ['UMean','pMean']

#main function
def main():
    #global renderView, caseName, availableCellArrays
    global UREF,LREF,CREF,FREF,fullCaseSetupDict
    #pvPostSetupDict = getPVSetup()

    print('''\n\t####\t\tEZ-CFD PARAVIEW POST-PROCESSING V1.0\t\t ####\n\n''')

    # Start timing the script execution
    begin = time.time()

    print('\tProcessing case: %s' % (caseName))

    
    #Getting reference values and geometries
    UREF,LREF,CREF,FREF,fullCaseSetupDict = getRef()

    #getting the boundaries
    geomKeys = list(getGeometry(fullCaseSetupDict).keys())
    geomList = []
    for geom in geomKeys:
        #remove internal domain surfaces if any
        if not 'IDOM' in geom:
            geomList.append(geom.split('.')[0])
    

    
    

    
    # Initialize the renderView
    print('\n\tGenerating the initial render view')
    renderView                           = GetRenderView()
    renderView.ViewSize                  = [4096,3000] # default
    renderView.CameraParallelProjection  = 1
    renderView.Background                = [1,1,1]
    renderView.OrientationAxesVisibility = 0
    renderView.CenterAxesVisibility      = 0
    renderView.CameraPosition            = [0,0,0]
    renderView.CameraFocalPoint          = [0,0,0]
    renderView.CameraParallelScale       = 1

    # Load in the data file from EnSight
    print("\n\tLoading data...")
    try:
        filename = '%s/EnSight/%s.case' % (casePath,caseName)
    except:
        sys.exit('ERROR! Unable to find EnSight File!')
    
    # read the ensight file
    reader   = EnSightReader(CaseFileName=filename)

    print('\tInitial data read, setting cell arrays')

    #getting available time steps in the simulation
    timesteps  = reader.TimestepValues
    print("\t\tTimesteps in Case:")
    for i in timesteps:
        print('\t\t\t' + str(i))

    print('\t\tUsing latest time step: %s' % (timesteps[-1]))
    renderView.ViewTime = timesteps[-1]

    #print available cell arrays
    #loads all available variables by default
    print('\t\tVariables in Case:')
    for var in reader.CellArrays:
        print('\t\t\t' + var)

    #if using predefined variables is selected...
    if USE_PRE_DEF_VARS == True:
        print('\t\tDEBUG_OPTION: Using Predefined Variable List!')
        reader.CellArrays = PRE_DEF_VAR_LIST

    print('\tCell Arrays set, updating pipeline')
    reader.UpdatePipeline()
    source   = CellDatatoPointData(Input=reader)
    print('\tCompleted cell to point data')
    source.UpdatePipeline()
    print('\tData Read Time (s): ' + str(round(time.time()-begin,3)) + '\n\n')

    #getting the boundaries
    
    print('\tGetting boundaries...')
    boundaries, selections = getBoundaries(source)
    include_patches = geomList
    selectors = includeSelectors(selections,include_patches)
    geomSurface  = ExtractBlock(Input = source)

    
def getPVSetup():
    '''
        Looks for the setup config file to know the settings
        - if the file cannot be found in CWD it will use the default setup file
        - will check to see if all the required variables exists, if not, it will
        populate based on the default config file
        
    '''
    pass
    pvPostSetupDict = 1
    return pvPostSetupDict

def getRef():
    '''
        gets all the reference values based on the fullCaseSetupDict file
    '''
    print('\n\tReading fullCaseSetupDict from case directory...')
    caseSetupConfig = configparser.ConfigParser()
    caseSetupConfig.optionxform = str
    #tries to read the caseSetup file, if it can't find it, it throws an error
    try:
        caseSetupConfig.read_file(open("%s/fullCaseSetupDict" % (casePath)))
    except Exception as e:
        print('\n\n\nERROR! fullCaseSetupDict invalid!')
        print(e)
        exit()

    caseSetupConfigSections = caseSetupConfig.sections()
    print('\n\t\tSections:')
    for sections in caseSetupConfigSections:
        print('\t\t\t' + sections)

    UREF = float(caseSetupConfig['BC_SETUP']['INLET_MAG'])
    LREF = float(caseSetupConfig['BC_SETUP']['REFLEN'])
    CREF = np.array(caseSetupConfig['BC_SETUP']['REFCOR'].split(' ')).astype('float')
    FREF = np.array(caseSetupConfig['BC_SETUP']['FRT_WH_CTR'].split(' ')).astype('float')

    print('\t\tReference Values:')
    print('\t\t\tUREF: %1.3g, LREF: %1.3g, CREF: [%1.3g %1.3g %1.3g] , FREF: [%1.3g %1.3g %1.3g]' % (UREF,LREF,CREF[0],CREF[1],CREF[2],FREF[0],FREF[1],FREF[2]))
    
    return UREF,LREF,CREF,FREF,caseSetupConfig

def getBoundaries(source):
    boundaries = dict()
    selections = dict()
    nDataSets = source.GetDataInformation().GetNumberOfDataSets()
    print("\n\t\tFinding Blocks:\n\t\t\tNumber of Blocks = {}".format(nDataSets)) # Used to be called blocks
    
    counter = 0
    for i in range(1,nDataSets,1):
        counter += 1
        print("\t\t\t\tFound Block Name {} ".format(source.GetDataInformation().GetBlockName(i)))
        boundaries[source.GetDataInformation().GetBlockName(i)] = counter
        selections['/Root/{}'.format(source.GetDataInformation().GetBlockName(i))] = counter
        counter += 1

    return boundaries, selections

def includeBoundaries(boundaries,include_patches):
    allRegions  = list(boundaries.keys())
    meshRegions = []
    print('\n\tIncluding patches containing: ',include_patches)
    for boundary in boundaries:
        for patch in include_patches:
            if patch.lower() in boundary.lower() and boundary not in meshRegions:
                meshRegions.append(boundary)

    print('\n\tUsing patches:')
    for meshRegion in meshRegions:
        print('\t\t' + meshRegion)


    print('\n')
    return meshRegions

def includeSelectors(selections,include_patches):
    selectors = []
    print('\t\t\tSelected boundaries for surface:')
    for selector in selections.keys():
        if any(x in selector for x in include_patches):
            selectors.append(selector)
            print('\t\t\t\t%s' % (selector))

    return selectors

def getGeometry(fullCaseSetupDict):    
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
    

    return geomDict

def geomToDict(geomDict,geometryList,geomColumnNames):
    
    try:
        print('\t\tGetting geometry...')
        geometryList = geometryList.split('\n')
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












