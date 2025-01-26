import os
import sys
import numpy as np
import pandas as pd
import time as time
import re
import configparser
import glob
import json

#import paraview modules
from paraview.simple import *
from paraview.numpy_support import vtk_to_numpy


#disable first views
paraview.simple._DisableFirstRenderCameraReset()

#getting case information
casePath = os.getcwd() #path of directory this script is being run in
caseName = casePath.split('/')[-1] #directory name
#templateLoc = '/home/openfoam/openFoam/scripts/caseSetup/caseSetup/postUtilities'
templateLoc = (os.path.dirname(os.path.realpath(__file__)))

ASPECT_RATIO = 16/9
YRES = 3000
XRES = ASPECT_RATIO * YRES
RESOLUTION = [int(XRES),int(YRES)]

#debug options
USE_PRE_DEF_VARS = False #using only predefined variables in list PRE_DEF_VAR_LIST
PRE_DEF_VAR_LIST = ['UnwMean','pMean']

#main function
def main():
    #global renderView, caseName, availableCellArrays
    global UREF,LREF,CREF,FREF,WREF,fullCaseSetupDict,renderView
    

    print('''\n\t####\t\tEZ-CFD PARAVIEW POST-PROCESSING V1.0\t\t ####\n\n''')

    # Start timing the script execution
    begin = time.time()

    print('\tProcessing case: %s' % (caseName))

    
    #Getting reference values and geometries
    UREF,LREF,WREF,CREF,FREF,fullCaseSetupDict = getRef()
    
    #Getting pvPostSetup file
    pvPostSetupDict = getPVSetup(templateLoc)
    if pvPostSetupDict['PV_POST_MAIN']['VARIABLE_DICT'] == 'default':
        varDict,viewsDict = getVariableDicts('%s/default/defaultVariables'% (templateLoc),pvPostSetupDict['PV_POST_MAIN']['CAMERA_VIEWS'])
    else:
        varDict,viewsDict = getVariableDicts(pvPostSetupDict['PV_POST_MAIN']['VARIABLE_DICT'],pvPostSetupDict['PV_POST_MAIN']['CAMERA_VIEWS'])
    
    
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
    renderView.ViewSize                  = RESOLUTION # default
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
    include_patches = geomList #including the geometry surfaces that are PID's
    #selecting the geometry surfaces
    selectors = includeSelectors(selections,include_patches)
    
    #extract the geometry surface for surface processing
    geomSurf                         = ExtractBlock(Input = source)
    geomSurf.Selectors               = selectors
    if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower() == 'half':
        print('\t\tFound case as half, reflecting surfaces...')
        geomSurface = Reflect(Input = geomSurf)
        geomSurface.Plane = 'Y'
        geomSurface.Center = CREF[1]
    else:
        geomSurface = geomSurf

    #geomSurface.UpdatePipeline()

    #extracting internal domain
    internalVol                      = ExtractBlock(Input = source)
    internalVol.Selectors            = '/Root/internalMesh'

    if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower() == 'half':
        print('\t\tFound case as half, reflecting volume...')
        internalVolume = Reflect(Input = internalVol)
        internalVolume.Plane = 'Y'
        internalVolume.Center = CREF[1]
    else:
        internalVolume = internalVol

    #internalVolume.UpdatePipeline()
    if pvPostSetupDict['PV_POST_MAIN']['SURFACE_IMG'].lower() == 'true':
        generateSurfaceContours(geomSurface,renderView,pvPostSetupDict['SURFACE']['VARIABLES'],pvPostSetupDict['SURFACE']['VIEWS'],varDict,viewsDict)

    if pvPostSetupDict['PV_POST_MAIN']['SLICE_IMG'].lower() == 'true':
        generateSlices(internalVolume,renderView,pvPostSetupDict['SLICE'],varDict,viewsDict)

def getVariableDicts(variablePaths, viewsPath):
    
    #getting default variables
    with open(variablePaths) as f: 
        varDict = f.read()
        varDict = varDict.replace('\\','\\\\')\
                         .replace('UREF',str(UREF))\
                         .replace('LREF',str(LREF))\
                         .replace('WREF',str(WREF)) #replaces certain key words with respective reference values

    varDict = json.loads(varDict)
    print('\t\tAvailable variables:')
    for key in varDict.keys():
        print('\t\t\t%s:'% (key))
        for postType in varDict[key].keys():
            print('\t\t\t\t%s:' % (postType))
            for variable in varDict[key][postType].keys():
                print('\t\t\t\t\t%s: %s' % (variable,varDict[key][postType][variable]))

    #if views path is not default it will use the path given to import the dict file
    #replaces the key words with values
    if viewsPath.lower() != 'default':
        with open(viewsPath) as f: 
            viewsDict = f.read()
            viewsDict = viewsDict.replace('UREF',str(UREF))\
                                 .replace('LREF',str(LREF))\
                                 .replace('WREF',str(WREF))\
                                 .replace('CREF',str(CREF))\
                                 .replace('FREF',str(FREF)) #replaces certain key words with respective reference values

        viewsDict = json.loads(viewsDict)
        print('\t\tAvailable views:')
        for key in viewsDict.keys():
            print('\t\t\t%s' % (key))
            for item in viewsDict[key].keys():
                if any(x in item for x in ['campos','focalpos','viewup']):
                    val = np.array(eval(viewsDict[key][item]))
                    viewsDict[key][item] = val
                else:
                    val = viewsDict[key][item]
                print('\t\t\t\t%s: %s' % (item,str(val)))
    else:
        viewsDict = generateDefaultViews(LREF,CREF,FREF,WREF)

    return varDict,viewsDict


def generateSurfaceContours(surfaceSource,renderView,surfaceVars,views,varDict,viewsDict):
    '''
        Generates surface contour plots, the variables are given in a csv file.

    '''
    beginSurface = time.time()
    print('\tGenerating surface contour plots...')
    if surfaceVars == 'default':
        surfaceVars = ['Geom','CpMean','UnwMean','UnwMeanX','UnwMeanY','UnwMeanY','UnwMeanZ','CpPrime2Mean','CfMean']
    if views == 'default':
        views = ['Rear','Front','Left','Right','Top','Bottom','FrontRight','FrontLeft','RearRight','RearLeft']


    print('\t\tRequested surface variables: %s' % (surfaceVars))

    availableVars = list(varDict['surfaceVariables'].keys())


    for variable in surfaceVars:
        if 'geom' in variable.lower():
            surfaceSourceDisplay = Show(surfaceSource,renderView,'UnstructuredGridRepresentation')
            surfaceSourceDisplay.Representation = 'Surface'
        else:

            #check that variable exists in varDict
            surfaceVarsAvailable = varDict['surfaceVariables'].keys()
            if not variable in surfaceVarsAvailable:
                print('\n\t\t%s not available, skipping...\n' % (variable))
                continue

            #calculating variable using variable equation
            calculator = Calculator(registrationName='calculator', Input=surfaceSource)
            calculator.Function = str(varDict['surfaceVariables'][variable]['equation'])
            calculator.ResultArrayName = variable
            
            surfaceSourceDisplay = Show(calculator,renderView,'UnstructuredGridRepresentation')
            surfaceSourceDisplay.Representation = 'Surface'
            surfaceSourceDisplay.BackfaceRepresentation = 'Cull Frontface'
            
            
            ColorBy(surfaceSourceDisplay,('POINTS',variable))
            LUT = GetColorTransferFunction(variable)
            PWF  = GetOpacityTransferFunction(variable)
            title = varDict['surfaceVariables'][variable]['label']
            varRange = np.array(varDict['surfaceVariables'][variable]['range'])
            tableValues = 15
            color = varDict['surfaceVariables'][variable]['color']

            LUT.NumberOfTableValues = tableValues #default
            LUT.RescaleTransferFunction(varRange[0],varRange[1])
            PWF.RescaleTransferFunction(varRange[0],varRange[1])
            LUT.ApplyPreset(color,True)
            PWF.ApplyPreset(color,True)
            renderView.Update()


            colorBar = GetScalarBar(LUT,renderView)
            colorBar.Title = title
            colorBar.TitleFontFamily = 'Times'
            colorBar.LabelFontFamily = 'Times'
            colorBar.ComponentTitle = ''
            colorBar.Orientation = 'Horizontal'
            colorBar.WindowLocation = 'Lower Center'
            if  'Cf' in variable:
                colorBar.LabelFormat = '%-1.3g'
                colorBar.RangeLabelFormat = '%-1.3g'
                colorBar.CustomLabels = np.linspace(varRange[0],varRange[1],3)         
            else:
                colorBar.LabelFormat = '%-1.1f'
                colorBar.RangeLabelFormat = '%-1.1f'
                colorBar.CustomLabels = np.linspace(varRange[0],varRange[1],5)         
            colorBar.ScalarBarLength = 0.3
            colorBar.ScalarBarThickness = 160
            colorBar.TitleFontSize = 160
            colorBar.LabelFontSize = 160
            colorBar.TitleColor = [0,0,0]
            colorBar.LabelColor = [0,0,0]
            colorBar.AddRangeLabels = 1
            colorBar.AutomaticLabelFormat = 0
            colorBar.DrawAnnotations = 0
            colorBar.DrawTickLabels = 1
            colorBar.UseCustomLabels = 1
               
                    

        HideUnusedScalarBars()  
        for view in views:
            renderView = setView(renderView,viewsDict[view])
            saveImages(renderView,caseName,variable,'Surface',view)
        
        try:
            Delete(colorBar)
        except:
            print('')
        Hide(surfaceSourceDisplay,renderView)
        Delete(surfaceSourceDisplay)
        try:
            Delete(calculator)
        except:
            print('')

    print('\tSurface Contour Generation Time (s): ' + str(round(time.time()-beginSurface,3)) + '\n\n')

def generateSlices(volumeSource,renderView,sliceDict,varDict,viewsDict):
    

    print('\tGenerating slices...')

    #set defaults
    if sliceDict['VARIABLES'].lower() == 'default':
        sliceVars = ['CpMean','CptMean','UMean','CpPrime2Mean','UMeanX','UMeanY','UMeanZ','vorticityMean','vorticityMeanX','vorticityMeanY','vorticityMeanZ']
    else:
        sliceVars = sliceDict['VARIABLES'].split()

    if sliceDict['VIEWS'].lower() == 'default':
        sliceViews = ['Front','Left','Top']
    else:
        sliceViews = sliceDict['VIEWS'].split()

    if sliceDict['NORMALS'].lower() =='default':
        normalsList = ['X','Y','Z']
    else:
        normalsList = sliceDict['NORMALS'].split()

    if sliceDict['NSLICES'].lower() =='default':
        nSliceList = [50,20,20]

        if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower() == 'half':
            nSliceList = [50,20,20]
        else:
            nSliceList = [50,40,20]
    else:
        nSliceList = sliceDict['NSLICES'].split()


    if sliceDict['SLICE_RANGE'].lower() =='default':
        frontAxleXLoc = float(FREF[0])
        frontAxleYLoc = float(FREF[1])
        frontAxleZLoc = float(FREF[2])
        if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower() == 'half':
            sliceRangeList = ['[%1.4f,%1.4f]' % (frontAxleXLoc + (-LREF*0.5),frontAxleXLoc+LREF*2),'[%1.4f,%1.4f]' % (frontAxleYLoc+(-WREF*0.75),frontAxleYLoc),'[%1.4f,%1.4f]' % (frontAxleZLoc,frontAxleZLoc + WREF*1)]
        else:
            sliceRangeList = ['[%1.4f,%1.4f]' % (frontAxleXLoc + (-LREF*0.5),frontAxleXLoc+LREF*2),'[%1.4f,%1.4f]' % (frontAxleYLoc+(-WREF*0.75),frontAxleYLoc+(WREF*0.75)),'[%1.4f,%1.4f]' % (frontAxleZLoc,frontAxleZLoc + WREF*1)]
    else:
        sliceRangeList = sliceDict['SLICE_RANGE'].split()


    print('\t\tRequested surface variables: %s' % (sliceVars))
    availableVars = list(varDict['sliceVariables'].keys())

    #check if inputs are of the same length
    if not len(normalsList) == len(nSliceList) == len(sliceRangeList) == len(sliceViews):
        sys.exit('ERROR! The inputs in pvPostSetup is not correct for [SLICE] -> number of inputs are not equal!')

    #start slice loop
    begin_slice = time.time()
    for normal, nslices, sliceRange,sliceView in zip(normalsList,nSliceList,sliceRangeList,sliceViews):
        sliceRange = np.array(sliceRange.replace('[','').replace(']','').split(',')).astype('float')
        sliceArray = np.round(np.linspace(sliceRange[0],sliceRange[1],int(nslices)),3)

        n = 0 # start counter
        for coord in sliceArray:
            slice1           = Slice(Input=volumeSource, SliceType="Plane" )
            slice1.SliceOffsetValues = 0
            if normal == 'X':
                slice1.SliceType.Normal  = [1,0,0]
                slice1.SliceType.Origin  = [coord,0,0]
            elif normal == 'Y':
                slice1.SliceType.Normal  = [0,1,0]
                slice1.SliceType.Origin  = [0,coord,0]
            elif normal == 'Z':
                slice1.SliceType.Normal  = [0,0,1]
                slice1.SliceType.Origin  = [0,0,coord]
            
            slice1.SliceType         = "Plane"
            slice1.Triangulatetheslice = 0
            #start looping through variables
            for variable in sliceVars:
                #if vorticity in internal volume do not mirror vectors

                #check if variables are available
                if not variable in availableVars:
                    print('\n\t\t%s not available, skipping...\n' % (variable))
                    continue

                #calculating variable using variable equation
                calculator = Calculator(registrationName='calculator', Input=slice1)
                calculator.Function = str(varDict['sliceVariables'][variable]['equation'])
                calculator.ResultArrayName = variable
                
                sliceSourceDisplay = Show(calculator,renderView,'UnstructuredGridRepresentation')
                sliceSourceDisplay.Representation = 'Surface'
                
                
                
                ColorBy(sliceSourceDisplay,('POINTS',variable))
                LUT = GetColorTransferFunction(variable)
                PWF  = GetOpacityTransferFunction(variable)
                title = varDict['sliceVariables'][variable]['label']
                varRange = np.array(varDict['sliceVariables'][variable]['range'])
                tableValues = 15
                color = varDict['sliceVariables'][variable]['color']

                LUT.NumberOfTableValues = tableValues #default
                LUT.RescaleTransferFunction(varRange[0],varRange[1])
                PWF.RescaleTransferFunction(varRange[0],varRange[1])
                LUT.ApplyPreset(color,True)
                PWF.ApplyPreset(color,True)
                renderView.Update()


                colorBar = GetScalarBar(LUT,renderView)
                colorBar.Title = title
                colorBar.TitleFontFamily = 'Times'
                colorBar.LabelFontFamily = 'Times'
                colorBar.ComponentTitle = ''
                colorBar.Orientation = 'Horizontal'
                colorBar.WindowLocation = 'Lower Center'
                if  'Cf' in variable:
                    colorBar.LabelFormat = '%-1.3g'
                    colorBar.RangeLabelFormat = '%-1.3g'
                    colorBar.CustomLabels = np.linspace(varRange[0],varRange[1],3)         
                else:
                    colorBar.LabelFormat = '%-1.1f'
                    colorBar.RangeLabelFormat = '%-1.1f'
                    colorBar.CustomLabels = np.linspace(varRange[0],varRange[1],5)         
                colorBar.ScalarBarLength = 0.3
                colorBar.ScalarBarThickness = 160
                colorBar.TitleFontSize = 160
                colorBar.LabelFontSize = 160
                colorBar.TitleColor = [0,0,0]
                colorBar.LabelColor = [0,0,0]
                colorBar.AddRangeLabels = 1
                colorBar.AutomaticLabelFormat = 0
                colorBar.DrawAnnotations = 0
                colorBar.DrawTickLabels = 1
                colorBar.UseCustomLabels = 1
                colorBar.BackgroundColor = [1,1,1,1]
                colorBar.BackgroundPadding = 10  
                       
                            

                HideUnusedScalarBars()  
                
                renderView = setView(renderView,viewsDict[sliceView])
                saveImages(renderView,caseName,variable,'slice',sliceView,normal=normal,position=coord,counter=n)
              
                
                try:
                    Delete(colorBar)
                except:
                    print('')
                Hide(sliceSourceDisplay,renderView)
                Delete(sliceSourceDisplay)
                Delete(calculator)

            n = n + 1


    
                

def setView(renderView,view):
    renderView.CameraParallelProjection = view['pp']
    renderView.CameraPosition           = view['campos']
    renderView.CameraFocalPoint         = view['focalpos']
    renderView.CameraViewUp             = view['viewup']
    renderView.CameraParallelScale      = view['parallelscale']
   
    return renderView

def saveImages(renderView,caseName,variable,imageType,view,normal=None,position=None,counter=None):
    #checking for image type, if not slice then omit last two variables in name
    if imageType.lower() != 'slice':
        fileName = '%s_%s_%s_%s.png' % (caseName,variable,imageType,view)
        dirName = '%s_%s' % (variable,imageType)
    else:
        fileName = '%s_%s_%s_%s_%s_%1.5g_%s.png' % (caseName,variable,imageType,view,normal,position,format(int(counter),"05"))
        dirName = '%s_%s' % (variable,imageType)


    #checking if required folders exists
    if os.path.isdir('postProcessing'):
        if os.path.isdir('postProcessing/images'):
            if not os.path.isdir('postProcessing/images/%s' % (dirName)):
                print('\t\tDirectory for %s not found, creating...' % (dirName))
                os.system('mkdir postProcessing/images/%s' % (dirName))
        else:
            print('\t\tDirectory for images not found, creating...')
            os.system('mkdir postProcessing/images')
            print('\t\tDirectory for %s not found, creating...' % (dirName))
            os.system('mkdir postProcessing/images/%s' % (dirName))
    else:
        print('\t\tpostProcessing directory not found, creating...')
        os.system('mkdir postProcessing')
        print('\t\tDirectory for images not found, creating...')
        os.system('mkdir postProcessing/images')
        print('\t\tDirectory for %s not found, creating...' % (dirName))
        os.system('mkdir postProcessing/images/%s' % (dirName))

    print('\t\t\tSaving screenshot: %s' % (fileName))
    SaveScreenshot(filename='postProcessing/images/' + dirName + '/' + fileName,viewOrLayout=renderView,OverrideColorPalette = 'WhiteBackground')
        

def generateDefaultViews(LREF,CREF,FREF,WREF):
    print('\t\tCreating default views...')

    #create point of focus for front, side and rear views
    focalYZ = np.array(CREF)
    WREF = float(WREF)
    

    A = (7*WREF)/(5*(1-(800/RESOLUTION[1])))
    a = 7*WREF/5
    dy = (A/2)-(a/2)


    Af = 5*(WREF/3)/(1-800/RESOLUTION[1])
    df = Af/2-(800*Af/RESOLUTION[1])
    df = df-(df*0.2)

    focalYZ[2] = df
    

    viewsDict = {'Rear':{'pp':1,'campos':focalYZ + [LREF*10,0,0],'focalpos':focalYZ,'viewup':[0,0,1],'parallelscale':Af/2.5},
                 'Front':{'pp':1,'campos':focalYZ + [LREF*-10,0,0],'focalpos':focalYZ,'viewup':[0,0,1],'parallelscale':Af/2.5},
                 'Left':{'pp':1,'campos':focalYZ + [LREF*0.2,LREF*-10,0],'focalpos':focalYZ + [LREF*0.2,0,0],'viewup':[0,0,1],'parallelscale':Af/2.2},
                 'Right':{'pp':1,'campos':focalYZ + [LREF*0.2,LREF*10,0],'focalpos':focalYZ + [LREF*0.2,0,0],'viewup':[0,0,1],'parallelscale':Af/2.2},
                 'Top':{'pp':1,'campos':focalYZ + [0,-dy,LREF*10],'focalpos':focalYZ + [0,-dy,0],'viewup':[0,1,0],'parallelscale':A/1.8},
                 'Bottom':{'pp':1,'campos':focalYZ + [0,dy,-LREF*10],'focalpos':focalYZ + [0,dy,0],'viewup':[0,-1,0],'parallelscale':A/1.8},
                 'FrontLeft':{'pp':1,'campos':focalYZ + [-LREF*5,-LREF*5,LREF*1.5],'focalpos':focalYZ + [CREF[0]*0,0,-df/2],'viewup':[0,0,1],'parallelscale':Af/2},
                 'FrontRight':{'pp':1,'campos':focalYZ + [-LREF*5,LREF*5,LREF*1.5],'focalpos':focalYZ + [CREF[0]*0,0,-df/2],'viewup':[0,0,1],'parallelscale':Af/2},
                 'RearLeft':{'pp':1,'campos':focalYZ + [LREF*5,-LREF*5,LREF*1.5],'focalpos':focalYZ + [CREF[0]*0,0,-df/2],'viewup':[0,0,1],'parallelscale':Af/2},
                 'RearRight':{'pp':1,'campos':focalYZ + [LREF*5,LREF*5,LREF*1.5],'focalpos':focalYZ + [CREF[0]*0,0,-df/2],'viewup':[0,0,1],'parallelscale':Af/2},
                 }

    for key in viewsDict.keys():
        print('\t\t\t%s' % (key))
        for item in viewsDict[key].keys():
                
            val = viewsDict[key][item]
            
            print('\t\t\t\t%s: %s' % (item,str(val)))

    return viewsDict




def getPVSetup(templateLoc):
    '''
        Looks for the setup config file to know the settings
        - if the file cannot be found in CWD it will use the default setup file
        - will check to see if all the required variables exists, if not, it will
        populate based on the default config file
        
    '''
    print('\tGetting pvPostSetup file...')
    pvSetupList = glob.glob('pvPostSetup')
    if len(pvSetupList) == 0:
        print('\t\tUnable to find pvPostSetup file! Using default settings...')
        pvSetupPath = '%s/%s/pvPostSetup' % (templateLoc,'default')
        

    pvPostSetupConfig = configparser.ConfigParser()
    pvPostSetupConfig.optionxform = str
    
    try:
        pvPostSetupConfig.read_file(open(pvSetupPath))
    except Exception as e:
        print('\n\n\nERROR! Unable to get the pvPostSetup!')
        print(e)
        sys.exit()

    print('\n\t\tpvPostSetup settings:')
    for key in pvPostSetupConfig.keys():
        if key == 'DEFAULT':
            continue
        else:        
            print('\t\t\t%s' % (key))
            for item in pvPostSetupConfig[key].keys():
                value = pvPostSetupConfig[key][item]
                print('\t\t\t\t%s: %s'% (item,value))

    
    pvPostSetupDict = pvPostSetupConfig
    
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
    WREF = float(caseSetupConfig['BC_SETUP']['REFWIDTH'])
    CREF = np.array(caseSetupConfig['BC_SETUP']['REFCOR'].split(' ')).astype('float')
    FREF = np.array(caseSetupConfig['BC_SETUP']['FRT_WH_CTR'].split(' ')).astype('float')

    print('\t\tReference Values:')
    print('\t\t\tUREF: %1.3g, LREF: %1.3g, WREF: %1.3g, CREF: [%1.3g %1.3g %1.3g] , FREF: [%1.3g %1.3g %1.3g]' % (UREF,LREF,WREF,CREF[0],CREF[1],CREF[2],FREF[0],FREF[1],FREF[2]))
    
    return UREF,LREF,WREF,CREF,FREF,caseSetupConfig

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
        selections['/Root/{}'.format(source.GetDataInformation().GetBlockName(i).replace('-',''))] = counter
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
        if any(x.replace('-','') in selector for x in include_patches):
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












