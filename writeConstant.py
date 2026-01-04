import os
import sys
import numpy as np
import subprocess as sp
import math
from utilities import *
from writeSystem import *

propertiesDict = {'air':{'nu':1.51e-5},
                  'custom':{'nu':1.51e-5}}


def writeOptions(templateLoc, geomDict,fullCaseSetupDict):
    print('\tWriting fvOptions...')
    fvOptionPath = '%s/defaultDicts/constant/fvOptions' % (templateLoc)
    localFvOptionPath = 'constant/fvOptions' 
    copyTemplateToCase(fvOptionPath,localFvOptionPath)
    
    optionList = []
    #write porous media options
    
    for geom in geomDict:
        if geom.startswith('POR'):
            
            geomName = geom.split('.')[0]
            print('\t\tWriting out porous media: %s' % (geomName))
            porString = '''GEOMNAME{type explicitPorositySource;active yes; explicitPorositySourceCoeffs {type DarcyForchheimer; selectionMode cellZone;cellZone GEOMNAME_INTERNAL; DarcyForchheimerCoeffs{d d [0 -2 0 0 0 0 0] (DCOEFFS 3e10 3e10);f f [0 -1 0 0 0 0 0] (FCOEFFS 1e5 1e5); coordinateSystem {type cartesian; origin  (0 0 0); coordinateRotation {type axesRotation;e1 (VEC1);e2 (VEC2);}}}}}\n'''
            
            string = porString.replace('GEOMNAME',geomName)
            for key in fullCaseSetupDict[geomName].keys():
                
                string = string.replace(key,' '.join(fullCaseSetupDict[geomName][key]))
            
            optionList.append(string)
    
    
    optionList = ''.join(optionList)
    
    search_and_replace(localFvOptionPath,'<OPTION_MEDIA>',optionList)

def writeMRFG(templateLoc, geomDict,fullCaseSetupDict):
    print('\tWriting MRF Regions...')
    MRFPropertiesPath = '%s/defaultDicts/constant/MRFProperties' % (templateLoc)
    localMRFPropertiesPath = 'constant/MRFProperties' 
    copyTemplateToCase(MRFPropertiesPath,localMRFPropertiesPath)
    
    mrfList = []
    #write mrf options
    
    for geom in geomDict:
        if geom.startswith('MRFG'):
            
            geomName = geom.split('.')[0]

            mrfRadDefault = fullCaseSetupDict[geomName]['MRF_RAD'][0].lower()
            mrfCenterDefault = fullCaseSetupDict[geomName]['MRF_CENTER'][0].lower()
            mrfAxisDefault = fullCaseSetupDict[geomName]['MRF_AXIS'][0].lower()

            if 'default' in [mrfRadDefault,mrfCenterDefault,mrfAxisDefault]:
                try:
                    vertices, faces = readGeomFile(geom)
                    xcenter,ycenter,zcenter, xaxis,yaxis,zaxis = find_wheel_axis(vertices,faces)
                    radius = zcenter - float(fullCaseSetupDict['BC_SETUP']['FRT_WH_CTR'][2]) #wheel contact patch intersection with ground
                except Exception as error:
                    print('ERROR! Unable to calculate coordinates for %s, please check that your geometry is valid or manually input coordinates!' % (geom))
                    print('%s' % (error))
                    sys.exit()
                whOrig = '%1.6f %1.6f %1.6f' % (xcenter, ycenter, zcenter) #sets the wheel origin based on bounding box
                whAxis = '%1.2g %1.2g %1.2g' % (xaxis,yaxis,zaxis)
                if fullCaseSetupDict[geom.split('.')[0]]['MRF_RAD'][0].lower() != 'default':
                        try:
                            radius = float(fullCaseSetupDict[geom.split('.')[0]]['WH_RAD'][0])
                        except Exception as error:
                            print('\t\t\t\tWARNING! Value given for radius in %s is not valid, using calculated value instead.' % (geom.split('.')[0]))
                            print('\t\t\t\t%s' % (error))
                if fullCaseSetupDict[geom.split('.')[0]]['MRF_AXIS'][0].lower() != 'default':
                    try:
                        whAxisX = float(fullCaseSetupDict[geom.split('.')[0]]['WH_AXIS'][0])
                        whAxisY = float(fullCaseSetupDict[geom.split('.')[0]]['WH_AXIS'][1])
                        whAxisZ = float(fullCaseSetupDict[geom.split('.')[0]]['WH_AXIS'][2])
                        
                        whAxis = '%1.6f %1.6f %1.6f' % (whAxisX, whAxisY, whAxisZ)
                    except:
                        print('\t\t\t\tWARNING! Value given for wheel axis in %s is not valid, using calculated value instead.' % (geom.split('.')[0]))
                if fullCaseSetupDict[geom.split('.')[0]]['MRF_CENTER'][0].lower() != 'default':
                    try:
                        whCentX = float(fullCaseSetupDict[geom.split('.')[0]]['WH_CENTER'][0])
                        whCentY = float(fullCaseSetupDict[geom.split('.')[0]]['WH_CENTER'][1])
                        whCentZ = float(fullCaseSetupDict[geom.split('.')[0]]['WH_CENTER'][2])
                        whOrig = '%1.6f %1.6f %1.6f' % (whCentX, whCentY, whCentZ)
                    except Exception as error:
                        print('\t\t\t\tWARNING! Value given for wheel center in %s is not valid, using calculated value instead.' % (geom.split('.')[0]))
                        print('\t\t\t\t%s' % (error))
                        try:
                            vertices, faces = readGeomFile(geom)
                            xcenter,ycenter,zcenter, xaxis,yaxis,zaxis = find_wheel_axis(vertices,faces)
                            radius = zcenter
                            #bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ = getBoundingBox(geom.replace('.gz',''))
                        
                        except Exception as error:
                            print('ERROR! Unable to calculate coordinates for %s, please check that your geometry is valid or manually input coordinates!' % (geom))
                            print('%s' % (error))
                            sys.exit()
                            continue
                        whOrig = '%1.6f %1.6f %1.6f' % (xcenter, ycenter, zcenter) #sets the wheel origin based on bounding box
            

            

            print('\t\t\t\tWheel Center: %s' % (whOrig))
            print('\t\t\t\tWheel Axis: %s' % (whAxis))
            if fullCaseSetupDict[geom.split('.')[0]]['ROT_MRF'][0].lower() == 'true':
                rotaVel = calcRotaVel(float(fullCaseSetupDict['BC_SETUP']['INLET_MAG'][0]),radius)
                print('\t\t\t\tCalculated Radius: %1.4f m' % (radius))
                print('\t\t\t\tCalculated Radial Velocity: %1.4f rad/s' % (rotaVel))
                whVel = str(rotaVel)
                print('\t\tWriting out MRF region: %s' % (geomName))
                mrfString = '''toint-GEOM_NAME{cellZone fluid-GEOM_NAME;active yes;origin (CENTER_POINT);axis (ROT_AXIS);omega constant ROT_VEL;nonRotatingPatches ();}\n'''
                string = mrfString.replace('GEOM_NAME',geomName).replace('CENTER_POINT',whOrig).replace('ROT_AXIS',whAxis).replace('ROT_VEL',whVel)
                mrfList.append(string)
            elif fullCaseSetupDict[geom.split('.')[0]]['ROT_MRF'][0].lower() == 'false':
                rotaVel = 0
                mrfList.append('')
            else:
                sys.exit('ERROR! [%s] -> ROT_MRF input invalid!' % (geom.split('.')[0]))
    
    
    mrfList = ''.join(mrfList)
    
    search_and_replace(localMRFPropertiesPath,'<MRF_ZONES>',mrfList)
    
    
def writeTransportProperties(templateLoc, fullCaseSetupDict):
    print('\tWriting transportProperties...')
    transportPropertiesPath = '%s/defaultDicts/constant/transportProperties' % (templateLoc)
    localTransportPropertiesPath = 'constant/transportProperties' 
    copyTemplateToCase(transportPropertiesPath,localTransportPropertiesPath)
    
    #writing viscosity
    material = fullCaseSetupDict['GLOBAL_MATERIAL']['MATERIAL_TYPES']
    
    if len(material) > 1:
        sys.exit('ERROR! Multiple materials specified in [GLOBAL_MATERIAL], only single phase simulations available at this time! Please enter only 1 material.')
    elif len(material) < 1:
        sys.exit('ERROR! No materials specified in [GLOBAL_MATERIAL]! Please specify atleast 1 material!')
    else:
        material = material[0]
        print('\t\tMaterial: %s' % (material))
        
    nu = fullCaseSetupDict['GLOBAL_MATERIAL']['VISCOSITY']
    
    
    #nu must be specified always, if using default, put default
    if len(nu) != 1:
        sys.exit('ERROR! Viscosity in [GLOBAL_MATERIAL] must be specified, if intending to use default, enter ''default''!')
    if nu == 'default':
        nu = propertiesDict[material]
        
    print('\t\tViscosity = %s' % (str(nu[0])))        
    search_and_replace(localTransportPropertiesPath,'<VISCOSITY>',str(nu[0]))
    
    
def writeTurbulenceProperties(templateLoc, fullCaseSetupDict):
    print('\tWriting turbulenceProperties...')
    turbulencePropertiesSolvePath = '%s/defaultDicts/constant/turbulencePropertiesSolve' % (templateLoc)
    turbulencePropertiesInitPath = '%s/defaultDicts/constant/turbulencePropertiesInit' % (templateLoc)
    localTurbulencePropertiesSolvePath = 'constant/turbulencePropertiesSolve' 
    localTurbulencePropertiesInitPath = 'constant/turbulencePropertiesInit' 
    copyTemplateToCase(turbulencePropertiesSolvePath,localTurbulencePropertiesSolvePath)
    copyTemplateToCase(turbulencePropertiesInitPath,localTurbulencePropertiesInitPath)
    
    #writing viscosity
    simtype = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_TYPE'][0]
    turbmodel = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['TURB_MODEL'][0]
    initType = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0]
    
    if not simtype in dictDict['solverType'].keys():
        sys.exit('ERROR! Solver type is not valid in [GLOBAL_SIM_CONTROL]! Available solver types: %s' % (list(dictDict['solverType'].keys())))
    elif simtype == 'steady':
        sim_type = 'RAS'
    elif simtype == 'transient':
        sim_type = 'LES'

    
    
    if turbmodel in dictDict['solverType'][simtype].keys():
        turb_model = dictDict['solverType'][simtype][turbmodel]
        print('\t\tTurbulence Model: %s' % (turb_model))
    else:
        sys.exit('ERROR! Turbulence model is not valid in [GLOBAL_SIM_CONTROL]! Available turbulence models: %s' % (list(dictDict['solverType'][simtype].keys())))
        
    
    if initType == 'steady':
        search_and_replace(localTurbulencePropertiesInitPath,'<SIMTYPE>','RAS')
        search_and_replace(localTurbulencePropertiesInitPath,'<TURBMODEL>',dictDict['solverType']['steady'][turbmodel])
        
    search_and_replace(localTurbulencePropertiesSolvePath,'<SIMTYPE>',sim_type)
    search_and_replace(localTurbulencePropertiesSolvePath,'<TURBMODEL>',turb_model)
    
    
   
        
    
    
    
    #nu must be specified always, if using default, put default
    

    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
            