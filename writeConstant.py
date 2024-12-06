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
    turbulencePropertiesPath = '%s/defaultDicts/constant/turbulenceProperties' % (templateLoc)
    localTurbulencePropertiesiesPath = 'constant/turbulenceProperties' 
    copyTemplateToCase(turbulencePropertiesPath,localTurbulencePropertiesiesPath)
    
    #writing viscosity
    simtype = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_TYPE'][0]
    turbmodel = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['TURB_MODEL'][0]
    
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
        
    
    
    search_and_replace(localTurbulencePropertiesiesPath,'<SIMTYPE>',sim_type)
    search_and_replace(localTurbulencePropertiesiesPath,'<TURBMODEL>',turb_model)
    
   
        
    
    
    
    #nu must be specified always, if using default, put default
    

    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
            