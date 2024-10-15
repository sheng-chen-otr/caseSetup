import os
import sys
import numpy as np
import subprocess as sp
import math
from utilities import *
from writeSystem import *


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
            