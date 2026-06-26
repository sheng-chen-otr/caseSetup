import os
import sys
import numpy as np
import subprocess as sp
import math
import re
from utilities import *


steadyTurb = {'sa':'SpalartAllmaras','kosst':'kOmegaSST'}

transientTurb = {'sa':'SpalartAllmarasIDDES','kosst':'kOmegaSSTIDDES'}

initDict = {'steady': {'potential':'controlDictPotential','none':''},
            'transient':{'potential':'controlDictPotential','steady':'controlDictSimpleInit','none':''}}
            
controlDictDict = {'steady':'controlDictSimple',
                   'transient':'controlDictPiso'}

exportDict = {'steady':'controlDictSimpleExport',
              'transient':'controlDictPisoExport'}

#SRF cornering controlDicts. steady -> SRFSimpleFoam, transient -> SRFPimpleFoam. createZeroDirectory
#picks the field-creation template from the controlDict `application`, so these drive both the solver
#and the Urel field creation.
controlDictCornerDict = {'steady':'controlDictSRFSimple',
                         'transient':'controlDictSRFPiso'}

exportCornerDict = {'steady':'controlDictSRFSimpleExport',
                    'transient':'controlDictSRFPisoExport'}

dictDict = {'solverType':{'steady':steadyTurb,'transient':transientTurb},
            'initType': initDict,
            'controlDictType':controlDictDict,
            'exportDictType':exportDict,
            'controlDictCornerType':controlDictCornerDict,
            'exportCornerType':exportCornerDict}

foList = ['averageFieldsDict','cptMeanDict','nearWallFieldsDict','wallShearStressDict','vorticityDict','QCriterionDict','yPlusDict','surfaceFieldAverage','surfaces']
coeffList = ['forceCoeffs','forceCoeffsExport','forceCoeffSetup']
prefixToIgnore = ['IDOM','SMP','REFX','REF','MRFG','POR'] #list of prefixes to not include in the boundary conditions for geometry as well as for forceCalculations

forceVecDict = {'drag':{'default':'1 0 0'},
                'lift':{'default':'0 0 1'},
                'pitch':{'default':'0 1 0'}}

def writeSurfaceFeatureExtract(templateLoc,geomDict,fullCaseSetupDict):
    print('\n\tWriting out to surfaceFeatureExtract...')
    print('\t\tGeometries to write:')
    
    templatePath = '%s/defaultDicts/system/surfaceFeatureExtractDict' % (templateLoc)
    localPath = 'system/surfaceFeatureExtractDict'
    copyTemplateToCase(templatePath,localPath)

        
    
    featExtractStringArray = []
    featExtractString = 'GEOM_NAME{extractionMethod extractFromSurface; includedAngle FEAT_ANGLE;subsetFeatures{nonManifoldEdges no; openEdges yes;}writeObj no;}'
    feat_angle = fullCaseSetupDict['GLOBAL_REFINEMENT']['FEAT_EDGE_ANGLE'][0]
    for geom in geomDict.keys():
        print('\t\t\t'+geom)
        featString = featExtractString.replace('GEOM_NAME',geom.replace('.gz',''))\
                                      .replace('FEAT_ANGLE',feat_angle)
        featExtractStringArray.append(featString)
    
    featExtractString = '\n'.join(featExtractStringArray)
    search_and_replace('system/surfaceFeatureExtractDict','<SURFACE_FEAT_EXTRACT>',featExtractString)
    

   
    
def writeControlDict(templateLoc, fullCaseSetupDict):
    print('\n\tWriting out controlDicts...')

    #get the control dict template paths
    setupDict = fullCaseSetupDict['GLOBAL_CONTROL']
    controlDictTemplatePath = "%s/defaultDicts/system" % (templateLoc)
    simType = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_TYPE'][0].lower()
    simInit = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0].lower()
    simTurb = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['TURB_MODEL'][0].lower()
    #SRF cornering is solved in the rotating frame: the relative velocity Urel is the primary field
    #and createZeroDirectory must emit it. A dedicated controlDict (application SRFSimpleFoam) drives
    #both the solver and the field creation, so cornering selects the SRF controlDict templates below.
    runCornering = ('CORNERING_SETUP' in fullCaseSetupDict and
                    fullCaseSetupDict['CORNERING_SETUP']['RUN_CORNERING'][0].lower() == 'true')
    if runCornering and simType not in dictDict['controlDictCornerType'].keys():
        print('\n\tERROR: SRF cornering is not supported for SIM_TYPE = %s.' % (simType))
        print('\tSupported cornering sim types are: %s' % (list(dictDict['controlDictCornerType'].keys())))
        sys.exit()
    if simType in dictDict['solverType'].keys():
        #sets up the controlDict for initialization run
        possibleInit = dictDict['initType'][simType].keys()
        
        #writes out controlDict for initialization method
        if runCornering:
            #cornering skips the absolute-frame initialisation entirely; writing controlDictPotential or
            #controlDictSimpleInit here would let the solve preamble swap controlDict and mislead
            #createZeroDirectory about the application, so no initialisation controlDict is written.
            print('\t\tCornering case: skipping initialisation controlDict')
        elif simInit in possibleInit:
            controlDictTemplate = dictDict['initType'][simType][simInit]
            initTemplatePath = '%s/%s' % (controlDictTemplatePath,controlDictTemplate)
            if controlDictTemplate != '':
                localControlDictPath = 'system/%s' % (controlDictTemplate)
                copyTemplateToCase(initTemplatePath,localControlDictPath)

                print('\t\tWriting controlDict for initialization: %s' % (simInit))
                for key in setupDict.keys():
                    search_and_replace(localControlDictPath, '<%s>' % (key), setupDict[key][0])
            else:
                print('\t\t\tNo initialization selected!')
        else:
            print('\n\t\tERROR: Unable to find initialization type in [GLOBAL_SIM_CONTROL] -> SIM_INIT = %s' % (simInit))
            print('\t\tPossible initialization methods for SIM_TYPE = %s is %s' % (simType,list(possibleInit)))
            print('\t\t\tIf trying to skip initialization, do SIM_INIT = none')
            sys.exit()
        #sets up the controlDict for actual run
        print('\t\tWriting controlDict for solver: %s' % (simType))
        if runCornering:
            controlDictTemplate = dictDict['controlDictCornerType'][simType]
            exportDictTemplate = dictDict['exportCornerType'][simType]
        else:
            controlDictTemplate = dictDict['controlDictType'][simType]
            exportDictTemplate = dictDict['exportDictType'][simType]
        solverPath = '%s/%s' % (controlDictTemplatePath,controlDictTemplate)
        exportPath = '%s/%s' % (controlDictTemplatePath,exportDictTemplate)
        localControlDictPath = 'system/%s' % (controlDictTemplate)
        localExportDictPath = 'system/controlDictExport'
        copyTemplateToCase(solverPath,localControlDictPath)
        copyTemplateToCase(exportPath,localExportDictPath)
        for key in setupDict.keys():
            search_and_replace(localControlDictPath, '<%s>' % (key), setupDict[key][0])
    else:
        print('\n\tERROR: Unable to find sim type in [GLOBAL_SIM_CONTROL] -> SIM_TYPE = %s' % (simType))
        print('\tPossible sim types are: %s' % (list(dictDict['solverType'].keys()))) 
        sys.exit()
    
def copyFunctionObjects(templateLoc,foList,geomDict,fullCaseSetupDict):
    print('\tCopying functionObjects...')
    for fo in foList:
        print('\t\t%s' % (fo))
        foLoc = '%s/defaultDicts/system/%s' % (templateLoc,fo)
        localPath = 'system/%s' % (fo)
        try:
            copyTemplateToCase(foLoc,localPath)
        except Exception as e:
            sys.exit('\t\t%s' % (e))
            
            continue
        for section in fullCaseSetupDict.keys():
            for var in fullCaseSetupDict[section].keys():
                search_and_replace(localPath, '<%s>' % (var), ' '.join(fullCaseSetupDict[section][var]))
        
        

def writeCorneringRelativeFields(fullCaseSetupDict):
    #SRF cornering solves the relative velocity Urel in the rotating frame. U is the absolute velocity
    #(Urel + omega x r) and is dominated by the solid-body rotation away from the rotation centre, so it
    #is not a meaningful field for post-processing a cornering case. Repoint every function object that
    #reads the absolute velocity / mean velocity at its relative-frame equivalent so all derived fields
    #(vorticity, Q, wall shear, yPlus, Cp/Cpt, surface samples, force coefficients) are computed from Urel.
    #Only the *input* velocity field references are swapped; result/output names (QMean, vorticityMean,
    #UnwMean, cpMean, ...) are kept so downstream sampling, pvPost and the report tooling are unaffected.
    #averageFieldsDict is deliberately left alone: it averages both U and Urel, so UrelMean/UrelPrime2Mean
    #already exist for these function objects to consume.
    runCornering = ('CORNERING_SETUP' in fullCaseSetupDict and
                    fullCaseSetupDict['CORNERING_SETUP']['RUN_CORNERING'][0].lower() == 'true')
    if not runCornering:
        return
    print('\n\tCornering case: repointing velocity-dependent function objects to Urel/UrelMean...')
    #UMean and the other derived velocity fields are never used as dictionary keywords, so they can be
    #swapped wherever they appear (input field references and sampled-field lists alike). Order longest
    #first so UMeanNear/UPrime2Mean are not partially matched by the bare UMean rule.
    meanSwaps = [(r'\bUPrime2Mean\b', 'UrelPrime2Mean'),
                 (r'\bUMeanNear\b', 'UrelMeanNear'),
                 (r'\bUMean\b', 'UrelMean')]
    #bare U is also used as a dictionary KEYWORD (e.g. the pressure / forceCoeffs `U  UMean;` entry), so
    #only swap bare U -> Urel in files where it appears purely as a field VALUE. The word-boundary regex
    #leaves Urel, UMean, UnwMean, UInf, magUInf untouched.
    valueFiles = ['controlDictSRFSimple', 'controlDictSRFPiso',
                  'vorticityDict', 'QCriterionDict', 'wallShearStressDict', 'yPlusDict',
                  'nearWallFieldsDict', 'surfaceFieldAverage', 'surfaces']
    keywordFiles = ['cptMeanDict', 'forceCoeffsExport']
    for fo in valueFiles + keywordFiles:
        localPath = 'system/%s' % (fo)
        if not os.path.exists(localPath):
            continue
        with open(localPath, 'r') as f:
            content = f.read()
        for pattern, repl in meanSwaps:
            content = re.sub(pattern, repl, content)
        if fo in valueFiles:
            content = re.sub(r'\bU\b', 'Urel', content)
        with open(localPath, 'w') as f:
            f.write(content)
        print('\t\t%s' % (fo))


def writeForceCoeff(templateLoc,geomDict,fullCaseSetupDict):

    #copying forceCoeff templates
    print('\tCopying forceCoeff templates...')
    for coeff in coeffList:
        print('\t\t%s' % (coeff))
        foLoc = '%s/defaultDicts/system/%s' % (templateLoc,coeff)
        localPath = 'system/%s' % (coeff)
        try:
            copyTemplateToCase(foLoc,localPath)
        except Exception as e:
            sys.exit('\t\t%s' % (e))
    #writing to the forceCoeffs and forceCoeffsExport
    allList = []
    geomList = []
    
    for geom in geomDict:
        if any(x in geom for x in prefixToIgnore):
            continue
        
        geomName = geom.split('.')[0]
        allList.append('''".*%s.*" ''' % (geomName))
        geomOut = """%s{patches (".*%s.*"); #include "forceCoeffSetup"}\n""" % (geomName,geomName)
        geomList.append(geomOut)
    
    allString = ''.join(allList)
    allLine = """all{patches (%s); #include "forceCoeffSetup"}\n""" % (allString)
    geomList.append(allLine)
    geomList = "".join(geomList)
    print("\tSetting up forceCoeffs...")
    try:
        search_and_replace("system/forceCoeffs", "<FORCE_PATCHES>",geomList)
    except Exception as e:
        sys.exit('ERROR! Unable to find forceCoeffs in system directory! %s' % (e))
        
        
    print("\tSetting up forceCoeffSetup and forceCoeffsExport...")
    try:
        search_and_replace("system/forceCoeffsExport" , "<FORCE_PATCHES>",allString)
        for key in fullCaseSetupDict['BC_SETUP'].keys():
            if 'drag' in key.lower() or 'lift' in key.lower() or 'pitch' in key.lower():
                if not '_VEC' in key:
                    continue
                forceType = key.split('_')[0].lower()
                print('\t\tSetting %s vector:' % (forceType.upper()))
                forceVec = fullCaseSetupDict['BC_SETUP'][key] #initial force vec, checking if single key word or a vector

                if forceType in forceVecDict.keys():
                    if forceVec[0].lower() in forceVecDict[forceType].keys():
                        forceKey = forceVec[0].lower()
                        forceVec = forceVecDict[forceType][forceVec[0].lower()]
                        print('\t\t\tusing %s: %s' % (forceKey,forceVec))
                        print('\t\t\ttransforming force vector to match domain...')
                        forceVec = transformForceVector(fullCaseSetupDict,forceVec)
                        print('\t\t\t\tupdated force vector: %s' % (forceVec))

                    else:
                        print('\t\t\tusing user input: %s' % (forceVec))
                        forceVec = ' '.join(forceVec)
                if isinstance(forceVec, list):
                    #forceType not in forceVecDict (e.g. a *_VEC key with an unrecognised prefix): forceVec
                    #is still the raw token list. search_and_replace needs a string, so join it here.
                    forceVec = ' '.join(forceVec)
                search_and_replace("system/forceCoeffsExport", "<%s>" % (key),forceVec) 
                search_and_replace("system/forceCoeffSetup", "<%s>" % (key),forceVec) 
            elif 'COR' in key:
                search_and_replace("system/forceCoeffsExport", "<%s>" % (key),' '.join(fullCaseSetupDict['BC_SETUP'][key])) 
                search_and_replace("system/forceCoeffSetup", "<%s>" % (key),' '.join(fullCaseSetupDict['BC_SETUP'][key])) 
            search_and_replace("system/forceCoeffsExport", "<%s>" % (key),fullCaseSetupDict['BC_SETUP'][key][0]) 
            search_and_replace("system/forceCoeffSetup", "<%s>" % (key),fullCaseSetupDict['BC_SETUP'][key][0])
    
    except Exception as e:
        sys.exit('ERROR! Unable to find forceCoeffsExport in system directory! %s' % (e))
    

    
def writeSurfaces(templateLoc,geomDict,fullCaseSetupDict):        
    #no need to copy surfacesTemplate due to it being copied by fo copy
    geomList = []
    geomOutList = []
    for geom in geomDict:
        if any(x in geom for x in prefixToIgnore):
            continue
        else:          
            geomName = geom.split('.')[0]
            geomList.append("""".*%s.*" """ % (geomName))
            geomOut = """\t%s{$patchSurface;patches (".*%s.*");}\n""" % (geomName,geomName)
            geomOutList.append(geomOut)

    geomString = ''.join(geomList)
    allLine = """\tall{$patchSurface;patches (%s);}\n""" % (geomString)
    geomOutList.append(allLine)
    geomOutList = "".join(geomOutList)

    print("\tSetting up surfaces...")
    search_and_replace("system/surfaces", "<SURFACE_PATCHES>",geomOutList)
    for var in fullCaseSetupDict['GLOBAL_POST_PRO'].keys():
        val = ''.join(fullCaseSetupDict['GLOBAL_POST_PRO'][var])
        search_and_replace("system/surfaces", "<%s>" % (var), val)
        
        
        

def writePostProSurfaceList(fullCaseSetupDict):
    #slices this would require updating the slices from the previous versions
    #might not work well with the pv-post output names
    #will run with defaults using track width and 
    print("\tWriting out surface list...")
    sliceList = ['XSLICE','YSLICE','ZSLICE']
    sliceDict = {}
    for sliceType in sliceList:
        sliceTypeName = "%s" % (sliceType)
        slicePositions = fullCaseSetupDict['GLOBAL_POST_PRO'][sliceType][0]
        if slicePositions == 'default':

            #calculate the default values for slices:
            #for x normals it's 1.5 times REFLEN ahead of the center and 3 times REFLEN behind
            #for y normals it's in 2/3 times REFLEN in -y and 2/3 times REFLEN in +y-max
            #for z normals it's from 0.001m to 2/3 times REFLEN up
            
            #getting the wheel base
            reflen = float(fullCaseSetupDict['BC_SETUP']['REFLEN'][0])
            #get the reference center rotation
            refcen = np.array(fullCaseSetupDict['BC_SETUP']['REFCOR']).astype(float)

            
            #calculating start and end for direction

            if 'X' in sliceType:
                start = refcen[0] - (reflen*1.5)
                end = refcen[0] + (reflen*3)
                ds = 0.1
            elif 'Y' in sliceType:
                start = refcen[1] - (reflen*(2/3))
                end = refcen[1] + (reflen*(2/3))
                ds = 0.1
            elif 'Z' in sliceType:
                start = refcen[2] + 0.001
                end = refcen[2] + (reflen*(2/3))
                ds = 0.1
                
            start = round(start,3)
            end = round(end,3)
            ds = round(ds,3)
            sliceDict[sliceType] = {}
            sliceDict[sliceType]['start'] = start
            sliceDict[sliceType]['end'] = end
            sliceDict[sliceType]['ds'] = ds
            sliceDict[sliceType]['nslice'] = int((end - start)/ds)
            
            print('\t\t%s' % (sliceType))
            print('\t\t\tStart: %s'% (start))
            print('\t\t\tEnd: %s' % (end))
            print('\t\t\tDs: %s' % (ds))
            print('\t\t\tnSlices: %s' % (int((end-start)/ds)))
        
  
        
        else:
            try:
                sliceDict[sliceType] = {}
                sliceDict[sliceType]['start'] = slidePositions[0]
                sliceDict[sliceType]['end'] = slidePositions[1]
                sliceDict[sliceType]['ds'] = slidePositions[2]
                sliceDict[sliceType]['nslice'] = int((slidePositions[1]-slidePositions[0])/slidePositions[2])
                
                print('\t%s' % (sliceType))
                print('\t\tStart: %s'% (slidePositions[0]))
                print('\t\tEnd: %s' % (slidePositions[1]))
                print('\t\tDs: %s' % (slidePositions[2]))
                print('\t\tnSlices: %s' % (int((slidePositions[1]-slidePositions[0])/slidePositions[2])))
            except:
                print('\t%s not requested or invalid input, skipping...' % (sliceType))
    
    with open('system/surfaceSetupList', 'w') as surfaceSetupList:
        surfacesList = []
        for sliceType in sliceDict.keys():
            sliceRange = np.linspace(sliceDict[sliceType]['start'],sliceDict[sliceType]['end'],sliceDict[sliceType]['nslice'])
            n = 0
            for i in sliceRange:
                if 'X' in sliceType:
                    basePoint = [i,0,0]
                    normalVector = [1,0,0]
                    normal = 'x'
                    #print(basePoint)
                elif 'Y' in sliceType:
                    basePoint = [0,i,0]
                    normalVector = [0,1,0]
                    normal = 'y'
                elif 'Z' in sliceType:
                    basePoint = [0,0,i]
                    normalVector = [0,0,1]
                    normal = 'z'
                stringToAdd = """%s_%sNormal_%04d{type cuttingPlane; planeType pointAndNormal; pointAndNormalDict{ point (%1.4f %1.4f %1.4f); normal (%s %s %s);} interpolate true;}\n""" % (sliceType,normal,float(n),basePoint[0],basePoint[1],basePoint[2],normalVector[0],normalVector[1],normalVector[2])
                n = n + 1
                surfacesList.append(stringToAdd)
        surfaceSetupList.write(''.join(surfacesList))
        
def writeSolution(templateLoc, fullCaseSetupDict):
    print('\tWriting fvSolution...')
    simType = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_TYPE'][0]
    templatePathRefLims = '%s/defaultDicts/system/refsAndLimitsIncludeDict' % (templateLoc) 
    templatePathPiso = '%s/defaultDicts/system/fvSolutionPiso' % (templateLoc) 
    templatePathSimple = '%s/defaultDicts/system/fvSolutionSimple' % (templateLoc)
    templatePathSRFSimple = '%s/defaultDicts/system/fvSolutionSRFSimple' % (templateLoc)
    templatePathSRFPiso = '%s/defaultDicts/system/fvSolutionSRFPiso' % (templateLoc)
    templatePathPotential = '%s/defaultDicts/system/fvSolutionPotential' % (templateLoc) 
    copyTemplateToCase(templatePathRefLims, 'system/refsAndLimitsIncludeDict')
    runCornering = ('CORNERING_SETUP' in fullCaseSetupDict and
                    fullCaseSetupDict['CORNERING_SETUP']['RUN_CORNERING'][0].lower() == 'true')
    if runCornering:
        #cornering solves Urel in the rotating frame; the SRF fvSolution adds the Urel solver and
        #relaxation entries. Initialisation is skipped, so no init fvSolution is written. Transient
        #cornering uses the PIMPLE-based SRF variant, steady uses the SIMPLE-based one.
        print('\t\tCornering case:')
        if simType.lower() == 'transient':
            print('\t\t\tCopying fvSolutionSRFPiso')
            copyTemplateToCase(templatePathSRFPiso, 'system/fvSolutionSRFPiso')
        else:
            print('\t\t\tCopying fvSolutionSRFSimple')
            copyTemplateToCase(templatePathSRFSimple, 'system/fvSolutionSRFSimple')
        return
    
    if simType.lower() == 'transient':
        #copying over the fvSchemesPiso
        print('\t\tIs transient case:')
        
        print('\t\t\tCopying fvSolutionPiso')
        
        copyTemplateToCase(templatePathPiso, 'system/fvSolutionPiso')       
        if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0].lower() == 'potential':
            print('\t\t\tCopying fvSolutionPotential')
            copyTemplateToCase(templatePathPotential, 'system/fvSolutionPotential')
        elif fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0].lower() == 'steady':   
            print('\t\t\tCopying fvSolutionSimple')
            copyTemplateToCase(templatePathSimple, 'system/fvSolutionSimple')
    elif simType.lower() == 'steady':
        print('\t\tIs steady case:')
        print('\t\t\tCopying fvSolutionSimple')

        copyTemplateToCase(templatePathSimple, 'system/fvSolutionSimple')
        
        if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0].lower() == 'potential':
            print('\t\t\tCopying fvSolutionPotential')
            copyTemplateToCase(templatePathPotential, 'system/fvSolutionPotential')
    else:
        sys.exit('ERROR! SIM_TYPE is invalid!')


def writeSchemes(templateLoc, fullCaseSetupDict):
    print('\tWriting fvSchemes...')
    simType = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_TYPE'][0]
    templatePathPiso = '%s/defaultDicts/system/fvSchemesPiso' % (templateLoc) 
    templatePathSimple = '%s/defaultDicts/system/fvSchemesSimple' % (templateLoc) 
    templatePathSRFSimple = '%s/defaultDicts/system/fvSchemesSRFSimple' % (templateLoc) 
    templatePathSRFPiso = '%s/defaultDicts/system/fvSchemesSRFPiso' % (templateLoc) 
    runCornering = ('CORNERING_SETUP' in fullCaseSetupDict and
                    fullCaseSetupDict['CORNERING_SETUP']['RUN_CORNERING'][0].lower() == 'true')
    if runCornering:
        #cornering needs Urel momentum schemes (div(phi,Urel), laplacian(nuEff,Urel), grad(Urel)).
        #Transient cornering uses the backward-ddt PIMPLE variant, steady the steadyState SIMPLE one.
        print('\t\tCornering case:')
        if simType.lower() == 'transient':
            print('\t\t\tCopying fvSchemesSRFPiso')
            copyTemplateToCase(templatePathSRFPiso, 'system/fvSchemesSRFPiso')
        else:
            print('\t\t\tCopying fvSchemesSRFSimple')
            copyTemplateToCase(templatePathSRFSimple, 'system/fvSchemesSRFSimple')
        return
    
    if simType.lower() == 'transient':
        #copying over the fvSchemesPiso
        print('\t\tIs transient case:')

        print('\t\t\tCopying fvSchemesPiso')
        
        copyTemplateToCase(templatePathPiso, 'system/fvSchemesPiso')
        
        if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0].lower() == 'steady':   
            print('\t\t\tCopying fvSchemesSimple')
            copyTemplateToCase(templatePathSimple, 'system/fvSchemesSimple')

    elif simType.lower() == 'steady':
        print('\t\tIs steady case:')
        print('\t\t\tCopying fvSchemesSimple')

        copyTemplateToCase(templatePathSimple, 'system/fvSchemesSimple')
        
    else:
        sys.exit('ERROR! SIM_TYPE is invalid!')
        
    

def writeBoundaries(templateLoc,geomDict,fullCaseSetupDict):
    print('\n\tWriting out caseProperties...')
    bcStrings = {'GEOM':'GEOM_NAME{category wall; type noSlip; patches ("GEOM_NAME.*"); options {wallFunction WALL_MODEL; motion stationary;} values {$:initialConditions;}}',
                 'ROTA':'GEOM_NAME{category wall; type noSlip; patches ("GEOM_NAME.*"); options {wallFunction WALL_MODEL; motion rotating;} values {type rotatingWallVelocity;origin (WH_ORIG); axis (WH_AXIS); rotVel WH_VEL; $:initialConditions;}}',
                 'MOVG':'GEOM_NAME{category wall; type noSlip; patches ("GEOM_NAME.*"); options {wallFunction WALL_MODEL; motion stationary;} values {$:initialConditions;}}'
                }
    initDict = {'U':{'default':'uniform (INLET_VECTOR);'},
                'p':{'default':'p uniform 0;','userinput':'p uniform INIT_P;','csname':'INIT_P'},
                'k':{'default':'k uniform 0.24;','userinput':'k uniform INIT_K;','csname':'INIT_K'},
                'omega':{'default':'omega uniform 1.78;','userinput':'omega uniform INIT_OMEGA;','csname':'INIT_OMEGA'},
                'nut':{'default':'nut uniform 0;','userinput':'nut uniform INIT_NUT;','csname':'INIT_NUT'},
                'nuTilda':{'default':'nuTilda uniform 0.05;','userinput':'nuTilda uniform INIT_NUTILDA;','csname':'INIT_NUTILDA'}
                }
    
    initStringArray = []
    #creating initial conditions
    print('\t\tCreating initial conditions:')
    inletMag = fullCaseSetupDict['BC_SETUP']['INLET_MAG']
    yaw = fullCaseSetupDict['BC_SETUP']['YAW']
    sim_sym = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM']
    pitch = float(fullCaseSetupDict['BC_SETUP']['DOMAIN_PITCH'][0])
    if sim_sym[0].lower() == 'half':
        inletVec = velVector(float(inletMag[0]),0,pitch)
    else:
        inletVec = velVector(float(inletMag[0]),float(yaw[0]),pitch)
    velString = 'U uniform (%s);' % (inletVec[0])
    initStringArray.append(velString)   
    #SRF cornering solves the relative velocity Urel, which createZeroDirectory generates as a field.
    #createZeroDirectory sets each field's internalField to ${:initialConditions.<fieldName>}, so the
    #initialConditions block must contain a Urel entry or 0/Urel generation aborts. Seed Urel with the
    #same freestream vector as U (a reasonable initial guess; the SRFFreestream/SRFVelocity BCs and the
    #solver correct it for the rotating frame).
    runCornering = ('CORNERING_SETUP' in fullCaseSetupDict and
                    fullCaseSetupDict['CORNERING_SETUP']['RUN_CORNERING'][0].lower() == 'true')
    if runCornering:
        urelString = 'Urel uniform (%s);' % (inletVec[0])
        initStringArray.append(urelString)
        #SRFFreestreamVelocity reads UInf as a single vector keyword (no 'uniform' prefix), unlike the
        #'value'/internalField fields which are uniform Fields. Expose a bare-vector UInf in initialConditions
        #so the srfFreestream boundary template can reference ${:VALUE.UInf} without a 'uniform' token.
        #UInf is the freestream in the INERTIAL frame: SRFFreestreamVelocity imposes a patch value of
        #(UInf rotated into the frame) - omega x r. A car cornering through still air has zero inertial-frame
        #freestream; all of its relative wind comes from the rotation (omega x r, which equals INLET_MAG at the
        #car's radius). Seeding UInf with the INLET_MAG vector double-counts the speed (~2x velocity, ~4x forces)
        #and piles the velocity/pressure maxima onto the car, so UInf must be the zero vector for cornering.
        uInfString = 'UInf (0 0 0);'
        initStringArray.append(uInfString)
    
    
    for key in initDict.keys():
        if key != 'U':
            initName = initDict[key]['csname']
            initVal = fullCaseSetupDict['GLOBAL_SIM_CONTROL'][initName][0]
            
            if initVal.lower() == 'default':
                initString = initDict[key]['default']
                initStringArray.append(initString)
                print('\t\t\t%s: %s' % (key,initVal)) 
            else:
                try:
                    #testing if the value is valid, if it's not floatable then throw error!
                    initVal = float(initVal)
                except:
                    print('ERROR! %s does not have a valid input! Must be float or int or default!' % (initName))
                    sys.exit()
                initString = initDict[key]['userinput'].replace(initName,str(initVal))
                initStringArray.append(initString)
                print('\t\t\t%s: %1.4f' % (key,initVal))
    initStringArray = '\n'.join(initStringArray)
    print('\t\tSetting up boundary conditions for geometry...')
    #set up the domain boundaries first
    domainWallDictFull = {'inlet':'inlet{category inlet; type subSonic; patches (".*x-min.*" ".*y-min.*"); options {flowSpecification fixedVelocity;} values {$:initialConditions;}}',
                          'ground':'ground{category wall; type noSlip; patches (".*z-min.*"); options {wallFunction highReynolds; motion GROUND_MOV;} values {$:initialConditions;}}',
                          'outlet':'outlet{category outlet; type subSonic; patches (".*x-max.*" ".*y-max.*"); options {returnFlow default;} values {$:initialConditions;}}',
                          'walls':'walls{category wall; type slip; patches (".*z-max.*"); values {$:initialConditions;}}'
                          }
    domainWallDictHalf = {'inlet':'inlet{category inlet; type subSonic; patches (".*x-min.*"); options {flowSpecification fixedVelocity;} values {$:initialConditions;}}',
                          'ground':'ground{category wall; type noSlip; patches (".*z-min.*"); options {wallFunction highReynolds; motion GROUND_MOV;} values {$:initialConditions;}}',
                          'outlet':'outlet{category outlet; type subSonic; patches (".*x-max.*"); options {returnFlow default;} values {$:initialConditions;}}',
                          'walls':'walls{category wall; type slip; patches (".*z-max.*" ".*y-min.*"); values {$:initialConditions;}}',
                          'symmetry':'symmetry{category symmetry; type symmetry; patches (".*y-max.*");}'
                          }
    #SRF cornering: whole domain rotates with the car. All vertical (x/y) walls become self-selecting
    #SRFFreestream patches, the ground becomes an SRFVelocity moving wall, and the roof stays slip. No
    #symmetry plane (the cornering flow field is not laterally symmetric, so the case must be run full).
    domainWallDictCorner = {'inlet':'inlet{category inlet; type freestream; patches (".*x-min.*" ".*x-max.*" ".*y-min.*" ".*y-max.*"); options {flowSpecification srfFreestream;} values {$:initialConditions;}}',
                            'ground':'ground{category wall; type noSlip; patches (".*z-min.*"); options {wallFunction highReynolds; motion srf;} values {$:initialConditions;}}',
                            'walls':'walls{category wall; type slip; patches (".*z-max.*"); values {$:initialConditions;}}'
                            }
    internalDomainDict = {'inlet':'NAME{category inlet; type subSonic; patches (PATCHNAME); options {flowSpecification fixedVelocity;} values {INITIAL_CONDITIONS}}',
                          'movingwall':'NAME{category wall; type noSlip; patches PATCHNAME); options {wallFunction highReynolds; motion moving;} values {INITIAL_CONDITIONS}}',
                          'outlet':'NAME{category outlet; type subSonic; patches (PATCHNAME); options {returnFlow default;} values {$:initialConditions;}}',
                          'wall':'NAME{category wall; type slip; patches (PATCHNAME); values {$:initialConditions;}}',
                          'noslipwall':'NAME{category wall; type noSlip; patches (PATCHNAME); options {wallFunction highReynolds; motion stationary;} values {$:initialConditions;}}',
                          'symmetry':'NAME{category symmetry; type symmetry; patches (PATCHNAME);}'
                          }
    internalDomain = False
    for geom in geomDict.keys():
        if 'IDOM' in geom:
            internalDomain = True
    runCornering = ('CORNERING_SETUP' in fullCaseSetupDict and
                    fullCaseSetupDict['CORNERING_SETUP']['RUN_CORNERING'][0].lower() == 'true')
    caseSetupStringArray  = []    
    if internalDomain == False:
        if runCornering:
            domainWallDict = domainWallDictCorner
        elif fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0].lower() == 'half':
            domainWallDict = domainWallDictHalf

        else:
            domainWallDict = domainWallDictFull
        
        
        for patch in domainWallDict.keys():
            print('\t\t\t%s' % (patch))
            
            patchString = domainWallDict[patch].replace('GROUND_MOV',fullCaseSetupDict['BC_SETUP']['GROUND'][0])
            caseSetupStringArray.append(patchString)
    else:
        domainWallDict = internalDomainDict
        for section in fullCaseSetupDict.keys():
            if 'IDOM' in section:
                internalPatchName = section
                patchType = fullCaseSetupDict[section]['PATCH_TYPE'][0]
                if patchType.lower() == 'inlet':
                    if fullCaseSetupDict[section]['INLET_VEL'][0].lower() == 'default':
                        initialConditions = '$:initialConditions;'
                    else:
                        inletVel = fullCaseSetupDict[section]['INLET_VEL']
                        inletVelX = inletVel[0]
                        inletVelY = inletVel[1]
                        inletVelZ = inletVel[2]
                        initialConditions = 'U uniform (%1.4f %1.4f %1.4f); p uniform 0; k uniform 0.24; omega uniform 1.78; nut uniform 0; nuTilda uniform 0.05;' % (inletVelX,inletVelY, inletVelZ)
                    patchString = domainWallDict['inlet'].replace('PATCHNAME',"%s.*" % (section))\
                                                         .replace('NAME','%s' % (section))\
                                                         .replace('INITIAL_CONDITIONS',initialConditions)
                elif patchType.lower() == 'movingwall':
                        if fullCaseSetupDict[section]['PATCH_MOV'][0].lower() == 'default':
                            initialConditions = '$:initialConditions;'
                        else:
                            patchVel = fullCaseSetupDict[section]['PATCH_MOV']
                            inletVelX = patchVel[0]
                            inletVelY = patchVel[1]
                            inletVelZ = patchVel[2]
                            initialConditions = 'U uniform (%1.4f %1.4f %1.4f); p uniform 0; k uniform 0.24; omega uniform 1.78; nut uniform 0; nuTilda uniform 0.05;' % (inletVelX,inletVelY, inletVelZ)
                        patchString = domainWallDict['movingwall'].replace('PATCHNAME',"%s.*" % (section))\
                                                                .replace('NAME',"%s" % (section))\
                                                                .replace('INITIAL_CONDITIONS',initialConditions)
                elif patchType.lower() == 'outlet':
                        patchString = domainWallDict['outlet'].replace('PATCHNAME',"%s.*" % (section))\
                                                              .replace('NAME',"%s" % (section))
                elif patchType.lower() == 'wall':
                        patchString = domainWallDict['wall'].replace('PATCHNAME',"%s.*" % (section))\
                                                            .replace('NAME',"%s" % (section))
                elif patchType.lower() == 'noslipwall':
                        patchString = domainWallDict['noslipwall'].replace('PATCHNAME',"%s.*" % (section))\
                                                            .replace('NAME',"%s" % (section))
                elif patchType.lower() == 'symmetry':
                        patchString = domainWallDict['symmetry'].replace('PATCHNAME',"%s.*" % (section))\
                                                                .replace('NAME',"%s" % (section))
                else:
                    print('ERROR! Invalid PATCH_TYPE in %s' % (section))
                    sys.exit('Valid patch types are: %s' % (domainWallDict.keys()))
                
                caseSetupStringArray.append(patchString)
                

    wallFunctionDict = {'high':'highReynolds',
                        'low':'lowReynolds',
                        'ally':'allReynolds'
                        }
    for geom in geomDict.keys():
        print('\t\t\t%s' % (geom.split('.')[0]))
        geomPrefix = geom.split('-')[0]
        geomWallModel = geomDict[geom]['wallmodel'] #put in the wall model for the geometry
        
        availableWallFunctions = wallFunctionDict.keys()
        try:
            geomWallModel = wallFunctionDict[geomWallModel.lower()]
        except:
            print('ERROR! Wall model not valid, available wall models:')
            for model in availableWallFunctions:
                print(model)
            sys.exit()
        whAxis = '0 1 0' #initialize the default values
        whVel = '' 
        whOrig = ''
        
        if any(x in geom for x in prefixToIgnore):
            continue

        if geomPrefix in bcStrings.keys():
        
            geomString = bcStrings[geomPrefix]
            #sets the axis for rotational bodies, bodies that require rotating wall conditions
            if geomPrefix == 'ROTA' and fullCaseSetupDict[geom.split('.')[0]]['ROT_WH'][0].lower() == 'true':
                #checks if any of the values are default, if any are default, requires that bounding calculations are done
                defaultRad = fullCaseSetupDict[geom.split('.')[0]]['WH_RAD'][0].lower() == 'default'
                defaultAxis = fullCaseSetupDict[geom.split('.')[0]]['WH_AXIS'][0].lower() == 'default'
                defaultOrig = fullCaseSetupDict[geom.split('.')[0]]['WH_CENTER'][0].lower() == 'default'
                #if any are true, then run the bounding box calculation
                if True in [defaultRad,defaultAxis,defaultOrig]:
                    
                    try:
                        vertices, faces = readGeomFile(geom)
                        xcenter,ycenter,zcenter, xaxis,yaxis,zaxis = find_wheel_axis(vertices,faces)
                        #effective rolling radius = axle to the GROUND PLANE at the contact patch below the axle.
                        #the tyre penetrates the floor and is cut at the ground by snappyHexMesh, so the rolling-wall
                        #contact is at the ground plane (not the tyre geometry). calcLoadedRadius intersects the vertical
                        #line through the axle with the heave-corrected, pitch/roll-tilted floor plane.
                        radius = calcLoadedRadius(xcenter, ycenter, zcenter, fullCaseSetupDict) #axle to ground plane below axle
                        #bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ = getBoundingBox(geom.replace('.gz',''))
                        
                    except Exception as error:
                        print('ERROR! Unable to calculate coordinates for %s, please check that your geometry is valid or manually input coordinates!' % (geom))
                        print('%s' % (error))
                        sys.exit()
                        continue
                            
                    #xcenter,ycenter,zcenter,radius = getRotaCoordinates(bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ)
                    whOrig = '%1.6f %1.6f %1.6f' % (xcenter, ycenter, zcenter) #sets the wheel origin based on bounding box
                    whAxis = '%1.2g %1.2g %1.2g' % (xaxis,yaxis,zaxis)
                if fullCaseSetupDict[geom.split('.')[0]]['WH_RAD'][0].lower() != 'default':
                    try:
                        radius = float(fullCaseSetupDict[geom.split('.')[0]]['WH_RAD'][0])
                    except Exception as error:
                        print('\t\t\t\tWARNING! Value given for radius in %s is not valid, using calculated value instead.' % (geom.split('.')[0]))
                        print('\t\t\t\t%s' % (error))
                if fullCaseSetupDict[geom.split('.')[0]]['WH_AXIS'][0].lower() != 'default':
                    try:
                        whAxisX = float(fullCaseSetupDict[geom.split('.')[0]]['WH_AXIS'][0])
                        whAxisY = float(fullCaseSetupDict[geom.split('.')[0]]['WH_AXIS'][1])
                        whAxisZ = float(fullCaseSetupDict[geom.split('.')[0]]['WH_AXIS'][2])
                        
                        whAxis = '%1.6f %1.6f %1.6f' % (whAxisX, whAxisY, whAxisZ)
                    except:
                        print('\t\t\t\tWARNING! Value given for wheel axis in %s is not valid, using calculated value instead.' % (geom.split('.')[0]))
                if fullCaseSetupDict[geom.split('.')[0]]['WH_CENTER'][0].lower() != 'default':
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
                            radius = calcLoadedRadius(xcenter, ycenter, zcenter, fullCaseSetupDict) #axle to ground plane below axle
                            #bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ = getBoundingBox(geom.replace('.gz',''))
                        
                        except Exception as error:
                            print('ERROR! Unable to calculate coordinates for %s, please check that your geometry is valid or manually input coordinates!' % (geom))
                            print('%s' % (error))
                            sys.exit()
                            continue
                        whOrig = '%1.6f %1.6f %1.6f' % (xcenter, ycenter, zcenter) #sets the wheel origin based on bounding box
   
                print('\t\t\t\tWheel Center: %s' % (whOrig))
                print('\t\t\t\tWheel Axis: %s' % (whAxis))
                if fullCaseSetupDict[geom.split('.')[0]]['ROT_WH'][0].lower() == 'true':
                    if runCornering:
                        #cornering: local ground speed differs per wheel. In the rotating frame the wheel
                        #rolls at V_local = |omega_corner| * (horizontal distance from the corner axis to the
                        #wheel centre), so inner wheels spin slower and outer wheels faster than INLET_MAG.
                        omegaSigned, cAxis, cCentre = corneringFrame(fullCaseSetupDict)
                        whOrigParts = whOrig.split()
                        dx = float(whOrigParts[0]) - cCentre[0]
                        dy = float(whOrigParts[1]) - cCentre[1]
                        localSpeed = abs(omegaSigned) * math.sqrt(dx * dx + dy * dy)
                        print('\t\t\t\tCornering Local Ground Speed: %1.4f m/s' % (localSpeed))
                        rotaVel = calcRotaVel(localSpeed,radius)
                    else:
                        rotaVel = calcRotaVel(float(inletMag[0]),radius)
                elif fullCaseSetupDict[geom.split('.')[0]]['ROT_WH'][0].lower() == 'false':
                    rotaVel = 0
                else:
                    sys.exit('ERROR! [%s] -> ROT_WH input invalid!' % (geom.split('.')[0]))
                print('\t\t\t\tCalculated Radius: %1.4f m' % (radius))
                print('\t\t\t\tCalculated Radial Velocity: %1.4f rad/s' % (rotaVel))
                whVel = str(rotaVel)
            
                geomString = geomString.replace('WALL_MODEL',geomWallModel).replace('WH_ORIG',whOrig).replace('WH_AXIS',whAxis).replace('WH_VEL',whVel).replace('GEOM_NAME',geom.split('.')[0])
                caseSetupStringArray.append(geomString)

        else:
        #if geom is not a rotating geometry then use the default no-slip patch
            geomString = bcStrings['GEOM']
            geomString = geomString.replace('WALL_MODEL',geomWallModel).replace('WH_ORIG',whOrig).replace('WH_AXIS',whAxis).replace('WH_VEL',whVel).replace('GEOM_NAME',geom.split('.')[0])
            
            caseSetupStringArray.append(geomString)
                
    caseSetupStrings = '\n'.join(caseSetupStringArray)
    
    templateLoc = '%s/defaultDicts/system/caseProperties' % (templateLoc)
    templateDest = 'system/caseProperties'
    copyTemplateToCase(templateLoc,templateDest)
    
    search_and_replace('system/caseProperties','<INITIAL_CONDITIONS>',initStringArray)
    search_and_replace('system/caseProperties','<BOUNDARY_CONDITIONS>',caseSetupStrings)

    
def writeDecomposeParDict(templateLoc,fullCaseSetupDict):
    print('\n\tWriting out decomposeParDict...')
    templatePath = '%s/defaultDicts/system/decomposeParDict' % (templateLoc)
    copyTemplateToCase(templatePath,'system/decomposeParDict')
    decomposition = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP']['DECOMP']
    nproc = np.prod(np.array(decomposition,dtype=int))
    print('\t\tNumber of processors: %1.0f' % (nproc))
    print('\t\tDecomposition: %s %s %s' % (decomposition[0],decomposition[1],decomposition[2]))
    search_and_replace('system/decomposeParDict','<NPROC>',str(nproc))
    search_and_replace('system/decomposeParDict','<DECOMPOSITION>','(%s %s %s)' % (decomposition[0],decomposition[1],decomposition[2]))
    

def createTopoSet(templateLoc,geomDict,fullCaseSetupDict):
    #must check for check for moving grounds, and rotating geometry
    pass
        
def writeBlockMesh(templateLoc, fullCaseSetupDict):
    print('\n\tWriting out blockMeshDict...')
    simSym = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0]
    blockMeshDictTemplatePath = "%s/defaultDicts/system/blockMeshDict" % (templateLoc)
    localBlockMeshPath = 'system/blockMeshDict'
    
    #copying the template from template location to case
    copyTemplateToCase(blockMeshDictTemplatePath,localBlockMeshPath)
    #getting the wheel center location projected to the ground
    #used to determine where the "center" of the domain should lie
    #the offset is used to shift the domain in X by amount specified
    whctx = fullCaseSetupDict['BC_SETUP']['FRT_WH_CTR'][0] #wheel front axle location projected to ground
    whcty = fullCaseSetupDict['BC_SETUP']['REFCOR'][1] #symmetry plane location
    whctz = fullCaseSetupDict['BC_SETUP']['REFCOR'][2] #wheel contact patch intersection with ground
    domOffset = fullCaseSetupDict['BC_SETUP']['DOMAIN_OFFSET'][0] #offset value of the domain
    xdom = fullCaseSetupDict['BC_SETUP']['DOMAIN_SIZE'][0] #domain size in X
    ydom = fullCaseSetupDict['BC_SETUP']['DOMAIN_SIZE'][1] #domain size in Y
    zdom = fullCaseSetupDict['BC_SETUP']['DOMAIN_SIZE'][2] #domain size in Z
    
    #calculating the domain locations
    minX = -1*float(xdom)
    minY = -1*float(ydom)/2 + float(whcty) #only keeps the lhs
    minZ = float(whctz)
    maxX = 0
    maxY = float(ydom)/2 + float(whcty)
    maxZ = float(zdom) + float(whctz)
       
    #if symmetry corrects the ymax to symmetry plane
    if simSym == 'half':
        print('\t\tSetting up case as: %s' % (simSym))
        maxY = 0 + float(whcty) #only keeps the lhs
    else:
        print('\t\tSetting up case as: %s' % (simSym))
    
    #correcting for offset
    minX = minX + float(xdom)*float(domOffset) + float(whctx)
    maxX = maxX + float(xdom)*float(domOffset) + float(whctx)
    
    print('\t\tDomain Bounding Box: MIN[%1.4f %1.4f %1.4f] MAX[%1.4f %1.4f %1.4f]' % (minX,minY,minZ,maxX,maxY,maxZ))
    print('\t\tDomain X Length: %1.2fm' % (maxX-minX))
    print('\t\tDomain Y Length: %1.2fm' % (maxY-minY))
    print('\t\tDomain Z Length: %1.2fm' % (maxZ-minZ))
    blockMeshVert = calcBlockMeshVert(minX,minY,minZ,maxX,maxY,maxZ)
    blockMeshBoundaryArray = []
    patchString = 'PATCHNAME {type PATCHTYPE;faces ((PATCHFACES));}'
    domainDict = {'y-max':{'type':{'full':'patch',
                                   'half':'symmetryPlane'},
                           'faces':'3 7 6 2'},
                  'y-min':{'type':{'full':'patch',
                                   'half':'wall'},
                           'faces':'1 5 4 0'},
                  'x-min':{'type':{'full':'patch',
                                   'half':'patch'},
                           'faces':'0 4 7 3'},
                  'x-max':{'type':{'full':'patch',
                                   'half':'patch'},
                           'faces':'2 6 5 1'},
                  'z-max':{'type':{'full':'wall',
                                   'half':'wall'},
                           'faces':'4 5 6 7'},
                  'z-min':{'type':{'full':'wall',
                                   'half':'wall'},
                           'faces':'0 3 2 1'},
                 }
    print('\t\tCreating boundary patches:')
    for key in domainDict:
        patchType = domainDict[key]['type'][simSym]
        patchFace = domainDict[key]['faces']
        patchLine = patchString.replace('PATCHNAME',key)\
                               .replace('PATCHTYPE',patchType)\
                               .replace('PATCHFACES',patchFace)\
                               .replace('\t','')
        print('\t\t\t%s: %s' % (key,patchType))

        blockMeshBoundaryArray.append(patchLine)
    
    #calculating number of blocks in each direction
    backGroundSize = float(fullCaseSetupDict['BC_SETUP']['BASE_CELL_SIZE'][0])
    nXblocks = int((maxX - minX)/backGroundSize)
    nYblocks = int((maxY - minY)/backGroundSize)
    nZblocks = int((maxZ - minZ)/backGroundSize)
    blockNumbers = '(%s %s %s)' % (nXblocks,nYblocks,nZblocks)
    #replace text in template
    blockMeshBoundary = '\n'.join(blockMeshBoundaryArray)
    search_and_replace('system/blockMeshDict','<BLOCKMESH_VERT>',blockMeshVert)
    search_and_replace('system/blockMeshDict','<BLOCKMESH_BOUNDARY>',blockMeshBoundary)
    search_and_replace('system/blockMeshDict','<BLOCK_NUMBERS>',blockNumbers)
    
    
def getBashOutput(command,searchType,searchVar):

    if type(command) == list:
        command = command
    elif type(command) == str:
        command = command.split(' ')
    else:
        sys.exit("ERROR: Invalid Type in bashOutput")
        
    bashOutput = sp.check_output(command)
    bashOutput = bashOutput.decode('utf-8')
    bashOutput = bashOutput.split('\n')
    outputLine = []
    for line in bashOutput:
        
        if searchType.lower() == 'start' and line.startswith(searchVar):
            outputLine.append(line)      
        elif searchType.lower() == 'end' and line.endswith(searchVar):
            outputLine.append(line)    
        elif searchType.lower() == 'contains' and searchVar in line:
            outputLine.append(line)    
        else:
            continue   
    return outputLine 
    
def calcBlockMeshVert(minX,minY,minZ,maxX,maxY,maxZ):
    vertArray = []
    verticies = np.array([[minX, minY, minZ], [maxX, minY, minZ], [maxX,maxY,minZ] ,[minX,maxY,minZ] ,[minX, minY, maxZ], [maxX, minY, maxZ], [maxX,maxY,maxZ], [minX,maxY,maxZ]])
    for i in verticies:
        vertString = '(%s %s %s)' % (i[0],i[1],i[2])
        vertArray.append(vertString)
    
    vertString = '\n'.join(vertArray)
    
    return vertString
    
    
