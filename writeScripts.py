import os
import sys
import numpy as np
import subprocess as sp
import math
import glob as glob
from utilities import *
import fileinput
import json
from writeSystem import *
from writeConstant import *


# clusterDict = {'slurm':{'header':   '''
                                    # #SBATCH --nodes=NUM_NODES
                                    # #SBATCH --ntasks-per-node=TASK_PER_NODE
                                    # #SBATCH --partition=PARTITION_NAME
                                    # #SBATCH --output=SCRIPT_NUMBER_%x_SCRIPT_TYPE_%j.out
                                    # #SBATCH --exclusive
                                    # #SBATCH --job-name=JOB_CODE_CASENAME_SYMMETRY
                                    # ''' ,
                        # 'meshing': {'NUM_NODES': 'NUM_NODES', 
                                    # 'TASK_PER_NODE':'TASK_PER_NODE', 
                                    # 'JOB_CODE':'JOB_CODE', 
                                    # 'CASE_NAME': 'CASENAME',
                                    # '_SYMMETRY' : 'SIM_SYM'
                                    # },
                        # 'solve':    {'NUM_NODES': 'NUM_NODES', 
                                    # 'TASK_PER_NODE':'TASK_PER_NODE', 
                                    # 'JOB_CODE':'JOB_CODE', 
                                    # 'CASE_NAME': 'CASENAME',
                                    # '_SYMMETRY' : 'SIM_SYM'
                                    # },
                        # 'export':   {'NUM_NODES': 'NUM_NODES', 
                                    # 'TASK_PER_NODE':'TASK_PER_NODE', 
                                    # 'JOB_CODE':'JOB_CODE', 
                                    # 'CASE_NAME': 'CASENAME',
                                    # '_SYMMETRY' : 'SIM_SYM'
                                    # },
                        # 'postPro':  {'NUM_NODES': 1, 
                                    # 'TASK_PER_NODE':8, 
                                    # 'JOB_CODE':'JOB_CODE', 
                                    # 'CASE_NAME': 'CASENAME',
                                    # '_SYMMETRY' : 'SIM_SYM'
                                    # }
                        # }   
                # }



def makeScripts(templateLoc,fullCaseSetupDict):

    clusterDictPath = '%s/defaultCluster/slurm/clusterDict' % (templateLoc)

    #getting default variables
    with open(clusterDictPath) as f: 
        clusterDict = f.read()

    clusterDict = json.loads(clusterDict)
    print('\t\tImporting cluster script values')
    meshingScriptArray = []
    solveScriptArray = []
    exportScriptArray = []
    postScriptArray = []
    postMeshScriptArray = []
    
    if 'ansa' in fullCaseSetupDict['GLOBAL_REFINEMENT']['TEMPLATE_TYPE'][0].lower():
        meshKey = 'ansaMesh'
    else:
        meshKey = 'snappyHexMesh'
    if meshKey in clusterDict.keys():
        for line in clusterDict[meshKey].keys():
            meshingScriptArray.append(clusterDict[meshKey][line])
        meshingScript = '\n'.join(meshingScriptArray)
        with open('meshingScript', 'w') as m:
            m.write(meshingScript)

    #srf cornering solves in the rotating frame with CORNER_SOLVER, so skip the
    #absolute-frame initialisers (they'd seed a bad field)
    runCornering = ('CORNERING_SETUP' in fullCaseSetupDict and
                    fullCaseSetupDict['CORNERING_SETUP']['RUN_CORNERING'][0].lower() == 'true')

    if 'solve' in clusterDict.keys():
        for line in clusterDict['solve'].keys():
            if 'initialize' in line:
                if runCornering:
                    continue
                if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0] == 'potential':
                    solveScriptArray.append(clusterDict['solve']['initialize']['initializePotential'])
                elif fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0] == 'steady':
                    solveScriptArray.append(clusterDict['solve']['initialize']['initializeSteady'])
            elif 'topoSet' in line:
                #ansaMesh can't snap the porous/MRF interface, so build those cellZones with topoSet
                #(system/topoSetDict, written by createTopoSet) before the solver runs.
                if hasTopoSetRegions(fullCaseSetupDict):
                    solveScriptArray.append(clusterDict['solve']['topoSet'])
            elif 'groundPatch' in line:
                #split the z-min ground into user belt/plate patches (topoSetDictGround +
                #createPatchDict) after decomposePar, before createZeroDirectory needs them
                if hasGroundZones(fullCaseSetupDict):
                    solveScriptArray.append(clusterDict['solve']['groundPatch'])
            elif 'solve' in line:
                if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_TYPE'][0].lower() == 'steady':
                    solveScriptArray.append(clusterDict['solve']['solve']['steady'])
                elif fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_TYPE'][0].lower() == 'transient':
                    solveScriptArray.append(clusterDict['solve']['solve']['transient'])

            else:
                solveScriptArray.append(clusterDict['solve'][line])
        solveScript = '\n'.join(solveScriptArray)
        if runCornering:
            solverApp = fullCaseSetupDict['CORNERING_SETUP']['CORNER_SOLVER'][0]
            simTypeLower = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_TYPE'][0].lower()
            if simTypeLower == 'transient':
                #transient cornering uses SRFPimpleFoam, auto-upgrade from the steady default
                if solverApp == 'SRFSimpleFoam':
                    solverApp = 'SRFPimpleFoam'
                #bare swaps so the controlDictPiso guards and cp lines both flip to SRF
                fileSwaps = [('system/controlDictPiso', 'system/controlDictSRFPiso'),
                            ('system/fvSolutionPiso',  'system/fvSolutionSRFPiso'),
                            ('system/fvSchemesPiso',   'system/fvSchemesSRFPiso')]
                execOld = 'pisoFoam -parallel >> log.pisoFoam'
                execNew = '%s -parallel >> log.pisoFoam' % (solverApp)
            else:
                #steady cornering swaps the Simple dicts for SRF. trailing space anchors the
                #filename so the Export/Init variants aren't matched
                fileSwaps = [('system/controlDictSimple ', 'system/controlDictSRFSimple '),
                            ('system/fvSolutionSimple ',  'system/fvSolutionSRFSimple '),
                            ('system/fvSchemesSimple ',   'system/fvSchemesSRFSimple ')]
                execOld = 'simpleFoam -parallel >> log.simpleFoam'
                execNew = '%s -parallel >> log.simpleFoam' % (solverApp)
            #replace the whole getControlDict chain with an explicit copy of the SRF controlDict
            #so meshing can't fall back to controlDictPotential. done before the swap loop so the
            #original string still matches
            srfControlDict = 'controlDictSRFPiso' if simTypeLower == 'transient' else 'controlDictSRFSimple'
            corneringGetControlDict = ('if [ -f system/%s ]; then cp system/%s system/controlDict; '
                                    'else exit 1; fi') % (srfControlDict, srfControlDict)
            meshingScript = meshingScript.replace(clusterDict[meshKey]['getControlDict'], corneringGetControlDict)
            for src, dst in fileSwaps:
                meshingScript = meshingScript.replace(src, dst)
                solveScript = solveScript.replace(src, dst)
            solveScript = solveScript.replace(execOld, execNew)
            #drop the potential-init preamble for cornering - a stale controlDictPotential would
            #clobber the SRF controlDict and make createZeroDirectory emit U instead of Urel
            potentialPreamble = (" if [[ -f 'system/controlDictPotential' ]]; then "
                                "cp system/controlDictPotential system/controlDict; "
                                "cp system/fvSolutionPotential system/fvSolution; fi;")
            solveScript = solveScript.replace(potentialPreamble, '')
        with open('solveScript', 'w') as m:
            m.write(solveScript)


    if 'export' in clusterDict.keys():
        for line in clusterDict['export'].keys():
            exportScriptArray.append(clusterDict['export'][line])
        exportScript = '\n'.join(exportScriptArray)

        if runCornering:
            #export must run the SRF solver (not pisoFoam) so -postProcess knows the frame, and
            #-fields must preload everything the FOs read (UrelMean/Urel aren't loaded from disk)
            exportScript = exportScript.replace('pisoFoam -postProcess',
                                                '%s -postProcess' % (solverApp))
            exportScript = exportScript.replace('(pMean UMean)', '(pMean UMean UrelMean Urel)')
        with open('exportScript', 'w') as m:
            m.write(exportScript)

    if 'post' in clusterDict.keys():
        for line in clusterDict['post'].keys():
            postScriptArray.append(clusterDict['post'][line])
        postScript = '\n'.join(postScriptArray)
        with open('postScript', 'w') as m:
            m.write(postScript)

    if 'postMesh' in clusterDict.keys():
        for line in clusterDict['postMesh'].keys():
            postMeshScriptArray.append(clusterDict['postMesh'][line])
        postMeshScript = '\n'.join(postMeshScriptArray)

        with open('postMeshScript', 'w') as m:
            m.write(postMeshScript)




   

    

    

    

    
    
    
    
def copyScripts(templateLoc, fullCaseSetupDict,case):

    scriptDict = {'meshingScript':'MS',
                  'solveScript':'SO',
                  'exportScript':'EX',
                  'postScript':'PO',
                  'postMeshScript':'PM'}

    print('\tWriting copyingScripts...')
    
    makeScripts(templateLoc,fullCaseSetupDict)
    for scriptName in scriptDict.keys():
        keyCode = scriptDict[scriptName]
        jobCode = fullCaseSetupDict['TITLES']['JOBCODE'][0]
        sym = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0]
        lineToFind = '#SBATCH --job-name='
        newLine = '#SBATCH --job-name=%s_%s%s_%s\n\n' % (jobCode,keyCode,case,sym)
        if os.path.isfile(scriptName):
            replace_line_in_file(scriptName,lineToFind,newLine)

def replace_line_in_file(file_path,line_to_find , replacement_line):
    # Read the file content into a list of lines
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Modify the desired line(s)
    for i, line in enumerate(lines):
        if line_to_find in line: # Check if the target string is in the line
            lines[i] = replacement_line + '\n' # Replace the line (add newline character)

    # Overwrite the file with the modified content
    with open(file_path, 'w') as file:
        file.writelines(lines)
           

    
