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


# def getClusterType(templateLoc, fullCaseSetupDict):
    
    # #replaceScript header







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
    for line in clusterDict[meshKey].keys():
        meshingScriptArray.append(clusterDict[meshKey][line])
    meshingScript = '\n'.join(meshingScriptArray)

    #SRF cornering is solved directly in the rotating frame with CORNER_SOLVER (default SRFSimpleFoam).
    #The standard initialisers (potentialFoam / simpleFoam) run in the absolute frame and would seed an
    #inconsistent field for an SRF solve, so initialisation is skipped entirely when cornering.
    runCornering = ('CORNERING_SETUP' in fullCaseSetupDict and
                    fullCaseSetupDict['CORNERING_SETUP']['RUN_CORNERING'][0].lower() == 'true')

    for line in clusterDict['solve'].keys():
        if 'initialize' in line:
            if runCornering:
                continue
            if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0] == 'potential':
                solveScriptArray.append(clusterDict['solve']['initialize']['initializePotential'])
            elif fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0] == 'steady':
                solveScriptArray.append(clusterDict['solve']['initialize']['initializeSteady'])
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
            #transient cornering runs SRFPimpleFoam. CORNER_SOLVER defaults to the steady solver name
            #(SRFSimpleFoam), so auto-upgrade it to the transient SRF solver when left at that default.
            if solverApp == 'SRFSimpleFoam':
                solverApp = 'SRFPimpleFoam'
            #bare (un-anchored) swaps so both the `[ -f system/controlDictPiso ]` guards and the
            #`cp system/controlDictPiso ...` lines flip to the SRF variant. Safe because the meshing and
            #solve scripts never reference controlDictPisoExport / fvSolutionPisoExport.
            fileSwaps = [('system/controlDictPiso', 'system/controlDictSRFPiso'),
                         ('system/fvSolutionPiso',  'system/fvSolutionSRFPiso'),
                         ('system/fvSchemesPiso',   'system/fvSchemesSRFPiso')]
            execOld = 'foamExec pisoFoam -parallel >> log.pisoFoam'
            execNew = 'foamExec %s -parallel >> log.pisoFoam' % (solverApp)
        else:
            #cornering swaps the steady "Simple" dictionaries for their SRF counterparts so that meshing
            #(getControlDict / getFvDict) and the solve step copy the SRF controlDict/fvSolution/fvSchemes,
            #and createZeroDirectory reads application = SRFSimpleFoam and emits the Urel field. The trailing
            #space anchors the filename so controlDictSimpleExport / controlDictSimpleInit are not matched.
            fileSwaps = [('system/controlDictSimple ', 'system/controlDictSRFSimple '),
                         ('system/fvSolutionSimple ',  'system/fvSolutionSRFSimple '),
                         ('system/fvSchemesSimple ',   'system/fvSchemesSRFSimple ')]
            execOld = 'foamExec simpleFoam -parallel >> log.simpleFoam'
            execNew = 'foamExec %s -parallel >> log.simpleFoam' % (solverApp)
        #Meshing's getControlDict selects the controlDict by existence and falls back to
        #controlDictPotential. For cornering only the SRF controlDict is ever written, and relying on the
        #substring swaps below to flip that chain is fragile: if it misses, meshing leaves controlDict as
        #controlDictPotential and the solve's createZeroDirectory then reads application = potentialFoam
        #instead of the SRF solver. Replace the whole getControlDict chain with an explicit copy of the SRF
        #controlDict so it can never fall back to controlDictPotential. Done before the swap loop so the
        #exact original string still matches; the SRF replacement contains no swap-target substrings.
        srfControlDict = 'controlDictSRFPiso' if simTypeLower == 'transient' else 'controlDictSRFSimple'
        corneringGetControlDict = ('if [ -f system/%s ]; then cp system/%s system/controlDict; '
                                   'else exit 1; fi') % (srfControlDict, srfControlDict)
        meshingScript = meshingScript.replace(clusterDict[meshKey]['getControlDict'], corneringGetControlDict)
        for src, dst in fileSwaps:
            meshingScript = meshingScript.replace(src, dst)
            solveScript = solveScript.replace(src, dst)
        solveScript = solveScript.replace(execOld, execNew)
        #The solve preamble copies controlDictPotential->controlDict (absolute-frame potential init) right
        #before createZeroDirectory reads the application. Cornering never writes controlDictPotential, but a
        #stale file from an earlier run would clobber the SRF controlDict and make createZeroDirectory emit U
        #instead of Urel, so drop that copy entirely for cornering cases.
        potentialPreamble = (" if [[ -f 'system/controlDictPotential' ]]; then "
                             "cp system/controlDictPotential system/controlDict; "
                             "cp system/fvSolutionPotential system/fvSolution; fi;")
        solveScript = solveScript.replace(potentialPreamble, '')

    for line in clusterDict['export'].keys():
        exportScriptArray.append(clusterDict['export'][line])
    exportScript = '\n'.join(exportScriptArray)

    if runCornering:
        #The export step template runs `pisoFoam -postProcess -latestTime -fields "(pMean UMean)"`. Two
        #coordinated changes are required for cornering:
        # 1. pisoFoam runs in the absolute frame and does not know about the SRF rotation, so it cannot
        #    relate the rotating-frame fields. The export must run the same SRF solver used for the solve
        #    (SRFSimpleFoam / SRFPimpleFoam, held in solverApp) so -postProcess understands the frame.
        # 2. The -fields list is the ONLY thing read into the objectRegistry, and the cornering function
        #    objects read both the absolute means (UMean) and the rotating-frame fields (UrelMean, Urel).
        #    Unless every consumed field is preloaded the FOs fail with e.g. "failed lookup of UrelMean
        #    (objectRegistry region0)", even though U/UMean and Urel/UrelMean all exist on disk. Preload all
        #    four (pMean UMean UrelMean Urel) so every lookup resolves regardless of frame.
        exportScript = exportScript.replace('foamExec pisoFoam -postProcess',
                                            'foamExec %s -postProcess' % (solverApp))
        exportScript = exportScript.replace('(pMean UMean)', '(pMean UMean UrelMean Urel)')

    for line in clusterDict['post'].keys():
        postScriptArray.append(clusterDict['post'][line])
    postScript = '\n'.join(postScriptArray)

    for line in clusterDict['postMesh'].keys():
        postMeshScriptArray.append(clusterDict['postMesh'][line])
    postMeshScript = '\n'.join(postMeshScriptArray)




    with open('meshingScript', 'w') as m:
        m.write(meshingScript)

    with open('solveScript', 'w') as m:
        m.write(solveScript)

    with open('exportScript', 'w') as m:
        m.write(exportScript)

    with open('postScript', 'w') as m:
        m.write(postScript)

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
           

    
