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
        for line in clusterDict['ansaMesh'].keys():
            meshingScriptArray.append(clusterDict['ansaMesh'][line])
    else:
        for line in clusterDict['snappyHexMesh'].keys():
            meshingScriptArray.append(clusterDict['snappyHexMesh'][line])
    meshingScript = '\n'.join(meshingScriptArray)

    for line in clusterDict['solve'].keys():
        if 'initialize' in line:
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

    for line in clusterDict['export'].keys():
        exportScriptArray.append(clusterDict['export'][line])
    exportScript = '\n'.join(exportScriptArray)

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
           

    
