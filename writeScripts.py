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
    
    for line in clusterDict['meshing'].keys():
        meshingScriptArray.append(clusterDict['meshing'][line])
    meshingScript = '\n'.join(meshingScriptArray)

    for line in clusterDict['solve'].keys():
        if 'initialize' in line:
            if fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0] == 'potential':
                solveScriptArray.append(clusterDict['solve']['initialize']['initializePotential'])
            elif fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_INIT'][0] == 'steady':
                solveScriptArray.append(clusterDict['solve']['initialize']['initializeSteady'])

        else:
            solveScriptArray.append(clusterDict['solve'][line])
    solveScript = '\n'.join(solveScriptArray)

    for line in clusterDict['export'].keys():
        exportScriptArray.append(clusterDict['export'][line])
    exportScript = '\n'.join(exportScriptArray)

    for line in clusterDict['post'].keys():
        postScriptArray.append(clusterDict['post'][line])
    postScript = '\n'.join(postScriptArray)




    with open('meshingScript', 'w') as m:
        m.write(meshingScript)

    with open('solveScript', 'w') as m:
        m.write(solveScript)

    with open('exportScript', 'w') as m:
        m.write(exportScript)

    with open('postScript', 'w') as m:
        m.write(postScript)
    
    
    
def copyScripts(templateLoc, fullCaseSetupDict,case):

    scriptDict = {'meshing':'MS',
                  'solve':'SO',
                  'export':'EX',
                  'postPro':'PP',
                  'postScript':'PO'}

    print('\tWriting copyingScripts...')
    scriptVersion = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP']['SCRIPT_VERSION'][0]
    copyScriptsPath = '%s/defaultCluster/slurm/scripts/%s/*Script' % (templateLoc,scriptVersion)
    
    makeScripts(templateLoc,fullCaseSetupDict)
    scriptsList = glob.glob(copyScriptsPath)
    for script in scriptsList:
        scriptName = script.split('/')[-1]
        #copyTemplateToCase(script,scriptName)
        if os.path.isfile(scriptName):
            for line in fileinput.  input(scriptName, inplace=True):
                if line.strip().startswith('#SBATCH --job-name='):
                    for key in scriptDict.keys():
                        if key in scriptName:
                            keyCode = scriptDict[key]
                            jobCode = fullCaseSetupDict['TITLES']['JOBCODE'][0]
                            sym = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0]
                            line = '#SBATCH --job-name=%s_%s%s_%s\n\n' % (jobCode,keyCode,case,sym)
                            
                sys.stdout.write(line)
    
           
        
    #print(scriptsList)
    
    
    
    
    