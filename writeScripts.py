import os
import sys
import numpy as np
import subprocess as sp
import math
import glob as glob
from utilities import *
import fileinput
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
    
    
def copyScripts(templateLoc, fullCaseSetupDict,case):

    scriptDict = {'meshing':'MS',
                  'solve':'SO',
                  'export':'EX',
                  'postPro':'PP',
                  'postScript':'PO'}

    print('\tWriting copyingScripts...')
    scriptVersion = fullCaseSetupDict['GLOBAL_COMPUTE_SETUP']['SCRIPT_VERSION'][0]
    copyScriptsPath = '%s/defaultCluster/slurm/scripts/%s/*Script' % (templateLoc,scriptVersion)
    
    
    scriptsList = glob.glob(copyScriptsPath)
    for script in scriptsList:
        scriptName = script.split('/')[-1]
        copyTemplateToCase(script,scriptName)
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
    
    
    
    
    