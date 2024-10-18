import os
import sys
import numpy as np
import subprocess as sp
import math
from utilities import *
from writeSystem import *
from writeConstant import *


clusterDict = {'slurm':{'header':   '''
                                    #SBATCH --nodes=NUM_NODES
                                    #SBATCH --ntasks-per-node=TASK_PER_NODE
                                    #SBATCH --partition=PARTITION_NAME
                                    #SBATCH --output=01_%x_meshing_%j.out
                                    #SBATCH --exclusive
                                    #SBATCH --job-name=JOB_CODE_CASENAME_SYMMETRY
                                    ''' ,
                        'meshing': {'NUM_NODES': 'NUM_NODES', 
                                    'TASK_PER_NODE':'TASK_PER_NODE', 
                                    'JOB_CODE':'JOB_CODE', 
                                    'CASE_NAME': 'CASENAME',
                                    '_SYMMETRY' : 'SIM_SYM'
                                    }
                        'solve':    {'NUM_NODES': 'NUM_NODES', 
                                    'TASK_PER_NODE':'TASK_PER_NODE', 
                                    'JOB_CODE':'JOB_CODE', 
                                    'CASE_NAME': 'CASENAME',
                                    '_SYMMETRY' : 'SIM_SYM'
                                    }
                        'export':   {'NUM_NODES': 'NUM_NODES', 
                                    'TASK_PER_NODE':'TASK_PER_NODE', 
                                    'JOB_CODE':'JOB_CODE', 
                                    'CASE_NAME': 'CASENAME',
                                    '_SYMMETRY' : 'SIM_SYM'
                                    }
                        'postPro':  {'NUM_NODES': '1', 
                                    'TASK_PER_NODE':'8', 
                                    'JOB_CODE':'JOB_CODE', 
                                    'CASE_NAME': 'CASENAME',
                                    '_SYMMETRY' : 'SIM_SYM'
                                    }
                        }   
                }


def getClusterType(templateLoc, fullCaseSetupDict):
    
    #replaceScript header
    
    
    
    
    