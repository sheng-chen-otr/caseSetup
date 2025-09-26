import os
import sys
import numpy as np
import re
import pandas as pd
import configparser
import argparse
import matplotlib.pyplot as plt
import scipy.stats as st
import glob
from plotForces import *
from summary import *
#from estimateStatisticalError import *
from forceConvergencePlot import *

# Set default matplotlib parameters
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)

def main():

    # Get case information
    global casePath, caseName, path, case,caseLoc,job
    casePath = os.getcwd() 
    caseName = casePath.split('/')[-1]
    path = os.path.split(casePath)[0]
    case = caseName
    job = os.path.basename(os.path.dirname(path))
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        prog='OpenFOAM Post Processing Tool',
        description='Process and analyze OpenFOAM simulation results'
    )
    
    # Add main operation mode arguments
    parser.add_argument('--summary', action='store_true', 
                       help='Generate case summary')
    parser.add_argument('--forces', action='store_true', 
                       help='Plot force coefficients')
    
    # Add force plotting specific arguments
    parser.add_argument('-p', '--plotData', default=['Cd','Cl','CoP'],
                       nargs='+', choices=['Cd','Cl','CoP','Cl(f)','Cl(r)','Cs(f)','Cs(r)'],
                       help='Data to plot when using --forces')
    parser.add_argument('-t', '--trial', default=[caseName], nargs='+',
                       help='Trials to plot (default: current trial)')
    parser.add_argument('-s', '--savePlots', action='store_true',
                       help='Save generated plots')
    parser.add_argument('--skipStats', action='store_true',
                       help='Skip calculating statistics')
    parser.add_argument('--yscaling', default='default',
                       help='Y-axis scaling for plots')
    parser.add_argument('--saveFormat', default='png',
                       choices=['png', 'eps', 'jpeg'],
                       help='Format for saved plots')
    parser.add_argument('--avgTime', type=float,
                       help='Manual averaging start time')
    parser.add_argument('--plotEveryOther', default=10, type=int,
                       help='Number of every other time steps to plot, reduces messiness in plot')

    args = parser.parse_args()
    
    if args.summary:
        generate_summary()

    if args.forces:
        casePathDict, caseLoc = setCasePaths(args.trial,casePath)
        casePathDict = getCaseData(casePathDict)
        casePathDict = makePandasArrays(args,casePathDict)
        plotData(args,caseLoc,casePathDict)
        
def generate_summary():
        
    caseSetupPath="%s/caseSetup" % (casePath)
    fullCaseSetupDict = configparser.ConfigParser()
    fullCaseSetupDict.optionxform = str
    fullCaseSetupDict.read_file(open(caseSetupPath))
    configSections = fullCaseSetupDict.sections()

    coeffFiles = getCoeffPaths(casePath, case)
    for part in coeffFiles:
        if part != 'all':
            avgData = averageCoeffs(fullCaseSetupDict,case,part,coeffFiles)
    avgData = averageCoeffs(fullCaseSetupDict,case,'all',coeffFiles)
    
    numCells,mesher,sym = cellCount(fullCaseSetupDict,casePath,case)
    inletMag,lastTime,yaw,movingGround,rotatingWheels,simType,turbModel = bcParser(fullCaseSetupDict,path,case)
    runDate,runTime,version,solver = getOfVersion(casePath)
    refArea = float(fullCaseSetupDict['BC_SETUP']['REFAREA'][0])
    #default datas
    rowNames = ['Job','Trial','Solver','Version','Run Date','Solve Time','Num. Cells','Mesher','Symmetry','Ref. Area (m^2)','Iterations','Simulation Type','Moving Ground','Rotating Wheels','Turbulence Model','Velocity','Yaw','Cd','Cl','Cl/Cd','%Front','Cd CI','Cl CI']
    data = [job,case,solver,version,runDate,runTime,numCells,mesher,sym.lower(),refArea,avgData['endTime'],simType.lower(),movingGround,rotatingWheels,turbModel,inletMag,yaw,avgData['cd'],avgData['cl'],avgData['cl/cd'],avgData['cop'],avgData['cd_ci'],avgData['cl_ci']]
    try:
        porousData = getPorousData()
        #adding porous data
        for key in porousData.keys():
            rowNames.append(str(key))
            data.append(str(porousData[key]))
    except:
        print('\tUnable to get porous media data, skipping...')
    summary = pd.DataFrame(columns=rowNames)
    summary.loc[-1] = data
    print("\n\n")
    for col in summary.columns:
        print('{:>100s}{:>30s}'.format(col,str(summary[col].values[0])))
        
    summary = summary.transpose()
    summary.to_csv("%s/%s/summary.csv"% (path,case),header=False)

main()
