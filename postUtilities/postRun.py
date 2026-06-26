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
from collections import OrderedDict
#from plotForces import *
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


def isCaseComplete(casePath):
    #A ride-height child is "complete" when its OpenFOAM solver finished, which the solver signals by
    #printing a standalone "End" line at the very end of its log. Checking that marker (rather than
    #comparing the latest time against controlDict endTime) is robust to convergence-based early stops
    #(residualControl), where a steady case legitimately ends before endTime. SRF cornering solves still
    #write log.simpleFoam (steady) / log.pisoFoam (transient), so the same log names cover cornering.
    for logName in ('log.simpleFoam', 'log.pisoFoam'):
        logPath = os.path.join(casePath, logName)
        if not os.path.isfile(logPath):
            continue
        try:
            with open(logPath, 'r') as solveLog:
                for line in solveLog:
                    if line.strip() == 'End':
                        return True
        except Exception:
            return False
        #log exists but never reached the End marker -> still running or crashed
        return False
    return False


def getCorneringInfo(fullCaseSetupDict, casePath, case):
    #Returns the cornering / per-corner descriptors for a single case as an ordered dict of
    #column-name -> value, ready to extend the summary rowNames/data.
    #  - Cornering flag, corner radius and direction come from THIS case's own caseSetup
    #    (createRideHeightCases injects the per-point CORNER_RADIUS / CORNER_DIR into each child).
    #  - The per-corner ride-height change and steer angle come from the parent's
    #    rideHeights_updated.csv (written by createRideHeightCases, one row per point keyed by
    #    'caseName'), which lives one directory above the child case.
    #Standalone / non-cornering cases fall back to N/A so the column set stays consistent.
    info = OrderedDict([
        ('Cornering', 'False'),
        ('Corner Radius (m)', 'N/A'),
        ('Corner Direction', 'N/A'),
        ('Steer Angle (deg)', 'N/A'),
        ('RH_FL', 'N/A'),
        ('RH_FR', 'N/A'),
        ('RH_RL', 'N/A'),
        ('RH_RR', 'N/A'),
    ])

    if fullCaseSetupDict.has_section('CORNERING_SETUP'):
        runCornering = fullCaseSetupDict['CORNERING_SETUP'].get('RUN_CORNERING', 'False').strip().lower() == 'true'
        info['Cornering'] = str(runCornering)
        if runCornering:
            info['Corner Radius (m)'] = fullCaseSetupDict['CORNERING_SETUP'].get('CORNER_RADIUS', '').strip() or 'N/A'
            info['Corner Direction'] = fullCaseSetupDict['CORNERING_SETUP'].get('CORNER_DIR', '').strip() or 'N/A'

    #per-corner ride-height change + steer from the parent rideHeights_updated.csv
    rhCsvPath = os.path.join(os.path.dirname(casePath), 'rideHeights_updated.csv')
    if os.path.isfile(rhCsvPath):
        try:
            rhMap = pd.read_csv(rhCsvPath)
            if 'caseName' in rhMap.columns:
                rowMatch = rhMap[rhMap['caseName'].astype(str) == case]
                if len(rowMatch) > 0:
                    row = rowMatch.iloc[0]
                    cornerCols = OrderedDict([
                        ('Ride Height Change FL', 'wheel_fl'),
                        ('Ride Height Change FR', 'wheel_fr'),
                        ('Ride Height Change RL', 'wheel_rl'),
                        ('Ride Height Change RR', 'wheel_rr'),
                    ])
                    for label, col in cornerCols.items():
                        if col in rhMap.columns:
                            info[label] = round(float(row[col]), 3)
                    for steerCand in ('steer', 'steer_deg', 'steer_angle'):
                        if steerCand in rhMap.columns:
                            info['Steer Angle (deg)'] = round(float(row[steerCand]), 3)
                            break
        except Exception as e:
            print('\tUnable to read cornering/ride-height info from %s: %s' % (rhCsvPath, e))

    return info


def generate_summary():
        
    caseSetupPath="%s/caseSetup" % (casePath)
    fullCaseSetupDict = configparser.ConfigParser()
    fullCaseSetupDict.optionxform = str
    fullCaseSetupDict.read_file(open(caseSetupPath))
    configSections = fullCaseSetupDict.sections()
    rideHeightSetup = fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RUN_RIDE_HEIGHT']
   

    if rideHeightSetup.lower() == 'true' and not '_' in os.path.basename(os.getcwd()):
        print('Detected ride height map, averaging map values!')
        runPoints = fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RUN_RH_POINTS']
        caseName = os.path.basename(os.getcwd())
        rhCases = []
        #rhAvgData = pd.DataFrame(columns=['CD','CL','CLF','CLR','CSF','CSR','CI-CD','CI-CL'])
        rhAvgData = pd.DataFrame()
        for point in runPoints:
            if point == ' ':
                continue
            rideHeightDir = "%s_%s" % (caseName,point)
            if not os.path.isdir(rideHeightDir):
                print('\tWARNING! Ride height point %s not found (%s), skipping from average.' % (point, rideHeightDir))
                continue
            if not isCaseComplete(os.path.join(os.getcwd(), rideHeightDir)):
                print('\tWARNING! Ride height point %s (%s) is not complete (no "End" in solve log), skipping from average.' % (point, rideHeightDir))
                continue
            rhCases.append(rideHeightDir)
        if len(rhCases) < 1:
            print('\tNo complete ride height points to average; skipping parent summary.')
            return
        porousDict = {}
        partsDict = {}
        for case in rhCases:
            rhPath = os.path.join(os.getcwd(),case)
            try:
                coeffFiles = getCoeffPaths(rhPath)
                for part in coeffFiles:
                   
                    if part != 'all':
                        partAverage,partAverageArray = averageCoeffs(fullCaseSetupDict,case,part,coeffFiles)
                        partAverage = pd.DataFrame([partAverage])
                        partsDict[case] = {} #prepare to collect the forces for parts
                        partsDict[case][part] = partAverage
                avgData,averagedArray = averageCoeffs(fullCaseSetupDict,case,'all',coeffFiles)
                numCells,mesher,sym = cellCount(fullCaseSetupDict,os.getcwd(),case)
                inletMag,lastTime,yaw,movingGround,rotatingWheels,simType,turbModel = bcParser(fullCaseSetupDict,os.getcwd(),case)
                runDate,runTime,version,solver = getOfVersion(rhPath)
                refArea = float(fullCaseSetupDict['BC_SETUP']['REFAREA'][0])
                caseAvgData = pd.DataFrame([avgData])
                #rhAvgData = pd.concat([rhAvgData,averagedArray],axis=0)
                rhAvgData = pd.concat([rhAvgData,caseAvgData],axis=0)
            except Exception as E:
                print('\t\tUnable to average: %s' % (case))
                #print(E)
            
            try:
                porousData = getPorousData(path,case)
            except Exception as e:
                print('\tUnable to get porous media data, skipping...')
                print(e)
            if len(porousData.keys()) < 1:
                porousDict = {}
            else:
                porousDict[case] = porousData
        #average all the part data
        avgPartDict = {}
        for part in partsDict[list(partsDict.keys())[0]].keys():
            avgPartDict[part] = pd.concat([partsDict[case][part] for case in partsDict.keys()]).mean(axis=0).round(3)
           
        
        

        
        rhMeans = rhAvgData.mean(axis=0).round(3)
        #average all the porous media data
        if len(porousDict.keys()) < 1:
            print('\t\tNo porous data to average!')
            avgPorous = None
        else:
            #get the names of all the porous media data
            porousKeys = porousDict[porousDict.keys()[0]].keys()
            porousArray = pd.DataFrame(columns=porousKeys)
            for case in porousDict.keys():
                casePorousArray = pd.DataFrame.from_dict(porousDict[case])
                porousArray = pd.concat([porousArray,casePorousArray])
            avgPorous = porousArray.mean(axis=0).round(3)

        #default datas
        rowNames = ['Job','Trial','Solver','Version','Run Date','Solve Time','Num. Cells','Mesher','Symmetry','Ref. Area (m^2)','Iterations','Simulation Type','Moving Ground','Rotating Wheels','Turbulence Model','Velocity','Yaw','Cd','Cl','Cl/Cd','%Front','Cd CI','Cl CI']
        data = [job,caseName,solver,version,'N/A','N/A','N/A',mesher,sym.lower(),refArea,'N/A',simType.lower(),movingGround,rotatingWheels,turbModel,inletMag,yaw,rhMeans['cd'],rhMeans['cl'],rhMeans['cl/cd'],rhMeans['cop'],rhMeans['cd_ci'],rhMeans['cl_ci']]
        if len(avgPartDict.keys()) > 0:
            for part in avgPartDict.keys():
                partVarDict = {'CL':'cl',
                               'CD':'cd'}
                
                for varkey in partVarDict.keys():
                    rowNames.append(part + ' ' + varkey)
                    data.append(avgPartDict[part][partVarDict[varkey]])
        if avgPorous != None:
            for col in avgPorous.index:
                rowNames.append(col)
            for val in avgPorous:
                data.append(val)

        summary = pd.DataFrame(columns=rowNames)
        summary.loc[-1] = data
        print("\n\n")
        for col in summary.columns:
            print('{:>100s}{:>30s}'.format(col,str(summary[col].values[0])))
            
        summary = summary.transpose()
        summary.to_csv("%s/summary.csv"% (casePath),header=False)
                
    else:
        case = os.path.basename(casePath)
        coeffFiles = getCoeffPaths(casePath)
        partsDict = {}
        for part in coeffFiles:
            if part != 'all':
                partsDict[part],avgDataArray = averageCoeffs(fullCaseSetupDict,case,part,coeffFiles)
        avgData,allDataArray = averageCoeffs(fullCaseSetupDict,case,'all',coeffFiles)
        
        numCells,mesher,sym = cellCount(fullCaseSetupDict,casePath,case)
        inletMag,lastTime,yaw,movingGround,rotatingWheels,simType,turbModel = bcParser(fullCaseSetupDict,path,case)
        runDate,runTime,version,solver = getOfVersion(casePath)
        refArea = float(fullCaseSetupDict['BC_SETUP']['REFAREA'][0])
    
        #default datas
        rowNames = ['Job','Trial','Solver','Version','Run Date','Solve Time','Num. Cells','Mesher','Symmetry','Ref. Area (m^2)','Iterations','Simulation Type','Moving Ground','Rotating Wheels','Turbulence Model','Velocity','Yaw','Cd','Cl','Cl/Cd','%Front','Cd CI','Cl CI']
        data = [job,case,solver,version,runDate,runTime,numCells,mesher,sym.lower(),refArea,avgData['endTime'],simType.lower(),movingGround,rotatingWheels,turbModel,inletMag,yaw,avgData['cd'],avgData['cl'],avgData['cl/cd'],avgData['cop'],avgData['cd_ci'],avgData['cl_ci']]
        #cornering + per-corner ride-height / steer descriptors for this case
        corneringInfo = getCorneringInfo(fullCaseSetupDict, casePath, case)
        for label, value in corneringInfo.items():
            rowNames.append(label)
            data.append(value)
        for part in partsDict.keys():
            partVarDict = {'CL':'cl',
                           'CD':'cd'}
            
            for varkey in partVarDict.keys():
                rowNames.append(part + ' ' + varkey)
                data.append(partsDict[part][partVarDict[varkey]])
        try:
            porousData = getPorousData(path,case)
            #adding porous data
            for key in porousData.keys():
                rowNames.append(str(key))
                data.append(str(porousData[key]))
        except Exception as e:
            print('\tUnable to get porous media data, skipping...')
            print(e)
        summary = pd.DataFrame(columns=rowNames)
        summary.loc[-1] = data
        print("\n\n")
        for col in summary.columns:
            print('{:>100s}{:>30s}'.format(col,str(summary[col].values[0])))
            
        summary = summary.transpose()
        summary.to_csv("%s/%s/summary.csv"% (path,case),header=False)

main()
