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
    parser.add_argument('--wingPlots', action='store_true',
                       help='Plot variable versus x from wing intersection CSV files')
    parser.add_argument('--wingCases', nargs='*', default=[],
                       help='Additional case names or paths to overlay on wing plots')
    parser.add_argument('--wingVariables', nargs='*', default=[],   
                       help='Wing variable columns to plot (default: all available)')
    
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

    if args.wingPlots:
        casePathDict, caseLoc = setCasePaths(args.trial,casePath)
        plotWingPressure(args, casePathDict, caseLoc)

    if not args.summary and not args.forces and not args.wingPlots:
        parser.print_help()


def plotWingPressure(args, casePathDict, caseLoc):
    """Plot per-wing, per-variable, per-y scatter overlays for pressure and profile coordinates."""
    overlayCaseMap = resolveOverlayCaseMap(casePathDict, args.wingCases, os.getcwd())
    if not overlayCaseMap:
        print('\tNo valid cases found for wing pressure plotting.')
        return

    colors = plt.get_cmap('tab10').colors
    markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '<', '>', '*']

    for case, caseInfo in casePathDict.items():
        casePath = caseInfo['path']
        outputDir = os.path.join(casePath, 'postProcessing', 'wingPressure')
        if not os.path.isdir(outputDir):
            print('\tNo wing pressure directory found for %s, skipping output location.' % case)
            continue

        for wingName in ('frontWing', 'rearWing'):
            caseWingData = {}
            variableTargets = set()
            for idxCase, (caseLabel, thisCasePath) in enumerate(overlayCaseMap.items()):
                variableFileMap = discoverWingVariableCsvs(thisCasePath, wingName)
                if not variableFileMap:
                    continue
                for variable in variableFileMap.keys():
                    variableTargets.add(variable)
                caseWingData[caseLabel] = {
                    'path': thisCasePath,
                    'index': idxCase,
                    'variableFiles': variableFileMap,
                }

            if not caseWingData or not variableTargets:
                print('\tNo %s wing-pressure CSV files found for selected cases, skipping.' % wingName)
                continue

            selectedVariables = selectWingPressureVariables(
                args.wingVariables,
                sorted(variableTargets),
                wingName,
            )
            if not selectedVariables:
                continue

            for variable in selectedVariables:
                yTargets = set()
                for caseData in caseWingData.values():
                    yMap = caseData['variableFiles'].get(variable, {})
                    for yKey in yMap.keys():
                        try:
                            yTargets.add(float(yKey))
                        except ValueError:
                            pass

                if not yTargets:
                    continue

                for yVal in sorted(yTargets):
                    yKey = format(float(yVal), '.8g')
                    fig, axes = plt.subplots(2, 1, figsize=[10, 16], frameon=True)
                    pressureAx, profileAx = axes
                    plottedPressure = 0
                    plottedProfile = 0

                    for caseLabel, caseData in caseWingData.items():
                        csvPath = caseData['variableFiles'].get(variable, {}).get(yKey)
                        if csvPath is None:
                            continue
                        try:
                            df = pd.read_csv(csvPath)
                        except Exception as error:
                            print('\tWARNING! Unable to read %s: %s' % (csvPath, error))
                            continue

                        color = colors[caseData['index'] % len(colors)]
                        marker = markers[caseData['index'] % len(markers)]

                        if {'x', variable}.issubset(df.columns):
                            pData = df.dropna(subset=['x', variable]).sort_values('x')
                            if not pData.empty:
                                pressureAx.scatter(pData['x'], pData[variable], s=18,
                                                   marker=marker, color=color,
                                                   label=caseLabel, alpha=0.9)
                                plottedPressure += 1

                        if {'x', 'z'}.issubset(df.columns):
                            profileData = df.dropna(subset=['x', 'z']).sort_values('x')
                            if not profileData.empty:
                                profileAx.scatter(profileData['x'], profileData['z'], s=18,
                                                  marker=marker, color=color,
                                                  label=caseLabel, alpha=0.9)
                                plottedProfile += 1

                    if not plottedPressure and not plottedProfile:
                        plt.close(fig)
                        continue

                    pressureAx.set_xlabel('x (m)')
                    pressureAx.set_ylabel(variable)
                    pressureAx.set_title('%s %s (y = %+.4g m)' % (wingName, variable, yVal))
                    pressureAx.grid(True, alpha=0.3)
                   
                    profileAx.set_xlabel('x (m)')
                    profileAx.set_ylabel('z (m)')
                    profileAx.set_title('%s profile (x-z points, y = %+.4g m)' % (wingName, yVal))
                    profileAx.grid(True, alpha=0.3)
                    profileAx.set_aspect('equal', adjustable='box')

                    if plottedPressure:
                        pressureAx.legend(loc='best', fontsize=8)
                    if plottedProfile:
                        profileAx.legend(loc='best', fontsize=8)

                    overlayTag = sanitizeOverlayTag(overlayCaseMap.keys())
                    outputPath = os.path.join(
                        outputDir,
                        '%s_%s_overlay_%s_y_%s.%s' %
                        (wingName, variable, overlayTag, yKey, args.saveFormat)
                    )
                    fig.suptitle('%s %s overlays at y = %+.4g m' % (wingName, variable, yVal), fontsize=11)
                    fig.tight_layout()
                    fig.savefig(outputPath, dpi=300, bbox_inches='tight')
                    plt.close(fig)
                    print('\tWrote %s' % outputPath)


def wingPressureYValue(csvPath):
    """Extract the y coordinate from a wing CSV filename for sorting."""
    stem = os.path.splitext(os.path.basename(csvPath))[0]
    try:
        return float(stem.rsplit('_', 1)[-1])
    except ValueError:
        return 0.0


def resolveOverlayCaseMap(primaryCasePathDict, overlayCaseArgs, cwd):
    """Resolve case labels to case paths for primary and overlay wing-pressure plots."""
    caseMap = OrderedDict()

    for caseLabel, caseInfo in primaryCasePathDict.items():
        casePath = os.path.abspath(caseInfo['path'])
        if os.path.isdir(casePath):
            caseMap[caseLabel] = casePath

    if not overlayCaseArgs:
        return caseMap

    parentPath = os.path.dirname(os.path.abspath(cwd))
    for token in overlayCaseArgs:
        tokenPath = os.path.expanduser(token)
        candidates = []
        if os.path.isabs(tokenPath):
            candidates.append(tokenPath)
        else:
            candidates.append(os.path.abspath(tokenPath))
            candidates.append(os.path.abspath(os.path.join(parentPath, tokenPath)))

        resolved = None
        for candidate in candidates:
            if os.path.isdir(candidate):
                resolved = candidate
                break

        if resolved is None:
            print('\tWARNING! Could not resolve overlay case %s; skipping.' % token)
            continue

        duplicate = False
        for existingPath in caseMap.values():
            if os.path.abspath(existingPath) == os.path.abspath(resolved):
                duplicate = True
                break
        if duplicate:
            continue

        overlayLabel = os.path.basename(os.path.abspath(resolved))
        label = overlayLabel
        counter = 2
        while label in caseMap:
            label = '%s_%d' % (overlayLabel, counter)
            counter += 1
        caseMap[label] = resolved

    return caseMap


def sanitizeOverlayTag(caseLabels):
    """Create a compact filename-safe tag from plotted case labels."""
    labels = list(caseLabels)
    if not labels:
        return 'none'
    if len(labels) == 1:
        base = labels[0]
    else:
        base = '%s_plus_%d' % (labels[0], len(labels) - 1)
    return re.sub(r'[^A-Za-z0-9_\-\.]+', '_', base)


def discoverWingVariableCsvs(casePath, wingName):
    """Return mapping of variable -> yKey -> csvPath for wing-pressure exports."""
    pressureDir = os.path.join(casePath, 'postProcessing', 'wingPressure')
    if not os.path.isdir(pressureDir):
        return {}

    csvFiles = glob.glob(os.path.join(pressureDir, '%s_*.csv' % wingName))
    variableMap = {}
    coordColumns = {'x', 'y', 'z'}

    for csvPath in csvFiles:
        stem = os.path.splitext(os.path.basename(csvPath))[0]
        prefix = '%s_' % wingName
        if not stem.startswith(prefix):
            continue
        tail = stem[len(prefix):]

        # New combined format: <wingName>_<y>.csv (variables are columns in the file).
        isCombined = False
        try:
            float(tail)
            isCombined = True
        except ValueError:
            isCombined = False

        if isCombined:
            try:
                headerOnly = pd.read_csv(csvPath, nrows=0)
            except Exception:
                continue
            yKey = format(wingPressureYValue(csvPath), '.8g')
            for column in headerOnly.columns:
                if column in coordColumns:
                    continue
                variableMap.setdefault(column, {})[yKey] = csvPath
            continue

        # Legacy format: <wingName>_<variable>_<y>.csv
        if '_' not in tail:
            continue
        variable, yToken = tail.rsplit('_', 1)
        if not variable or not yToken:
            continue
        yKey = format(wingPressureYValue(csvPath), '.8g')
        variableMap.setdefault(variable, {})[yKey] = csvPath
    return variableMap


def selectWingPressureVariables(requestedVariables, availableVariables, wingName):
    """Resolve requested wing variables; default to all available when none requested."""
    if not requestedVariables:
        return availableVariables

    requested = []
    for token in requestedVariables:
        for part in token.replace(',', ' ').split():
            if part and part not in requested:
                requested.append(part)

    selected = [variable for variable in requested if variable in availableVariables]
    missing = [variable for variable in requested if variable not in availableVariables]
    if missing:
        print('\tWARNING! %s missing requested wing variables: %s' %
              (wingName, ', '.join(missing)))

    if not selected:
        print('\tWARNING! No requested variables available for %s. Available: %s' %
              (wingName, ', '.join(availableVariables)))
    return selected


def isCaseComplete(casePath):
    #a case is done when the solver wrote a standalone "End" line in its log. more robust
    #than comparing latest time vs endTime (residualControl can stop steady early)
    for logName in ('log.simpleFoam', 'log.pisoFoam', 'log.SRFSimpleFoam', 'log.SRFPimpleFoam'):
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
    #cornering/per-corner descriptors for a case as an ordered dict to extend the summary.
    #flag/radius/dir from this case's caseSetup, ride-height change + steer from the parent
    #rideHeights_updated.csv. non-cornering falls back to N/A
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
                        ('RH_FL', 'wheel_fl'),
                        ('RH_FR', 'wheel_fr'),
                        ('RH_RL', 'wheel_rl'),
                        ('RH_RR', 'wheel_rr'),
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


def discoverRideHeightChildCases(parentPath, parentCaseName):
    """Return child case directory names matching parentCaseName_# pattern."""
    pattern = re.compile(r'^%s_\d+$' % re.escape(parentCaseName))
    children = []
    try:
        for entry in os.listdir(parentPath):
            if not pattern.match(entry):
                continue
            if os.path.isdir(os.path.join(parentPath, entry)):
                children.append(entry)
    except Exception:
        return []
    return sorted(children)


def readChildSummaryCsv(summaryPath):
    """Read summary.csv written as key,value rows and return a dict."""
    try:
        table = pd.read_csv(summaryPath, header=None)
    except Exception:
        return None
    if table.shape[1] < 2:
        return None
    keys = table.iloc[:, 0].astype(str).str.strip()
    vals = table.iloc[:, 1]
    return dict(zip(keys, vals))


def generate_summary():
    caseSetupPath = "%s/fullCaseSetupDict" % (casePath)
    fullCaseSetupDict = configparser.ConfigParser()
    fullCaseSetupDict.optionxform = str
    fullCaseSetupDict.read_file(open(caseSetupPath))

    parentCaseName = os.path.basename(os.getcwd())
    childCases = discoverRideHeightChildCases(os.getcwd(), parentCaseName)

    if len(childCases) > 0:
        print('Detected ride height map by child directories, averaging child summaries!')
        childRows = []
        for child in childCases:
            childPath = os.path.join(os.getcwd(), child)
            if not isCaseComplete(childPath):
                print('\tWARNING! Child case %s is incomplete, skipping.' % child)
                continue

            summaryPath = os.path.join(childPath, 'summary.csv')
            if not os.path.isfile(summaryPath):
                print('\tWARNING! Child case %s missing summary.csv, skipping.' % child)
                continue

            summaryDict = readChildSummaryCsv(summaryPath)
            if not summaryDict:
                print('\tWARNING! Child case %s has unreadable summary.csv, skipping.' % child)
                continue

            childRows.append(summaryDict)

        if len(childRows) < 1:
            print('\tNo valid child summaries found; skipping parent summary.')
            return

        summaryFrame = pd.DataFrame(childRows)
        numericFrame = summaryFrame.apply(pd.to_numeric, errors='coerce')
        meanNumeric = numericFrame.mean(axis=0, skipna=True)
        first = childRows[0]

        rowNames = ['Job','Trial','Solver','Version','Run Date','Solve Time','Num. Cells','Mesher','Symmetry','Ref. Area (m^2)','Iterations','Simulation Type','Moving Ground','Rotating Wheels','Turbulence Model','Velocity','Yaw','Cd','Cl','Cl/Cd','%Front','Cd CI','Cl CI']
        data = [
            job,
            parentCaseName,
            first.get('Solver', 'N/A'),
            first.get('Version', 'N/A'),
            'N/A',
            'N/A',
            'N/A',
            first.get('Mesher', 'N/A'),
            str(first.get('Symmetry', 'N/A')).lower(),
            meanNumeric.get('Ref. Area (m^2)', np.nan),
            meanNumeric.get('Iterations', np.nan),
            str(first.get('Simulation Type', 'N/A')).lower(),
            first.get('Moving Ground', 'N/A'),
            first.get('Rotating Wheels', 'N/A'),
            first.get('Turbulence Model', 'N/A'),
            meanNumeric.get('Velocity', np.nan),
            meanNumeric.get('Yaw', np.nan),
            meanNumeric.get('Cd', np.nan),
            meanNumeric.get('Cl', np.nan),
            meanNumeric.get('Cl/Cd', np.nan),
            meanNumeric.get('%Front', np.nan),
            meanNumeric.get('Cd CI', np.nan),
            meanNumeric.get('Cl CI', np.nan),
        ]

        baseSet = set(rowNames)
        for column in meanNumeric.index:
            if column in baseSet:
                continue
            value = meanNumeric[column]
            if pd.isna(value):
                continue
            rowNames.append(column)
            data.append(value)

        summary = pd.DataFrame(columns=rowNames)
        summary.loc[-1] = data
        print("\n\n")
        for col in summary.columns:
            print('{:>100s}{:>30s}'.format(col, str(summary[col].values[0])))

        summary = summary.transpose()
        summary.to_csv("%s/summary.csv" % (casePath), header=False)
        return

    case = os.path.basename(casePath)
    coeffFiles = getCoeffPaths(casePath)
    partsDict = {}
    for part in coeffFiles:
        if part != 'all':
            partsDict[part], avgDataArray = averageCoeffs(fullCaseSetupDict, case, part, coeffFiles)
    avgData, allDataArray = averageCoeffs(fullCaseSetupDict, case, 'all', coeffFiles)

    numCells, mesher, sym = cellCount(fullCaseSetupDict, casePath, case)
    inletMag, lastTime, yaw, movingGround, rotatingWheels, simType, turbModel = bcParser(fullCaseSetupDict, path, case)
    runDate, runTime, version, solver = getOfVersion(casePath)
    refArea = float(fullCaseSetupDict['BC_SETUP']['REFAREA'][0])

    rowNames = ['Job','Trial','Solver','Version','Run Date','Solve Time','Num. Cells','Mesher','Symmetry','Ref. Area (m^2)','Iterations','Simulation Type','Moving Ground','Rotating Wheels','Turbulence Model','Velocity','Yaw','Cd','Cl','Cl/Cd','%Front','Cd CI','Cl CI']
    data = [job, case, solver, version, runDate, runTime, numCells, mesher, sym.lower(), refArea, avgData['endTime'], simType.lower(), movingGround, rotatingWheels, turbModel, inletMag, yaw, avgData['cd'], avgData['cl'], avgData['cl/cd'], avgData['cop'], avgData['cd_ci'], avgData['cl_ci']]

    corneringInfo = getCorneringInfo(fullCaseSetupDict, casePath, case)
    for label, value in corneringInfo.items():
        rowNames.append(label)
        data.append(value)

    for part in partsDict.keys():
        partVarDict = {'CL': 'cl', 'CD': 'cd'}
        for varkey in partVarDict.keys():
            rowNames.append(part + ' ' + varkey)
            data.append(partsDict[part][partVarDict[varkey]])

    try:
        porousData = getPorousData(path, case)
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
        print('{:>100s}{:>30s}'.format(col, str(summary[col].values[0])))

    summary = summary.transpose()
    summary.to_csv("%s/%s/summary.csv" % (path, case), header=False)

main()
