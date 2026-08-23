import os
import sys
import numpy as np
import re
import pandas as pd
import configparser
import argparse
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.interpolate import griddata
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
    parser.add_argument('--sensitivityPlots', action='store_true',
                       help='Plot ride-height/yaw/cornering sensitivity sweeps for a ride-height mapping parent case')
    parser.add_argument('--includeSideForce', action='store_true',
                       help='Include Cs(f)/Cs(r) in sensitivity sweep plots (off by default)')

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

    if args.sensitivityPlots:
        plotRideHeightSensitivity(casePath, includeSideForce=args.includeSideForce)

    if not args.summary and not args.forces and not args.wingPlots and not args.sensitivityPlots:
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


def formatSummaryNumericValues(summaryDf, decimals=3):
    """Format all numeric summary fields to a fixed decimal precision."""
    if summaryDf.empty:
        return summaryDf

    rowIdx = summaryDf.index[0]
    for col in summaryDf.columns:
        try:
            numericValue = float(summaryDf.at[rowIdx, col])
            if np.isfinite(numericValue):
                summaryDf.at[rowIdx, col] = ('%%.%df' % decimals) % numericValue
        except Exception:
            pass
    return summaryDf


#ride-height map sensitivity sweeps: groups of ride-height-CSV columns that can each act as a
#single "swept" variable. 'cols' lists the underlying CSV columns that belong to the group;
#a group counts as "changing" if its columns vary within the rows being considered (columns of
#the same group moving together, e.g. fl==fr, still count as ONE variable). Because the
#ride-height map can be a full grid (front x rear x yaw x ...), a group is plotted as a sweep for
#each FIXED-COMBO subset of the other groups (incl. STEER_ONLY_GROUP) where that group varies
#and at least MIN_SWEEP_POINTS rows are available -- i.e. one figure per single-variable sweep,
#even if several such sweeps of the same group exist at different fixed values of everything else.
MIN_SWEEP_POINTS = 3
RIDE_HEIGHT_SWEEP_GROUPS = OrderedDict([
    ('Front Ride Height', {'cols': ['fl', 'fr'], 'xLabel': 'Front Ride Height (avg fl/fr)'}),
    ('Rear Ride Height', {'cols': ['rl', 'rr'], 'xLabel': 'Rear Ride Height (avg rl/rr)'}),
    ('Yaw', {'cols': ['yaw'], 'xLabel': 'Yaw (deg)'}),
    ('Corner', {'cols': ['corner_radius', 'corner_dir'], 'xLabel': 'Corner Radius (m)'}),
])
#constancy-only: never plotted as its own sweep, but must stay constant for any of the
#groups above to be considered a clean single-variable sweep
RIDE_HEIGHT_STEER_GROUP = {'cols': ['steer', 'steer_deg', 'steer_angle']}

#default force/moment coefficients plotted for every sensitivity sweep; Cs(f)/Cs(r) are only
#added when includeSideForce=True (see plotRideHeightSensitivity)
DEFAULT_SWEEP_METRICS = ['Cd', 'Cl', 'Cl(f)', 'Cl(r)']
OPTIONAL_SWEEP_METRICS = ['Cs(f)', 'Cs(r)']


def loadRideHeightMap(casePath):
    """Load the parent case's rideHeights_updated.csv, or None if not present/unreadable."""
    rhPath = os.path.join(casePath, 'rideHeights_updated.csv')
    if not os.path.isfile(rhPath):
        return None
    try:
        return pd.read_csv(rhPath)
    except Exception as e:
        print('\tUnable to read %s: %s' % (rhPath, e))
        return None


def buildSweepDataset(casePath, rhMap):
    """Merge each ride-height-map row with its child case's summary.csv metrics.

    Skips child cases that are incomplete, missing summary.csv, or unreadable (mirrors
    the child-skip behavior used when averaging the parent summary).
    """
    if 'caseName' not in rhMap.columns:
        print('\tWARNING! rideHeights_updated.csv has no caseName column; skipping sensitivity plots.')
        return None

    metricKeys = DEFAULT_SWEEP_METRICS + OPTIONAL_SWEEP_METRICS
    rows = []
    for _, row in rhMap.iterrows():
        caseName = str(row.get('caseName', '')).strip()
        if not caseName:
            continue
        childPath = os.path.join(casePath, caseName)
        if not isCaseComplete(childPath):
            print('\tWARNING! Child case %s is incomplete, skipping for sensitivity plots.' % caseName)
            continue

        summaryPath = os.path.join(childPath, 'summary.csv')
        if not os.path.isfile(summaryPath):
            print('\tWARNING! Child case %s missing summary.csv, skipping for sensitivity plots.' % caseName)
            continue

        summaryDict = readChildSummaryCsv(summaryPath)
        if not summaryDict:
            print('\tWARNING! Child case %s has unreadable summary.csv, skipping for sensitivity plots.' % caseName)
            continue

        entry = row.to_dict()
        for metric in metricKeys:
            try:
                entry[metric] = float(summaryDict.get(metric, np.nan))
            except Exception:
                entry[metric] = np.nan
        rows.append(entry)

    if len(rows) < 2:
        print('\tNot enough complete child cases with data to build sensitivity plots.')
        return None

    return pd.DataFrame(rows)


def detectRideHeightSweeps(df, minPoints=MIN_SWEEP_POINTS):
    """Return (sweeps, activeGroups).

    activeGroups maps sweep-group name -> {'cols', 'xLabel'} for every group whose columns
    are present in df.

    sweeps is a list of {'name': groupName, 'df': subsetDataFrame} entries. Each entry is a
    fixed-combo subset of df -- rows where every OTHER group (including the steer-only
    constancy group) has a single constant value -- within which the named group's columns
    vary across at least `minPoints` rows. Because the ride-height map can be a full grid
    (e.g. front x rear combinations), the SAME group can produce several sweep entries at
    different fixed values of the other variables (front sweep with rear held at 0, another
    front sweep with rear held at -0.01, etc); each is returned separately so it can be
    plotted on its own figure.
    """
    activeGroups = OrderedDict()
    for name, spec in RIDE_HEIGHT_SWEEP_GROUPS.items():
        cols = [c for c in spec['cols'] if c in df.columns]
        if cols:
            activeGroups[name] = {'cols': cols, 'xLabel': spec['xLabel']}

    steerCols = [c for c in RIDE_HEIGHT_STEER_GROUP['cols'] if c in df.columns]

    def effectiveValue(frame, cols):
        #columns within a group (e.g. fl/fr) move together for a pure ride-height sweep;
        #average them into one representative value per row so both columns changing IN STEP
        #still counts as a single variable (fl == fr -> treated as one axis).
        return frame[cols].mean(axis=1).round(9)

    sweeps = []
    for name, spec in activeGroups.items():
        otherCols = []
        for otherName, otherSpec in activeGroups.items():
            if otherName != name:
                otherCols.extend(otherSpec['cols'])
        otherCols.extend(steerCols)
        otherCols = [c for c in dict.fromkeys(otherCols) if c in df.columns]

        if otherCols:
            groupKey = df[otherCols].round(9).apply(tuple, axis=1)
        else:
            groupKey = pd.Series(0, index=df.index)

        for _, subIdx in df.groupby(groupKey).groups.items():
            subDf = df.loc[subIdx].reset_index(drop=True)
            if len(subDf) < minPoints:
                continue
            if effectiveValue(subDf, spec['cols']).nunique(dropna=True) <= 1:
                continue
            sweeps.append({'name': name, 'df': subDf})

    return sweeps, activeGroups


def _sweepFixedParts(sweepName, df, activeGroups):
    """List of 'col=value' strings for every OTHER group/steer column held constant in df."""
    fixedParts = []
    for otherName, otherSpec in activeGroups.items():
        if otherName == sweepName:
            continue
        for col in otherSpec['cols']:
            if df[col].nunique(dropna=True) <= 1:
                fixedParts.append('%s=%s' % (col, df[col].iloc[0]))

    steerCols = [c for c in RIDE_HEIGHT_STEER_GROUP['cols'] if c in df.columns]
    for col in steerCols:
        if df[col].nunique(dropna=True) <= 1:
            fixedParts.append('%s=%s' % (col, df[col].iloc[0]))
            break
    return fixedParts


def buildSweepTitle(sweepName, df, activeGroups):
    """Short title naming just the sweep; fixed configuration is shown in a side panel instead
    (see plotRideHeightSweep) so long lists of held-constant values don't overrun the title."""
    return '%s Sensitivity Sweep' % sweepName


def getSweepXValues(sweepName, groupSpec, df):
    """Return (xValues array, xLabel) for the swept variable."""
    if sweepName == 'Front Ride Height' and {'fl', 'fr'}.issubset(df.columns):
        return df[['fl', 'fr']].mean(axis=1).to_numpy(), groupSpec['xLabel']
    if sweepName == 'Rear Ride Height' and {'rl', 'rr'}.issubset(df.columns):
        return df[['rl', 'rr']].mean(axis=1).to_numpy(), groupSpec['xLabel']
    if sweepName == 'Corner':
        if 'corner_radius' in df.columns and df['corner_radius'].nunique(dropna=True) > 1:
            return df['corner_radius'].to_numpy(), 'Corner Radius (m)'
        if 'corner_dir' in df.columns:
            return df['corner_dir'].astype(str).to_numpy(), 'Corner Direction'
    col = groupSpec['cols'][0]
    return df[col].to_numpy(), groupSpec.get('xLabel', col)


def plotRideHeightSweep(df, sweepName, groupSpec, activeGroups, casePath, includeSideForce=False,
                         fileSuffix=''):
    """Plot one figure (scatter + polynomial curve fit, sorted by x) for a single-variable sweep.

    The swept variable's held-constant siblings are listed in a side panel (not the title) so
    long configuration lists don't overrun the figure width.
    """
    xValues, xLabel = getSweepXValues(sweepName, groupSpec, df)

    metrics = list(DEFAULT_SWEEP_METRICS)
    if includeSideForce:
        metrics += OPTIONAL_SWEEP_METRICS
    metrics = [m for m in metrics if m in df.columns and df[m].notna().any()]
    if not metrics:
        print('\tNo usable force coefficient columns for %s sweep, skipping plot.' % sweepName)
        return

    xArr = np.asarray(xValues)
    isNumericX = np.issubdtype(xArr.dtype, np.number)
    sortOrder = np.argsort(xArr) if isNumericX else np.argsort(xArr.astype(str))

    fig, axes = plt.subplots(len(metrics), 1, figsize=(7, 3 * len(metrics)), sharex=True, squeeze=False)
    axes = axes[:, 0]

    for ax, metric in zip(axes, metrics):
        yValues = df[metric].to_numpy()
        xSorted = xArr[sortOrder]
        ySorted = yValues[sortOrder]
        ax.scatter(xSorted, ySorted, marker='o', zorder=3)

        if isNumericX:
            validMask = np.isfinite(xSorted.astype(float)) & np.isfinite(ySorted.astype(float))
            nValid = int(np.count_nonzero(validMask))
            #degree scales with available points but stays low-order (avoid overfitting a
            #handful of ride-height points); need at least degree+1 points to fit.
            degree = min(3, nValid - 1) if nValid > 1 else 0
            if degree >= 1 and np.unique(xSorted[validMask]).size > degree:
                coeffs = np.polyfit(xSorted[validMask].astype(float), ySorted[validMask].astype(float), degree)
                xFit = np.linspace(xSorted[validMask].min(), xSorted[validMask].max(), 100)
                ax.plot(xFit, np.polyval(coeffs, xFit), linestyle='-', zorder=2)
            else:
                ax.plot(xSorted, ySorted, linestyle='-', zorder=2)
        else:
            ax.plot(xSorted, ySorted, linestyle='-', zorder=2)

        ax.set_ylabel(metric)
        ax.grid(True)

    axes[-1].set_xlabel(xLabel)
    fig.suptitle(buildSweepTitle(sweepName, df, activeGroups))

    fixedParts = _sweepFixedParts(sweepName, df, activeGroups)
    if fixedParts:
        sideText = 'Fixed configuration:\n' + '\n'.join(fixedParts)
        fig.subplots_adjust(right=0.72)
        fig.text(0.75, 0.5, sideText, va='center', ha='left', fontsize=8,
                  bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray'))
    else:
        fig.tight_layout()

    outputDir = os.path.join(casePath, 'postProcessing', 'sensitivityPlots')
    os.makedirs(outputDir, exist_ok=True)
    outFile = os.path.join(outputDir, '%s%s_sweep.png' % (sweepName.replace(' ', ''), fileSuffix))
    fig.savefig(outFile, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('\tSaved %s sensitivity sweep: %s' % (sweepName, outFile))


#minimum number of rows / unique front-rear combinations required to build a FRH-RRH contour
MIN_CONTOUR_POINTS = 4
#other columns besides front/rear ride height that must be held constant for a contour group
#(roll itself is enforced separately -- only roll==0 rows are considered at all)
FRH_RRH_CONTOUR_CONSTANT_GROUPS = ['Yaw', 'Corner']


def buildFrhRrhContourGroups(df, activeGroups, minPoints=MIN_CONTOUR_POINTS):
    """Return a list of fixed-combo subsets of df suitable for a FRH-vs-RRH contour plot.

    Requires 'roll' == 0 (within tolerance) and both front (fl/fr) and rear (rl/rr) ride
    height columns present and varying with at least 2 unique values each. Every other
    tracked group (Yaw, Corner) and the steer-only constancy group must be held constant
    within each returned subset, matching the same fixed-combo logic used for 1D sweeps.
    """
    if 'roll' not in df.columns:
        return []
    if not {'fl', 'fr', 'rl', 'rr'}.issubset(df.columns):
        return []

    rollFiltered = df[np.isclose(df['roll'].astype(float), 0.0, atol=1e-6)]
    if rollFiltered.empty:
        return []

    steerCols = [c for c in RIDE_HEIGHT_STEER_GROUP['cols'] if c in rollFiltered.columns]
    otherCols = []
    for groupName in FRH_RRH_CONTOUR_CONSTANT_GROUPS:
        spec = activeGroups.get(groupName)
        if spec:
            otherCols.extend(spec['cols'])
    otherCols.extend(steerCols)
    otherCols = [c for c in dict.fromkeys(otherCols) if c in rollFiltered.columns]

    if otherCols:
        groupKey = rollFiltered[otherCols].round(9).apply(tuple, axis=1)
    else:
        groupKey = pd.Series(0, index=rollFiltered.index)

    contourGroups = []
    for _, subIdx in rollFiltered.groupby(groupKey).groups.items():
        subDf = rollFiltered.loc[subIdx].reset_index(drop=True)
        if len(subDf) < minPoints:
            continue
        frontVal = subDf[['fl', 'fr']].mean(axis=1).round(9)
        rearVal = subDf[['rl', 'rr']].mean(axis=1).round(9)
        if frontVal.nunique() < 2 or rearVal.nunique() < 2:
            continue
        contourGroups.append(subDf)

    return contourGroups


def plotFrhRrhContour(subDf, activeGroups, casePath, includeSideForce=False, fileSuffix=''):
    """Plot a Front-Ride-Height (x) vs Rear-Ride-Height (y) contour, cubic-interpolated, one
    subplot per force coefficient metric, for a fixed-combo subset with roll held at 0."""
    frontVal = subDf[['fl', 'fr']].mean(axis=1).to_numpy(dtype=float)
    rearVal = subDf[['rl', 'rr']].mean(axis=1).to_numpy(dtype=float)

    metrics = list(DEFAULT_SWEEP_METRICS)
    if includeSideForce:
        metrics += OPTIONAL_SWEEP_METRICS
    metrics = [m for m in metrics if m in subDf.columns and subDf[m].notna().any()]
    if not metrics:
        print('\tNo usable force coefficient columns for FRH-RRH contour, skipping plot.')
        return

    gridX, gridY = np.meshgrid(
        np.linspace(frontVal.min(), frontVal.max(), 100),
        np.linspace(rearVal.min(), rearVal.max(), 100),
    )

    fig, axes = plt.subplots(1, len(metrics), figsize=(6 * len(metrics), 5), squeeze=False)
    axes = axes[0, :]

    for ax, metric in zip(axes, metrics):
        zValues = subDf[metric].to_numpy(dtype=float)
        gridZ = griddata((frontVal, rearVal), zValues, (gridX, gridY), method='cubic')

        contourf = ax.contourf(gridX, gridY, gridZ, levels=20, cmap='viridis')
        ax.contour(gridX, gridY, gridZ, levels=20, colors='black', linewidths=0.4, alpha=0.5)
        ax.scatter(frontVal, rearVal, c='white', edgecolors='black', s=25, zorder=3)
        fig.colorbar(contourf, ax=ax, label=metric)
        ax.set_xlabel('Front Ride Height (avg fl/fr)')
        ax.set_ylabel('Rear Ride Height (avg rl/rr)')
        ax.set_title(metric)

    fig.suptitle('Front vs Rear Ride Height Sensitivity Contour (roll=0)')

    fixedParts = []
    for groupName in FRH_RRH_CONTOUR_CONSTANT_GROUPS:
        spec = activeGroups.get(groupName)
        if not spec:
            continue
        for col in spec['cols']:
            if subDf[col].nunique(dropna=True) <= 1:
                fixedParts.append('%s=%s' % (col, subDf[col].iloc[0]))
    fixedParts.append('roll=0')
    steerCols = [c for c in RIDE_HEIGHT_STEER_GROUP['cols'] if c in subDf.columns]
    for col in steerCols:
        if subDf[col].nunique(dropna=True) <= 1:
            fixedParts.append('%s=%s' % (col, subDf[col].iloc[0]))
            break

    if fixedParts:
        sideText = 'Fixed configuration:\n' + '\n'.join(fixedParts)
        fig.subplots_adjust(top=0.8)
        fig.text(0.5, 0.90, sideText, va='center', ha='center', fontsize=8,
                  bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray'))

    outputDir = os.path.join(casePath, 'postProcessing', 'sensitivityPlots')
    os.makedirs(outputDir, exist_ok=True)
    outFile = os.path.join(outputDir, 'FrontRearRideHeight_contour%s.png' % fileSuffix)
    fig.savefig(outFile, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('\tSaved FRH-RRH sensitivity contour: %s' % outFile)


def plotFrhRrhContours(df, casePath, includeSideForce=False):
    """Detect and plot FRH-vs-RRH contour(s) for fixed-combo subsets with roll held at 0."""
    _, activeGroups = detectRideHeightSweeps(df)
    contourGroups = buildFrhRrhContourGroups(df, activeGroups)
    if not contourGroups:
        print('\tNo roll=0 Front/Rear Ride Height grid detected; skipping FRH-RRH contour plot.')
        return

    usedSuffixes = {}
    for subDf in contourGroups:
        fixedParts = []
        for groupName in FRH_RRH_CONTOUR_CONSTANT_GROUPS:
            spec = activeGroups.get(groupName)
            if not spec:
                continue
            for col in spec['cols']:
                if subDf[col].nunique(dropna=True) <= 1:
                    fixedParts.append('%s%s' % (col, subDf[col].iloc[0]))
        slug = '_'.join(fixedParts).replace(' ', '')
        fileSuffix = ('_%s' % slug) if slug else ''

        usedSuffixes[fileSuffix] = usedSuffixes.get(fileSuffix, 0) + 1
        if usedSuffixes[fileSuffix] > 1:
            fileSuffix = '%s_%d' % (fileSuffix, usedSuffixes[fileSuffix])

        plotFrhRrhContour(subDf, activeGroups, casePath, includeSideForce=includeSideForce,
                           fileSuffix=fileSuffix)


def plotRideHeightSensitivity(casePath, includeSideForce=False):
    """Detect and plot single-variable sensitivity sweeps for a ride-height mapping parent case.

    Only runs for parent cases (those with child dirs matching caseName_#); does nothing for
    plain single cases or when run from inside a child case.
    """
    parentCaseName = os.path.basename(casePath)
    childCases = discoverRideHeightChildCases(casePath, parentCaseName)
    if not childCases:
        print('\tNo ride-height child cases detected; skipping sensitivity plots.')
        return

    rhMap = loadRideHeightMap(casePath)
    if rhMap is None:
        print('\tNo rideHeights_updated.csv found; skipping sensitivity plots.')
        return

    df = buildSweepDataset(casePath, rhMap)
    if df is None:
        return

    sweeps, activeGroups = detectRideHeightSweeps(df)
    if not sweeps:
        print('\tNo single-variable sweeps detected among ride-height child cases.')
    else:
        usedSuffixes = {}
        for sweep in sweeps:
            sweepName = sweep['name']
            subDf = sweep['df']

            fixedParts = _sweepFixedParts(sweepName, subDf, activeGroups)
            slug = '_'.join(p.replace('=', '') for p in fixedParts).replace(' ', '')
            fileSuffix = ('_%s' % slug) if slug else ''

            #disambiguate on the rare chance two sub-sweeps of the same group produce the same slug
            key = (sweepName, fileSuffix)
            usedSuffixes[key] = usedSuffixes.get(key, 0) + 1
            if usedSuffixes[key] > 1:
                fileSuffix = '%s_%d' % (fileSuffix, usedSuffixes[key])

            plotRideHeightSweep(subDf, sweepName, activeGroups[sweepName], activeGroups, casePath,
                                 includeSideForce=includeSideForce, fileSuffix=fileSuffix)

    plotFrhRrhContours(df, casePath, includeSideForce=includeSideForce)



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

        rowNames = ['Job','Trial','Solver','Version','Run Date','Solve Time','Num. Cells','Mesher','Symmetry','Ref. Area (m^2)','Iterations','Simulation Type','Moving Ground','Rotating Wheels','Turbulence Model','Velocity','Yaw','Cd','Cl','Cl/Cd','%Front','Cd CI','Cl CI','Cl(f)','Cl(r)','Cs(f)','Cs(r)']
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
            meanNumeric.get('Cl(f)', meanNumeric.get('clf', np.nan)),
            meanNumeric.get('Cl(r)', meanNumeric.get('clr', np.nan)),
            meanNumeric.get('Cs(f)', meanNumeric.get('csf', np.nan)),
            meanNumeric.get('Cs(r)', meanNumeric.get('csr', np.nan)),
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
        summary = formatSummaryNumericValues(summary, decimals=3)
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

    rowNames = ['Job','Trial','Solver','Version','Run Date','Solve Time','Num. Cells','Mesher','Symmetry','Ref. Area (m^2)','Iterations','Simulation Type','Moving Ground','Rotating Wheels','Turbulence Model','Velocity','Yaw','Cd','Cl','Cl/Cd','%Front','Cd CI','Cl CI','Cl(f)','Cl(r)','Cs(f)','Cs(r)']
    data = [job, case, solver, version, runDate, runTime, numCells, mesher, sym.lower(), refArea, avgData['endTime'], simType.lower(), movingGround, rotatingWheels, turbModel, inletMag, yaw, avgData['cd'], avgData['cl'], avgData['cl/cd'], avgData['cop'], avgData['cd_ci'], avgData['cl_ci'], avgData.get('clf', np.nan), avgData.get('clr', np.nan), avgData.get('csf', np.nan), avgData.get('csr', np.nan)]

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
    summary = formatSummaryNumericValues(summary, decimals=3)
    print("\n\n")
    for col in summary.columns:
        print('{:>100s}{:>30s}'.format(col, str(summary[col].values[0])))

    summary = summary.transpose()
    summary.to_csv("%s/%s/summary.csv" % (path, case), header=False)

main()
