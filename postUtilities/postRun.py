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
from scipy import ndimage
import glob
from collections import OrderedDict
from datetime import date
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
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
    parser.add_argument('--sensitivityCases', nargs='+', default=[],
                       help='Other ride-height mapping parent case paths to overlay on the sensitivity sweep plots for comparison')
    parser.add_argument('--pptReport', action='store_true',
                       help='Generate a PowerPoint (.pptx) report summarizing the trial(s) in --trial')
    parser.add_argument('--addMovies', action='store_true',
                       help='Generate and embed slice movies (requires ffmpeg on PATH) when building a --pptReport; off by default since it can be slow')

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
        plotRideHeightSensitivity(casePath, includeSideForce=args.includeSideForce,
                                   compareCasePaths=args.sensitivityCases)

    if args.pptReport:
        generate_ppt_report(args)

    if not args.summary and not args.forces and not args.wingPlots and not args.sensitivityPlots and not args.pptReport:
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


def plotRideHeightSweep(entries, sweepName, groupSpec, activeGroups, casePath, includeSideForce=False,
                         fileSuffix=''):
    """Plot one figure (scatter + polynomial curve fit, sorted by x) for a single-variable sweep.

    `entries` is a list of (caseLabel, df) pairs. A single-case sweep passes one entry with
    caseLabel=None (no legend shown); a multi-case comparison passes one entry per matched
    parent case (same sweep type + same fixed configuration), each drawn in its own color with
    a legend so sensitivities can be compared across ride-height maps.

    The swept variable's held-constant siblings are listed in a side panel (not the title) so
    long configuration lists don't overrun the figure width.
    """
    metrics = list(DEFAULT_SWEEP_METRICS)
    if includeSideForce:
        metrics += OPTIONAL_SWEEP_METRICS
    metrics = [m for m in metrics if any(m in df.columns and df[m].notna().any() for _, df in entries)]
    if not metrics:
        print('\tNo usable force coefficient columns for %s sweep, skipping plot.' % sweepName)
        return

    fig, axes = plt.subplots(len(metrics), 1, figsize=(7, 3 * len(metrics)), sharex=True, squeeze=False)
    axes = axes[:, 0]

    colorCycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0'])
    showLegend = len(entries) > 1

    xLabel = groupSpec.get('xLabel', sweepName)
    for entryIdx, (caseLabel, df) in enumerate(entries):
        if metrics and not any(m in df.columns and df[m].notna().any() for m in metrics):
            continue
        xValues, entryXLabel = getSweepXValues(sweepName, groupSpec, df)
        xLabel = entryXLabel
        xArr = np.asarray(xValues)
        isNumericX = np.issubdtype(xArr.dtype, np.number)
        sortOrder = np.argsort(xArr) if isNumericX else np.argsort(xArr.astype(str))
        color = colorCycle[entryIdx % len(colorCycle)]

        for ax, metric in zip(axes, metrics):
            if metric not in df.columns:
                continue
            yValues = df[metric].to_numpy()
            xSorted = xArr[sortOrder]
            ySorted = yValues[sortOrder]
            ax.scatter(xSorted, ySorted, marker='o', zorder=3, color=color,
                       label=caseLabel if showLegend else None)

            if isNumericX:
                validMask = np.isfinite(xSorted.astype(float)) & np.isfinite(ySorted.astype(float))
                nValid = int(np.count_nonzero(validMask))
                #degree scales with available points but stays low-order (avoid overfitting a
                #handful of ride-height points); need at least degree+1 points to fit.
                degree = min(3, nValid - 1) if nValid > 1 else 0
                if degree >= 1 and np.unique(xSorted[validMask]).size > degree:
                    coeffs = np.polyfit(xSorted[validMask].astype(float), ySorted[validMask].astype(float), degree)
                    xFit = np.linspace(xSorted[validMask].min(), xSorted[validMask].max(), 100)
                    ax.plot(xFit, np.polyval(coeffs, xFit), linestyle='-', zorder=2, color=color)
                else:
                    ax.plot(xSorted, ySorted, linestyle='-', zorder=2, color=color)
            else:
                ax.plot(xSorted, ySorted, linestyle='-', zorder=2, color=color)

    for ax, metric in zip(axes, metrics):
        ax.set_ylabel(metric)
        ax.grid(True)

    if showLegend:
        axes[0].legend(fontsize=8)

    axes[-1].set_xlabel(xLabel)
    titleDf = entries[0][1]
    fig.suptitle(buildSweepTitle(sweepName, titleDf, activeGroups))

    fixedParts = _sweepFixedParts(sweepName, titleDf, activeGroups)
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
    fig.savefig(outFile, dpi=300, bbox_inches='tight')
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

    nCols = 2
    nRows = int(np.ceil(len(metrics) / float(nCols)))
    fig, axes = plt.subplots(nRows, nCols, figsize=(7.5 * nCols, 6 * nRows), squeeze=False)
    flatAxes = axes.flatten()
    fig.subplots_adjust(wspace=0.55, hspace=0.45)

    for ax, metric in zip(flatAxes, metrics):
        zValues = subDf[metric].to_numpy(dtype=float)
        gridZ = griddata((frontVal, rearVal), zValues, (gridX, gridY), method='cubic')

        contourf = ax.contourf(gridX, gridY, gridZ, levels=20, cmap='viridis')
        ax.contour(gridX, gridY, gridZ, levels=20, colors='black', linewidths=0.4, alpha=0.5)
        ax.scatter(frontVal, rearVal, c='white', edgecolors='black', s=25, zorder=3)
        fig.colorbar(contourf, ax=ax, label=metric, pad=0.03)
        ax.set_xlabel('Front Ride Height (avg fl/fr)')
        ax.set_ylabel('Rear Ride Height (avg rl/rr)')
        ax.set_title(metric, pad=10)

    #hide any unused axes (e.g. only 3 metrics in a 2x2 grid)
    for ax in flatAxes[len(metrics):]:
        ax.set_visible(False)

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
    fig.savefig(outFile, dpi=300, bbox_inches='tight')
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


def loadCaseSweepDataset(path):
    """Build the ride-height sweep dataset for a parent case directory (used for both the
    primary case and any --sensitivityCases comparison cases). Returns None if `path` is not a
    valid ride-height mapping parent (no matching child dirs, or no rideHeights_updated.csv)."""
    path = os.path.abspath(path)
    parentCaseName = os.path.basename(path)
    childCases = discoverRideHeightChildCases(path, parentCaseName)
    if not childCases:
        print('\tNo ride-height child cases detected in %s; skipping.' % path)
        return None

    rhMap = loadRideHeightMap(path)
    if rhMap is None:
        print('\tNo rideHeights_updated.csv found in %s; skipping.' % path)
        return None

    return buildSweepDataset(path, rhMap)


def plotRideHeightSensitivity(casePath, includeSideForce=False, compareCasePaths=None):
    """Detect and plot single-variable sensitivity sweeps for a ride-height mapping parent case.

    Only runs for parent cases (those with child dirs matching caseName_#); does nothing for
    plain single cases or when run from inside a child case.

    If `compareCasePaths` is given, sub-sweeps from those other parent ride-height cases are
    overlaid on the same figure whenever they match the primary case's sweep type AND fixed
    configuration (e.g. both are Front Ride Height sweeps with Rear Ride Height=-0.01), so
    sensitivities can be compared across cases. FRH-RRH contour plots remain single-case.
    """
    df = loadCaseSweepDataset(casePath)
    if df is None:
        print('\tSkipping sensitivity plots for %s.' % casePath)
        return

    caseLabel = os.path.basename(os.path.normpath(casePath))
    sweeps, activeGroups = detectRideHeightSweeps(df)

    #collect comparison cases' sweeps under the same (sweepName, fixedConfig) matching used
    #for the primary case's own sub-sweeps, so all matching sub-sweeps across cases share a plot
    matchedSweeps = OrderedDict()
    for sweep in sweeps:
        sweepName = sweep['name']
        subDf = sweep['df']
        fixedKey = frozenset(_sweepFixedParts(sweepName, subDf, activeGroups))
        matchedSweeps.setdefault((sweepName, fixedKey), []).append((caseLabel, subDf))

    for comparePath in (compareCasePaths or []):
        compareDf = loadCaseSweepDataset(comparePath)
        if compareDf is None:
            continue
        compareLabel = os.path.basename(os.path.normpath(comparePath))
        compareSweeps, compareActiveGroups = detectRideHeightSweeps(compareDf)
        for sweep in compareSweeps:
            sweepName = sweep['name']
            subDf = sweep['df']
            fixedKey = frozenset(_sweepFixedParts(sweepName, subDf, compareActiveGroups))
            key = (sweepName, fixedKey)
            if key in matchedSweeps:
                matchedSweeps[key].append((compareLabel, subDf))

    if not matchedSweeps:
        print('\tNo single-variable sweeps detected among ride-height child cases.')
    else:
        usedSuffixes = {}
        for (sweepName, fixedKey), entries in matchedSweeps.items():
            slug = '_'.join(p.replace('=', '') for p in sorted(fixedKey)).replace(' ', '')
            fileSuffix = ('_%s' % slug) if slug else ''

            #disambiguate on the rare chance two sub-sweeps of the same group produce the same slug
            key = (sweepName, fileSuffix)
            usedSuffixes[key] = usedSuffixes.get(key, 0) + 1
            if usedSuffixes[key] > 1:
                fileSuffix = '%s_%d' % (fileSuffix, usedSuffixes[key])

            plotRideHeightSweep(entries, sweepName, activeGroups[sweepName], activeGroups, casePath,
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


#### PPT Report Generation ####
#field names (matching summary.csv row labels written by generate_summary) shown on the
#trial-setup/BC table and the results table; fields missing from a case's summary.csv are
#shown as 'N/A' rather than being skipped, so tables stay aligned across trials
PPT_INFO_FIELDS = ['Simulation Type', 'Velocity', 'Turbulence Model', 'Moving Ground',
                   'Rotating Wheels', 'Yaw', 'Symmetry', 'Iterations', 'Num. Cells']
PPT_RESULTS_FIELDS = ['Cd', 'Cl', 'Cl/Cd', '%Front', 'Cd CI', 'Cl CI',
                      'Cl(f)', 'Cl(r)', 'Cs(f)', 'Cs(r)']


def _pptAddTitleTextbox(slide, prs, text):
    box = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), prs.slide_width - Inches(1.0), Inches(0.8))
    tf = box.text_frame
    tf.text = text
    tf.paragraphs[0].font.size = Pt(28)
    tf.paragraphs[0].font.bold = True
    return box


def _pptStyleTable(table, nRows, nCols):
    for j in range(nCols):
        cell = table.cell(0, j)
        cell.fill.solid()
        cell.fill.fore_color.rgb = RGBColor(0x30, 0x30, 0x30)
        for para in cell.text_frame.paragraphs:
            para.font.size = Pt(12)
            para.font.bold = True
            para.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    for i in range(1, nRows):
        for j in range(nCols):
            for para in table.cell(i, j).text_frame.paragraphs:
                para.font.size = Pt(11)


def buildPptTitleSlide(prs, job, caseArray):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    box = slide.shapes.add_textbox(Inches(1.0), Inches(2.5), prs.slide_width - Inches(2.0), Inches(2.0))
    tf = box.text_frame
    tf.text = '%s - Post Processing Report' % (job)
    tf.paragraphs[0].font.size = Pt(40)
    tf.paragraphs[0].font.bold = True
    trialsPara = tf.add_paragraph()
    trialsPara.text = ' | '.join(caseArray)
    trialsPara.font.size = Pt(20)
    datePara = tf.add_paragraph()
    datePara.text = str(date.today())
    datePara.font.size = Pt(14)
    return slide


def buildPptFieldTableSlide(prs, title, caseArray, caseSummaries, fields):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _pptAddTitleTextbox(slide, prs, title)

    availableFields = [f for f in fields if any(f in caseSummaries.get(c, {}) for c in caseArray)]
    if not availableFields:
        return slide

    nRows = len(availableFields) + 1
    nCols = len(caseArray) + 1
    tableShape = slide.shapes.add_table(nRows, nCols, Inches(0.5), Inches(1.3),
                                         prs.slide_width - Inches(1.0), Inches(0.4) * nRows)
    table = tableShape.table
    table.cell(0, 0).text = 'Metric'
    for j, caseName in enumerate(caseArray):
        table.cell(0, j + 1).text = caseName
    for i, field in enumerate(availableFields):
        table.cell(i + 1, 0).text = field
        for j, caseName in enumerate(caseArray):
            table.cell(i + 1, j + 1).text = str(caseSummaries.get(caseName, {}).get(field, 'N/A'))
    _pptStyleTable(table, nRows, nCols)
    return slide


def buildPptDeltaTableSlide(prs, caseArray, caseSummaries, fields):
    refCase = caseArray[0]
    compareCases = caseArray[1:]
    if not compareCases:
        return None

    availableFields = [f for f in fields if f in caseSummaries.get(refCase, {})]
    if not availableFields:
        return None

    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _pptAddTitleTextbox(slide, prs, 'Results Delta to %s' % (refCase))

    nRows = len(availableFields) + 1
    nCols = len(compareCases) + 1
    tableShape = slide.shapes.add_table(nRows, nCols, Inches(0.5), Inches(1.3),
                                         prs.slide_width - Inches(1.0), Inches(0.4) * nRows)
    table = tableShape.table
    table.cell(0, 0).text = 'Metric'
    for j, caseName in enumerate(compareCases):
        table.cell(0, j + 1).text = '%s - %s' % (caseName, refCase)
    for i, field in enumerate(availableFields):
        table.cell(i + 1, 0).text = field
        refValue = caseSummaries[refCase].get(field, None)
        for j, caseName in enumerate(compareCases):
            compareValue = caseSummaries.get(caseName, {}).get(field, None)
            try:
                text = '%0.3f' % (float(compareValue) - float(refValue))
            except (TypeError, ValueError):
                text = 'N/A'
            table.cell(i + 1, j + 1).text = text
    _pptStyleTable(table, nRows, nCols)
    return slide


def _fitAndCenterPicture(shape, boxLeft, boxTop, boxWidth, boxHeight):
    """Rescales a picture/movie shape (preserving its current aspect ratio) to fit within the
    given box, then centers it in that box both horizontally and vertically."""
    scale = min(boxWidth / shape.width, boxHeight / shape.height)
    shape.width = int(shape.width * scale)
    shape.height = int(shape.height * scale)
    shape.left = int(boxLeft + (boxWidth - shape.width) / 2)
    shape.top = int(boxTop + (boxHeight - shape.height) / 2)


def addPptImageSlide(prs, title, imagePath):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _pptAddTitleTextbox(slide, prs, title)
    pictureTop = Inches(1.2)
    boxWidth = prs.slide_width - Inches(2.0)
    boxHeight = prs.slide_height - pictureTop - Inches(0.3)
    picture = slide.shapes.add_picture(imagePath, left=Inches(1.0), top=pictureTop, height=boxHeight)
    _fitAndCenterPicture(picture, Inches(1.0), pictureTop, boxWidth, boxHeight)
    return slide


def addPptSideBySideImageSlide(prs, title, leftImagePath, rightImagePath):
    """Adds a slide with two images placed side by side, each scaled to fit half the slide
    width while preserving its own aspect ratio (so neither image is squished), and centered
    (both horizontally and vertically) within its half."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _pptAddTitleTextbox(slide, prs, title)
    pictureTop = Inches(1.2)
    margin = Inches(0.4)
    gap = Inches(0.3)
    availableHeight = prs.slide_height - pictureTop - Inches(0.3)
    halfWidth = (prs.slide_width - 2 * margin - gap) / 2

    for i, imagePath in enumerate((leftImagePath, rightImagePath)):
        left = margin + i * (halfWidth + gap)
        picture = slide.shapes.add_picture(imagePath, left=left, top=pictureTop, width=halfWidth)
        _fitAndCenterPicture(picture, left, pictureTop, halfWidth, availableHeight)
    return slide


#CFD render images are written by pvPost.py's saveImages() as:
#   postProcessing/images/<variable>_<imageType>/<caseName>_<variable>_<imageType>_<view>.jpeg
#(see pvPost.py saveImages()). Note the .jpeg extension (not .png) and that the iso-surface
#variable names are renamed from their vtp filenames ('isoCtp.vtp' -> 'Cpt', 'isoQ.vtp' -> 'Q').
PPT_IMAGE_GROUPS = [
    ('Geom', 'Surface', 'Geometry'),
    ('CpMean', 'Surface', 'Cp Mean'),
    ('CfMean', 'Surface', 'Cf Mean'),
    ('CpPrime2Mean', 'Surface', 'Cp RMS'),
    ('Cpt', 'isoSurface', 'Cpt = 0 Iso-Surface'),
    ('Q', 'isoSurface', 'Q-Criterion Iso-Surface'),
]
PPT_IMAGE_VIEWS = ['Front', 'FrontLeft', 'Left', 'Bottom', 'RearLeft', 'Rear']


def findPvPostImage(trialPath, variable, imageType, caseName, view):
    """Locate a pvPost.py-rendered image, matching its current dirName/fileName convention
    exactly (dirName='<variable>_<imageType>', fileName='<caseName>_<variable>_<imageType>_<view>.jpeg).
    Returns None if not found."""
    dirName = '%s_%s' % (variable, imageType)
    fileName = '%s_%s_%s_%s.jpeg' % (caseName, variable, imageType, view)
    imagePath = os.path.join(trialPath, 'postProcessing', 'images', dirName, fileName)
    return imagePath if os.path.isfile(imagePath) else None


def addPvPostImageSlides(prs, path, caseArray):
    """Add one slide per found (variable, imageType, view) render for each trial, using
    pvPost.py's current image naming. Silently skips any combination that isn't found (a
    trial may not have run every view/variable)."""
    print('\tAdding CFD render image slides...')
    anyFound = False
    for variable, imageType, label in PPT_IMAGE_GROUPS:
        for view in PPT_IMAGE_VIEWS:
            for trial in caseArray:
                #pvPost.py names its own output files using only the case directory's own
                #basename (it has no knowledge of any ride-height parent prefix), so an expanded
                #child trial like 'parentTrial/childName' must be looked up by 'childName' even
                #though its containing directory is path/parentTrial/childName.
                imagePath = findPvPostImage(os.path.join(path, trial), variable, imageType,
                                             os.path.basename(trial), view)
                if imagePath:
                    anyFound = True
                    addPptImageSlide(prs, '%s - %s - %s' % (label, trial, view), imagePath)
    if not anyFound:
        print('\tNo CFD render images found in postProcessing/images (run pvPost.py first), skipping.')


#pvPost.py's generateSlices() writes each frame of a slice sweep to:
#   postProcessing/images/<variable>_slice/<caseName>_<variable>_slice_<view>_<normal>_<position>_<counter>.jpeg
#(see pvPost.py generateSlices()/saveImages()). We stitch each (variable, normal, view) sweep for a
#trial into an .mp4 via ffmpeg (same approach as createMovies.py) so it can be embedded as a
#playable movie in the PPT deck.
#Views/normals match pvPost.py generateSlices()'s default sliceViews/normalsList pairing.
PPT_SLICE_MOVIE_VIEWS = [
    ('X', 'Front'),
    ('Y', 'Left'),
    ('Y', 'LeftForward'),
    ('Y', 'LeftBack'),
    ('Z', 'Top'),
]
PPT_SLICE_MOVIE_VARS = ['CpMean', 'UMean', 'CptMean']
PPT_SLICE_MOVIE_FPS = 2


def _extractSliceCounter(filename):
    numbers = re.findall(r'\d+', os.path.basename(filename))
    return int(numbers[-1]) if numbers else 0


def generateSliceMovie(trialPath, caseName, variable, normal, view):
    """Stitch a slice sweep's frames into an .mp4 via ffmpeg (assumes ffmpeg is on PATH).
    Returns (moviePath, posterImagePath) on success, or None if no matching frames were found
    or ffmpeg failed."""
    imagesDir = os.path.join(trialPath, 'postProcessing', 'images', '%s_slice' % (variable))
    pattern = os.path.join(imagesDir, '%s_%s_slice_%s_%s_*.jpeg' % (caseName, variable, view, normal))
    images = sorted(glob.glob(pattern), key=_extractSliceCounter)
    if not images:
        return None

    prefix = '%s_%s_slice_%s_%s' % (caseName, variable, view, normal)
    listFile = os.path.join(imagesDir, '%s_images.txt' % (prefix))
    with open(listFile, 'w') as f:
        for image in images:
            f.write("file '%s'\n" % (image))

    moviePath = os.path.join(imagesDir, '%s.mp4' % (prefix))
    #NOTE: -framerate MUST be an input option (before -i) so the concat demuxer displays every
    #listed frame for 1/framerate seconds. Putting the rate after -i (as an output -r) makes
    #ffmpeg resample/decimate the sequence against its default 25fps input assumption, silently
    #dropping most of the slice frames instead of including all of them.
    cmd = ("ffmpeg -y -framerate %s -f concat -safe 0 -i '%s' -vf scale=1920:1080 "
           "-c:v libx264 -pix_fmt yuv420p '%s' >> log.pptReport" % (PPT_SLICE_MOVIE_FPS, listFile, moviePath))
    ret = os.system(cmd)
    if ret != 0 or not os.path.isfile(moviePath):
        print('\tWARNING! ffmpeg failed to build slice movie for %s (see log.pptReport)' % (prefix))
        return None
    return moviePath, images[0]


def addPptMovieSlide(prs, title, moviePath, posterImagePath):
    """Embeds a movie sized to preserve the poster image's aspect ratio (add_movie itself has no
    concept of natural size, unlike add_picture, so the poster frame's real pixel dimensions are
    used to compute it), then centers it within the slide's picture box."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    _pptAddTitleTextbox(slide, prs, title)
    pictureTop = Inches(1.2)
    boxWidth = prs.slide_width - Inches(2.0)
    boxHeight = prs.slide_height - pictureTop - Inches(0.3)

    posterHeightPx, posterWidthPx = plt.imread(posterImagePath).shape[:2]
    aspect = posterHeightPx / posterWidthPx
    movieWidth = boxWidth
    movieHeight = int(movieWidth * aspect)
    movie = slide.shapes.add_movie(moviePath, left=Inches(1.0), top=pictureTop,
                                    width=movieWidth, height=movieHeight,
                                    poster_frame_image=posterImagePath)
    _fitAndCenterPicture(movie, Inches(1.0), pictureTop, boxWidth, boxHeight)
    return slide


def addSliceMovieSlides(prs, path, caseArray):
    """Build and embed one slice-sweep movie slide per (variable, view) for each trial, covering
    every default slice view for that variable, so they play back directly in the slideshow.
    Requires ffmpeg on PATH. Silently skips any variable/view/trial combination with no matching
    slice frames."""
    print('\tGenerating slice movies (requires ffmpeg on PATH)...')
    anyFound = False
    for variable in PPT_SLICE_MOVIE_VARS:
        for normal, view in PPT_SLICE_MOVIE_VIEWS:
            for trial in caseArray:
                #see addPvPostImageSlides: pvPost.py's own slice frame filenames use only the
                #case directory's basename, not any ride-height parent prefix.
                result = generateSliceMovie(os.path.join(path, trial), os.path.basename(trial),
                                             variable, normal, view)
                if result:
                    anyFound = True
                    moviePath, posterImagePath = result
                    title = '%s - %s Slices - %s' % (variable, view, trial)
                    addPptMovieSlide(prs, title, moviePath, posterImagePath)
    if not anyFound:
        print('\tNo slice images found for movie generation, skipping.')


def buildForceHistoryImages(args, caseArray):
    """Generate (or refresh) the forceHistory_<var> plots for caseArray, reusing the same
    setCasePaths/getCaseData/makePandasArrays/plotData pipeline as --forces. Returns the
    directory the images were saved to."""
    localArgs = argparse.Namespace(**vars(args))
    localArgs.trial = caseArray
    casePathDict, localCaseLoc = setCasePaths(caseArray, casePath)
    casePathDict = getCaseData(casePathDict)
    casePathDict = makePandasArrays(localArgs, casePathDict)
    plotData(localArgs, localCaseLoc, casePathDict)
    if localCaseLoc.lower() == 'outtrial':
        return os.path.join(path, caseArray[0])
    return casePath


#--- Binned force plot (ported from legacy forceBinPlot.py, x-axis units per binPlotForces_v2_0.py) ---
#pvPost.py's Left view camera looks from +Y toward the origin with viewup=+Z (see
#default/defaultViews), and the Front/Rear cameras confirm +X is the front of the vehicle.
#cross(forward=-Y, up=+Z) = -X, i.e. screen-right = -X, so the vehicle's front (+X) renders on
#the LEFT edge of the Left-view image and the rear on the right edge. Combined with bin index
#0 = front (per user), the leftmost opaque (body) pixel column maps to the smallest x co-ord
#and the rightmost maps to the largest.
PPT_BIN_PLOT_WHITE_THRESH = 235  # 0-255; pixels with all channels >= this are background
PPT_BIN_PLOT_BLACK_THRESH = 60   # 0-255; pixels with all channels <= this are title text
PPT_BIN_PLOT_IMG_ALPHA_SINGLE = 0.5
PPT_BIN_PLOT_IMG_ALPHA_MULTI = 0.35


def parseBinXCoords(coeffFile):
    """Parses the '# x co-ords :' header line out of a binForceCoeffs .dat file (as
    binPlotForces_v2_0.py's importBinData() does) to get each bin's real x position in meters.
    Returns a numpy array, or None if the header line wasn't found."""
    with open(coeffFile, 'r') as f:
        for line in f:
            if 'x co-ords' in line:
                xCoords = line.replace('#', '').replace(':', '').replace('x', '').replace('co-ords', '')
                return np.array([float(v) for v in xCoords.split()])
    return None


def loadBinForceCoeffs(fullCaseSetupDict, path, case):
    """Reads postProcessing/binForceCoeffs/<latestTime>/forceCoeffBin.dat for a case, replicating
    the original forceBinPlot.py's per-bin coefficient/force extraction, plus the real per-bin x
    co-ordinates (in meters) parsed the same way as binPlotForces_v2_0.py. Returns a dict with
    xCoeffs/yCoeffs/zCoeffs/xForce/yForce/zForce/xCoords (100-element arrays), or None if no bin
    force data (or no x co-ords header) is available for this case."""
    scale = 2 if 'half' in case else 1
    surfaceLoc = glob.glob(os.path.join(path, case, 'postProcessing', 'binForceCoeffs', '*'))
    if not surfaceLoc:
        return None
    lastTime = os.path.basename(surfaceLoc[0])
    coeffFile = os.path.join(path, case, 'postProcessing', 'binForceCoeffs', lastTime, 'forceCoeffBin.dat')
    if not os.path.isfile(coeffFile):
        return None

    xCoords = parseBinXCoords(coeffFile)
    if xCoords is None:
        print('\tWARNING! Could not find "x co-ords" header in %s, skipping.' % (coeffFile))
        return None

    inletMag = bcParser(fullCaseSetupDict, path, case)[0]
    forceMultiplier = (float(inletMag) ** 2) * 0.5

    coeffs = np.loadtxt(coeffFile, dtype='float', comments='#', delimiter=None, skiprows=10)
    coeffs = coeffs[1:]
    n = 100
    xCoeffs = np.zeros(n)
    yCoeffs = np.zeros(n)
    zCoeffs = np.zeros(n)
    for i in range(n):
        start = 9 * i
        binVals = coeffs[start:start + 3] * scale
        xCoeffs[i], yCoeffs[i], zCoeffs[i] = binVals[0], binVals[1], binVals[2] * -1

    return {
        'xCoeffs': xCoeffs, 'yCoeffs': yCoeffs, 'zCoeffs': zCoeffs,
        'xForce': xCoeffs * forceMultiplier * 1.225,
        'yForce': yCoeffs * forceMultiplier * 1.225,
        'zForce': zCoeffs * forceMultiplier * 1.225,
        'xCoords': xCoords,
    }


def loadVehicleSideImageRGBA(imagePath, whiteThresh=PPT_BIN_PLOT_WHITE_THRESH, blackThresh=PPT_BIN_PLOT_BLACK_THRESH):
    """Loads a pvPost.py Geom_Surface Left-view render and makes the white background and black
    title text transparent, leaving only the grey vehicle body opaque. Anti-aliased edge pixels
    (e.g. faint outline remnants of the removed title text, or the body's own anti-aliased
    silhouette edge) can be neither near-white nor near-black and so would otherwise survive the
    threshold and get included in the body's bounding box, leaving a faint gap between the
    cropped image's edge and the actual vehicle silhouette. To avoid that, only the single
    largest connected blob of kept pixels (the vehicle body itself) is retained; any other
    disconnected leftover speckle (e.g. text-edge anti-aliasing) is discarded. Returns the RGBA
    array cropped to the columns spanning the vehicle body (so column 0 = vehicle front, last
    column = vehicle rear), or None if the image has no visible body left after masking."""
    img = plt.imread(imagePath)
    if img.dtype == np.uint8:
        img = img.astype(float) / 255.0
    rgb = img[:, :, :3]
    isWhite = np.all(rgb >= whiteThresh / 255.0, axis=2)
    isBlack = np.all(rgb <= blackThresh / 255.0, axis=2)
    keep = ~(isWhite | isBlack)

    labeled, numBlobs = ndimage.label(keep)
    if numBlobs == 0:
        return None
    blobSizes = ndimage.sum(keep, labeled, index=range(1, numBlobs + 1))
    largestBlob = int(np.argmax(blobSizes)) + 1
    keep = labeled == largestBlob

    bodyCols = np.nonzero(keep.any(axis=0))[0]
    if bodyCols.size == 0:
        return None
    firstCol, lastCol = int(bodyCols[0]), int(bodyCols[-1])

    rgba = np.zeros((img.shape[0], img.shape[1], 4), dtype=float)
    rgba[:, :, :3] = rgb
    rgba[:, :, 3] = np.where(keep, 1.0, 0.0)
    return rgba[:, firstCol:lastCol + 1, :]


def _plotBinForceComponent(ax, caseBinData, coeffKey, label, colors):
    """Plots one binned coefficient component (e.g. 'zCoeffs'/Cl or 'xCoeffs'/Cd) for every case
    onto ax using each case's real per-bin x co-ordinates (in meters, parsed from the
    binForceCoeffs header, per binPlotForces_v2_0.py's importBinData()). Overlays each case's
    vehicle silhouette scaled to that same case's x co-ord extent (as binPlotForces_v2_0.py
    does), and sets the axes box's physical aspect ratio to match the (undistorted) vehicle
    image so the car isn't squished/stretched. The shared x-axis limit spans the min/max x
    co-ord across all cases."""
    maxCoeff = max(float(np.max(np.abs(d[coeffKey]))) for d in caseBinData.values())
    yTop = maxCoeff * 1.2 if maxCoeff > 0 else 1.0
    imgAlpha = PPT_BIN_PLOT_IMG_ALPHA_SINGLE if len(caseBinData) == 1 else PPT_BIN_PLOT_IMG_ALPHA_MULTI

    xMin = min(float(np.min(d['xCoords'])) for d in caseBinData.values())
    xMax = max(float(np.max(d['xCoords'])) for d in caseBinData.values())

    imgAspect = None
    for i, (case, binData) in enumerate(caseBinData.items()):
        color = colors[i % len(colors)]
        #case may be a compound 'parentTrial/childName' identifier (expanded ride-height
        #child); the parent prefix is redundant in the legend since it's shared across all
        #plotted cases, so just show the case's own basename.
        ax.plot(binData['xCoords'], binData[coeffKey], '-', linewidth=1, color=color, label=os.path.basename(case))

        imagePath = findPvPostImage(os.path.join(path, case), 'Geom', 'Surface', os.path.basename(case), 'Left')
        if not imagePath:
            print('\tWARNING! No left-view geometry image found for %s, skipping image overlay.' % (case))
            continue
        rgba = loadVehicleSideImageRGBA(imagePath)
        if rgba is None:
            print('\tWARNING! Could not isolate vehicle body in %s, skipping image overlay.' % (imagePath))
            continue
        #each case's image is scaled to that case's own x co-ord extent (matching
        #binPlotForces_v2_0.py), so cases with slightly different bin extents each still line up
        #their own vehicle silhouette against their own data.
        caseXMin = float(np.min(binData['xCoords']))
        caseXMax = float(np.max(binData['xCoords']))
        #vehicle image is already cropped to exactly the body extent, so mapping its full width
        #to [caseXMin, caseXMax] and its full height to [0, yTop] means the axes box's physical
        #aspect ratio (set below) purely determines whether it's distorted, independent of yTop.
        if imgAspect is None:
            imgAspect = rgba.shape[0] / rgba.shape[1]
        ax.imshow(rgba, extent=[caseXMin, caseXMax, 0, yTop], aspect='auto', alpha=imgAlpha, zorder=0)

    ax.set_xlim(xMin, xMax)
    ax.set_ylim(0, yTop)
    if imgAspect is not None:
        ax.set_box_aspect(imgAspect)
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Coefficient')
    ax.set_title(label)
    ax.legend()


def _plotBinForceDeltaComponent(ax, baseCase, caseBinData, coeffKey, label, colors):
    """Plots, for every case except baseCase, the delta of one binned coefficient component
    (compareCase - baseCase) against distance (m). Since each case has its own real bin x
    co-ords, each compare case's coefficient values are linearly interpolated (np.interp) onto
    baseCase's x co-ords before differencing, so the delta is defined at baseCase's bin
    positions. Overlays baseCase's own vehicle silhouette (as the common reference geometry)
    scaled to baseCase's x co-ord extent, with the axes box's aspect ratio matched to the image.
    """
    baseBinData = caseBinData[baseCase]
    baseX = baseBinData['xCoords']
    baseVals = baseBinData[coeffKey]
    compareCases = [c for c in caseBinData if c != baseCase]

    deltas = {}
    for case in compareCases:
        binData = caseBinData[case]
        interpVals = np.interp(baseX, binData['xCoords'], binData[coeffKey])
        deltas[case] = interpVals - baseVals

    if not deltas:
        ax.set_title(label)
        return

    maxDelta = max(float(np.max(np.abs(d))) for d in deltas.values())
    yLim = maxDelta * 1.2 if maxDelta > 0 else 1.0
    imgAlpha = PPT_BIN_PLOT_IMG_ALPHA_SINGLE if len(deltas) == 1 else PPT_BIN_PLOT_IMG_ALPHA_MULTI

    #colors offset by 1 so a compare case keeps the same color it used in the absolute plot
    #(where baseCase took colors[0]).
    for i, case in enumerate(compareCases):
        ax.plot(baseX, deltas[case], '-', linewidth=1, color=colors[(i + 1) % len(colors)],
                label='%s - %s' % (os.path.basename(case), os.path.basename(baseCase)))

    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')

    imgAspect = None
    imagePath = findPvPostImage(os.path.join(path, baseCase), 'Geom', 'Surface', os.path.basename(baseCase), 'Left')
    if imagePath:
        rgba = loadVehicleSideImageRGBA(imagePath)
        if rgba is not None:
            imgAspect = rgba.shape[0] / rgba.shape[1]
            baseXMin, baseXMax = float(np.min(baseX)), float(np.max(baseX))
            ax.imshow(rgba, extent=[baseXMin, baseXMax, -yLim, yLim], aspect='auto', alpha=imgAlpha, zorder=0)
        else:
            print('\tWARNING! Could not isolate vehicle body in %s, skipping image overlay.' % (imagePath))
    else:
        print('\tWARNING! No left-view geometry image found for %s, skipping image overlay.' % (baseCase))

    ax.set_xlim(float(np.min(baseX)), float(np.max(baseX)))
    ax.set_ylim(-yLim, yLim)
    if imgAspect is not None:
        ax.set_box_aspect(imgAspect)
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Coefficient Delta')
    ax.set_title(label)
    ax.legend()


def buildBinForcePlots(path, caseArray, outputDir):
    """Builds separate binned Cl (downforce) and Cd (drag) coefficient-vs-distance(m) figures
    (one PNG each) for caseArray, each overlaying every case's own left-view vehicle silhouette
    (background/text removed) scaled to that case's real bin x co-ordinate extent (parsed from
    binForceCoeffs, per binPlotForces_v2_0.py), with the axes box's aspect ratio matched to the
    vehicle image so it isn't squished/stretched. When multiple cases are given, each case's plot
    line and image use a distinct color and the images are drawn at reduced opacity so they can
    be visually stacked/compared. Also builds Cl/Cd delta-vs-baseline figures (baseline =
    caseArray[0], matching the convention used by the Results Delta table) when more than one
    case has bin data, each compare case's coefficients interpolated onto the baseline's x
    co-ords before differencing. Also writes a trial<case>_binForces.csv per case, like the
    legacy script. Returns a dict with 'clPlotPath'/'cdPlotPath' (always present) and
    'clDeltaPlotPath'/'cdDeltaPlotPath' (present only if there were >= 2 cases with bin data), or
    None if no case had bin force data."""
    caseBinData = {}
    for case in caseArray:
        caseSetupPath = os.path.join(path, case, 'fullCaseSetupDict')
        if not os.path.isfile(caseSetupPath):
            print('\tWARNING! %s missing fullCaseSetupDict, skipping bin plot.' % (case))
            continue
        fullCaseSetupDict = configparser.ConfigParser()
        fullCaseSetupDict.optionxform = str
        fullCaseSetupDict.read_file(open(caseSetupPath))

        binData = loadBinForceCoeffs(fullCaseSetupDict, path, case)
        if binData is None:
            print('\tWARNING! %s has no binForceCoeffs data, skipping.' % (case))
            continue
        caseBinData[case] = binData

    if not caseBinData:
        return None

    for case, binData in caseBinData.items():
        allForces = np.vstack((binData['xCoords'], binData['xCoeffs'], binData['yCoeffs'],
                                binData['zCoeffs'], binData['xForce'], binData['yForce'],
                                binData['zForce'])).T
        #case may be a compound 'parentTrial/childName' identifier (expanded ride-height child);
        #sanitize to a flat filename-safe tag since this is a NEW output file, not a lookup of
        #something pvPost.py already wrote.
        caseTag = case.replace('/', '_')
        np.savetxt(os.path.join(outputDir, 'trial%s_binForces.csv' % (caseTag)), allForces, delimiter=',',
                   header='xCoords,xCoeffs,yCoeffs,zCoeffs,xForce,yForce,zForce')

    colors = plt.cm.tab10.colors
    reportName = buildCaseTag(list(caseBinData.keys()))

    figCl, axCl = plt.subplots()
    _plotBinForceComponent(axCl, caseBinData, 'zCoeffs', 'Binned Cl (Downforce)', colors)
    figCl.tight_layout()
    clPlotPath = os.path.join(outputDir, '%s_binnedCl.png' % (reportName))
    figCl.savefig(clPlotPath, dpi=300)
    plt.close(figCl)

    figCd, axCd = plt.subplots()
    _plotBinForceComponent(axCd, caseBinData, 'xCoeffs', 'Binned Cd (Drag)', colors)
    figCd.tight_layout()
    cdPlotPath = os.path.join(outputDir, '%s_binnedCd.png' % (reportName))
    figCd.savefig(cdPlotPath, dpi=300)
    plt.close(figCd)

    result = {'clPlotPath': clPlotPath, 'cdPlotPath': cdPlotPath}

    if len(caseBinData) > 1:
        baseCase = next(c for c in caseArray if c in caseBinData)

        figClDelta, axClDelta = plt.subplots()
        _plotBinForceDeltaComponent(axClDelta, baseCase, caseBinData, 'zCoeffs',
                                     'Binned Cl Delta vs %s' % (baseCase), colors)
        figClDelta.tight_layout()
        clDeltaPlotPath = os.path.join(outputDir, '%s_binnedClDelta.png' % (reportName))
        figClDelta.savefig(clDeltaPlotPath, dpi=300)
        plt.close(figClDelta)

        figCdDelta, axCdDelta = plt.subplots()
        _plotBinForceDeltaComponent(axCdDelta, baseCase, caseBinData, 'xCoeffs',
                                     'Binned Cd Delta vs %s' % (baseCase), colors)
        figCdDelta.tight_layout()
        cdDeltaPlotPath = os.path.join(outputDir, '%s_binnedCdDelta.png' % (reportName))
        figCdDelta.savefig(cdDeltaPlotPath, dpi=300)
        plt.close(figCdDelta)

        result['clDeltaPlotPath'] = clDeltaPlotPath
        result['cdDeltaPlotPath'] = cdDeltaPlotPath

    return result


def generate_ppt_report(args):
    print('\n\tGenerating PowerPoint report...')

    #ride-height mapping parent cases (any --trial entry with child_# dirs, whether that's the
    #cwd default or an explicitly-passed -t) are expanded into their individual child cases so
    #each grid point gets full CFD render/bin-plot/etc. slides, same as manually listing every
    #child case with -t. Child case directories live NESTED inside the parent
    #(path/parentTrial/childName), not directly under path, so each child is tracked as the
    #compound relative path 'parentTrial/childName' -- every downstream os.path.join(path, trial)
    #call then resolves to the correct nested location. The parent name itself is kept aside
    #(not added to caseArray) so its own summary.csv (the parent-level average, written by
    #--summary) can still drive the "Ride Height Map Averages" table + sensitivity sweep plots
    #below.
    rideHeightParents = []
    expandedTrials = []
    for trial in args.trial:
        children = discoverRideHeightChildCases(os.path.join(path, trial), trial)
        if children:
            rideHeightParents.append(trial)
            print('\tDetected ride-height mapping case %s, expanding into %d child case(s).' %
                  (trial, len(children)))
            expandedTrials.extend(os.path.join(trial, child) for child in children)
        else:
            expandedTrials.append(trial)

    caseSummaries = {}
    for trial in expandedTrials:
        summaryPath = os.path.join(path, trial, 'summary.csv')
        if not os.path.isfile(summaryPath):
            print('\tWARNING! %s is missing summary.csv (run --summary for it first), skipping.' % (trial))
            continue
        summaryDict = readChildSummaryCsv(summaryPath)
        if not summaryDict:
            print('\tWARNING! Unable to read summary.csv for %s, skipping.' % (trial))
            continue
        caseSummaries[trial] = summaryDict

    caseArray = list(caseSummaries.keys())
    if len(caseArray) < 1:
        sys.exit('ERROR! No valid case summaries found, cannot build PPT report. Run --summary for each trial first.')

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    buildPptTitleSlide(prs, job, caseArray)
    buildPptFieldTableSlide(prs, 'Trial Setup and Boundary Conditions', caseArray, caseSummaries, PPT_INFO_FIELDS)
    buildPptFieldTableSlide(prs, 'Results', caseArray, caseSummaries, PPT_RESULTS_FIELDS)
    if len(caseArray) > 1:
        buildPptDeltaTableSlide(prs, caseArray, caseSummaries, PPT_RESULTS_FIELDS)

    print('\tGenerating force history plots...')
    try:
        forceImageDir = buildForceHistoryImages(args, caseArray)
        #matches forceConvergencePlot.py's plotData(), which collapses compound
        #'parentTrial/childName' case identifiers (expanded ride-height children) into a
        #short '<parent>_rhmap' tag before joining them into the saved filename.
        trialTag = buildCaseTag(caseArray)
        for var in args.plotData:
            imagePath = os.path.join(forceImageDir, '%s_forceHistory_%s.%s' %
                                      (trialTag, var, args.saveFormat))
            if os.path.isfile(imagePath):
                addPptImageSlide(prs, 'Force History - %s' % (var), imagePath)
            else:
                print('\tWARNING! Could not find %s, skipping slide.' % (imagePath))
    except Exception as e:
        print('\tWARNING! Unable to generate force history plots: %s' % (e))

    binPlots = buildBinForcePlots(path, caseArray, casePath)
    if binPlots:
        addPptSideBySideImageSlide(prs, 'Binned Forces', binPlots['clPlotPath'], binPlots['cdPlotPath'])
        if 'clDeltaPlotPath' in binPlots:
            addPptSideBySideImageSlide(prs, 'Binned Force Deltas', binPlots['clDeltaPlotPath'],
                                        binPlots['cdDeltaPlotPath'])
    else:
        print('\tNo binForceCoeffs data found for any trial, skipping binned force plots.')

    addPvPostImageSlides(prs, path, caseArray)
    if args.addMovies:
        addSliceMovieSlides(prs, path, caseArray)
    else:
        print('\tSkipping slice movie generation (pass --addMovies to include).')

    

    #ride height: any ride-height mapping parent detected above (see expansion at the top of
    #this function) gets a map-averages table (each parent's own already-averaged summary.csv)
    #plus per-point sweep comparison plots (reusing plotRideHeightSensitivity's cross-case sweep
    #matching)
    if rideHeightParents:
        parentSummaries = {}
        for parentTrial in rideHeightParents:
            summaryPath = os.path.join(path, parentTrial, 'summary.csv')
            if not os.path.isfile(summaryPath):
                print('\tWARNING! %s missing its own (parent) summary.csv (run --summary on it), '
                      'skipping from Ride Height Map Averages.' % (parentTrial))
                continue
            summaryDict = readChildSummaryCsv(summaryPath)
            if summaryDict:
                parentSummaries[parentTrial] = summaryDict

        if parentSummaries:
            rideHeightTrials = list(parentSummaries.keys())
            print('\tDetected ride-height mapping case(s): %s' % (', '.join(rideHeightTrials)))
            buildPptFieldTableSlide(prs, 'Ride Height Map Averages', rideHeightTrials, parentSummaries,
                                     PPT_RESULTS_FIELDS)

            primaryPath = os.path.join(path, rideHeightTrials[0])
            comparePaths = [os.path.join(path, t) for t in rideHeightTrials[1:]]
            sensitivityDir = os.path.join(primaryPath, 'postProcessing', 'sensitivityPlots')
            plotRideHeightSensitivity(primaryPath, includeSideForce=args.includeSideForce,
                                       compareCasePaths=comparePaths)
            for plotFile in sorted(glob.glob(os.path.join(sensitivityDir, '*.png'))):
                title = 'Ride Height - %s' % (os.path.splitext(os.path.basename(plotFile))[0])
                addPptImageSlide(prs, title, plotFile)

    reportName = buildCaseTag(caseArray)
    outputPath = os.path.join(casePath, '%s_report_%s.pptx' % (reportName, date.today()))
    prs.save(outputPath)
    print('\tSaved PPT report to %s' % (outputPath))


main()
