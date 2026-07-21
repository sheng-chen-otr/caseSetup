import os
import sys
import numpy as np
import subprocess as sp
import math
import gzip
import re
import json
import csv
import multiprocessing
import pandas as pd
import configparser
import struct
import shutil
import subprocess
from datetime import datetime
from copy import deepcopy

updateCaseSetupFlag = False



def transformForceVector(fullCaseSetupDict,forceVec):
    domain_pitch = math.radians(float(fullCaseSetupDict['BC_SETUP']['DOMAIN_PITCH'][0]))
    domain_roll = math.radians(float(fullCaseSetupDict['BC_SETUP']['DOMAIN_ROLL'][0]))
    forceVec = forceVec.split(' ')
    forceVec = x_rotation(np.array(forceVec).astype(float),domain_roll)
    forceVec = y_rotation(np.array(forceVec).astype(float),domain_pitch)
    return " ".join(str(round(x,6)) for x in forceVec)



def velVector(inletMag,yaw,pitch):
    initVel = [float(inletMag), 0, 0]
    pitch = math.radians(pitch)
    yaw = math.radians(yaw)
    initDragVec = [1,0,0]
    initLiftVec = [0,0,1]
    #transform inlet vector for pitch
    pitchTransVel = y_rotation(initVel,pitch)
    #transform inlet vector for yaw
    yawTransVel = z_rotation(pitchTransVel,yaw)
    
    #" ".join(str(x) for x in xs)
    #transform drag and lift vectors for pitch
    dragVec = y_rotation(initDragVec,pitch)
    liftVec = y_rotation(initLiftVec,pitch)
    return " ".join(str(x) for x in yawTransVel), " ".join(str(x) for x in dragVec)," ".join(str(x) for x in liftVec)


def _truthy(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ['true', '1', 'yes', 'on']


def _parse_keyword_list(raw, toLower=True):
    if raw is None:
        return []
    if isinstance(raw, list):
        tokens = []
        for item in raw:
            if isinstance(item, str):
                tokens.extend(item.replace(',', ' ').split())
            else:
                tokens.append(str(item))
    else:
        tokens = str(raw).replace(',', ' ').split()
    cleaned = [str(t).strip() for t in tokens if str(t).strip() != '']
    if toLower:
        return [t.lower() for t in cleaned]
    return cleaned


def _resolve_existing_path(pathCandidate):
    if os.path.exists(pathCandidate):
        return pathCandidate
    execDir = os.path.dirname(os.path.realpath(__file__))
    fallback = os.path.join(execDir, pathCandidate)
    if os.path.exists(fallback):
        return fallback
    fallbackRh = os.path.join(execDir, 'rideHeightUtils', pathCandidate)
    if os.path.exists(fallbackRh):
        return fallbackRh
    return None


def _scale_suspension_hardpoints(obj, scale):
    if isinstance(obj, dict):
        return {k: _scale_suspension_hardpoints(v, scale) for k, v in obj.items()}
    if isinstance(obj, list):
        if len(obj) == 3 and all(isinstance(x, (int, float)) for x in obj):
            return [float(x) * scale for x in obj]
        return [_scale_suspension_hardpoints(v, scale) for v in obj]
    return obj


def _parse_cfg_vec3(raw, keyName):
    vals = str(raw).replace(',', ' ').split()
    if len(vals) != 3:
        raise ValueError('%s must contain exactly 3 numeric values' % (keyName))
    return [float(vals[0]), float(vals[1]), float(vals[2])]


def _load_suspension_hardpoints_cfg(hardpointPath):
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read_file(open(hardpointPath))

    cornerSections = {
        'fl': 'SUSP_FL',
        'fr': 'SUSP_FR',
        'rl': 'SUSP_RL',
        'rr': 'SUSP_RR',
    }

    keyMap = {
        'wheel_center_static': 'WHEEL_CENTER_STATIC',
        'uca_f_inner': 'UCA_F_INNER',
        'uca_r_inner': 'UCA_R_INNER',
        'lca_f_inner': 'LCA_F_INNER',
        'lca_r_inner': 'LCA_R_INNER',
        'tie_inner': 'TIE_INNER',
        'uca_outer_static': 'UCA_OUTER_STATIC',
        'lca_outer_static': 'LCA_OUTER_STATIC',
        'tie_outer_static': 'TIE_OUTER_STATIC',
        'pushrod_outer_static': 'PUSHROD_OUTER_STATIC',
        'wheel_axis_local': 'WHEEL_AXIS_LOCAL',
        'wheel_forward_local': 'WHEEL_FORWARD_LOCAL',
        'rack_axis': 'RACK_AXIS',
    }

    rockerKeyMap = {
        'pivot': 'ROCKER_PIVOT',
        'axis': 'ROCKER_AXIS',
        'pushrod_joint_ref': 'ROCKER_PUSHROD_JOINT_REF',
        'damper_joint_ref': 'ROCKER_DAMPER_JOINT_REF',
        'damper_chassis': 'ROCKER_DAMPER_CHASSIS',
    }

    corners = {}
    cornerPidKeywords = {'fl': [], 'fr': [], 'rl': [], 'rr': []}

    for corner, sec in cornerSections.items():
        if not cfg.has_section(sec):
            continue

        c = {}
        for outKey, cfgKey in keyMap.items():
            if cfg.has_option(sec, cfgKey):
                c[outKey] = _parse_cfg_vec3(cfg.get(sec, cfgKey), '%s.%s' % (sec, cfgKey))

        rocker = {}
        for outKey, cfgKey in rockerKeyMap.items():
            if not cfg.has_option(sec, cfgKey):
                raise ValueError('Missing required key %s in section %s' % (cfgKey, sec))
            rocker[outKey] = _parse_cfg_vec3(cfg.get(sec, cfgKey), '%s.%s' % (sec, cfgKey))
        c['rocker'] = rocker

        for required in ['wheel_center_static','uca_f_inner','uca_r_inner','lca_f_inner','lca_r_inner','tie_inner','uca_outer_static','lca_outer_static','tie_outer_static','pushrod_outer_static']:
            if required not in c:
                raise ValueError('Missing required key for %s: %s' % (sec, required))

        # rack travel direction for steering. optional, defaults to vehicle-lateral (y)
        if 'rack_axis' not in c:
            c['rack_axis'] = [0.0, 1.0, 0.0]

        corners[corner] = c

        for key, val in cfg.items(sec):
            keyUp = key.upper()
            if 'PID' not in keyUp:
                continue
            v = str(val).strip()
            if v != '':
                cornerPidKeywords[corner].append(v)

    # optional component sections: merge PID keywords by CORNER field
    for sec in cfg.sections():
        if sec in cornerSections.values():
            continue
        if not cfg.has_option(sec, 'CORNER'):
            continue
        c = cfg.get(sec, 'CORNER').strip().lower()
        if c not in cornerPidKeywords:
            continue
        for key, val in cfg.items(sec):
            keyUp = key.upper()
            if 'PID' not in keyUp:
                continue
            v = str(val).strip()
            if v != '':
                cornerPidKeywords[c].append(v)

    return corners, cornerPidKeywords


def _load_suspension_kinematics_setup(fullCaseSetupDict):
    # read from [RIDE_HEIGHT_SETUP] so no custom section handling is needed
    section = None
    sectionName = None
    if 'RIDE_HEIGHT_SETUP' in fullCaseSetupDict and 'USE_KINEMATIC_SOLVER' in fullCaseSetupDict['RIDE_HEIGHT_SETUP']:
        section = fullCaseSetupDict['RIDE_HEIGHT_SETUP']
        sectionName = 'RIDE_HEIGHT_SETUP'

    if section is None:
        return None

    if not _truthy(section.get('USE_KINEMATIC_SOLVER', ['false'])[0]):
        return None

    hardpointFile = section.get('HARDPOINT_FILE', [''])[0]
    if hardpointFile == '':
        hardpointFile = section.get('SUSPENSION_HARDPOINT_FILE', [''])[0]
    if hardpointFile == '':
        print('\t\tWARNING! [%s] enabled but HARDPOINT_FILE is empty. Skipping kinematic solver.' % (sectionName))
        return None

    hardpointPath = _resolve_existing_path(hardpointFile)
    if hardpointPath is None:
        print('\t\tWARNING! Suspension hardpoint file cannot be found: %s. Skipping kinematic solver.' % (hardpointFile))
        return None

    if hardpointPath.lower().endswith('.cfg') or hardpointPath.lower().endswith('.ini'):
        try:
            corners, cornerPidKeywords = _load_suspension_hardpoints_cfg(hardpointPath)
        except Exception as e:
            print('\t\tWARNING! Invalid suspension CFG format: %s. Skipping kinematic solver.' % (e))
            return None
    else:
        print('\t\tWARNING! Unsupported hardpoint file type (CFG required): %s. Skipping kinematic solver.' % (hardpointPath))
        return None

    if not isinstance(corners, dict) or len(corners) == 0:
        print('\t\tWARNING! Suspension hardpoint file has invalid format. Expected corner dictionary. Skipping kinematic solver.')
        return None

    scaleVal = section.get('HARDPOINT_SCALE', [''])[0]
    if scaleVal == '':
        scaleVal = section.get('SUSPENSION_HARDPOINT_SCALE', ['1.0'])[0]
    scale = float(scaleVal)
    if abs(scale - 1.0) > 1e-15:
        corners = _scale_suspension_hardpoints(corners, scale)

    return {
        'corners': corners,
        'corner_pid_keywords': cornerPidKeywords,
        'source': hardpointPath,
    }


def _append_suspension_kinematics(rideHeights, fullCaseSetupDict):
    setup = _load_suspension_kinematics_setup(fullCaseSetupDict)
    if setup is None:
        return rideHeights

    try:
        from rideHeightUtils.doubleWishbonePushrodKinematics import DoubleWishbonePushrodSolver
    except Exception as e:
        print('\t\tWARNING! Could not import DoubleWishbonePushrodSolver, skipping kinematics: %s' % (e))
        return rideHeights

    print('\tApplying suspension kinematic solver using hardpoints: %s' % (setup['source']))

    cornerToTravelCol = {
        'fl': 'wheel_fl',
        'fr': 'wheel_fr',
        'rl': 'wheel_rl',
        'rr': 'wheel_rr',
    }

    nPoints = len(rideHeights)

    # optional front steer (FL road-wheel angle, deg). FL/FR share a rack so the FL
    # rack travel is reused for FR, letting FR angle (Ackermann) emerge
    steerCol = None
    for cand in ('steer', 'steer_deg', 'steer_angle'):
        if cand in rideHeights.columns:
            steerCol = cand
            break
    if steerCol is not None:
        frontSteer = [float(x) for x in rideHeights[steerCol].values]
        hasSteer = any(abs(s) > 1e-12 for s in frontSteer)
    else:
        frontSteer = [0.0] * nPoints
        hasSteer = False

    referenceCorner = 'fl'

    def _solve_corner_series(solver, travelList, steerList, steerMode):
        out = []
        state = None
        for dz, st in zip(travelList, steerList):
            solved = solver.solve_for_wheel_travel(
                dz, steer=st, steer_mode=steerMode, initial_state=state
            )
            out.append(solved)
            state = solved['state']
        return out

    frontRackTravel = None  # rack travel per point determined from the reference wheel

    for corner in ('fl', 'fr', 'rl', 'rr'):
        if corner not in setup['corners']:
            continue

        try:
            solver = DoubleWishbonePushrodSolver(setup['corners'][corner])
        except Exception as e:
            print('\t\tWARNING! Failed to initialize kinematic solver for corner %s: %s' % (corner, e))
            continue

        travel = [float(x) for x in rideHeights[cornerToTravelCol[corner]].values]

        try:
            if corner == referenceCorner and hasSteer:
                solved = _solve_corner_series(solver, travel, frontSteer, 'angle')
                frontRackTravel = [row['rack_travel'] for row in solved]
            elif corner == 'fr' and hasSteer and frontRackTravel is not None:
                solved = _solve_corner_series(solver, travel, frontRackTravel, 'rack')
            else:
                solved = _solve_corner_series(solver, travel, [0.0] * nPoints, 'none')
        except Exception as e:
            print('\t\tWARNING! Failed to solve kinematics for corner %s: %s' % (corner, e))
            continue

        rideHeights['%s_camber_deg' % (corner)] = [row['camber_deg'] for row in solved]
        rideHeights['%s_toe_deg' % (corner)] = [row['toe_deg'] for row in solved]
        rideHeights['%s_damper_delta' % (corner)] = [row['damper_delta'] for row in solved]
        rideHeights['%s_rocker_deg' % (corner)] = [row['rocker_theta_deg'] for row in solved]
        rideHeights['%s_kin_residual' % (corner)] = [row['residual_norm'] for row in solved]
        rideHeights['%s_rack_travel' % (corner)] = [row['rack_travel'] for row in solved]

        rideHeights['%s_wc_x' % (corner)] = [row['wheel_center'][0] for row in solved]
        rideHeights['%s_wc_y' % (corner)] = [row['wheel_center'][1] for row in solved]
        rideHeights['%s_wc_z' % (corner)] = [row['wheel_center'][2] for row in solved]

        rideHeights['%s_rvec_x' % (corner)] = [row['state'].rvec[0] for row in solved]
        rideHeights['%s_rvec_y' % (corner)] = [row['state'].rvec[1] for row in solved]
        rideHeights['%s_rvec_z' % (corner)] = [row['state'].rvec[2] for row in solved]

    return rideHeights


def calculateRideHeights(fullCaseSetupDict):
    '''
    Calculate ride heights and update dataframe with pitch, roll angles and wheel movements.
    
    Since CFD is in vehicle coordinate system (tunnel/wheels move around car body):
    - Positive pitch: front goes up, rear goes down
    - Positive roll: right side goes up, left side goes down
    - Heave: vertical translation of entire car
    
    :param fullCaseSetupDict: full case setup dictionary from caseSetupDict file
    :return: rideHeights dataframe with added columns: pitch, roll, wheel_fl, wheel_fr, wheel_rl, wheel_rr
    '''
    print('\tCalculating ride heights...')
    #get reference values
    REFCOR = fullCaseSetupDict['BC_SETUP']['REFCOR']
    REFWIDTH = float(fullCaseSetupDict['BC_SETUP']['REFWIDTH'][0])
    REFLEN = float(fullCaseSetupDict['BC_SETUP']['REFLEN'][0])
    


    #get rideheight file
    RIDE_HEIGHT_FILE = fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RIDE_HEIGHT_FILE'][0]
    if os.path.exists(RIDE_HEIGHT_FILE):
        print('\t\tReading ride height file from: %s' % RIDE_HEIGHT_FILE)
        rideHeightFilePath = RIDE_HEIGHT_FILE
    else:
        execDir = os.path.dirname(os.path.realpath(__file__))
        rideHeightFilePath = os.path.join(execDir,'rideHeightUtils',RIDE_HEIGHT_FILE)
        print('\t\tReading ride height file from: %s' % rideHeightFilePath)
    if not os.path.exists(rideHeightFilePath):
        sys.exit('ERROR! Ride height file %s cannot be found!' % (rideHeightFilePath))
    rideHeights = pd.read_csv(rideHeightFilePath)

    #converting to meters if units are mm
    if fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RH_UNIT'][0].lower() == 'mm':
        rideHeights[['fl','fr','rl','rr']] = rideHeights[['fl','fr','rl','rr']].multiply(0.001)
    elif fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RH_UNIT'][0].lower() == 'm':
        rideHeights[['fl','fr','rl','rr']] = rideHeights[['fl','fr','rl','rr']].multiply(1)

    #get ride height settings
    INIT_RH = fullCaseSetupDict['RIDE_HEIGHT_SETUP']['INIT_RH'][0]
    if INIT_RH == '':
        print('\t\tNo initial ride height provided, assuming delta positions at wheel ground points.')
        initFL = 0
        initFR = 0
        initRL = 0
        initRR = 0
    elif len(INIT_RH) == 4:
        initFL = float(INIT_RH[0])
        initFR = float(INIT_RH[1])
        initRL = float(INIT_RH[2])
        initRR = float(INIT_RH[3])
    else:
        sys.exit('ERROR! Invalid initial ride height provided, please provide 4 values for FL, FR, RL, RR or leave blank to assume ride height points are wheel movements!')

    if fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RH_REF_WIDTH'][0] == '':
        bw = [REFWIDTH,REFWIDTH]
    else:
        bw = [float(fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RH_REF_WIDTH'][0]),float(fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RH_REF_WIDTH'][0])]

    if fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RH_REF_LEN'][0] == '':
        bl = REFLEN
    else:        
        bl = float(fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RH_REF_LEN'][0])

    # Initialize new columns in dataframe
    rideHeights['pitch'] = 0.0
    rideHeights['roll'] = 0.0
    rideHeights['heave'] = 0.0
    rideHeights['wheel_fl'] = 0.0
    rideHeights['wheel_fr'] = 0.0
    rideHeights['wheel_rl'] = 0.0
    rideHeights['wheel_rr'] = 0.0
    rideHeights['tunnel_pitch'] = 0.0
    rideHeights['tunnel_roll'] = 0.0
    rideHeights['tunnel_heave'] = 0.0

    for idx, row in rideHeights.iterrows():
        # Calculate roll angles at front and rear
        fr_roll_angle, fr_dz_center = calculateRollAngles(row['fl'],row['fr'],bw[0])
        rr_roll_angle, rr_dz_center = calculateRollAngles(row['rl'],row['rr'],bw[1])
        
        #check if front and rear roll angles are the same
        drollAngle = abs(fr_roll_angle - rr_roll_angle)
        if drollAngle > 0.05:
            sys.exit('Error! Front and rear roll angles mis-match for point: %s' % (row['point']))
        
        # Calculate pitch angle from front and rear center heights
        pitch_angle, dz_wb_center = calculatePitchAngles(fr_dz_center, rr_dz_center, bl, REFCOR)
        
        # Store pitch and roll angles
        rideHeights.loc[idx, 'pitch'] = pitch_angle
        rideHeights.loc[idx, 'roll'] = fr_roll_angle  # Use front roll angle since they should match
        rideHeights.loc[idx, 'heave'] = dz_wb_center  # Heave is center displacement
        
        # wheel movements in tunnel coords - in tunnel frame the reference moves
        # relative to the wheels (inverse of the vehicle-frame motion)
        wheel_fl, wheel_fr, wheel_rl, wheel_rr, tunnel_pitch, tunnel_roll, tunnel_dz = calculateWheelMovements(
            row['fl'], row['fr'], row['rl'], row['rr'],
            pitch_angle, fr_roll_angle, dz_wb_center
        )
        
        rideHeights.loc[idx, 'wheel_fl'] = wheel_fl
        rideHeights.loc[idx, 'wheel_fr'] = wheel_fr
        rideHeights.loc[idx, 'wheel_rl'] = wheel_rl
        rideHeights.loc[idx, 'wheel_rr'] = wheel_rr
        rideHeights.loc[idx, 'tunnel_pitch'] = tunnel_pitch
        rideHeights.loc[idx, 'tunnel_roll'] = tunnel_roll
        rideHeights.loc[idx, 'tunnel_heave'] = tunnel_dz

    rideHeights = _append_suspension_kinematics(rideHeights, fullCaseSetupDict)

    print('\nUpdated rideHeights dataframe:')
    print(rideHeights)
    
    return rideHeights

def createRideHeightCases(rideHeights, fullCaseSetupDict):
    '''
    Creates the cases for each ride height point by copying the base case and modifying the necessary files with the new ride height information.

    The cases are created inside each case directory with the name of the ride height point. 
    For example, if the base case is "baseCase" and there is a ride height point 1, then the new case will be created as "baseCase_1".
    The caseSetup file is updated with the new ride height information for each case, the DOMAIN_PITCH and DOMAIN_ROLL are updated with the calculated pitch and roll angles for each ride height point, and the wheel movements are updated in the appropriate files 
    
    :param case: Description
    :param rideHeights: Description
    :param fullCaseSetupDict: Description
    '''
    baseCaseDir = os.path.join(os.getcwd())
    baseCaseName = os.path.basename(baseCaseDir)
    print('\t\tCreating ride height cases based on base case: %s' % (baseCaseName))
    rideHeights['caseName'] = ''
    updatedRideHeights = pd.DataFrame(columns=rideHeights.columns)
    baseTunnelPitch = float(fullCaseSetupDict['BC_SETUP']['DOMAIN_PITCH'][0])
    baseTunnelRoll = float(fullCaseSetupDict['BC_SETUP']['DOMAIN_ROLL'][0])
    #baseREFCOR = float(fullCaseSetupDict['BC_SETUP']['REFCOR'][0])
    currentREFCOR = fullCaseSetupDict['BC_SETUP']['REFCOR']
    SIM_SYM = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM']

    # srf cornering: ride-height map drives the per-point corner radius/direction from
    # optional columns, falling back to base [CORNERING_SETUP] when missing
    runCornering = ('CORNERING_SETUP' in fullCaseSetupDict and
                    fullCaseSetupDict['CORNERING_SETUP']['RUN_CORNERING'][0].lower() == 'true')
    cornerRadiusCol = None
    cornerDirCol = None
    if runCornering:
        for cand in ('corner_radius', 'cornerradius', 'corner_r', 'radius'):
            if cand in rideHeights.columns:
                cornerRadiusCol = cand
                break
        for cand in ('corner_dir', 'cornerdir', 'turn_dir', 'direction'):
            if cand in rideHeights.columns:
                cornerDirCol = cand
                break
        if cornerRadiusCol is None:
            print("\t\tWARNING! RUN_CORNERING is True but no corner radius column (corner_radius) was "
                  "found in the ride height file; all ride-height cases will use the base "
                  "CORNER_RADIUS = %s." % (fullCaseSetupDict['CORNERING_SETUP']['CORNER_RADIUS'][0]))
        else:
            print("\t\tCornering enabled: per-point corner radius taken from ride height column '%s'."
                  % (cornerRadiusCol))

    for idx, row in rideHeights.iterrows():
        
        if fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RUN_RH_POINTS'][0] == '':
            print('\t\tRide Height Point: %s' % (int(row['point'])))
        elif str(int(row['point'])) not in fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RUN_RH_POINTS']:
            continue
        
        if row['yaw'] != 0 and SIM_SYM[0].lower() == 'half':
            print('\t\t\tSkipping ride height point %s since it has non-zero yaw and simulation is set to half symmetry...' % (int(row['point'])))
            continue
        elif row['roll'] != 0 and SIM_SYM[0].lower() == 'half':
            print('\t\t\tSkipping ride height point %s since it has non-zero roll and simulation is set to half symmetry...' % (int(row['point'])))
            continue
        newCaseName = "%s_%s" % (baseCaseName, int(row['point']))
        newCaseDir = os.path.join(os.getcwd(), newCaseName)
        if not os.path.exists(newCaseDir):
            os.makedirs(newCaseDir)
        
        #write new caseSetup file with updated pitch and roll angles for domain
        newCaseSetupPath = os.path.join(newCaseDir,'caseSetup')
        updatedFullCaseSetupDict = deepcopy(fullCaseSetupDict)
        updatedFullCaseSetupDict['BC_SETUP']['DOMAIN_PITCH'] = baseTunnelPitch + float(row['tunnel_pitch'])
        updatedFullCaseSetupDict['BC_SETUP']['DOMAIN_ROLL'] = baseTunnelRoll + float(row['tunnel_roll'])
        
        updatedFullCaseSetupDict['BC_SETUP']['REFCOR'] = [currentREFCOR[0], currentREFCOR[1], str(float(currentREFCOR[2]) + row['tunnel_heave'])]

        #per-point corner radius / direction from the ride-height map (cornering only)
        if runCornering:
            if cornerRadiusCol is not None:
                rawRadius = row[cornerRadiusCol]
                try:
                    radiusVal = float(rawRadius)
                    if radiusVal <= 0:
                        raise ValueError('corner radius must be positive')
                    updatedFullCaseSetupDict['CORNERING_SETUP']['CORNER_RADIUS'] = [str(radiusVal)]
                except Exception:
                    print('\t\t\tWARNING! Point %s has invalid corner radius %r; using base '
                          'CORNER_RADIUS = %s.' % (int(row['point']), rawRadius,
                          fullCaseSetupDict['CORNERING_SETUP']['CORNER_RADIUS'][0]))
            if cornerDirCol is not None:
                rawDir = str(row[cornerDirCol]).strip().lower()
                if rawDir in ('left', 'right'):
                    updatedFullCaseSetupDict['CORNERING_SETUP']['CORNER_DIR'] = [rawDir]
                elif rawDir not in ('', 'nan'):
                    print('\t\t\tWARNING! Point %s has invalid corner direction %r; using base '
                          'CORNER_DIR = %s.' % (int(row['point']), rawDir,
                          fullCaseSetupDict['CORNERING_SETUP']['CORNER_DIR'][0]))
            #fail fast with a clear per-point message if the resolved radius is too small for the box
            checkCorneringDomain(updatedFullCaseSetupDict)

        writeToRHCaseSetup(updatedFullCaseSetupDict,newCaseSetupPath)
        rideHeights.loc[idx, 'caseName'] = newCaseName
        updatedRideHeights = pd.concat([updatedRideHeights, rideHeights.loc[idx].to_frame().T], ignore_index=True)


        print('\t\t\tPoint %s - pitch: %s, roll: %s, heave: %s' % (int(row['point']), row['tunnel_pitch'], row['tunnel_roll'], row['tunnel_heave']))
    
    #write out the updated rideHeights dataframe with case names for each point
    rideHeightsOutputPath = os.path.join(os.getcwd(),'rideHeights_updated.csv')
    updatedRideHeights.to_csv(rideHeightsOutputPath, index=False)

    return updatedRideHeights
def _categorize_pid_keywords_by_component(cornerPidKeywords):
    '''
    Categorize PID keywords by component type based on keyword content.
    Returns dict mapping corner -> {component_type -> [keywords]}
    '''
    compTypes = ['UCA', 'LCA', 'ROCKER', 'PUSHROD', 'DAMPER', 'WHEEL', 'TIE']
    # accept common abbreviations too. aliases are matched as substrings of the upper-cased keyword
    compAliases = {
        'UCA': ['UCA'],
        'LCA': ['LCA'],
        'ROCKER': ['ROCKER', 'ROCK', 'RKR'],
        'PUSHROD': ['PUSHROD', 'PROD', 'PSHRD', 'PUSH', 'PR'],
        'DAMPER': ['DAMPER', 'DAMP', 'SHOCK', 'STRUT'],
        'WHEEL': ['WHEEL', 'WHL'],
        'TIE': ['TIE', 'TIEROD', 'TROD'],
    }
    categorized = {}
    
    for corner, keywords in cornerPidKeywords.items():
        categorized[corner] = {ct: [] for ct in compTypes}
        for kw in keywords:
            kw_upper = kw.upper()
            for ct in compTypes:
                if any(alias in kw_upper for alias in compAliases[ct]):
                    categorized[corner][ct].append(kw)
                    break
    
    return categorized


def _validate_component_transforms(compTransforms, corner):
    #warn if any transform exceeds sane bounds. returns a list of warnings
    warnings = []
    
    for compType, transform in compTransforms.items():
        # rotation magnitude
        rvec = np.array(transform.get('rotation_rvec', [0, 0, 0]))
        rot_angle_rad = np.linalg.norm(rvec)
        rot_angle_deg = np.degrees(rot_angle_rad)
        
        # sane bounds per component
        if compType == 'WHEEL':
            max_rot = 15.0  # camber/toe
        elif compType in ['UCA', 'LCA']:
            max_rot = 20.0  # arm rotation
        elif compType == 'TIE':
            max_rot = 25.0  # tie rod
        elif compType == 'ROCKER':
            max_rot = 30.0  # rocker
        else:
            max_rot = 25.0  # default
        
        if rot_angle_deg > max_rot:
            warnings.append(f'{compType} ({corner}): rotation {rot_angle_deg:.2f}° exceeds typical bound {max_rot}°')
        
        # translation magnitude
        trans = np.array(transform.get('translation', [0, 0, 0]))
        trans_mag = np.linalg.norm(trans)
        if trans_mag > 0.5:  # >50cm looks wrong
            warnings.append(f'{compType} ({corner}): translation {trans_mag:.4f}m exceeds typical bound 0.5m')
    
    return warnings


def _write_component_transforms_log(childName, corner, row, compTransforms):
    #dump per-component transforms to csv for inspection
    logPath = os.path.join(childName, f'component_transforms_{corner}.csv')
    
    # make sure parent dir exists
    os.makedirs(os.path.dirname(logPath), exist_ok=True)
    
    rows = []
    for compType, transform in compTransforms.items():
        rvec = transform.get('rotation_rvec', [0, 0, 0])
        rot_mag = np.linalg.norm(rvec)
        trans = transform.get('translation', [0, 0, 0])
        pivot = transform.get('pivot', [0, 0, 0])
        morph = transform.get('morph', None)
        damper_delta = transform.get('damper_delta', np.nan)
        damper_len_static = transform.get('damper_length_static', np.nan)
        damper_len_current = transform.get('damper_length_current', np.nan)
        damper_morph_factor = transform.get('damper_morph_factor', np.nan)
        
        rows.append({
            'component': compType,
            'translation_x': trans[0],
            'translation_y': trans[1],
            'translation_z': trans[2],
            'translation_mag': np.linalg.norm(trans),
            'rotation_rvec_x': rvec[0],
            'rotation_rvec_y': rvec[1],
            'rotation_rvec_z': rvec[2],
            'rotation_angle_rad': rot_mag,
            'rotation_angle_deg': np.degrees(rot_mag),
            'pivot_x': pivot[0],
            'pivot_y': pivot[1],
            'pivot_z': pivot[2],
            'morph_mode': '' if morph is None else morph.get('mode', ''),
            'damper_delta': damper_delta,
            'damper_length_static': damper_len_static,
            'damper_length_current': damper_len_current,
            'damper_morph_factor': damper_morph_factor,
        })
    
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(logPath, index=False)


def _plot_suspension_linkage(kinSetup, corner, row, compTransforms, childName, rideHeights):
    #plot the linkage before/after to eyeball the kinematics, save a png
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
    except ImportError:
        print(f'\t\t\tWARNING: matplotlib not available, skipping linkage visualization for {corner}')
        return
    
    hardpoints = kinSetup['corners'][corner]

    def _rotate_point_about_pivot(point, pivot, rvec):
        theta = np.linalg.norm(rvec)
        if theta <= 1e-15:
            return point.copy()
        k = rvec / theta
        v = point - pivot
        v_rot = (
            v * np.cos(theta)
            + np.cross(k, v) * np.sin(theta)
            + k * np.dot(k, v) * (1.0 - np.cos(theta))
        )
        return pivot + v_rot
    
    # Prepare figure
    fig = plt.figure(figsize=(14, 6))
    
    # ===== BEFORE (left plot) =====
    ax1 = fig.add_subplot(121, projection='3d')
    
    # Static configuration
    colors = {'UCA': 'red', 'LCA': 'blue', 'ROCKER': 'green', 'TIE': 'orange', 'WHEEL': 'black'}
    
    # UCA
    if 'uca_f_inner' in hardpoints:
        uca_fi = np.array(hardpoints['uca_f_inner'])
        uca_ri = np.array(hardpoints['uca_r_inner'])
        uca_o = np.array(hardpoints['uca_outer_static'])
        ax1.plot([uca_fi[0], uca_ri[0]], [uca_fi[1], uca_ri[1]], [uca_fi[2], uca_ri[2]], 'r-', linewidth=2, label='UCA chassis')
        ax1.plot([uca_ri[0], uca_o[0]], [uca_ri[1], uca_o[1]], [uca_ri[2], uca_o[2]], 'r--', linewidth=1.5, label='UCA outer')
        ax1.scatter(*uca_fi, color='red', s=50)
        ax1.scatter(*uca_ri, color='red', s=50)
        ax1.scatter(*uca_o, color='red', s=100, marker='s')
    
    # LCA
    if 'lca_f_inner' in hardpoints:
        lca_fi = np.array(hardpoints['lca_f_inner'])
        lca_ri = np.array(hardpoints['lca_r_inner'])
        lca_o = np.array(hardpoints['lca_outer_static'])
        ax1.plot([lca_fi[0], lca_ri[0]], [lca_fi[1], lca_ri[1]], [lca_fi[2], lca_ri[2]], 'b-', linewidth=2, label='LCA chassis')
        ax1.plot([lca_ri[0], lca_o[0]], [lca_ri[1], lca_o[1]], [lca_ri[2], lca_o[2]], 'b--', linewidth=1.5, label='LCA outer')
        ax1.scatter(*lca_fi, color='blue', s=50)
        ax1.scatter(*lca_ri, color='blue', s=50)
        ax1.scatter(*lca_o, color='blue', s=100, marker='s')
    
    # ROCKER
    if 'rocker' in hardpoints:
        rocker_pivot = np.array(hardpoints['rocker']['pivot'])
        rocker_axis = np.array(hardpoints['rocker']['axis']) / np.linalg.norm(hardpoints['rocker']['axis'])
        ax1.scatter(*rocker_pivot, color='green', s=150, marker='^', label='Rocker pivot')
        # Draw rocker axis
        rocker_axis_end = rocker_pivot + rocker_axis * 0.2
        ax1.plot([rocker_pivot[0], rocker_axis_end[0]], [rocker_pivot[1], rocker_axis_end[1]], [rocker_pivot[2], rocker_axis_end[2]], 'g--', linewidth=2)
    
    # TIE
    if 'tie_inner' in hardpoints:
        tie_i = np.array(hardpoints['tie_inner'])
        tie_o = np.array(hardpoints['tie_outer_static'])
        ax1.plot([tie_i[0], tie_o[0]], [tie_i[1], tie_o[1]], [tie_i[2], tie_o[2]], 'orange', linewidth=2, label='Tie rod')
        ax1.scatter(*tie_i, color='orange', s=50)
        ax1.scatter(*tie_o, color='orange', s=100, marker='s')

    # PUSHROD
    if 'rocker' in hardpoints and 'pushrod_outer_static' in hardpoints:
        pr_i = np.array(hardpoints['rocker']['pushrod_joint_ref'])
        pr_o = np.array(hardpoints['pushrod_outer_static'])
        ax1.plot([pr_i[0], pr_o[0]], [pr_i[1], pr_o[1]], [pr_i[2], pr_o[2]], 'purple', linewidth=2, label='Pushrod')

    # DAMPER
    if 'rocker' in hardpoints:
        dm_i = np.array(hardpoints['rocker']['damper_joint_ref'])
        dm_o = np.array(hardpoints['rocker']['damper_chassis'])
        ax1.plot([dm_i[0], dm_o[0]], [dm_i[1], dm_o[1]], [dm_i[2], dm_o[2]], 'brown', linewidth=2, label='Damper')
    
    # WHEEL
    if 'wheel_center_static' in hardpoints:
        wc = np.array(hardpoints['wheel_center_static'])
        ax1.scatter(*wc, color='black', s=200, marker='o', label='Wheel center')
    
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_title(f'Suspension Linkage - STATIC ({corner.upper()})')
    ax1.legend(fontsize=8)
    
    # ===== AFTER (right plot) =====
    ax2 = fig.add_subplot(122, projection='3d')
    
    # Current configuration (apply transforms)
    if 'WHEEL' in compTransforms:
        wheel_trans = np.array(compTransforms['WHEEL']['translation'])
    else:
        wheel_trans = np.array([0, 0, 0])
    
    # UCA with rotation
    if 'uca_f_inner' in hardpoints and 'UCA' in compTransforms:
        uca_fi_curr = np.array(hardpoints['uca_f_inner'])  # Chassis point fixed
        uca_ri_curr = np.array(hardpoints['uca_r_inner'])
        uca_o_curr = np.array(hardpoints['uca_outer_static']) + wheel_trans
        
        ax2.plot([uca_fi_curr[0], uca_ri_curr[0]], [uca_fi_curr[1], uca_ri_curr[1]], [uca_fi_curr[2], uca_ri_curr[2]], 'r-', linewidth=2)
        ax2.plot([uca_ri_curr[0], uca_o_curr[0]], [uca_ri_curr[1], uca_o_curr[1]], [uca_ri_curr[2], uca_o_curr[2]], 'r--', linewidth=1.5)
        ax2.scatter(*uca_fi_curr, color='red', s=50)
        ax2.scatter(*uca_ri_curr, color='red', s=50)
        ax2.scatter(*uca_o_curr, color='red', s=100, marker='s')
    
    # LCA with rotation
    if 'lca_f_inner' in hardpoints and 'LCA' in compTransforms:
        lca_fi_curr = np.array(hardpoints['lca_f_inner'])
        lca_ri_curr = np.array(hardpoints['lca_r_inner'])
        lca_o_curr = np.array(hardpoints['lca_outer_static']) + wheel_trans
        
        ax2.plot([lca_fi_curr[0], lca_ri_curr[0]], [lca_fi_curr[1], lca_ri_curr[1]], [lca_fi_curr[2], lca_ri_curr[2]], 'b-', linewidth=2)
        ax2.plot([lca_ri_curr[0], lca_o_curr[0]], [lca_ri_curr[1], lca_o_curr[1]], [lca_ri_curr[2], lca_o_curr[2]], 'b--', linewidth=1.5)
        ax2.scatter(*lca_fi_curr, color='blue', s=50)
        ax2.scatter(*lca_ri_curr, color='blue', s=50)
        ax2.scatter(*lca_o_curr, color='blue', s=100, marker='s')
    
    # ROCKER with rotation
    if 'rocker' in hardpoints and 'ROCKER' in compTransforms:
        rocker_pivot = np.array(hardpoints['rocker']['pivot'])
        rocker_axis = np.array(hardpoints['rocker']['axis']) / np.linalg.norm(hardpoints['rocker']['axis'])
        ax2.scatter(*rocker_pivot, color='green', s=150, marker='^')
        rocker_axis_end = rocker_pivot + rocker_axis * 0.2
        ax2.plot([rocker_pivot[0], rocker_axis_end[0]], [rocker_pivot[1], rocker_axis_end[1]], [rocker_pivot[2], rocker_axis_end[2]], 'g--', linewidth=2)
    
    # TIE with rotation
    if 'tie_inner' in hardpoints and 'TIE' in compTransforms:
        tie_i = np.array(hardpoints['tie_inner'])
        tie_o_curr = np.array(hardpoints['tie_outer_static']) + wheel_trans
        ax2.plot([tie_i[0], tie_o_curr[0]], [tie_i[1], tie_o_curr[1]], [tie_i[2], tie_o_curr[2]], 'orange', linewidth=2)
        ax2.scatter(*tie_i, color='orange', s=50)
        ax2.scatter(*tie_o_curr, color='orange', s=100, marker='s')

    # PUSHROD current
    if 'rocker' in hardpoints and 'pushrod_outer_static' in hardpoints and 'ROCKER' in compTransforms:
        rocker_pivot = np.array(hardpoints['rocker']['pivot'])
        rocker_rvec = np.array(compTransforms['ROCKER']['rotation_rvec'])
        pr_i_static = np.array(hardpoints['rocker']['pushrod_joint_ref'])
        pr_i_curr = _rotate_point_about_pivot(pr_i_static, rocker_pivot, rocker_rvec)
        pr_o_curr = np.array(hardpoints['pushrod_outer_static']) + wheel_trans
        ax2.plot([pr_i_curr[0], pr_o_curr[0]], [pr_i_curr[1], pr_o_curr[1]], [pr_i_curr[2], pr_o_curr[2]], 'purple', linewidth=2)

    # DAMPER current
    if 'rocker' in hardpoints and 'ROCKER' in compTransforms:
        rocker_pivot = np.array(hardpoints['rocker']['pivot'])
        rocker_rvec = np.array(compTransforms['ROCKER']['rotation_rvec'])
        dm_i_static = np.array(hardpoints['rocker']['damper_joint_ref'])
        dm_i_curr = _rotate_point_about_pivot(dm_i_static, rocker_pivot, rocker_rvec)
        dm_o = np.array(hardpoints['rocker']['damper_chassis'])
        ax2.plot([dm_i_curr[0], dm_o[0]], [dm_i_curr[1], dm_o[1]], [dm_i_curr[2], dm_o[2]], 'brown', linewidth=2)
    
    # WHEEL with rotation
    if 'wheel_center_static' in hardpoints and 'WHEEL' in compTransforms:
        wc_curr = np.array(hardpoints['wheel_center_static']) + wheel_trans
        ax2.scatter(*wc_curr, color='black', s=200, marker='o')
    
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.set_title(f'Suspension Linkage - CURRENT ({corner.upper()})')
    
    plt.tight_layout()
    
    # Save image
    imagePath = os.path.join(childName, f'linkage_validation_{corner}.png')
    # Ensure parent directory exists
    os.makedirs(os.path.dirname(imagePath), exist_ok=True)
    plt.savefig(imagePath, dpi=100, bbox_inches='tight')
    print(f'\t\t\t\tSaved linkage plot: {imagePath}')
    plt.close()


def _create_linkage_gifs(baseCaseName, points, corners=('fl', 'fr', 'rl', 'rr'), frameDurationMs=450):
    #build a gif per corner from the per-point linkage pngs
    try:
        from PIL import Image
    except ImportError:
        print('\t\t\tWARNING: Pillow not available, skipping GIF creation.')
        return

    for corner in corners:
        framePaths = []
        for point in points:
            p = os.path.join(f'{baseCaseName}_{int(point)}', f'linkage_validation_{corner}.png')
            if os.path.exists(p):
                framePaths.append(p)

        if len(framePaths) == 0:
            continue

        # keep only the first and last frames
        if len(framePaths) > 1:
            framePaths = [framePaths[0], framePaths[-1]]

        images = [Image.open(fp).convert('P', palette=Image.ADAPTIVE) for fp in framePaths]
        outPath = os.path.join(os.getcwd(), f'linkage_animation_{corner}.gif')
        images[0].save(
            outPath,
            save_all=True,
            append_images=images[1:],
            duration=frameDurationMs,
            loop=0,
            optimize=False,
        )
        for img in images:
            img.close()
        print(f'\t\t\tSaved linkage animation GIF: {outPath}')


def _compute_component_transforms(kinSetup, corner, row, rideHeights):
    #per-component transforms from the solver state. returns component -> {translation, rotation_rvec, pivot}
    compTransforms = {}
    hardpoints = kinSetup['corners'][corner]

    def _rvec_from_vectors(v_static, v_current):
        v_static_norm = v_static / (np.linalg.norm(v_static) + 1e-15)
        v_current_norm = v_current / (np.linalg.norm(v_current) + 1e-15)
        rot_axis = np.cross(v_static_norm, v_current_norm)
        rot_axis_norm = np.linalg.norm(rot_axis)
        cos_angle = np.dot(v_static_norm, v_current_norm)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        rot_angle = np.arccos(cos_angle)
        if rot_axis_norm > 1e-15:
            return (rot_axis / rot_axis_norm) * rot_angle
        return np.array([0.0, 0.0, 0.0])

    def _rotate_point_about_pivot(point, pivot, rvec):
        theta = np.linalg.norm(rvec)
        if theta <= 1e-15:
            return point.copy()
        k = rvec / theta
        v = point - pivot
        v_rot = (
            v * np.cos(theta)
            + np.cross(k, v) * np.sin(theta)
            + k * np.dot(k, v) * (1.0 - np.cos(theta))
        )
        return pivot + v_rot

    def _rvec_about_axis_to_target(axisPoint, axisDir, p_from, p_to):
        #rvec about a fixed hinge line (through axisPoint, dir axisDir) that swings p_from
        #onto p_to. only perpendicular components count, so points on the line stay fixed
        axis = np.asarray(axisDir, dtype=np.float64)
        n = np.linalg.norm(axis)
        if n <= 1e-15:
            return np.array([0.0, 0.0, 0.0])
        axis = axis / n
        ap0 = np.asarray(axisPoint, dtype=np.float64)
        a = np.asarray(p_from, dtype=np.float64) - ap0
        b = np.asarray(p_to, dtype=np.float64) - ap0
        a_perp = a - np.dot(a, axis) * axis
        b_perp = b - np.dot(b, axis) * axis
        na = np.linalg.norm(a_perp)
        nb = np.linalg.norm(b_perp)
        if na <= 1e-15 or nb <= 1e-15:
            return np.array([0.0, 0.0, 0.0])
        a_hat = a_perp / na
        b_hat = b_perp / nb
        cos_t = np.clip(np.dot(a_hat, b_hat), -1.0, 1.0)
        sin_t = np.dot(np.cross(a_hat, b_hat), axis)
        theta = math.atan2(sin_t, cos_t)
        return axis * theta

    # wheel: follows wheel center, rotates with camber/toe from the solver
    wcx, wcy, wcz = '%s_wc_x' % corner, '%s_wc_y' % corner, '%s_wc_z' % corner
    if all(k in rideHeights.columns for k in [wcx, wcy, wcz]):
        staticWC = np.array(hardpoints['wheel_center_static'], dtype=np.float64)
        currentWC = np.array([float(row[wcx]), float(row[wcy]), float(row[wcz])], dtype=np.float64)
        wheelTranslation = (currentWC - staticWC).tolist()
        
        # wheel rotation from the solver (camber + toe)
        rvec_x_col, rvec_y_col, rvec_z_col = '%s_rvec_x' % corner, '%s_rvec_y' % corner, '%s_rvec_z' % corner
        if all(k in rideHeights.columns for k in [rvec_x_col, rvec_y_col, rvec_z_col]):
            wheelRvec = [float(row[rvec_x_col]), float(row[rvec_y_col]), float(row[rvec_z_col])]
        else:
            wheelRvec = [0.0, 0.0, 0.0]
        
        compTransforms['WHEEL'] = {
            'translation': wheelTranslation,
            'rotation_rvec': wheelRvec,
            'pivot': staticWC.tolist(),
        }
    
    # UCA & LCA: rigid A-arms that only rotate about the line through their two chassis
    # bushings. bushings stay fixed, outer ball joint swings to follow the wheel pose
    for compType in ['UCA', 'LCA']:
        inner_f_key = '%s_f_inner' % compType.lower()
        inner_r_key = '%s_r_inner' % compType.lower()
        outer_key = '%s_outer_static' % compType.lower()
        
        if all(k in hardpoints for k in [inner_f_key, inner_r_key, outer_key]):
            # static configuration
            innerF = np.array(hardpoints[inner_f_key], dtype=np.float64)
            innerR = np.array(hardpoints[inner_r_key], dtype=np.float64)
            outerStatic = np.array(hardpoints[outer_key], dtype=np.float64)
            
            # new outer ball-joint position from the full wheel pose
            if 'WHEEL' in compTransforms:
                wheelPivot = np.array(compTransforms['WHEEL']['pivot'], dtype=np.float64)
                wheelRvecArr = np.array(compTransforms['WHEEL']['rotation_rvec'], dtype=np.float64)
                wheelTrans = np.array(compTransforms['WHEEL']['translation'], dtype=np.float64)
                outerNew = _rotate_point_about_pivot(outerStatic, wheelPivot, wheelRvecArr) + wheelTrans
            else:
                outerNew = outerStatic
            
            # hinge axis is the line through the two bushings; rotate about it to carry
            # the outer joint from static to new position
            hingeAxis = innerR - innerF
            rvec = _rvec_about_axis_to_target(innerF, hingeAxis, outerStatic, outerNew)
            
            # pivot is a bushing on the hinge axis, so both bushings stay fixed
            compTransforms[compType] = {
                'translation': [0.0, 0.0, 0.0],
                'rotation_rvec': rvec.tolist(),
                'pivot': innerF.tolist(),
            }
    
    # ROCKER: rotates about ROCKER_PIVOT
    if 'rocker' in hardpoints:
        rocker_pivot = np.array(hardpoints['rocker']['pivot'], dtype=np.float64)
        
        # rocker angle from the solver
        rocker_deg_col = '%s_rocker_deg' % corner
        if rocker_deg_col in rideHeights.columns:
            rocker_angle = float(row[rocker_deg_col])
            
            # rocker axis
            rocker_axis = np.array(hardpoints['rocker']['axis'], dtype=np.float64)
            rocker_axis = rocker_axis / np.linalg.norm(rocker_axis)  # normalize
            
            # angle -> rotation vector (axis * radians)
            rvec = rocker_axis * np.radians(rocker_angle)
            
            compTransforms['ROCKER'] = {
                'translation': [0.0, 0.0, 0.0],
                'rotation_rvec': rvec.tolist(),
                'pivot': rocker_pivot.tolist(),
            }
    
    # PUSHROD: one end on rocker, one end on wheel side. solve from endpoint motion
    if 'ROCKER' in compTransforms and 'pushrod_outer_static' in hardpoints and 'rocker' in hardpoints and 'WHEEL' in compTransforms:
        pushrodInnerStatic = np.array(hardpoints['rocker']['pushrod_joint_ref'], dtype=np.float64)
        pushrodOuterStatic = np.array(hardpoints['pushrod_outer_static'], dtype=np.float64)

        rockerPivot = np.array(hardpoints['rocker']['pivot'], dtype=np.float64)
        rockerRvec = np.array(compTransforms['ROCKER']['rotation_rvec'], dtype=np.float64)
        rockerTrans = np.array(compTransforms['ROCKER']['translation'], dtype=np.float64)
        pushrodInnerCurrent = _rotate_point_about_pivot(pushrodInnerStatic, rockerPivot, rockerRvec) + rockerTrans

        # pushrod outer end is bolted to the UCA, so it follows the UCA hinge motion, not
        # the upright. using the upright pose would add steer and a spurious fore/aft shift
        if 'UCA' in compTransforms:
            ucaPivot = np.array(compTransforms['UCA']['pivot'], dtype=np.float64)
            ucaRvec = np.array(compTransforms['UCA']['rotation_rvec'], dtype=np.float64)
            ucaTrans = np.array(compTransforms['UCA']['translation'], dtype=np.float64)
            pushrodOuterCurrent = _rotate_point_about_pivot(pushrodOuterStatic, ucaPivot, ucaRvec) + ucaTrans
        else:
            wheelPivot = np.array(compTransforms['WHEEL']['pivot'], dtype=np.float64)
            wheelRvec = np.array(compTransforms['WHEEL']['rotation_rvec'], dtype=np.float64)
            wheelTrans = np.array(compTransforms['WHEEL']['translation'], dtype=np.float64)
            pushrodOuterCurrent = _rotate_point_about_pivot(pushrodOuterStatic, wheelPivot, wheelRvec) + wheelTrans

        v_static = pushrodOuterStatic - pushrodInnerStatic
        v_current = pushrodOuterCurrent - pushrodInnerCurrent
        pushrodRvec = _rvec_from_vectors(v_static, v_current)

        compTransforms['PUSHROD'] = {
            'translation': (pushrodInnerCurrent - pushrodInnerStatic).tolist(),
            'rotation_rvec': pushrodRvec.tolist(),
            'pivot': pushrodInnerStatic.tolist(),
        }

    # DAMPER: chassis end fixed, rocker end follows the rocker; apply stroke morph from solver delta
    if 'ROCKER' in compTransforms and 'rocker' in hardpoints:
        damperJointStatic = np.array(hardpoints['rocker']['damper_joint_ref'], dtype=np.float64)
        damperChassis = np.array(hardpoints['rocker']['damper_chassis'], dtype=np.float64)

        rockerPivot = np.array(hardpoints['rocker']['pivot'], dtype=np.float64)
        rockerRvec = np.array(compTransforms['ROCKER']['rotation_rvec'], dtype=np.float64)
        rockerTrans = np.array(compTransforms['ROCKER']['translation'], dtype=np.float64)
        damperJointCurrent = _rotate_point_about_pivot(damperJointStatic, rockerPivot, rockerRvec) + rockerTrans

        v_static = damperJointStatic - damperChassis
        v_current = damperJointCurrent - damperChassis
        damperRvec = _rvec_from_vectors(v_static, v_current)

        static_len = float(np.linalg.norm(v_static))
        damper_delta_col = '%s_damper_delta' % corner
        if damper_delta_col in rideHeights.columns:
            damper_delta = float(row[damper_delta_col])
        else:
            damper_delta = float(np.linalg.norm(v_current) - static_len)
        current_len = static_len + damper_delta
        if static_len > 1e-15:
            morph_factor = max(current_len / static_len, 1e-6)
        else:
            morph_factor = 1.0

        compTransforms['DAMPER'] = {
            'translation': [0.0, 0.0, 0.0],
            'rotation_rvec': damperRvec.tolist(),
            'pivot': damperChassis.tolist(),
            'morph': {
                'mode': 'vector_scale',
                'vector': v_static.tolist(),
                'factor': morph_factor,
                'origin': damperChassis.tolist(),
                'limits': {'min': 0.0, 'max': static_len},
            },
            'damper_delta': damper_delta,
            'damper_length_static': static_len,
            'damper_length_current': current_len,
            'damper_morph_factor': morph_factor,
        }
    
    # TIE: inner (rack end) slides with rack travel, outer is rigid on the upright and
    # follows the wheel pose. tie length is preserved, so rebuild as a rotation about the inner end
    if 'tie_inner' in hardpoints and 'tie_outer_static' in hardpoints and 'WHEEL' in compTransforms:
        tieInnerStatic = np.array(hardpoints['tie_inner'], dtype=np.float64)
        tieOuterStatic = np.array(hardpoints['tie_outer_static'], dtype=np.float64)

        staticWC = np.array(compTransforms['WHEEL']['pivot'], dtype=np.float64)
        wheelRvec = np.array(compTransforms['WHEEL']['rotation_rvec'], dtype=np.float64)
        wheelTrans = np.array(compTransforms['WHEEL']['translation'], dtype=np.float64)

        # tie outer is on the upright -> apply full wheel pose
        tieOuterCurrent = _rotate_point_about_pivot(tieOuterStatic, staticWC, wheelRvec) + wheelTrans

        # tie inner (rack end) slides along the rack axis by the solved rack travel
        rackAxis = np.array(hardpoints.get('rack_axis', [0.0, 1.0, 0.0]), dtype=np.float64)
        rackNorm = np.linalg.norm(rackAxis)
        if rackNorm > 1e-15:
            rackAxis = rackAxis / rackNorm
        else:
            rackAxis = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        rack_col = '%s_rack_travel' % corner
        rackTravel = float(row[rack_col]) if rack_col in rideHeights.columns else 0.0
        tieInnerCurrent = tieInnerStatic + rackTravel * rackAxis

        v_static = tieOuterStatic - tieInnerStatic
        v_current = tieOuterCurrent - tieInnerCurrent
        rvec = _rvec_from_vectors(v_static, v_current)

        compTransforms['TIE'] = {
            'translation': (tieInnerCurrent - tieInnerStatic).tolist(),
            'rotation_rvec': rvec.tolist(),
            'pivot': tieInnerStatic.tolist(),
        }

    return compTransforms


def transformGeom(fullCaseSetupDict,rideHeights,geomDict):
    baseCaseDir = os.getcwd()
    baseCaseName = os.path.basename(baseCaseDir)
    #base domain tilt + ground height, used to build each point's ground-zone transform
    baseDomainPitch = float(fullCaseSetupDict['BC_SETUP']['DOMAIN_PITCH'][0])
    baseDomainRoll = float(fullCaseSetupDict['BC_SETUP']['DOMAIN_ROLL'][0])
    baseREFCOR = fullCaseSetupDict['BC_SETUP']['REFCOR']
    if 'ansa' in fullCaseSetupDict['GLOBAL_REFINEMENT']['TEMPLATE_TYPE'][0].lower():
        USE_ANSA = True
    else:
        USE_ANSA = False

    kinSetup = _load_suspension_kinematics_setup(fullCaseSetupDict)
    rideHeightSetup = fullCaseSetupDict.get('RIDE_HEIGHT_SETUP', {})
    requireSuspensionPidMatch = _truthy(rideHeightSetup.get('SUSP_REQUIRE_ALL_PIDS_MATCHED', ['false'])[0])

    # corner keywords from the hardpoint cfg, then optional caseSetup override/extension
    cornerPidKeywords = {'fl': [], 'fr': [], 'rl': [], 'rr': []}
    if kinSetup is not None and 'corner_pid_keywords' in kinSetup:
        for c in cornerPidKeywords.keys():
            cornerPidKeywords[c].extend(list(kinSetup['corner_pid_keywords'].get(c, [])))

    # de-duplicate while preserving order
    for c in cornerPidKeywords.keys():
        dedup = []
        seen = set()
        for k in cornerPidKeywords[c]:
            if k not in seen:
                dedup.append(k)
                seen.add(k)
        cornerPidKeywords[c] = dedup

    # categorize keywords by component type
    categorizedKeywords = _categorize_pid_keywords_by_component(cornerPidKeywords)

    # copies/transforms all geometries into each child case directory
    processedPoints = []
    for idx, row in rideHeights.iterrows():
        point = int(row['point'])
        childName = baseCaseName + '_' + str(point)
        transformedGeoms = set()
        processedPoints.append(point)

        print('\t\t\tTransforming geometries for point %s' % str(point))

        # build per-component, per-corner transforms from kinematics
        pidTransformDict = {}
        
        for corner in ['fr', 'fl', 'rr', 'rl']:
            if kinSetup is None or corner not in kinSetup['corners']:
                continue
            
            # compute transforms for each component type
            compTransforms = _compute_component_transforms(kinSetup, corner, row, rideHeights)
            
            # validate transforms
            validation_warnings = _validate_component_transforms(compTransforms, corner)
            if validation_warnings:
                for warn in validation_warnings:
                    print(f'\t\t\t\tWARNING: {warn}')
            
            # write transform log
            _write_component_transforms_log(childName, corner, row, compTransforms)
            
            # plot linkage before/after
            _plot_suspension_linkage(kinSetup, corner, row, compTransforms, childName, rideHeights)
            
            # map each component's keywords to its transform
            for compType, compTransform in compTransforms.items():
                keywords = categorizedKeywords[corner].get(compType, [])
                if len(keywords) == 0:
                    continue
                
                regex = '|'.join(re.escape(k) for k in keywords)
                pidTransformDict[regex] = compTransform
                print('\t\t\t\t%s (%s) transform: translation=[%1.6f, %1.6f, %1.6f]' % 
                      (compType, corner, compTransform['translation'][0], 
                       compTransform['translation'][1], compTransform['translation'][2]))

        for geom in geomDict.keys():
            geomFilePath = 'constant/triSurface/%s' % (geom)
            outputPath = os.path.join(childName, 'constant', 'triSurface', geom)

            # try PID-aware suspension transform first
            if len(pidTransformDict) > 0:
                try:
                    pidResult = transformGeometryByPIDRegex(
                        inputFile=geomFilePath,
                        outputFile=outputPath,
                        pidTransformDict=pidTransformDict,
                        pivot=None,
                        requireAllPIDsMatched=requireSuspensionPidMatch,
                        logBasePath=outputPath + '.susp_pid_transform'
                    )
                    if pidResult.get('matched_pid_count', 0) > 0:
                        transformedGeoms.add(geom)
                        continue
                except Exception as e:
                    # skip non-OBJ/STL or binary STL on the PID path, fall back to filename
                    print('\t\t\t\tPID scan skipped for %s: %s' % (geom, e))

        # legacy filename fallback for wheel files not PID-transformed (translation only)
        if kinSetup is not None:
            for corner in ['fr', 'fl', 'rr', 'rl']:
                if corner not in kinSetup['corners']:
                    continue
                compTransforms = _compute_component_transforms(kinSetup, corner, row, rideHeights)
                if 'WHEEL' in compTransforms:
                    wheelTransform = compTransforms['WHEEL']
                    tVec = wheelTransform['translation']
                    wheelRvec = wheelTransform.get('rotation_rvec', None)
                    wheelPivot = wheelTransform.get('pivot', None)
                    for geom in geomDict.keys():
                        if geom in transformedGeoms:
                            continue
                        geomFilePath = 'constant/triSurface/%s' % (geom)
                        if corner in geom.lower():
                            transformParts(geom, geomFilePath, tVec, childName,
                                           rotation_rvec=wheelRvec, pivot=wheelPivot)
                            transformedGeoms.add(geom)

        # copy all untouched files
        for geom in geomDict.keys():
            if geom in transformedGeoms:
                continue
            geomFilePath = 'constant/triSurface/%s' % (geom)
            geomChildPath = '%s/constant/triSurface/%s' % (childName,geom)
            if geom.startswith('GRND'):
                #ground zone follows the domain: heave up to the new REFCOR height, then tilt about
                #it (roll x then pitch y), matching the blockMesh transform for this point
                childPitch = math.radians(baseDomainPitch + float(row['tunnel_pitch']))
                childRoll = math.radians(baseDomainRoll + float(row['tunnel_roll']))
                childPivot = [float(baseREFCOR[0]), float(baseREFCOR[1]), float(baseREFCOR[2]) + float(row['tunnel_heave'])]
                tmp = geomChildPath + '.grndtmp'
                print('\t\t\tTransforming ground zone: %s' % (geom))
                transformGeometryPreservePID(geomFilePath,outputFile=geomChildPath,translation=[0,0,float(row['tunnel_heave'])])
                transformGeometryPreservePID(geomChildPath,outputFile=tmp,rotation_rvec=[childRoll,0,0],pivot=childPivot)
                transformGeometryPreservePID(tmp,outputFile=geomChildPath,rotation_rvec=[0,childPitch,0],pivot=childPivot)
                os.remove(tmp)
                continue
            print('\t\t\tCopying non-moving geometry: %s' % (geom))
            copyWithMkdir(geomFilePath, geomChildPath)

    # build corner animation gifs from the per-point linkage plots
    if len(processedPoints) > 0:
        _create_linkage_gifs(baseCaseName, processedPoints, corners=('fl', 'fr', 'rl', 'rr'))

    return rideHeights



def transformBlockMesh(fullCaseSetupDict):
    print('\t\t\tTransforming blockMesh based on domain rotations specified...')
    #get blockMeshDict file
    blockMeshDictPath = os.path.join(os.getcwd(),'system','blockMeshDict')
    if not os.path.exists(blockMeshDictPath):
        sys.exit('ERROR! blockMeshDict file cannot be found at %s!' % (blockMeshDictPath))
    domain_roll = float(fullCaseSetupDict['BC_SETUP']['DOMAIN_ROLL'][0])
    domain_pitch = float(fullCaseSetupDict['BC_SETUP']['DOMAIN_PITCH'][0])
    print('\t\t\t\tDomain Roll: %s, Domain Pitch: %s' % (domain_roll, domain_pitch))

    transformBlockmeshPoints(blockMeshDictPath,outputFile=os.path.join(os.getcwd(),'system','blockMeshDict'),rotation={'x':domain_roll, 'y':domain_pitch, 'z':0},REFCOR=fullCaseSetupDict['BC_SETUP']['REFCOR'])

    
   
def transformGroundGeom(geomDict,fullCaseSetupDict):
    #tilt ground-zone surfaces to follow the domain (roll about x then pitch about y, about
    #REFCOR), same convention as transformBlockMesh. no translation, ground moves with the mesh
    grndGeoms = [g for g in geomDict.keys() if g.startswith('GRND')]
    if len(grndGeoms) == 0:
        return
    domain_roll = math.radians(float(fullCaseSetupDict['BC_SETUP']['DOMAIN_ROLL'][0]))
    domain_pitch = math.radians(float(fullCaseSetupDict['BC_SETUP']['DOMAIN_PITCH'][0]))
    REFCOR = fullCaseSetupDict['BC_SETUP']['REFCOR']
    pivot = [float(REFCOR[0]), float(REFCOR[1]), float(REFCOR[2])]
    print('\t\t\tTilting ground zone surfaces to match domain (roll %s, pitch %s)...' %
          (math.degrees(domain_roll), math.degrees(domain_pitch)))
    for geom in grndGeoms:
        localPath = os.path.join('constant','triSurface',geom)
        if not os.path.exists(localPath):
            print('\t\t\t\tWARNING! %s not found, skipping tilt.' % (geom))
            continue
        #read the untilted source (master via symlink, or the local copy in a child case)
        srcPath = os.path.realpath(localPath)
        d = os.path.dirname(localPath)
        tmp1 = os.path.join(d,'_grndtilt1_'+geom)
        tmp2 = os.path.join(d,'_grndtilt2_'+geom)
        #roll about x then pitch about y, both about REFCOR (matches transformBlockmeshPoints)
        transformGeometryPreservePID(srcPath,outputFile=tmp1,rotation_rvec=[domain_roll,0,0],pivot=pivot)
        transformGeometryPreservePID(tmp1,outputFile=tmp2,rotation_rvec=[0,domain_pitch,0],pivot=pivot)
        if os.path.islink(localPath):
            os.remove(localPath)
        os.replace(tmp2,localPath)
        if os.path.exists(tmp1):
            os.remove(tmp1)
        print('\t\t\t\t%s tilted about REFCOR.' % (geom))


def transformParts(geom,geomFilePath, transformAmount,childCaseName, rotation_rvec=None, pivot=None):
    if isinstance(transformAmount, (list, tuple, np.ndarray)):
        tVec = [float(transformAmount[0]), float(transformAmount[1]), float(transformAmount[2])]
    else:
        tVec = [0.0, 0.0, float(transformAmount)]

    if rotation_rvec is not None and float(np.linalg.norm(rotation_rvec)) > 1e-15:
        print('\t\t\t\tMoving %s by translation [%s, %s, %s] and rotation_rvec [%s, %s, %s]' %
              (geom, tVec[0], tVec[1], tVec[2], rotation_rvec[0], rotation_rvec[1], rotation_rvec[2]))
    else:
        print('\t\t\t\tMoving %s by translation [%s, %s, %s]' % (geom,tVec[0],tVec[1],tVec[2]))
    # dummy transform
    try:
        outputPath = os.path.join(childCaseName,'constant','triSurface',geom)
        if not 'dummy' in geom:
            transformGeometryPreservePID(geomFilePath,outputFile=outputPath,translation=tVec,
                                         rotation_rvec=rotation_rvec,pivot=pivot)
        else:
            print('')
        print('\t\t\t\t\ttransform OK!')
    except Exception as E:
        print(E)
        sys.exit('\t\t\t\t\ttransform ERROR!')

    
    #copyWithMkdir(geomFilePath, outputPath) #copy the geometry to the new case directory with the same name, but transformed position
def transformGeometryPreservePID(inputFile, outputFile, rotation=None, translation=None, scale=None, morphing_dict=None, rotation_rvec=None, pivot=None):
    '''
    Transform geometry files (OBJ or STL, with or without .gz compression) while preserving
    PID naming (groups in OBJ, solid names in STL).
    
    For OBJ files: preserves 'g <name>' and 'o <name>' groupings
    For STL files: preserves 'solid <name>' / 'endsolid <name>' multi-solid structure
                   (binary STLs only have a single solid name in their header)
    
    :param inputFile: Path to input geometry file (.obj, .stl, .obj.gz, .stl.gz)
    :param outputFile: Path to output geometry file (same supported extensions)
    :param rotation: Dict with keys 'x', 'y', 'z' for rotation angles in degrees
    :param translation: List/array [dx, dy, dz] for translation offsets
    :param scale: Scalar or [sx, sy, sz] for scaling
    :param morphing_dict: Dict mapping vertex indices to displacement vectors
    :param rotation_rvec: Optional axis-angle rotation vector [rx, ry, rz] in radians,
        applied about ``pivot`` before the euler ``rotation`` and ``translation``
    :param pivot: Optional [px, py, pz] pivot point for ``rotation_rvec`` (defaults to origin)
    :return: Path to output file
    '''
    print(f'\tTransforming {inputFile} -> {outputFile} (preserving PIDs)...')
    
    # Determine file format
    lowerName = inputFile.lower()
    isGz = lowerName.endswith('.gz')
    if '.obj' in lowerName:
        fmt = 'obj'
    elif '.stl' in lowerName:
        fmt = 'stl'
    else:
        sys.exit('ERROR! Input file must be .obj, .stl, .obj.gz, or .stl.gz')
    
    # Open input (text or binary)
    def openInput(mode='rt'):
        if isGz:
            return gzip.open(inputFile, mode)
        else:
            return open(inputFile, mode)
    
    # Ensure output directory exists
    def ensureOutputDir():
        outDir = os.path.dirname(outputFile)
        if outDir and not os.path.exists(outDir):
            os.makedirs(outDir)
    
    outIsGz = outputFile.lower().endswith('.gz')
    
    # precompute axis-angle rotation matrix + pivot for rotation_rvec
    rvecMatrix = None
    if rotation_rvec is not None and float(np.linalg.norm(rotation_rvec)) > 1e-15:
        rvecMatrix = _rodrigues_matrix_from_rvec(rotation_rvec)
    pivotArr = np.array(pivot, dtype=np.float64) if pivot is not None else None
    
    # ---- Build transformation function for a vertex ----
    def transformVertex(v, idx=None):
        v = np.array(v, dtype=np.float64)
        if morphing_dict is not None and idx is not None and idx in morphing_dict:
            v = v + np.array(morphing_dict[idx])
        if scale is not None:
            if isinstance(scale, (int, float)):
                v = v * scale
            else:
                v = v * np.array(scale)
        if rvecMatrix is not None:
            if pivotArr is not None:
                v = rvecMatrix @ (v - pivotArr) + pivotArr
            else:
                v = rvecMatrix @ v
        if rotation is not None:
            rx = rotation.get('x', 0)
            ry = rotation.get('y', 0)
            rz = rotation.get('z', 0)
            if rx != 0:
                v = x_rotation(v, math.radians(rx))
            if ry != 0:
                v = y_rotation(v, math.radians(ry))
            if rz != 0:
                v = z_rotation(v, math.radians(rz))
        if translation is not None:
            v = v + np.array(translation)
        return v
    
    # ---- OBJ handling ----
    if fmt == 'obj':
        with openInput('rt') as f:
            lines = f.readlines()
        
        outLines = []
        vertexIdx = 0
        for line in lines:
            stripped = line.strip()
            if stripped.startswith('v '):
                parts = stripped.split()
                v = [float(parts[1]), float(parts[2]), float(parts[3])]
                vNew = transformVertex(v, vertexIdx)
                outLines.append(f'v {vNew[0]:.6f} {vNew[1]:.6f} {vNew[2]:.6f}\n')
                vertexIdx += 1
            else:
                # Preserve everything else: g, o, f, vn, vt, comments, mtllib, usemtl
                outLines.append(line if line.endswith('\n') else line + '\n')
        
        ensureOutputDir()
        if outIsGz:
            with gzip.open(outputFile, 'wt') as f:
                f.writelines(outLines)
        else:
            with open(outputFile, 'wt') as f:
                f.writelines(outLines)
        
        print(f'\t\tTransformed {vertexIdx} vertices, preserved OBJ groups.')
        return outputFile
    
    # ---- STL handling ----
    if fmt == 'stl':
        # Detect ASCII vs binary
        with openInput('rb') as f:
            head = f.read(512)
        try:
            headStr = head.decode('utf-8', errors='ignore').lstrip().lower()
            asciiLikely = headStr.startswith('solid') and 'facet' in headStr
        except Exception:
            asciiLikely = False
        
        if asciiLikely:
            # ASCII STL — preserve solid/endsolid names and structure
            with openInput('rt') as f:
                lines = f.readlines()
            
            outLines = []
            triBuffer = []  # collect 3 vertices for current facet
            
            for line in lines:
                stripped = line.strip()
                low = stripped.lower()
                
                if low.startswith('solid') or low.startswith('endsolid') or \
                   low.startswith('outer loop') or low.startswith('endloop') or \
                   low.startswith('endfacet'):
                    outLines.append(line if line.endswith('\n') else line + '\n')
                    if low.startswith('endfacet'):
                        triBuffer = []
                
                elif low.startswith('facet normal'):
                    # Hold a placeholder; will be replaced after the 3 vertices
                    triBuffer = []
                    outLines.append(('__FACET_NORMAL_PLACEHOLDER__',))
                
                elif low.startswith('vertex'):
                    parts = stripped.split()
                    v = [float(parts[1]), float(parts[2]), float(parts[3])]
                    vNew = transformVertex(v)
                    triBuffer.append(vNew)
                    outLines.append(f'      vertex {vNew[0]:.6e} {vNew[1]:.6e} {vNew[2]:.6e}\n')
                    
                    if len(triBuffer) == 3:
                        edge1 = triBuffer[1] - triBuffer[0]
                        edge2 = triBuffer[2] - triBuffer[0]
                        n = np.cross(edge1, edge2)
                        nLen = np.linalg.norm(n)
                        if nLen > 0:
                            n = n / nLen
                        else:
                            n = np.array([0.0, 0.0, 0.0])
                        for i in range(len(outLines) - 1, -1, -1):
                            if isinstance(outLines[i], tuple) and outLines[i][0] == '__FACET_NORMAL_PLACEHOLDER__':
                                outLines[i] = f'  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n'
                                break
                else:
                    outLines.append(line if line.endswith('\n') else line + '\n')
            
            # Safety: replace any unfilled placeholders
            outLines = [l if not isinstance(l, tuple)
                        else '  facet normal 0.000000e+00 0.000000e+00 0.000000e+00\n'
                        for l in outLines]
            
            ensureOutputDir()
            if outIsGz:
                with gzip.open(outputFile, 'wt') as f:
                    f.writelines(outLines)
            else:
                with open(outputFile, 'wt') as f:
                    f.writelines(outLines)
            
            print('\t\tTransformed ASCII STL, preserved solid names.')
            return outputFile
        
        else:
            # binary STL - keep the 80-byte header (holds the solid name)
            with openInput('rb') as f:
                data = f.read()
            
            header = data[:80]
            numTri = struct.unpack('<I', data[80:84])[0]
            offset = 84
            
            outBuf = bytearray()
            outBuf += header
            outBuf += struct.pack('<I', numTri)
            
            for _ in range(numTri):
                # skip stored normal — recompute
                offset += 12
                v1 = np.array(struct.unpack('<fff', data[offset:offset+12])); offset += 12
                v2 = np.array(struct.unpack('<fff', data[offset:offset+12])); offset += 12
                v3 = np.array(struct.unpack('<fff', data[offset:offset+12])); offset += 12
                attr = data[offset:offset+2]; offset += 2
                
                v1 = transformVertex(v1)
                v2 = transformVertex(v2)
                v3 = transformVertex(v3)
                
                n = np.cross(v2 - v1, v3 - v1)
                nLen = np.linalg.norm(n)
                if nLen > 0:
                    n = n / nLen
                else:
                    n = np.array([0.0, 0.0, 0.0])
                
                outBuf += struct.pack('<fff', *n)
                outBuf += struct.pack('<fff', *v1)
                outBuf += struct.pack('<fff', *v2)
                outBuf += struct.pack('<fff', *v3)
                outBuf += attr
            
            ensureOutputDir()
            if outIsGz:
                with gzip.open(outputFile, 'wb') as f:
                    f.write(outBuf)
            else:
                with open(outputFile, 'wb') as f:
                    f.write(outBuf)
            
            print(f'\t\tTransformed binary STL ({numTri} triangles), preserved header/solid name.')
            return outputFile


def _open_text_file_maybe_gz(path, mode='rt', compresslevel=1):
    if path.lower().endswith('.gz'):
        # compresslevel=1 (fastest) for writes - these gz outputs are intermediate
        # meshing inputs so write speed matters more than size. reads ignore it
        if 'w' in mode or 'a' in mode or 'x' in mode:
            return gzip.open(path, mode, compresslevel=compresslevel)
        return gzip.open(path, mode)
    return open(path, mode)


def _is_ascii_stl_path(inputFile):
    with _open_text_file_maybe_gz(inputFile, 'rb') as f:
        head = f.read(1024)
    headStr = head.decode('utf-8', errors='ignore').lstrip().lower()
    return headStr.startswith('solid') and ('facet' in headStr or 'endsolid' in headStr)


def _build_pid_regex_mapping(pidNames, pidTransformDict):
    pidToConfig = {}
    unmatched = []

    for pid in pidNames:
        matches = []
        for pattern, cfg in pidTransformDict.items():
            if re.search(pattern, pid):
                matches.append((pattern, cfg))

        if len(matches) == 0:
            unmatched.append(pid)
        elif len(matches) > 1:
            patterns = [m[0] for m in matches]
            raise ValueError(f'PID "{pid}" matches multiple regex patterns: {patterns}')
        else:
            pidToConfig[pid] = matches[0][1]

    return pidToConfig, unmatched


def _canonical_transform_signature(cfg):
    # stable signature to detect conflicting transforms on shared OBJ vertices
    if cfg is None:
        return json.dumps({}, sort_keys=True)
    return json.dumps(cfg, sort_keys=True)


def _rodrigues_matrix_from_rvec(rvec):
    rv = np.array(rvec, dtype=np.float64)
    theta = np.linalg.norm(rv)
    if theta < 1e-15:
        return np.eye(3)
    k = rv / theta
    kx, ky, kz = k
    K = np.array(
        [[0.0, -kz, ky], [kz, 0.0, -kx], [-ky, kx, 0.0]],
        dtype=np.float64,
    )
    return np.eye(3) + math.sin(theta) * K + (1.0 - math.cos(theta)) * (K @ K)


def _build_vertex_transformer(cfg, pivot):
    if cfg is None:
        cfg = {}

    translation = cfg.get('translation', [0.0, 0.0, 0.0])
    rotation = cfg.get('rotation', {'x': 0.0, 'y': 0.0, 'z': 0.0})
    rotationRvec = cfg.get('rotation_rvec', None)
    morphing_dict = cfg.get('morphing_dict', None)
    morph = cfg.get('morph', None)
    
    # pivot from config if present, else the parameter (per-regex override)
    pivot = cfg.get('pivot', pivot)
    # a None pivot must default to origin - np.array(None) is NaN and corrupts the block
    if pivot is None:
        pivot = [0.0, 0.0, 0.0]

    tx, ty, tz = translation
    rx = float(rotation.get('x', 0.0))
    ry = float(rotation.get('y', 0.0))
    rz = float(rotation.get('z', 0.0))

    p = np.array(pivot, dtype=np.float64)
    d = np.array([tx, ty, tz], dtype=np.float64)

    Rr = None
    if rotationRvec is not None:
        Rr = _rodrigues_matrix_from_rvec(rotationRvec)

    axisToIdx = {'x': 0, 'y': 1, 'z': 2}

    def transform_vertex(v, global_idx=None):
        vv = np.array(v, dtype=np.float64)

        if morph is not None:
            mode = morph.get('mode', '').lower()

            if mode == 'displacement':
                disp = np.array(morph.get('vector', [0.0, 0.0, 0.0]), dtype=np.float64)
                vv = vv + disp

            elif mode == 'axis_scale':
                axis = morph.get('axis', 'z').lower()
                if axis not in axisToIdx:
                    raise ValueError(f'Invalid morph axis: {axis}. Expected one of x/y/z.')

                factor = float(morph.get('factor', 1.0))
                origin = np.array(morph.get('origin', pivot), dtype=np.float64)
                limits = morph.get('limits', None)
                axisIdx = axisToIdx[axis]

                inRange = True
                if limits is not None:
                    minVal = limits.get('min', None)
                    maxVal = limits.get('max', None)
                    coord = vv[axisIdx]
                    if minVal is not None and coord < float(minVal):
                        inRange = False
                    if maxVal is not None and coord > float(maxVal):
                        inRange = False

                if inRange:
                    vv[axisIdx] = origin[axisIdx] + (vv[axisIdx] - origin[axisIdx]) * factor

            elif mode == 'vector_stroke':
                axisVector = np.array(morph.get('vector', [0.0, 0.0, 1.0]), dtype=np.float64)
                vnorm = np.linalg.norm(axisVector)
                if vnorm == 0.0:
                    raise ValueError('vector_stroke requires a non-zero morph vector.')
                u = axisVector / vnorm

                origin = np.array(morph.get('origin', pivot), dtype=np.float64)
                distance = float(morph.get('distance', 0.0))
                # positive distance compresses toward origin, negative extends away
                targetLength = float(morph.get('target_length', 0.0))
                limits = morph.get('limits', None)

                rel = vv - origin
                s = float(np.dot(rel, u))
                radial = rel - s * u

                inRange = True
                if limits is not None:
                    minVal = limits.get('min', None)
                    maxVal = limits.get('max', None)
                    if minVal is not None and s < float(minVal):
                        inRange = False
                    if maxVal is not None and s > float(maxVal):
                        inRange = False

                if inRange:
                    # Compress by explicit distance along the selected axis segment.
                    newS = s - distance
                    if targetLength > 0.0:
                        newS = max(0.0, min(targetLength, newS))
                    vv = origin + radial + newS * u

            elif mode == 'vector_scale':
                axisVector = np.array(morph.get('vector', [0.0, 0.0, 1.0]), dtype=np.float64)
                vnorm = np.linalg.norm(axisVector)
                if vnorm == 0.0:
                    raise ValueError('vector_scale requires a non-zero morph vector.')
                u = axisVector / vnorm

                origin = np.array(morph.get('origin', pivot), dtype=np.float64)
                factor = float(morph.get('factor', 1.0))
                limits = morph.get('limits', None)

                rel = vv - origin
                s = float(np.dot(rel, u))
                radial = rel - s * u

                inRange = True
                if limits is not None:
                    minVal = limits.get('min', None)
                    maxVal = limits.get('max', None)
                    if minVal is not None and s < float(minVal):
                        inRange = False
                    if maxVal is not None and s > float(maxVal):
                        inRange = False

                if inRange:
                    vv = origin + radial + (s * factor) * u

            elif mode != '':
                raise ValueError(
                    f'Unsupported morph mode: {mode}. Supported modes: displacement, axis_scale, vector_stroke, vector_scale.'
                )

        if morphing_dict is not None and global_idx is not None and global_idx in morphing_dict:
            vv = vv + np.array(morphing_dict[global_idx], dtype=np.float64)

        vv = vv - p
        if Rr is not None:
            vv = Rr @ vv
        if rx != 0.0:
            vv = x_rotation(vv, math.radians(rx))
        if ry != 0.0:
            vv = y_rotation(vv, math.radians(ry))
        if rz != 0.0:
            vv = z_rotation(vv, math.radians(rz))
        vv = vv + p
        vv = vv + d
        return vv

    return transform_vertex


def _build_block_rigid_transform(cfg, pivot):
    '''
    Build the (R, pivot, translation) for a pure rigid transform (rotation about pivot + translation),
    matching the per-vertex composition order in _build_vertex_transformer:
        v' = Rz @ Ry @ Rx @ Rrvec @ (v - pivot) + pivot + translation.
    Returns (R, p, d) as numpy arrays. Only valid when cfg has no morph / morphing_dict.
    '''
    if cfg is None:
        cfg = {}

    translation = cfg.get('translation', [0.0, 0.0, 0.0])
    rotation = cfg.get('rotation', {'x': 0.0, 'y': 0.0, 'z': 0.0})
    rotationRvec = cfg.get('rotation_rvec', None)
    pivot = cfg.get('pivot', pivot)
    if pivot is None:
        pivot = [0.0, 0.0, 0.0]

    rx = math.radians(float(rotation.get('x', 0.0)))
    ry = math.radians(float(rotation.get('y', 0.0)))
    rz = math.radians(float(rotation.get('z', 0.0)))

    R = np.eye(3)
    if rotationRvec is not None:
        R = _rodrigues_matrix_from_rvec(rotationRvec) @ R
    if rx != 0.0:
        cx, sx = math.cos(rx), math.sin(rx)
        R = np.array([[1.0, 0.0, 0.0], [0.0, cx, -sx], [0.0, sx, cx]]) @ R
    if ry != 0.0:
        cy, sy = math.cos(ry), math.sin(ry)
        R = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]]) @ R
    if rz != 0.0:
        cz, sz = math.cos(rz), math.sin(rz)
        R = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]]) @ R

    p = np.array(pivot, dtype=np.float64)
    tx, ty, tz = translation
    d = np.array([tx, ty, tz], dtype=np.float64)
    return R, p, d


def _transform_vertex_block(cfg, pivot, block, globalIndices):
    '''
    Transform an (M, 3) block of vertices that all share the same transform config.

    For pure rigid transforms (no morph / morphing_dict) the whole block is transformed with a single
    vectorized matrix multiply. When per-vertex morphing is present, falls back to the per-vertex
    transformer (built once for the block) so morph/morphing_dict semantics are preserved exactly.
    '''
    hasMorph = cfg is not None and (cfg.get('morph', None) is not None or cfg.get('morphing_dict', None) is not None)
    if not hasMorph:
        R, p, d = _build_block_rigid_transform(cfg, pivot)
        return (block - p) @ R.T + p + d

    fn = _build_vertex_transformer(cfg, pivot)
    out = np.empty_like(block)
    for k in range(block.shape[0]):
        out[k] = fn(block[k], int(globalIndices[k]))
    return out


def _write_pid_transform_log(logBasePath, logData):
    jsonPath = f'{logBasePath}.json'
    csvPath = f'{logBasePath}.csv'

    outDir = os.path.dirname(jsonPath)
    if outDir and not os.path.exists(outDir):
        os.makedirs(outDir)

    with open(jsonPath, 'w') as jf:
        json.dump(logData, jf, indent=2)

    with open(csvPath, 'w', newline='') as cf:
        writer = csv.DictWriter(
            cf,
            fieldnames=['pid', 'regex', 'translation', 'rotation', 'morph', 'vertex_count', 'warnings']
        )
        writer.writeheader()
        for row in logData['pid_rows']:
            writer.writerow(row)


def transformGeometryByPIDRegex(
    inputFile,
    outputFile,
    pidTransformDict,
    pivot,
    requireAllPIDsMatched=True,
    logBasePath=None
):
    '''
    Unpack PID groups from OBJ or ASCII STL, transform by PID regex, and repack to the same
    format/compression while preserving line ordering and indexing.

    PID mapping:
    - OBJ: current `g` or `o` name (most recent encountered name).
    - ASCII STL: `solid <name>`.

    Regex behavior is case-sensitive by design.

    :param inputFile: Input geometry path (.obj, .obj.gz, .stl, .stl.gz)
    :param outputFile: Output geometry path (must keep same format and compression as input)
    :param pidTransformDict: Dict[regexPattern] -> {
        'translation': [dx, dy, dz],
        'rotation': {'x': deg, 'y': deg, 'z': deg},
        'rotation_rvec': [rx, ry, rz] axis-angle rotation vector in radians (optional),
        'morphing_dict': {globalVertexIndex: [dx, dy, dz]} (optional),
        'morph': {
            'mode': 'axis_scale' or 'displacement' or 'vector_stroke' or 'vector_scale',
            # axis_scale:
            'axis': 'x'/'y'/'z',
            'factor': 1.0,
            'origin': [x, y, z],
            'limits': {'min': float, 'max': float} (optional),
            # displacement:
            'vector': [dx, dy, dz],
            # vector_stroke:
            'vector': [vx, vy, vz],
            'distance': float,
            'origin': [x, y, z],
            'target_length': float (optional),
            # vector_scale:
            'vector': [vx, vy, vz],
            'factor': float,
            'origin': [x, y, z],
            'limits': {'min': float, 'max': float} (optional)
        } (optional)
    }
    :param pivot: [x, y, z] fixed pivot for extrinsic XYZ rotations
    :param requireAllPIDsMatched: If True, error when any PID is not matched by regex map
    :param logBasePath: Base path for JSON/CSV log output (without extension)
    :return: Dict with output path and warning metadata
    '''
    if not isinstance(pidTransformDict, dict) or len(pidTransformDict) == 0:
        raise ValueError('pidTransformDict must be a non-empty dictionary of regex -> transform config.')

    if outputFile is None:
        raise ValueError('outputFile must be provided.')

    inLower = inputFile.lower()
    outLower = outputFile.lower()

    inFmt = 'obj' if '.obj' in inLower else ('stl' if '.stl' in inLower else None)
    outFmt = 'obj' if '.obj' in outLower else ('stl' if '.stl' in outLower else None)
    if inFmt is None or outFmt is None:
        raise ValueError('Only OBJ/STL files are supported.')

    inGz = inLower.endswith('.gz')
    outGz = outLower.endswith('.gz')
    if inFmt != outFmt or inGz != outGz:
        raise ValueError('outputFile must preserve input format and compression exactly.')

    if inFmt == 'stl' and not _is_ascii_stl_path(inputFile):
        raise ValueError('Only ASCII STL is supported for PID-preserving regex transforms.')

    warnings = []
    timestamp = datetime.utcnow().isoformat() + 'Z'

    if inFmt == 'obj':
        with _open_text_file_maybe_gz(inputFile, 'rt') as f:
            lines = f.readlines()

        vertices = []
        currentPid = '__DEFAULT__'
        pidVertexRefs = {}
        pidVertexRefs[currentPid] = set()

        # hot loop over millions of lines: work on the raw line and cache the current
        # group's vertex-ref set to skip per-line strip() and setdefault() lookups
        currentRefSet = pidVertexRefs[currentPid]
        for line in lines:
            if line.startswith('v '):
                parts = line.split()
                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif line.startswith('f '):
                nVerts = len(vertices)
                for ref in line.split()[1:]:
                    if '/' in ref:
                        ref = ref.split('/', 1)[0]
                    idx = int(ref)
                    vIdx = idx - 1 if idx > 0 else nVerts + idx
                    if vIdx >= 0:
                        currentRefSet.add(vIdx)
            elif line.startswith('g ') or line.startswith('o '):
                toks = line.split(maxsplit=1)
                currentPid = toks[1].strip() if len(toks) > 1 else '__EMPTY__'
                currentRefSet = pidVertexRefs.get(currentPid)
                if currentRefSet is None:
                    currentRefSet = set()
                    pidVertexRefs[currentPid] = currentRefSet

        pidNames = [pid for pid, refs in pidVertexRefs.items() if len(refs) > 0]
        pidToCfg, unmatched = _build_pid_regex_mapping(pidNames, pidTransformDict)
        if requireAllPIDsMatched and len(unmatched) > 0:
            raise ValueError(f'Unmatched OBJ PIDs: {unmatched}')
        if len(unmatched) > 0:
            warnings.append(f'Unmatched OBJ PIDs: {unmatched}')

        vertexToSignature = {}
        vertexToCfg = {}
        vertexToPid = {}
        for pid in pidNames:
            cfg = pidToCfg.get(pid, None)
            # unmatched PIDs have no transform - leave them identity. a signature would
            # route them through the rigid transform with pivot=None -> NaN, corrupting static parts
            if cfg is None:
                continue
            sig = _canonical_transform_signature(cfg)
            for vIdx in pidVertexRefs[pid]:
                if vIdx in vertexToSignature and vertexToSignature[vIdx] != sig:
                    prevPid = vertexToPid.get(vIdx, '<unknown>')
                    raise ValueError(
                        f'Vertex {vIdx} is shared by multiple PIDs with different transforms. '
                        'Cannot preserve original indexing without splitting vertices. '
                        f'Offending PID pair: "{prevPid}" [transform={vertexToSignature[vIdx]}] '
                        f'vs "{pid}" [transform={sig}].'
                    )
                vertexToSignature[vIdx] = sig
                vertexToCfg[vIdx] = cfg
                vertexToPid[vIdx] = pid

        # vectorized transform: group vertices by signature and apply each transform to the
        # whole block at once (one numpy matmul per group). unreferenced vertices keep identity
        transformedArr = np.array(vertices, dtype=np.float64) if len(vertices) > 0 else np.zeros((0, 3), dtype=np.float64)

        sigToIndices = {}
        for vIdx, sig in vertexToSignature.items():
            sigToIndices.setdefault(sig, []).append(vIdx)

        for sig, idxList in sigToIndices.items():
            idxArr = np.array(idxList, dtype=np.int64)
            cfg = vertexToCfg[idxList[0]]
            transformedArr[idxArr] = _transform_vertex_block(cfg, pivot, transformedArr[idxArr], idxArr)

        outLines = []
        vCounter = 0
        for line in lines:
            if line.startswith('v '):
                vv = transformedArr[vCounter]
                outLines.append(f'v {vv[0]:.6f} {vv[1]:.6f} {vv[2]:.6f}\n')
                vCounter += 1
            else:
                outLines.append(line if line.endswith('\n') else line + '\n')

        outDir = os.path.dirname(outputFile)
        if outDir and not os.path.exists(outDir):
            os.makedirs(outDir)
        with _open_text_file_maybe_gz(outputFile, 'wt') as f:
            f.writelines(outLines)

        pidRows = []
        for pid in pidNames:
            matchedPattern = None
            for pattern in pidTransformDict.keys():
                if re.search(pattern, pid):
                    matchedPattern = pattern
                    break
            cfg = pidToCfg.get(pid, {})
            pidRows.append({
                'pid': pid,
                'regex': matchedPattern,
                'translation': cfg.get('translation', [0.0, 0.0, 0.0]) if cfg else [0.0, 0.0, 0.0],
                'rotation': cfg.get('rotation', {'x': 0.0, 'y': 0.0, 'z': 0.0}) if cfg else {'x': 0.0, 'y': 0.0, 'z': 0.0},
                'morph': cfg.get('morph', None) if cfg else None,
                'vertex_count': len(pidVertexRefs.get(pid, set())),
                'warnings': '; '.join(warnings)
            })

    else:
        with _open_text_file_maybe_gz(inputFile, 'rt') as f:
            lines = f.readlines()

        pidNames = []
        currentPid = None
        for line in lines:
            stripped = line.strip()
            low = stripped.lower()
            if low.startswith('solid'):
                toks = stripped.split(maxsplit=1)
                currentPid = toks[1] if len(toks) > 1 else '__EMPTY__'
                pidNames.append(currentPid)

        pidToCfg, unmatched = _build_pid_regex_mapping(pidNames, pidTransformDict)
        if requireAllPIDsMatched and len(unmatched) > 0:
            raise ValueError(f'Unmatched STL PIDs: {unmatched}')
        if len(unmatched) > 0:
            warnings.append(f'Unmatched STL PIDs: {unmatched}')

        outLines = []
        currentPid = None
        triBuffer = []
        facetPlaceholderIdx = None
        globalVertexIdx = 0
        pidVertexCounter = {pid: 0 for pid in pidNames}
        # cache one transformer per PID instead of rebuilding the closure per vertex
        pidTransformerCache = {}

        for line in lines:
            stripped = line.strip()
            low = stripped.lower()

            if low.startswith('solid'):
                toks = stripped.split(maxsplit=1)
                currentPid = toks[1] if len(toks) > 1 else '__EMPTY__'
                outLines.append(line if line.endswith('\n') else line + '\n')
                continue

            if low.startswith('facet normal'):
                triBuffer = []
                facetPlaceholderIdx = len(outLines)
                outLines.append('__FACET_NORMAL_PLACEHOLDER__\n')
                continue

            if low.startswith('vertex'):
                if currentPid is None:
                    raise ValueError('Encountered STL vertex line outside of a solid block.')
                parts = stripped.split()
                v = np.array([float(parts[1]), float(parts[2]), float(parts[3])], dtype=np.float64)
                fn = pidTransformerCache.get(currentPid)
                if fn is None:
                    fn = _build_vertex_transformer(pidToCfg.get(currentPid, None), pivot)
                    pidTransformerCache[currentPid] = fn
                vNew = fn(v, globalVertexIdx)
                globalVertexIdx += 1
                pidVertexCounter[currentPid] = pidVertexCounter.get(currentPid, 0) + 1
                triBuffer.append(vNew)
                outLines.append(f'      vertex {vNew[0]:.6e} {vNew[1]:.6e} {vNew[2]:.6e}\n')

                if len(triBuffer) == 3 and facetPlaceholderIdx is not None:
                    edge1 = triBuffer[1] - triBuffer[0]
                    edge2 = triBuffer[2] - triBuffer[0]
                    n = np.cross(edge1, edge2)
                    nLen = np.linalg.norm(n)
                    if nLen > 0:
                        n = n / nLen
                    else:
                        n = np.array([0.0, 0.0, 0.0])
                    outLines[facetPlaceholderIdx] = (
                        f'  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n'
                    )
                continue

            outLines.append(line if line.endswith('\n') else line + '\n')

        outDir = os.path.dirname(outputFile)
        if outDir and not os.path.exists(outDir):
            os.makedirs(outDir)
        with _open_text_file_maybe_gz(outputFile, 'wt') as f:
            f.writelines(outLines)

        pidRows = []
        for pid in pidNames:
            matchedPattern = None
            for pattern in pidTransformDict.keys():
                if re.search(pattern, pid):
                    matchedPattern = pattern
                    break
            cfg = pidToCfg.get(pid, {})
            pidRows.append({
                'pid': pid,
                'regex': matchedPattern,
                'translation': cfg.get('translation', [0.0, 0.0, 0.0]) if cfg else [0.0, 0.0, 0.0],
                'rotation': cfg.get('rotation', {'x': 0.0, 'y': 0.0, 'z': 0.0}) if cfg else {'x': 0.0, 'y': 0.0, 'z': 0.0},
                'morph': cfg.get('morph', None) if cfg else None,
                'vertex_count': pidVertexCounter.get(pid, 0),
                'warnings': '; '.join(warnings)
            })

    if logBasePath is None:
        logBasePath = outputFile + '.pid_transform_log'

    matchedPidCount = sum(1 for row in pidRows if row.get('regex') is not None)

    logData = {
        'timestamp': timestamp,
        'input_file': inputFile,
        'output_file': outputFile,
        'pivot': list(pivot) if pivot is not None else None,
        'warnings': warnings,
        'matched_pid_count': matchedPidCount,
        'pid_rows': pidRows
    }
    _write_pid_transform_log(logBasePath, logData)

    return {
        'output_file': outputFile,
        'warnings': warnings,
        'matched_pid_count': matchedPidCount,
        'log_json': f'{logBasePath}.json',
        'log_csv': f'{logBasePath}.csv'
    }

def transformGeometry(inputFile, outputFile=None, rotation=None, translation=None, scale=None, morphing_dict=None):
    '''
    Transform geometry files (OBJ or STL, with or without .gz compression).
    
    Supports:
    - Rotation: around x, y, z axes (in degrees)
    - Translation: offset in x, y, z directions
    - Scaling: uniform or non-uniform scaling
    - Morphing: vertex-specific displacements based on a morphing dictionary
    
    :param inputFile: Path to input geometry file (OBJ or STL, can be .gz compressed)
    :param outputFile: Path to output geometry file. If None, returns vertices and faces only
    :param rotation: Dict with keys 'x', 'y', 'z' for rotation angles in degrees, e.g. {'x': 0, 'y': 0, 'z': 0}
    :param translation: List/array [dx, dy, dz] for translation offsets
    :param scale: Scalar for uniform scaling, or list/array [sx, sy, sz] for non-uniform scaling
    :param morphing_dict: Dict mapping vertex indices to displacement vectors, e.g. {0: [0.1, 0, 0], 5: [0, 0.2, 0]}
    :return: Tuple of (vertices, faces) - transformed geometry
    '''
    
    # Load the geometry file
    print(f'\tLoading geometry from {inputFile}...')
    vertices, faces = readGeomFile(inputFile)
    vertices = np.array(vertices, dtype=np.float64)
    
    # Apply morphing first (vertex-specific displacements)
    if morphing_dict is not None:
        print('\tApplying morphing deformations...')
        for vertex_idx, displacement in morphing_dict.items():
            if vertex_idx < len(vertices):
                vertices[vertex_idx] += np.array(displacement)
            else:
                print(f'\t\tWarning: Vertex index {vertex_idx} out of range')
    
    # Apply scaling
    if scale is not None:
        print('\tApplying scaling...')
        if isinstance(scale, (int, float)):
            # Uniform scaling
            vertices *= scale
        else:
            # Non-uniform scaling [sx, sy, sz]
            vertices *= np.array(scale)
    
    # Apply rotation (order: x, y, z)
    if rotation is not None:
        print('\tApplying rotations...')
        rot_x = rotation.get('x', 0)
        rot_y = rotation.get('y', 0)
        rot_z = rotation.get('z', 0)
        
        if rot_x != 0:
            rot_x_rad = math.radians(rot_x)
            for i in range(len(vertices)):
                vertices[i] = x_rotation(vertices[i], rot_x_rad)
        
        if rot_y != 0:
            rot_y_rad = math.radians(rot_y)
            for i in range(len(vertices)):
                vertices[i] = y_rotation(vertices[i], rot_y_rad)
        
        if rot_z != 0:
            rot_z_rad = math.radians(rot_z)
            for i in range(len(vertices)):
                vertices[i] = z_rotation(vertices[i], rot_z_rad)
    
    # Apply translation
    if translation is not None:
        print('\tApplying translation...')
        vertices += np.array(translation)
    
    # Write output file if specified
    if outputFile is not None:
        print(f'\tWriting transformed geometry to {outputFile}...')
        writeGeometryFile(outputFile, vertices, faces)
    
    return vertices, faces
    
def runRHCaseSetup(rideHeights, fullCaseSetupDict):
    # Get the path to caseSetup.py
    execDir = os.path.dirname(os.path.realpath(__file__))
    caseSetupScript = os.path.join(execDir, 'caseSetup.py')
    
    # Get venv python path
    venv_python = os.path.join(execDir, '.venv', 'bin', 'python3')
    if not os.path.exists(venv_python):
        venv_python = 'python3'  # Fallback to system python
    
    for idx, row in rideHeights.iterrows():
        point = (int(row['point']))
        childName = os.path.basename(os.getcwd()) + '_' + str(point)
        childCasePath = os.path.join(os.getcwd(), childName)
        print('\t\tRunning caseSetup for case: %s' % (childName))
        
        # Run caseSetup in the child directory
        cmd = [venv_python, caseSetupScript, '-s', 'otr', '--rideHeightMode']
        result = subprocess.run(cmd, cwd=childCasePath)
        
        if result.returncode != 0:
            print(f'\t\t\tWARNING: caseSetup failed for {childName}')
        else:
            print(f'\t\t\tcaseSetup completed successfully for {childName}')
    
# def runCaseSetup(caseSetupPath, setupType='otr', venv_python=None):
#     """Run caseSetup.py in a subprocess."""
#     if venv_python is None:
#         venv_python = 'python3'
    
#     cmd = [venv_python, 'caseSetup.py', '-s', setupType, '-rideHeightMode']
#     result = subprocess.run(cmd, cwd=os.path.dirname(caseSetupPath))
#     return result.returncode

def writeToRHCaseSetup(writeCaseSetupDict,newCasePath):
    writeConfig = configparser.ConfigParser()
    writeConfig.optionxform = str
    for module in writeCaseSetupDict.keys():
        writeConfig.add_section(module)
        for key in writeCaseSetupDict[module].keys():
            #print(" ".join(list(writeCaseSetupDict[module][key])))
            try:
                writeConfig.set(module,key," ".join(list(writeCaseSetupDict[module][key])))
            except:
                writeConfig.set(module,key,str(writeCaseSetupDict[module][key]))

    with open(newCasePath,'w') as caseSetupFile:
                writeConfig.write(caseSetupFile)
    
    if updateCaseSetupFlag == True:
        sys.exit('\nWARNING: caseSetup has been updated with default values, please check caseSetup and rerun.')
                

def calculateRollAngles(dl,dr,bw):
    alpha = math.degrees(math.asin((dr-dl)/bw))
    d_center = (dl+dr)/2
    return round(alpha,3), round(d_center,4)

def calculatePitchAngles(dzf,dzr,bl,REFCOR):
    theta = math.degrees(math.asin((dzf-dzr)/bl))
    dz_wb_center = ((dzf+dzr)/2)-float(REFCOR[2])

    return round(theta,3), round(dz_wb_center,4)




def calculateWheelMovements(z_fl, z_fr, z_rl, z_rr, pitch_deg, roll_deg, heave):
    '''
    Calculate wheel movements in tunnel coordinate system.
    
    In the CFD coordinate system, the tunnel (and wheels) move around the car body.
    Since the tunnel/wheels move around the car, the wheel movements are the opposite
    of the ride height values.
    
    REFCOR: Reference center of rotation at the wheelbase center (y=0) projected to ground
    
    Convention:
    - Positive pitch: front goes up, rear goes down
    - Positive roll: right side goes up, left side goes down  
    - Heave: vertical translation of entire car
    
    The wheel vertical positions (z_fl, z_fr, z_rl, z_rr) from the ride height file
    represent car geometry. In the tunnel frame, wheel movements are opposite.
    
    :param z_fl, z_fr, z_rl, z_rr: vertical displacements of each wheel (car frame)
    :param pitch_deg: pitch angle in degrees (positive = front up)
    :param roll_deg: roll angle in degrees (positive = right up)
    :param heave: vertical heave displacement at reference center
    :param bw: wheel track [front_track, rear_track]
    :param bl: wheelbase length
    :param REFCOR: reference corner coordinates [x, y, z] - wheelbase center at y=0 on ground
    :return: wheel movements for FL, FR, RL, RR in tunnel frame
    '''
    
    # Wheel movements are opposite direction to ride height (tunnel moves around car)
    wheel_fl = -z_fl
    wheel_fr = -z_fr
    wheel_rl = -z_rl
    wheel_rr = -z_rr

    # calculate tunnel movement due to pitch (rotation around y-axis)
    tunnel_pitch = -pitch_deg
    tunnel_roll = -roll_deg
    # calculate tunnel vertical movement before applying roll and pitch
    tunnel_dz = -heave

    


    
    
    return wheel_fl, wheel_fr, wheel_rl, wheel_rr, tunnel_pitch, tunnel_roll, tunnel_dz




def getBoundingBoxOBJ(geomFile):
    objPath = 'constant/triSurface/%s' % (geomFile)
    
    
    #check that obj file is indeed an obj file
    
    if not geomFile.split('.')[-1].lower() == 'obj' or geomFile.split('.')[-2].lower() == 'obj':
        sys.exit('ERROR! Cannot attempt to import non-OBJ files using getBoundingBoxOBJ function!')
        
    #read in the file from the path given
    if geomFile.split('.')[-1].lower() == 'gz':
        import gzip
        import shutil
        with gzip.open(objPath, 'rt') as objFile:
            #initializing the arrays
            xCoords = []
            yCoords = []
            zCoords = []
            
            for line in objFile:
                if line.startswith('v'):
                    xCoords.append(line.split(' ')[1])
                    yCoords.append(line.split(' ')[2])
                    zCoords.append(line.split(' ')[3])
    else:
        with open(objPath, 'rt') as objFile:
            #initializing the arrays
            xCoords = []
            yCoords = []
            zCoords = []
            
            for line in objFile:
                
                if line.startswith('v'):
                    xCoords.append(line.split(' ')[1])
                    yCoords.append(line.split(' ')[2])
                    zCoords.append(line.split(' ')[3])
    
    xCoords = np.array(xCoords)
    yCoords = np.array(yCoords)
    zCoords = np.array(zCoords)
    
    #getting bounds
    
    bbminX = min(xCoords)
    bbmaxX = max(xCoords)
    bbminY = min(yCoords)
    bbmaxY = max(yCoords)
    bbminZ = min(zCoords)
    bbmaxZ = max(zCoords)
    sys.exit()
    
    return bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ
        
def getBoundingBoxSTL(geomFile):
    
    stlPath = 'constant/triSurface/%s' % (geomFile)
    
    #check that obj file is indeed an obj file
    if not geomFile.split('.')[-1].lower() == 'stl' or geomFile.split('.')[-2].lower() == 'stl':
        sys.exit('ERROR! Cannot attempt to import non-STL files using getBoundingBoxSTL function!')
        
    #read in the file from the path given
    if geomFile.split('.')[-1].lower() == 'gz':
        import gzip
        import shutil
        with gzip.open(stlPath, 'rt') as stlFile:
            #initializing the arrays
            xCoords = []
            yCoords = []
            zCoords = []
            
            for line in stlFile:
                if 'vertex' in line:
                    line.replace('  ',',').replace(' ',',')
                    xCoords.append(line.split(',')[1])
                    yCoords.append(line.split(',')[2])
                    zCoords.append(line.split(',')[3])
    else:
        with open(stlPath, 'rt') as stlFile:
            #initializing the arrays
            xCoords = []
            yCoords = []
            zCoords = []
            
            for line in stlFile:
                if 'vertex' in line:
                    line.replace('  ',',').replace(' ',',')
                    xCoords.append(line.split(',')[1])
                    yCoords.append(line.split(',')[2])
                    zCoords.append(line.split(',')[3])
    
    xCoords = np.array(xCoords)
    yCoords = np.array(yCoords)
    zCoords = np.array(zCoords)
    
    #getting bounds
    
    bbminX = min(xCoords)
    bbmaxX = max(xCoords)
    bbminY = min(yCoords)
    bbmaxY = max(yCoords)
    bbminZ = min(zCoords)
    bbmaxZ = max(zCoords)
    
    sys.exit()
    
    
    return bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ   

def getBoundingBoxPv(geomFile):
    geomPath = 'constant/triSurface/%s' % (geomFile)
    paraview.simple._DisableFirstRenderCameraReset()
    
    if '.obj' in geomFile or '.OBJ' in geomFile:
        geomReader = WavefrontOBJReader(registrationName='geomFile',FileName = geomPath)
    elif '.stl' in geomFile or '.STL' in geomFile:
        geomReader = STLReader(registrationName='geomFile',FileNames = geomPath)
        
    renderView1 = GetActiveViewOrCreate('RenderView')
    geomDisplay = Show(geomReader,renderView1)
    bounds = geomReader.GetDataInformation().GetBounds()
    bbminX = bounds[0]
    bbmaxX = bounds[1]
    bbminY = bounds[2]
    bbmaxY = bounds[3]
    bbminZ = bounds[4]
    bbmaxZ = bounds[5]
    
    del geomDisplay
    del renderView1
    del geomReader
    return bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ

def getRotaCoordinates(bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ):   
    xcenter = (bbminX + bbmaxX)/2
    ycenter = (bbminY + bbmaxY)/2
    zcenter = (bbminZ + bbmaxZ)/2
    
    radius = (bbmaxZ - bbminZ)/2
    return round(xcenter,6),round(ycenter,6),round(zcenter,6),round(radius,6)

def calcRotaVel(inletMag,radius):
    rotaVel = inletMag/radius
    
    return rotaVel

def calcLoadedRadius(xcenter, ycenter, zcenter, fullCaseSetupDict):
    '''
    Effective rolling radius for the rotatingWallVelocity BC = vertical distance from the wheel center
    (axle) down to the GROUND PLANE at the contact patch directly beneath the axle.

    The radius is measured to the ground plane, not to the tyre geometry. The tyre mesh penetrates the
    ground and is cut at the floor by snappyHexMesh, so the actual rolling-wall contact patch lies on the
    ground plane (z_ground), not at the lowest tyre vertex (which sits below the road). The BC sets the
    tyre surface velocity at that contact point to the road speed, so omega = V / R with R = z_axle -
    z_ground. This correctly makes R shrink when the suspension compresses (axle drops, penetration
    increases) and grow when it extends, so the effective rolling radius genuinely differs between
    ride-height points.

    The tunnel floor passes through REFCOR (heave-corrected per ride-height child case) and is rotated
    about REFCOR by roll (about x) then pitch (about y), matching transformBlockmeshPoints. Its normal is
        n = Ry(pitch) @ Rx(roll) @ (0,0,1) = (sin(pitch)cos(roll), -sin(roll), cos(pitch)cos(roll)).
    The ground height under the axle is found by intersecting the vertical line through (xcenter, ycenter)
    with that plane:
        z_ground = REFCOR_z - (n_x*(xcenter-REFCOR_x) + n_y*(ycenter-REFCOR_y)) / n_z
    which reduces to z_ground = REFCOR_z for zero pitch/roll.

    :param xcenter, ycenter, zcenter: wheel center / axle coordinates (m)
    :param fullCaseSetupDict: case setup dict providing BC_SETUP REFCOR, DOMAIN_PITCH, DOMAIN_ROLL
    :return: effective rolling radius (m)
    '''
    refcor = fullCaseSetupDict['BC_SETUP']['REFCOR']
    refX = float(refcor[0])
    refY = float(refcor[1])
    refZ = float(refcor[2])
    pitch = math.radians(float(fullCaseSetupDict['BC_SETUP']['DOMAIN_PITCH'][0]))
    roll = math.radians(float(fullCaseSetupDict['BC_SETUP']['DOMAIN_ROLL'][0]))

    nx = math.sin(pitch) * math.cos(roll)
    ny = -math.sin(roll)
    nz = math.cos(pitch) * math.cos(roll)

    #ground z directly below the axle (vertical line through xcenter, ycenter intersecting the floor plane)
    groundZ = refZ - (nx * (xcenter - refX) + ny * (ycenter - refY)) / nz
    radius = zcenter - groundZ

    return radius

def corneringAxis(fullCaseSetupDict):
    '''
    Vertical (tilted) rotation axis for SRF cornering, identical to the ride-height ground normal so the
    single rotating frame stays aligned with the floor when the domain is pitched/rolled:
        n = Ry(pitch) @ Rx(roll) @ (0,0,1) = (sin(pitch)cos(roll), -sin(roll), cos(pitch)cos(roll))

    :param fullCaseSetupDict: case setup dict providing BC_SETUP DOMAIN_PITCH, DOMAIN_ROLL
    :return: unit numpy axis vector
    '''
    pitch = math.radians(float(fullCaseSetupDict['BC_SETUP']['DOMAIN_PITCH'][0]))
    roll = math.radians(float(fullCaseSetupDict['BC_SETUP']['DOMAIN_ROLL'][0]))
    axis = np.array([math.sin(pitch) * math.cos(roll),
                     -math.sin(roll),
                     math.cos(pitch) * math.cos(roll)])
    return axis / np.linalg.norm(axis)

def corneringFrame(fullCaseSetupDict):
    '''
    Single-rotating-frame (SRF) parameters for a steady curved-path (cornering) simulation. The whole
    domain is solved in a frame rotating at omega about a vertical axis through the corner centre, so the
    car follows a circular path of radius CORNER_RADIUS at speed INLET_MAG:

        |omega| = INLET_MAG / CORNER_RADIUS                     (rad/s)

    Sign convention (car nose points toward the inlet at -x, freestream travels +x; up is +z):
        CORNER_DIR = left  -> centre on the car's left  (-y), yaw about +axis -> positive omega
        CORNER_DIR = right -> centre on the car's right (+y), yaw about +axis -> negative omega

    The corner centre is measured from REFCOR (the reference centre of rotation at the wheelbase centre,
    y=0, projected to ground) with a lateral +/- CORNER_RADIUS offset in y (the dominant lateral term;
    pitch/roll tilt of the lateral offset is neglected as those angles are small). CORNER_CENTER may be
    set explicitly (three floats) to override the derived centre. The rotation axis follows the tilted
    ground normal (see corneringAxis) so SRF stays aligned with ride-height attitude.

    :param fullCaseSetupDict: case setup dict (BC_SETUP, CORNERING_SETUP)
    :return: (omegaSigned [rad/s], axis [unit numpy], centre [numpy x y z])
    '''
    inletMag = float(fullCaseSetupDict['BC_SETUP']['INLET_MAG'][0])
    radius = float(fullCaseSetupDict['CORNERING_SETUP']['CORNER_RADIUS'][0])
    cornerDir = str(fullCaseSetupDict['CORNERING_SETUP']['CORNER_DIR'][0]).strip().lower()

    if radius <= 0:
        sys.exit('ERROR! [CORNERING_SETUP] -> CORNER_RADIUS must be a positive number when RUN_CORNERING is True!')

    omegaMag = inletMag / radius

    if cornerDir == 'left':
        sign = 1.0
        lateral = -1.0
    elif cornerDir == 'right':
        sign = -1.0
        lateral = 1.0
    else:
        sys.exit("ERROR! [CORNERING_SETUP] -> CORNER_DIR must be 'left' or 'right'!")

    omegaSigned = sign * omegaMag
    axis = corneringAxis(fullCaseSetupDict)

    cornerCenter = fullCaseSetupDict['CORNERING_SETUP']['CORNER_CENTER']
    if str(cornerCenter[0]).strip().lower() != 'default':
        try:
            centre = np.array([float(cornerCenter[0]), float(cornerCenter[1]), float(cornerCenter[2])])
        except Exception as error:
            print("\t\tWARNING! [CORNERING_SETUP] -> CORNER_CENTER is not three valid floats, using derived centre instead.")
            print('\t\t%s' % (error))
            cornerCenter = ['default']
    if str(cornerCenter[0]).strip().lower() == 'default':
        refCor = fullCaseSetupDict['BC_SETUP']['REFCOR']
        centre = np.array([float(refCor[0]),
                           float(refCor[1]) + lateral * radius,
                           float(refCor[2])])

    return omegaSigned, axis, centre

def checkCorneringDomain(fullCaseSetupDict):
    '''
    Validity guard for SRF cornering. The rectangular box approximates an annular sector only when the
    domain's lateral half-width is small relative to the turn radius; otherwise inner/outer walls span
    very different path radii and the box is a poor model of the curved tunnel.

        require  CORNER_RADIUS >= CORNER_CLEARANCE_FACTOR * (DOMAIN_SIZE_y / 2)

    Never auto-resizes; on failure reports the minimum admissible radius and exits.

    :param fullCaseSetupDict: case setup dict (BC_SETUP, CORNERING_SETUP)
    :return: minimum admissible radius (m)
    '''
    try:
        radius = float(fullCaseSetupDict['CORNERING_SETUP']['CORNER_RADIUS'][0])
    except (ValueError, IndexError):
        sys.exit('ERROR! [CORNERING_SETUP] -> CORNER_RADIUS must be a positive number when '
                 'RUN_CORNERING is True (set it directly, or provide a per-point corner_radius '
                 'column in the ride height file).')
    factor = float(fullCaseSetupDict['CORNERING_SETUP']['CORNER_CLEARANCE_FACTOR'][0])
    ydom = float(fullCaseSetupDict['BC_SETUP']['DOMAIN_SIZE'][1])
    halfWidth = ydom / 2.0
    rMin = factor * halfWidth
    if radius < rMin:
        sys.exit('ERROR! [CORNERING_SETUP] -> CORNER_RADIUS=%1.3fm is too small for this domain '
                 '(lateral half-width=%1.3fm, clearance factor=%1.2f). Minimum admissible radius is '
                 '%1.3fm. Increase CORNER_RADIUS, reduce DOMAIN_SIZE y, or lower CORNER_CLEARANCE_FACTOR.'
                 % (radius, halfWidth, factor, rMin))
    return rMin

def getBoundingBox(geomFile):
    command = 'surfaceCheck -outputThreshold 0 constant/triSurface/%s' % (geomFile)
    searchType = 'contains'
    searchVar = 'Bounding'
    boundingBox = getBashOutput(command,searchType,searchVar)
    boundingBox = boundingBox[0].replace('Bounding Box : ','')\
                                .replace('(','')\
                                .replace(')','')\
                                .split(' ')
    
    if len(boundingBox) < 6:
        sys.exit('ERROR! %s is not a valid geometry, unable to get bounding box!')
    else:    
        bbminX = float(boundingBox[0])
        bbminY = float(boundingBox[1])
        bbminZ = float(boundingBox[2])
        bbmaxX = float(boundingBox[3])
        bbmaxY = float(boundingBox[4])
        bbmaxZ = float(boundingBox[5])
    
    
        return bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ
    

    

     
def getBashOutput(command,searchType,searchVar):

    if type(command) == list:
        command = command
    elif type(command) == str:
        command = command.split(' ')
    else:
        sys.exit("ERROR: Invalid Type in bashOutput")
        
    bashOutput = sp.check_output(command)
    bashOutput = bashOutput.decode('utf-8')
    bashOutput = bashOutput.split('\n')
    outputLine = []
    for line in bashOutput:
        
        if searchType.lower() == 'start' and line.startswith(searchVar):
            outputLine.append(line)      
        elif searchType.lower() == 'end' and line.endswith(searchVar):
            outputLine.append(line)    
        elif searchType.lower() == 'contains' and searchVar in line:
            outputLine.append(line)    
        else:
            continue   
    return outputLine 
    
def stripExt(string):
    string = string.split('.')
    return string[0]
    
    
def copyTemplateToCase(templatePath,templateDest):
    #check if the source path is valid
    if not os.path.exists(templatePath):
        sys.exit('\n\t\tERROR: %s cannot be found!' % (templatePath))
    try:
        os.system('cp %s %s' % (templatePath,templateDest))
    except Exception as e:
        print('ERROR: Unable to copy template to case!')
        print('\n\t\t\t'+ e )
        sys.exit()
def getGeomPID(geometry):
    command = """surfaceSplitByPatch %s""" % (geometry)
    searchType = 'contains'
    searchVar = 'Zone'
    pid = getBashOutput(command,searchType,searchVar)
    pidArray = []
    
    for i in pid:
        pidArray.append(i.split('"')[1])   

    return pidArray

    
def search_and_replace(file_path, search_word, replace_word):
   with open(file_path, 'r') as file:
      file_contents = file.read()
      updated_contents = file_contents.replace(search_word, replace_word)
   with open(file_path, 'w') as file:
      file.write(updated_contents)  

def x_rotation(vector,theta):
    """Rotates 3-D vector around x-axis"""
    R = np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0, np.sin(theta), np.cos(theta)]])
    return np.dot(R,vector)

def y_rotation(vector,theta):
    """Rotates 3-D vector around y-axis"""
    R = np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta), 0, np.cos(theta)]])
    return np.dot(R,vector)

def z_rotation(vector,theta):
    """Rotates 3-D vector around z-axis"""
    R = np.array([[np.cos(theta), -np.sin(theta),0],[np.sin(theta), np.cos(theta),0],[0,0,1]])
    return np.dot(R,vector)   

def load_obj(f):
    vertices = []
    faces = []
    #with open(filename, 'r') as f:
    for line in f:
        if line.startswith('v '):
            vertices.append([float(x) for x in line.split()[1:]])
        elif line.startswith('f '):
            face = []
            for v in line.split()[1:]:  # Handle faces with more than 3 vertices
                face.append(int(v.split('/')[0]) - 1) #-1 because .obj indices start from 1
            if len(face)==3:
                faces.append(face)
            elif len(face) > 3: # Triangulate if needed
                for i in range(1, len(face)-1):
                    faces.append([face[0], face[i], face[i+1]])

    return np.array(vertices), np.array(faces)


def is_ascii_stl(f):
    """Checks if an STL file is likely ASCII format."""
    try:
        
        # Check the first few lines for "solid" keyword (ASCII STL starts with "solid")
        for _ in range(5):  # Check the first 5 lines
            line = f.readline()
            if "solid" in line.lower():  # Case-insensitive check
                return True
            elif line.strip() == "":  # Skip empty lines
                continue
            else: # Contains data that is not "solid"
                return False # Likely binary if not ASCII
        return False  # If "solid" not found in first 5 lines, likely binary
    except UnicodeDecodeError:
        return False  # If UnicodeDecodeError, it's likely binary

def calculate_moi_tensor(vertices, faces):
    """
    Calculates the moment of inertia tensor for a 3D object defined by vertices and faces.

    Args:
        vertices: A NumPy array of shape (n, 3) representing the 3D coordinates of the vertices.
        faces: A NumPy array of shape (m, 3) representing the vertex indices of each triangular face.

    Returns:
        A 3x3 NumPy array representing the moment of inertia tensor.  Returns None if the input is invalid or if a calculation error occurs.
    """

    if not isinstance(vertices, np.ndarray) or vertices.shape[1] != 3:
        print("Error: vertices must be a NumPy array of shape (n, 3).")
        return None

    if not isinstance(faces, np.ndarray) or faces.shape[1] != 3:
        print("Error: faces must be a NumPy array of shape (m, 3).")
        return None

    try:  # Catch potential errors like invalid indices in faces
        moi = np.zeros((3, 3))

        for face in faces:
            v1 = vertices[face[0]]
            v2 = vertices[face[1]]
            v3 = vertices[face[2]]

            # Calculate the volume of the tetrahedron formed by the triangle and the origin
            tet_volume = np.dot(np.cross(v1, v2), v3) / 6.0

            # Calculate integrals for the tetrahedron's inertia contribution
            # Using formulas derived from  ∫∫∫ xᵢxⱼ dV over the tetrahedron.
            #  (See https://en.wikipedia.org/wiki/List_of_moments_of_inertia or similar resource).

            for i in range(3):
                for j in range(3):
                    integral_term = 0
                    for k in range(3):  # Iterate over vertices of the face
                        integral_term += v1[k] + v2[k] + v3[k] if i==j and i==k else 0  # Diagonal terms are different
                        integral_term += v1[k] + v2[k]       if i==j and i!=k and k==0 else 0
                        integral_term += v1[k] + v3[k]       if i==j and i!=k and k==1 else 0
                        integral_term += v2[k] + v3[k]       if i==j and i!=k and k==2 else 0
                        integral_term += v1[k]               if i!=j and (k==i or k==j) else 0
                    moi[i, j] += tet_volume * integral_term / 20.0


        # The MOI tensor calculated above represents the inertia w.r.t. the origin.
        #  For many applications, you'll want the inertia tensor relative to the center of mass.
        total_mass = np.sum(moi.diagonal()) # Mass if uniform density
        center_of_mass = np.sum(vertices, axis=0) / vertices.shape[0]  # Approximate CoM

        # Apply the parallel axis theorem:
        moi -= total_mass * (np.outer(center_of_mass, center_of_mass) * np.eye(3) - np.outer(center_of_mass, center_of_mass))
        
        return moi
    except (IndexError, TypeError) as e:
        print(f"An error occurred during calculation: {e}")
        return None  # Indicate failure

def load_binary_stl(f):
    """Loads vertices and faces from an STL file (binary format)."""
    vertices = []
    faces = []

    f.seek(80)

    # Read the number of triangles
    num_triangles = struct.unpack('<I', f.read(4))[0]

    for _ in range(num_triangles):
        # Read the normal vector (we don't need it for MOI)
        f.read(12)

        # Read the vertices of the triangle
        triangle_vertices = []
        for _ in range(3):
            vertex = struct.unpack('<fff', f.read(12))
            triangle_vertices.append(vertex)
        vertices.extend(triangle_vertices) # Add the unique vertices
        
        # Define the face by the indices of the vertices
        faces.append([len(vertices)-3, len(vertices)-2, len(vertices)-1]) # Add a new face

        # Read the attribute byte count (2 bytes)
        f.read(2)

    return np.array(vertices), np.array(faces)

def load_ascii_stl(f):
    vertices = []
    faces = []
    vertex_index = 0  # Keep track of vertex indices for faces

    for line in f:
        if "vertex" in line:
            match = re.findall(r"vertex\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)", line)
            if match:
                vertices.append([float(x) for x in match[0]])
                vertex_index += 1  # Increment the vertex index
        elif "endfacet" in line: # Add face after processing all vertices of a facet
            faces.append([vertex_index-3, vertex_index-2, vertex_index-1]) # Define the face using last 3 vertex indices

    return np.array(vertices), np.array(faces)
    
def process_chunk(chunk, centroid):
    """Processes a chunk of vertices."""
    centered_chunk = chunk - centroid
    return np.dot(centered_chunk.T, centered_chunk)  # Calculate partial covariance

def find_wheel_axis(vertices,faces, num_processes = None):
    """
    Finds the rotation axis of a car wheel from an STL file.

    Args:
        filename: Path to the STL file.

    Returns:
        A NumPy array representing the unit vector of the wheel's rotation axis.
        Returns None if the file cannot be loaded or the axis cannot be determined.
    """
    try:
        print('\t\t\t\t\tCalculating rotational axis...')

        # 2. Calculate the best-fit plane (using SVD):
        centroid = np.mean(vertices, axis=0)
        centered_vertices = vertices - centroid
        num_vertices = len(vertices)

        # Parallel processing:
        if num_processes is None:
            num_processes = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=num_processes)
        chunk_size = num_vertices // num_processes
        results = []
        for i in range(num_processes):
            start = i * chunk_size
            end = (i + 1) * chunk_size if i < num_processes - 1 else num_vertices # last chunk could have more elements
            chunk = vertices[start:end]
            results.append(pool.apply_async(process_chunk, (chunk, centroid)))
        pool.close()
        pool.join()

        # Combine results:
        covariance_matrix = np.sum([r.get() for r in results], axis=0)

        # Incremental SVD on the covariance matrix:
        _, _, V = np.linalg.svd(covariance_matrix) # SVD is now performed on a much smaller matrix
        normal_vector = V[-1, :]
        axis = normal_vector / np.linalg.norm(normal_vector)
        return centroid[0], centroid[1], centroid[2], axis[0], -1*np.abs(axis[1]), axis[2]

    except Exception as e:
        print(f"Error: Could not determine wheel axis. {e}")
        return None

def incremental_svd(points, k=3):
    """Performs incremental SVD to find the k largest singular values/vectors."""
    A = np.zeros((k, k))
    for i, point in enumerate(points):
        A = A + np.outer(point, point) # Update the matrix incrementally
        if (i+1) % k == 0: # Perform SVD every few data points to make it more stable
            U, S, V = np.linalg.svd(A)
            A = np.diag(S).dot(V) # Update the covariance matrix

    U, S, V = np.linalg.svd(A)
    return U, S, V

def readGeomFile(fileName):
    
    filePath = 'constant/triSurface/%s' % (fileName)
    
    if '.gz' in filePath:
        print('\t\t\t\tUnzipping gz file...')
        try:
            file = gzip.open(filePath, 'rt')
        except:
            file = gzip.open(filePath, 'rb')
    else:
        file = open(filePath,'rt')

    if '.obj' in filePath:
        print('\t\t\t\tLoading obj file...')
        vertices, faces = load_obj(file)
    elif '.stl' in filePath:
        if is_ascii_stl(file):
            print('\t\t\t\tLoading ascii stl file...')
            vertices, faces = load_ascii_stl(file)
        elif is_ascii_stl(file) == False:
            print('\t\t\t\tLoading binary stl file...')
            vertices, faces = load_binary_stl(file)

    return vertices,faces


def calculate_planar_surface_geometry(vertices, faces, planarity_tolerance=1.0e-5):
    """Return area-weighted geometry data for a planar triangular surface.

    The normal is obtained from the least-variance PCA direction, so the
    result is independent of inconsistent OBJ/STL face winding.  Face areas
    are used for the centroid and equivalent circular diameter.  The returned
    normal has no meaningful sign; callers must orient it using a target point
    or an explicit user vector.
    """
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=int)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or len(vertices) < 3:
        raise ValueError('surface must contain at least three 3-D vertices')
    if faces.ndim != 2 or faces.shape[1] != 3 or len(faces) == 0:
        raise ValueError('surface must contain triangular faces')
    if np.any(faces < 0) or np.any(faces >= len(vertices)):
        raise ValueError('surface contains an invalid face index')

    tri = vertices[faces]
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    double_area = np.linalg.norm(cross, axis=1)
    valid = double_area > np.finfo(float).eps
    if not np.any(valid):
        raise ValueError('surface has zero area')

    tri = tri[valid]
    areas = 0.5 * double_area[valid]
    face_centres = np.mean(tri, axis=1)
    area_total = float(np.sum(areas))
    centre = np.sum(face_centres * areas[:, None], axis=0) / area_total

    centred = vertices - centre
    _, singular_values, vh = np.linalg.svd(centred, full_matrices=False)
    normal = vh[-1]
    normal_norm = np.linalg.norm(normal)
    if normal_norm <= np.finfo(float).eps:
        raise ValueError('unable to determine a surface normal')
    normal = normal / normal_norm

    scale = max(float(np.max(np.linalg.norm(centred, axis=1))), 1.0)
    planarity_error = float(singular_values[-1] / scale)
    if planarity_error > planarity_tolerance:
        raise ValueError(
            'surface is not planar (relative error %.6g > %.6g)' %
            (planarity_error, planarity_tolerance)
        )

    diameter = 2.0 * np.sqrt(area_total / np.pi)
    return {
        'center': centre,
        'normal': normal,
        'area': area_total,
        'diameter': float(diameter),
        'planarity_error': planarity_error,
    }


def calculate_surface_centroid(vertices, faces):
    """Return an area-weighted centroid for any triangulated surface."""
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=int)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or len(vertices) < 3:
        raise ValueError('surface must contain at least three 3-D vertices')
    if faces.ndim != 2 or faces.shape[1] != 3 or len(faces) == 0:
        raise ValueError('surface must contain triangular faces')
    tri = vertices[faces]
    areas = 0.5 * np.linalg.norm(
        np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1
    )
    valid = areas > np.finfo(float).eps
    if not np.any(valid):
        raise ValueError('surface has zero area')
    areas = areas[valid]
    centres = np.mean(tri[valid], axis=1)
    return np.sum(centres * areas[:, None], axis=0) / np.sum(areas)


def orient_surface_normal(normal, source_center, target_center):
    """Orient a surface normal from source_center toward target_center."""
    normal = np.asarray(normal, dtype=float)
    source_center = np.asarray(source_center, dtype=float)
    target_center = np.asarray(target_center, dtype=float)
    normal_norm = np.linalg.norm(normal)
    target_vector = target_center - source_center
    target_norm = np.linalg.norm(target_vector)
    if normal_norm <= np.finfo(float).eps:
        raise ValueError('surface normal is zero')
    if target_norm <= np.finfo(float).eps:
        raise ValueError('target geometry centre coincides with source centre')
    normal = normal / normal_norm
    target_vector = target_vector / target_norm
    if abs(float(np.dot(normal, target_vector))) <= 1.0e-10:
        raise ValueError('target geometry lies in the source surface plane')
    if np.dot(normal, target_vector) < 0.0:
        normal = -normal
    return normal





def writeGeometryFile(outputFile, vertices, faces):
    '''
    Write geometry to OBJ or STL file (with optional .gz compression).
    
    :param outputFile: Path to output file (determines format from extension)
    :param vertices: Array of vertex coordinates
    :param faces: Array of face vertex indices
    '''
    
    # Handle .gz compression
    if outputFile.endswith('.gz'):
        # Determine actual format from filename
        if '.obj' in outputFile:
            format_type = 'obj'
            actual_file = outputFile.replace('.gz', '')
        elif '.stl' in outputFile:
            format_type = 'stl'
            actual_file = outputFile.replace('.gz', '')
        else:
            sys.exit('ERROR! Output file must be .obj.gz or .stl.gz')
        
        # Write to temporary file first
        if format_type == 'obj':
            writeOBJFile(actual_file, vertices, faces)
        else:
            writeSTLFile(actual_file, vertices, faces, binary=True)
        
        # Compress to .gz
        print(f'\t\tCompressing to {outputFile}...')
        with open(actual_file, 'rb') as f_in:
            with gzip.open(outputFile, 'wb') as f_out:
                f_out.writelines(f_in)
        os.remove(actual_file)
    else:
        # Write uncompressed
        if '.obj' in outputFile:
            writeOBJFile(outputFile, vertices, faces)
        elif '.stl' in outputFile:
            # Auto-detect binary vs ASCII from user preference or default to binary
            writeSTLFile(outputFile, vertices, faces, binary=True)
        else:
            sys.exit('ERROR! Output file must be .obj or .stl format')


def writeOBJFile(filename, vertices, faces):
    '''Write vertices and faces to OBJ file.'''
    print(f'\t\tWriting OBJ file: {filename}')
    with open(filename, 'w') as f:
        # Write vertices
        for vertex in vertices:
            f.write(f'v {vertex[0]:.6f} {vertex[1]:.6f} {vertex[2]:.6f}\n')
        
        # Write faces (OBJ uses 1-based indexing)
        for face in faces:
            f.write(f'f {face[0]+1} {face[1]+1} {face[2]+1}\n')


def writeSTLFile(filename, vertices, faces, binary=True):
    '''Write vertices and faces to STL file (ASCII or binary).'''
    if binary:
        print(f'\t\tWriting binary STL file: {filename}')
        with open(filename, 'wb') as f:
            # Write 80-byte header
            header = b'Binary STL file created by utilities.py'.ljust(80, b'\0')
            f.write(header)
            
            # Write number of triangles
            f.write(struct.pack('<I', len(faces)))
            
            # Write each triangle
            for face in faces:
                v1 = vertices[face[0]]
                v2 = vertices[face[1]]
                v3 = vertices[face[2]]
                
                # Calculate normal
                edge1 = v2 - v1
                edge2 = v3 - v1
                normal = np.cross(edge1, edge2)
                norm_length = np.linalg.norm(normal)
                if norm_length > 0:
                    normal = normal / norm_length
                else:
                    normal = np.array([0, 0, 0])
                
                # Write normal
                f.write(struct.pack('<fff', *normal))
                
                # Write vertices
                f.write(struct.pack('<fff', *v1))
                f.write(struct.pack('<fff', *v2))
                f.write(struct.pack('<fff', *v3))
                
                # Write attribute byte count
                f.write(struct.pack('<H', 0))
    else:
        print(f'\t\tWriting ASCII STL file: {filename}')
        with open(filename, 'w') as f:
            f.write('solid geometry\n')
            
            for face in faces:
                v1 = vertices[face[0]]
                v2 = vertices[face[1]]
                v3 = vertices[face[2]]
                
                # Calculate normal
                edge1 = v2 - v1
                edge2 = v3 - v1
                normal = np.cross(edge1, edge2)
                norm_length = np.linalg.norm(normal)
                if norm_length > 0:
                    normal = normal / norm_length
                else:
                    normal = np.array([0, 0, 0])
                
                f.write(f'  facet normal {normal[0]:.6e} {normal[1]:.6e} {normal[2]:.6e}\n')
                f.write('    outer loop\n')
                f.write(f'      vertex {v1[0]:.6e} {v1[1]:.6e} {v1[2]:.6e}\n')
                f.write(f'      vertex {v2[0]:.6e} {v2[1]:.6e} {v2[2]:.6e}\n')
                f.write(f'      vertex {v3[0]:.6e} {v3[1]:.6e} {v3[2]:.6e}\n')
                f.write('    endloop\n')
                f.write('  endfacet\n')
            
            f.write('endsolid geometry\n')


def copyWithMkdir(src, dest):
    """Copy file and create destination directory if needed."""
    dest_dir = os.path.dirname(dest)
    if dest_dir and not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    shutil.copy2(src, dest)


def transformBlockmeshPoints(inputFile, outputFile=None, rotation=None, translation=None, scale=None, REFCOR=None):
    '''
    Transform blockmesh dictionary points by rotation, translation, and/or scaling.
    
    Reads a OpenFOAM blockMeshDict file, transforms all vertex points, and optionally
    writes to a new file.
    
    :param inputFile: Path to input blockMeshDict file
    :param outputFile: Path to output blockMeshDict file. If None, only returns transformed points
    :param rotation: Dict with keys 'x', 'y', 'z' for rotation angles in degrees
    :param translation: List/array [dx, dy, dz] for translation offsets
    :param scale: Scalar for uniform scaling, or list/array [sx, sy, sz] for non-uniform scaling
    :param REFCOR: Point [x, y, z] to rotate about. If None, rotates about origin (0, 0, 0)
    :return: Tuple of (vertices, new_content) - transformed vertices array and full modified blockMeshDict
    '''
    
    print(f'\t\t\tReading blockMeshDict from {inputFile}...')
    
    with open(inputFile, 'r') as f:
        content = f.read()
    
    # Extract vertices section using regex
    import re
    vertices_match = re.search(r'vertices\s*\((.*?)\);', content, re.DOTALL)
    
    if not vertices_match:
        sys.exit('ERROR! Could not find vertices section in blockMeshDict')
    
    vertices_text = vertices_match.group(1)
    
    # Parse vertices - format is typically: (x y z)
    vertex_pattern = r'\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)'
    vertices_list = re.findall(vertex_pattern, vertices_text)
    
    if not vertices_list:
        sys.exit('ERROR! Could not parse vertices from blockMeshDict')
    
    print(f'\tFound {len(vertices_list)} vertices')
    
    # Convert to numpy array for transformations
    vertices = np.array([[float(v[0]), float(v[1]), float(v[2])] for v in vertices_list], dtype=np.float64)
    
    # Set rotation center point
    if REFCOR is None:
        rotation_center = np.array([0.0, 0.0, 0.0])
    else:
        rotation_center = np.array([float(REFCOR[0]), float(REFCOR[1]), float(REFCOR[2])])
        print(f'\tRotating about REFCOR: {rotation_center}')
    
    # Apply scaling
    if scale is not None:
        print('\tApplying scaling...')
        if isinstance(scale, (int, float)):
            vertices *= scale
        else:
            vertices *= np.array(scale)
    
    # Apply rotation about REFCOR (order: x, y, z)
    if rotation is not None:
        print('\t\t\tApplying rotations about reference point...')
        rot_x = rotation.get('x', 0)
        rot_y = rotation.get('y', 0)
        rot_z = rotation.get('z', 0)
        
        # Translate vertices to origin, rotate, then translate back
        vertices -= rotation_center
        
        if rot_x != 0:
            rot_x_rad = math.radians(rot_x)
            for i in range(len(vertices)):
                vertices[i] = x_rotation(vertices[i], rot_x_rad)
        
        if rot_y != 0:
            rot_y_rad = math.radians(rot_y)
            for i in range(len(vertices)):
                vertices[i] = y_rotation(vertices[i], rot_y_rad)
        
        if rot_z != 0:
            rot_z_rad = math.radians(rot_z)
            for i in range(len(vertices)):
                vertices[i] = z_rotation(vertices[i], rot_z_rad)
        
        # Translate back to original reference point
        vertices += rotation_center
    
    # Apply translation
    if translation is not None:
        print('\tApplying translation...')
        vertices += np.array(translation)
    
    # Generate new vertices section
    new_vertices_text = 'vertices\n(\n'
    for vertex in vertices:
        new_vertices_text += f'    ({vertex[0]:.6e} {vertex[1]:.6e} {vertex[2]:.6e})\n'
    new_vertices_text += ');'
    
    # Replace vertices section in content
    new_content = re.sub(
        r'vertices\s*\((.*?)\);',
        new_vertices_text,
        content,
        flags=re.DOTALL
    )
    
    # Write output file if specified
    if outputFile is not None:
        output_dir = os.path.dirname(outputFile)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        print(f'\t\t\tWriting transformed blockMeshDict to {outputFile}...')
        with open(outputFile, 'w') as f:
            f.write(new_content)
    
    return vertices, new_content
