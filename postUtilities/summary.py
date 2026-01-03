import os
import numpy as np
import fileinput
import sys
import re
import math
import argparse as argparse
import glob as glob
import pandas as pd
import scipy.stats as st
import configparser
from estimateStatisticalError import *

print("#### CASE SUMMARY ####")
# parser = argparse.ArgumentParser(prog='CASE SUMMARY',description='Summarizes the case into one csv file.')

# args = parser.parse_args()

# path = os.path.split(os.getcwd())[0]
# case = os.path.split(os.getcwd())[1]
# job = os.path.basename(os.path.dirname(path))

# if not os.getcwd().split('/')[-2].lower() == 'cases':
#     sys.exit('ERROR! Not executed in a trial directory!')
# print('Reading caseSetup file...')
# caseSetupPath="%s/%s/caseSetup" % (path,case)
# fullCaseSetupDict = configparser.ConfigParser()
# fullCaseSetupDict.optionxform = str
# fullCaseSetupDict.read_file(open(caseSetupPath))
# configSections = fullCaseSetupDict.sections()

def main():
    
    coeffFiles = getCoeffPaths()
    for part in coeffFiles:
        if part != 'all':
            avgData = averageCoeffs(case,part,coeffFiles)
    avgData = averageCoeffs(case,'all',coeffFiles)
    numCells,mesher,sym = cellCount()
    inletMag,lastTime,yaw,movingGround,rotatingWheels,simType,turbModel = bcParser(path,case)
    runDate,runTime,version,solver = getOfVersion()
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

    


def getGeometry(fullCaseSetupDict):    
    geomDict = {}
    geomColumnNames = ['scale','refinement','layers','expansion','wallmodel']
    #getting geometry list from geometry section
    try:
        geometryList = fullCaseSetupDict['GEOMETRY']['GEOM']
    except:
        print('ERROR! [GEOMETRY] section invalid!')
        exit()
    #breaking down geometry list into geometry dict
    geomDict = geomToDict(geomDict,geometryList,geomColumnNames)
    return geomDict
           

def geomToDict(geomDict,geometryList,geomColumnNames):
    
    try:
        print('\t\tGetting geometry...')
        geometryList = geometryList.split('\n')
        print(geometryList)
        nGeoms = len(geometryList)
        print('\t\t\tFound: %1.0f' % (nGeoms))
    except:
        print('ERROR! Geometry input invalid!')
        
    #splitting geoms into columns
    for geometry in geometryList:
        geometry = geometry.split(',')
        print(geometry)
        
        if len(geometry) != len(geomColumnNames) + 1:
            print('ERROR! Geometry input invalid! Insufficient options input for: %s' % (geometry[0]))
            exit()
        else:
            geometryName = geometry[0]
            if geometryName in geomDict.keys():
                print('ERROR! Duplicate geometry name found: %s' % (geometryName))
                exit()
            geomDict[geometryName] = {}
            n = 0
            while n < (len(geomColumnNames)):
                geomCol = geomColumnNames[n]
                geomDict[geometryName][geomCol] = geometry[n+1]
                n = n + 1
    
    
    return geomDict

    
def getGeomArray(geometry):
    global geomArray,geomTypeArray,GEOM
    geomArray = []
    scaleArray = []
    refArray = []
    layerArray = []
    expRatioArray = []
    blSetup = []
    geomTypeArray = []
    GEOM = geometry
    GEOM = GEOM.splitlines()
    if len(GEOM) < 2:
        SINGLEGEOM = True
    else:
        SINGLEGEOM = False
    GEOM = np.genfromtxt(GEOM, delimiter=',', dtype=str)
    for row in GEOM:
        if SINGLEGEOM == True:
            row = GEOM
        geomArray.append(row[0].split('.')[0])
        geomTypeArray.append("." + ".".join(row[0].split('.')[1:]))
        scaleArray.append(row[1])
        refArray.append(row[2])
        layerArray.append(row[3])
        expRatioArray.append(row[4])
        blSetup.append(row[5])
        if SINGLEGEOM == True:
            break

def magnitude(vector):
    return math.sqrt(sum(pow(element, 2) for element in vector))
    
def bcParser(fullCaseSetupDict,path,case):
    print("Getting boundary conditions...")
    with open('%s/%s/system/controlDict' % (path,case)) as controlDict:
        lines = controlDict.readlines()
    time = []
    for line in lines:
        if line.startswith("endTime"):
            time.append(re.search('[0-9]+',line).group(0))
            lastTime = time[0]
        elif line.startswith("application"):
            application = re.search('application\s* (.+?);',line).group(1)
        with open('%s/%s/system/caseProperties' % (path,case)) as caseProperties:
            lines = caseProperties.readlines()
            wheelRotation = False
            for line in lines:
                if "U" in line:
                    inletVel = np.array(re.findall('[0-9]+',line))
                    n=0
                    inletVel = map(float,inletVel)
                    inletVel = list(inletVel)
                elif "rotating" in line:
                    wheelRotation = True
                elif "ground" in line:
                    groundFlag = 1
                    wheelFlag = 0
                elif "*wh*" in line:
                    wheelFlag = 1
                    groundFlag = 0
                   
    inletMag = fullCaseSetupDict['BC_SETUP']['INLET_MAG']
    turbModel = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['TURB_MODEL']
    simType = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_TYPE']
    movingGround = fullCaseSetupDict['BC_SETUP']['GROUND']
    if 'half' in case or simType.lower() == 'half':
        yaw = 0
    else:
        yaw = fullCaseSetupDict['BC_SETUP']['YAW'][0]
        
    
    return inletMag,lastTime,yaw,movingGround,wheelRotation,simType,turbModel

def getCoeffPaths(path,case):
    print("Getting force coefficients...")
    #check if postProcessing dir exists
    postProPath = "%s/postProcessing" % (path)
    coeffFiles = {}
    excludeParts = ['allExports','binForceCoeffs','yPlusMean','wallShearStress','images','.DS_Store']
    if os.path.isdir(postProPath):
        for part in os.listdir(postProPath):
            if part in excludeParts or any(x in part for x in excludeParts):
                continue

            else:
                partPath = "%s/%s" % (postProPath,part)
                
                #get list of times in part
                timeList = os.listdir(partPath)
                if len(timeList) > 0:
                    if len(glob.glob("%s/*/coefficient*.dat" % (partPath))) > 0:
                        coeffFiles[part] = {}
                    else:
                        continue
                    for time in timeList:
                        coeffFilePath = glob.glob("%s/%s/coefficient*.dat" % (partPath,time))
                        if len(coeffFilePath)>0:                            
                            coeffFiles[part][time] = coeffFilePath[-1]
                else:
                    continue 

    return coeffFiles
def averageCoeffs(fullCaseSetupDict,case,part,coeffFiles):
    print("\tAveraging force coefficients for %s..." % (part))
    simType = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0]
    avgStart = float(fullCaseSetupDict['GLOBAL_CONTROL']['AVGSTART'][0])
    dataHeader = 'time,cd,cdf,cdr,cl,clf,clr,cmpitch,cmroll,cmyaw,cs,csf,csr'
    variables = ['cd','cdf','cdr','cl','clf','clr','cmpitch','cmroll','cmyaw','cs','csf','csr']
    dataHeader = dataHeader.split(",")
    coeffs = pd.DataFrame(columns=dataHeader)

    for time in coeffFiles[part].keys():
        timeCoeffs = pd.read_csv(coeffFiles[part][time],skiprows=13,delim_whitespace=True,names = dataHeader)
        coeffs = pd.concat([coeffs,timeCoeffs],ignore_index=True,axis=0)
    endTime = coeffs['time'].iloc[-1]
    avgStartRows = coeffs[coeffs['time'] >= avgStart]
    dt = coeffs['time'].iloc[-1] - coeffs['time'].iloc[-2]
    averagedData = {}
    for var in variables:
        data = avgStartRows[var]
        results = estimate_statistical_error(data,dt)
        
        if 'half' in case or fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower() == 'half':
            averagedData[var] = round(results['total_mean']*2,3)
        else:
            averagedData[var] = round(results['total_mean'],3)
        averagedData[var+'_ci'] = round(np.abs(results['mean_95_confidence_interval'][1] - results['mean_95_confidence_interval'][0])/2,4)

    

    cop = round(((averagedData['clf'])/(averagedData['cl'])) * 100,2)
    clcd = round(averagedData['cl']/averagedData['cd'],3)
    averagedData['endTime'] = endTime
    averagedData['cop'] = cop
    averagedData['cl/cd'] = clcd
    avgs = [averagedData['cd'],averagedData['cl'],averagedData['clf'],averagedData['clr'],averagedData['csf'],averagedData['csr'],averagedData['cd_ci'],averagedData['cl_ci']]
    np.savetxt("trial%s_AVG_%s_coeff.csv" % (case, part), avgs, delimiter=",",header="Time,CD,CL,CLF,CLR,CSF,CSR,CI-CD,CI-CL")
    return averagedData
        

def cellCount(fullCaseSetupDict,path,case):
    print("Getting cell counts...")
    checkmesh = "%s/log.checkMesh" % (path)
    logfidelity = "%s/log.fidelityMesh" % (path)
    logsnappy = "%s/log.snappyHexMesh" % (path)
    if os.path.isfile(logfidelity):
        mesher = "Fidelity Hexpress"
    elif os.path.isfile(logsnappy):
        mesher = "snappyHexMesh"
    else:
        mesher = "N/A"
    if "half" in case or fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower() == 'half':
        sym = "Half"
    elif "corver" in case:
        sym = "Corner"
    else:
        sym = "Full"
    if os.path.isfile(checkmesh):
        checkmeshfile = open(checkmesh, "r")
        for line in checkmeshfile:
            if "cells:" in line:
                numCells = line.split()[1]
    else:
        numCells = "N/A"
    return numCells, mesher, sym
def getOfVersion(path):
    print("Getting OF version...")
    logsimplefoam = "%s/log.simpleFoam" % (path)
    logpisofoam = "%s/log.pisoFoam" % (path)
    if os.path.isfile(logpisofoam):
        logpath = logpisofoam
        solver = "pisoFoam"
    elif os.path.isfile(logsimplefoam):
        logpath = logsimplefoam
        solver = "simpleFoam"
    else:
        logpath = False
        solver = "N/A"
    if not logpath == False:
        checksolvefile = open(logpath,"r")
        for line in checksolvefile:
            if "Date" in line:
                runDate = ' '.join(line.split()[2:])
                #print(runDate)
            elif "ExecutionTime" in line:
                runTime = round(float(line.split()[2])/3600,2)
            elif "Build" in line:
                version = line.split()[2]
    else:
        runDate = "N/A"
        runTime = "N/A"
        version = "N/A"
    return runDate, runTime,version,solver
#get the information required for radiator flow rates   
def getPorousData():
    print("Getting porous data...")
    postProcessingPath = "%s/%s/postProcessing" % (path,case)
    if os.path.isdir("postProcessing"):
        ppFileList = []
        for dir in os.listdir(postProcessingPath):
            if os.path.isdir("%s/%s" % (postProcessingPath,dir)) and dir.startswith('POR'):
                surfVal = glob.glob("%s/%s/*/surfaceFieldValue.dat*" % (postProcessingPath,dir),recursive=True)
                ppFileList.extend(surfVal)
        porousData = {}
        if len(ppFileList) < 1:
            return porousData
        for file in ppFileList:
            surfaceName = file.split('/')[-3]
            porousData[surfaceName + ' Volume Flow Rate (m^3/s)'] = []
            porousData[surfaceName + ' Velocity (m/s)'] = []
            porousData[surfaceName + ' Area (m^2)'] = []
            surfaceVals = open(file,'r').readlines()
            surfArea = 0
            surfVel = 0
            volFlowRate = 0
            n = 1
            surfaceFieldAvg = pd.read_csv(file,skiprows=4,delimiter='\t')
            for line in surfaceVals:
                #print(line)
                if 'Area' in line:
                    surfArea = float(line.split(':')[1])
                n = n + 1
            print(surfaceFieldAvg.columns)
            surfVel = surfaceFieldAvg['areaNormalAverage(UMean)'][0]
            volFlowRate = surfArea*surfVel
            porousData[surfaceName + ' Volume Flow Rate (m^3/s)'] = str(round(volFlowRate,3))
            porousData[surfaceName + ' Velocity (m/s)'] = str(round(surfVel,3))
            porousData[surfaceName + ' Area (m^2)'] = str(round(surfArea,3))
    else:
        porousData = {}
    return porousData
    






# main()