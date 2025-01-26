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

print("#### EZ-CFD CASE SUMMARY ####")
parser = argparse.ArgumentParser(prog='EZ-CFD CASE SUMMARY',description='Summarizes the case into one csv file.')

args = parser.parse_args()

path = os.path.split(os.getcwd())[0]
case = os.path.split(os.getcwd())[1]
job = os.path.basename(os.path.dirname(path))
print('Reading caseSetup file...')
caseSetupPath="%s/%s/caseSetup" % (path,case)
fullCaseSetupDict = configparser.ConfigParser()
fullCaseSetupDict.optionxform = str
fullCaseSetupDict.read_file(open(caseSetupPath))
configSections = fullCaseSetupDict.sections()




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
# def assignVar(variable):
#    if variable in varList:
#       varIndex = np.where(varList == variable)
#       variableValue = varVal[varIndex[0]]
#       return variableValue[0]
#    else:
#       variableValue = False
#       return variableValue
def magnitude(vector):
    return math.sqrt(sum(pow(element, 2) for element in vector))
    
def bcParser(path,case):
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
    if 'half' in case or simType.lower() == 'half':
        yaw = 0
    else:
        yaw = fullCaseSetupDict['BC_SETUP']['YAW'][0]
        
    
    return inletMag,lastTime,yaw,wheelRotation,simType,turbModel

def getCoeffPaths():
    print("Getting force coefficients...")
    #check if postProcessing dir exists
    postProPath = "%s/%s/postProcessing" % (path,case)
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
def averageCoeffs(case,part,coeffFiles):
    print("\tAveraging force coefficients for %s..." % (part))
    simType = fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'][0]
    avgStart = float(fullCaseSetupDict['GLOBAL_CONTROL']['AVGSTART'][0])
    dataHeader = 'time,cd,cdf,cdr,cl,clf,clr,cmpitch,cmroll,cmyaw,cs,csf,csr'
    dataHeader = dataHeader.split(",")
    coeffs = pd.DataFrame(columns=dataHeader)

    for time in coeffFiles[part].keys():
        timeCoeffs = pd.read_csv(coeffFiles[part][time],skiprows=13,delim_whitespace=True,names = dataHeader)
        coeffs = pd.concat([coeffs,timeCoeffs],ignore_index=True,axis=0)
    endTime = coeffs['time'].iloc[-1]
    avgStartRows = coeffs[coeffs['time'] >= avgStart]
    avgRow = avgStartRows.mean(axis=0)
    confInt = avgStartRows.apply(lambda x: st.t.interval(0.95, len(x)-1, loc=np.mean(x), scale=st.sem(x)), axis=0)
    confInt = confInt.diff().iloc[1,:]
    if 'half' in case or fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower() == 'half':
        cd = round(avgRow['cd'] * 2,4)
        cd_ci = round(confInt['cd'] * 2,4)
        cl = round(avgRow['cl'] * 2,4)
        cl_ci = round(confInt['cl'] * 2,4)
        clf = round(avgRow['clf'] * 2,4)
        clr = round(avgRow['clr'] * 2,4)
        csf = round(avgRow['csf'] * 2,4)
        csr = round(avgRow['csr'] * 2,4)
        cop = round(((avgRow['clf'])/(avgRow['cl'])) * 100,2)
    else:
        cd = round(avgRow['cd'],4)
        cd_ci = round(confInt['cd'],4)
        cl = round(avgRow['cl'],4)
        cl_ci = round(confInt['cl'],4)
        clf = round(avgRow['clf'],4)
        clr = round(avgRow['clr'],4)
        csf = round(avgRow['csf'],4)
        csr = round(avgRow['csr'],4)
        cop = round(((avgRow['clf'])/(avgRow['cl'])) * 100,2)
    clcd = round(cl/cd,3)
    
        
    averagedData = pd.DataFrame(columns=['endTime','cd','cl','clf','clr','csf','csr','cd_ci','cl_ci','cop','cl/cd'])
    averagedData['endTime'] = [endTime]
    averagedData['cd'] = [cd]
    averagedData['cl'] = [cl]
    averagedData['clf'] = [clf]
    averagedData['clr'] = [clr]
    averagedData['cop'] = [cop]
    averagedData['cl/cd'] = [clcd]
    averagedData['cd_ci'] = [cd_ci]
    averagedData['cl_ci'] = [cl_ci]
    averagedData['csf'] = [csf]
    averagedData['csr'] = [csr]
    return averagedData
        

def cellCount():
    print("Getting cell counts...")
    checkmesh = "%s/%s/log.checkMesh" % (path,case)
    logfidelity = "%s/%s/log.fidelityMesh" % (path,case)
    logsnappy = "%s/%s/log.snappyHexMesh" % (path,case)
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
def getOfVersion():
    print("Getting OF version...")
    logsimplefoam = "%s/%s/log.simpleFoam" % (path,case)
    logpisofoam = "%s/%s/log.pisoFoam" % (path,case)
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
    




#geomDict = getGeometry(fullCaseSetupDict)  
#geomList

getPorousData()
coeffFiles = getCoeffPaths()
avgData = averageCoeffs(case,'all',coeffFiles)
numCells,mesher,sym = cellCount()
inletMag,lastTime,yaw,wheelRotation,simType,turbModel = bcParser(path,case)
runDate,runTime,version,solver = getOfVersion()
refArea = float(fullCaseSetupDict['BC_SETUP']['REFAREA'][0])
#default datas
rowNames = ['Job','Trial','Solver','Version','Run Date','Solve Time','Num. Cells','Mesher','Symmetry','Ref. Area (m^2)','Iterations','Simulation Type','Moving Ground','Turbulence Model','Velocity','Yaw','Cd','Cl','Cl/Cd','%Front','Cd CI','Cl CI']
data = [job,case,solver,version,runDate,runTime,numCells,mesher,sym.lower(),refArea,avgData['endTime'].values[0],simType.lower(),wheelRotation,turbModel,inletMag,yaw,avgData['cd'].values[0],avgData['cl'].values[0],avgData['cl/cd'].values[0],avgData['cop'].values[0],avgData['cd_ci'].values[0],avgData['cl_ci'].values[0]]
porousData = getPorousData()
#adding porous data
for key in porousData.keys():
    rowNames.append(str(key))
    data.append(str(porousData[key]))
summary = pd.DataFrame(columns=rowNames)
summary.loc[-1] = data
print("\n\n")
for col in summary.columns:
    print('{:>100s}{:>30s}'.format(col,str(summary[col].values[0])))
    
summary = summary.transpose()
summary.to_csv("%s/%s/summary.csv"% (path,case),header=False)
