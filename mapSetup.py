import sys
import fileinput
import os
import getopt
import re
import glob
import argparse as argparse
import configparser
import numpy as np
import pandas as pd
import math
from glob import glob
path = os.path.split(os.getcwd())[0]
case = os.path.split(os.getcwd())[1]
jobPath = os.path.abspath(os.path.join(path,os.pardir))

archivePath = "/media/openfoam/archive/jobs"
scriptLoc = "/home/openfoam/openFoam/scripts"
templateLoc = "%s/templates" % (os.path.dirname(os.path.realpath(__file__)))
print(templateLoc)
if not os.path.isdir(archivePath):
    sys.exit("ERROR: Archive path incorrect!")

if not os.path.isdir(templateLoc):
    sys.exit("ERROR: Template path incorrect!")




parser = argparse.ArgumentParser(prog='mapSetup-v1.0',description='Sets up and executes ride hight and roll maps using a donor case.')

                    
parser.add_argument("-c","--case", action="store_true",
                    help='Creates case based on values given in the config file.')
parser.add_argument("-m","--mesh", action="store_true",
                    help='Sets up meshing dicts and links geometry files. Will copy template case files if the trial folder is empty.')
parser.add_argument('-s',"--fixScripts", action="store_true", 
                    help='Runs script fixer.')
parser.add_argument('-d',"--controlDict", action="store_true", 
                    help='Writes only controlDict for solver.')

args = parser.parse_args()
MESH = args.mesh
CASE = args.case
FIXSCRIPT = args.fixScripts
CONTROLDICT = args.controlDict



################## 
#### PREAMBLE ####
print("#### MAP SETUP V1.0 ####")
##################
caseSetupPath="%s/%s/caseSetup.cfg" % (path,case)
config = configparser.ConfigParser()
config.optionxform = str
config.read_file(open(caseSetupPath))
configSections = config.sections()
variableSet = {}
varList = []
varVal = []

for section in configSections:
   variableSet[section] = np.array(config.items(section)) 
   if not "GEOMETRY" in section:
      print("  Reading %s values:" % (section))
   
   for item in np.array(config.items(section)):
      if not "GEOMETRY" in section:
         print("     %s: %s" % (item[0],item[1]))
      varList = np.append(varList,item[0])
      varVal = np.append(varVal,item[1])

with open("%s/%s/caseSetup.cfg" % (path,case), 'w') as configfile:
  config.write(configfile)

def main():
    global runMapList
    getGeomArray(assignVar('GEOM'))
    checkMapSetup()
    getDonorGeometry()
    runMapList = createMapArray()
    transformGeom()


def checkMapSetup():
    
    sectionName = "MAP_SETUP"
    mapType = ["Ride_Height","Roll"]
    rideHeightData = ["Map_Trial_Start","Map_Type","frMin_frMax_frIncrement","rrMin_rrMax_rrIncrement","lhsMin_lhsMax_lhsIncrement","rhsMin_rhsMax_rhsIncrement","rollCenter"]
    defaultRideHeightVal = ["001_half","Ride_Height","-0.025 0.025 5","-0.025 0.025 5","-0.025 0.025 2","-0.025 0.025 2","1.2 0 0.01"]

    checkSetupFlag = False
    if sectionName in config.sections():
        if assignVar("Map_Type") != None:
            
            #check if ride height input is given
    
            for rhData,rhVal in zip(rideHeightData,defaultRideHeightVal):               
                if rhData in dict(config.items(sectionName)).keys():
                    if rhData not in rideHeightData[0:2]:
                        #checking to make sure the values are correct
                        for val in assignVar(rhData).split(" "):
                            try: 
                                val = float(val)      
                            except:
                                sys.exit("Error in MAP_SETUP, please make sure inputs are integers or float!")
                else:
                    #if a rhData is missing from the config file, write default values
                    checkSetupFlag = True
                    config[sectionName][rhData] = rhVal
        else:
            config[sectionName]["Map_Type"] = "Ride_Height"
            for rhData,rhVal in zip(rideHeightData,defaultRideHeightVal):               
                if rhData in dict(config.items(sectionName)).keys():
                    if rhData not in rideHeightData[0:2]:
                        #checking to make sure the values are correct
                        for val in assignVar(rhData).split(" "):
                            try: 
                                val = float(val)      
                            except:
                                sys.exit("Error in MAP_SETUP, please make sure inputs are integers or float!")
                else:
                    #if a rhData is missing from the config file, write default values
                    checkSetupFlag = True
                    config[sectionName][rhData] = rhVal
            checkSetupFlag = True
        
                    
                   
    else: 
        config[sectionName] = {}
        for rhData,rhVal in zip(rideHeightData,defaultRideHeightVal):
            config[sectionName][rhData] = rhVal
        checkSetupFlag = True
            
            
    if checkSetupFlag == True:
        writeCfgCheck()
   
    
    
    
def writeCfgCheck():
    with open("%s/%s/caseSetup.cfg" % (path,case), 'w') as configfile:
        config.write(configfile)
        
    sys.exit("ERROR: caseSetup.cfg needs to be updated due to new data added, please check!")
def assignVar(variable):
   
   if variable in varList:
      varIndex = np.where(varList == variable)
      variableValue = varVal[varIndex[0]]
      return variableValue[0]
   else:
      variableValue = None
      return variableValue
      
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
        #geomArray.append(row[0].replace('.stl',''))
        geomArray.append(row[0].split('.')[0])
        geomTypeArray.append("." + ".".join(row[0].split('.')[1:]))
        scaleArray.append(row[1])
        refArray.append(row[2])
        layerArray.append(row[3])
        expRatioArray.append(row[4])
        blSetup.append(row[5])
        if SINGLEGEOM == True:
            break
            
def getGeomPID(geometry):
    command = """surfaceSplitByPatch %s""" % (geometry)
    searchType = 'contains'
    searchVar = 'Zone'
    pid = getBashOutput(command,searchType,searchVar)
    pidArray = []
    
    for i in pid:
        pidArray.append(i.split('"')[1])   

    return pidArray

    
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



def getDonorGeometry():
    global donorGeomDict,donorGeomList
    donorGeomDict = {}
    for geomName,geomType in zip(geomArray,geomTypeArray):
        #initializing donorGeomDict with donor geometry names, the list will contained all the names of the transformed geometry
        donorGeomDict[geomName + geomType] = {}
        donorGeomDict[geomName + geomType]["fileNames"] = {}
        donorGeomDict[geomName + geomType]["fileExtension"] = geomType
    donorGeomList = donorGeomDict.keys()
    
def createMapCases(mapType,mapArray,mapMatrix,frhArray):
    #check if start case number is available and if all case names are available
    nRuns =  len(list(mapArray['frhM'].values())) * len(list(mapArray['rrhM'].values()))
    runStart = assignVar("Map_Trial_Start")
    runStartNum = int(runStart.split("_")[0])
    runMapList = ["{0:03}".format(i) for i in range(runStartNum,runStartNum+nRuns+1)] #creating trial name list
    try:
        runStartType = runStart.split("_")[1]
        runMapList = [s + "_" + runStartType for s in runMapList] #add suffix if exists
    except:
        pass
    
    
    #checking if there's conflicts
    caseList = os.listdir(jobPath + "/CASES")
    jobNumber = jobPath.split('/')[-1]
    archiveList = os.listdir(archivePath + "/" + jobNumber + "/" + "CASES")
    for i in archiveList:
        caseList.append(i)
    #print(caseList)
    for case in runMapList:
        if case in caseList:
            sys.exit("#### ERROR #### Map case sequence conflicts with existing cases!")
        else:
            pass
    
    print("\tWill create the following cases:")
    for case in runMapList:
        print("\t\t%s" %(case))   
    
    return runMapList
    
    
    
def createMapArray():
    print("\n\n\n")
    print("Creating map values...")
    global frh,rrh, pitchArray,rollAng,nRuns
    mapArray = {}
    
    #calculate the roll and pitch angles using the wheel center data
    frVals = np.asarray(assignVar("frMin_frMax_frIncrement").split(" "),dtype = float)
    rrVals = np.asarray(assignVar("rrMin_rrMax_rrIncrement").split(" "),dtype = float)
    FR_LHS_CENTER = np.asarray(assignVar("FR_LHS_CENTER").split(" "),dtype = float)
    RR_LHS_CENTER = np.asarray(assignVar("RR_LHS_CENTER").split(" "),dtype = float)
    
    frArray = np.linspace(frVals[0],frVals[1],int(frVals[2]))
    rrArray = np.linspace(rrVals[0],rrVals[1],int(rrVals[2]))
    mapArray['frh'] = {}
    mapArray['rrh'] = {}
    mapArray['frhM'] = {}
    mapArray['rrhM'] = {}
    for h,i in zip(frArray,range(len(frArray))):
        mapArray['frh'][i] = round(h,6)
        mapArray['frhM'][i] = int(h*1000)
    for h,i in zip(rrArray,range(len(rrArray))):
        mapArray['rrh'][i] = round(h,6)
        mapArray['rrhM'][i] = int(h*1000)
    
    try:
        frlhscent = float(assignVar("FR_LHS_CENTER").split(" ")[0])
        rrlhscent = float(assignVar("RR_LHS_CENTER").split(" ")[0])
        wheelBase = rrlhscent - frlhscent
    except:
        sys.exit("#### ERROR #### Wheel center coordinates invalid!")
    
    #front ride height is the row and rear ride height is the columns of the matrix
    
    pitchArray = np.zeros((int(frVals[2]),int(rrVals[2])))#init pitch angle array
    frct = 0
    while frct < int(frVals[2]):
        rrct = 0
        while rrct < int(rrVals[2]):
            pitchArray[frct][rrct] = math.degrees(math.asin((mapArray['frh'][frct]-mapArray['rrh'][rrct])/wheelBase))
            rrct = rrct + 1
        frct = frct + 1
    
    frct = 0
    while frct < int(frVals[2]):
        rrct = 0
        while rrct < int(rrVals[2]):
            pitchArray[frct][rrct] = math.degrees(math.asin((mapArray['frh'][frct]-mapArray['rrh'][rrct])/wheelBase))
            rrct = rrct + 1
        frct = frct + 1
    
    nRuns = len(list(mapArray['frhM'].values())) * len(list(mapArray['rrhM'].values()))
  
    if assignVar("Map_Type").lower() == 'ride_height':
        print("\tMap Type: Ride Height")
        print("\tFront Ride Height Changes (mm): %s" % (list(mapArray['frhM'].values())))
        print("\tRear Ride Height Changes (mm): %s" % (list(mapArray['rrhM'].values())))
        print("\tTotal Number of Runs: %s" % (nRuns))
        runMapList = createMapCases(assignVar("Map_Type").lower(),mapArray,pitchArray,mapArray['frh'])
        return runMapList
    elif assignVar("Map_Type").lower() == "roll":
        print("\tMap Type: Roll")
        sys.exit("####ERROR#### Roll not currently implemented, check back later...")
    else:
        sys.exit("####ERROR#### Incorrect map type!")
    

def matrix_multiply(*matrices):
    if len(matrices) == 1:
        return matrices
    else:
        try:
            m_other = matrix_multiply(*matrices[1:])
            return np.matmul(matrices[0], m_other)
        except:
            print(matrices[0])
            print(m_other)
            raise

def vectorRotation(p, transformVector , theta): #theta in degrees
    p = [[pp] for pp in p + [1]]
    x1 = [0,0,0]
    x1, y1, z1 = x1
    x2, y2, z2 = transformVector

    U = [x2-x1, y2-y1, z2-z1]
    U = np.array(U) / np.sqrt(np.dot(U,U))
    a,b,c = U
    d = np.sqrt(b**2 + c**2)

    T = [[1,0,0,-x1],[0,1,0,-y1],[0,0,1,-z1],[0,0,0,1]]
    T_inv = [[1,0,0,x1],[0,1,0,y1],[0,0,1,z1],[0,0,0,1]]

    R_x = [[1,0,0,0],[0,c/d,-b/d,0],[0,b/d,c/d,0],[0,0,0,1]]
    R_x_inv = [[1,0,0,0],[0,c/d,b/d,0],[0,-b/d,c/d,0],[0,0,0,1]]

    R_y = [[d,0,-a,0],[0,1,0,0],[a,0,d,0],[0,0,0,1]]
    R_y_inv = [[d,0,a,0],[0,1,0,0],[-a,0,d,0],[0,0,0,1]]

    ct = np.cos(math.radians(theta))
    st = np.sin(math.radians(theta))
    R_z = [[ct,st,0,0],[-st,ct,0,0],[0,0,1,0],[0,0,0,1]]

    p2 = matrix_multiply(T_inv, R_x_inv, R_y_inv, R_z, R_y, R_x, T, p)
    return p2[0][:3]      
    
    
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
        
    
    
def transformGeom():
    geometryList = assignVar("GEOM")
    geometryList = geometryList.splitlines()
    geometryList = np.genfromtxt(geometryList, delimiter=',', dtype=str)
    geometryToTransform = []
    for geometryLine in geometryList:
        if "-wh-" in geometryLine[0]:
            pass
        else:
            geometryToTransform.append(geometryLine[0])
    
    #get the geometries that are required to have a vector change
    vectorDataPrefix = {'POR':['VEC','POINT']}
    vectorsToTransform = {}
    #vectorsToTransform['DONOR'] = {}
    print("Getting input points/vectors...")
    for section in configSections:
        for key in vectorDataPrefix.keys():
            if key in section:
                print("\t%s" % (section))
                vectorsToTransform[section] = {}
                for subsection in config[section].keys():
                    for subkey in vectorDataPrefix[key]:
                        if subkey in subsection:
                            vectorsToTransform[section][subsection] = config[section][subsection]
                            print("\t\t" + subsection + " : " + vectorsToTransform[section][subsection])
    
    print("Generating new vectors for each run...")
    transformedVectors = {}
    
    for key in vectorsToTransform.keys():
        transformedVectors[key] = {}
        for run in runMapList:
            transformedVectors[key][run] = {}
            for prefKey in vectorDataPrefix.keys():
                if key.startswith(prefKey):
                    pass
            
    
    
    print(vectorsToTransform)
                       


main()                    
        
    







