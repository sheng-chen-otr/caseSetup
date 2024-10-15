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
import subprocess as sp


path = os.path.split(os.getcwd())[0]
case = os.path.split(os.getcwd())[1]
jobPath = os.path.abspath(os.path.join(path,os.pardir))
scriptLoc = "/home/openfoam/openFoam/scripts"
templateLoc = "%s/templates" % (os.path.dirname(os.path.realpath(__file__)))
print('Template path: ' + templateLoc)

if not os.path.isdir(templateLoc):
    sys.exit("ERROR: Template path incorrect!")

parser = argparse.ArgumentParser(prog='caseSetup-v3.0',description='Set us the case based on settings written out in the caseSetup.cfg')
                    
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
#### FUNCTIONS ####

def writeCfgCheck():
    with open("%s/%s/caseSetup.cfg" % (path,case), 'w') as configfile:
        config.write(configfile)
        
    sys.exit("ERROR: caseSetup.cfg needs to be updated due to new data added, please check!")

def checkDefaults():
    defaultDict = {}
    sections = ['TITLES','TEMPLATE','GEOMETRY','CONTROLSETUP','FLUIDPROP','FORCESETUP','WHEEL_SETUP','POSTPROSETUP']
    #initialize sections
    for section in sections:
        defaultDict[section] = {}
    
def writePostProSurfaceList():
    #slices this would require updating the slices from the previous versions
    #might not work well with the pv-post output names
    print("Writing out surface list...")
    sliceList = ['XSLICE','YSLICE','ZSLICE']
    sliceDict = {}
    for sliceType in sliceList:
        sliceTypeName = "%s" % (sliceType)
        try:
            sliceVals = np.array(assignVar(sliceTypeName).split(' ')).astype('float')
            sliceDict[sliceType] = {}
            sliceDict[sliceType]['start'] = sliceVals[0]
            sliceDict[sliceType]['end'] = sliceVals[1]
            sliceDict[sliceType]['ds'] = sliceVals[2]
            sliceDict[sliceType]['nslice'] = int((sliceVals[1]-sliceVals[0])/sliceVals[2])
            
            print('\t%s' % (sliceType))
            print('\t\tStart: %s'% (sliceVals[0]))
            print('\t\tEnd: %s' % (sliceVals[1]))
            print('\t\tDs: %s' % (sliceVals[2]))
            print('\t\tnSlices: %s' % (int((sliceVals[1]-sliceVals[0])/sliceVals[2])))
        except:
            print('\t%s not requested or invalid input, skipping...' % (sliceType))
    
    with open('system/surfaceSetupList', 'w') as surfaceSetupList:
        surfacesList = []
        for sliceType in sliceDict.keys():
            sliceRange = np.linspace(sliceDict[sliceType]['start'],sliceDict[sliceType]['end'],sliceDict[sliceType]['nslice'])
            n = 0
            for i in sliceRange:
                if 'X' in sliceType:
                    basePoint = [i,0,0]
                    normalVector = [1,0,0]
                    normal = 'x'
                    #print(basePoint)
                elif 'Y' in sliceType:
                    basePoint = [0,i,0]
                    normalVector = [0,1,0]
                    normal = 'y'
                elif 'Z' in sliceType:
                    basePoint = [0,0,i]
                    normalVector = [0,0,1]
                    normal = 'z'
                stringToAdd = """%s_%sNormal_%04d{type cuttingPlane; planeType pointAndNormal; pointAndNormalDict{ basePoint (%1.4f %1.4f %1.4f); normalVector (%s %s %s);} interpolate true;}\n""" % (sliceType,normal,float(n),basePoint[0],basePoint[1],basePoint[2],normalVector[0],normalVector[1],normalVector[2])
                n = n + 1
                surfacesList.append(stringToAdd)
        surfaceSetupList.write(''.join(surfacesList))
            

    
    
    
    
def checkCorner():
    
    requiredData = ["CORNER_RAD"]
    if assignVar("SYMMETRY") == "corner":
        print("Found case as corner, setting up corner case...")
        for i in requiredData:
            if assignVar(i) != None:
                print("\t%s: %0.4f" % (i,float(assignVar(i))))
            else:
                checkSetupFlag = True
                print("\tCorner Data not found, adding to caseSetup.cfg")
                config["FORCESETUP"]["CORNER_RAD"] = "16"
                writeCfgCheck()
                return checkSetupFlag
                
                
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

def checkForPorousData():
    
    requiredData = ["_DCOEFFS","_FCOEFFS","_VEC1","_VEC2","_POINT"]
    defaultData = ['100','100','1 0 0','0 1 0','0 0 0']
    getGeomArray(assignVar('GEOM'))
    print(" Porous Media Check...")
    for i in geomArray:
        
        if i.startswith("POR"):
            
            sectionName = i
            print("\t%s" % (sectionName))
            if sectionName in config.sections():
                
                for item in config.items(sectionName):
                    if any(x in item[0] for x in requiredData):
                        if item[1] != '':
                            print("\t\t%s: %s OK!" % (item[0],item[1]))
                        else:
                            print("\t\t%s: %s FAIL!" % (item[0],item[1]))
                            checkSetupFlag = True
                            return checkSetupFlag
                
            else:
                checkSetupFlag = True
                print("\t\tPorous Media Input Not Found, Adding %s to caseSetup..." % (i))
                config["%s" % (i)] = {}
                for data,val in zip(requiredData,defaultData):
                    config["%s" % (i)]["%s%s" % (i,data)] = val
                 
                writeCfgCheck()
    
def search_and_replace(file_path, search_word, replace_word):
   with open(file_path, 'r') as file:
      file_contents = file.read()
      updated_contents = file_contents.replace(search_word, replace_word)
   with open(file_path, 'w') as file:
      file.write(updated_contents)
def assignVar(variable):
   
   if variable in varList:
      varIndex = np.where(varList == variable)
      variableValue = varVal[varIndex[0]]
      return variableValue[0]
   else:
      variableValue = None
      return variableValue
def setupPorousMedia():
    porousArray = []
    porStringArray = []
    fieldAverageArray = []
    print("Setting up porous media...")
    print("     Found the following porous media:")
    for i in geomArray:
        if i.startswith("POR"):
            porousArray.append(i)
            print("             %s" % (i))
    
    if len(porousArray) < 1:
        print("             No porous media found...")
        if os.path.isfile("%s/%s/constant/fvOptions" % (path,case)):
            search_and_replace("%s/%s/constant/fvOptions" % (path,case), "<POROUS_MEDIA>"," ")
        
        if os.path.isfile("%s/%s/system/surfaceFieldAverage" % (path,case)):
            search_and_replace("%s/%s/system/surfaceFieldAverage" % (path,case), "<POROUS_MEDIA>"," ")
            
        
    else:
        print("Setting up porous media...")
        for i in porousArray:
            porousString = """
%s
{
type explicitPorositySource;
active yes;
explicitPorositySourceCoeffs {type DarcyForchheimer; selectionMode cellZone;cellZone %s_INTERNAL;
DarcyForchheimerCoeffs{d d [0 -2 0 0 0 0 0] (%s 3e10 3e10);f f [0 -1 0 0 0 0 0] (%s 1e5 1e5); coordinateSystem {type cartesian; origin  (0 0 0); coordinateRotation {type axesRotation;e1 (%s);e2 (%s);}}}}}\n
                            """ % (i,i,assignVar("%s_DCOEFFS" % (i)),assignVar("%s_FCOEFFS" % (i)),assignVar("%s_VEC1" % (i)),assignVar("%s_VEC2" % (i)))
            porIndex = geomArray.index(i)
            porousType = geomTypeArray[porIndex]
            porousFile = "%s%s" % (i,porousType)
            porousFilePath = "%s/%s/constant/triSurface/%s%s" % (path,case,i,porousType)
            print("\tGetting surface PIDs for %s..." % (porousFile))
            pid = getGeomPID(porousFilePath)
            print("\t\tFound the following PID's:")
            for surface in pid:
                print("\t\t\t%s" % (surface.split('.')[0].split('/')[-1]))
                fieldAverageString = """%s{type surfaceFieldValue; surfaceFormat vtk; libs (fieldFunctionObjects);fields (U UMean);operation  areaNormalAverage;regionType sampledSurface;name sampledSurface; sampledSurfaceDict{type sampledTriSurfaceMesh; surface %s; source cells; interpolate true;} writeFields true;writeToFile true;writeControl timeStep;}\n""" % (surface.split('.')[0].split('/')[-1],surface.split('/')[-1])
                fieldAverageArray.append(fieldAverageString)
            
            porStringArray.append(porousString)
        porStringArray = "".join(porStringArray)
        fieldAverageArray = "".join(fieldAverageArray)
        if os.path.isfile("%s/%s/constant/fvOptions" % (path,case)):
            print("\tWriting out porous media...")
            search_and_replace("%s/%s/constant/fvOptions" % (path,case), "<POROUS_MEDIA>",porStringArray)
        else:
            sys.exit("      ERROR! fvOptions not valid!")
        
        if os.path.isfile("%s/%s/system/surfaceFieldAverage" % (path,case)):
            print("\tWriting out porous media to surfaceFieldAverage...")
            if len(fieldAverageArray) < 1:
                print("\tWARNING: No porous surfaces written out to surfaceFieldAverage, check files or caseSetup!")
                search_and_replace("%s/%s/system/surfaceFieldAverage" % (path,case), "<POROUS_MEDIA>",'')
            elif len(fieldAverageArray) > 0:
                search_and_replace("%s/%s/system/surfaceFieldAverage" % (path,case), "<POROUS_MEDIA>",fieldAverageArray)
            else:
                print("\tWARNING: No porous surfaces written out to surfaceFieldAverage, check files or caseSetup!")
                search_and_replace("%s/%s/system/surfaceFieldAverage" % (path,case), "<POROUS_MEDIA>",'')
    
    
def getGeomArray(geometry):
    global geomArray,blSetup,geomTypeArray
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
        if ".stl" in row[0]:
            geomTypeArray.append(".stl")
        elif ".obj" in row[0]:
            geomTypeArray.append(".obj")
        elif ".nas" in row[0]:
            geomTypeArray.append(".nas")
        # if ".gz" in row[0]:
            # geomTypeArray.append("." + ".".join(row[0].split('.')[1:]).replace(".gz",""))
        # else:
            # geomTypeArray.append(".".join(row[0].split('.')[1:]))
        #geomTypeArray.append("." + ".".join(row[0].split('.')[1:]).replace(".gz",""))
        scaleArray.append(row[1])
        refArray.append(row[2])
        layerArray.append(row[3])
        expRatioArray.append(row[4])
        blSetup.append(row[5])
        if SINGLEGEOM == True:
            break
    
def addSnappyGeometry(geometry,template):
   global geomArray,blSetup,GEOM,SINGLEGEOM
   TEMPLATE = template
   GEOM = geometry
   
   GEOM = GEOM.splitlines()
   if len(GEOM) < 2:
      SINGLEGEOM = True
   else:
      SINGLEGEOM = False
   GEOM = np.genfromtxt(GEOM, delimiter=',', dtype=str)
   
   print("\n\n\n\nUsing %s for snappyHexMeshDict..." % (TEMPLATE))
   print("\n\n\n")
   print("Getting geometries...")
   
   if SINGLEGEOM == True:
      
      print("     " + GEOM[0])
      print("         Scale:" + GEOM[1])
      print("         Ref:" + GEOM[2])
      print("         nLayers:" + GEOM[3])
      print("         expR:" + GEOM[4])
      print("         BL:" + GEOM[5] + "\n")
   else:
      for row in GEOM:
          print("     " + row[0])
          print("         Scale:" + row[1])
          print("         Ref:" + row[2])
          print("         nLayers:" + row[3])
          print("         expR:" + row[4])
          print("         BL:" + row[5] + "\n")
   geomArray = []
   scaleArray = []
   refArray = []
   layerArray = []
   expRatioArray = []
   blSetup = []
   geomTypeArray = []
   
   if os.path.isdir("constant/triSurface"):
      os.system("rm -r constant/triSurface")
      os.system("mkdir constant/triSurface")
   for geom in GEOM:
      if SINGLEGEOM == True:
         geom = GEOM
            
      if not os.path.isfile("%s/02_reference/MSH/%s" % (jobPath,geom[0])):
         sys.exit("#### ERROR ####\n   TriSurface: %s not found in MSH folder!" % (geom[0]))
      if os.path.isfile("constant/triSurface/%s" % (geom[0])):
         print("     %s is already linked, skipping..." % (geom[0]))
         continue
      print("     Creating symbolic link for %s in triSurface" % (geom[0]))
      cmd = "ln -s ../../../../02_reference/MSH/%s constant/triSurface/%s" % (geom[0],geom[0])
      
      os.system(cmd)
      if SINGLEGEOM == True:
         break
   for row in GEOM:
       if SINGLEGEOM == True:
         row = GEOM
       #geomArray.append(row[0].replace('.stl',''))
       geomArray.append(row[0].split('.')[0])
       if ".stl" in row[0]:
           geomTypeArray.append(".stl")
       elif ".obj" in row[0]:
           geomTypeArray.append(".obj")
       elif ".nas" in row[0]:
           geomTypeArray.append(".nas")
       #geomTypeArray.append("." + ".".join(row[0].split('.')[1:]))
       scaleArray.append(row[1])
       refArray.append(row[2])
       layerArray.append(row[3])
       expRatioArray.append(row[4])
       blSetup.append(row[5])
       if SINGLEGEOM == True:
         break
    
   GEOM_STRING=[]
   FEAT_STRING=[]
   REF_STRING=[]
   LAYER_STRING=[]
   REGION_STRING=[]
   nGeom = range(len(geomArray))
   for n in nGeom:
       #print(geomArray[n])
       geomString=""
       featureString=""
       refinementString=""
       layerString=""
       if geomArray[n].startswith("POR"):
            geomString = """%s{type distributedTriSurfaceMesh;scale %s;file "%s%s";}\n""" % (geomArray[n],scaleArray[n],geomArray[n],geomTypeArray[n])
       elif geomArray[n].startswith("REF"):
            geomString = """%s{type triSurfaceMesh;scale %s;file "%s%s";}\n""" % (geomArray[n],scaleArray[n],geomArray[n],geomTypeArray[n])   
       else:
            geomString = """%s{type distributedTriSurfaceMesh;scale %s;file "%s%s";regions{".*";}}\n""" % (geomArray[n],scaleArray[n],geomArray[n],geomTypeArray[n])
       
       if geomArray[n].startswith("REF"):
            pass
       else:
            featureString = """{file "%s.eMesh";scale %s;level 10;}\n""" % (geomArray[n],scaleArray[n])
       
       if geomArray[n].startswith("POR"):
            refinementString = """%s{level (10 10);faceZone %s;cellZone %s_INTERNAL;cellZoneInside insidePoint;insidePoint (%s);}\n""" % (geomArray[n], geomArray[n],geomArray[n],assignVar("%s_POINT" % (geomArray[n])))
       
       elif geomArray[n].startswith("REF"):
            pass
       else:
            refinementString = """%s{level (%s %s);regions{#include"snappyRefinementDict"}}\n""" % (geomArray[n], refArray[n],refArray[n])
       
       if geomArray[n].startswith("REF"):
            pass
       else:
            layerString = """".*%s.*"{nSurfaceLayers %s;expansionRatio %s;}\n""" % (geomArray[n],layerArray[n],expRatioArray[n])
       
       if geomArray[n].startswith("REF"):
            regionRefString = """%s{mode inside;levels ((1E15 %s));}\n""" % (geomArray[n], refArray[n])
       GEOM_STRING.append(geomString)
       FEAT_STRING.append(featureString)
       REF_STRING.append(refinementString)
       LAYER_STRING.append(layerString)
       try:
            REGION_STRING.append(regionRefString)
       except:
            pass
       
   GEOM_STRING="".join(GEOM_STRING)
   FEAT_STRING="".join(FEAT_STRING)
   REF_STRING="".join(REF_STRING)
   LAYER_STRING="".join(LAYER_STRING)
   REGION_STRING="".join(REGION_STRING)
       
   print("Copying %s to trial folder..." % (TEMPLATE))
   if os.path.isfile("/home/openfoam/openFoam/templates/snappyTemplates/%s" % (TEMPLATE)):
       os.system("cp /home/openfoam/openFoam/templates/snappyTemplates/%s %s/%s/system/snappyHexMeshDict" % (TEMPLATE,path,case))
   else:
       sys.exit("Template not found! Please select valid template!")
           
   print("         Writing out geometry...")
   search_and_replace("%s/%s/system/snappyHexMeshDict" % (path,case), "<GEOMETRY>",GEOM_STRING)
   print("         Writing out featureEdges...")
   search_and_replace("%s/%s/system/snappyHexMeshDict" % (path,case), "<FEATURE_EDGE>",FEAT_STRING)
   print("         Writing out refinementSurfaces...")
   search_and_replace("%s/%s/system/snappyHexMeshDict" % (path,case), "<REFINEMENT_SURFACES>",REF_STRING)
   print("         Writing out refinementRegions...")
   search_and_replace("%s/%s/system/snappyHexMeshDict" % (path,case), "<REFINEMENT_REGIONS>",REGION_STRING)
   print("         Writing out layers...")
   search_and_replace("%s/%s/system/snappyHexMeshDict" % (path,case), "<LAYERS>",LAYER_STRING)
   #print("\n\n\n\n\n")
   print("#### Finished creating snappyHexMeshDict! ####")
   print("\n\n\n\n\n")
def getCaseFiles():
   #global geomArray,blSetup
   simtype = assignVar('SIMTYPE')
   sym = assignVar('SYMMETRY')
   common = assignVar('COMMON_TEMPLATE')
   TEMPLATE = assignVar('TEMPLATE')
   GEOM = assignVar('GEOM')
   
   geomArray = []
   geomTypeArray = []
   blSetup = []

   if not os.path.isdir("%s/%s" % (templateLoc,simtype)):
      sys.exit("#### ERROR #####\nSIMTYPE: %s not valid! Please check your caseSetup.cfg file!" % (simtype))
   elif not os.path.isdir("%s/%s/%s" % (templateLoc,simtype,sym)):
      sys.exit("#### ERROR #####\nSYMMETRY: %s not valid! Please check your caseSetup.cfg file!" % (sym))
   else:
      print("Copying files for:")
      print("  SIMTYPE: %s" % (simtype))
      print("  SYMMETRY: %s" % (sym))
      os.system("cp -r %s/%s/%s/* %s/%s/." % (templateLoc,simtype,sym,path,case))
      os.system("cp -r %s/%s/%s/* %s/%s/." % (templateLoc,simtype,common,path,case))

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
def setDecomposition():
    
    print("Setting decomposeParDict for: ")
    print("    nCores: %s" % (assignVar("NUMCORES")))
    print("    hierarchicalDecomp: %s" % (assignVar("DECOMP")))

    for i in varList:
      search_and_replace("%s/%s/system/decomposeParDict" % (path,case), "<%s>" % (i),"%s" % (assignVar(i)))

def velVector(inletMag,yaw,pitch):
    initVel = [float(inletMag), 0, 0]
    pitch = math.radians(-1*pitch)
    yaw = math.radians(yaw)
    initDragVec = [1,0,0]
    initLiftVec = [0,0,1]
    #transform inlet vector for yaw
    yawTransVel = z_rotation(initVel,yaw)
    #transform inlet vector for pitch
    pitchTransVel = y_rotation(yawTransVel,pitch)
    #" ".join(str(x) for x in xs)
    #transform drag and lift vectors for pitch
    dragVec = y_rotation(initDragVec,pitch)
    liftVec = y_rotation(initLiftVec,pitch)
    return " ".join(str(x) for x in pitchTransVel), " ".join(str(x) for x in dragVec)," ".join(str(x) for x in liftVec)
def wheelRotVel(inletMag,radius):
   rotVel = float(inletMag)/float(radius)
   return rotVel
def appendRunScript(scriptName,command,appendAfterCmd):
    files = glob.glob("%s*" % (scriptName))

    for fileName in files:
        newFile = []
        repeat = False
        
        with open(fileName,'r') as file:
        
            #check if the file already contains the same command, if so break
            for line in file:
                if command in line:
                    print('\t\tCommand already in file, skipping!')
                    repeat = True

            
        with open(fileName,'r') as file:
            #this time if above command is not met
            if repeat == False:    
                for line in file:
                    newFile.append(line)
                    if appendAfterCmd in line:
                        newFile.append(command)
                fileWrite = open(fileName,'w')
                fileWrite.writelines(newFile)
                
            
def checkPitchType():
    #get transform tunnel variable boolean
    
    try: 
        transformTunnel = assignVar('TRANSFORM_TUNNEL').capitalize()
        if transformTunnel.lower() == 'true':
            transformTunnel = True
        elif transformTunnel.lower() == 'false':
            transformTunnel = False
        else:
            transformTunnel = False
            print('\n WARNING: TRANSFORM_TUNNEL value invalid, using default False!\n')
    except:
        transformTunnel = False #defaults to false if not given

    #if tunnel move is not enabled, then check for wheel transform, if not set, default to false!
    if transformTunnel == False:
        
        try:
            transformWheel == assignVar('TRANSFORM_WHEEL').capitalize()
            if transformWheel.lower() == 'true':
                transformWheel = True
            elif transformWheel.lower() == 'false':
                transformWheel = False
            else:
                transformWheel = False
                print('\n WARNING: TRANSFORM_WHEEL value invalid, using default False!\n')
        except:
            transformWheel = False
    #if tunnel move is enabled, then check for wheel transform, if not set default to true!
    elif transformTunnel == False:
        try:
            transformWheel == assignVar('TRANSFORM_WHEEL').capitalize()
            if transformWheel.lower() == 'true':
                transformWheel = True
            elif transformWheel.lower() == 'false':
                transformWheel = False
            else:
                transformWheel = False
                print('\n WARNING: TRANSFORM_WHEEL value invalid, using default False!\n')
        except:
            transformWheel = True #defaults to true!
    else:
        try:
            transformWheel == assignVar('TRANSFORM_WHEEL').capitalize()
            if transformWheel.lower() == 'true':
                transformWheel = True
            elif transformWheel.lower() == 'false':
                transformWheel = False
            else:
                transformWheel = False
                print('\n WARNING: TRANSFORM_WHEEL value invalid, using default False!\n')
        except:
            transformWheel = True #defaults to true!
        
    return transformTunnel, transformWheel
    
        
        
    
def writeCaseProp():
   print("Setting up caseProperties:")
   inletMag = float(assignVar('INLETMAG'))
   yaw = float(assignVar('YAW'))
   try: 
        pitch = float(assignVar('PITCH'))
   except:
        pitch = 0
        
   if transformTunnel == False:
        pitch = 0
   if "half" in case:
      print("   Found case as half case, using Yaw = 0...")
      yaw = 0
   elif "half" in assignVar("SYMMETRY"):
      print("   Found case as half case, using Yaw = 0...")
      yaw = 0
   elif "corner" in assignVar("SYMMETRY"):
      print("   Found case as corner case, using Yaw = 0...")
      yaw = 0
      cornerRad = float(assignVar("CORNER_RAD"))

      domainRotVel = inletMag/cornerRad
      domainRotVelRad = domainRotVel
      domainRotVel = domainRotVel*(60/(2*math.pi))
      corPoint = assignVar("REFCOR").split(" ")[0]
      domainOrigin = "%s %s 0" % (corPoint,-1*cornerRad)

      #domainRotVector = "0 0 %s" % (domainRotVel)
      search_and_replace("%s/%s/constant/SRFProperties" % (path,case), "<DOM_ROT>","%s" % (domainRotVel))
      search_and_replace("%s/%s/constant/SRFProperties" % (path,case), "<DOM_ORIG>","%s" % (domainOrigin))
      


   if not "corner" in assignVar("SYMMETRY"):
      inletVec,dragVec,liftVec = velVector(inletMag,yaw,pitch)
      groundCond = assignVar('GROUND')
      if groundCond not in ['moving','stationary']:
         sys.exit("INCORRECT GROUND CONDITION!")
      print("  Inlet Velocity: %s" % (inletMag))
      print("  Yaw: %s" % (yaw))
      print("  Pitch: %s" % (pitch))
      print("  Inlet Velocity Vector: %s\n" % (inletVec))
      print("  Ground Motion: %s\n" % (groundCond))
      search_and_replace("%s/%s/system/caseProperties" % (path,case), "<INLET_VELOCITY>","%s" % (inletVec))
      bcList = []
      for geom,bl in zip(geomArray,blSetup):
         if "REF" in geom:
            continue
         elif "POR" in geom:
            continue
            
         if bl == "high":
            wfunc = "highReynolds"
         elif bl == "low":
            wfunc = "lowReynolds"
         else:
               sys.exit("#### ERROR #####\nWall Function option for %s not valid!" % (geom))
         print("  setting up %s with %s wall function" % (geom,wfunc))
         if groundCond == "moving":
            
            if not "-wh-" in geom:
               
               wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion stationary;}values{$:initialConditions;}}\n""" % (geom,geom,wfunc)
               bcList.append(wallBc)
               
            elif all(i in geom for i in ["fr","wh","lhs"]):
               center = assignVar("FR_LHS_CENTER")
               axis = assignVar("FR_LHS_AXIS")
               rad = assignVar("FR_RAD")
               rotvel = wheelRotVel(inletMag,rad)
               
               wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion rotating;}values{type rotatingWallVelocity;origin (%s);axis (%s);rotVel -%s; $:initialConditions;}}\n""" % (geom,geom,wfunc,center,axis,rotvel)
               bcList.append(wallBc)
            elif all(i in geom for i in ["fr","wh","rhs"]):
               center = assignVar("FR_RHS_CENTER")
               axis = assignVar("FR_RHS_AXIS")
               rad = assignVar("FR_RAD")
               rotvel = wheelRotVel(inletMag,rad)
               
               wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion rotating;}values{type rotatingWallVelocity;origin (%s);axis (%s);rotVel -%s; $:initialConditions;}}\n""" % (geom,geom,wfunc,center,axis,rotvel)
               if not "half" in [assignVar("SYMMETRY"),case]: 
                  bcList.append(wallBc)
            elif all(i in geom for i in ["rr","wh","lhs"]):
               center = assignVar("RR_LHS_CENTER")
               axis = assignVar("RR_LHS_AXIS")
               rad = assignVar("RR_RAD")
               rotvel = wheelRotVel(inletMag,rad)
              
               wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion rotating;}values{type rotatingWallVelocity;origin (%s);axis (%s);rotVel -%s; $:initialConditions;}}\n""" % (geom,geom,wfunc,center,axis,rotvel)
               bcList.append(wallBc)
            elif all(i in geom for i in ["rr","wh","rhs"]):
               center = assignVar("RR_RHS_CENTER")
               axis = assignVar("RR_RHS_AXIS")
               rad = assignVar("RR_RAD")
               rotvel = wheelRotVel(inletMag,rad)
               
               
               wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion rotating;}values{type rotatingWallVelocity;origin (%s);axis (%s);rotVel -%s; $:initialConditions;}}\n""" % (geom,geom,wfunc,center,axis,rotvel)
               if not "half" in [assignVar("SYMMETRY"),case]:
                  bcList.append(wallBc)
         else:
            
            if "half" in [assignVar("SYMMETRY"),case] and "rhs" in geom: 
               continue

            wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion stationary;}values{$:initialConditions;}}\n""" % (geom,geom,wfunc)
            bcList.append(wallBc)
      bcList="".join(bcList)
      search_and_replace("%s/%s/system/caseProperties" % (path,case), "<PATCHBC>","%s" % (bcList))
      groundBc = """ ground{category wall;type noSlip;patches (z-min); options {wallFunction highReynolds; motion %s;} values{$:initialConditions;}}\n """ % (groundCond)
      search_and_replace("%s/%s/system/caseProperties" % (path,case), "<GROUND_PATCH>","%s" % (groundBc))
      print("Writing out to cpMeanDict...")
      search_and_replace("%s/%s/system/cpMeanDict" % (path,case), "<INLET_VELOCITY>","%s" % (inletVec))
      print("Writing out to ctpMeanDict...")
      search_and_replace("%s/%s/system/ctpMeanDict" % (path,case), "<INLET_VELOCITY>","%s" % (inletVec))
   else:
      inletVec,vec = velVector(inletMag,yaw)
      groundCond = assignVar('GROUND')
      if groundCond not in ['moving','stationary']:
         sys.exit("INCORRECT GROUND CONDITION!")
      print("  Inlet Velocity: %s" % (inletMag))
      print("  Corner Radius: %s" % (cornerRad))
      print("  Domain Rotation Speed: %s rpm" % (domainRotVel))
      print("  Ground Motion: %s" % (groundCond))
      search_and_replace("%s/%s/system/caseProperties" % (path,case), "<INLET_VELOCITY>", inletVec)
      bcList = []
      for geom,bl in zip(geomArray,blSetup):
         if bl == "high":
            wfunc = "highReynolds"
         elif bl == "low":
            wfunc = "lowReynolds"
         else:
               sys.exit("#### ERROR #####\nWall Function option for %s not valid!" % (geom))
         print("  setting up %s with %s wall function" % (geom,wfunc))
         if groundCond == "moving":
            if not "-wh-" in geom:
               
               wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion stationary;}values{$:initialConditions;}}\n""" % (geom,geom,wfunc)
               bcList.append(wallBc)
            elif all(i in geom for i in ["fr","wh","lhs"]):
               center = assignVar("FR_LHS_CENTER")
               axis = assignVar("FR_LHS_AXIS")
               rad = assignVar("FR_RAD")
               rotvel = wheelRotVel(inletMag,rad)
               
               wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion rotating;}values{type rotatingWallVelocity;origin (%s);axis (%s);rotVel -%s; $:initialConditions;}}\n""" % (geom,geom,wfunc,center,axis,rotvel)
               bcList.append(wallBc)
            elif all(i in geom for i in ["fr","wh","rhs"]):
               center = assignVar("FR_RHS_CENTER")
               axis = assignVar("FR_RHS_AXIS")
               rad = assignVar("FR_RAD")
               rotvel = wheelRotVel(inletMag,rad)
               
               wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion rotating;}values{type rotatingWallVelocity;origin (%s);axis (%s);rotVel -%s; $:initialConditions;}}\n""" % (geom,geom,wfunc,center,axis,rotvel)
               if not "half" in [assignVar("SYMMETRY"),case]: 
                  bcList.append(wallBc)
            elif all(i in geom for i in ["rr","wh","lhs"]):
               center = assignVar("RR_LHS_CENTER")
               axis = assignVar("RR_LHS_AXIS")
               rad = assignVar("RR_RAD")
               rotvel = wheelRotVel(inletMag,rad)
              
               wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion rotating;}values{type rotatingWallVelocity;origin (%s);axis (%s);rotVel -%s; $:initialConditions;}}\n""" % (geom,geom,wfunc,center,axis,rotvel)
               bcList.append(wallBc)
            elif all(i in geom for i in ["rr","wh","rhs"]):
               center = assignVar("RR_RHS_CENTER")
               axis = assignVar("RR_RHS_AXIS")
               rad = assignVar("RR_RAD")
               rotvel = wheelRotVel(inletMag,rad)
               
               
               wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion rotating;}values{type rotatingWallVelocity;origin (%s);axis (%s);rotVel -%s; $:initialConditions;}}\n""" % (geom,geom,wfunc,center,axis,rotvel)
               if not "half" in [assignVar("SYMMETRY"),case]:
                  bcList.append(wallBc)
         else:
            
            if "half" in [assignVar("SYMMETRY"),case] and "rhs" in geom: 
               continue

            wallBc = """   %s{category wall;type noSlip;patches (".*%s.*");options {wallFunction %s;motion stationary;}values{$:initialConditions;}}\n""" % (geom,geom,wfunc)
            bcList.append(wallBc)
      bcList="".join(bcList)
      search_and_replace("%s/%s/system/caseProperties" % (path,case), "<PATCHBC>","%s" % (bcList))
      groundBc = """ground{category wall;type noSlip;patches (z-min);options {wallFunction highReynolds;motion rotating;}values{type rotatingWallVelocity;origin (%s);axis (0 0 -1);rotVel %s; $:initialConditions;}}\n""" % (domainOrigin,domainRotVelRad)
      search_and_replace("%s/%s/system/caseProperties" % (path,case), "<GROUND_PATCH>","%s" % (groundBc))
      print("Writing out to cpMeanDict...")
      search_and_replace("%s/%s/system/cpMeanDict" % (path,case), "<INLET_VELOCITY>","%s" % (inletVec))
      print("Writing out to ctpMeanDict...")
      search_and_replace("%s/%s/system/ctpMeanDict" % (path,case), "<INLET_VELOCITY>","%s" % (inletVec))

def writeForceCoeffSetup():
   inletMag = float(assignVar('INLETMAG'))
   yaw = float(assignVar('YAW'))
   try:
        pitch = float(assignVar('PITCH'))
   except:
        pitch = 0
   if transformTunnel != True:
        pitch = 0
   if "half" in case:
      print("   Found case as half case, using Yaw = 0!")
      yaw = 0
   elif "half" in assignVar("SYMMETRY"):
      print("   Found case as half case, using Yaw = 0!")
      yaw = 0
   inletVec,dragVec,liftVec = velVector(inletMag,yaw,pitch)
   search_and_replace("%s/%s/system/forceCoeffSetup" % (path,case), "<DRAGVEC>","%s" % (dragVec))
   search_and_replace("%s/%s/system/forceCoeffSetup" % (path,case), "<LIFTVEC>","%s" % (liftVec))
   search_and_replace("%s/%s/system/forceCoeffsExport" % (path,case), "<DRAGVEC>","%s" % (dragVec))
   search_and_replace("%s/%s/system/forceCoeffsExport" % (path,case), "<LIFTVEC>","%s" % (liftVec))

   print("Setting up forceCoeffSetup...")
   for i in varList:
      search_and_replace("%s/%s/system/forceCoeffSetup" % (path,case), "<%s>" % (i),"%s" % (assignVar(i)))
   
   print("Setting up forceCoeffExport...")
   for i in varList:
      search_and_replace("%s/%s/system/forceCoeffsExport" % (path,case), "<%s>" % (i),"%s" % (assignVar(i)))
def writeForceCoeff():
   geomList = []
   geomOutList = []
   for geom in geomArray:
      geomList.append("""".*%s.*" """ % (geom))
      geomOut = """%s{patches (".*%s.*"); #include "forceCoeffSetup"}\n""" % (geom,geom)
      geomOutList.append(geomOut)
   
   geomString = ''.join(geomList)
   allLine = """all{patches (%s); #include "forceCoeffSetup"}\n""" % (geomString)
   geomOutList.append(allLine)
   geomOutList = "".join(geomOutList)
   print("Setting up forceCoeffs...")
   search_and_replace("%s/%s/system/forceCoeffs" % (path,case), "<FORCE_PATCHES>",geomOutList)
def writePostProDict():
   geomList = []
   for geom in geomArray:
      geomList.append("""".*%s.*" """ % (geom))
   geomString = ''.join(geomList)
   
   print("Setting up forceCoeffsExport...")
   search_and_replace("%s/%s/system/forceCoeffsExport" % (path,case), "<ALL_PATCHES>",geomString)
   print("Setting up nearWallFieldsDict...")
   search_and_replace("%s/%s/system/nearWallFieldsDict" % (path,case), "<ALL_PATCHES>",geomString)
   print("Setting up wallShearStressDict...")
   search_and_replace("%s/%s/system/wallShearStressDict" % (path,case), "<ALL_PATCHES>",geomString)
def writeConstant():
   constantList = os.listdir("%s/%s/constant/" % (path,case))
   print("Setting up constant dicts:")
   for const in constantList:
      print("  %s..." % (const))
      for i in varList:
         if const not in ["triSurface","polyMesh","extendedFeatureEdgeMesh"]:
            search_and_replace("%s/%s/constant/%s" % (path,case,const), "<%s>" % (i),"%s" % (assignVar(i)))
def writeToSurfaceFeature(path,case,GEOM):
   
   print("Writing out surfaceFeatureExtractDict...")
   open('%s/%s/system/surfaceFeatureExtractDict' % (path,case),'w').close() #clearing contents of surfaceFeatureExtract
   header = """FoamFile
   {
       version     2.0;
       format      ascii;
       class       dictionary;
       object      surfaceFeatureExtractDict;
   }
   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // """
   textArray = [header]
   for geom in GEOM:
      
      if SINGLEGEOM == True:
        geom = GEOM
      
        
     
      if geom[0].startswith("REF"):
            continue
            
      else:
          print("     Writing out %s to surfaceFeatureExtractDict" % (geom[0]))
          textArray.append(geom[0].replace(".gz",""))
          textArray.append('{')
          textArray.append('   #include "extractSetupDict"')
          textArray.append('}')
      if SINGLEGEOM == True:
        break
   with open('%s/%s/system/surfaceFeatureExtractDict' % (path,case),'w') as surfFile:
      surfFile.write('\n'.join(textArray))
def writeControlDict():
   print("Setting up controlDict...")
   simtype = assignVar('SIMTYPE')
   sym = assignVar('SYMMETRY')
   common = assignVar('COMMON_TEMPLATE')
   os.system("cp -r %s/%s/%s/system/controlDict* %s/%s/system/." % (templateLoc,simtype,common,path,case))
   for i in varList:
      
      search_and_replace("%s/%s/system/controlDictSimpleFoam" % (path,case), "<%s>" % (i),"%s" % (assignVar(i)))
      search_and_replace("%s/%s/system/controlDict" % (path,case), "<%s>" % (i),"%s" % (assignVar(i)))
      if os.path.isfile("%s/%s/system/controlDictPisoFoam" % (path,case)):
         search_and_replace("%s/%s/system/controlDictPisoFoam" % (path,case), "<%s>" % (i),"%s" % (assignVar(i)))
      if os.path.isfile("%s/%s/system/controlDictSRFSimpleFoam" % (path,case)):
         search_and_replace("%s/%s/system/controlDictSRFSimpleFoam" % (path,case), "<%s>" % (i),"%s" % (assignVar(i)))
def writeAverageFieldsDict():
   print("Setting up averageFieldsDict...")
   avgList = ['averageFieldsDict','averageFieldsDictPisoFoam','averageFieldsDictSimpleFoam']
   for j in avgList:
      if os.path.isfile("%s/%s/system/%s" % (path,case,j)):
         for i in varList:
            search_and_replace("%s/%s/system/%s" % (path,case,j), "<%s>" % (i),"%s" % (assignVar(i)))
def getConfigVariables(caseSetupPath):
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
    return config,variableSet,varList,varVal,configSections      

def writeSurfaces():
   geomList = []
   geomOutList = []
   for geom in geomArray:
      geomList.append("""".*%s.*" """ % (geom))
      geomOut = """%s{$patchSurface;patches (".*%s.*");}\n""" % (geom,geom)
      geomOutList.append(geomOut)
   
   geomString = ''.join(geomList)
   allLine = """all{$patchSurface;patches (%s);}\n""" % (geomString)
   geomOutList.append(allLine)
   geomOutList = "".join(geomOutList)
   print("Setting up surfaces...")
   search_and_replace("%s/%s/system/surfaces" % (path,case), "<SURFACE_PATCHES>",geomOutList)
   for i in varList:
      search_and_replace("%s/%s/system/surfaces" % (path,case), "<%s>" % (i),"%s" % (assignVar(i)))
def fixScripts(job):
    print("Changing run script trial names...")
    #################	Meshing Script	########################
    mesh=glob('meshingScript*')
    for meshName in mesh:
        if os.path.isfile("%s" % (meshName)):
            for line in fileinput.input("%s" % (meshName), inplace=True):
                if line.strip().startswith('#SBATCH --job-name='):
                    line=("#SBATCH --job-name=%s_MS%s \n" % (job,case))
                sys.stdout.write(line)
            print("     Meshing Script for %s_MS%s in %s" % (job,case,meshName))
    mesh=glob('fidelityMeshing*')
    for meshName in mesh:
        if os.path.isfile("%s" % (meshName)):
            for line in fileinput.input("%s" % (meshName), inplace=True):
                if line.strip().startswith('#SBATCH --job-name='):
                    line=("#SBATCH --job-name=%s_MS%s \n" % (job,case))
                sys.stdout.write(line)
            print("     Fidelity Meshing Script for %s_MS%s in %s" % (job,case,meshName))
    #################	Solve Script	########################
    solve=glob('solveScript-v*')
    for solveName in solve:
        if os.path.isfile("%s" % (solveName)):
            for line in fileinput.input("%s" % (solveName), inplace=True):
                if line.strip().startswith('#SBATCH --job-name='):
                    line=("#SBATCH --job-name=%s_SO%s \n" % (job,case))
                elif line.strip().startswith('RESTART='):
                    if assignVar("RESTART") in ["no","No","NO"]:
                        line=("RESTART=FALSE\n")
                    elif assignVar("RESTART") in ["yes","Yes","YES"]:
                        line=("RESTART=TRUE\n")
                    else:
                        sys.exit("RESTART variable has incorrect entry: " + assignVar("RESTART"))
                sys.stdout.write(line)
            print("     Solve Script for %s_SO%s in %s" % (job,case,solveName))
        
    solveSA=glob('solveScriptSA-v*')
    for solveSAName in solveSA:
        if os.path.isfile("%s" % (solveSAName)):
            for line in fileinput.input("%s" % (solveSAName), inplace=True):
                if line.strip().startswith('#SBATCH --job-name='):
                    line=("#SBATCH --job-name=%s_SO%s \n" % (job,case))
                elif line.strip().startswith('RESTART='):
                    if assignVar("RESTART") in ["no","No","NO"]:
                        line=("RESTART=FALSE\n")
                    elif assignVar("RESTART") in ["yes","Yes","YES"]:
                        line=("RESTART=TRUE\n")
                    else:
                        sys.exit("RESTART variable has incorrect entry: " + assignVar("RESTART"))
                sys.stdout.write(line)
            print("     Solve Script for %s_SO%s in %s" % (job,case,solveSAName))   
    solveTr=glob('solveScript*Transient*')
    for solveTrName in solveTr:
        if os.path.isfile('%s' % (solveTrName)):
            for line in fileinput.input("%s" % (solveTrName), inplace=True):
                if line.strip().startswith('#SBATCH --job-name='):
                    line=("#SBATCH --job-name=%s_TR%s \n" % (job,case))
                elif line.strip().startswith('RESTART='):
                    if assignVar("RESTART") in ["no","No","NO"]:
                        line=("RESTART=FALSE\n")
                    elif assignVar("RESTART") in ["yes","Yes","YES"]:
                        line=("RESTART=TRUE\n")
                    else:
                        sys.exit("RESTART variable has incorrect entry: " + assignVar("RESTART"))
                sys.stdout.write(line)
            print("     Solve Script Transient for %s_SO%s in %s" % (job,case,solveTrName)) 
             
    #################	Export Script	########################
    export=glob('exportScript*')
    for exportName in export:
        if os.path.isfile('%s' % (exportName)):
            for line in fileinput.input("%s" % (exportName), inplace=True):
                if line.strip().startswith('#SBATCH --job-name='):
                    line=("#SBATCH --job-name=%s_EX%s \n" % (job,case))
                sys.stdout.write(line)
            print("     Export Script for %s_EX%s in %s" % (job,case,exportName))
    #################	Post Pro Script	########################
    export=glob('postProScript-v*')
    for exportName in export:
        if os.path.isfile('%s' % (exportName)):
            for line in fileinput.input("%s" % (exportName), inplace=True):
                if line.strip().startswith('#SBATCH --job-name='):
                    line=("#SBATCH --job-name=%s_PP%s \n" % (job,case))
                sys.stdout.write(line)
            print("     Post Pro Script for %s_PP%s in %s" % (job,case,exportName))
   ################# Post Script   ########################
    export=glob('postScript-v*')
    for exportName in export:
        if os.path.isfile('%s' % (exportName)):
            for line in fileinput.input("%s" % (exportName), inplace=True):
                if line.strip().startswith('#SBATCH --job-name='):
                    line=("#SBATCH --job-name=%s_PO%s \n" % (job,case))
                sys.stdout.write(line)
            print("     Post Script for %s_PO%s in %s" % (job,case,exportName))

################## 
#### PREAMBLE ####
print("#### CASE SETUP V-3.0 ####")
##################

caseSetupPath="%s/%s/caseSetup.cfg" % (path,case)
config,variableSet,varList,varVal,configSections = getConfigVariables(caseSetupPath)         
checkCorner()          
checkForPorousData()
transformTunnel,transformWheel = checkPitchType()
        
        
        
with open("%s/%s/caseSetup.cfg" % (path,case), 'w') as configfile:
  config.write(configfile)




if CASE == True or os.path.isdir("%s/%s/constant" % (path,case)) == False or os.path.isdir("%s/%s/system" % (path,case)) == False:
    
    getCaseFiles()
    writeCaseProp()
    writeForceCoeffSetup()
    writeForceCoeff()
    writeConstant()
    writePostProDict()
    writeControlDict()
    writeAverageFieldsDict()
    writeSurfaces()
    setDecomposition()
    writePostProSurfaceList()
    
    if not MESH == True:
        getGeomArray(assignVar('GEOM'))
        setupPorousMedia()
        
if MESH == True:
    addSnappyGeometry(assignVar('GEOM'),assignVar('TEMPLATE'))
    writeToSurfaceFeature(path,case,GEOM)
    setupPorousMedia()    
    
  
if FIXSCRIPT == True:
    fixScripts(assignVar("JOBNAME"))
if CONTROLDICT == True:
    writeControlDict()

print("#### CASE SETUP COMPLETE! ####")
