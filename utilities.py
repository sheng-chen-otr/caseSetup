import os
import sys
import numpy as np
import subprocess as sp
import math
try:
    from paraview.simple import *
except:
    print('ERROR! Unable to import paraview modules, skipping dependant operations!')

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