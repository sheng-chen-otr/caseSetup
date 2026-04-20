import os
import sys
import numpy as np
import subprocess as sp
import math
import gzip
import re
import multiprocessing
import pandas as pd
import configparser

updateCaseSetupFlag = False





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
    execDir = os.path.dirname(os.path.realpath(__file__))
    rideHeightFilePath = os.path.join(execDir,'rideHeightUtils',RIDE_HEIGHT_FILE)
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
        
        # Calculate wheel movements in tunnel coordinate system
        # In vehicle coordinates: wheel positions relative to reference
        # In tunnel coordinates: reference moves relative to wheels (inverse transformation)
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
    for idx, row in rideHeights.iterrows():
        if row['point'] not in fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RUN_RH_POINTS']:
            continue
        newCaseName = "%s_%s" % (baseCaseName, int(row['point']))
        newCaseDir = os.path.join(os.getcwd(), newCaseName)
        if not os.path.exists(newCaseDir):
            os.makedirs(newCaseDir)
        
        #write new caseSetup file with updated pitch and roll angles for domain
        newCaseSetupPath = os.path.join(newCaseDir,'caseSetup')
        updatedFullCaseSetupDict = fullCaseSetupDict.copy()
        updatedFullCaseSetupDict['BC_SETUP']['DOMAIN_PITCH'] = float(row['tunnel_pitch'])
        updatedFullCaseSetupDict['BC_SETUP']['DOMAIN_ROLL'] = float(row['tunnel_roll'])
        currentREFCOR = updatedFullCaseSetupDict['BC_SETUP']['REFCOR']
        updatedFullCaseSetupDict['BC_SETUP']['REFCOR'] = [currentREFCOR[0], currentREFCOR[1], str(float(currentREFCOR[2]) + row['tunnel_heave'])]
        writeToRHCaseSetup(updatedFullCaseSetupDict,newCaseSetupPath)


        print('\t\tCreated case: %s with updated pitch: %s, roll: %s, heave: %s' % (newCaseName, row['tunnel_pitch'], row['tunnel_roll'], row['tunnel_heave']))
    

        
        


def writeToRHCaseSetup(writeCaseSetupDict,newCasePath,output='default'):
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
    
    try:
        #caseSetupConfig.read_file(open("%s/%s/caseSetup" % (path,case)))
        if output == 'default':
            with open(newCasePath,'w') as caseSetupFile:
                writeConfig.write(caseSetupFile)
        else:
            with open(newCasePath,'w') as caseSetupFile:
                writeConfig.write(caseSetupFile)
    except:
        print('\nERROR! Unable to write caseSetup!')
        
    if updateCaseSetupFlag == True:
        sys.exit('\nWARNING: caseSetup has been updated with default values, please check caseSetup and rerun.')
                

def calculateRollAngles(dl,dr,bw):
    alpha = math.degrees(math.asin((dl-dr)/bw))
    d_center = (dl+dr)/2
    return alpha, d_center

def calculatePitchAngles(dzf,dzr,bl,REFCOR):
    theta = math.degrees(math.asin((dzf-dzr)/bl))
    dz_wb_center = ((dzf+dzr)/2)-float(REFCOR[2])

    return theta, dz_wb_center




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