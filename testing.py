import numpy as np
import shutil

def main():

    geomFile = 'ROTA-fr-wh-lhs.stl'
    #getBoundingBoxOBJ(geomFile)
    getBoundingBoxSTL(geomFile)
    
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
                    line = line.replace('  ',',').replace(' ',',')
                    xCoords.append(line.split(',')[-3])
                    yCoords.append(line.split(',')[-2])
                    zCoords.append(line.split(',')[-1])
    else:
        with open(stlPath, 'rt') as stlFile:
            #initializing the arrays
            xCoords = []
            yCoords = []
            zCoords = []
            
            for line in stlFile:
                
                if 'vertex' in line:
                    line = line.replace('  ',',').replace(' ',',')
                    print(line.split(',')[-3])
                    xCoords.append(line.split(',')[-3])
                    yCoords.append(line.split(',')[-2])
                    zCoords.append(line.split(',')[-1])
    
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
    
    print(bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ)  
    sys.exit()
    
    
    return bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ  
    
    

main()