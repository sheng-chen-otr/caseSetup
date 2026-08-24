import sys
import os
import pandas as pd
import argparse as argparse
import matplotlib.pyplot as plt
from itertools import repeat
import numpy as np
from PIL import Image
import cv2
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
from bcParser_v1_0 import bcParser

path = os.path.split(os.getcwd())[0]
case = os.path.split(os.getcwd())[1]

 

parser = argparse.ArgumentParser(prog='Force Development Plotter',description='Will plot force development and forward average values if sufficient iterations available.\
                                                                    If there are insufficient time steps for averaging, the forward average line will be equal\
                                                                    to the raw force integral line.')
parser.add_argument("-t","--caseNames", metavar='trial',nargs='+',
                    help='Additional trials to be plotted, if doing current trial, leave blank!')
                    
parser.add_argument('-a',"--forwardAvg", metavar='avgStart', type=float, 
                    help='Time to start forward averaging')

parser.add_argument("-s","--save", action="store_true",
                    help='Saves all plots (regardless if -p options are used) to current trial.')

parser.add_argument("-n","--noShow", action="store_true",
                    help='Does not show plots.')    

parser.add_argument("-i","--noImage", action="store_true",
                    help='Does not show geometry images for all trials.')                 
                    
args = parser.parse_args()


save = args.save
noShow = args.noShow
noImage = args.noImage

    

if args.caseNames is None:
    caseNames = [case]
else:
    caseNames = args.caseNames
    caseNames = [case] + caseNames


#check if all the cases listed can be found can be found, and checking that it's type is the same

def findCase(caseName):
    global inletMag,lastTime,yaw,wheelRotation,simType,turbModel
    if os.path.isdir("%s/%s" % (path,caseName)):
        print("     Found trial%s..." % (caseName))
        [inletMag,lastTime,yaw,wheelRotation,simType,turbModel] = bcParser(path,caseName)
        print("         Simulation Type: %s\n" % (simType))
        
    else:
        sys.exit("  Trial%s not found! Check your entry!\n" % (caseName))
    
def check(lst):
    # Use the repeat function to generate an iterator that returns the first element of the list repeated len(lst) times
    repeated = list(repeat(lst[0], len(lst)))
    # Compare the repeated list to the original list and return the result
    if not repeated == lst:
        sys.exit("Simulation type not consistent! Exiting!")

def getGeomImage(path,caseName):
    imagePath = "%s/%s/postProcessing/images/Geom_Surface/%s_Geom_Surface_Left.png" % (path,caseName,caseName)
    if not os.path.isfile(imagePath):
        #sys.exit("     Cannot find geometry image for trial%s, exiting..." % (caseName))
        print("     Cannot find geometry image for trial%s, skipping for this trial..." % (caseName))
        SKIPIMAGE = True
        ROI = []
        return ROI, SKIPIMAGE
    elif noImage == True:
        #sys.exit("     Cannot find geometry image for trial%s, exiting..." % (caseName))
        print("     No Image option selected, will not show geometry image...")
        SKIPIMAGE = True
        ROI = []
        return ROI, SKIPIMAGE
    else:
        ROI = convertImage(imagePath)
        SKIPIMAGE = False
        return ROI, SKIPIMAGE
        
    
    
def convertImage(imagePath):
    # Load image, convert to grayscale, Gaussian blur, Otsu's threshold
    image = cv2.imread(imagePath)
    original = image.copy()
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    #gray = image
    #blur = cv2.GaussianBlur(gray, (3,3), 0)
    thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)[1]

    # Obtain bounding rectangle and extract ROI
    x,y,w,h = cv2.boundingRect(thresh)
    cv2.rectangle(image, (x, y), (x + w, y + h), (36,255,12), 2)
    ROI = original[y:y+h, x:x+w]

    # Add alpha channel
    b,g,r = cv2.split(ROI)
    alpha = np.ones(b.shape, dtype=b.dtype) *100
    ROI = cv2.merge([b,g,r,alpha])

    #cv2.imshow('thresh', thresh)
    #cv2.imshow('image', image)
    #cv2.imshow('ROI', ROI)
    #cv2.waitKey()
    return ROI
     

def importBinData(path,caseName):
    xCoords = []
    forceHeader = []
    forceCoeffs = []
    [inletMag,lastTime,yaw,wheelRotation,simType,turbModel] = bcParser(path,caseName)
    exportTime=lastTime
    if not os.path.isdir("%s/%s/postProcessing/binForceCoeffs/" % (path,caseName)):
        sys.exit("Bin force directory not found! The run might not be complete...")
    timeList = os.listdir("%s/%s/postProcessing/binForceCoeffs/" % (path,caseName))
    if not lastTime in timeList:
        exportTime = max(timeList)
        #print(exportTime)
        print("End time specified in controlDict does not match available times. Using lastest available time %s" % (exportTime))
    
   
    if not os.path.isdir("%s/%s/postProcessing/binForceCoeffs/%s/" % (path,caseName,exportTime)):
        sys.exit("Bin forces export might not have been run successfully for trial %s, please check your log files!" % (caseName))
    
    coeffList = os.listdir("%s/%s/postProcessing/binForceCoeffs/%s/" % (path,caseName,exportTime))
    
    if len(coeffList) < 1:
        sys.exit("Bin forces export might not have been run successfully for trial %s, please check your log files!" % (caseName))
    
    for coeff in coeffList:
        filePath = "%s/%s/postProcessing/binForceCoeffs/%s/%s" % (path,caseName,exportTime,coeff)
        
        with open(r"%s/%s/postProcessing/binForceCoeffs/%s/%s" % (path,caseName,exportTime,coeff), 'r') as fp:
            lines = len(fp.readlines())
            
            
        if lines > 10 and lastTime in coeff:
            binCoeff = coeff
            break
        elif lines > 10:
            binCoeff = coeff
    

    if binCoeff == "":
        sys.exit("None of the coefficients files in trial %s have enough data!" % (caseName))
                
                
    
        
    with open(r"%s/%s/postProcessing/binForceCoeffs/%s/%s" % (path,caseName,exportTime,binCoeff), 'r') as fp:

        for line in fp:
            if "x co-ords" in line:
                xCoords =  line.replace("#","").replace(":","").replace("x","").replace("co-ords","")
                xCoords = xCoords.split()
                xCoords = [float(x) for x in xCoords]
            elif "# Time" in line:
                forceHeader = line.replace("#","").replace("Time","")
                forceHeader = forceHeader.split()
                #print(len(forceHeader))
            elif exportTime in line:
                #print(lastTime)
                forceCoeffs = line.replace(exportTime,"")
                
                forceCoeffs = forceCoeffs.split()
                #print(len(forceCoeffs))
                forceCoeffs = [float(x) for x in forceCoeffs]
                
               
               
    dataFrameZ = pd.DataFrame(forceCoeffs, index = forceHeader,columns = ['cl'])
    dataFrameX = pd.DataFrame(forceCoeffs, index = forceHeader,columns = ['cd'])
    
    dataFrameZ = dataFrameZ[dataFrameZ.index.str.contains("internal|patch|_x|_y") == False]   
    dataFrameX = dataFrameX[dataFrameX.index.str.contains("internal|patch|_z|_y") == False]
    
    dataFrameZ["xCoords"] = xCoords
    dataFrameX["xCoords"] = xCoords
    
    
    dataFrameZ = dataFrameZ.set_index(['xCoords'])
    dataFrameX = dataFrameX.set_index(['xCoords'])
    
    return dataFrameZ, dataFrameX
    
simulationArray = []
lastTimeArray = []    
for i in caseNames:
    print("Checking for trial%s...\n" % (i))
    findCase(i)
    simulationArray = simulationArray + [simType]
    lastTimeArray = lastTimeArray + [lastTime]
    
#check(simulationArray)


allClBin = pd.DataFrame()
allClBinMin = []
allClBinMax = []

allCdBin = pd.DataFrame()
allCdBinMin = []
allCdBinMax = []

allxcoordBinMin = []
allxcoordBinMax = []



#plt.style.use('bmh')
clbinfig = plt.figure(figsize=(10, 5))
clbinax = clbinfig.add_subplot(111, frame_on=True)
clbinax.set_title("Cl Development")

cdbinfig = plt.figure(figsize=(10, 5))
cdbinax = cdbinfig.add_subplot(111, frame_on=True)
cdbinax.set_title("Cd Development")

#plotting the bin data first
for i in caseNames:
    dataFrameZ,dataFrameX = importBinData(path,i)
    allClBin[i+" - distance"] = dataFrameZ.index
    if "half" in i:
        allClBin[i] = dataFrameZ.cl.to_numpy()*2
        clMin = dataFrameZ.cl.to_numpy().min()*2
        clMax = dataFrameZ.cl.to_numpy().max()*2
        #print(allClBin.to_string())
    else:
        allClBin[i] = dataFrameZ.cl.to_numpy()
        clMin = dataFrameZ.cl.to_numpy().min()
        clMax = dataFrameZ.cl.to_numpy().max()
    
    allClBinMin.append(clMin)
    allClBinMax.append(clMax)
    
    allCdBin[i+" - distance"] = dataFrameX.index
    if "half" in i:
        allCdBin[i] = dataFrameX.cd.to_numpy()*2
        cdMin = dataFrameX.cd.to_numpy().min()*2
        cdMax = dataFrameX.cd.to_numpy().max()*2
    else:
        allCdBin[i] = dataFrameX.cd.to_numpy()
        cdMin = dataFrameX.cd.to_numpy().min()
        cdMax = dataFrameX.cd.to_numpy().max()
        
    allCdBinMin.append(cdMin)
    allCdBinMax.append(cdMax)
        
    xMin = dataFrameZ.index.to_numpy().min()
    xMax = dataFrameZ.index.to_numpy().max()
    allxcoordBinMin.append(xMin)
    allxcoordBinMax.append(xMax)
    
    allCdBin.plot(x=i+" - distance",y=i,ax=cdbinax)
    allClBin.plot(x=i+" - distance",y=i,ax=clbinax)
    

#displaying the figures    
count = 0    
for i in caseNames:
    ROI,SKIPIMAGE = getGeomImage(path,i)
    #if there is no image found, skip this this part
    if SKIPIMAGE == True:
        continue
        
    h,w,c = ROI.shape
    r = h/w
    xMin = allxcoordBinMin[count]
    xMax = allxcoordBinMax[count]
    clMax = max(allClBinMax)
    clMin = min(allClBinMin)
    
    plotR = (xMax-xMin)/(clMax-clMin)
    rFactor = r/plotR
    clbinax.imshow(ROI,origin='upper',extent = [xMin,xMax,clMin,clMin+(xMax-xMin)*r])
    cdbinax.imshow(ROI,origin='upper',extent = [xMin,xMax,cdMin,cdMin+(xMax-xMin)*r])
    count += 1
   
   
    
    


#clbinax.autoscale()
clbinax.margins(0,0)
clbinMin = min(allClBinMin)
clbinMax = max(allClBinMax)
posBinMin = min(allxcoordBinMin)
postBinMax = max(allxcoordBinMax)
clbinax.set_xlim([posBinMin,postBinMax])

clbinax.set_xlabel('Distance (m)')
clbinax.set_ylabel('Cl')
clbinax.legend(loc='upper left')
clbinax.grid(visible=True,which='both',axis='both')
clbinax.xaxis.set_minor_locator(AutoMinorLocator())
clbinax.yaxis.set_minor_locator(AutoMinorLocator())
clbinax.tick_params(which='both', width=2)
clbinax.tick_params(which='major', length=7)
clbinax.tick_params(which='minor', length=4, color='b')



cdbinax.autoscale()
cdbinax.margins(0,0)
cdbinMin = min(allCdBinMin)
cdbinMax = max(allCdBinMax)
posBinMin = min(allxcoordBinMin)
postBinMax = max(allxcoordBinMax)
cdbinax.set_xlim([posBinMin,postBinMax])

cdbinax.set_xlabel('Distance (m)')
cdbinax.set_ylabel('Cd')
cdbinax.legend(loc='upper left')
cdbinax.grid(visible=True,which='both',axis='both')
cdbinax.xaxis.set_minor_locator(AutoMinorLocator())
cdbinax.yaxis.set_minor_locator(AutoMinorLocator())
cdbinax.tick_params(which='both', width=2)
cdbinax.tick_params(which='major', length=7)
cdbinax.tick_params(which='minor', length=4, color='b')


plotList = [cdbinfig,clbinfig]
plotNames = ['cd','cl']
if save == True:
    for figure,plotName in zip(plotList,plotNames):
        figure.savefig("%s-development_%s.png" % (plotName,'_'.join(caseNames)),dpi=300,bbox_inches='tight')


if not noShow == True:
    plt.show()















