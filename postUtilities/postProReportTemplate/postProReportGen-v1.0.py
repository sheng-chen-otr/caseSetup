import os
import sys
import re



path = os.path.split(os.getcwd())[0] 		
case = os.path.split(os.getcwd())[1]
job = os.path.dirname(path)







postProImages = [".".join(f.split(".")[:-1]) for f in os.listdir("%s/%s/postProcessing/images/" % (path,case))]


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]



def getCpImage(postProImages,path,case):
	cPList = []
	for image in postProImages:

		if "_cP_" in image and not "Slice" in image:
			cPList.append(image)
	
	cPListSorted = sorted(cPList)



	return cPListSorted


def getWssImage(postProImages,path,case):
	wssList = []
	for image in postProImages:

		if "_Wss_" in image and not "Slice" in image:
			wssList.append(image)
	
	wssListSorted = sorted(wssList)

	return wssListSorted


def getCtpXSlice(postProImages,path,case):
	ctpXSliceList = []
	for image in postProImages:

		if "_ctp_xSlice" in image:
			ctpXSliceList.append(image)
	
	ctpXSliceList.sort(key=natural_keys)

	return ctpXSliceList



cPListSorted = getCpImage(postProImages,path,case)
wssListSorted = getWssImage(postProImages,path,case)
ctpXSliceList = getCtpXSlice(postProImages,path,case)



def makeCpReport(cPListSorted):

	imageList = cPListSorted

	numcp = len(imageList)

	texText = []

	n = 0

	for image in imageList:
		template = open("postProReport/cP_plots_template.tex","r")

		if n == 0:
			nextLine = imageList[n + 1]
			prevLine = imageList[numcp - 1]
			
		elif n == numcp - 1:
			nextLine = imageList[0]
			prevLine = imageList[n - 1]
			
		else:
			nextLine = imageList[n + 1]
			prevLine = imageList[n - 1]
			


		for line in template:
			
			if "TRIAL" in line:
				reLine = line.replace("TRIAL",case)
				texText.append(reLine)

			elif "IMAGEPATH" in line:
				reLine = line.replace("IMAGEPATH","%s/%s/postProcessing/images/%s.png" % (path,case,image))
				texText.append(reLine)
			
			elif "IMAGELABEL" in line:
				reLine = line.replace("IMAGELABEL",image)
				texText.append(reLine)

			elif "NEXTLABEL" in line:
				reLine = line.replace("NEXTLABEL",nextLine)
				texText.append(reLine)

			elif "PREVLABEL" in line:
				reLine = line.replace("PREVLABEL",prevLine)
				texText.append(reLine)
			else:
				reLine = line
				texText.append(reLine)



		n = n + 1


	with open('%s/%s/postProReport/cP_plots.tex' % (path,case),'w') as texFile:
		texFile.write('\n'.join(texText))




makeCpReport(cPListSorted)












