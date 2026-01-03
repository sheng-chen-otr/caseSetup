import os
import sys
import re
import numpy as np


path = os.path.split(os.getcwd())[0] 		
case = os.path.split(os.getcwd())[1]
job = os.path.split(os.path.dirname(path))[1]
scriptPath = os.path.dirname(os.path.abspath(__file__))
comp = False

templatePath = os.path.join(scriptPath,'postProReportTemplate')
logoPath = os.path.join(scriptPath,'postProReportTemplate','logos','otrLogo_greyOnWhite.png')


if len(sys.argv) > 1:
	case2 = sys.argv[1]
	comp = True


if not os.path.isdir("%s/%s/postProReport" % (path,case)):
	print("No postProReport directory found in trial%s, making..." % (case))
	os.system('mkdir %s/%s/postProReport' % (path,case))
	print("Copying postProReport_template to postProReport directory...")
	os.system("cp %s/postProReport_template.tex %s/%s/postProReport/postProReport.tex" % (templatePath,path,case))

	# if comp == True:
	# 	print("No postProReport directory found in trial%s, making..." % (case2))
	# 	os.system('mkdir %s/%s/postProReport' % (path,case2))
	# 	print("Copying postProReport_template to postProReport directory...")
	# 	os.system("cp %s/postProReport_template.tex %s/%s/postProReport/postProReport.tex" % (templatePath,path,case2))

else:
	print("postProReport directory found trial%s,, proceeding..." % (case))
	print("Copying postProReport_template to postProReport directory...")
	os.system("cp %s/postProReport_template.tex %s/%s/postProReport/postProReport.tex" % (templatePath,path,case))

	# if comp == True:
	# 	print("postProReport directory found trial%s, proceeding..." % (case2))
	# 	#print("Copying postProReport_template to postProReport directory...")
	# 	#os.system("cp %s/postProReport_template.tex %s/%s/postProReport/postProReport.tex" % (templatePath,path,case2))



def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


objectTypes = ['Geom','CpMean','CptMean','isoQ','Wss','UMeanNear']
objectLists = ['Geom','Cp','IsoCtp','IsoQ','Wss','UMeanNear']
objectArray = np.empty((10,6),dtype=object)
objectArray2 = np.empty((10,6),dtype=object)
objectListSlice = ['cp_xSlice','cp_ySlice','cp_zSlice','ctp_xSlice','ctp_ySlice','ctp_zSlice','ctp_LIC_xSlice','ctp_LIC_ySlice','ctp_LIC_zSlice','percentU_xSlice','percentU_ySlice','percentU_zSlice']


def titleText(path,job,case,logoPath):
	texTitle = []
	mainTemplate = open("%s/%s/postProReport/postProReport.tex" % (path,case),"r")
	for line in mainTemplate:

		
		if "TITLE" in line:
			if comp == True:
				reLine = line.replace("TITLE",("\\verb|%s - Trial%s and Trial%s|" % (job,case,case2)))

			else:
				reLine = line.replace("TITLE",("\\verb|%s - Trial%s|" % (job,case)))

			texTitle.append(reLine)
		elif "LOGO" in line:
			reLine = line.replace("LOGO",(logoPath))
			texTitle.append(reLine)
		else:
			reLine = line
			texTitle.append(reLine)

	with open('%s/%s/postProReport/postProReport.tex' % (path,case),'w') as texFile:
		texFile.write('\n'.join(texTitle))






def getImage(path,case,objectType):
	List = []
	postProImages = [".".join(f.split(".")[:-1]) for f in os.listdir("%s/%s/postProcessing/images/" % (path,case))]
	for image in postProImages:

		if objectType in image and not "Slice" in image:
			List.append(image)
	
	listSorted = sorted(List)
	listSorted = np.array(listSorted)




	return listSorted

#titleText(path,job,case,logoPath)

n = 0
for obj in objectTypes:

	imList = getImage(path,case,obj)
	#print(imList)
	objectArray[:,n] = imList
	if comp == True:
		imList2 = getImage(path,case2,obj)
		objectArray2[:,n] = imList2

	n = n + 1


objectItemLen = len(objectArray[0,:])

def makeReport(postProList,comp,postProList2,objectItemCount):

	imageList = postProList[:,objectItemCount]
	num1 = len(imageList)

	if comp == True:
		imageList2 = postProList2[:,objectItemCount]
		num2 = len(imageList2)

		if num1 != num2:
			sys.exit("The number of images don't match in both trials!")


	texText = []
	texText.append("\\section{%s Plots}" % (objectLists[objectItemCount]))
	#texText.append("\\subsection{\\url{Trial_%s}}" % (case))

	n = 0

	for image in imageList:
		template = open("%s/plots_template.tex" % (templatePath),"r")

		if n == 0:
			nextLine = imageList[n + 1]
			prevLine = imageList[num1 - 1]
			
		elif n == num1 - 1:
			nextLine = imageList[0]
			prevLine = imageList[n - 1]
			
		else:
			nextLine = imageList[n + 1]
			prevLine = imageList[n - 1]
			


		for line in template:
			
			if "TRIAL" in line:
				reLine = line.replace("TRIAL",case)
				reLine = reLine.replace("PLOTNAME",objectLists[objectItemCount])
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
				if comp == True:
					nextTrial = "\\hyperlink{%s}{\\beamergotobutton{NEXT TRIAL}}" % (imageList2[n])
					texText.append(nextTrial)
			else:
				reLine = line
				texText.append(reLine)



		n = n + 1

	n = 0

	if comp == True:
		#texText.append("\\subsection{\\url{Trial_%s}}" % (case2))
		for image in imageList2:
			template = open("%s/plots_template.tex" % (templatePath),"r")

			if n == 0:
				nextLine = imageList2[n + 1]
				prevLine = imageList2[num2 - 1]
				
			elif n == num2 - 1:
				nextLine = imageList2[0]
				prevLine = imageList2[n - 1]
				
			else:
				nextLine = imageList2[n + 1]
				prevLine = imageList2[n - 1]
				


			for line in template:
				
				if "TRIAL" in line:
					reLine = line.replace("TRIAL",case2)
					print(n)
					reLine = reLine.replace("PLOTNAME",objectLists[objectItemCount])
					texText.append(reLine)
					

				elif "IMAGEPATH" in line:
					reLine = line.replace("IMAGEPATH","%s/%s/postProcessing/images/%s.png" % (path,case2,image))
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
					if comp == True:
						nextTrial = "\\hyperlink{%s}{\\beamergotobutton{NEXT TRIAL}}" % (imageList[n])
						texText.append(nextTrial)
				else:
					reLine = line
					texText.append(reLine)



			n = n + 1


	with open('%s/%s/postProReport/%s_plots.tex' % (path,case,objectLists[objectItemCount]),'w') as texFile:
		texFile.write('\n'.join(texText))


for n in range(objectItemLen):
	
	makeReport(objectArray,comp,objectArray2,n)


def makeWssreport(wssListSorted):

	imageList = wssListSorted

	numcp = len(imageList)

	texText = []
	texText.append("\\section{WSS Plots}")

	n = 0

	for image in imageList:
		template = open("%s/wss_plots_template.tex" % (templatePath),"r")

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

	
	with open('%s/%s/postProReport/wss_plots.tex' % (path,case),'w') as texFile:
		texFile.write('\n'.join(texText))







print("Making pdf using pdfLatex...")
#print("pdflatex %s/%s/postProReport/postProReport.tex -output-directory=%s/%s/postProReport" % (path,case,path,case))
os.system("cd %s/%s/postProReport;pdflatex %s/%s/postProReport/postProReport.tex -output-directory=%s/%s/postProReport" % (path,case,path,case,path,case))

if comp == True:
	os.system("mv postProReport/postProReport.pdf trial%s_trial%s_postProReport.pdf" % (case,case2))
else:
	os.system("mv postProReport/postProReport.pdf trial%s_postProReport.pdf" % (case))
print("Done!")












