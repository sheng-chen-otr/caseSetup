import os
import sys
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import argparse
import glob
import configparser
import scipy.stats as st
from scipy.optimize import curve_fit
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)


'''
	Plots the integral forces and allows for saving it or comparing against other cases

'''

#getting case information
casePath = os.getcwd() #path of directory this script is being run in
caseName = casePath.split('/')[-1] #directory name


parser = argparse.ArgumentParser(
                    prog='Force Plotter V1.0',
                    description='Plots the integral force development in your case',
                    epilog='Read the documentation for more information on how to use this')
parser.add_argument('-p','--plotData',default=['Cd','Cl','CoP'],choices = ['Cd','Cl','CoP','Cl(f)','Cl(r)','Cs(f)','Cs(r)'])
parser.add_argument('-t','--trial',default = [caseName],nargs='+',help='selects the trials to plot, if plotting current ONLY trial do not use this flag')
parser.add_argument('-s','--savePlots',action='store_true',help='saves all plots')
parser.add_argument('--skipStats',action='store_true',help='skips calculating statistics')
parser.add_argument('--saveFormat',default = 'png',help = 'sets the format for the saved plots, by default [png]',
					choices = ['png','eps','jpeg']) 
parser.add_argument('--avgTime',help = 'manually set averaging start time',type=float) 


args = parser.parse_args()


def main():
	global caseLoc

	casePathDict, caseLoc = setCasePaths(args.trial)

	casePathDict = getCaseData(casePathDict)
	casePathDict = makePandasArrays(casePathDict)
	plotData(casePathDict)



def setCasePaths(inputCaseList):
	casePathDict = {}
	if casePath.split('/')[-2].lower() == 'cases':
		print('\tExecuting in trial directory...')
		searchPath = os.path.join('/',*casePath.split('/')[:-1])
		caseLoc = 'inTrial'
	elif casePath.split('/')[-1].lower() == 'cases':
		print('\tExecuting in CASES directory...')
		searchPath = os.path.join('/',*casePath.split('/'))
		caseLoc = 'outTrial'
	else:
		sys.exit('ERROR! Unable to plot, script not executed in CASES directory or in trial directory!')
	
	print('\tLooking for cases...')
	for case in inputCaseList:
		if os.path.isdir(searchPath + '/' + case):
			casePathDict[case] = {}
			casePathDict[case]['path'] = searchPath + '/' + case
			print('\t\tFound %s' % (case))
		else:
			print('\t\tUnable to find %s, skipping...' % (case))

	if len(list(casePathDict.keys())) < 1:
		sys.exit('ERROR! Unable to find any requested cases!')
	
	return casePathDict, caseLoc
	

def getCaseData(casePathDict):
	'''
		Looking for case data, finds if the case has been restarted or not, if it has it will add the restart to it
		If there exists multiple coefficient.dat files it will use the latest one only
		The different times will be appended together starting from the oldest to the newest time

	'''
	global caseSetupConfig
	

	for case in casePathDict.keys():
		casePath = casePathDict[case]['path']
		UREF,LREF,WREF,CREF,FREF,AVGT,caseSetupConfig = getRef(casePath)
		postProPath = os.path.join(casePath,'postProcessing/all/*')

		#checking to see if there was a restart at a different time, will join data if so

		forceTimeList = glob.glob(postProPath)
		forceTimeList.sort(key=os.path.getmtime)
		casePathDict[case]['caseTimes'] = {}
		casePathDict[case]['avgStart'] = AVGT


		for caseTime in forceTimeList:
			coefficientsList = glob.glob(os.path.join(caseTime,'coefficient*.dat'))
			coefficientsList.sort(key=os.path.getmtime)
			casePathDict[case]['caseTimes'][caseTime.split('/')[-1]] = coefficientsList[-1]


	return casePathDict

			
def makePandasArrays(casePathDict):

	print('\n\tConsolidating data in each case...')


	for case in casePathDict.keys():
		print('\t\tCase times in %s: ' % (case))
		avgt = casePathDict[case]['avgStart']
		
		if len(casePathDict[case]['caseTimes'].keys()) < 0:
			sys.exit('ERROR! No data availalbe for %s' % (case))
	
		n = 0
		for caseTime in casePathDict[case]['caseTimes'].keys():
			print('\t\t\t%s' % (caseTime))
			
			if n == 0:
				forceData = pd.read_csv(casePathDict[case]['caseTimes'][caseTime],delimiter='\t',skiprows = 12)
			else:
				newForceData = pd.read_csv(casePathDict[case]['caseTimes'][caseTime],delimiter='\t',skiprows = 12)
				forceData = pd.concat((forceData,newForceData),axis=0,ignore_index = True)

			


			
			
				
			
			n = n + 1

		#make new columns
		newColumn = []
		for column in forceData.columns:
			newColumn.append(column.replace('\t','').replace(' ',''))
		forceData.columns = newColumn
		
		casePathDict[case]['data'] = forceData
		#calculate center of pressure
		forceData['CoP'] = forceData['Cl(f)'].div(forceData['Cl']).mul(100)
		maxCurrentTime = max(forceData['#Time'])



		if args.avgTime != None:
			avgt = args.avgTime
			
		if maxCurrentTime < avgt:
			avgForceData = forceData
		else:
			avgForceData = forceData[forceData['#Time']>=avgt].copy()
			
		time = avgForceData['#Time'].to_numpy()
		for var in args.plotData:
			if args.skipStats == False:
				data = avgForceData[var].to_numpy()		
				print('\t\t\tCalculating statistics for %s' % (var))
				if len(time) > 100:
					sampleMean,ci,errorPercent = calcCI(time, data,100)


					if sampleMean == 1:
						continue
						print('\t\t\tInsufficient averaging time to calculate statistics...')
					else:
						avgForceData[var+'SampleMean'] = sampleMean 
						avgForceData[var+'CI'] = ci 
						avgForceData[var+'Error'] = errorPercent
				else:
					print('In-sufficient time to calculate statistics!')
			
			casePathDict[case]['avgData'] = avgForceData

	return casePathDict
		

def plotData(casePathDict):
	labelDict = {'Cl':{'label':'$C_L$','title':'$C_L$'},
				 'Cd':{'label':'$C_D$','title':'$C_D$'},
				 'Cl(f)':{'label':'$C_{L,f}$','title':'$C_{L,f}$'},
				 'Cl(r)':{'label':'$C_{L,r}$','title':'$C_{L,r}$'},
				 'Cs(f)':{'label':'$C_{S,f}$','title':'$C_{S,f}$'},
				 'Cs(r)':{'label':'$C_{S,r}$','title':'$C_{S,r}$'},
				 'CoP':{'label':'Percent Front','title':'Center of Pressure'},
				}
	
	for var in args.plotData:
		fig1,ax = plt.subplots(figsize=[7,4],frameon=True)
		
		caseMeans = []
		for case in casePathDict.keys():
			rawData = casePathDict[case]['data']
			avgData = casePathDict[case]['avgData']
			data = rawData[['#Time',var]]
			if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
				caseMeans.append(np.mean(avgData[var][-20:])) #using this to calculate the range for scaling axis
			else:
				caseMeans.append(np.mean(avgData[var][-20:])*2) #using this to calculate the range for scaling axis
			if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
				ax.plot(rawData['#Time'],rawData[var]*2,label=case,alpha=0.5)
			else:
				ax.plot(rawData['#Time'],rawData[var],label=case,alpha=0.5)
			if args.skipStats == False:
				#try:
				for statvar in ['SampleMean','CI']:
					if statvar == 'SampleMean':
						statvarlabel = 'Fwd. Avg'
						if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
							ax.plot(avgData['#Time'],avgData[var+statvar]*2,label=case + '(%s)'% (statvarlabel),linestyle = '--')
						else:
							ax.plot(avgData['#Time'],avgData[var+statvar],label=case + '(%s)'% (statvarlabel),linestyle = '--')
					# elif statvar == 'CI':
					# 	statvarlabel = statvar
					# 	ax.fill_between(avgData['#Time'],avgData[var + 'SampleMean'] - avgData[var+statvar],avgData[var + 'SampleMean'] + avgData[var+statvar],label=case + '(%s)'% (statvarlabel),alpha = 0.5)
						
				#except:
					#print('WARNING! Unable to plot statistics data!')
						

		caseMeans = np.array(caseMeans)
		ax.set_ylim(min(caseMeans)-(np.mean(caseMeans)*0.1),max(caseMeans)+(np.mean(caseMeans)*0.1))
		ax.set_ylabel(labelDict[var]['label'])
		ax.set_title(labelDict[var]['title'])
		ax.set_xlabel('Time (s)')
		plt.legend()
		plt.savefig('%s_forceHistory_%s.%s' % ('_'.join(args.trial),var,args.saveFormat),dpi = 300,bbox_inches='tight')



def calcCI(time,data,nchunks = 100):
	if len(data) < nchunks:
		
		return 1,1,1,1,1

	#calculating the 95% confidence interval for the data to gauge convergence
	ciArrayMean = []
	errorArrayMean = []
	sampleMeanArrayMean = []
	lowerBoundArrayMean = []
	upperBoundArrayMean = []
	ciArrayMean = []
	errorPercentArrayMean = []

	ciArrayStd = []
	errorArrayStd = []
	sampleStdArrayStd = []
	lowerBoundArrayStd = []
	upperBoundArrayStd = []
	ciArrayStd = []
	errorPercentArrayStd = []

	totalMean = np.mean(data)
	totalStd = np.std(data)
	#break data into chunks, about 30 chunks
	for n in range(len(data)):
		if n < 1:
			samples = data[:n+1]
			std = data[:n+1]
			errorChunk = np.sqrt(np.mean(data[:n+1]-totalMean)**2)/totalMean
		else:
			samples = []
			std = []
			errorChunk = []
			#chunking data
			chunks = np.array_split(data[:n],nchunks)
			#calculate mean of each chunk
			for chunk in chunks:
				if len(chunk) > 0:
					samples.append(float(np.mean(chunk)))
					errorChunk.append(np.sqrt(np.mean((np.mean(chunk)-totalMean)**2))/totalMean)
					std.append(float(np.std(chunk)))
			samples = np.array(samples)
			std = np.array(std)
			errorChunk = np.array(errorChunk)

		#errorMean = np.sqrt(np.mean((samples-totalMean)**2))/totalMean
		errorMean = errorChunk
		
		errorStd = np.sqrt(np.mean((std-totalStd)**2))/totalStd

		sampleMean = np.mean(samples)
		sampleStd = np.std(samples)

		errorArrayMean.append(errorMean)
		lowerBoundMean = totalMean/(1+2*errorMean)
		upperBoundMean = totalMean/(1-2*errorMean)
		lowerBoundArrayMean.append(lowerBoundMean)
		upperBoundArrayMean.append(upperBoundMean)
		ciArrayMean.append((upperBoundMean-lowerBoundMean)/2)
		errorPercentArrayMean.append(abs(100*((upperBoundMean-lowerBoundMean)/2)/totalMean))
		sampleMeanArrayMean.append(sampleMean)

		errorArrayStd.append(errorStd)
		lowerBoundStd = totalMean/(1+2*errorStd)
		upperBoundStd = totalMean/(1-2*errorStd)
		lowerBoundArrayStd.append(lowerBoundStd)
		upperBoundArrayStd.append(upperBoundStd)
		ciArrayStd.append((upperBoundMean-lowerBoundStd)/2)
		errorPercentArrayStd.append(abs(100*((upperBoundStd-lowerBoundStd)/2)/totalStd))



	

	return sampleMeanArrayMean,ciArrayMean,errorPercentArrayMean






def getRef(casePath):
    '''
        gets all the reference values based on the fullCaseSetupDict file
    '''
    print('\n\tReading fullCaseSetupDict from case directory...')
    caseSetupConfig = configparser.ConfigParser()
    caseSetupConfig.optionxform = str
    #tries to read the caseSetup file, if it can't find it, it throws an error
    try:
        caseSetupConfig.read_file(open(casePath+'/caseSetup'))
    except Exception as e:
        print('\n\n\nERROR! fullCaseSetupDict invalid!')
        print(e)
        exit()

    caseSetupConfigSections = caseSetupConfig.sections()    

    UREF = float(caseSetupConfig['BC_SETUP']['INLET_MAG'])
    LREF = float(caseSetupConfig['BC_SETUP']['REFLEN'])
    WREF = float(caseSetupConfig['BC_SETUP']['REFWIDTH'])
    CREF = np.array(caseSetupConfig['BC_SETUP']['REFCOR'].split(' ')).astype('float')
    FREF = np.array(caseSetupConfig['BC_SETUP']['FRT_WH_CTR'].split(' ')).astype('float')
    AVGT = float(caseSetupConfig['GLOBAL_CONTROL']['AVGSTART'])

    #print('\t\tReference Values:')
    #print('\t\t\tUREF: %1.3g, LREF: %1.3g, WREF: %1.3g, CREF: [%1.3g %1.3g %1.3g] , FREF: [%1.3g %1.3g %1.3g]' % (UREF,LREF,WREF,CREF[0],CREF[1],CREF[2],FREF[0],FREF[1],FREF[2]))
    
    return UREF,LREF,WREF,CREF,FREF,AVGT,caseSetupConfig




























main()