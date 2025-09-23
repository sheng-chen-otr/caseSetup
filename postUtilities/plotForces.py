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
from estimateStatisticalError import *
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)


def setCasePaths(inputCaseList,casePath):
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

			
def makePandasArrays(args,casePathDict):

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
		avgData = pd.DataFrame()
		for var in args.plotData:
			if args.skipStats == False:
				data = avgForceData[var].to_numpy()		
				print('\t\t\tCalculating statistics for %s' % (var))
				dt = time[-1]-time[-2]
				if len(time) > 10:
					statTime, avgDataArray,ciUpperArray,ciLowerArray  = calculateAcrossTime(time,data,dt=dt)
					avgData[var+'SampleMean'] = avgDataArray
					avgData[var+'UpperCI'] = ciUpperArray
					avgData[var+'LowerCI'] = ciLowerArray
					avgData[var+'StatTime'] = statTime

				else:
					print('In-sufficient time to calculate statistics!')
			
			casePathDict[case]['avgData'] = avgData

	return casePathDict
		

def plotData(args,caseLoc,casePathDict):
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
		caseMins = []
		caseMax = []
		for case in casePathDict.keys():
			rawData = casePathDict[case]['data']
			avgData = casePathDict[case]['avgData']
			data = rawData[['#Time',var]]
			avgStart = casePathDict[case]['avgStart']
			
			if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
				caseMeans.append(np.mean(data[var][data['#Time']>=avgStart])*2) #using this to calculate the range for scaling axis
				#caseMax.append(max(data[var][data['#Time']>=avgStart])*2)
				#caseMins.append(min(data[var][data['#Time']>=avgStart])*2)
			else:
				caseMeans.append(np.mean(data[var][data['#Time']>=avgStart])) #using this to calculate the range for scaling axis\
				#caseMax.append(max(data[var][data['#Time']>=avgStart]))
				#caseMins.append(min(data[var][data['#Time']>=avgStart]))
			if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
				ax.plot(rawData['#Time'],rawData[var]*2,label=case,alpha=0.5)
			else:
				ax.plot(rawData['#Time'],rawData[var],label=case,alpha=0.5)
			if args.skipStats == False:
				try:
					for statvar in ['SampleMean','CI']:
						if statvar == 'SampleMean':
							statvarlabel = 'Fwd. Avg'
							if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
								ax.plot(avgData[var+'StatTime'],avgData[var+statvar]*2,label=case + '(%s)'% (statvarlabel),linestyle = '--')
							else:
								ax.plot(avgData[var+'StatTime'],avgData[var+statvar],label=case + '(%s)'% (statvarlabel),linestyle = '--')
						elif statvar == 'CI':
							statvarlabel = statvar
							ci = (avgData[var+'UpperCI']-avgData[var + 'LowerCI'])/2
							
							if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
								ax.fill_between(avgData[var+'StatTime'],avgData[var+'SampleMean']*2+ci,
											avgData[var+'SampleMean']*2-ci,label=case + '(%s)'% (statvarlabel),alpha = 0.5)
							else:
								ax.fill_between(avgData[var+'StatTime'],avgData[var+'SampleMean']+ci,
											avgData[var+'SampleMean']-ci,label=case + '(%s)'% (statvarlabel),alpha = 0.5)
				except Exception as error:
					print('WARNING! Unable to plot statistics data!')
					print(error)
						

		caseMeans = np.array(caseMeans)
		if args.yscaling == 'default':
			ax.set_ylim(min(caseMeans)-(np.mean(caseMeans)*0.1),max(caseMeans)+(np.mean(caseMeans)*0.1))
		elif args.yscaling == 'matDefault':
			print('')
		# elif args.yscaling == 'bounds':
		# 	ax.set_ylim(min(caseMins)-(min(caseMins)*0.1),max(caseMax)+(max(caseMax)*0.5))
		
		
		ax.set_ylabel(labelDict[var]['label'])
		ax.set_title(labelDict[var]['title'])
		ax.set_xlabel('Time (s)')
		plt.legend()
		if caseLoc.lower() == 'outtrial':
			plt.savefig('%s/%s_forceHistory_%s.%s' % (list(casePathDict.keys())[0],'_'.join(args.trial),var,args.saveFormat),dpi = 300,bbox_inches='tight')
		else :
			plt.savefig('%s_forceHistory_%s.%s' % ('_'.join(args.trial),var,args.saveFormat),dpi = 300,bbox_inches='tight')




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
























