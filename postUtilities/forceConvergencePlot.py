# %%
import os
import sys
import numpy as np
import re
import pandas as pd
import configparser
import argparse
import matplotlib.pyplot as plt
import scipy
import scipy.stats as st
import glob
from estimateStatisticalError import *



plt.style.use(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'pltStyles.txt')) #latex styling

# #### Paths ####
# cdHistoryPath = '/Users/zchen147/Library/CloudStorage/OneDrive-azureford/Documents/researchMasc/data/sheng_wheel_data/loaded/dataFiles/r0_undef/cd-history.csv'
# clHistoryPath = '/Users/zchen147/Library/CloudStorage/OneDrive-azureford/Documents/researchMasc/data/sheng_wheel_data/loaded/dataFiles/r0_undef/cl-history.csv'
# savePath = '/Users/zchen147/Library/CloudStorage/OneDrive-azureford/Documents/researchMasc/data/sheng_wheel_data/loaded/dataFiles/r0_undef'
# #### End Paths ####

# #### Settings ####
# avgStartTime = 0.324

# #### End Settings ####

# def main():
#     cdData = calculateStatistics(cdHistoryPath,forceName='cd')
#     clData = calculateStatistics(clHistoryPath,forceName='cl')
#     plotData(cdData,'$C_{D}$',ylim = [0.9, 1.15],text='Avg. Start', textpos=[avgStartTime+0.01, 1.1], fontsize=14)
#     plotData(clData,'$C_{L}$',ylim = [0.65, 0.85],text='Avg. Start', textpos=[avgStartTime+0.01, 0.8], fontsize=14)
    
    
#     plt.show()


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

        avgData = {}
        rawTime = forceData['#Time'].to_numpy()
        for var in args.plotData:
            if args.skipStats == False:
                data = avgForceData[var].to_numpy()		
                print('\t\t\tCalculating statistics for %s' % (var))
                dt = time[-1]-time[-2]
                if len(time) > 10:
                    avgData = calculateStatistics(avgData,time,data,var)
   
                else:
                    print('In-sufficient time to calculate statistics!')
            
            casePathDict[case]['avgData'] = avgData

    return casePathDict

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

def plotData(args,caseLoc,casePathDict):
    labelDict = {'Cl':{'label':'$C_L$','title':'$C_L$'},
                    'Cd':{'label':'$C_D$','title':'$C_D$'},
                    'Cl(f)':{'label':'$C_{L,f}$','title':'$C_{L,f}$'},
                    'Cl(r)':{'label':'$C_{L,r}$','title':'$C_{L,r}$'},
                    'Cs(f)':{'label':'$C_{S,f}$','title':'$C_{S,f}$'},
                    'Cs(r)':{'label':'$C_{S,r}$','title':'$C_{S,r}$'},
                    'CoP':{'label':'Percent Front','title':'Center of Pressure'}
                }

    for var in args.plotData:

        fig1,ax = plt.subplots(figsize=[10,6],frameon=True)

        

        caseEndTimes = []
        caseMeans = []
        caseMins = []
        caseMax = []
        for case in casePathDict.keys():
            rawData = casePathDict[case]['data']
            avgData = casePathDict[case]['avgData']
            # for key in avgData.keys():
            #     print(key)
            data = rawData[['#Time',var]]
            avgStart = casePathDict[case]['avgStart']
            caseEndTimes.append(max(data['#Time']))
            if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
                if var == 'CoP':
                    factor = 1
                else:
                    factor = 2
                caseMeans.append(np.mean(data[var][data['#Time']>=avgStart])*factor) #using this to calculate the range for scaling axis
                
            else:
                caseMeans.append(np.mean(data[var][data['#Time']>=avgStart])) #using this to calculate the range for scaling axis\
                
            if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
                if var == 'CoP':
                    factor = 1
                else:
                    factor = 2
                ax.plot(rawData['#Time'],rawData[var]*factor,label=case,alpha=0.5)
            else:
                ax.plot(rawData['#Time'],rawData[var],label=case,alpha=0.5)
            
            if args.skipStats == False:
                for statvar in ['SampleMean', 'CI']:
                    if statvar == 'SampleMean':
                        print('\tPlotting forward averages for %s' % (var))
                        statvarlabel = 'Fwd. Avg'
                        if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
                            if var == 'CoP':
                                factor = 1
                            else:
                                factor = 2
                                
                            ax.plot(avgData[var + 'StatTime'], avgData[var + statvar] * factor, label=case + '(%s)' % (statvarlabel), linestyle='--')
                        else:
                            ax.plot(avgData[var + 'StatTime'], avgData[var + statvar], label=case + '(%s)' % (statvarlabel), linestyle='--')
                    elif statvar == 'CI':
                        statvarlabel = statvar
                        #ci = (avgData[var + 'UpperCI'] - avgData[var + 'LowerCI']) / 2
                        if 'half' in caseSetupConfig['GLOBAL_SIM_CONTROL']['SIM_SYM'].lower():
                            if var == 'CoP':
                                factor = 1
                            else:
                                factor = 2
                            try:
                                ax.errorbar(avgData[var + 'ciTime'],
                                avgData[var + 'ciMean']*factor,
                                yerr=avgData[var + 'ciRaw']*factor,
                                fmt='none',capsize=3,linewidth=1)
                            except Exception as error: 
                                print('\t\tWARNING! Unable to plot fwd avg data...')
                                print('\t\t\t%s' % (error))
                                
                        else:
                            try:
                                ax.errorbar(avgData[var + 'ciTime'],
                                avgData[var + 'ciMean'],
                                yerr=avgData[var + 'ciRaw'],
                                fmt='none',capsize=3,linewidth=1)
                                
                            except Exception as error: 
                                print('\t\tWARNING! Unable to plot fwd avg data...')
                                print('\t\t\t%s' % (error))
                                
                        
                  
        caseMeans = np.array(caseMeans)
        
        if args.yscaling == 'default':
            ax.set_ylim(min(caseMeans)-(np.mean(caseMeans)*0.1),max(caseMeans)+(np.mean(caseMeans)*0.1))
        elif args.yscaling == 'matDefault':
            print('')
        else:
            try:
                yscaling = float(args.yscaling)
                ax.set_ylim(min(caseMeans)-(np.mean(caseMeans)*yscaling),max(caseMeans)+(np.mean(caseMeans)*yscaling))
            except Exception as error:
                print(error)
                print('WARNING! Y-axis scaling factor incorrect, maybe not a float or int? using default scaling!')
                ax.set_ylim(min(caseMeans)-(np.mean(caseMeans)*0.1),max(caseMeans)+(np.mean(caseMeans)*0.1))
        
        ax.set_xlim(0,max(caseEndTimes))
        ax.set_ylabel(labelDict[var]['label'])
        ax.set_title(labelDict[var]['title'])
        ax.set_xlabel('Time (s)')

        if len(casePathDict.keys()) < 2:
            fig1.legend(loc='upper center',fontsize=10,
                        
                        )
        else: 
            fig1.legend(loc='upper center',fontsize=10,
                        ncols = len(casePathDict.keys()),
                        ) #kinda wierd that it doesn't like to do 1 col ... so if less than 2 cases it will not have a ncol setting

        if caseLoc.lower() == 'outtrial':
            plt.savefig('%s/%s_forceHistory_%s.%s' % (list(casePathDict.keys())[0],'_'.join(args.trial),var,args.saveFormat),dpi = 300,
                        bbox_inches='tight'
                        )
        else:
            plt.savefig('%s_forceHistory_%s.%s' % ('_'.join(args.trial),var,args.saveFormat),dpi = 300,bbox_inches='tight')

def calculateStatistics(avgData,time,data,forceName,contCiResolution = 100):
        
    avgT, avgVal, std = forward_average(time,data)
    
    ciTime = []
    ciArrayUpper = []
    ciArrayLower = []
    ciRawArray = []
    ciMeanArray = []
    
    blockSpace = np.linspace(0,len(avgT),int(len(avgT)/contCiResolution)).astype(int)
    
    for i in blockSpace:
        if i == 0:
            continue
        ciTime.append(time[i-1])
        mean = data[:i].mean()
        ciMeanArray.append(mean)
        ci = calc_confidence_interval(data[:i-1],i_samp = calc_indept_samples(data[:i-1]))
        
        # estCi = estimate_statistical_error(data[:i-1],dt=time[1]-time[0])['mean_95_confidence_interval']
        # ci = abs(estCi[1]-estCi[0])/2
        ciArrayUpper.append(mean+ci)
        ciArrayLower.append(mean-ci)
        ciRawArray.append(ci)

    
    avgData[forceName +'rawData'] = data
    avgData[forceName + 'StatTime'] = np.array(avgT)
    avgData[forceName + 'SampleMean'] = np.array(avgVal)
    avgData[forceName +'ciTime'] = np.array(ciTime)
    avgData[forceName + 'UpperCI'] = np.array(ciArrayUpper)
    avgData[forceName + 'LowerCI'] = np.array(ciArrayLower)
    avgData[forceName +'ciRaw'] = np.array(ciRawArray)
    avgData[forceName +'ciMean'] = np.array(ciMeanArray)

    print('\t\t\tFinal CI: ' + str(ci))
    
    outputData = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in avgData.items() ]))
    outputData.to_csv(path_or_buf=os.path.join(forceName + '_statistics.csv'),na_rep='')
    
    return avgData
    

def get_moving_average(force, ma_window):
    chunk = int(0.5 * ma_window)
    force_ma = np.zeros_like(force)
    for i in range(len(force) - 1, -1, -1):
        f_window = min(len(force - 1) - i, chunk)
        r_window = min(i, chunk)
        force_ma[i] = np.average(force[i - r_window:i + f_window])
    return force_ma


def forward_average(time, arr, start_time=None):
    """Compute an array showing evolution of the mean as more samples are
    included.

    Returns
    -------
    time : numpy.ndarray
        Time array over which averaging was performed.
    ave_arr : numpy.ndarray
        Fordward-averaged array.
    """
    time_t = np.array(time)
    arr = np.array(arr)
    if start_time is not None:
        arr = arr[time >= start_time]
        time_t = time[time >= start_time]
    t_samp = 1
    try:
        t_samp = calc_indept_samples(arr)
    except:
        print("error in calculating independent samples, assuming n_samp = 1")
    i_sample = np.linspace(0, t_samp, len(arr))
    interval = 1.0 - (1.0 - 0.95) / 2.0
    std_arr = scipy.stats.t.ppf(interval, t_samp) * np.array(pd.Series(arr).expanding(2).std(ddof=0)) / np.sqrt(
        i_sample)
    ave_arr = np.array(pd.Series(arr).expanding(0).mean())
    return time_t, ave_arr, std_arr

def smooth(data):
    """Smooth data with a 3rd order Butterworth filter."""
    b, a = scipy.signal.butter(3, 7 / len(data))
    fdata = scipy.signal.filtfilt(b, a, data)
    return fdata

def two_tailed_t_score(confidence_interval, n):
    """
    Calculate the two-tailed t-score with the given confidence
    interval (0 - 1.0) and number of samples
    """

    if confidence_interval is not None:
        interval = 1.0 - (1.0 - confidence_interval) / 2.0
        return scipy.stats.t.ppf(interval, n - 1)
    else:
        return 1.0

def get_rr_var(x, n_bins, batch_overlap=0.75):
    """
    #Get the estimated variance of using repeat reference method

    Keyword arguments:
    x -- array to generate autocorrelation function of
    n_bins -- number of bins to divide the signal into
    batch_overlap -- amount to overlap the batches by
    """

    batch_len = int(len(x) / n_bins)
    batch_overlap = int(0.75 * batch_len)
    cumsum = np.cumsum(np.insert(x, 0, 0))
    mavg = (cumsum[batch_len:] - cumsum[:-batch_len]) / float(batch_len)
    step = batch_len - batch_overlap
    return np.var(mavg[::step])


def get_mean_var(x, n_divs_max=40):
    """
    Use method of Mockett to get the estimated variance in the mean

    Keyword arguments:
    x -- array to get variance of
    n_divs_max -- maximum number of divisions to use to find B_min
    """
    x = x - np.mean(x)
    n_divs_max = min(n_divs_max, len(x))
    rr_var = np.zeros(n_divs_max)
    var = np.var(x)
    for i in range(4, n_divs_max):
        t = float(len(x)) / i
        rr_var[i] = var / (2 * t * get_rr_var(x, i))
    try:
        b = np.min(rr_var[4:])
    except:
        return -1
    return var/(2*b*len(x))


def calc_indept_samples(values):
    """
    Compute number of independent samples
    """
    var = get_mean_var(np.copy(values))
    
    return np.var(values) / var


def calc_mean_std_err(values):
    """Calculate the standard error of the mean of a time series."""
    values = np.array(values)
    if len(values) > 1:
        n = calc_indept_samples(values)
        
        err = values.std() / np.sqrt(n)
    else:
        err = np.nan
    return err


def calc_confidence_interval(values, interval=0.95, i_samp=1):
    """Calculate the confidence interval of the mean of a time series."""
    values = np.array(values)
    if len(values) > 1:
        err = (values.std() / np.sqrt(i_samp)) * two_tailed_t_score(interval, i_samp)
    else:
        err = np.nan
    return err



# %%
