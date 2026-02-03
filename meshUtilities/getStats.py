import os 
import sys
import glob
import numpy as np
import csv
import pandas as pd


def main():
    trialDirs = glob.glob('*')
    summaryList = pd.DataFrame()
    for trialDir in trialDirs:
        if os.path.isdir(trialDir):
            if os.path.exists(os.path.join(trialDir,'summary.csv')):
                summary = pd.read_csv(os.path.join(trialDir,'summary.csv'),header=None,delimiter=',',index_col=0)
                # if len(summaryList) == 0:
                #     summaryList = pd.DataFrame(index=summary.iloc[:,0])
                if summary.loc['Mesher'].isnull().values.any():
                    summary.loc['Mesher'] = 'ansaMesh'  
                summaryList = pd.concat([summaryList,summary.iloc[:,0]],axis=1)
    summaryList.transpose().to_csv('summaryData.csv',index=False,header=True)   

main()
    