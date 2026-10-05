"""
Copyright (c) 2026 Zhi Sheng Chen

Licensed under the PolyForm Noncommercial License 1.0.0 (the "License").
You may not use this file except in compliance with the License. Any
commercial use requires a separate commercial license from the copyright
holder.

You may obtain a copy of the License at:
    https://polyformproject.org/licenses/noncommercial/1.0.0

Or see the LICENSE file in the root of this repository.

This software is provided "as is", without warranty of any kind, express or
implied. See the License for the specific language governing permissions
and limitations.
"""

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
    