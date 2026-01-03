import os
import numpy as np
import fileinput
import sys
import datetime
import time
import re
import math
import openpyxl 
import argparse as argparse
import glob as glob
import pandas as pd
import scipy.stats as st
import configparser
from pathlib import Path

#### WARNING ####
# not currently working

path = os.path.split(os.getcwd())[0]
case = os.path.split(os.getcwd())[1]
job = os.path.basename(os.path.dirname(path))

if not os.getcwd().split('/')[-2].lower() == 'cases':
    sys.exit('ERROR! Not executed in a trial directory!')
print('Reading caseSetup file...')
caseSetupPath="%s/%s/caseSetup" % (path,case)
fullCaseSetupDict = configparser.ConfigParser()
fullCaseSetupDict.optionxform = str
fullCaseSetupDict.read_file(open(caseSetupPath))
configSections = fullCaseSetupDict.sections()

timeOut = 20


def main():
    getTrialsList(jobPath=Path(os.getcwd()).parents[1],job=job)



def getTrialsList(jobPath,job):
    
    trialsListPath = str(jobPath)+'/*' + job +'_Trials_List.xlsx'
    xlFiles = glob.glob(str(jobPath)+'/*' + job +'_Trials_List.xlsx')
    startCheck = datetime.datetime.now()
    while any('~$' in x for x in xlFiles):
        xlFiles = glob.glob(str(jobPath)+'/*' + job +'_Trials_List.xlsx')
        print('\t Waiting... Excel File is currently open...') 

        time.sleep(5)
        if datetime.datetime.now() - startCheck > datetime.timedelta(seconds=timeOut):
            print('Timed out at: %s' % (timeOut))
            print('Unable to open excel file due to file being open, update trials list manually or re-run program at a later time!')
            sys.exit()

    
    

    print('\tWriting data to trials list...')
        
    


main()