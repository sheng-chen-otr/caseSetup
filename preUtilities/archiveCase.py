import os
import argparse
import sys
import shutil
import glob


parser = argparse.ArgumentParser(
                    prog='Archive Case V1.0',
                    description='Archives cases to archive drive',
                    epilog='Can provide compression of cases if requested or run in batch mode on queue')

parser.add_argument('-t','--trial',required=True,nargs='+',help='selects the start and end values, will archive values in that range, includes any cases with "_half", do not include suffixes!')
parser.add_argument('--batch',action='store_true',help='submits the cases into queue and runs it in batch mode')
parser.add_argument('--exclude',required=False,nargs='+',help='case numbers to exclude, do not include _half or _full or any suffix')




args = parser.parse_args()
print("###### Archive Case ######")
path = os.path.split(os.getcwd())[0]
case = os.path.split(os.getcwd())[1]
casesDir = os.getcwd()
def main():
    verifyLocation(case)
    getCases(args.trial)
def getCases(trials):
    try:
        trialStart = int(trials[0])
        trialEnd = int(trials[1])
    except:
        sys.exit('ERROR! Input trials must be integers and not strings or float!')
      
        
    testRange = range(trialStart,trialEnd)
    trialList = []
    for file in testRange:
        fileList = glob.glob('%03d*' % (file))
        if len(fileList) == 0:
            print('\tNo cases with name %03d' % (file))
        elif len(fileList) > 1:
            sys.exit('ERROR! More than 1 case with the same case number!')
        else:
            print('\t%03d found...' % (file))
            trialList.append(fileList[0])

            

def verifyLocation(case):
    if not case == 'CASES':
        sys.exit('ERROR! Must be run in the CASES directory!')
    


main()