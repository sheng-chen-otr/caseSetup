import os
import sys
import glob
import gzip
import argparse


parser = argparse.ArgumentParser(
                    prog='Archive Compress',
                    description='Reduces the size of a job in the archive',
                    epilog='Helps reduce the size of the archived data, especially data that is older. This must be run in the top level of the job folder.')


parser.add_argument('--compress_case',
                    action='store_true', help='compresses each case')
parser.add_argument('--clean_case',
                    action='store_true',
                    help='proceeds to clean up and remove postProcessing directories and any stls or extendedFeatureEdgeMesh directories')
parser.add_argument('--compress_geom',
                    action='store_true',
                    help='compresses the files in MSH folder')

args = parser.parse_args()

def checkLocation():
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd,'CASES')) and 'archive' in cwd:
        print('Working in %s' % (os.path.basename(cwd)))
    else:
        print('ERROR! Not run in a job directory or CASES directory does not exist!')

def removeCaseInternals():
    caseList = glob.glob(os.path.join('CASES','*'))
    for case in caseList:
        for postDir in glob.glob(os.path.join(case,'*')):
            removeList = ['postProcessing','fidelityMesh_']
            dirName = os.path.basename(postDir)
            if any(x in dirName for x in removeList):
                print('Removing ' + postDir)
                #os.system('rm -r %s' % (postDir))
        for trisurface in glob.glob(os.path.join(case,'constant','triSurface','*')):
            print('Removing ' + trisurface)
            os.system('rm -r %s' % (trisurface))
        for edgemesh in glob.glob(os.path.join(case,'constant','extendedFeatureEdgeMesh')):
            print('Removing ' + edgemesh)
            os.system('rm -r %s' % (edgemesh))
        for polyMesh in glob.glob(os.path.join(case,'constant','polyMesh')):
            print('Removing ' + polyMesh)
            os.system('rm -r %s' % (polyMesh))
        if args.compress_case:
            if 'tar.gz' in os.path.basename(case):
                continue
            else:
                print('compressing %s' % (os.path.basename(case)))
                os.system('tar czf %s.tar.gz %s' % (case,case))
                print('deleting %s' % (os.path.basename(case)))
                os.system('rm -r %s' % (case))


def compressGeometry():
    stlList = glob.glob(os.path.join('02_reference','MSH','*'))
    for stl in stlList:
        stlName = os.path.basename(stl)
        extList = ['.stl','.obj']
        if any(x in stlName for x in extList):
            if stlName.endswith('.gz'):
                print('skipping' + stlName)
            else:
                print('compressing ' + stlName)
                compression = compressFile(stl)
                if compression == 0:
                    print('compression successful, removing %s' % (stlName))
                    os.system('rm %s' % (stl))
                    


def compressFile(filePath):
    import shutil
    import time
    try:
        time_start = time.time()
        with open(filePath, 'rb') as f_in:
            with gzip.open(filePath + '.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                f_out.close()
            f_in.close()
        time_end = time.time() - time_start
        print('\t\tFile compressed in ',time_end)
        return 0
    except Exception as e:
        print('ERROR! File compression was not successful!')
        print('e')
        sys.exit()
        return 1

def main():
    checkLocation()
    if args.clean_case:
        removeCaseInternals()
    
    if args.compress_geom:
        compressGeometry()

main()