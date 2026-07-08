import os
import pathlib as Path
import shutil
import configparser
import argparse



#script to archive openfoam cases


# Parse command line arguments
parser = argparse.ArgumentParser(
    prog='OpenFOAM Archival Tool',
    description='Archive OpenFOAM case data'
)
parser.add_argument('-t', '--trial', nargs='+',
                   help='Trial range to archive')
parser.add_argument('-c', '--copy', action='store_true',
                   help='Copy data to archive')
parser.add_argument('-n', '--noCompress', action='store_true',
                   help='Do not compress data when archiving, compresses by default')
parser.add_argument('-r', '--retrieve', action='store_true',
                   help='Retrieve data from archive')
parser.add_argument('--compress', action='store_true',
                   help='Compress the selected cases already in the archive location. Uncompressed cases are compressed to .tar.gz; already-compressed cases are skipped.')
parser.add_argument('--overwrite', action='store_true',
                   help='Merge into existing archived/retrieved data: refresh files in place and add missing ones, but never delete existing files. Compressed .tar.gz archives are left untouched.')   
parser.add_argument('--keepEnsight', action='store_true',
                   help='Keep Ensight files when archiving')
parser.add_argument('--archiveRefData', action='store_true',
                   help='Also copy the job reference data (02_reference, e.g. MSH) to the archive. Reference data is copied, never removed.')

parser.add_argument('--batchMode', action='store_true',
                   help='Run in batch mode without user prompts.')


args = parser.parse_args()

defaultArchivePath = '/media/openfoam/archive/jobs'

def main():
    job,parentPath = getJob()

    if args.trial is None:
        # reference data is job-level (not tied to any trial), so it can be archived on
        # its own without a trial range - just run the standalone reference archive
        if args.archiveRefData:
            summarizeReferenceData(job, parentPath)
            if not args.batchMode:
                if not confirmArchive('Proceeding will copy the reference data above to '
                                    'the archive location. The local copy is kept.'):
                    print("Archive cancelled by user. No changes made.")
                    return
                archiveReferenceData(job, parentPath)
            return
        print("Please specify a trial range using the -t or --trial option.")
        return
    if len(args.trial) == 1:
        trialRange = [int(args.trial[0])]
    elif len(args.trial) == 2:
        trialRange = range(int(args.trial[0]), int(args.trial[1]) + 1)
    else:
        print("Please specify either a single trial or a range of trials (two values).")
        return

    # retrieve pulls cases back out of the archive into the local job; it's the inverse
    # of the -c archive move and doesn't touch anything already present unless --overwrite
    if args.retrieve:
        retrieveCases(job, trialRange)
        return

    # compress only operates on data already in the archive: uncompressed archived cases
    # are tarred up, cases already compressed are left alone
    if args.compress:
        compressArchivedCases(job, trialRange)
        return

    if len(args.trial) == 2:
        print(f"Archiving trials {args.trial[0]}-{args.trial[1]}")

    # collect every case (trial or ride-height point) that will be archived
    casesToProcess, rideHeightParents = collectCases(trialRange)
    if len(casesToProcess) == 0:
        print("No valid cases found to archive.")
        return

    # only copying to the archive location (-c) is gated behind a status summary and
    # confirmation, since that's the operation the user reviews before committing to
    if args.copy:
        summarizeCases(job, casesToProcess)
        if args.archiveRefData:
            summarizeReferenceData(job, parentPath)
        if not args.batchMode:
            if not confirmArchive():
                print("Archive cancelled by user. No changes made.")
                return

    for casePath in casesToProcess:
        print('\tArchiving %s...' % (casePath))
        if args.copy:
            source = casePath
            dest = os.path.join(defaultArchivePath, job, 'CASES', casePath)
            copied = copyData(job, parentPath, source, dest)
            #-c is a move: once the case is safely in the archive, remove the original
            #entirely. if the copy was skipped (already archived, no --overwrite) leave
            #the original untouched to avoid deleting data that wasn't refreshed
            if copied:
                print('\t\tRemoving original case %s...' % (casePath))
                shutil.rmtree(casePath)
        else:
            clearCase(casePath)

    # ride-height trials keep shared data (e.g. caseSetup) directly in the trial dir
    # alongside the point cases; archive that too and drop the emptied trial dir
    if args.copy:
        for trialName in rideHeightParents:
            archiveTrialParent(job, trialName)

    # reference data (02_reference) is shared job input, so it is only copied to the
    # archive - never removed locally like the cases are
    if args.copy and args.archiveRefData:
        archiveReferenceData(job, parentPath)


def archiveTrialParent(job, trialName):
    '''
    Archives the shared, parent-level data that lives directly in a ride-height trial
    directory (e.g. caseSetup) into the archived trial directory, then removes the
    local trial directory once it no longer holds any point cases.
    '''
    if not os.path.exists(trialName):
        return
    archiveTrialDir = os.path.join(defaultArchivePath, job, 'CASES', trialName)
    os.makedirs(archiveTrialDir, exist_ok=True)
    #points still present locally were skipped (already archived, no --overwrite); leave them
    remainingPoints = collectRideHeightPoints(trialName)
    for item in os.listdir(trialName):
        if item in remainingPoints:
            continue
        src = os.path.join(trialName, item)
        dst = os.path.join(archiveTrialDir, item)
        if os.path.exists(dst) and not args.overwrite:
            print(f"\t\tParent data {os.path.join(trialName, item)} already archived. Skipping.")
            continue
        print(f"\t\tArchiving parent data {os.path.join(trialName, item)}...")
        if os.path.isdir(src):
            #merge into any existing archived dir: refresh/add files, keep extras
            shutil.copytree(src, dst, dirs_exist_ok=True)
            shutil.rmtree(src)
        else:
            shutil.copy2(src, dst)
            os.remove(src)
    #only remove the trial dir once it's empty (all points and parent data archived)
    if len(os.listdir(trialName)) == 0:
        print(f"\t\tRemoving empty trial directory {trialName}...")
        os.rmdir(trialName)
    else:
        print(f"\t\tTrial directory {trialName} still has cases; leaving in place.")


def retrieveCases(job, trialRange):
    '''
    Copies/uncompresses cases from the archive location back into the local job.
    Cases already present locally are skipped unless --overwrite is given.
    '''
    casesToRetrieve, rideHeightParents = collectArchivedCases(job, trialRange)
    if len(casesToRetrieve) == 0:
        print("No archived cases found to retrieve.")
        return
    for casePath in casesToRetrieve:
        print('\tRetrieving %s...' % (casePath))
        retrieveData(job, casePath)
    #bring back the shared parent-level data for ride-height trials too
    for trialName in rideHeightParents:
        retrieveTrialParent(job, trialName)


def retrieveTrialParent(job, trialName):
    '''
    Restores the shared, parent-level data of a ride-height trial (e.g. caseSetup)
    from the archive back into the local trial directory. Existing local files are
    skipped unless --overwrite is given.
    '''
    archiveTrialDir = os.path.join(defaultArchivePath, job, 'CASES', trialName)
    if not os.path.isdir(archiveTrialDir):
        return
    os.makedirs(trialName, exist_ok=True)
    for item in os.listdir(archiveTrialDir):
        #skip the point cases themselves - those are restored by retrieveData
        if item.startswith(trialName + '_'):
            continue
        src = os.path.join(archiveTrialDir, item)
        dst = os.path.join(trialName, item)
        if os.path.exists(dst) and not args.overwrite:
            print(f"\t\tParent data {dst} already exists locally. Skipping.")
            continue
        print(f"\t\tRetrieving parent data {dst}...")
        if os.path.isdir(src):
            #merge into any existing local dir: refresh/add files, keep extras
            shutil.copytree(src, dst, dirs_exist_ok=True)
        else:
            shutil.copy2(src, dst)


def collectArchivedCases(job, trialRange):
    '''
    Expands a trial range into a flat list of case paths that exist in the archive,
    resolving ride-height trials into their individual ride-height point cases.
    Mirrors collectCases() but looks in the archive location instead of the cwd.
    '''
    casesRoot = os.path.join(defaultArchivePath, job, 'CASES')
    if not os.path.exists(casesRoot):
        print(f"No archive found for job {job} at {casesRoot}.")
        return []
    cases = []
    rideHeightParents = []
    for trial in trialRange:
        trialName = f"{trial:03d}"
        trialDir = os.path.join(casesRoot, trialName)
        trialArchive = trialDir + '.tar.gz'
        if os.path.isdir(trialDir):
            #a ride-height container holds NNN_M cases (dirs or .tar.gz); a plain
            #uncompressed case holds OpenFOAM folders, so has no NNN_ entries
            points = [item for item in os.listdir(trialDir)
                      if item.startswith(trialName + '_')]
            if len(points) > 0:
                for point in sorted(points):
                    pointName = stripArchiveExt(point)
                    cases.append(os.path.join(trialName, pointName))
                rideHeightParents.append(trialName)
            else:
                cases.append(trialName)
        elif os.path.exists(trialArchive):
            cases.append(trialName)
        else:
            print(f"Case {trialName} not found in archive. Skipping.")
    return cases, rideHeightParents


def stripArchiveExt(name):
    return name[:-len('.tar.gz')] if name.endswith('.tar.gz') else name


def compressArchivedCases(job, trialRange):
    '''
    Compresses the selected cases that already live in the archive location. An
    uncompressed archived case (a directory) is compressed to <case>.tar.gz and the
    directory removed; a case that is already compressed is skipped.
    '''
    result = collectArchivedCases(job, trialRange)
    if not result:
        return
    cases, rideHeightParents = result
    casesRoot = os.path.join(defaultArchivePath, job, 'CASES')

    # figure out which cases actually need compressing before touching anything
    toCompress = []
    for casePath in cases:
        archivePath = os.path.join(casesRoot, casePath)
        if os.path.exists(archivePath + '.tar.gz'):
            print(f"\t{casePath} already compressed. Skipping.")
        elif os.path.isdir(archivePath):
            toCompress.append(casePath)
        else:
            print(f"\t{casePath} not found in archive. Skipping.")

    if len(toCompress) == 0:
        print("No uncompressed archived cases to compress.")
        return

    for casePath in toCompress:
        archivePath = os.path.join(casesRoot, casePath)
        print(f"\tCompressing {casePath}...")
        shutil.make_archive(archivePath, 'gztar', archivePath)
        shutil.rmtree(archivePath)
    print("Compression complete.")


def retrieveData(job, casePath):
    '''
    Restores a single case from the archive back into the local job location.
    Returns True if the case was restored, False if skipped or missing.
    '''
    source = os.path.join(defaultArchivePath, job, 'CASES', casePath)
    dest = casePath  # relative to the CASES cwd, matches the archive layout

    # skip cases already present locally unless the user forces an overwrite
    if os.path.exists(dest) and not args.overwrite:
        print(f"Case {dest} already exists locally. Skipping retrieval.")
        return False
    if os.path.exists(dest) and args.overwrite:
        print(f"Updating existing local case {dest} (existing files kept, missing ones added).")

    # ride-height cases live under a trial subdir that must exist first
    parent = os.path.dirname(dest)
    if parent:
        os.makedirs(parent, exist_ok=True)

    if os.path.isdir(source):
        print(f"\t\tCopying uncompressed case from archive...")
        #merge into any existing local case: refresh/add files, keep extras
        shutil.copytree(source, dest, dirs_exist_ok=True)
    elif os.path.exists(source + '.tar.gz'):
        print(f"\t\tUncompressing case from archive...")
        os.makedirs(dest, exist_ok=True)
        shutil.unpack_archive(source + '.tar.gz', dest)
    else:
        print(f"Archived data for {casePath} not found. Skipping.")
        return False
    return True



def collectCases(trialRange):
    '''
    Expands a trial range into a flat list of case paths, resolving ride-height
    trials into their individual ride-height point directories.
    '''
    cases = []
    rideHeightParents = []
    for trial in trialRange:
        trialName = f"{trial:03d}"
        if not os.path.exists(trialName):
            print(f"Case path {trialName} does not exist. Skipping.")
            continue
        rideHeightPoints = collectRideHeightPoints(trialName)
        if len(rideHeightPoints) > 0:
            for point in sorted(rideHeightPoints):
                casePath = os.path.join(trialName, point)
                if not os.path.exists(casePath):
                    print(f"Case path {casePath} does not exist. Skipping.")
                    continue
                cases.append(casePath)
            rideHeightParents.append(trialName)
        else:
            cases.append(trialName)
    return cases, rideHeightParents


#column layout shared by the summary header and rows
_summaryRowFmt = '{:<20}{:>10}{:>7}{:>6}{:>9}{:>9}{:>6}{:>10}'


def summarizeCases(job, casesToProcess):
    '''
    Prints a per-case status summary of everything that is about to be archived,
    so the user can review before confirming.
    '''
    print('\n' + '=' * 77)
    print('  ARCHIVE SUMMARY - job %s' % (job))
    print('=' * 77)
    print(_summaryRowFmt.format('CASE', 'SIZE', 'PROCS', 'MESH',
                                'RESULTS', 'ENSIGHT', 'ANSA', 'ARCHIVED'))
    print('-' * 77)
    for casePath in casesToProcess:
        status = getCaseStatus(job, casePath)
        print(_summaryRowFmt.format(status['case'], status['size'], status['procs'],
                                    status['mesh'], status['results'], status['ensight'],
                                    status['ansa'], status['archived']))
    print('=' * 77)


#name of the job reference-data folder (shared input, e.g. meshes)
referenceDirName = '02_reference'


def summarizeReferenceData(job, parentPath):
    '''
    Prints a summary of the job reference data that will be copied to the archive,
    breaking out the MSH folder contents so the user can review before confirming.
    '''
    refPath = os.path.join(parentPath, referenceDirName)
    print('\n' + '=' * 77)
    print('  REFERENCE DATA SUMMARY - %s (copied, not removed)' % (referenceDirName))
    print('=' * 77)
    if not os.path.exists(refPath):
        print('  No %s directory found for job %s.' % (referenceDirName, job))
        print('=' * 77)
        return
    archived = isReferenceArchived(job)
    print('  Location : %s' % (refPath))
    print('  Size     : %s' % (humanSize(getDirSize(refPath))))
    print('  Archived : %s' % ('yes' if archived else 'no'))

    mshPath = os.path.join(refPath, 'MSH')
    print('-' * 77)
    if os.path.isdir(mshPath):
        mshFiles = sorted(f for f in os.listdir(mshPath)
                          if os.path.isfile(os.path.join(mshPath, f)))
        print('  MSH folder (%s files, %s):' % (len(mshFiles), humanSize(getDirSize(mshPath))))
        for f in mshFiles:
            fp = os.path.join(mshPath, f)
            print('    {:<40}{:>12}'.format(f, humanSize(os.path.getsize(fp))))
    else:
        print('  MSH folder: not found')
    print('=' * 77)


def isReferenceArchived(job):
    archiveRefPath = os.path.join(defaultArchivePath, job, referenceDirName)
    return os.path.exists(archiveRefPath)


def archiveReferenceData(job, parentPath):
    '''
    Copies the job reference data (02_reference) into the archive. Reference data is
    shared job input, so it is only copied - never removed from the local job.
    '''
    refPath = os.path.join(parentPath, referenceDirName)
    if not os.path.exists(refPath):
        print('\tNo %s directory to archive for job %s.' % (referenceDirName, job))
        return
    dest = os.path.join(defaultArchivePath, job, referenceDirName)
    print('\tArchiving reference data (%s)...' % (referenceDirName))
    if os.path.exists(dest) and not args.overwrite:
        print('\t\tReference data already exists in archive. Skipping copy.')
        return
    if os.path.exists(dest) and args.overwrite:
        print('\t\tUpdating existing reference data in archive (existing files kept, missing ones added).')
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    #merge into any existing archived reference data: refresh/add files, keep extras
    shutil.copytree(refPath, dest, dirs_exist_ok=True)
    print('\t\tCopied %s to archive.' % (referenceDirName))


def getCaseStatus(job, casePath):
    '''
    Gathers the status of a single case: on-disk size and which artefacts are
    present (processor dirs, mesh, results, EnSight, ANSA) plus whether it is
    already present in the archive location.
    '''
    contents = os.listdir(casePath)
    processors = [d for d in contents if d.startswith('processor')]
    hasMesh = os.path.exists(os.path.join(casePath, 'constant/polyMesh'))
    hasEnsight = any(f.endswith('EnSight') for f in contents)
    hasAnsa = any(f.endswith('.ansa.gz') for f in contents)
    hasResults = os.path.exists(os.path.join(casePath, 'postProcessing')) or \
                 any(isTimeDir(os.path.join(casePath, d)) for d in contents)
    return {
        'case': casePath,
        'size': humanSize(getDirSize(casePath)),
        'procs': len(processors),
        'mesh': 'yes' if hasMesh else 'no',
        'results': 'yes' if hasResults else 'no',
        'ensight': 'yes' if hasEnsight else 'no',
        'ansa': 'yes' if hasAnsa else 'no',
        'archived': 'yes' if isArchived(job, casePath) else 'no',
    }


def isTimeDir(path):
    #an OpenFOAM results time directory is a directory named with a positive number
    if not os.path.isdir(path):
        return False
    try:
        return float(os.path.basename(path)) > 0
    except ValueError:
        return False


def isArchived(job, casePath):
    #matches copyData's archive layout: <archive>/<job>/CASES/<case>[.tar.gz]
    archiveCasePath = os.path.join(defaultArchivePath, job, 'CASES', casePath)
    return os.path.exists(archiveCasePath) or os.path.exists(archiveCasePath + '.tar.gz')


def getDirSize(path):
    total = 0
    for root, dirs, files in os.walk(path):
        for f in files:
            fp = os.path.join(root, f)
            if not os.path.islink(fp):
                try:
                    total += os.path.getsize(fp)
                except OSError:
                    pass
    return total


def humanSize(nbytes):
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if nbytes < 1024 or unit == 'TB':
            return f"{nbytes:.1f}{unit}"
        nbytes /= 1024


def confirmArchive(message='Proceeding will copy the cases above to the archive '
                   'location and then remove the originals.'):
    print('\n' + message)
    while True:
        response = input('Proceed with archiving? [y/N]: ').strip().lower()
        if response in ('y', 'yes'):
            return True
        if response in ('', 'n', 'no'):
            return False
        print("Please answer 'y' or 'n'.")


def collectRideHeightPoints(trialName):
    # Check if the trial directory exists
    if not os.path.exists(trialName):
        print(f"ERROR! Trial directory {trialName} does not exist.")
        return []

    # List all items in the trial directory
    caseContents = os.listdir(trialName)

    # Filter for directories that match the ride height point pattern
    rideHeightPoints = [item for item in caseContents if item.startswith(trialName + "_") and os.path.isdir(os.path.join(trialName, item))]

    return rideHeightPoints

        
def copyData(job,jobPath,source,dest):
    '''
    Copies data from source to archive location. If case exists in archive, and no overwrite flag given, the case will skip copying.
    If it's the first time the job is archived, then the job will be created. If the archive reference data flag is given, then the reference data will be copied to the archive location.
    '''

    # check if job exists in archive location
    archiveJobPath = os.path.join(defaultArchivePath, job)
    if not os.path.exists(archiveJobPath):
        print(f"Job {job} does not exist in archive. Creating job directory.")
        os.makedirs(archiveJobPath)
        # make the CASES directory
        casesPath = os.path.join(archiveJobPath, "CASES")
        os.makedirs(casesPath)
    else:
        # check if CASES directory exists, if not create it
        casesPath = os.path.join(archiveJobPath, "CASES")
        if not os.path.exists(casesPath):
            print(f"Creating CASES directory for job {job}.")
            os.makedirs(casesPath)
    
    # ensure the destination parent exists (ride-height cases live under trial subdirs)
    os.makedirs(os.path.dirname(dest), exist_ok=True)

    if args.noCompress:
        # uncompressed archive: the case lives as a directory
        if os.path.exists(dest):
            if not args.overwrite:
                print(f"Data for {dest} already exists in archive. Skipping copy.")
                return False
            #overwrite merges: refresh files in place, add missing ones, keep extras
            print(f"Updating existing archive data for {dest} (existing files kept, missing ones added).")
        #a compressed archive of the same case would be stale; leave it, but warn
        if os.path.exists(dest + '.tar.gz'):
            print(f"Note: a compressed archive {dest}.tar.gz also exists and is left untouched.")
        shutil.copytree(source, dest, dirs_exist_ok=True)
    else:
        # compressed archive: a tarball can't be merged, so an existing one is left as-is
        if os.path.exists(dest + '.tar.gz'):
            print(f"Compressed archive for {dest} already exists. Skipping (cannot merge into a tarball).")
            return False
        if os.path.exists(dest):
            print(f"Uncompressed archive data for {dest} already exists. Skipping to avoid mixing formats.")
            return False
        shutil.make_archive(dest, 'gztar', source)
    return True



def getJob():
    '''
    Returns the job name based on the current working directory.
    '''
    #get parent directory name of current working directory

    cwd = os.getcwd()
    cwdName = os.path.basename(cwd)
    # check if case run in the right directory
    if cwdName != "CASES":
        return print('ERROR! Please run in CASES directory!')
    parent = os.path.dirname(cwd)

    job = os.path.basename(parent)
    return job,parent



def getFullCaseSetupDict(casePath):
    caseSetupPath="%s/fullCaseSetupDict" % (casePath)
    fullCaseSetupDict = configparser.ConfigParser()
    fullCaseSetupDict.optionxform = str
    fullCaseSetupDict.read_file(open(caseSetupPath))
    configSections = fullCaseSetupDict.sections()
    return fullCaseSetupDict

def clearCase(casePath):
    # clear the processor directories and .ansa files and ensight files
    print('\t\tRemoving processor directories...')
    processorDirs = [d for d in os.listdir(casePath) if d.startswith('processor')]
    for processorDir in processorDirs:
        shutil.rmtree(os.path.join(casePath, processorDir))
    polyMeshPath = os.path.join(casePath,'constant/polyMesh')
    if os.path.exists(polyMeshPath):
        shutil.rmtree(polyMeshPath)
        
    print('\t\tRemoving ANSA files...')
    ansaFiles = [f for f in os.listdir(casePath) if f.endswith('.ansa.gz')]
    for ansaFile in ansaFiles:
        os.remove(os.path.join(casePath, ansaFile))
    if not args.keepEnsight:
        print('\t\tRemoving Ensight files...')
        ensightFiles = [f for f in os.listdir(casePath) if f.endswith('EnSight')]
        for ensightFile in ensightFiles:
            shutil.rmtree(os.path.join(casePath, ensightFile))




main()