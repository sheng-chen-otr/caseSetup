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
import shutil
import subprocess
import argparse


#script to archive openfoam cases
#
#archive layout, per case (and per ride-height point):
#   <archive>/<job>/CASES/<case>/           plain copy of everything small (caseSetup,
#                                           system, constant without polyMesh, logs,
#                                           summary.csv, postProcessing, images, ...)
#   <archive>/<job>/CASES/<case>/EnSight.tar.zst   final averaged EnSight fields
#   <archive>/<job>/02_reference/           job reference data, copied once per job
#processor dirs, constant/polyMesh and ANSA files are not archived: EnSight holds the
#recomposed results and the mesh can be rebuilt from the setup + reference geometry.
#cases archived by older versions as <case>.tar.gz are left untouched.


# Parse command line arguments
parser = argparse.ArgumentParser(
    prog='OpenFOAM Archival Tool',
    description='Archive OpenFOAM case data. With no mode flag (-c, -r, --compress) the '
                'selected local cases are cleaned in place instead: processor dirs, '
                'constant/polyMesh and ANSA files are deleted (EnSight is kept unless '
                '--dropEnsight). Cases with no EnSight output are left untouched.'
)
parser.add_argument('-t', '--trial', nargs='+',
                   help='Trial or trial range (two values) to process')
parser.add_argument('-c', '--copy', action='store_true',
                   help='Move cases to the archive: setup, logs and post-processing are copied '
                        'as plain files, EnSight is compressed to EnSight.tar.zst, then the '
                        'local case is removed. Job reference data is copied once per job.')
parser.add_argument('-r', '--retrieve', action='store_true',
                   help='Retrieve the plain part of archived cases (setup, logs, summary)')
parser.add_argument('--fields', action='store_true',
                   help='With -r, also unpack the archived EnSight fields')
parser.add_argument('--compress', action='store_true',
                   help='Convert old uncompressed archived cases to the current layout: compress '
                        'EnSight and delete archived processor dirs, polyMesh and ANSA files. '
                        'Old .tar.gz archives are skipped.')
parser.add_argument('--overwrite', action='store_true',
                   help='Update cases that already exist at the destination: plain files are '
                        'merged (refreshed/added, never deleted), EnSight is replaced as a whole. '
                        'Old .tar.gz archives are never touched.')
parser.add_argument('--dropEnsight', action='store_true',
                   help='When cleaning in place, also delete EnSight (removes all field results)')
parser.add_argument('--archiveRefData', action='store_true',
                   help='Copy the job reference data (02_reference) to the archive. Runs on its '
                        'own when no -t is given; -c always does this automatically.')
parser.add_argument('--batchMode', action='store_true',
                   help='Run in batch mode without user prompts.')
#EnSight is now kept by default when cleaning; flag accepted so older commands still run
parser.add_argument('--keepEnsight', action='store_true', help=argparse.SUPPRESS)


args = parser.parse_args()

defaultArchivePath = '/media/openfoam/archive/jobs'

#name of the job reference-data folder (shared input, e.g. meshes)
referenceDirName = '02_reference'
#compressed EnSight inside each archived case
ensightTarName = 'EnSight.tar.zst'
partialExt = '.partial'
legacyExt = '.tar.gz'


def main():
    job, parentPath = getJob()

    modeCount = sum([args.copy, args.retrieve, args.compress])
    if modeCount > 1:
        sys.exit('ERROR! Use only one of -c, -r and --compress at a time.')
    if args.fields and not args.retrieve:
        sys.exit('ERROR! --fields only applies to -r.')
    if args.dropEnsight and modeCount > 0:
        sys.exit('ERROR! --dropEnsight only applies when cleaning in place (no -c, -r or --compress).')

    if args.trial is None:
        # reference data is job-level (not tied to any trial), so it can be archived on
        # its own without a trial range
        if args.archiveRefData:
            summarizeReferenceData(job, parentPath)
            if not args.batchMode and not confirm('Proceeding will copy the reference data above to '
                                                  'the archive location. The local copy is kept.'):
                print("Archive cancelled by user. No changes made.")
                return
            archiveReferenceData(job, parentPath)
            return
        sys.exit("Please specify a trial range using the -t or --trial option.")

    try:
        trialNums = [int(t) for t in args.trial]
    except ValueError:
        sys.exit('ERROR! Trials must be numbers, e.g. -t 5 or -t 5 12.')
    if len(trialNums) == 1:
        trialRange = trialNums
    elif len(trialNums) == 2:
        trialRange = range(trialNums[0], trialNums[1] + 1)
    else:
        sys.exit("Please specify either a single trial or a range of trials (two values).")

    if args.retrieve:
        if args.fields:
            requireZstd()
        retrieveCases(job, trialRange)
        return

    if args.compress:
        requireZstd()
        compressArchivedCases(job, trialRange)
        return

    if args.copy:
        requireZstd()

    casesToProcess, rideHeightParents = collectCases(trialRange)
    if len(casesToProcess) == 0:
        print("No valid cases found.")
        return

    if args.copy:
        archiveCases(job, parentPath, casesToProcess, rideHeightParents)
    else:
        clearCases(casesToProcess)


#------------------------------------------------------------------------------
# case contents
#------------------------------------------------------------------------------

def classifyPath(parts):
    '''
    Sorts a path inside a case (given as a list of components relative to the case
    root) into 'plain' (archived as-is), 'ensight' (compressed) or 'dropped'
    (not archived: processor dirs, polyMesh, ANSA files).
    '''
    top = parts[0]
    if top.startswith('processor') or top.endswith('.ansa.gz'):
        return 'dropped'
    if top.endswith('EnSight'):
        return 'ensight'
    if top == 'constant' and len(parts) > 1 and parts[1] == 'polyMesh':
        return 'dropped'
    return 'plain'


def plainOnlyIgnore(caseRoot):
    '''copytree ignore function that copies only the 'plain' part of a case.'''
    caseRoot = os.path.abspath(caseRoot)

    def ignore(directory, names):
        rel = os.path.relpath(os.path.abspath(directory), caseRoot)
        parts = [] if rel == '.' else rel.split(os.sep)
        return {n for n in names if classifyPath(parts + [n]) != 'plain'}
    return ignore


def ensightItems(casePath):
    '''Top-level EnSight output (normally a single EnSight dir) of a case.'''
    return sorted(item for item in os.listdir(casePath) if item.endswith('EnSight'))


def hasEnsightDir(casePath):
    return any(os.path.isdir(os.path.join(casePath, item)) for item in ensightItems(casePath))


def getCaseSizes(casePath):
    '''Walks a case once and returns the byte totals of its plain/ensight/dropped parts.'''
    sizes = {'plain': 0, 'ensight': 0, 'dropped': 0}
    for root, dirs, files in os.walk(casePath):
        rel = os.path.relpath(root, casePath)
        parts = [] if rel == '.' else rel.split(os.sep)
        for f in files:
            fp = os.path.join(root, f)
            if os.path.islink(fp):
                continue
            try:
                sizes[classifyPath(parts + [f])] += os.path.getsize(fp)
            except OSError:
                pass
    return sizes


def removePath(path):
    #EnSight etc. are normally dirs, but don't crash if one is a plain file
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    else:
        os.remove(path)


#------------------------------------------------------------------------------
# zstd tarballs
#------------------------------------------------------------------------------

def requireZstd():
    missing = [tool for tool in ('tar', 'zstd') if shutil.which(tool) is None]
    if missing:
        sys.exit('ERROR! %s not found on PATH. zstd is required to compress/uncompress '
                 'EnSight data - install it (e.g. apt install zstd) and retry.' % (', '.join(missing)))


def makeZstTar(srcDir, members, outPath):
    '''
    Tars `members` (relative to srcDir) into outPath with multithreaded zstd. Written to
    a .partial file, verified by listing it back, then renamed into place, so outPath
    is only ever a complete archive. Raises on any failure.
    '''
    partial = outPath + partialExt
    if os.path.exists(partial):
        os.remove(partial)
    try:
        subprocess.run(['tar', '-I', 'zstd -T0', '-cf', partial, '-C', srcDir] + members,
                       check=True)
        #reading the whole archive back checks both the zstd stream and the tar structure
        subprocess.run(['tar', '-I', 'zstd', '-tf', partial],
                       check=True, stdout=subprocess.DEVNULL)
        os.replace(partial, outPath)
    except BaseException:
        if os.path.exists(partial):
            os.remove(partial)
        raise


def unpackEnsight(tarPath, destCase):
    '''
    Unpacks an EnSight.tar.zst into destCase, replacing any existing EnSight output.
    Unpacks to a temp dir first so a failed unpack leaves the local EnSight untouched.
    '''
    tmp = os.path.join(destCase, '.EnSight.unpack')
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)
    try:
        subprocess.run(['tar', '-I', 'zstd', '-xf', tarPath, '-C', tmp], check=True)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    for item in os.listdir(tmp):
        target = os.path.join(destCase, item)
        if os.path.lexists(target):
            removePath(target)
        os.replace(os.path.join(tmp, item), target)
    os.rmdir(tmp)


#------------------------------------------------------------------------------
# archive (-c)
#------------------------------------------------------------------------------

def casesRoot(job):
    return os.path.join(defaultArchivePath, job, 'CASES')


def isLegacyDir(archivePath):
    '''An uncompressed archive from an older version still holds the heavy data.'''
    for item in os.listdir(archivePath):
        if classifyPath([item]) != 'plain' and item != ensightTarName and not item.endswith(partialExt):
            return True
    return os.path.isdir(os.path.join(archivePath, 'constant', 'polyMesh'))


def archiveState(job, casePath, localHasEnsight):
    '''
    State of a case in the archive: 'none', 'legacyTar' (old .tar.gz), 'legacyDir'
    (old uncompressed dir), 'incomplete' (plain part copied but EnSight tarball
    missing, e.g. an interrupted run) or 'complete'.
    '''
    dest = os.path.join(casesRoot(job), casePath)
    if os.path.exists(dest + legacyExt):
        return 'legacyTar'
    if not os.path.isdir(dest):
        return 'none'
    if isLegacyDir(dest):
        return 'legacyDir'
    if localHasEnsight and not os.path.exists(os.path.join(dest, ensightTarName)):
        return 'incomplete'
    return 'complete'


def archiveAction(state):
    if state == 'none':
        return 'new'
    if state == 'incomplete':
        return 'resume'
    if state == 'complete':
        return 'overwrite' if args.overwrite else 'skip: archived'
    if state == 'legacyTar':
        return 'skip: old .tar.gz'
    return 'skip: run --compress'


def archiveCases(job, parentPath, casesToProcess, rideHeightParents):
    actions = {}
    rows = []
    for casePath in casesToProcess:
        state = archiveState(job, casePath, hasEnsightDir(casePath))
        actions[casePath] = archiveAction(state)
        sizes = getCaseSizes(casePath)
        total = sizes['plain'] + sizes['ensight'] + sizes['dropped']
        rows.append((casePath, humanSize(total), humanSize(sizes['plain']),
                     humanSize(sizes['ensight']), countProcessors(casePath), actions[casePath]))
    printSummary('ARCHIVE SUMMARY - job %s' % (job),
                 ('CASE', 'TOTAL', 'PLAIN', 'ENSIGHT', 'PROCS', 'ACTION'), rows,
                 'PLAIN is copied as-is, ENSIGHT is compressed; processor dirs, polyMesh and '
                 'ANSA files are not archived.')
    for casePath in casesToProcess:
        if not actions[casePath].startswith('skip') and not hasEnsightDir(casePath):
            print('  WARNING! %s has no EnSight output - only its plain files will be archived.' % (casePath))
    summarizeReferenceData(job, parentPath)

    if not args.batchMode and not confirm('Proceeding will copy the cases above to the archive '
                                          'location and then remove the local originals.'):
        print("Archive cancelled by user. No changes made.")
        return

    for casePath in casesToProcess:
        action = actions[casePath]
        if action.startswith('skip'):
            print('\tSkipping %s (%s).' % (casePath, action[len('skip: '):]))
            continue
        print('\tArchiving %s...' % (casePath))
        try:
            archiveCase(job, casePath)
        except (subprocess.CalledProcessError, OSError) as e:
            print('\t\tERROR! Archiving %s failed (%s). Local case kept.' % (casePath, e))
            continue
        #-c is a move: the case is safely in the archive, so remove the original
        print('\t\tRemoving original case %s...' % (casePath))
        shutil.rmtree(casePath)

    # ride-height trials keep shared data (e.g. caseSetup) directly in the trial dir
    # alongside the point cases; archive that too and drop the emptied trial dir
    for trialName in rideHeightParents:
        archiveTrialParent(job, trialName)

    # reference data (02_reference) is shared job input, so it is only copied to the
    # archive - never removed locally like the cases are
    archiveReferenceData(job, parentPath)


def archiveCase(job, casePath):
    '''
    Copies the plain part of a case into the archive (merging into an existing archived
    case) and compresses its EnSight output to EnSight.tar.zst. Raises on failure.
    '''
    dest = os.path.join(casesRoot(job), casePath)
    os.makedirs(dest, exist_ok=True)
    print('\t\tCopying setup, logs and post-processing...')
    shutil.copytree(casePath, dest, ignore=plainOnlyIgnore(casePath), dirs_exist_ok=True)

    tarPath = os.path.join(dest, ensightTarName)
    if hasEnsightDir(casePath):
        print('\t\tCompressing EnSight...')
        makeZstTar(casePath, ensightItems(casePath), tarPath)
    elif os.path.exists(tarPath):
        print('\t\tNo local EnSight output; keeping the EnSight already in the archive.')
    else:
        print('\t\tNo EnSight output to archive.')


def archiveTrialParent(job, trialName):
    '''
    Archives the shared, parent-level data that lives directly in a ride-height trial
    directory (e.g. caseSetup) into the archived trial directory. Local parent data is
    only removed once no point cases remain locally, so a partly archived trial keeps
    its setup alongside the points that were skipped.
    '''
    if not os.path.exists(trialName):
        return
    archiveTrialDir = os.path.join(casesRoot(job), trialName)
    os.makedirs(archiveTrialDir, exist_ok=True)
    remainingPoints = collectRideHeightPoints(trialName)
    for item in os.listdir(trialName):
        if item in remainingPoints:
            continue
        src = os.path.join(trialName, item)
        dst = os.path.join(archiveTrialDir, item)
        if os.path.exists(dst) and not args.overwrite:
            print(f"\t\tParent data {src} already archived. Skipping copy.")
        else:
            print(f"\t\tArchiving parent data {src}...")
            if os.path.isdir(src):
                #merge into any existing archived dir: refresh/add files, keep extras
                shutil.copytree(src, dst, dirs_exist_ok=True)
            else:
                shutil.copy2(src, dst)
        if not remainingPoints:
            removePath(src)
    if len(os.listdir(trialName)) == 0:
        print(f"\t\tRemoving empty trial directory {trialName}...")
        os.rmdir(trialName)
    else:
        print(f"\t\tTrial directory {trialName} still has cases; leaving its parent data in place.")


#------------------------------------------------------------------------------
# reference data
#------------------------------------------------------------------------------

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
    if archived and not args.overwrite:
        action = 'skip (already archived)'
    elif archived:
        action = 'overwrite (merge)'
    else:
        action = 'copy'
    print('  Location : %s' % (refPath))
    print('  Size     : %s' % (humanSize(getDirSize(refPath))))
    print('  Archived : %s' % ('yes' if archived else 'no'))
    print('  Action   : %s' % (action))

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
    return os.path.exists(os.path.join(defaultArchivePath, job, referenceDirName))


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
    if os.path.exists(dest) and not args.overwrite:
        print('\tReference data already in archive. Skipping copy.')
        return
    print('\tArchiving reference data (%s)...' % (referenceDirName))
    if os.path.exists(dest):
        print('\t\tUpdating existing reference data in archive (existing files kept, missing ones added).')
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    #merge into any existing archived reference data: refresh/add files, keep extras
    shutil.copytree(refPath, dest, dirs_exist_ok=True)
    print('\t\tCopied %s to archive.' % (referenceDirName))


#------------------------------------------------------------------------------
# retrieve (-r)
#------------------------------------------------------------------------------

def retrieveCases(job, trialRange):
    '''
    Copies cases from the archive back into the local job: the plain part always,
    EnSight with --fields. Old .tar.gz archives are unpacked whole.
    '''
    casesToRetrieve, rideHeightParents = collectArchivedCases(job, trialRange)
    if len(casesToRetrieve) == 0:
        print("No archived cases found to retrieve.")
        return
    for casePath in casesToRetrieve:
        print('\tRetrieving %s...' % (casePath))
        try:
            retrieveData(job, casePath)
        except (subprocess.CalledProcessError, OSError) as e:
            print('\t\tERROR! Retrieving %s failed (%s).' % (casePath, e))
    #bring back the shared parent-level data for ride-height trials too
    for trialName in rideHeightParents:
        retrieveTrialParent(job, trialName)


def retrieveData(job, casePath):
    '''
    Restores a single case from the archive back into the local job location.
    Without --overwrite an existing local case is left alone, except that --fields
    still unpacks EnSight into it when it has none yet.
    '''
    source = os.path.join(casesRoot(job), casePath)
    dest = casePath  # relative to the CASES cwd, matches the archive layout
    exists = os.path.exists(dest)

    # ride-height cases live under a trial subdir that must exist first
    parent = os.path.dirname(dest)
    if parent:
        os.makedirs(parent, exist_ok=True)

    if os.path.exists(source + legacyExt):
        if exists and not args.overwrite:
            print(f"\t\tCase {dest} already exists locally. Skipping.")
            return
        print('\t\tUnpacking old .tar.gz archive (full case)...')
        os.makedirs(dest, exist_ok=True)
        shutil.unpack_archive(source + legacyExt, dest)
        return

    if not os.path.isdir(source):
        print(f"\t\tArchived data for {casePath} not found. Skipping.")
        return

    if exists and not args.overwrite:
        print(f"\t\tCase {dest} already exists locally. Skipping plain files.")
    else:
        if exists:
            print(f"\t\tUpdating existing local case {dest} (existing files kept, missing ones added).")
        print('\t\tCopying setup, logs and post-processing...')
        shutil.copytree(source, dest, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns(ensightTarName, '*' + partialExt))

    if args.fields:
        tarPath = os.path.join(source, ensightTarName)
        if not os.path.exists(tarPath):
            print('\t\tNo archived EnSight for %s.' % (casePath))
        elif hasEnsightDir(dest) and not args.overwrite:
            print('\t\tEnSight already exists locally. Skipping (use --overwrite to replace).')
        else:
            print('\t\tUnpacking EnSight...')
            unpackEnsight(tarPath, dest)


def retrieveTrialParent(job, trialName):
    '''
    Restores the shared, parent-level data of a ride-height trial (e.g. caseSetup)
    from the archive back into the local trial directory. Existing local files are
    skipped unless --overwrite is given.
    '''
    archiveTrialDir = os.path.join(casesRoot(job), trialName)
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
    root = casesRoot(job)
    if not os.path.exists(root):
        print(f"No archive found for job {job} at {root}.")
        return [], []
    cases = []
    rideHeightParents = []
    for trial in trialRange:
        trialName = f"{trial:03d}"
        trialDir = os.path.join(root, trialName)
        if os.path.isdir(trialDir):
            #a ride-height container holds NNN_M cases (dirs or .tar.gz); a plain case doesn't
            points = {stripArchiveExt(item) for item in os.listdir(trialDir)
                      if item.startswith(trialName + '_')}
            if len(points) > 0:
                for point in sorted(points):
                    cases.append(os.path.join(trialName, point))
                rideHeightParents.append(trialName)
            else:
                cases.append(trialName)
        elif os.path.exists(trialDir + legacyExt):
            cases.append(trialName)
        else:
            print(f"Case {trialName} not found in archive. Skipping.")
    return cases, rideHeightParents


def stripArchiveExt(name):
    return name[:-len(legacyExt)] if name.endswith(legacyExt) else name


#------------------------------------------------------------------------------
# convert old uncompressed archives (--compress)
#------------------------------------------------------------------------------

def compressArchivedCases(job, trialRange):
    '''
    Converts old uncompressed archived cases to the current layout: EnSight is
    compressed to EnSight.tar.zst and processor dirs, polyMesh and ANSA files are
    deleted from the archive. Old .tar.gz archives and cases already in the current
    layout are skipped.
    '''
    cases, rideHeightParents = collectArchivedCases(job, trialRange)
    root = casesRoot(job)

    toConvert = []
    for casePath in cases:
        archivePath = os.path.join(root, casePath)
        if os.path.exists(archivePath + legacyExt):
            print(f"\t{casePath} is an old .tar.gz archive. Leaving as is.")
        elif not os.path.isdir(archivePath):
            print(f"\t{casePath} not found in archive. Skipping.")
        elif not isLegacyDir(archivePath):
            print(f"\t{casePath} already in the current layout. Skipping.")
        else:
            toConvert.append(casePath)

    if len(toConvert) == 0:
        print("No old uncompressed archived cases to convert.")
        return

    rows = []
    for casePath in toConvert:
        archivePath = os.path.join(root, casePath)
        sizes = getCaseSizes(archivePath)
        total = sizes['plain'] + sizes['ensight'] + sizes['dropped']
        rows.append((casePath, humanSize(total), humanSize(sizes['plain']),
                     humanSize(sizes['ensight']), countProcessors(archivePath), 'convert'))
    printSummary('ARCHIVE CONVERSION - job %s' % (job),
                 ('CASE', 'TOTAL', 'PLAIN', 'ENSIGHT', 'PROCS', 'ACTION'), rows,
                 'ENSIGHT is compressed; processor dirs, polyMesh and ANSA files are deleted '
                 'from the archive.')
    if not args.batchMode and not confirm('Proceeding will compress EnSight and delete processor dirs, '
                                          'polyMesh and ANSA files from the archived cases above.'):
        print("Conversion cancelled by user. No changes made.")
        return

    for casePath in toConvert:
        archivePath = os.path.join(root, casePath)
        print(f"\tConverting {casePath}...")
        try:
            items = [i for i in ensightItems(archivePath) if i != ensightTarName and not i.endswith(partialExt)]
            if items:
                print('\t\tCompressing EnSight...')
                makeZstTar(archivePath, items, os.path.join(archivePath, ensightTarName))
                for item in items:
                    removePath(os.path.join(archivePath, item))
            print('\t\tRemoving processor dirs, polyMesh and ANSA files...')
            removeDropped(archivePath)
        except (subprocess.CalledProcessError, OSError) as e:
            print('\t\tERROR! Converting %s failed (%s). Archived data left as is.' % (casePath, e))
    print("Conversion complete.")


def removeDropped(casePath):
    for item in os.listdir(casePath):
        if classifyPath([item]) == 'dropped':
            removePath(os.path.join(casePath, item))
    polyMeshPath = os.path.join(casePath, 'constant', 'polyMesh')
    if os.path.exists(polyMeshPath):
        shutil.rmtree(polyMeshPath)


#------------------------------------------------------------------------------
# clean in place (no mode flag)
#------------------------------------------------------------------------------

def clearCases(casesToProcess):
    actions = {}
    rows = []
    for casePath in casesToProcess:
        sizes = getCaseSizes(casePath)
        total = sizes['plain'] + sizes['ensight'] + sizes['dropped']
        if not hasEnsightDir(casePath):
            actions[casePath] = 'skip: no EnSight'
            freed = 0
        else:
            actions[casePath] = 'clean'
            freed = sizes['dropped'] + (sizes['ensight'] if args.dropEnsight else 0)
        rows.append((casePath, humanSize(total), humanSize(freed),
                     humanSize(sizes['ensight']), countProcessors(casePath), actions[casePath]))
    removed = 'processor dirs, polyMesh and ANSA files'
    if args.dropEnsight:
        removed += ' and EnSight'
    printSummary('CLEAN SUMMARY (local cases)',
                 ('CASE', 'TOTAL', 'FREED', 'ENSIGHT', 'PROCS', 'ACTION'), rows,
                 'Deletes %s. Cases with no EnSight output are left untouched.' % (removed))

    message = 'Proceeding will permanently delete the %s of the cases above.' % (removed)
    if args.dropEnsight:
        message += ' Without EnSight and processor dirs, no field results will remain.'
    if not args.batchMode and not confirm(message, 'Proceed with cleaning? [y/N]: '):
        print("Clean cancelled by user. No changes made.")
        return

    for casePath in casesToProcess:
        if actions[casePath] != 'clean':
            print('\tSkipping %s (no EnSight output; run may have failed).' % (casePath))
            continue
        print('\tCleaning %s...' % (casePath))
        clearCase(casePath)


def clearCase(casePath):
    print('\t\tRemoving processor dirs, polyMesh and ANSA files...')
    removeDropped(casePath)
    if args.dropEnsight:
        print('\t\tRemoving EnSight files...')
        for item in ensightItems(casePath):
            removePath(os.path.join(casePath, item))


#------------------------------------------------------------------------------
# shared helpers
#------------------------------------------------------------------------------

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
                cases.append(os.path.join(trialName, point))
            rideHeightParents.append(trialName)
        else:
            cases.append(trialName)
    return cases, rideHeightParents


def collectRideHeightPoints(trialName):
    if not os.path.exists(trialName):
        print(f"ERROR! Trial directory {trialName} does not exist.")
        return []
    return [item for item in os.listdir(trialName)
            if item.startswith(trialName + "_") and os.path.isdir(os.path.join(trialName, item))]


def countProcessors(casePath):
    return len([d for d in os.listdir(casePath) if d.startswith('processor')])


#column layout shared by the summary header and rows
_summaryRowFmt = '{:<20}{:>10}{:>10}{:>10}{:>7}{:>20}'


def printSummary(title, header, rows, note):
    print('\n' + '=' * 77)
    print('  ' + title)
    print('=' * 77)
    print(_summaryRowFmt.format(*header))
    print('-' * 77)
    for row in rows:
        print(_summaryRowFmt.format(*row))
    print('-' * 77)
    print('  ' + note)
    print('=' * 77)


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


def confirm(message, prompt='Proceed with archiving? [y/N]: '):
    print('\n' + message)
    while True:
        response = input(prompt).strip().lower()
        if response in ('y', 'yes'):
            return True
        if response in ('', 'n', 'no'):
            return False
        print("Please answer 'y' or 'n'.")


def getJob():
    '''
    Returns the job name and job path based on the current working directory.
    '''
    cwd = os.getcwd()
    if os.path.basename(cwd) != "CASES":
        sys.exit('ERROR! Please run in CASES directory!')
    parent = os.path.dirname(cwd)
    return os.path.basename(parent), parent


if __name__ == '__main__':
    main()
