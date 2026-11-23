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

#internal flow / Stirling engine case setup, reads the internalSetup file in the case directory
#phase 3: closed volume with WALL geometry and MOVE bodies (displacementInterpolation), rhoPimpleFoam

import os
import sys
import shutil
import argparse

import setupCore
from setupCore import __version__
import writeInternal

#available addon keywords (REGN, SMP come in a later phase)
addonKeyWords = ['WALL', 'MOVE', 'FAN', 'POR', 'REGN', 'SMP']
supportedPrefixes = ['WALL', 'MOVE', 'FAN', 'POR']
geomColumnNames = ['scale', 'refinement', 'layers', 'expansion', 'wallmodel']
configName = 'internalSetup'


def getExecDir():
    #directory of the tool (or of the executable when built with PyInstaller)
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.realpath(__file__))


def parseArgs(argv=None):
    parser = argparse.ArgumentParser(prog='internalSetup-v%s' % __version__,
                                     description='Sets up an internal flow (Stirling engine) case from the internalSetup file')
    parser.add_argument('-s', '--setup', default='internal',
                        help='Identifies which setup templates to use. Available setups are: internal, internal-otr')
    parser.add_argument('-d', '--controlDict', action='store_true',
                        help='Writes only controlDict for the solver.')
    parser.add_argument('--new', action='store_true',
                        help='Writes a new internalSetup, will overwrite what is currently there if it exists!')
    parser.add_argument('--modules', action='store_true',
                        help='Shows all possible modules.')
    return parser.parse_args(argv)


def findFromSections(writeCaseSetupDict, geomNames, pidDict=None):
    #MOVE sections without their own file (FROM = <geometry>, like a piston crown in the chamber)
    #FROM has to point at a geometry still in GEOM, otherwise it's left over from a removed part
    #with pidDict: ZONE_PIDS also has to match a PID of that geometry, otherwise the PID was removed
    sections = []
    for s in writeCaseSetupDict:
        source = writeCaseSetupDict[s].get('FROM', ['self'])[0]
        if not s.startswith('MOVE') or s in geomNames or source not in geomNames:
            continue
        keywords = [k for k in writeCaseSetupDict[s].get('ZONE_PIDS', ['']) if k != '']
        if pidDict is not None and keywords and source in pidDict:
            if not any(kw.lower() in pid.lower() for kw in keywords for pid in pidDict[source]):
                print('\t\t[%s] ZONE_PIDS %s matches no PID of %s any more, section removed' % (s, ' '.join(keywords), source))
                continue
        sections.append(s)
    return sections


def findPidBodies(pidDict):
    #MOVE-<name> PIDs in files that are not MOVE files: one body per PID, {section: (geometry, pid)}
    bodies = {}
    for geomName, pids in pidDict.items():
        if geomName.startswith('MOVE'):
            continue
        for pid in pids:
            if not pid.upper().startswith('MOVE-') or len(pid) <= len('MOVE-'):
                continue
            section = 'MOVE-%s' % (pid[len('MOVE-'):])
            if section in bodies:
                sys.exit('ERROR! PID %s is in %s and %s: MOVE- PIDs need unique names'
                         % (pid, bodies[section][0], geomName))
            bodies[section] = (geomName, pid)
    return bodies


def addPidBodies(writeCaseSetupDict, pidBodies, geomNames, moveItems):
    #new [MOVE-<name>] section for each MOVE- PID: FROM, ZONE_PIDS and MOTION_COLUMN filled, the rest from defaultMove
    #an existing section is kept as it is (only missing keys are added later)
    updateFlag = False
    for section, (geomName, pid) in pidBodies.items():
        if section in geomNames:
            print('\t\tWARNING: PID %s in %s has the same name as the geometry %s, no section added for it'
                  % (pid, geomName, section))
            continue
        if section in writeCaseSetupDict:
            continue
        values = {'FROM': [geomName], 'ZONE_PIDS': [pid], 'MOTION_COLUMN': [section[len('MOVE-'):]]}
        writeCaseSetupDict[section] = dict((key, values.get(key, [val])) for key, val in moveItems)
        print('\n\tAdded [%s] for the PID %s in %s (FROM, ZONE_PIDS, MOTION_COLUMN set), please check it'
              % (section, pid, geomName))
        updateFlag = True
    return updateFlag


def copyPostSetup(execDir, caseDir):
    #copying pvPostInternalSetup into the case, an existing copy is kept
    dest = os.path.join(caseDir, 'pvPostInternalSetup')
    if os.path.isfile(dest):
        print('\n\tpvPostInternalSetup exists in the case, kept')
        return
    shutil.copyfile(os.path.join(execDir, 'postUtilities', 'default', 'pvPostInternalSetup'), dest)
    print('\n\tCopied pvPostInternalSetup (image panels and videos for postScript) into the case')


def checkGeomPrefixes(geomDict):
    for geom in geomDict:
        prefix = geom.split('-')[0]
        if prefix not in addonKeyWords:
            sys.exit('ERROR! %s has no internal prefix, name geometry files %s-<name>.stl'
                     % (geom, '/'.join(addonKeyWords)))
        if prefix not in supportedPrefixes:
            sys.exit('ERROR! %s: %s geometry is not supported yet (available: %s)'
                     % (geom, prefix, '/'.join(supportedPrefixes)))


def main(argv=None):
    args = parseArgs(argv)
    titleText = '''\t##############################\n\t\tinternalSetup-v%s\t\n\t##############################''' % __version__
    print(titleText)

    execDir = getExecDir()
    templateLoc = os.path.join(execDir, 'setupTemplates', args.setup)
    if not os.path.isdir(templateLoc):
        sys.exit('ERROR! %s is not a valid setup template location!' % (templateLoc))
    print('\n\tUsing %s setup template in path: %s' % (args.setup, templateLoc))

    #the tool runs from inside the case: <job>/CASES/<case>, geometry in <job>/02_reference/MSH
    caseDir = os.getcwd()
    jobPath = os.path.abspath(os.path.join(caseDir, os.pardir, os.pardir))
    mshDir = os.path.join(jobPath, '02_reference', 'MSH')

    defaultDict = setupCore.getGlobalDefaults(templateLoc, args.modules)
    if args.new:
        setupCore.writeNewCaseSetup(defaultDict, caseDir, configName)
        return

    caseSetupDict, writeCaseSetupDict, fullCaseSetupDict, updateFlag = setupCore.getCaseSetup(
        defaultDict, caseDir, addonKeyWords, configName)

    if args.controlDict:
        writeInternal.writeControlDictInternal(templateLoc, fullCaseSetupDict)
        return

    #geometry: GEOM line -> geomDict, WALL sections filled from defaultWall (group lines kept)
    try:
        geometryList = fullCaseSetupDict['GEOMETRY']['GEOM']
    except KeyError:
        sys.exit('ERROR! [GEOMETRY] section invalid!')
    geomDict = setupCore.geomToDict({}, geometryList, geomColumnNames)
    checkGeomPrefixes(geomDict)
    for prefix, template, header in (('WALL', 'defaultWall', 'WALL_SETUP'), ('MOVE', 'defaultMove', 'MOVE_SETUP'),
                                     ('FAN', 'defaultFan', 'FAN_SETUP'), ('POR', 'defaultPorous', 'POR_SETUP')):
        writeCaseSetupDict, flag = setupCore.checkAddonDefaults(
            geomDict, writeCaseSetupDict, templateLoc, prefix, template, header, keepExtra=True)
        updateFlag = updateFlag or flag
    geomNames = [setupCore.stripExt(geom) for geom in geomDict]
    #MOVE- PIDs in WALL files become their own MOVE sections (FROM = that file)
    pidDict = writeInternal.readGeomPIDs(geomDict, mshDir)
    moveItems = setupCore.readAddonTemplate(templateLoc, 'defaultMove', 'MOVE_SETUP')
    if addPidBodies(writeCaseSetupDict, findPidBodies(pidDict), geomNames, moveItems):
        updateFlag = True
    fromSections = findFromSections(writeCaseSetupDict, geomNames, pidDict)
    for section in fromSections:
        if setupCore.fillAddonSection(writeCaseSetupDict, section, moveItems, keepExtra=True):
            updateFlag = True
    setupCore.cleanUpCaseSetup(geomDict, writeCaseSetupDict, fullCaseSetupDict, addonKeyWords, fromSections)
    setupCore.writeToCaseSetup(writeCaseSetupDict, caseDir, configName, updateFlag)

    #PIDs, mesh motion and BC groups, checked before anything is written to the case
    writeInternal.getFlow(fullCaseSetupDict)
    motionInfo = writeInternal.buildMotion(geomDict, fullCaseSetupDict, pidDict, mshDir, caseDir)
    if motionInfo and writeInternal.isSteady(fullCaseSetupDict):
        sys.exit('ERROR! %s is a steady solver: MOVE bodies need rhoPimpleFoam or pimpleFoam'
                 % (fullCaseSetupDict['GLOBAL_SIM_CONTROL']['SOLVER'][0]))
    movingPids = motionInfo['movingPids'] if motionInfo else {}
    wallGroups = writeInternal.getWallGroups(geomDict, fullCaseSetupDict, pidDict, movingPids)

    setupCore.linkGeomFiles(geomDict, caseDir, jobPath)
    if motionInfo:
        writeInternal.writeMotionGeometry(motionInfo, mshDir)
        writeInternal.writeRemeshMeshes(motionInfo, mshDir)
    writeInternal.writeBlockMeshInternal(templateLoc, geomDict, fullCaseSetupDict)
    writeInternal.writeSnappyInternal(templateLoc, geomDict, fullCaseSetupDict, motionInfo)
    writeInternal.writeCasePropertiesInternal(templateLoc, geomDict, fullCaseSetupDict, wallGroups)
    writeInternal.writeThermoInternal(templateLoc, fullCaseSetupDict)
    fanNames = writeInternal.writeFanBaffles(geomDict, fullCaseSetupDict)
    porous = writeInternal.getPorous(geomDict, fullCaseSetupDict)
    porOwners = motionInfo['porOwners'] if motionInfo else {}
    writeInternal.writeFvOptionsInternal(porous, fullCaseSetupDict, porOwners)
    warmStart = writeInternal.getWarmStart(fullCaseSetupDict, porous)
    writeInternal.writeTurbulenceInternal(templateLoc, fullCaseSetupDict)
    writeInternal.writeControlDictInternal(templateLoc, fullCaseSetupDict)
    writeInternal.writeNumericsInternal(templateLoc, fullCaseSetupDict, writeInternal.hasPressureBC(wallGroups))
    writeInternal.writeFunctionObjectsInternal(wallGroups, fullCaseSetupDict, motionInfo, porous)
    motionZones = writeInternal.writeMotionDicts(templateLoc, motionInfo, fullCaseSetupDict) if motionInfo else None
    writeInternal.makeScriptsInternal(templateLoc, fullCaseSetupDict, execDir, motionZones, fanNames, motionInfo,
                                      warmStart)

    copyPostSetup(execDir, caseDir)
    fullCaseSetupDict = writeInternal.addSetupInfo(fullCaseSetupDict, wallGroups, motionInfo)
    setupCore.writeToCaseSetup(fullCaseSetupDict, caseDir, 'fullCaseSetupDict')
    print('\nCase setup completed successfully! runScript does mesh + solve + export in one job')
    if motionInfo:
        print('Mesh motion dry run first: ./runScript --meshOnly then ./motionCheckScript (cluster: ./submitScript --check)')
        if motionInfo.get('segments'):
            print('Remeshing: runScript also makes the meshes in meshes/seg*, the solve maps onto the next one at each '
                  'segment end (system/remeshSegments)')
    print('Run: ./runScript then ./postScript locally, or ./submitScript on the cluster (from the case directory in $SCRATCH, chains postScript)')


if __name__ == '__main__':
    main()
