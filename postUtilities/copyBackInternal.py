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

#copies the run data of internal cases back to a Google Drive folder: <folder>/<JOB>/<case>/ with the case's own
#sub-paths (postProcessing/<fo>/<time>/*.dat keep their folders). logs, Slurm .out files, setup files, summary and
#plots, function object .dat files, the system folder, images and videos. only new or changed files are uploaded (md5 against Drive)
#usage: python3 copyBackInternal.py [case ...] [--dryRun] [--noImages] [--config file] [--folderId id]
#
#settings (arguments override environment variables, which override the file):
#  <CASESETUP_CONFIG or ~/.config/caseSetup>/drive.ini (or --config)   [drive] folderID = ... (id or folder URL), credentials = ... (optional)
#  INTERNAL_DRIVE_FOLDER, INTERNAL_DRIVE_CREDENTIALS
#  credentials default: ~/.config/caseSetup/internal_cfd_credentials.json (do not upload this to git!)
#the folder has to be on a Shared Drive (service accounts have no storage of their own) and shared with the service
#account as a Content manager / Editor. needs google-api-python-client and google-auth

import os
import re
import sys
import glob
import hashlib
import argparse
import configparser

#settings folder: CASESETUP_CONFIG (set by postScript from POST_CONFIG, e.g. a folder on a disk the compute nodes
#share), else ~/.config/caseSetup
CONFIG_DIR = os.getenv('CASESETUP_CONFIG') or os.path.join('~', '.config', 'caseSetup')
CONFIG_FILE = os.path.join(CONFIG_DIR, 'drive.ini')
DEFAULT_CREDENTIALS = os.path.join(CONFIG_DIR, 'internal_cfd_credentials.json')
FOLDER_TYPE = 'application/vnd.google-apps.folder'

#files of a case to copy back, relative to the case directory
CASE_FILES = ['log.*', '*.out', 'internalSetup', 'fullCaseSetupDict', 'pvPostInternalSetup', 'summary.csv',
 '*.png',
 os.path.join('constant','triSurface','*.stl*',),
 os.path.join('constant','triSurface','*.obj*',)]
DATA_FILES = [os.path.join('postProcessing', '**', '*.dat'),
              os.path.join('system', '**', '*')]
IMAGE_FILES = [os.path.join('postProcessing', 'images', '**', '*.png'),
               os.path.join('postProcessing', 'meshImages', '**', '*.png')]
MOVIE_FILES = [os.path.join('postProcessing', 'movies', '*.mp4'), os.path.join('postProcessing', 'meshMovies', '*.mp4')]


def parseArgs(argv=None):
    parser = argparse.ArgumentParser(prog='copyBackInternal', description='Copy internal run data back to Google Drive')
    parser.add_argument('cases', nargs='*', default=['.'], help='case directories (default: current directory)')
    parser.add_argument('--folderId', default=None, help='Drive folder id or URL (overrides INTERNAL_DRIVE_FOLDER / drive.ini)')
    parser.add_argument('--config', default=None, help='settings file (default: ~/.config/caseSetup/drive.ini)')
    parser.add_argument('--noImages', action='store_true', help='skip the png images (videos, plots and data still go)')
    parser.add_argument('--dryRun', action='store_true', help='show what would be uploaded, change nothing')
    return parser.parse_args(argv)


#----- settings

def cleanFolderId(value):
    #the id after /folders/ of a folder URL; quotes and a pasted URL are accepted
    value = str(value).strip().strip('\'"').strip()
    match = re.search(r'/folders/([A-Za-z0-9_-]+)', value)
    return match.group(1) if match else value


def getSettings(folderId=None, configFile=None):
    configPath = os.path.expanduser(os.path.expandvars(configFile or CONFIG_FILE))
    if configFile and not os.path.isfile(configPath):
        sys.exit('ERROR! Config file %s not found' % (configPath))
    config = configparser.ConfigParser()
    if os.path.isfile(configPath):
        config.read(configPath)
    fromFile = config['drive'] if config.has_section('drive') else {}
    settings = {
        'folderID': cleanFolderId(folderId or os.getenv('INTERNAL_DRIVE_FOLDER') or fromFile.get('folderID', '')),
        'credentials': os.path.expanduser(os.path.expandvars(
            os.getenv('INTERNAL_DRIVE_CREDENTIALS') or fromFile.get('credentials', DEFAULT_CREDENTIALS))),
        'configPath': configPath}
    if not settings['folderID']:
        sys.exit('ERROR! No Drive folder: set INTERNAL_DRIVE_FOLDER, --folderId or folderID in %s (found the file: %s)'
                 % (configPath, 'yes' if os.path.isfile(configPath) else 'no'))
    return settings


#----- local files

def jobAndCase(caseDir):
    #<job>/CASES/<case> -> (job, case); a case outside CASES goes under its own name
    caseDir = os.path.abspath(caseDir)
    case = os.path.basename(caseDir)
    if os.path.basename(os.path.dirname(caseDir)) == 'CASES':
        return os.path.basename(os.path.dirname(os.path.dirname(caseDir))), case
    return '', case


def collectFiles(caseDir, images=True):
    #relative paths (forward slashes) of the files to copy back
    patterns = CASE_FILES + DATA_FILES + MOVIE_FILES + (IMAGE_FILES if images else [])
    found = set()
    for pattern in patterns:
        for path in glob.glob(os.path.join(caseDir, pattern), recursive=True):
            if os.path.isfile(path):
                found.add(os.path.relpath(path, caseDir).replace('\\', '/'))
    return sorted(found)


def md5(path):
    digest = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def planUploads(local, remote):
    #local {relpath: md5}, remote {relpath: md5 or None} -> (new, changed, unchanged) relpaths
    new = sorted(path for path in local if path not in remote)
    changed = sorted(path for path in local if path in remote and remote[path] != local[path])
    unchanged = sorted(path for path in local if path in remote and remote[path] == local[path])
    return new, changed, unchanged


#----- Drive

class Drive:
    #folders and files under one Drive folder, listings cached per folder
    def __init__(self, settings):
        try:
            from google.oauth2 import service_account
            from googleapiclient.discovery import build
        except ImportError:
            sys.exit('ERROR! google-api-python-client / google-auth missing: pip install --user google-api-python-client google-auth')
        if not os.path.isfile(settings['credentials']):
            sys.exit('ERROR! Credentials file %s not found: set credentials in %s or INTERNAL_DRIVE_CREDENTIALS'
                     % (settings['credentials'], settings['configPath']))
        creds = service_account.Credentials.from_service_account_file(
            settings['credentials'], scopes=['https://www.googleapis.com/auth/drive'])
        self.email = creds.service_account_email
        self.service = build('drive', 'v3', credentials=creds, cache_discovery=False)
        self.root = settings['folderID']
        self.listings = {}

    def listFolder(self, folderId):
        #{name: {id, md5Checksum, mimeType}} of a folder, all pages
        if folderId in self.listings:
            return self.listings[folderId]
        items, token = {}, None
        while True:
            result = self.service.files().list(
                q="'%s' in parents and trashed = false" % (folderId),
                fields='nextPageToken, files(id, name, md5Checksum, mimeType)', pageSize=1000, pageToken=token,
                includeItemsFromAllDrives=True, supportsAllDrives=True).execute()
            for item in result.get('files', []):
                items[item['name']] = item
            token = result.get('nextPageToken')
            if not token:
                break
        self.listings[folderId] = items
        return items

    def folder(self, parentId, name, create=True):
        #id of a sub-folder, made if missing (None when create is False and it does not exist)
        item = self.listFolder(parentId).get(name)
        if item and item['mimeType'] == FOLDER_TYPE:
            return item['id']
        if not create:
            return None
        made = self.service.files().create(body={'name': name, 'mimeType': FOLDER_TYPE, 'parents': [parentId]},
                                           fields='id', supportsAllDrives=True).execute()
        self.listings[parentId][name] = {'id': made['id'], 'mimeType': FOLDER_TYPE}
        self.listings[made['id']] = {}
        return made['id']

    def folderOf(self, baseId, relDir, create=True):
        folderId = baseId
        for part in [p for p in relDir.split('/') if p]:
            folderId = self.folder(folderId, part, create)
            if folderId is None:
                return None
        return folderId

    def remoteMd5(self, baseId, relpaths):
        #{relpath: md5} of the files already on Drive (folders that do not exist yet are not made)
        remote = {}
        for relpath in relpaths:
            relDir, name = os.path.split(relpath)
            folderId = self.folderOf(baseId, relDir, create=False) if baseId else None
            if folderId is None:
                continue
            item = self.listFolder(folderId).get(name)
            if item and item['mimeType'] != FOLDER_TYPE:
                remote[relpath] = item.get('md5Checksum')
        return remote

    def upload(self, baseId, relpath, path):
        from googleapiclient.http import MediaFileUpload
        relDir, name = os.path.split(relpath)
        folderId = self.folderOf(baseId, relDir)
        media = MediaFileUpload(path, resumable=os.path.getsize(path) > (5 << 20))
        existing = self.listFolder(folderId).get(name)
        if existing and existing['mimeType'] != FOLDER_TYPE:
            self.service.files().update(fileId=existing['id'], media_body=media, supportsAllDrives=True).execute()
        else:
            made = self.service.files().create(body={'name': name, 'parents': [folderId]}, media_body=media,
                                               fields='id', supportsAllDrives=True).execute()
            self.listings[folderId][name] = {'id': made['id'], 'mimeType': 'file'}


def copyCase(drive, caseDir, images=True, dryRun=False):
    caseDir = os.path.abspath(caseDir)
    job, case = jobAndCase(caseDir)
    relCase = '/'.join(p for p in [job, case] if p)
    files = collectFiles(caseDir, images)
    if not files:
        print('\t%s: nothing to copy' % (relCase))
        return 0
    local = dict((relpath, md5(os.path.join(caseDir, relpath))) for relpath in files)
    baseId = drive.folderOf(drive.root, relCase, create=False)
    remote = drive.remoteMd5(baseId, files)
    new, changed, unchanged = planUploads(local, remote)
    size = sum(os.path.getsize(os.path.join(caseDir, p)) for p in new + changed)
    print('\t%s: %d files, %d new, %d changed, %d unchanged (%.1f MB to upload)'
          % (relCase, len(files), len(new), len(changed), len(unchanged), size/1e6))
    if dryRun:
        for relpath in new + changed:
            print('\t\t%s %s' % ('new    ' if relpath in new else 'changed', relpath))
        return 0
    baseId = drive.folderOf(drive.root, relCase)
    for k, relpath in enumerate(new + changed, start=1):
        drive.upload(baseId, relpath, os.path.join(caseDir, relpath))
        if k % 50 == 0 or k == len(new) + len(changed):
            print('\t\t%d/%d uploaded' % (k, len(new) + len(changed)))
    return len(new) + len(changed)


def main(argv=None):
    args = parseArgs(argv)
    settings = getSettings(args.folderId, args.config)
    drive = Drive(settings)
    print('\n\tCopying back to Drive folder %s as %s...' % (settings['folderID'], drive.email))
    try:
        from googleapiclient.errors import HttpError
    except ImportError:
        HttpError = Exception
    try:
        drive.listFolder(drive.root)
        for case in args.cases:
            copyCase(drive, case, images=not args.noImages, dryRun=args.dryRun)
    except HttpError as error:
        text = str(error)
        if 'storageQuotaExceeded' in text or 'do not have storage quota' in text:
            sys.exit('ERROR! Service accounts have no storage: put the folder on a Shared Drive.\n%s' % (text))
        if '404' in text:
            sys.exit('ERROR! Folder %s not found: check the id and share it with %s.\n%s'
                     % (settings['folderID'], drive.email, text))
        if '403' in text:
            sys.exit('ERROR! Permission denied (or the Drive API is not enabled for the key\'s project): share the '
                     'folder with %s as Content manager / Editor.\n%s' % (drive.email, text))
        raise


if __name__ == '__main__':
    main()
