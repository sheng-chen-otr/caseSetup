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

#pushes internal (Stirling) run summaries to one database gsheet: one row per run, found by
#JOB + NAME (added at the bottom if missing). header names in row 1; a summary key without a column gets a
#new column at the end. a run already in the sheet is overwritten (empty values clear old ones, unless
#--keepEmpty); COMMENTS is only written with --comment, columns not in the summary are left alone. reads summary.csv in the case (run postInternal.py first)
#usage: python3 gsheetInternal.py [case ...] [--comment "text"] [--dryRun] [--config file] [--sheetId id] [--worksheet name]
#
#settings (arguments override environment variables, which override the file):
#  <CASESETUP_CONFIG or ~/.config/caseSetup>/gsheet.ini (or --config)   [gsheet] sheetID = ..., worksheet = ... (optional, default first
#                                                 sheet), credentials = ... (optional)
#  INTERNAL_GSHEET_ID, INTERNAL_GSHEET_WORKSHEET, INTERNAL_GSHEET_CREDENTIALS
#  credentials default: ~/.config/caseSetup/internal_cfd_credentials.json (do not upload this to git!)

import argparse
import configparser
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import postInternal

#settings folder: CASESETUP_CONFIG (set by postScript from POST_CONFIG, e.g. a folder on a disk the compute nodes
#share), else ~/.config/caseSetup
CONFIG_DIR = os.getenv('CASESETUP_CONFIG') or os.path.join('~', '.config', 'caseSetup')
CONFIG_FILE = os.path.join(CONFIG_DIR, 'gsheet.ini')
DEFAULT_CREDENTIALS = os.path.join(CONFIG_DIR, 'internal_cfd_credentials.json')
HEADER_ROW = 1
KEY_COLUMNS = ['JOB', 'NAME']


def getSettings(sheetId=None, worksheet=None, configFile=None):
    #sheet id, worksheet name and credentials path from the arguments, environment or config file
    #configFile: --config path, else ~/.config/caseSetup/gsheet.ini (a given file has to exist)
    configPath = os.path.expanduser(os.path.expandvars(configFile or CONFIG_FILE))
    if configFile and not os.path.isfile(configPath):
        sys.exit('ERROR! Config file %s not found' % (configPath))
    config = configparser.ConfigParser()
    if os.path.isfile(configPath):
        config.read(configPath)
    fromFile = config['gsheet'] if config.has_section('gsheet') else {}
    settings = {
        'sheetID': sheetId or os.getenv('INTERNAL_GSHEET_ID') or fromFile.get('sheetID', ''),
        'worksheet': worksheet or os.getenv('INTERNAL_GSHEET_WORKSHEET') or fromFile.get('worksheet', ''),
        'credentials': os.path.expanduser(os.path.expandvars(
            os.getenv('INTERNAL_GSHEET_CREDENTIALS') or fromFile.get('credentials', DEFAULT_CREDENTIALS)))}
    settings['sheetID'] = cleanSheetId(settings['sheetID'])
    settings['worksheet'] = settings['worksheet'].strip().strip('\'"')
    settings['configPath'] = configPath
    if not settings['sheetID']:
        sys.exit('ERROR! No sheet id: set INTERNAL_GSHEET_ID, --sheetId or sheetID in %s (found the file: %s)'
                 % (configPath, 'yes' if os.path.isfile(configPath) else 'no'))
    return settings


def cleanSheetId(value):
    #the id between /d/ and the next / of the sheet URL; quotes and a pasted full URL are accepted
    value = str(value).strip().strip('\'"').strip()
    match = re.search(r'/spreadsheets/d/([A-Za-z0-9_-]+)', value)
    return match.group(1) if match else value


def openWorksheet(settings):
    try:
        import gspread
        from google.oauth2.service_account import Credentials
    except ImportError:
        sys.exit('ERROR! gspread / google-auth missing: pip install --user gspread google-auth')
    if not os.path.isfile(settings['credentials']):
        sys.exit('ERROR! Credentials file %s not found: set credentials in %s or INTERNAL_GSHEET_CREDENTIALS'
                 % (settings['credentials'], settings.get('configPath', CONFIG_FILE)))
    try:
        creds = Credentials.from_service_account_file(settings['credentials'],
                                                      scopes=['https://www.googleapis.com/auth/spreadsheets'])
    except ValueError as error:
        sys.exit('ERROR! %s is not a service account key (json from Google Cloud > IAM > Service accounts > Keys): %s'
                 % (settings['credentials'], error))
    print('\tservice account %s, sheet %s' % (creds.service_account_email, settings['sheetID']))
    try:
        spreadsheet = gspread.authorize(creds).open_by_key(settings['sheetID'])
        if settings['worksheet']:
            return spreadsheet.worksheet(settings['worksheet'])
        return spreadsheet.sheet1
    except gspread.exceptions.SpreadsheetNotFound:
        sys.exit('ERROR! Sheet %s not found: check the id (between /d/ and /edit in the URL) and share the sheet with '
                 '%s as Editor' % (settings['sheetID'], creds.service_account_email))
    except gspread.exceptions.WorksheetNotFound:
        sys.exit('ERROR! No tab "%s" in the sheet: check worksheet (tab names, case and spaces count)' % (settings['worksheet']))
    except gspread.exceptions.APIError as error:
        text = str(error)
        if 'has not been used' in text or 'is disabled' in text or 'SERVICE_DISABLED' in text:
            sys.exit('ERROR! The Google Sheets API is not enabled for the key\'s Google Cloud project: enable it under '
                     'APIs & Services, wait a minute, retry.\n%s' % (text))
        if 'PERMISSION_DENIED' in text or '403' in text:
            sys.exit('ERROR! Permission denied: share the sheet with %s as Editor.\n%s' % (creds.service_account_email, text))
        if '404' in text:
            sys.exit('ERROR! Sheet %s not found: check the id and the sharing with %s.\n%s'
                     % (settings['sheetID'], creds.service_account_email, text))
        raise


def normalise(name):
    return str(name).strip().upper()


def columnLetter(index):
    #0-based column index -> A, B, .., AA
    letters, n = '', index + 1
    while n:
        n, r = divmod(n - 1, 26)
        letters = chr(65 + r) + letters
    return letters


def cellValue(value):
    #numbers stay numbers in the sheet, nan/inf are left empty
    if isinstance(value, float) and (value != value or value in (float('inf'), float('-inf'))):
        return ''
    return value


def sameKey(a, b):
    #JOB / NAME match: spaces and case ignored, and 003 = 3 (a number the sheet stored without its leading zeros)
    a, b = str(a).strip(), str(b).strip()
    if a.lower() == b.lower():
        return True
    try:
        return float(a) == float(b)
    except ValueError:
        return False


def planUpdate(headers, rows, summary, comment=None, keepEmpty=False):
    #pure logic, no API: (new headers, row number or None for a new row, {column index: value}, added columns)
    #headers: header row values, rows: sheet rows below the header (lists of strings)
    #an existing run is overwritten: every summary column gets the new value, empty/nan clears the old one
    #(keepEmpty: leave old values where the new one is empty); COMMENTS and columns not in the summary are kept
    headers = list(headers)
    while headers and str(headers[-1]).strip() == '':
        headers.pop()
    if not any(str(h).strip() for h in headers):
        headers = list(postInternal.FIXED_KEYS)
    nOld = len(headers)
    index = {normalise(h): k for k, h in enumerate(headers) if str(h).strip()}
    values = dict(summary)
    values.pop('COMMENTS', None)
    if comment is not None:
        values['COMMENTS'] = comment
    for key in list(values):
        if normalise(key) not in index:
            if values[key] == '' or values[key] is None:
                continue
            index[normalise(key)] = len(headers)
            headers.append(key)
    for key in KEY_COLUMNS:
        if normalise(key) not in index:
            index[normalise(key)] = len(headers)
            headers.append(key)

    #row of this run: JOB and NAME match (the first one if the sheet has duplicates)
    jobCol, nameCol = index['JOB'], index['NAME']
    rowNumber = None
    for k, row in enumerate(rows):
        cells = list(row) + [''] * (len(headers) - len(row))
        if sameKey(cells[jobCol], values.get('JOB', '')) and sameKey(cells[nameCol], values.get('NAME', '')):
            rowNumber = HEADER_ROW + 1 + k
            break
    updates = {}
    for key, value in values.items():
        value = cellValue(value)
        if value is None:
            value = ''
        if value == '' and (rowNumber is None or keepEmpty or normalise(key) not in index):
            continue
        updates[index[normalise(key)]] = value
    #identity columns as text, so 003 stays 003
    for key in KEY_COLUMNS:
        if values.get(key, '') != '':
            updates[index[key]] = str(values[key])
    return headers, rowNumber, updates, headers[nOld:]


def pushSummary(worksheet, summary, comment=None, dryRun=False, keepEmpty=False):
    allValues = worksheet.get_all_values()
    headers = allValues[HEADER_ROW - 1] if len(allValues) >= HEADER_ROW else []
    rows = allValues[HEADER_ROW:]
    newHeaders, rowNumber, updates, added = planUpdate(headers, rows, summary, comment, keepEmpty)
    if rowNumber is None:
        #first row below the last non-empty row
        used = [k for k, row in enumerate(rows) if any(str(c).strip() for c in row)]
        rowNumber = HEADER_ROW + 1 + (used[-1] + 1 if used else 0)
        action = 'new row %d' % (rowNumber)
    else:
        action = 'row %d (overwritten)' % (rowNumber)
    print('\t%s / %s -> %s, %d values%s' % (summary.get('JOB', ''), summary.get('NAME', ''), action, len(updates),
                                           ', new columns: %s' % (' '.join(added)) if added else ''))
    if dryRun:
        for col in sorted(updates):
            print('\t\t%-4s %-30s %s' % (columnLetter(col), newHeaders[col], updates[col]))
        return
    if len(newHeaders) > worksheet.col_count:
        worksheet.add_cols(len(newHeaders) - worksheet.col_count)
    if rowNumber > worksheet.row_count:
        worksheet.add_rows(rowNumber - worksheet.row_count)
    batch = []
    if newHeaders != headers[:len(newHeaders)]:
        batch.append({'range': 'A%d:%s%d' % (HEADER_ROW, columnLetter(len(newHeaders) - 1), HEADER_ROW),
                      'values': [newHeaders]})
    for col, value in sorted(updates.items()):
        batch.append({'range': '%s%d' % (columnLetter(col), rowNumber), 'values': [[value]]})
    worksheet.batch_update(batch, value_input_option='RAW')


def main(argv=None):
    parser = argparse.ArgumentParser(description='Push internal run summaries to the database gsheet')
    parser.add_argument('cases', nargs='*', default=['.'], help='case directories (default: current directory)')
    parser.add_argument('--comment', default=None, help='text for the COMMENTS column (otherwise left as it is)')
    parser.add_argument('--sheetId', default=None, help='sheet id (overrides INTERNAL_GSHEET_ID / gsheet.ini)')
    parser.add_argument('--worksheet', default=None, help='worksheet name (default: first sheet)')
    parser.add_argument('--config', default=None, help='path to the config file (default: ~/.config/caseSetup/gsheet.ini)')
    parser.add_argument('--dryRun', action='store_true', help='show what would be written, change nothing')
    parser.add_argument('--keepEmpty', action='store_true',
                        help='for a run already in the sheet, keep old values where the new summary is empty')
    args = parser.parse_args(argv)
    settings = getSettings(args.sheetId, args.worksheet, args.config)
    worksheet = openWorksheet(settings)
    print('\n\tPushing to "%s"...' % (worksheet.title))
    for case in args.cases:
        pushSummary(worksheet, postInternal.readSummary(case), args.comment, args.dryRun, args.keepEmpty)


if __name__ == '__main__':
    main()
