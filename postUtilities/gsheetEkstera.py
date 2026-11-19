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

#pushes internal (Ekstera / Stirling) run summaries to one database gsheet: one row per run, found by
#JOB + NAME (added at the bottom if missing). header names in row 1; a summary key without a column gets a
#new column at the end. empty values never overwrite what is in the sheet, COMMENTS is only written with
#--comment. reads postInternal/summary.csv (run postInternal.py first)
#usage: python3 gsheetEkstera.py [case ...] [--comment "text"] [--dryRun]
#
#settings (environment variables override the file):
#  ~/.config/ekstera/gsheet.ini     [gsheet] sheetID = ..., worksheet = ... (optional, default first sheet),
#                                   credentials = ... (optional)
#  EKSTERA_GSHEET_ID, EKSTERA_GSHEET_WORKSHEET, EKSTERA_GSHEET_CREDENTIALS
#  credentials default: ~/.config/ekstera/credentials.json (do not upload this to git!)

import argparse
import configparser
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import postInternal

CONFIG_FILE = os.path.join('~', '.config', 'ekstera', 'gsheet.ini')
DEFAULT_CREDENTIALS = os.path.join('~', '.config', 'ekstera', 'credentials.json')
HEADER_ROW = 1
KEY_COLUMNS = ['JOB', 'NAME']


def getSettings(sheetId=None, worksheet=None):
    #sheet id, worksheet name and credentials path from the arguments, environment or config file
    config = configparser.ConfigParser()
    configPath = os.path.expanduser(CONFIG_FILE)
    if os.path.isfile(configPath):
        config.read(configPath)
    fromFile = config['gsheet'] if config.has_section('gsheet') else {}
    settings = {
        'sheetID': sheetId or os.getenv('EKSTERA_GSHEET_ID') or fromFile.get('sheetID', ''),
        'worksheet': worksheet or os.getenv('EKSTERA_GSHEET_WORKSHEET') or fromFile.get('worksheet', ''),
        'credentials': os.path.expanduser(os.path.expandvars(
            os.getenv('EKSTERA_GSHEET_CREDENTIALS') or fromFile.get('credentials', DEFAULT_CREDENTIALS)))}
    if not settings['sheetID']:
        sys.exit('ERROR! No Ekstera sheet id: set EKSTERA_GSHEET_ID or sheetID in %s' % (CONFIG_FILE))
    return settings


def openWorksheet(settings):
    import gspread
    from google.oauth2.service_account import Credentials
    if not os.path.isfile(settings['credentials']):
        sys.exit('ERROR! Credentials file %s not found: set EKSTERA_GSHEET_CREDENTIALS' % (settings['credentials']))
    creds = Credentials.from_service_account_file(settings['credentials'],
                                                  scopes=['https://www.googleapis.com/auth/spreadsheets'])
    spreadsheet = gspread.authorize(creds).open_by_key(settings['sheetID'])
    if settings['worksheet']:
        return spreadsheet.worksheet(settings['worksheet'])
    return spreadsheet.sheet1


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


def planUpdate(headers, rows, summary, comment=None):
    #pure logic, no API: (new headers, row number or None for a new row, {column index: value}, added columns)
    #headers: header row values, rows: sheet rows below the header (lists of strings)
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

    #row of this run: JOB and NAME match
    jobCol, nameCol = index['JOB'], index['NAME']
    rowNumber = None
    for k, row in enumerate(rows):
        cells = list(row) + [''] * (len(headers) - len(row))
        if str(cells[jobCol]).strip() == str(values.get('JOB', '')).strip() and \
                str(cells[nameCol]).strip() == str(values.get('NAME', '')).strip():
            rowNumber = HEADER_ROW + 1 + k
            break
    updates = {}
    for key, value in values.items():
        value = cellValue(value)
        if value == '' or value is None:
            continue
        updates[index[normalise(key)]] = value
    return headers, rowNumber, updates, headers[nOld:]


def pushSummary(worksheet, summary, comment=None, dryRun=False):
    allValues = worksheet.get_all_values()
    headers = allValues[HEADER_ROW - 1] if len(allValues) >= HEADER_ROW else []
    rows = allValues[HEADER_ROW:]
    newHeaders, rowNumber, updates, added = planUpdate(headers, rows, summary, comment)
    if rowNumber is None:
        #first row below the last non-empty row
        used = [k for k, row in enumerate(rows) if any(str(c).strip() for c in row)]
        rowNumber = HEADER_ROW + 1 + (used[-1] + 1 if used else 0)
        action = 'new row %d' % (rowNumber)
    else:
        action = 'row %d' % (rowNumber)
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
    parser = argparse.ArgumentParser(description='Push internal run summaries to the Ekstera database gsheet')
    parser.add_argument('cases', nargs='*', default=['.'], help='case directories (default: current directory)')
    parser.add_argument('--comment', default=None, help='text for the COMMENTS column (otherwise left as it is)')
    parser.add_argument('--sheetId', default=None, help='sheet id (overrides EKSTERA_GSHEET_ID / gsheet.ini)')
    parser.add_argument('--worksheet', default=None, help='worksheet name (default: first sheet)')
    parser.add_argument('--dryRun', action='store_true', help='show what would be written, change nothing')
    args = parser.parse_args(argv)
    settings = getSettings(args.sheetId, args.worksheet)
    worksheet = openWorksheet(settings)
    print('\n\tPushing to "%s"...' % (worksheet.title))
    for case in args.cases:
        pushSummary(worksheet, postInternal.readSummary(case), args.comment, args.dryRun)


if __name__ == '__main__':
    main()
