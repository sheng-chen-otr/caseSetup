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

import argparse
import configparser
import os
import re
import sys
import subprocess
from pathlib import Path
from summary import *
import gspread
from google.oauth2.service_account import Credentials

# sheet layout. columns are found by name, not letter: HEADER_ROW holds the column names
# (Trial, Date Run, Description, ...) and GROUP_ROW the merged group labels above them
# (Boundary Conditions, RAW, FW RAW, ...). a column is written as 'Header', or as
# ('Group', 'Header') when the header name is used more than once in the sheet
HEADER_ROW = 7
GROUP_ROW = HEADER_ROW - 1
TRIAL_COLUMN = 'Trial'

PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_CREDENTIALS_FILE = PROJECT_DIR / "credentials.json"

_WORKSHEET_CACHE = {}
_HEADER_CACHE = {}


def get_credentials_path():
    env_path = os.getenv("GOOGLE_APPLICATION_CREDENTIALS")
    if env_path:
        p = Path(env_path).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(
                f"GOOGLE_APPLICATION_CREDENTIALS points to missing file: {p}"
            )
        return str(p)

    if DEFAULT_CREDENTIALS_FILE.exists():
        return str(DEFAULT_CREDENTIALS_FILE)

    raise FileNotFoundError(
        "No credentials found. Set GOOGLE_APPLICATION_CREDENTIALS or add "
        "credentials.json in this project folder."
    )


def getSheetId(gsheetConfigPath,jobName):
    gsheetConfig = configparser.ConfigParser()
    gsheetConfig.optionxform = str
    try:
        print('Importing gsheet config...')
        gsheetConfig.read_file(open(gsheetConfigPath))
    except Exception as e:
        print('ERROR! Unable to import %s' % (gsheetConfigPath))
        sys.exit(e)

    return gsheetConfig[jobName]['sheetID']

def getWorksheetName(path):
    # Get the absolute path of the current script
    current_path = os.path.abspath(path)

    # Move two levels up
    two_up_path = os.path.abspath(os.path.join(current_path, os.pardir, os.pardir))

    # Return only the directory name (last folder in the path)
    return os.path.basename(two_up_path)


def get_worksheet(sheet_id, worksheet_name,jobnumber):
    # opened once per run and reused, every open costs Google API requests
    key = (sheet_id, worksheet_name)
    if key in _WORKSHEET_CACHE:
        return _WORKSHEET_CACHE[key]

    if not sheet_id or sheet_id == "PASTE_YOUR_GOOGLE_SHEET_ID":
        raise ValueError("Please set SHEET_ID in the script before running.")

    creds_path = get_credentials_path()
    scopes = ["https://www.googleapis.com/auth/spreadsheets"]
    creds = Credentials.from_service_account_file(creds_path, scopes=scopes)
    client = gspread.authorize(creds)
    spreadsheet = client.open_by_key(sheet_id)
    spreadsheetName = spreadsheet.title
    if not jobnumber in spreadsheetName:
        sys.exit('ERROR! Google sheets title does not match job code, please check!')

    worksheet = spreadsheet.worksheet(worksheet_name) if worksheet_name else spreadsheet.sheet1
    _WORKSHEET_CACHE[key] = worksheet
    return worksheet

def get_trial_number(case_name=None, case_path=None):
    """
    Return trial number as a string (keeps leading zeros).
    Examples handled: '003', 'trial003', 'trial_003', 'case_003_half'.
    """
    if case_name is None:
        base_path = case_path if case_path else os.getcwd()
        case_name = os.path.basename(os.path.normpath(base_path))

    case_name = str(case_name).strip()
    if case_name == "":
        raise ValueError("Case name is empty; cannot determine trial number.")

    # Case name itself is only digits (e.g. '003')
    if case_name.isdigit():
        return case_name

    # Common pattern: trial003 / trial_003 / trial-003
    trial_match = re.search(r"trial[_-]?(\d+)", case_name, flags=re.IGNORECASE)
    if trial_match:
        return trial_match.group(1)

    # Fallback: first number group in the case name
    number_match = re.search(r"(\d+)", case_name)
    if number_match:
        return number_match.group(1)

    raise ValueError(f"No numeric trial number found in case name: {case_name}")


def read_summary_csv(case_path):
    """Read summary.csv (key,value rows) and return a dictionary."""
    summary_path = Path(case_path) / "summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(f"summary.csv not found at: {summary_path}")
    summary_df = pd.read_csv(summary_path, header=None)
    if summary_df.shape[1] < 2:
        raise ValueError(f"summary.csv is invalid at: {summary_path}")
    keys = summary_df.iloc[:, 0].astype(str).str.strip()
    vals = summary_df.iloc[:, 1]
    return dict(zip(keys, vals))


def load_case_setup(case_path):
    case_setup_path = Path(case_path) / "caseSetup"
    if not case_setup_path.exists():
        return None
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read_file(open(case_setup_path))
    return cfg


def parse_float(value, default=""):
    try:
        return float(value)
    except Exception:
        return default


def get_part_coeff(summary_dict, prefix, coeff_key):
    suffix = f" {coeff_key.upper()}"
    for key, value in summary_dict.items():
        if key.lower().startswith(prefix.lower()) and key.endswith(suffix):
            return value
    return ""


def _first(value, default=""):
    """Return the first element of a list/tuple (as produced by caseSetup.py's
    in-memory caseSetup dict), or the value itself (as returned by configparser
    when reading a caseSetup file back off disk)."""
    if isinstance(value, (list, tuple)):
        return value[0] if len(value) else default
    if value is None:
        return default
    return value


def _full_text(value):
    """Whole text of a caseSetup value. The in-memory dict splits values into words, a
    caseSetup file read with configparser gives the plain string."""
    if isinstance(value, (list, tuple)):
        return ' '.join(str(v) for v in value).strip()
    if value is None:
        return ''
    return str(value).strip()


def _setup_value(case_setup_dict, section, key):
    try:
        return case_setup_dict[section][key]
    except Exception:
        return None


# ---------------------------------------------------------------------- values to write

def build_sheet_values_from_summary(summary_dict, case_setup_dict):
    """Result columns from summary.csv, as {column: value}."""
    symmetry = str(summary_dict.get("Symmetry", "")).strip().lower()
    is_half = "TRUE" if symmetry == "half" else "FALSE"

    density = 1.225
    wheel_base = ""
    if case_setup_dict is not None:
        density = parse_float(_first(_setup_value(case_setup_dict, 'GLOBAL_MATERIAL', 'DENSITY')), 1.225)
        wheel_base = parse_float(_first(_setup_value(case_setup_dict, 'BC_SETUP', 'REFLEN')), "")

    return {
        'Date Run': summary_dict.get("Run Date", ""),
        'Solve Time': summary_dict.get("Solve Time", ""),
        'Cell Count': summary_dict.get("Num. Cells", ""),
        'Mesher': summary_dict.get("Mesher", ""),
        'Is Half?': is_half,
        'Speed (m/s^2)': summary_dict.get("Velocity", ""),
        'Yaw': summary_dict.get("Yaw", ""),
        'Density': density,
        'Frontal Area': summary_dict.get("Ref. Area (m^2)", ""),
        'Wheel Base': wheel_base,
        ('RAW', 'CD'): summary_dict.get("Cd", ""),
        ('RAW', 'CL'): summary_dict.get("Cl", ""),
        ('RAW', 'CL(F)'): summary_dict.get("Cl(f)", ""),
        ('RAW', 'CL(R)'): summary_dict.get("Cl(r)", ""),
        ('RAW', 'CS(F)'): summary_dict.get("Cs(f)", ""),
        ('RAW', 'CS(R)'): summary_dict.get("Cs(r)", ""),
        ('RAW', 'CI-CD'): summary_dict.get("Cd CI", ""),
        ('RAW', 'CI-CL'): summary_dict.get("Cl CI", ""),
        ('FW RAW', 'CD'): get_part_coeff(summary_dict, "fw", "CD"),
        ('FW RAW', 'CL'): get_part_coeff(summary_dict, "fw", "CL"),
        ('RW RAW', 'CD'): get_part_coeff(summary_dict, "rw", "CD"),
        ('RW RAW', 'CL'): get_part_coeff(summary_dict, "rw", "CL"),
    }


def build_setup_time_values(case_setup_dict):
    """
    Setup-time-known columns (symmetry, velocity, yaw, density, reference area,
    wheelbase) from a case's caseSetup dict, as {column: value}. Result columns (run
    date, solve time, cell count, Cd, Cl, etc.) are intentionally left out since
    they are not known until the case has actually run.
    """
    sim_sym = str(_first(_setup_value(case_setup_dict, 'GLOBAL_SIM_CONTROL', 'SIM_SYM'), '')).strip().lower()
    is_half = "TRUE" if sim_sym == 'half' else "FALSE"

    velocity = _first(_setup_value(case_setup_dict, 'BC_SETUP', 'INLET_MAG'))
    yaw = 0 if sim_sym == 'half' else _first(_setup_value(case_setup_dict, 'BC_SETUP', 'YAW'))
    density = parse_float(_first(_setup_value(case_setup_dict, 'GLOBAL_MATERIAL', 'DENSITY')), 1.225)
    ref_area = _first(_setup_value(case_setup_dict, 'BC_SETUP', 'REFAREA'))
    wheel_base = parse_float(_first(_setup_value(case_setup_dict, 'BC_SETUP', 'REFLEN')), "")

    return {
        'Is Half?': is_half,
        'Speed (m/s^2)': velocity,
        'Yaw': yaw,
        'Density': density,
        'Frontal Area': ref_area,
        'Wheel Base': wheel_base,
    }


def build_title_values(case_setup_dict, point=None):
    """Description = CASENAME (plus ' - point N' for a ride-height child), Comments =
    CASEDESCRP. A blank value is left out so the cell keeps what is already there."""
    values = {}
    case_name = _full_text(_setup_value(case_setup_dict, 'TITLES', 'CASENAME'))
    case_descrp = _full_text(_setup_value(case_setup_dict, 'TITLES', 'CASEDESCRP'))
    if case_name:
        values['Description'] = case_name if point is None else '%s - point %s' % (case_name, point)
    if case_descrp:
        values['Comments'] = case_descrp
    return values


def _ride_height_unit_scale(case_setup_dict):
    # ride-height values -> mm
    unit = str(_first(_setup_value(case_setup_dict, 'RIDE_HEIGHT_SETUP', 'RH_UNIT'), 'mm')).strip().lower()
    return 1000.0 if unit == 'm' else 1.0


def _init_rh_mm(case_setup_dict):
    """INIT_RH (FL FR RL RR) in mm, or None when it is blank or not 4 numbers."""
    raw = _setup_value(case_setup_dict, 'RIDE_HEIGHT_SETUP', 'INIT_RH')
    tokens = _full_text(raw).split()
    if len(tokens) != 4:
        return None
    try:
        scale = _ride_height_unit_scale(case_setup_dict)
        return [float(t) * scale for t in tokens]
    except ValueError:
        return None


def build_ride_height_values(case_setup_dict, case_name, parent_case_path=None):
    """
    Ride Height Front/Rear (mm) and Roll (deg).
      ride-height child: from the parent's rideHeights_updated.csv row for this child.
        Front = (fl + fr)/2, Rear = (rl + rr)/2. These are target ride heights when
        INIT_RH is set, otherwise the change from the supplied geometry
      any other case: the INIT_RH averages when set, otherwise 0. Roll 0
    """
    if parent_case_path is not None:
        rh_csv = Path(parent_case_path) / 'rideHeights_updated.csv'
        try:
            rh_map = pd.read_csv(rh_csv)
            match = rh_map[rh_map['caseName'].astype(str) == str(case_name)]
            if len(match) == 0:
                raise ValueError('no row for %s' % (case_name))
            row = match.iloc[0]
            # fl/fr/rl/rr in rideHeights_updated.csv are already in m
            front = (float(row['fl']) + float(row['fr'])) / 2.0 * 1000.0
            rear = (float(row['rl']) + float(row['rr'])) / 2.0 * 1000.0
            roll = float(row['roll']) if 'roll' in rh_map.columns else 0.0
            return {
                ('Ride Height', 'Front'): round(front, 3),
                ('Ride Height', 'Rear'): round(rear, 3),
                'Roll': round(roll, 3),
            }
        except Exception as e:
            print('\t\tWARNING! Unable to read ride heights for %s from %s: %s' % (case_name, rh_csv, e))
            return {}

    init_rh = _init_rh_mm(case_setup_dict) if case_setup_dict is not None else None
    if init_rh is None:
        front = rear = 0
    else:
        front = round((init_rh[0] + init_rh[1]) / 2.0, 3)
        rear = round((init_rh[2] + init_rh[3]) / 2.0, 3)
    return {
        ('Ride Height', 'Front'): front,
        ('Ride Height', 'Rear'): rear,
        'Roll': 0,
    }


def build_case_identity_values(case_setup_dict, case_path):
    """Description, Comments, ride heights and roll for a case, working out whether it is
    a ride-height child (parentName_N) from its folder."""
    case_path = Path(case_path).resolve()
    parent_case_path = detect_parent_case_for_child(case_path)
    point = None
    if parent_case_path is not None:
        point = int(re.match(r'^.*_(\d+)$', case_path.name).group(1))
    values = {}
    if case_setup_dict is not None:
        values.update(build_title_values(case_setup_dict, point))
    values.update(build_ride_height_values(case_setup_dict, case_path.name, parent_case_path))
    return values


# ---------------------------------------------------------------------- sheet access

def _norm(text):
    return ' '.join(str(text).split()).lower()


def _column_letter(index):
    """0-based column index -> sheet letter (0 -> A, 26 -> AA)."""
    letters = ''
    n = index + 1
    while n:
        n, r = divmod(n - 1, 26)
        letters = chr(65 + r) + letters
    return letters


def _column_label(column):
    return '%s > %s' % column if isinstance(column, tuple) else column


def read_sheet_headers(worksheet):
    """[(group, header), ...] per column, normalised. A merged group label is only stored in
    its first cell, so each label is carried right until the next one."""
    key = worksheet.id
    if key in _HEADER_CACHE:
        return _HEADER_CACHE[key]
    groups = worksheet.row_values(GROUP_ROW)
    headers = worksheet.row_values(HEADER_ROW)
    if not any(str(h).strip() for h in headers):
        raise ValueError('header row %d of "%s" is empty. Check HEADER_ROW in gsheetExport.py.'
                         % (HEADER_ROW, worksheet.title))
    out = []
    group = ''
    for i in range(len(headers)):
        if i < len(groups) and str(groups[i]).strip():
            group = _norm(groups[i])
        out.append((group, _norm(headers[i])))
    _HEADER_CACHE[key] = out
    return out


def find_column(worksheet, column):
    """Sheet letter for a column given as 'Header' or ('Group', 'Header'). Returns None when
    it is not in the sheet. Raises ValueError when it matches more than one column."""
    headers = read_sheet_headers(worksheet)
    if isinstance(column, tuple):
        group, name = _norm(column[0]), _norm(column[1])
        matches = [i for i, (g, h) in enumerate(headers) if g == group and h == name]
    else:
        name = _norm(column)
        matches = [i for i, (_g, h) in enumerate(headers) if h == name]
    if len(matches) == 0:
        return None
    if len(matches) > 1:
        letters = ', '.join(_column_letter(i) for i in matches)
        fix = ('fix the group labels in row %d' % GROUP_ROW if isinstance(column, tuple) else
               'give it as (group, header) in gsheetExport.py')
        raise ValueError('column "%s" matches more than one column in rows %d-%d (%s), %s.'
                         % (_column_label(column), GROUP_ROW, HEADER_ROW, letters, fix))
    return _column_letter(matches[0])


def get_or_create_trial_row(worksheet, trial_name):
    """
    Find trial_name in the Trial column and return its row.
    If missing, write it to the next available row and return that row.
    """
    trial_name = str(trial_name).strip()
    if trial_name == "":
        raise ValueError("trial name is empty.")

    letter = find_column(worksheet, TRIAL_COLUMN)
    if letter is None:
        raise ValueError('no "%s" column in row %d of "%s". Check HEADER_ROW in gsheetExport.py.'
                         % (TRIAL_COLUMN, HEADER_ROW, worksheet.title))
    col_index = gspread.utils.a1_to_rowcol('%s1' % letter)[1]
    column_values = worksheet.col_values(col_index)

    for row_index, value in enumerate(column_values, start=1):
        if row_index > HEADER_ROW and str(value).strip() == trial_name:
            return row_index

    target_row = max(len(column_values), HEADER_ROW) + 1
    worksheet.update(
        range_name=f"{letter}{target_row}",
        values=[[trial_name]],
        value_input_option="RAW",
    )
    return target_row


def write_row(worksheet, target_row, values_by_column):
    """Write {column: value} into target_row in one request. Columns not found in the sheet
    are skipped with a warning. Returns the list of columns written."""
    if target_row <= HEADER_ROW:
        raise ValueError("target_row must be below the header row (%d)." % (HEADER_ROW))

    updates, written, missing = [], [], []
    for column, value in values_by_column.items():
        letter = find_column(worksheet, column)
        if letter is None:
            missing.append(_column_label(column))
            continue
        updates.append({'range': f"{letter}{target_row}", 'values': [[str(value)]]})
        written.append('%s (%s)' % (_column_label(column), letter))

    if missing:
        print('\t\tWARNING! These columns were not found in rows %d-%d of "%s" and were skipped: %s'
              % (GROUP_ROW, HEADER_ROW, worksheet.title, ', '.join(missing)))
    if updates:
        worksheet.batch_update(updates, value_input_option="RAW")
    return written


# ---------------------------------------------------------------------- prep (before a run)

def resolve_job_context(case_path):
    """Resolve (jobPath, jobName, worksheetName, sheetId) for a given case path,
    using the same job-root inference rules as main()."""
    case_path = Path(case_path).resolve()
    case_name = case_path.name

    if case_path.parent.name == 'CASES':
        jobPath = case_path.parent.parent
    elif case_path.parent.parent.name == 'CASES':
        jobPath = case_path.parent.parent.parent
    elif '_' in case_name:
        jobPath = case_path.parents[2]
    else:
        jobPath = case_path.parents[1]

    jobName = jobPath.name
    worksheet_name = '%s - Trials List' % (jobName)
    gsheet_id = getSheetId(
        os.path.join(jobPath, '02_reference', 'GSheet', '%s.gsheet' % jobName),
        jobName=jobName,
    )
    return jobPath, jobName, worksheet_name, gsheet_id


def prep_single_case_row(sheet_id, worksheet_name, case_path, job_name, case_setup_dict):
    """Ensure a row exists for the case (reserving the next row if missing) and
    fill in only the setup-time-known columns, leaving result columns untouched."""
    case_path = Path(case_path).resolve()
    case_name = case_path.name
    if case_setup_dict is None:
        print(f"\t\tWARNING! Unable to read caseSetup for {case_name}; skipping gsheet row prep.")
        return
    worksheet = get_worksheet(sheet_id, worksheet_name, job_name)
    target_row = get_or_create_trial_row(worksheet, case_name)
    values = build_setup_time_values(case_setup_dict)
    values.update(build_case_identity_values(case_setup_dict, case_path))
    write_row(worksheet, target_row, values)
    print(f"\t\tPrepared gsheet row {target_row} for {case_name}.")


def find_ride_height_children(parent_case_path):
    """Return [(pointNumber, childPath), ...] sorted by ascending point number for
    <parentName>_<n> child directories directly under the parent case directory."""
    parent_case_path = Path(parent_case_path)
    parent_name = parent_case_path.name
    pattern = re.compile(r'^%s_(\d+)$' % re.escape(parent_name))
    children = []
    for child in parent_case_path.iterdir():
        if not child.is_dir():
            continue
        m = pattern.match(child.name)
        if m:
            children.append((int(m.group(1)), child))
    children.sort(key=lambda item: item[0])
    return children


def prep_ride_height_rows(case_path, fullCaseSetupDict, jobName, worksheet_name, gsheet_id):
    """Ensure the parent row exists first, then the child rows in ascending
    ride-height-point order, for a ride-height mapping parent case."""
    case_path = Path(case_path).resolve()

    #parent row first
    prep_single_case_row(gsheet_id, worksheet_name, case_path, jobName, fullCaseSetupDict)

    #then children, in ascending ride-height-point order
    children = find_ride_height_children(case_path)
    if not children:
        print("\t\tWARNING! RUN_RIDE_HEIGHT is True but no child cases were found under: %s" % (case_path))
        return
    for pointNumber, childPath in children:
        childCfg = load_case_setup(childPath)
        prep_single_case_row(gsheet_id, worksheet_name, childPath, jobName, childCfg)


def prep_gsheet_rows(case_path, fullCaseSetupDict=None):
    """
    Ensure this case's gsheet row exists (and, for a ride-height mapping parent,
    every child's row as well, in parent-then-ascending-point-number order),
    pre-filled with the setup-time values already known from caseSetup. Intended
    to run right after caseSetup, before a case is submitted. Any failure (missing
    credentials, sheet not found, network error, etc.) is caught and printed as a
    warning only, so callers (e.g. caseSetup.py) can continue regardless.
    """
    case_path = Path(case_path).resolve()
    try:
        jobPath, jobName, worksheet_name, gsheet_id = resolve_job_context(case_path)
        print('\tFound gsheet id: %s' % (gsheet_id))

        if fullCaseSetupDict is None:
            fullCaseSetupDict = load_case_setup(case_path)

        is_ride_height_parent = False
        if fullCaseSetupDict is not None and 'RIDE_HEIGHT_SETUP' in fullCaseSetupDict:
            try:
                is_ride_height_parent = str(_first(
                    fullCaseSetupDict['RIDE_HEIGHT_SETUP']['RUN_RIDE_HEIGHT'], ''
                )).strip().lower() == 'true'
            except Exception:
                is_ride_height_parent = False

        if is_ride_height_parent:
            prep_ride_height_rows(case_path, fullCaseSetupDict, jobName, worksheet_name, gsheet_id)
        else:
            prep_single_case_row(gsheet_id, worksheet_name, case_path, jobName, fullCaseSetupDict)
    except Exception as e:
        print('\tWARNING! Unable to prepare gsheet row(s) for %s: %s' % (case_path.name, e))


# ---------------------------------------------------------------------- push (after a run)

def detect_parent_case_for_child(case_path):
    """Return parent case path if case is likely parentName_# under parentName folder."""
    case_path = Path(case_path).resolve()
    case_name = case_path.name
    m = re.match(r'^(.*)_\d+$', case_name)
    if not m:
        return None
    parent_name = m.group(1)

    immediate_parent = case_path.parent
    if immediate_parent.name == parent_name:
        return immediate_parent

    sibling_parent = case_path.parent / parent_name
    if sibling_parent.exists() and sibling_parent.is_dir():
        return sibling_parent

    return None


def regenerate_case_summary(case_path):
    post_run_script = Path(__file__).resolve().parent / 'postRun.py'
    if not post_run_script.exists():
        raise FileNotFoundError(f"postRun.py not found at: {post_run_script}")
    subprocess.run(
        [sys.executable, str(post_run_script), '--summary'],
        cwd=str(case_path),
        check=True,
    )


def push_case_summary(sheet_id, worksheet_name, case_path, job_name):
    case_path = Path(case_path).resolve()
    case_name = case_path.name
    summary_dict = read_summary_csv(case_path)
    case_setup_dict = load_case_setup(case_path)
    values = build_sheet_values_from_summary(summary_dict, case_setup_dict)
    values.update(build_case_identity_values(case_setup_dict, case_path))

    worksheet = get_worksheet(sheet_id, worksheet_name, job_name)
    target_row = get_or_create_trial_row(worksheet, case_name)
    written = write_row(worksheet, target_row, values)
    print(
        f"Done. Trial {case_name} mapped to row {target_row}. "
        f"Wrote {len(written)} value(s): {', '.join(written)}."
    )


def main():
    global jobName
    parser = argparse.ArgumentParser(
        description='Push case summary results to a Google Sheet, or reserve/pre-fill row(s) before a case runs.'
    )
    parser.add_argument(
        '--prep', action='store_true',
        help='Reserve gsheet row(s) for this case (and its ride-height children, if any, in parent-then-child '
             'order) and fill in the setup-time values already known from caseSetup, without touching result '
             'columns. Intended to run right after caseSetup, before the case is submitted.'
    )
    args = parser.parse_args()

    case_path = Path(os.getcwd()).resolve()
    case_name = case_path.name
    #checking if run in a trial directory
    # if os.path.basename(os.path.split(case_path)[0]) != 'CASES':
    #     sys.exit('ERROR! Please run in a trial directory!')

    if args.prep:
        prep_gsheet_rows(case_path)
        return

    jobPath, jobName, WORKSHEET_NAME, gsheetID = resolve_job_context(case_path)
    print('\tFound gsheet id: %s' % (gsheetID))

    parent_case_path = detect_parent_case_for_child(case_path)
    if parent_case_path is not None:
        print(f"Detected ride-height child case. Regenerating parent summary at: {parent_case_path}")
        regenerate_case_summary(parent_case_path)

    push_case_summary(gsheetID, WORKSHEET_NAME, case_path, jobName)

    if parent_case_path is not None:
        push_case_summary(gsheetID, WORKSHEET_NAME, parent_case_path, jobName)

if __name__ == "__main__":
    main()
