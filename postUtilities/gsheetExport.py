import configparser
import os
import re
import sys
import subprocess
from pathlib import Path
from summary import *
import gspread
from google.oauth2.service_account import Credentials

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

    return spreadsheet.worksheet(worksheet_name) if worksheet_name else spreadsheet.sheet1

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


def build_sheet_values_from_summary(summary_dict, case_setup_dict):
    symmetry = str(summary_dict.get("Symmetry", "")).strip().lower()
    is_half = "TRUE" if symmetry == "half" else "FALSE"

    density = 1.225
    wheel_base = ""
    if case_setup_dict is not None:
        try:
            density = parse_float(case_setup_dict['GLOBAL_MATERIAL']['DENSITY'], 1.225)
        except Exception:
            density = 1.225
        try:
            wheel_base = parse_float(case_setup_dict['BC_SETUP']['REFLEN'], "")
        except Exception:
            wheel_base = ""

    fw_cd = get_part_coeff(summary_dict, "fw", "CD")
    fw_cl = get_part_coeff(summary_dict, "fw", "CL")
    rw_cd = get_part_coeff(summary_dict, "rw", "CD")
    rw_cl = get_part_coeff(summary_dict, "rw", "CL")

    return [
        summary_dict.get("Run Date", ""),
        summary_dict.get("Solve Time", ""),
        summary_dict.get("Num. Cells", ""),
        summary_dict.get("Mesher", ""),
        is_half,
        summary_dict.get("Velocity", ""),
        summary_dict.get("Yaw", ""),
        density,
        "0",
        "0",
        summary_dict.get("Ref. Area (m^2)", ""),
        wheel_base,
        summary_dict.get("Cd", ""),
        summary_dict.get("Cl", ""),
        summary_dict.get("Cl(f)", ""),
        summary_dict.get("Cl(r)", ""),
        summary_dict.get("Cs(f)", ""),
        summary_dict.get("Cs(r)", ""),
        summary_dict.get("Cd CI", ""),
        summary_dict.get("Cl CI", ""),
        fw_cd,
        fw_cl,
        rw_cd,
        rw_cl,
    ]


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
    data_to_write = build_sheet_values_from_summary(summary_dict, case_setup_dict)
    target_row = get_or_create_trial_row(sheet_id, worksheet_name, case_name, job_name)
    write_to_sheet_cells(
        sheet_id,
        data_to_write,
        worksheet_name,
        target_row,
        TARGET_COLUMNS,
        jobName=job_name,
    )
    print(
        f"Done. Trial {case_name} mapped to row {target_row}. "
        f"Wrote {len(data_to_write)} value(s) in columns: "
        f"{', '.join(col.upper() for col in TARGET_COLUMNS)}."
    )

TARGET_COLUMNS = ["B","F", "G", "H", "M","N", "O", "P", "Q", "R", "S", "T", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK", "AL", "AM", "AN", "AO"]
PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_CREDENTIALS_FILE = PROJECT_DIR / "credentials.json"

def get_or_create_trial_row(sheet_id, worksheet_name, trial_number,jobNumber):
    """
    Find trial_number in column A and return its row.
    If missing, write it to the next available row in column A and return that row.
    """
    trial_number = str(trial_number).strip()
    if trial_number == "":
        raise ValueError("trial_number is empty.")

    worksheet = get_worksheet(sheet_id, worksheet_name,jobNumber)
    column_a = worksheet.col_values(1)

    for row_index, value in enumerate(column_a, start=1):
        if str(value).strip() == trial_number:
            return row_index

    target_row = len(column_a) + 1
    worksheet.update(
        range_name=f"A{target_row}",
        values=[[trial_number]],
        value_input_option="RAW",
    )
    return target_row


def write_to_sheet_cells(
    sheet_id,
    values,
    worksheet_name,
    target_row,
    target_columns,
    jobName
):
    if not sheet_id or sheet_id == "PASTE_YOUR_GOOGLE_SHEET_ID":
        raise ValueError("Please set SHEET_ID in the script before running.")
    if not values:
        raise ValueError("DATA_TO_WRITE is empty.")
    if target_row < 1:
        raise ValueError("target_row must be >= 1.")
    if len(values) != len(target_columns):
        raise ValueError("DATA_TO_WRITE and TARGET_COLUMNS must have the same length.")
    if not all(col and col.isalpha() for col in target_columns):
        raise ValueError("TARGET_COLUMNS must contain only column letters (e.g. A, B, AA).")

    worksheet = get_worksheet(sheet_id, worksheet_name,jobName)

    for col, value in zip(target_columns, values):
        cell = f"{col.upper()}{target_row}"
        worksheet.update(
            range_name=cell,
            values=[[str(value)]],
            value_input_option="RAW",
        )


def main():
    global jobName
    case_path = Path(os.getcwd()).resolve()
    case_name = case_path.name
    #checking if run in a trial directory
    # if os.path.basename(os.path.split(case_path)[0]) != 'CASES':
    #     sys.exit('ERROR! Please run in a trial directory!')
    
    if case_path.parent.name == 'CASES':
        jobPath = case_path.parent.parent
    elif case_path.parent.parent.name == 'CASES':
        jobPath = case_path.parent.parent.parent
    elif '_' in case_name:
        jobPath = Path(os.path.abspath(os.path.join(os.getcwd(), "../../../")))
    else:
        jobPath = Path(os.path.abspath(os.path.join(os.getcwd(), "../../")))
    jobName = jobPath.name
    #jobName = getWorksheetName(case_path)
    WORKSHEET_NAME = '%s - Trials List' % (jobName)
    #jobRoot = os.path.abspath(os.path.join(path, os.pardir))
    
    gsheetID = getSheetId(os.path.join(jobPath, '02_reference', 'GSheet', '%s.gsheet' % jobName), jobName=jobName)
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
