import configparser
import os
import re
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

SHEET_ID = "1yQ7-aeIgvioDSCrVY6mq3w_ZqJRzT0cZpypf8M3vYhU"
WORKSHEET_NAME = "250001 - Trials List"

def get_worksheet(sheet_id, worksheet_name):
    if not sheet_id or sheet_id == "PASTE_YOUR_GOOGLE_SHEET_ID":
        raise ValueError("Please set SHEET_ID in the script before running.")

    creds_path = get_credentials_path()
    scopes = ["https://www.googleapis.com/auth/spreadsheets"]
    creds = Credentials.from_service_account_file(creds_path, scopes=scopes)
    client = gspread.authorize(creds)

    spreadsheet = client.open_by_key(sheet_id)
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

TARGET_COLUMNS = ["B","F", "G", "H", "N", "O", "P", "Q", "R", "S", "T", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK", "AL", "AM", "AN", "AO"]
PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_CREDENTIALS_FILE = PROJECT_DIR / "credentials.json"

def get_or_create_trial_row(sheet_id, worksheet_name, trial_number):
    """
    Find trial_number in column A and return its row.
    If missing, write it to the next available row in column A and return that row.
    """
    trial_number = str(trial_number).strip()
    if trial_number == "":
        raise ValueError("trial_number is empty.")

    worksheet = get_worksheet(sheet_id, worksheet_name)
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

    worksheet = get_worksheet(sheet_id, worksheet_name)

    for col, value in zip(target_columns, values):
        cell = f"{col.upper()}{target_row}"
        worksheet.update(
            range_name=cell,
            values=[[str(value)]],
            value_input_option="RAW",
        )


def main():
    case_path = os.getcwd()
    case_name = os.path.basename(os.path.normpath(case_path))
    path = os.path.split(case_path)[0]

    case_setup_path = Path(case_path) / "caseSetup"
    if not case_setup_path.exists():
        raise FileNotFoundError(f"caseSetup file not found at: {case_setup_path}")

    full_case_setup_dict = configparser.ConfigParser()
    full_case_setup_dict.optionxform = str
    full_case_setup_dict.read_file(open(case_setup_path))

    coeff_files = getCoeffPaths(case_path, case_name)
    avg_data = averageCoeffs(full_case_setup_dict, case_name, "all", coeff_files)
    num_cells, mesher, _ = cellCount(full_case_setup_dict, case_path, case_name)
    inlet_mag, _, yaw, _, _, _, _ = bcParser(full_case_setup_dict, path, case_name)
    run_date, run_time, _, _ = getOfVersion(case_path)

    fw_cd = ""
    fw_cl = ""
    rw_cd = ""
    rw_cl = ""
    for part in coeff_files:
        part_lower = part.lower()
        if part_lower.startswith("fw") and fw_cd == "":
            fw_data = averageCoeffs(full_case_setup_dict, case_name, part, coeff_files)
            fw_cd = fw_data["cd"]
            fw_cl = fw_data["cl"]
        elif part_lower.startswith("rw") and rw_cd == "":
            rw_data = averageCoeffs(full_case_setup_dict, case_name, part, coeff_files)
            rw_cd = rw_data["cd"]
            rw_cl = rw_data["cl"]
        if fw_cd != "" and rw_cd != "":
            break

    data_to_write = [
        run_date,
        run_time,
        num_cells,
        mesher,
        inlet_mag,
        yaw,
        "1.225",
        "0",
        "0",
        "1",
        "1.55",
        avg_data["cd"],
        avg_data["cl"],
        avg_data["clf"],
        avg_data["clr"],
        avg_data["csf"],
        avg_data["csr"],
        avg_data["cd_ci"],
        avg_data["cl_ci"],
        fw_cd,
        fw_cl,
        rw_cd,
        rw_cl,
    ]

    trial_number = get_trial_number()
    target_row = get_or_create_trial_row(SHEET_ID, WORKSHEET_NAME, trial_number)

    write_to_sheet_cells(
        SHEET_ID,
        data_to_write,
        WORKSHEET_NAME,
        target_row,
        TARGET_COLUMNS,
    )
    print(
        f"Done. Trial {trial_number} mapped to row {target_row}. "
        f"Wrote {len(data_to_write)} value(s) in columns: "
        f"{', '.join(col.upper() for col in TARGET_COLUMNS)}."
    )

if __name__ == "__main__":
    main()