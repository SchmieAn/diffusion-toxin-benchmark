"""
Analyze structural similarity results from a ProteinMPNN pipeline run.
The same script was also used for the ESM Inverse Folding pipeline runs.
Therefore solely the path has to change to the ESM Inverse Folding directory.

Reads a TM-score CSV file and reports the fraction of entries meeting a
user-defined TM-score threshold, both overall and filtered to one selected
reference toxin (REFERENCE_FILE, e.g. 1KBA or 4K09). Reference toxins are only
available from run 6 onwards; for earlier runs only the overall evaluation is
reported.

Usage:
    python calc_structure_similarity.py
    -> prompts for a TM-score threshold (float)

Expected CSV columns:
    tm1             (float, 0–1)
    reference_file  (str, e.g. 'pdb1kba.ent'; only required for runs >= 6)

Note:
    This script was used for both the ProteinMPNN and the ESM Inverse Folding
    pipeline runs. Only the path(s) and the reference toxin in the Configuration
    section below need to be changed to point at the respective pipeline's output
    directory and the reference toxin of the analyzed run.
"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Directory containing the TM-score CSV file to analyze.
BASE_DIR = "/path_to_your_base_directory/Pipeline_ESM"

# Run number of the pipeline run to analyze.
RUN_NUMBER = 15

# Name of the input CSV file — depends on the run number.
FILENAME_TEMPLATE = "TM_score_{run_number}_filtered.csv"

# Reference toxin (value in the 'reference_file' column) for the toxin-specific
# evaluation. Runs 6-10 used 4K09 ("pdb4k09.ent"), runs 11-15 used 1KBA ("pdb1kba.ent").
# Reference toxins are only available from run 6 onwards; for earlier runs this
# value is ignored and the toxin-specific evaluation is skipped.
REFERENCE_FILE = "pdb1kba.ent"

# ---------------------------------------------------------------------------


import os
import pandas as pd


def analyze_run_data(run_number, base_dir, tm_threshold, reference_file):
    file_path = os.path.join(base_dir, FILENAME_TEMPLATE.format(run_number=run_number))
    print(f"Reading file: {file_path}")

    if not os.path.exists(file_path):
        print(f"Error: file '{file_path}' not found.")
        return

    try:
        df = pd.read_csv(file_path)

        if 'tm1' not in df.columns:
            print("Error: Column 'tm1' does not exist.")
            print(f"Found columns: {list(df.columns)}")
            return

        total_rows = len(df)
        if total_rows == 0:
            print("The CSV file is empty.")
            return

        # Overall threshold filter
        filtered_rows = len(df[df['tm1'] >= tm_threshold])
        percentage = (filtered_rows / total_rows) * 100

        print("\n" + "=" * 40)
        print(f"Results for run {run_number}")
        print("=" * 40)
        print(f"Total rows:               {total_rows}")
        print(f"Filtered rows:            {filtered_rows}  (tm1 >= {tm_threshold})")
        print(f"Percentage:               {percentage:.3f}%")

        # Filter to the selected reference toxin (only available from run 6 onwards)
        if run_number >= 6:
            tox_df = df[df['reference_file'] == reference_file]
            tox_filtered_rows = len(tox_df)
            if tox_filtered_rows == 0:
                print(f"No rows found for reference toxin '{reference_file}'.")
            else:
                tox_thresh_filtered = len(tox_df[tox_df['tm1'] >= tm_threshold])
                percentage_tox_only = (tox_thresh_filtered / tox_filtered_rows) * 100
                print(f"Reference toxin:          {reference_file}")
                print(f"Toxin reference rows:     {tox_filtered_rows}")
                print(f"Tox rows over threshold:  {tox_thresh_filtered}")
                print(f"Percentage tox rows:      {percentage_tox_only:.3f}%")
        else:
            print("No reference toxin evaluation (reference toxins exist from run 6 onwards).")
        print("=" * 40)

    except Exception as e:
        print(f"An error occurred while processing the file: {e}")


if __name__ == "__main__":
    try:
        tm_threshold = float(input("Threshold for TM-score: "))
        analyze_run_data(RUN_NUMBER, BASE_DIR, tm_threshold, REFERENCE_FILE)
    except ValueError:
        print("Error: Please enter an integer or float number.")
