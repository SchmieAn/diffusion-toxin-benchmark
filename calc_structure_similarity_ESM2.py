"""
Analyze structural similarity results from an ESM2 pipeline run.

Reads the TM-score CSV file given in CSV_FILENAME (e.g.
'TM_score_1KBA_fixed_AS_filtered.csv') and reports the fraction of entries
meeting a user-defined TM-score threshold, both overall and filtered to one
selected reference toxin (REFERENCE_FILE).

The ESM2 pipeline was run with two reference toxins, 4K09 and 1KBA. To analyze
the respective run, set CSV_FILENAME and REFERENCE_FILE accordingly
('pdb4k09.ent' or 'pdb1kba.ent').

Usage:
    python calc_structure_similarity_ESM2.py
    -> prompts for a TM-score threshold (float)

Expected CSV columns:
    tm1             (float, 0–1)
    reference_file  (str, e.g. 'pdb1kba.ent')
"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Directory containing the TM-score CSV file to analyze.
BASE_DIR = "/path_to_your_base_directory/ESM2_filtered"

# Name of the input CSV file.
CSV_FILENAME = "TM_score_1KBA_fixed_AS_filtered.csv"

# Reference toxin (value in the 'reference_file' column) for the toxin-specific
# evaluation: "pdb1kba.ent" (1KBA) or "pdb4k09.ent" (4K09).
REFERENCE_FILE = "pdb1kba.ent"

# ---------------------------------------------------------------------------


import os
import pandas as pd


def analyze_run_data(base_dir, tm_threshold, reference_file):
    file_path = os.path.join(base_dir, CSV_FILENAME)
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

        # Filter to the selected reference toxin
        tox_df = df[df['reference_file'] == reference_file]
        tox_filtered_rows = len(tox_df)

        print("\n" + "=" * 40)
        print(f"Total rows:               {total_rows}")
        print(f"Filtered rows:            {filtered_rows}  (tm1 >= {tm_threshold})")
        print(f"Percentage:               {percentage:.3f}%")
        if tox_filtered_rows == 0:
            print(f"No rows found for reference toxin '{reference_file}'.")
        else:
            tox_thresh_filtered = len(tox_df[tox_df['tm1'] >= tm_threshold])
            percentage_tox_only = (tox_thresh_filtered / tox_filtered_rows) * 100
            print(f"Reference toxin:          {reference_file}")
            print(f"Toxin reference rows:     {tox_filtered_rows}")
            print(f"Tox rows over threshold:  {tox_thresh_filtered}")
            print(f"Percentage tox rows:      {percentage_tox_only:.3f}%")
        print("=" * 40)

    except Exception as e:
        print(f"An error occurred while processing the file: {e}")


if __name__ == "__main__":
    try:
        tm_threshold = float(input("Threshold for TM-score: "))
        analyze_run_data(BASE_DIR, tm_threshold, REFERENCE_FILE)
    except ValueError:
        print("Error: Please enter an integer or float number.")
