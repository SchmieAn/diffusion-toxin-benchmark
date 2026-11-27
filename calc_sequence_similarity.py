"""
Analyze sequence similarity results from a ProteinMPNN pipeline run.
The same script was also used for the ESM Inverse Folding pipeline runs.
Therefore solely the path has to change to the ESM Inverse Folding directory.

Reads a CSV file containing pairwise sequence comparisons and reports
the fraction of entries meeting a user-defined percent identity threshold.

Usage:
    python calc_sequence_similarity.py
    -> prompts for a Percent_Identity threshold (float)

Expected CSV columns:
    Percent_Identity  (float, 0–100)

"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Directory containing the sequence comparison result CSV files.
BASE_DIR = "/path_to_your_base_directory/Pipeline_PMPNN_filtered"

# Run number to analyze. The script looks for the file:
#     {BASE_DIR}/sequence_comparison_results_run{RUN_NUMBER}.csv
RUN_NUMBER = 15

# Filename pattern for the result CSV files. "{run_number}" is replaced with RUN_NUMBER.
FILENAME_TEMPLATE = "sequence_comparison_results_run{run_number}.csv"

# ---------------------------------------------------------------------------


import os
import pandas as pd


def analyze_run_data(run_number, base_dir, identity_threshold):
    file_path = os.path.join(base_dir, FILENAME_TEMPLATE.format(run_number=run_number))
    print(f"Reading file: {file_path}")

    if not os.path.exists(file_path):
        print(f"Error: file '{file_path}' not found.")
        return

    try:
        df = pd.read_csv(file_path)

        if 'Percent_Identity' not in df.columns:
            print("Error: Column 'Percent_Identity' does not exist.")
            print(f"Found columns: {list(df.columns)}")
            return

        total_rows = len(df)
        if total_rows == 0:
            print("The CSV file is empty.")
            return

        filtered_rows = len(df[df['Percent_Identity'] >= identity_threshold])
        percentage = (filtered_rows / total_rows) * 100

        print("\n" + "=" * 40)
        print(f"Results for run {run_number}")
        print("=" * 40)
        print(f"Total rows:           {total_rows}")
        print(f"Filtered rows:        {filtered_rows}  (Percent_Identity >= {identity_threshold}%)")
        print(f"Percentage:           {percentage:.3f}%")
        print("=" * 40)

    except Exception as e:
        print(f"An error occurred while processing the file: {e}")


if __name__ == "__main__":
    try:
        identity_threshold = float(input("Threshold for Percent_Identity: "))
        analyze_run_data(RUN_NUMBER, BASE_DIR, identity_threshold)
    except ValueError:
        print("Error: Please enter an integer or float number.")
