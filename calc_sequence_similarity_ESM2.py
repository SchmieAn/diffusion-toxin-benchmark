"""
Analyze sequence similarity results from an ESM2 pipeline run.

Reads a CSV file here for example 'sequence_comparison_1KBA_fixed_AS_filtered.csv' and
reports the fraction of entries meeting a user-defined percent identity threshold.

Usage:
    python calc_sequence_similarity_ESM2.py
    -> prompts for a Percent_Identity threshold (float)

Expected CSV columns:
    Percent_Identity  (float, 0–100)
"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Directory containing the sequence comparison result CSV file.
BASE_DIR = "/path_to_your_base_directory/ESM2_filtered"

# Name of the result CSV file to analyze.
CSV_FILENAME = "sequence_comparison_1KBA_fixed_AS_filtered.csv"

# ---------------------------------------------------------------------------


import os
import pandas as pd


def analyze_run_data(base_dir, identity_threshold):
    file_path = os.path.join(base_dir, CSV_FILENAME)
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
        analyze_run_data(BASE_DIR, identity_threshold)
    except ValueError:
        print("Error: Please enter an integer or float number.")
