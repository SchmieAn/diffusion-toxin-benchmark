"""
Calculate TM-scores between filtered AlphaFold - ESM2 outputs and reference proteins.

Searches recursively for all .pdb files in the specific folder,
compares them against each reference protein, and writes results to a CSV file.

Usage:
    python calc_tm_alphafold_output_ESM2.py

Output CSV columns:
    generated_file, reference_file, tm1, tm2, rmsd, seqid, align_len, len1, len2
"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Path to the TMalign binary.
TMALIGN_PATH = "/path_to_your_base_directory/TM-score-stuff/TMalign/TMalign"

# Directory containing the ESM2-filtered PDB files to compare.
INPUT_DIR = "/path_to_your_base_directory/ESM2_filtered/1KBA_fixed_AS"

# Text file listing the reference protein PDB files, one per line.
REFERENCE_PROTEINS_FILE = "/path_to_your_base_directory/lists/reference_proteins.txt"

# Directory where the output CSV is written.
OUTPUT_DIR = "/path_to_your_base_directory/ESM2_filtered"

# Name of the output CSV file.
OUTPUT_FILENAME = "TM_score_1KBA_fixed_AS_filtered.csv"

# ---------------------------------------------------------------------------


import os
import subprocess
import glob
import csv
from pathlib import Path


def TMalign(gen_file, ref_file):
    """Run TMalign on two PDB files and return parsed scores as a dict."""
    if not os.path.exists(gen_file) or not os.path.exists(ref_file):
        raise FileNotFoundError(f'File {gen_file} or {ref_file} does not exist.')

    result = subprocess.run(
        f'{TMALIGN_PATH} {gen_file} {ref_file}',
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8'
    )

    if result.stderr:
        print(f"TMalign error: {gen_file} vs {ref_file}")
        print(result.stderr)
        return None

    lines = result.stdout.split('\n')
    return {
        'len1':      int(lines[9].split()[3]),
        'len2':      int(lines[10].split()[3]),
        'align_len': int(lines[12].split()[2][:-1]),
        'rmsd':      float(lines[12].split()[4][:-1]),
        'seqid':     float(lines[12].split()[-1]),
        'tm1':       float(lines[13].split()[1]),
        'tm2':       float(lines[14].split()[1]),
    }


if __name__ == "__main__":
    base_folder = Path(INPUT_DIR)

    with open(REFERENCE_PROTEINS_FILE) as f:
        reference_files = [line.strip().lower() for line in f if line.strip()]

    output_file = Path(OUTPUT_DIR) / OUTPUT_FILENAME

    fieldnames = ['generated_file', 'reference_file', 'tm1', 'tm2', 'rmsd', 'seqid', 'align_len', 'len1', 'len2']

    with open(output_file, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        # Recursively find all PDB files in the ESM2 output folder
        generated_files = glob.glob(os.path.join(base_folder, "**", "*.pdb"), recursive=True)

        for gen_file in generated_files:
            for ref_file in reference_files:
                print(f"Calculating TMscore: {gen_file} vs {ref_file}")
                result = TMalign(gen_file, ref_file)

                if result is not None:
                    writer.writerow({
                        'generated_file': Path(gen_file).name,
                        'reference_file': Path(ref_file).name,
                        **result
                    })
                else:
                    print(f"Warning: comparison failed for {gen_file} vs {ref_file}")
