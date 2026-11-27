"""
Split ProteinMPNN & ESM Inverse Folding FASTA output files into individual sequence files for AlphaFold input.

For each run folder listed in `SPECIFIC_FOLDERS`, the script reads all FASTA files
from the <run>/seqs/ subdirectory and writes each of the two designed sequences
as a separate .fa file into a corresponding <run>_alphafold/ output folder.

Expected FASTA line layout (0-indexed):
    0: header of sequence 0
    1: (unused)
    2: (unused)
    3: sequence 0
    4: header of sequence 1
    5: sequence 1

Usage:
    Set `INPUT_PATH` to the ProteinMPNN output directory and list the run folders
    to process in `SPECIFIC_FOLDERS`, then run:
        python fasta_file_split_cluster.py
"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# output directory containing the run folders to process
INPUT_PATH = "/path_to_your_base_directory/ProteinMPNN_Output/"

# AlphaFold input folders are written back next to the output
OUTPUT_PATH = INPUT_PATH

# Only these run folders will be processed
SPECIFIC_FOLDERS = ["run6", "run7"]

# ---------------------------------------------------------------------------


import os


for run_folder in SPECIFIC_FOLDERS:
    run_path = os.path.join(INPUT_PATH, run_folder)

    if not os.path.isdir(run_path):
        continue

    # Create the corresponding AlphaFold input folder
    output_run_path = os.path.join(OUTPUT_PATH, f"{run_folder}_alphafold")
    os.makedirs(output_run_path, exist_ok=True)

    seqs_path = os.path.join(run_path, "seqs")

    if not os.path.isdir(seqs_path):
        print(f"{seqs_path} not found.")
        continue

    for seq_file in os.listdir(seqs_path):
        seq_path = os.path.join(seqs_path, seq_file)
        base_name = os.path.splitext(seq_file)[0]

        with open(seq_path, "r", encoding="utf-8") as infile:
            lines = infile.readlines()

        # Write the two designed sequences as separate .fa files
        for i, (suffix, seq_line) in enumerate([("_00", lines[3]), ("_01", lines[5])]):
            out_path = os.path.join(output_run_path, f"{base_name}{suffix}.fa")
            with open(out_path, "w", encoding="utf-8") as out:
                out.write(f">{base_name}{suffix}\n")
                out.write(seq_line)
