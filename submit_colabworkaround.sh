#!/bin/bash
# Submit ColabFold prediction jobs for all .fa files across run folders 1–15.
#
# For each run folder (<BASE_DIR>/run<N>_alphafold), every .fa file is submitted
# as a separate SLURM job via run_colabfold.job. Output is written to a subfolder
# named <file>.fa-output/ next to the input file.
#
# Usage:
#   bash submit_colabworkaround.sh

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

BASE_DIR="/path_to_your_base_directory/RF-diff-generated/output"                  # Parent directory holding the run<N>_alphafold folders
JOB_SCRIPT="/path_to_your_base_directory/colabfold_workaround/run_colabfold.job"  # SLURM job script submitted per .fa file
# ---------------------------------------------------------------------------

for i in {1..15}; do
    cd "${BASE_DIR}/run${i}_alphafold" || { echo "Folder run${i}_alphafold not found, skipping."; continue; }

    for file in *.fa; do
        # Skip if no .fa files exist in this folder
        [ -e "$file" ] || continue

        mkdir -p "${file}-output"
        sbatch "$JOB_SCRIPT" \
            -i "$file" \
            -o "${file}-output"
    done

    cd - > /dev/null
done
