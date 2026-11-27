#!/bin/bash
# Run ProteinMPNN sequence design on all PDB files in a given folder.
#
# Usage: bash run_protein_mpnn.sh <folder_with_pdbs> <output_dir>

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

MPNN="/path_to_your_base_directory/ProteinMPNN"   # Path to the ProteinMPNN installation
NUM_SEQ_PER_TARGET=2                              # Number of sequences to design per input structure
SAMPLING_TEMP="0.1"                               # Sampling temperature
SAVE_PROBS=1                                      # Save per-position amino-acid probabilities (1 = yes)
# ---------------------------------------------------------------------------

# Input: folder containing PDB files
folder_with_pdbs=$1

# Output: directory for results (created if it does not exist)
output_dir=$2
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# Intermediate file: parsed chain information
path_for_parsed_chains="$output_dir/parsed_pdbs.jsonl"

# Parse all chains from the input PDB files
python "$MPNN/helper_scripts/parse_multiple_chains.py" \
    --input_path="$folder_with_pdbs" \
    --output_path="$path_for_parsed_chains"

# Run ProteinMPNN to design sequences
python "$MPNN/protein_mpnn_run.py" \
    --jsonl_path "$path_for_parsed_chains" \
    --out_folder "$output_dir" \
    --num_seq_per_target "$NUM_SEQ_PER_TARGET" \
    --sampling_temp "$SAMPLING_TEMP" \
    --save_probs "$SAVE_PROBS"
