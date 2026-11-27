#!/bin/bash
# Run RFdiffusion motif-scaffolding for a binding site from 8ZF1.
#
# Contig layout: [20-35 / A77-91 / 70-85]
#   - 20–35 randomly sampled residues (N-terminal flank)
#   - residues 77–91 on chain A of the input PDB (fixed motif)
#   - 70–85 randomly sampled residues (C-terminal flank)
#
# Generates 2000 designs.

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

PYTHON="/path_to_your_base_directory/envs/SE3nv/bin/python"                 # Python interpreter from the RFdiffusion (SE3nv) conda environment
RFDIFFUSION_SCRIPT="$HOME/RFdiffusion/scripts/run_inference.py"             # RFdiffusion inference script
OUTPUT_PREFIX="/path_to_your_base_directory/generated_proteins/run4/v2"     # Prefix (directory + name) for the generated designs
INPUT_PDB="/path_to_your_base_directory/TM-score-stuff/8zf1.pdb"            # Input PDB providing the motif
CONTIGS='[20-35/A77-91/70-85]'                                             # Contig map passed to RFdiffusion
NUM_DESIGNS=2000                                                            # Number of designs to generate
# ---------------------------------------------------------------------------

"$PYTHON" "$RFDIFFUSION_SCRIPT" \
    inference.output_prefix="$OUTPUT_PREFIX" \
    inference.input_pdb="$INPUT_PDB" \
    "contigmap.contigs=$CONTIGS" \
    inference.num_designs="$NUM_DESIGNS"
