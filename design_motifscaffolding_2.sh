#!/bin/bash
# Run RFdiffusion motif-scaffolding for a binding site from 4K09.
#
# Note: This script was used for both 4K09 and 1KBA as motif sources for the design
#       generation. The settings below correspond to one of the 4K09 runs (here 46% sequence were
#	kept as a foundation. The 81%, 71%, 63% and 54% as well as the parameters
#       used for the 1KBA runs are listed in the README.
#
# Contig layout (4K09 settings): [A15-53 / A57-61 / A67-68 / 35-60]
#   - residues 15–53 on chain A of the input PDB (fixed motif segment 1, 39 residues)
#   - residues 57–61 on chain A of the input PDB (fixed motif segment 2, 5 residues)
#   - residues 67–68 on chain A of the input PDB (fixed motif segment 3, 2 residues)
#   - 35–60 randomly sampled residues (C-terminal flank)
#
# Preserved sequence (4K09 settings): 46 residues of the input PDB are kept fixed.
# Together with the 35–60 sampled residues, each design has a total length of
# 81–106 residues, i.e. about 46 % of the designed sequence is preserved
# from the input protein
#
# Generates 500 designs.

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

PYTHON="/path_to_your_base_directory/envs/SE3nv/bin/python"                     # Python interpreter from the RFdiffusion (SE3nv) conda environment
RFDIFFUSION_SCRIPT="$HOME/software/RFdiffusion/scripts/run_inference.py"        # RFdiffusion inference script
OUTPUT_PREFIX="/path_to_your_base_directory/generated_proteins/run10/tox46p"    # Prefix (directory + name) for the generated designs
INPUT_PDB="/path_to_your_base_directory/comparison_results/4k09.pdb"            # Input PDB providing the motif
CONTIGS='[A15-53/A57-61/A67-68/35-60]'                                         # Contig map passed to RFdiffusion
NUM_DESIGNS=500                                                                 # Number of designs to generate
# ---------------------------------------------------------------------------

"$PYTHON" "$RFDIFFUSION_SCRIPT" \
    inference.output_prefix="$OUTPUT_PREFIX" \
    inference.input_pdb="$INPUT_PDB" \
    "contigmap.contigs=$CONTIGS" \
    inference.num_designs="$NUM_DESIGNS"
