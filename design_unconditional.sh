#!/bin/bash
# Run RFdiffusion unconditional protein design (no input structure).
#
# Designs proteins of 120–150 residues in length (randomly sampled per design).
# Generates 1000 designs.

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

PYTHON="/path_to_your_base_directory/envs/SE3nv/bin/python"                 # Python interpreter from the RFdiffusion (SE3nv) conda environment
RFDIFFUSION_SCRIPT="$HOME/RFdiffusion/scripts/run_inference.py"             # RFdiffusion inference script
OUTPUT_PREFIX="/path_to_your_base_directory/generated_proteins/run3/v1"     # Prefix (directory + name) for the generated designs
CONTIGS='[120-150]'                                                        # Contig map (design length range)
NUM_DESIGNS=1000                                                            # Number of designs to generate
# ---------------------------------------------------------------------------

"$PYTHON" "$RFDIFFUSION_SCRIPT" \
    inference.output_prefix="$OUTPUT_PREFIX" \
    "contigmap.contigs=$CONTIGS" \
    inference.num_designs="$NUM_DESIGNS"
