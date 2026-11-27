"""
Create one pLDDT vs. TM-score scatter plot per reference toxin for a pipeline run.

For the selected run, the script reads the filtered TM-score CSV file and the
filtered ColabFold models (.pdb files). The mean pLDDT of each generated model is
computed from the B-factor column of its CA atoms and matched to the model's
TM-score ('tm1') via the 'generated_file' column. For every reference toxin found
in the 'reference_file' column, a scatter plot (mean pLDDT on the x-axis, TM-score
on the y-axis) is created that shows all generated designs compared against that
toxin, the wildtype toxin itself (TM-score = 1.0, pLDDT = 100) and the structural
similarity (StS) threshold at TM-score = 0.6. Each plot is saved as a separate PDF.

Expected CSV columns:
    generated_file  (str, file name of the generated model, e.g. '*.pdb')
    reference_file  (str, e.g. 'pdb1kba.ent')
    tm1             (float, 0–1)
"""

# ---------------------------------------------------------------------------
# Note: This script was used for both the ProteinMPNN and the ESM Inverse Folding
#       pipeline runs. Only the paths and the model label (MODEL_LABEL) in the
#       Configuration section below need to be changed to point at the respective
#       pipeline's output directory and to label the plots correctly.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Run number of the ESM pipeline run to process
RUN_NUMBER = 6

# Custom filename suffix for the output plots (e.g. "_1kba", or "" for none)
NAMING_SUFFIX = "_4k09"

# Model label shown in the plot legend (e.g. "ESM + CF" for the ESM Inverse Folding
# pipeline or "PMPNN + CF" for the ProteinMPNN pipeline; CF = ColabFold)
MODEL_LABEL = "ESM + CF"

# Base directory holding the filtered TM-score CSV file
BASE_DIR = "/path_to_your_base_directory/comparison_results/Pipeline_ESM"

# Filtered TM-score CSV filename (depends on the run number)
CSV_FILENAME = "TM_score_{run_number}_filtered.csv"

# Directory tree containing the filtered .pdb model files (depends on the run number)
FILTERED_PDB_DIR = "/path_to_your_base_directory/ESM_filtered/run{run_number}_filtered"

# Output directory for the individual per-toxin plots (depends on the run number)
OUTPUT_DIR = "/path_to_your_base_directory/plots/pLDDT_vs_TM_per_toxin_ESM/run{run_number}"

# ---------------------------------------------------------------------------


import os
import glob
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 1. CONFIGURATION & PATHS
# ==========================================
# Input file paths
csvpath_alphafold = os.path.join(BASE_DIR, CSV_FILENAME.format(run_number=RUN_NUMBER))
filtered_pdb_dir = Path(FILTERED_PDB_DIR.format(run_number=RUN_NUMBER))

# Output directory for the individual toxin plots
output_dir = Path(OUTPUT_DIR.format(run_number=RUN_NUMBER))
output_dir.mkdir(parents=True, exist_ok=True)

# Plot styling configuration
plt.rcParams['axes.titlesize'] = 22
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['legend.title_fontsize'] = 18

# ==========================================
# 2. HELPER FUNCTIONS & PDB DISCOVERY
# ==========================================
def get_mean_plddt(pdb_path):
    """Extracts mean pLDDT score (stored in B-factor column) for CA atoms from a PDB file."""
    plddt_scores = []
    if not os.path.exists(pdb_path):
        return None
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    b_factor = float(line[60:66].strip())
                    plddt_scores.append(b_factor)
    return sum(plddt_scores) / len(plddt_scores) if plddt_scores else None

print("Scanning directory tree for PDB files...")
# Map filename to its full file path in memory for fast lookup
pdb_path_map = {file.name: file for file in filtered_pdb_dir.rglob("*.pdb")}
print(f"-> {len(pdb_path_map)} PDB files found.")

# ==========================================
# 3. DATA LOADING
# ==========================================
print("Loading CSV file...")
df_alphafold = pd.read_csv(csvpath_alphafold, low_memory=False)

# Fast lookup for mean pLDDT extraction
print("Calculating mean pLDDT for all models...")
def get_plddt_fast(fname):
    pdb_path = pdb_path_map.get(fname)
    return get_mean_plddt(pdb_path) if pdb_path else None

df_alphafold['plddt'] = df_alphafold['generated_file'].apply(get_plddt_fast)

# ==========================================
# 4. GENERATE PLOTS PER TOXIN
# ==========================================
# Automatically discover all unique toxins / reference PDBs
unique_toxins = sorted(df_alphafold['reference_file'].dropna().unique())
print(f"\nDiscovered toxins/references ({len(unique_toxins)}): {unique_toxins}")

for idx, toxin in enumerate(unique_toxins, start=1):
    toxin_name = toxin.replace('.ent', '').replace('.pdb', '').replace('pdb', '').upper()
    print(f"[{idx}/{len(unique_toxins)}] Creating plot for toxin: {toxin_name}...")

    # Filter dataset for current toxin
    df_toxin = df_alphafold[df_alphafold['reference_file'] == toxin]

    if df_toxin.empty:
        print(f"  -> No data available for {toxin}, skipping...")
        continue

    plt.figure(figsize=(9, 6.5))

    # Scatter plot for generated designs
    sns.scatterplot(
        data=df_toxin,
        x='plddt',
        y='tm1',
        color='tab:blue',
        label=f'{MODEL_LABEL} vs {toxin_name}',
        alpha=0.75,
        s=45
    )

    # Add Reference/Base Toxin point (Self-comparison: TM-score = 1.0, pLDDT = 100)
    ref_x = 100.0
    ref_y = 1.0
    plt.scatter(
        ref_x,
        ref_y,
        color='firebrick',
        s=90,
        zorder=5,
        label=f'Wildtype ({toxin_name})'
    )

    # Annotate the reference point clearly without overlap
    plt.annotate(
        f" {toxin_name}",
        (ref_x, ref_y),
        textcoords="offset points",
        xytext=(-15, -18),  # Shift text slightly down and left so it stays inside plot
        ha='right',
        fontsize=14,
        fontweight='bold',
        color='firebrick'
    )

# Threshold lines & labels
    plt.axhline(y=0.6, color="green", linestyle="--", linewidth=1.5, label="StS Threshold")

    plt.title(f'Protein structural similarity (StS) comparison: {toxin_name}')
    plt.xlabel('mean pLDDT')
    plt.ylabel('TM-score')

    # Standardize axis limits (x-axis up to 102 to provide margin for label)
    plt.ylim(0.0, 1.05)
    plt.xlim(69, 102)

    plt.legend(loc='lower right')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    # Construct modular filename using the configured suffix
    save_filename = output_dir / f"scatter_plddt_vs_tm_run{RUN_NUMBER}_{toxin_name}{NAMING_SUFFIX}.pdf"
    plt.savefig(save_filename, format="pdf")
    plt.close()  # Free memory resources

print(f"\nDone! All {len(unique_toxins)} plots saved successfully in:")
print(output_dir)
