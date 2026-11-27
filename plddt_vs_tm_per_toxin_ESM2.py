"""
Create one pLDDT vs. TM-score scatter plot per reference toxin for the ESM2 pipeline.

The script reads the filtered TM-score CSV file and the filtered ColabFold models
(.pdb files, searched recursively in all subfolders of FILTERED_PDB_DIR). The mean
pLDDT of each generated model is computed from the B-factor column of its CA atoms
(rescaled to 0-100 if stored as 0-1) and matched to the model's TM-score ('tm1')
via the 'generated_file' column. For every reference toxin found in the
'reference_file' column (here 4K09 and 1KBA), a scatter plot (mean pLDDT on the
x-axis, TM-score on the y-axis) is created that shows all ESM2-generated designs
compared against that toxin, the wildtype toxin itself (TM-score = 1.0,
pLDDT = 100) and the structural similarity (StS) threshold at TM-score = 0.6.
Each plot is saved as a separate PDF.

Expected CSV columns:
    generated_file  (str, file name of the generated model, e.g. '*.pdb')
    reference_file  (str, e.g. 'pdb4k09.ent', 'pdb1kba.ent')
    tm1             (float, 0–1)
"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Custom filename suffix for the output plots
NAMING_SUFFIX = "_4k09"

# Base directory holding the filtered TM-score CSV file
BASE_DIR = "/path_to_your_base_directory/comparison_results/ESM2_filtered"

# Filtered TM-score CSV filename
CSV_FILENAME = "TM_score_4K09_fixed_AS_filtered.csv"

# Directory containing the subfolders with .pdb files
FILTERED_PDB_DIR = "/path_to_your_base_directory/comparison_results/ESM2_filtered/4K09_fixed_AS"

# Output directory for the plots
OUTPUT_DIR = "/path_to_your_base_directory/plots/pLDDT_vs_TM_ESM2_fixed_AS"

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
# Input file paths matching your actual directory structure
csvpath_alphafold = os.path.join(BASE_DIR, CSV_FILENAME)

# Directory containing the subfolders with .pdb files
filtered_pdb_dir = Path(FILTERED_PDB_DIR)

# Output directory for the plots
output_dir = Path(OUTPUT_DIR)
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

print("Scanning directory tree for PDB files in subfolders...")
# Recursive lookup handles all *.fasta-output subdirectories automatically
pdb_path_map = {file.name: file for file in filtered_pdb_dir.rglob("*.pdb")}
print(f"-> {len(pdb_path_map)} PDB files found.")

# ==========================================
# 3. DATA LOADING & ENRICHMENT
# ==========================================
print("Loading CSV file...")
df_alphafold = pd.read_csv(csvpath_alphafold, low_memory=False)

print("Calculating mean pLDDT for all models...")
def get_plddt_fast(fname):
    if pd.isna(fname):
        return None
    # Strip any directory components in case 'generated_file' contains subpaths
    clean_fname = os.path.basename(str(fname))
    pdb_path = pdb_path_map.get(clean_fname)
    return get_mean_plddt(pdb_path) if pdb_path else None

df_alphafold['plddt'] = df_alphafold['generated_file'].apply(get_plddt_fast)

# Standardize pLDDT scale to 0-100 if scores are in 0.0-1.0 format
if df_alphafold['plddt'].dropna().max() <= 1.0:
    df_alphafold['plddt'] = df_alphafold['plddt'] * 100

# ==========================================
# 4. GENERATE PLOTS PER TOXIN
# ==========================================
unique_toxins = sorted(df_alphafold['reference_file'].dropna().unique())
print(f"\nDiscovered toxins/references ({len(unique_toxins)}): {unique_toxins}")

for idx, toxin in enumerate(unique_toxins, start=1):
    toxin_name = toxin.replace('.ent', '').replace('.pdb', '').replace('pdb', '').upper()
    print(f"[{idx}/{len(unique_toxins)}] Creating plot for toxin: {toxin_name}...")

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
        label=f'ESM2 + CF vs {toxin_name}',
        alpha=0.75,
        s=45
    )

    # Reference / Wildtype point
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

    plt.annotate(
        f" {toxin_name}",
        (ref_x, ref_y),
        textcoords="offset points",
        xytext=(-15, -18),
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

    plt.ylim(0.0, 1.05)
    plt.xlim(69, 102)

    plt.legend(loc='lower right')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    # Save plot
    save_filename = output_dir / f"scatter_plddt_vs_tm_{toxin_name}{NAMING_SUFFIX}.pdf"
    plt.savefig(save_filename, format="pdf")
    plt.close()

print(f"\nDone! All {len(unique_toxins)} plots saved successfully in:")
print(output_dir)
