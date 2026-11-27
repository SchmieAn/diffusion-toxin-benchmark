# Protein Design & Comparison Pipeline

Scripts for generating protein backbones, designing sequences, predicting structures,
and quantitatively comparing the generated proteins against a set of reference toxins.

Three parallel design pipelines are included and compared throughout:

- **PMPNN pipeline** - backbones from RFdiffusion, sequences from **ProteinMPNN**.
- **ESM Inverse Folding pipeline** - backbones from RFdiffusion, sequences from **ESM Inverse Folding**.
- **ESM-2 pipeline** - sequences from **ESM2** (files suffixed `_ESM2`).

Most analysis scripts exist in both variants; the `_ESM2` version is the ESM-2-pipeline
counterpart of the base (PMPNN & ESM Inverse Folding) script and differs mainly in its input/output paths.

> **Important: scripts shared across pipelines:** The analysis and plotting
> scripts/notebooks **without** `ESM2` in their filename were used for **both** the
> ProteinMPNN and the ESM Inverse Folding pipeline runs. The same code applies to
> either pipeline, only the path(s) in the Configuration section need to be changed
> to point at the respective pipeline's output directory. As a result, the path a
> given script currently contains simply reflects whichever pipeline it was last run
> on; switch it in the Configuration block to analyze the other. These shared scripts
> carry a `Note:` in their header stating this. (The `start_esm_temp_*.py` and
> `run_protein_mpnn.sh` generation scripts are the exception. Each is specific to one
> pipeline's sequence-design tool and is not interchangeable.)

> **Note on reproduction:** These scripts were written for a specific study and
> hardcode a particular directory layout. They are shared for transparency and reuse,
> **not** as a turnkey tool. You will need to adapt paths and input data to your own setup.

---

## Before you run anything: configuration

Every script has a clearly marked **Configuration** block that gathers all the values
you need to edit: directory paths, run numbers, input/output filenames, and key
parameters. In the Python scripts it looks like this and sits near the top of the file:

```python
# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

BASE_DIR = "/path_to_your_base_directory/Pipeline_PMPNN_filtered"
RUN_NUMBER = 15
FILENAME_TEMPLATE = "sequence_comparison_results_run{run_number}.csv"

# ---------------------------------------------------------------------------
```

Shell/job scripts have the same block (as shell variables) below the shebang and any
`#SBATCH` directives; notebooks have it in a Configuration cell near the top.

All machine-specific paths have been anonymised to `/path_to_your_base_directory/...`.
**Replace this prefix with your own base directory** before running. You generally do
not need to edit anything outside the Configuration block.

---

## Requirements

- Python 3.x with `pandas`, `numpy`, `matplotlib`, `biopython` (for sequence/structure parsing)
- [**RFdiffusion**](https://github.com/RosettaCommons/RFdiffusion) - backbone generation
- [**ProteinMPNN**](https://github.com/dauparas/ProteinMPNN) - sequence design (PMPNN pipeline)
- [**ESM Inverse Folding**] (https://github.com/Biohub/esm-archive/blob/main/examples/inverse_folding/README.md) - sequence design (ESM Inverse Folding pipeline)
- [**ESM-2**] (https://github.com/Biohub/esm-archive#esmfold) - sequence design (ESM-2 pipeline)
- [**ColabFold**](https://github.com/sokrypton/ColabFold) - structure prediction
- [**TM-align**](https://zhanggroup.org/TM-align/) - structural comparison (path set via `TMALIGN_PATH`)

> Please adjust and pin exact versions to match your environment.

---

## Pipeline overview

```
  RFdiffusion            ProteinMPNN / ESM Inverse Folding                  ColabFold                  Analysis & plots
 (backbones)      →                (sequences)                →       (predicted structures)   →   (similarity, TM, pLDDT)
```

---

## Scripts by stage

### 1. Backbone generation (RFdiffusion)
| Script | Purpose |
|---|---|
| `design_motifscaffolding.sh` | Motif-scaffolding design run |
| `design_motifscaffolding_2.sh` | Motif-scaffolding design (toxin motif/target) |
| `design_unconditional.sh` | Unconditional backbone generation |

### 2. Sequence design
| Script | Purpose |
|---|---|
| `run_protein_mpnn.sh` | Run ProteinMPNN on generated backbones (PMPNN pipeline) |
| `submit_esm_header.sh` | Run ESM Inverse Folding on generated backbones (ESM Inverse Folding pipeline) |
| `start_esm_temp_fixed_as.py` | ESM-2 — **iterative** masked sampling, fixed anchor residues |
| `start_esm_temp_lose_as.py` | ESM-2 — **single-shot** sampling, loose anchor residues |

### 3. Structure prediction (ColabFold / AlphaFold)
| Script | Purpose |
|---|---|
| `run_colabfold.job` | SLURM job to run ColabFold structure prediction |
| `submit_colabworkaround.sh` | Submit/orchestrate ColabFold jobs over generated outputs |

### 4. FASTA preparation
| Script | Purpose |
|---|---|
| `fasta_file_split_cluster.py` | Split clustered FASTA files into per-design sequence files |

### 5. Sequence comparison
| Script | Purpose |
|---|---|
| `sequence_comparison_proteins.py` / `_ESM2.py` | Compute pairwise percent-identity of designed vs reference sequences → CSV |
| `calc_sequence_similarity.py` / `_ESM2.py` | Report the fraction of pairs meeting a percent-identity threshold (prompts for the threshold) |

### 6. Structure comparison (TM-score)
| Script | Purpose |
|---|---|
| `calc_tm.py` | Run TM-align over generated structures vs references → TM-score CSV |
| `calc_tm_alphafold_output.py` / `_ESM2.py` | TM-align on ColabFold-predicted structures → TM-score CSV |
| `calc_structure_similarity.py` / `_ESM2.py` | Report the fraction of structures meeting a TM-score threshold (prompts for the threshold) |

### 7. Plots & summary analysis
| Script | Purpose |
|---|---|
| `plddt_vs_tm_per_toxin.py` / `_ESM2.py` | Scatter of pLDDT vs TM-score, per reference toxin |
| `sequence_similarity_plots.ipynb` / `_ESM2.ipynb` | Sequence-similarity distribution plots |
| `structure_comparison_plddt_tm.ipynb` / `_ESM2.ipynb` | pLDDT / TM-score comparison plots |
| `structure_comparison_scatter.ipynb` / `_ESM2.ipynb` | Structure-comparison scatter plots |
| `alpha-helix_content_calculation.ipynb` | Secondary-structure (α-helix content) analysis |

---

## Typical workflow

1. Generate backbones with the `design_*.sh` scripts (RFdiffusion).
2. Design sequences with `run_protein_mpnn.sh` / `submit_esm_header.sh` **or** use `start_esm_temp_*.py` scripts to generate sequences without diffusion.
3. Split/cluster sequences with `fasta_file_split_cluster.py`.
4. Predict structures with `run_colabfold.job` / `submit_colabworkaround.sh`.
5. Compute comparisons: `sequence_comparison_proteins*.py` and `calc_tm*.py`.
6. Summarise with `calc_sequence_similarity*.py` / `calc_structure_similarity*.py`
   (each prompts for a threshold) and generate plots from the notebooks.

---

## Citation

If you use these scripts, please cite:  

**Diffusion algorithms in protein design rarely propose known toxins: a benchmark on sequence and structure level**

## License
MIT

## Contact

Anne Schmieder — anne.schmieder@tu-dresden.de
