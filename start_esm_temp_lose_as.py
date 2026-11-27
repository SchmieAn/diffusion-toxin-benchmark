"""
Generate protein sequences with ESM-2 by single-shot sampling ("loose" variant).

The input sequences (see 'data') contain <mask> tokens at the positions to be
designed. The pretrained ESM-2 model (esm2_t33_650M_UR50D) computes the logits
for the whole sequence in a single forward pass, and NUM_SAMPLES sequences are
drawn by sampling an amino acid at every position from the softmax distribution
scaled by the sampling temperature. In contrast to the iterative script
(start_esm_temp_fixed_as.py), the non-masked residues are sampled as well and are
therefore not guaranteed to remain fixed.

Each generated sequence is saved as a FASTA file together with a contact map
(PNG) predicted by ESM-2 in OUTPUT_FOLDER.
"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Sampling temperature (analogous to ProteinMPNN)
TEMPERATURE = 0.1

# Number of sequences to generate per protein
NUM_SAMPLES = 5

# Output directory for the generated FASTA files and contact maps
OUTPUT_FOLDER = "/path_to_your_base_directory/results/loose_AS"

# ---------------------------------------------------------------------------


import torch
import os
import esm
import torch.nn.functional as F
import matplotlib.pyplot as plt

# ==============================================
# 1. MODEL & PARAMETER CONFIGURATION
# ==============================================
print("Loading ESM-2 model...")
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # Disables dropout for deterministic inference

os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# ==============================================
# 2. INPUT DATA
# ==============================================
# Sequences with <mask> tokens at positions to be designed
data = [
    ("protein3", "K A <mask> <mask> <mask> I S Q D I R L"),
]

batch_labels, batch_strs, batch_tokens = batch_converter(data)
batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

# ==============================================
# 3. MODEL FORWARD PASS (FOR SAMPLING)
# ==============================================
print("Computing logits for sampling...")
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=False)

logits = results["logits"]
probs = F.softmax(logits / TEMPERATURE, dim=-1)
batch_size, seq_len, alphabet_size = probs.shape
probs_flat = probs.view(-1, alphabet_size)

# Collect all sampled sequences: { protein_index: [seq1, seq2, ...] }
all_predicted_samples = {i: [] for i in range(batch_size)}

for sample_idx in range(NUM_SAMPLES):
    sampled_tokens_flat = torch.multinomial(probs_flat, num_samples=1)
    predicted_tokens = sampled_tokens_flat.view(batch_size, seq_len)

    for i, (label, original_seq) in enumerate(data):
        tokens_len = batch_lens[i]
        pred_indices = predicted_tokens[i, 1 : tokens_len - 1].tolist()
        predicted_seq = "".join([alphabet.get_tok(idx) for idx in pred_indices])
        all_predicted_samples[i].append(predicted_seq)

# ==============================================
# 4. SAVE FASTAs AND CONTACT MAPS PER SAMPLE
# ==============================================
print(f"\nGenerating FASTA files and contact maps for all {NUM_SAMPLES} samples...")

for i, (label, original_seq) in enumerate(data):
    sample_sequences = all_predicted_samples[i]

    # Prepare ESM-format batch for contact map computation
    esm_sample_data = [(f"{label}_sample_{idx+1}", seq) for idx, seq in enumerate(sample_sequences)]
    sample_labels, sample_strs, sample_tokens = batch_converter(esm_sample_data)
    sample_lens = (sample_tokens != alphabet.padding_idx).sum(1)

    with torch.no_grad():
        sample_results = model(sample_tokens, return_contacts=True)

    for sample_idx, (sample_label, predicted_seq) in enumerate(esm_sample_data):

        # Save FASTA file
        individual_fasta_path = os.path.join(OUTPUT_FOLDER, f"{sample_label}.fasta")
        with open(individual_fasta_path, "w") as fasta_file:
            fasta_file.write(f">{sample_label}\n")
            fasta_file.write(f"{predicted_seq}\n")

        # Save contact map
        attention_contacts = sample_results["contacts"][sample_idx]
        seq_len = len(predicted_seq)

        plt.matshow(attention_contacts[:seq_len, :seq_len])
        plt.title(f"{sample_label} Contact Map")
        plt.xlabel("Residue Index")
        plt.ylabel("Residue Index")

        plot_filename = os.path.join(OUTPUT_FOLDER, f"{sample_label}_contact_map.png")
        plt.savefig(plot_filename, bbox_inches='tight')
        plt.close()

    print(f" -> {NUM_SAMPLES} FASTAs and contact maps saved for: {label}")

print("\n[DONE] Single-shot generation complete.")
