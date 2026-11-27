"""
Generate protein sequences with ESM-2 by iterative masked sampling while keeping
the specified residues fixed.

The input sequences (see 'data') contain fixed residues given explicitly and
<mask> tokens at the positions to be designed. For each of the NUM_SAMPLES
samples, the masked positions are filled in one at a time with the pretrained
ESM-2 model (esm2_t33_650M_UR50D): in every iteration the model logits are
converted to probabilities using the sampling temperature, an amino acid is
sampled for every remaining mask position, and the position whose sampled amino
acid has the highest probability is filled in. This is repeated until no <mask>
token is left. The fixed residues are never changed.

Each completed sequence is saved as a FASTA file together with a contact map
(PNG) predicted by ESM-2 in OUTPUT_FOLDER.
"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Sampling temperature (analogous to ProteinMPNN)
TEMPERATURE = 0.1

# Number of sequences to generate per protein
NUM_SAMPLES = 20

# Output directory for the generated FASTA files and contact maps
OUTPUT_FOLDER = "/path_to_your_base_directory/results/fixed_AS"

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
model.eval()

os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# ==============================================
# 2. INPUT DATA
# ==============================================
# Fixed residues are given explicitly; <mask> tokens mark positions to be designed
data = [
    ("4K09_81p_25", "S L F E L G K M I L Q E T G K N P A K S Y G A Y G C N C G V L G R G K P K D A T D R C C Y V H K C C Y K K L T G C D Q K K D R Y S Y S W K D K T I V C G E N N P C L K E L C E C D <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask> <mask>"),
]

# Collect final sequences: { label: [seq1, seq2, ...] }
final_generated_sequences = {label: [] for label, _ in data}

# ==============================================
# 3. ITERATIVE MASKED SAMPLING
# ==============================================
# At each step, the most confident mask position is filled in first (greedy by confidence).
print(f"Starting iterative sampling ({NUM_SAMPLES} samples per protein, T={TEMPERATURE})...")

for label, original_seq in data:
    for sample_idx in range(NUM_SAMPLES):
        current_tokens_list = original_seq.split()

        while "<mask>" in current_tokens_list:
            temp_seq_str = " ".join(current_tokens_list)
            _, _, batch_tokens = batch_converter([(label, temp_seq_str)])

            with torch.no_grad():
                results = model(batch_tokens, repr_layers=[33])

            logits = results["logits"][0]
            probs = F.softmax(logits / TEMPERATURE, dim=-1)

            # Find the mask position with the highest sampling confidence
            best_mask_idx_in_seq = None
            best_mask_sampled_amino = None
            highest_confidence = -1.0

            for seq_pos, token in enumerate(current_tokens_list):
                if token == "<mask>":
                    token_idx_in_model = seq_pos + 1  # +1 for BOS token
                    pos_probs = probs[token_idx_in_model]

                    sampled_idx = torch.multinomial(pos_probs, num_samples=1).item()
                    sampled_amino = alphabet.get_tok(sampled_idx)
                    confidence = pos_probs[sampled_idx].item()

                    if confidence > highest_confidence:
                        highest_confidence = confidence
                        best_mask_idx_in_seq = seq_pos
                        best_mask_sampled_amino = sampled_amino

            current_tokens_list[best_mask_idx_in_seq] = best_mask_sampled_amino

        predicted_seq = "".join(current_tokens_list)
        final_generated_sequences[label].append(predicted_seq)

# ==============================================
# 4. SAVE FASTAs AND CONTACT MAPS PER SAMPLE
# ==============================================
print("\nGenerating FASTA files and contact maps for each sample...")

for label, _ in data:
    sample_sequences = final_generated_sequences[label]

    # Prepare ESM-format batch for contact map computation
    esm_sample_data = [(f"{label}_sample_{idx+1}", seq) for idx, seq in enumerate(sample_sequences)]
    sample_labels, sample_strs, sample_tokens = batch_converter(esm_sample_data)

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

        filename = os.path.join(OUTPUT_FOLDER, f"{sample_label}_contact_map.png")
        plt.savefig(filename, bbox_inches='tight')
        plt.close()

    print(f" -> {NUM_SAMPLES} FASTAs and contact maps saved for: {label}")

print("\n[DONE] Iterative generation complete.")
