"""
Compare generated protein sequences against toxin reference sequences using global pairwise alignment.

Loads toxin sequences from PDB files via SEQRES records, then iterates over generated
protein PDB files one at a time (memory-efficient), computes BLOSUM62-based global
alignments against each toxin, and writes results directly to a CSV file.

Usage:
    Set the paths in configuration below and run:
        python sequence_comparison_proteins.py

Output CSV columns:
    Generated_Protein_ID, Toxin_ID, Percent_Identity, Alignment_Score,
    Number_of_Matches, Length_gen_prot, Length_toxin

Note:
    This script was used for both the ProteinMPNN and the ESM Inverse Folding
    pipeline runs. Only the path(s) in the Configuration section below need to be
    changed to point at the respective pipeline's output directory.
"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Run number of the filtered run to process.
RUN_NUMBER = 5

# Directory containing the run-specific filtered generated-protein folders.
GENERATED_BASE_DIR = "/path_to_your_base_directory/PMPNN_filtered"

# Directory containing the toxin reference PDB files.
TOX_FOLDER = "/path_to_your_base_directory/toxine_pdb"

# Directory where the output CSV is written.
OUTPUT_DIR = "/path_to_your_base_directory/Pipeline_PMPNN_filtered"

# Name of the output CSV file — depends on the run number.
OUTPUT_FILENAME_TEMPLATE = "sequence_comparison_results_run{run_number}.csv"

# ---------------------------------------------------------------------------


from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align import substitution_matrices
import glob
import os
import csv


def load_toxin_sequences(tox_folder):
    """Load the first chain sequence from each PDB file in tox_folder.

    Returns a list of (sequence_str, record_id) tuples.
    """
    tox_files = glob.glob(os.path.join(tox_folder, "**", "*.pdb"), recursive=True)
    extracted_seq_id_pairs = []

    for tox_file in tox_files:
        print(f"Processing toxin file: {tox_file}")
        try:
            # 'pdb-seqres' reads SEQRES records and returns one-letter sequences
            for record in SeqIO.parse(tox_file, "pdb-seqres"):
                extracted_seq_id_pairs.append((str(record.seq), record.id))
                break  # Extract only the first chain per file
        except Exception as e:
            print(f"Error parsing {tox_file}: {e}")

    return extracted_seq_id_pairs


def calculate_pairwise_similarity(seq1_str, seq2_str):
    """Compute global pairwise alignment between two protein sequences.

    Uses BLOSUM62 with gap-open=-10, gap-extend=-0.5.

    Returns:
        tuple: (percent_identity, alignment_score, num_matches)
        Returns (0.0, 0.0, 0) if either sequence is empty or alignment fails.
    """
    if not seq1_str or not seq2_str:
        return 0.0, 0.0, 0

    matrix = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.globalds(seq1_str, seq2_str, matrix, -10, -0.5)

    if not alignments:
        return 0.0, 0.0, 0

    aligned_seq1, aligned_seq2, score, begin, end = alignments[0]

    matches = 0
    alignment_length = 0
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == b:
            matches += 1
        if a != '-' or b != '-':  # Count non-gap columns
            alignment_length += 1

    percent_identity = (matches / alignment_length * 100) if alignment_length > 0 else 0.0
    return percent_identity, score, matches


def compare_all_sequences_optimized_for_memory(generated_proteins_folder, tox_seq_id_pairs, output_csv_path):
    """Compare all generated proteins against all toxin sequences, writing results to CSV.

    Processes one generated PDB file at a time to avoid high memory usage.
    """
    print("\nStarting sequence comparison, writing results directly to CSV...")
    parser = PDBParser(PERMISSIVE=True)
    gen_files = glob.glob(os.path.join(generated_proteins_folder, "**", "*.pdb"), recursive=True)

    total_comparisons = len(gen_files) * len(tox_seq_id_pairs)
    current_comparison_count = 0

    fieldnames = [
        "Generated_Protein_ID", "Toxin_ID", "Percent_Identity",
        "Alignment_Score", "Number_of_Matches", "Length_gen_prot", "Length_toxin"
    ]

    with open(output_csv_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for gen_file_path in gen_files:
            current_gen_seq = ""
            current_gen_id = os.path.basename(gen_file_path).replace(".pdb", "")

            try:
                structure = parser.get_structure(current_gen_id, gen_file_path)

                for model in structure:
                    for chain in model:
                        chain_sequence_parts = []
                        for residue in chain:
                            if residue.get_id()[0] == ' ' and residue.resname != 'HOH':
                                try:
                                    chain_sequence_parts.append(seq1(residue.resname))
                                except KeyError:
                                    chain_sequence_parts.append('X')
                        if chain_sequence_parts:
                            current_gen_seq = "".join(chain_sequence_parts)
                            break  # First chain only
                    if current_gen_seq:
                        break  # First model only

                if not current_gen_seq:
                    print(f"Warning: No valid sequence found in {gen_file_path}. Skipping.")
                    continue

                length_gen_prot = len(current_gen_seq)

                for tox_seq, tox_id in tox_seq_id_pairs:
                    current_comparison_count += 1
                    if current_comparison_count % 1000 == 0:
                        print(f"Processed {current_comparison_count}/{total_comparisons} comparisons.")

                    percent_identity, alignment_score, num_matches = calculate_pairwise_similarity(
                        current_gen_seq, tox_seq
                    )

                    writer.writerow({
                        "Generated_Protein_ID": current_gen_id,
                        "Toxin_ID": tox_id,
                        "Percent_Identity": percent_identity,
                        "Alignment_Score": alignment_score,
                        "Number_of_Matches": num_matches,
                        "Length_gen_prot": length_gen_prot,
                        "Length_toxin": len(tox_seq)
                    })

            except Exception as e:
                print(f"Error processing {gen_file_path}: {e}")

    print(f"Comparison complete. Results saved to: {output_csv_path}")


def main():
    print("Program started")

    generated_proteins_folder = os.path.join(GENERATED_BASE_DIR, f"run{RUN_NUMBER}_filtered")
    tox_folder = TOX_FOLDER
    output_csv_path = os.path.join(OUTPUT_DIR, OUTPUT_FILENAME_TEMPLATE.format(run_number=RUN_NUMBER))

    toxin_sequences_with_id = load_toxin_sequences(tox_folder)
    print(f"Total toxin sequences loaded: {len(toxin_sequences_with_id)}")

    compare_all_sequences_optimized_for_memory(generated_proteins_folder, toxin_sequences_with_id, output_csv_path)

    print("\nProgram finished successfully.")


if __name__ == "__main__":
    main()
