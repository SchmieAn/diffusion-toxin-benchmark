"""
Download toxin PDB files from the RCSB and extract Cα atoms for each structure.

For each PDB code in tox_list, the file is downloaded to save_toxdir, parsed,
and all Cα atoms are extracted (with disordered atoms resolved to their default position).

Usage:
    python protein_compare.py
"""

# ---------------------------------------------------------------------------
# Configuration — edit these values to match your setup.
# ---------------------------------------------------------------------------

# Directory where downloaded toxin PDB files are saved.
SAVE_TOXDIR = "/path_to_your_base_directory/proteins"

# ---------------------------------------------------------------------------


from Bio.PDB import PDBList, PDBParser

# 7QFC excluded: no PDB file available
tox_list = [
    "1FCU", "2MLT", "3PLC", "4OM5", "2H8U", "3VTS", "1KBA", "6RC7", "4LFT", "1UG4",
    "1VYC", "2H5F", "1DTK", "4K09", "1ANS", "1VIB", "1QDP", "2I61", "2ASC", "2ATB"
]

pdbl = PDBList()
parser = PDBParser(QUIET=True)

for code in tox_list:
    pdbl.retrieve_pdb_file(code, pdir=SAVE_TOXDIR, file_format="pdb")

    structure_path = f"{SAVE_TOXDIR}/pdb{code.lower()}.ent"
    structure = parser.get_structure(code, structure_path)

    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    atom = residue["CA"]
                    # Disordered atoms: select the default (most occupied) position
                    if atom.is_disordered():
                        atom = atom.selected_child
                    ca_atoms.append(atom)

    print(f"{code}: {len(ca_atoms)} Cα atoms")
