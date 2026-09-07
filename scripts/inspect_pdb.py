from Bio.PDB import PDBParser

# Location of our test structure
pdb_file = "data/pdb/9NFU.pdb"

# Read the PDB structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("9NFU", pdb_file)

print("Structure: 9NFU")
print("--------------------")

# Find every chain in the structure
for model in structure:
    for chain in model:
        residues = [
            residue
            for residue in chain
            if residue.id[0] == " "
        ]

        print(
            f"Chain {chain.id}: "
            f"{len(residues)} amino-acid residues"
        )
