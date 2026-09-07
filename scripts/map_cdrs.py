from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

# Load structure
pdb_file = "data/pdb/9NFU.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("9NFU", pdb_file)

# Antibody/scFv is Chain C
chain = structure[0]["C"]

sequence = ""
pdb_residue_numbers = []

for residue in chain:
    if residue.id[0] != " ":
        continue

    try:
        aa = seq1(residue.resname)
    except:
        aa = "X"

    sequence += aa
    pdb_residue_numbers.append(residue.id[1])

print("9NFU Chain C sequence:")
print(sequence)
print()
print("Sequence length:", len(sequence))

print("\nFirst 20 sequence-to-PDB mappings:")
for i in range(min(20, len(sequence))):
    print(
        f"Sequence position {i+1}: "
        f"{sequence[i]} -> PDB residue {pdb_residue_numbers[i]}"
    )
