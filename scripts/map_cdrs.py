from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import csv

# Load structure
pdb_file = "data/pdb/9NFU.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("9NFU", pdb_file)

# Antibody/scFv chain
chain = structure[0]["C"]

sequence = ""
mapping = []

# Build sequence and residue mapping
for residue in chain:
    if residue.id[0] != " ":
        continue

    try:
        aa = seq1(residue.resname)
    except:
        aa = "X"

    sequence += aa

    mapping.append({
        "sequence_position": len(sequence),
        "amino_acid": aa,
        "pdb_residue_number": residue.id[1]
    })

# Temporary heavy/light split
for row in mapping:
    pos = row["sequence_position"]

    if pos <= 95:
        row["domain"] = "heavy"
    else:
        row["domain"] = "light"

print("9NFU Chain C sequence:")
print(sequence)
print("\nSequence length:", len(sequence))

output_file = "data/processed/9NFU_residue_mapping.csv"

with open(output_file, "w", newline="") as csvfile:

    fieldnames = [
        "sequence_position",
        "amino_acid",
        "pdb_residue_number",
        "domain"
    ]

    writer = csv.DictWriter(
        csvfile,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(mapping)

print(f"\nMapping saved to: {output_file}")
