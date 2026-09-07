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

# CDR sequences identified with AbNumber/ANARCII
cdr_sequences = {
    "H1": "GFSIKNTY",
    "H2": "IWPANGKT",
    "H3": "SRQLDPYNLYGNDV",
    "L1": "SSSY",
    "L2": "RNS",
    "L3": "STMNNDGNLV"
}

# Find the CDR ranges in the full sequence
cdr_ranges = {}

for cdr_name, cdr_seq in cdr_sequences.items():
    start_index = sequence.find(cdr_seq)

    if start_index == -1:
        print(f"WARNING: {cdr_name} not found")
        continue

    start_pos = start_index + 1
    end_pos = start_pos + len(cdr_seq) - 1

    cdr_ranges[cdr_name] = (start_pos, end_pos)

# Label every residue
for row in mapping:
    row["cdr"] = "framework"

    pos = row["sequence_position"]

    for cdr_name, (start_pos, end_pos) in cdr_ranges.items():
        if start_pos <= pos <= end_pos:
            row["cdr"] = cdr_name
            break

# Print basic information
print("9NFU Chain C sequence:")
print(sequence)
print("\nSequence length:", len(sequence))

print("\nCDR ranges:")
for cdr_name, (start_pos, end_pos) in cdr_ranges.items():
    print(f"{cdr_name}: sequence positions {start_pos}-{end_pos}")

# Save mapping
output_file = "data/processed/9NFU_residue_mapping.csv"

with open(output_file, "w", newline="") as csvfile:

    fieldnames = [
        "sequence_position",
        "amino_acid",
        "pdb_residue_number",
        "domain",
        "cdr"
    ]

    writer = csv.DictWriter(
        csvfile,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(mapping)

print(f"\nMapping saved to: {output_file}")
