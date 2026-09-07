from Bio.PDB import PDBParser, NeighborSearch

# Load the 9NFU structure
pdb_file = "data/pdb/9NFU.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("9NFU", pdb_file)

# Get first structural model
model = structure[0]

# 9NFU chain mapping
antigen = model["A"]
antibody = model["C"]

# Collect antigen atoms
antigen_atoms = list(antigen.get_atoms())
neighbor_search = NeighborSearch(antigen_atoms)

# Temporary contact-distance cutoff in Angstroms
cutoff = 5.0

contact_residues = []

# Search antibody residues for antigen contacts
for residue in antibody:

    # Ignore water and non-standard residues
    if residue.id[0] != " ":
        continue

    minimum_distance = float("inf")
    is_contact = False

    for atom in residue:
        nearby_atoms = neighbor_search.search(
            atom.coord,
            cutoff
        )

        if nearby_atoms:
            is_contact = True

            for antigen_atom in nearby_atoms:
                distance = atom - antigen_atom
                minimum_distance = min(
                    minimum_distance,
                    distance
                )

    if is_contact:
        contact_residues.append(
            (
                residue.resname,
                residue.id[1],
                minimum_distance
            )
        )

# Display results
print(f"Structure: 9NFU")
print(f"Antigen chain: A")
print(f"Antibody chain: C")
print(f"Contact cutoff: {cutoff} Å")
print(
    f"Number of antibody contact residues: "
    f"{len(contact_residues)}"
)
print()

for residue_name, residue_number, distance in contact_residues:
    print(
        f"Chain C | {residue_name} {residue_number} | "
        f"minimum distance = {distance:.2f} Å"
    )
import pandas as pd

# Convert detected contacts into a dataframe
contacts_df = pd.DataFrame(
    contact_residues,
    columns=[
        "amino_acid_3letter",
        "pdb_residue_number",
        "minimum_distance_A"
    ]
)

# Load sequence/PDB/CDR mapping
mapping_df = pd.read_csv(
    "data/processed/9NFU_residue_mapping.csv"
)

# Match structural contacts to their sequence and CDR positions
results = mapping_df.merge(
    contacts_df,
    on="pdb_residue_number",
    how="inner"
)

# Keep only contacts located within CDRs
cdr_contacts = results[
    results["cdr"] != "framework"
].copy()

# Save results
output_file = "data/processed/9NFU_cdr_contacts.csv"

cdr_contacts.to_csv(
    output_file,
    index=False
)

print("\n--------------------")
print("CDR CONTACT SUMMARY")
print("--------------------")

print(
    cdr_contacts[
        [
            "sequence_position",
            "amino_acid",
            "pdb_residue_number",
            "domain",
            "cdr",
            "minimum_distance_A"
        ]
    ].to_string(index=False)
)

print()
print(f"Total antibody contacts: {len(results)}")
print(f"Total CDR contacts: {len(cdr_contacts)}")
print(f"Saved to: {output_file}")
