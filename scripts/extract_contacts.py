from Bio.PDB import PDBParser, NeighborSearch
from Bio.SeqUtils import seq1
import pandas as pd
import sys

CUTOFF = 5.0

STRUCTURES = {
    "9NFU": {
        "antibody_chain": "C",
        "antigen_chains": ["A"],
        "mapping_file": "data/processed/9NFU_residue_mapping.csv",
        "output_file": "data/processed/9NFU_cdr_contacts.csv",
    },
    "9NH7": {
        "antibody_chain": "E",
        "antigen_chains": ["B", "H"],
        "cdr_ranges": {
            "H1": range(26, 34),
            "H2": range(51, 59),
            "H3": range(97, 112),
        },
        "output_file": "data/processed/9NH7_cdr_contacts.csv",
    },
}

# Read structure ID from command line
structure_id = "9NFU"

for arg in sys.argv[1:]:
    candidate = arg.upper()
    if candidate in STRUCTURES:
        structure_id = candidate
        break

config = STRUCTURES[structure_id]

# Load PDB
pdb_file = f"data/pdb/{structure_id}.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure(structure_id, pdb_file)
model = structure[0]

antibody = model[config["antibody_chain"]]

# Collect atoms from all antigen chains
antigen_atoms = []

for chain_id in config["antigen_chains"]:
    antigen_atoms.extend(
        list(model[chain_id].get_atoms())
    )

neighbor_search = NeighborSearch(antigen_atoms)

# Find antibody residues within cutoff
contacts = []

for residue in antibody:

    if residue.id[0] != " ":
        continue

    minimum_distance = float("inf")
    closest_antigen_chain = None
    closest_antigen_residue = None

    for atom in residue.get_atoms():

        nearby_atoms = neighbor_search.search(
            atom.coord,
            CUTOFF,
            level="A"
        )

        for antigen_atom in nearby_atoms:

            distance = atom - antigen_atom

            if distance < minimum_distance:
                minimum_distance = distance

                antigen_residue = antigen_atom.get_parent()
                antigen_chain = antigen_residue.get_parent()

                closest_antigen_chain = antigen_chain.id
                closest_antigen_residue = antigen_residue.id[1]

    if minimum_distance < float("inf"):

        contacts.append({
            "antibody_chain": config["antibody_chain"],
            "pdb_residue_number": residue.id[1],
            "amino_acid": seq1(residue.resname),
            "antigen_chain": closest_antigen_chain,
            "closest_antigen_residue": closest_antigen_residue,
            "minimum_distance_A": round(float(minimum_distance), 2),
        })

contacts_df = pd.DataFrame(contacts)

# Assign CDR labels
if structure_id == "9NFU":

    mapping_df = pd.read_csv(
        config["mapping_file"]
    )

    results = mapping_df.merge(
        contacts_df,
        on="pdb_residue_number",
        how="inner"
    )

    cdr_contacts = results[
        results["cdr"] != "framework"
    ].copy()

else:

    def assign_cdr(pdb_number):

        for cdr_name, residue_range in config["cdr_ranges"].items():

            if pdb_number in residue_range:
                return cdr_name

        return "framework"

    contacts_df["cdr"] = (
        contacts_df["pdb_residue_number"]
        .apply(assign_cdr)
    )

    cdr_contacts = contacts_df[
        contacts_df["cdr"] != "framework"
    ].copy()

# Save
cdr_contacts.to_csv(
    config["output_file"],
    index=False
)

# Summary
print("--------------------")
print("CONTACT SUMMARY")
print("--------------------")
print(f"Structure: {structure_id}")
print(f"Antibody chain: {config['antibody_chain']}")
print(f"Antigen chains: {config['antigen_chains']}")
print(f"Working cutoff: {CUTOFF} Å")
print(f"Total antibody contacts: {len(contacts_df)}")
print(f"Total CDR contacts: {len(cdr_contacts)}")

print("\nContacts by CDR:")
print(
    cdr_contacts["cdr"]
    .value_counts()
    .sort_index()
)

print(f"\nSaved to: {config['output_file']}")
