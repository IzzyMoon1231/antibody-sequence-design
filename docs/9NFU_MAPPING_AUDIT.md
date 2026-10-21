# 9NFU mapping and contact-method verification (review in progress)

## Verified from repository source

- `scripts/extract_contacts.py` uses `CUTOFF = 5.0` and selects antibody chain C against antigen chain A for 9NFU.
- The algorithm searches **all atoms** of the antigen against antibody-residue atoms, records the closest hit, and selects residues with a distance within the 5.0 Å neighborhood. This is an operational working cutoff, **not a verified Bennett et al. cutoff**.
- `scripts/map_cdrs.py` builds chain-C sequential positions from the PDB and finds six hard-coded CDR peptide substrings, then writes `data/processed/9NFU_residue_mapping.csv`.
- `data/processed/9NFU_cdr_contacts.csv` contains 14 historical CDR-contact rows (H1/H2/H3/L1/L2); these have **not** been independently regenerated.
- `results/9NFU/9NFU_IMGT_residue_mapping.csv` is an existing historical mapping output, but it is **not** the missing `data/processed/9NFU_residue_mapping.csv`.

## Blocking issues before treating mapping as authoritative

1. **Unverified numbering:** `map_cdrs.py` uses peptide-substring search rather than a saved residue-by-residue IMGT assignment with provenance. Verify each loop boundary against a pinned AbNumber/ANARCII version and numbering output.
2. **Hard-coded domain split:** `sequence_position <= 95` is labeled heavy and the rest light; this may classify linker residues as light and must not be treated as validated domain annotation.
3. **Non-unique join:** `extract_contacts.py` merges only on integer `pdb_residue_number`, ignoring insertion code and chain identity. Preserve model, chain, hetero flag, PDB residue number and insertion code, and enforce a unique mapping key.
4. **Missing file:** The contact extractor references `data/processed/9NFU_residue_mapping.csv`, absent on main. Do not invent mapping rows from saved contact tables.
5. **Methodology not reconciled:** Confirm whether the paper's definition uses all-heavy-atom contacts, hydrogen exclusion, minimum interatomic distance, and the same cutoff; sensitivity-test reasonable cutoffs once authoritative protocol is located.
6. **Original predictions not traced:** The archived PDBs and finalized FASTA input provenance have not been established as the origin of saved structural metrics. No paper claim should imply that these metrics were independently regenerated.

## Reproduction checklist

- [ ] Parse the actual deposited 9NFU PDB, selecting model and chain explicitly; record source SHA256 and insertion codes.
- [ ] Reconstruct chain-C amino-acid sequence and compare it to the sequence assumed by the CDR numbering source.
- [ ] Generate per-residue mapping with unique chain/model/number/insertion-code keys and auditable numbering source.
- [ ] Check all six CDR boundaries and scFv linker explicitly; fail on ambiguous substring matches.
- [ ] Rerun contact extraction and compare the resulting 14 historical CDR-contact rows, recording any differences.
- [ ] Reconcile contact threshold with published methods; do not label it authoritative prematurely.
- [ ] Preserve the historical outputs and document any new runs separately.

## Pull-request review order

1. PR #1: serializer and storage tests; CPU synthetic fixture is not biological inference validation.
2. PR #2: provenance manifest, input FASTA/source hashes, reference/map authority; unresolved status must remain unresolved.
3. PR #3: eligibility checks; contract validation must not be described as statistical or biological validation.

**Status:** This file is an audit checklist, not a validated residue mapping or a completed scientific reproduction.
