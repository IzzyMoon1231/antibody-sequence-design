# Baseline ESMFold sequence spot-check — 2026-10-20

Source: `sanithuheen/ProteinMPNN-Ablang-antibody-design` default branch.

Compared the sequence of unique C-alpha ATOM residues (chain + PDB author residue ID + insertion code) in each archived PDB to the FASTA record at the corresponding zero-based index. Residues were translated from PDB three-letter amino-acid codes to one-letter codes.

- **9NFU**: indices 000, 001, 025, 050, 075, 099; each PDB exactly matched its corresponding finalized FASTA sequence, length 247 amino acids.
- **9NH7**: indices 000, 001, 025, 050, 075, 099; each PDB exactly matched its corresponding finalized FASTA sequence, length 115 amino acids.
- **12/12 sampled comparisons passed.** This is **not** a full audit of all 200 baseline designs.
- The source FASTA files each contain 100 records.
- FASTA blob hashes: 9NFU `e0027963a950dfc397ec4497a49c7fdb5b3760a9`; 9NH7 `872fb3160ef966125081d4b304caf8a9e957b6d4`.

To verify all designs, run `python scripts/audit_esmfold_sequence_provenance.py /path/to/collaborator/repo --output esmfold_sequence_audit.csv` from the user repository checkout. This checks all six method/target cohorts and emits one row per design plus missing/extra files. Review the generated CSV and keep it as an artifact.

**Not established by this check:** whether the baseline designs were generated from the intended original PDB input, whether residue mapping/CDR assignments are biologically correct, whether the predicted structures are reliable, whether Edward's independent baseline validator passes, or whether OASIS/ablation reference inputs are correct. Those are separate checks.
