# 2026-10-19 structural evaluation handoff

## Status and immutable source references

Collaborator repository: https://github.com/sanithuheen/ProteinMPNN-Ablang-antibody-design

- RMSD CDR CSV Git blob: `33182e23e4f8242d46cca96144a3a0aeeb9807d2`
- RMSD framework CSV Git blob: `c0f2014940ce96ec4ab1175cf0fccf7cb8900835`
- Baseline finalized 9NFU FASTA at inspection: `e0027963a950dfc397ec4497a49c7fdb5b3760a9`
- Baseline finalized 9NH7 FASTA at inspection: `872fb3160ef966125081d4b304caf8a9e957b6d4`

Collaborator reports the baseline was generated from an unintended subset of the original PDB and is regenerating the baseline FASTAs. **The above baseline files and any predictions or metrics generated from them are provisional.** Do not overwrite or silently mix baseline generations.

## Results inspection

The current CDR CSV has 2,700 rows, matching 100 designs per method per target per loop: 9NH7 (3 methods × 3 loops × 100) and 9NFU (3 methods × 6 loops × 100).

Selected current arithmetic means (Å; **not** updated baseline claims):

| Target | Loop | Baseline | Ensemble | AbMPNN |
| --- | --- | ---: | ---: | ---: |
| 9NFU | H1 | 1.265 | 1.092 | 1.092 |
| 9NFU | H2 | 1.420 | 1.857 | 1.436 |
| 9NFU | H3 | 2.953 | 2.983 | 1.923 |
| 9NFU | L1 | 1.316 | 1.313 | 1.268 |
| 9NFU | L2 | 1.893 | 2.335 | 2.133 |
| 9NFU | L3 | 0.974 | 1.700 | 1.263 |
| 9NH7 | H1 | 2.578 | 2.626 | 2.587 |
| 9NH7 | H2 | 0.801 | 1.296 | 1.110 |
| 9NH7 | H3 | 0.906 | 2.916 | 4.030 |

Five 9NH7 AbMPNN designs carry `framework_outlier=1` in each loop (15 flagged CDR rows, not 15 independent designs). These are not to be silently excluded. Inspect the framework QC records and report sensitivity analysis with and without flagged designs, with exclusions clearly labeled.

## Next mandatory checks

1. After corrected baseline FASTAs land, record new blob hashes, sequence counts, and sequence-level differences from the above hashes.
2. Run `scripts/audit_esmfold_sequence_provenance.py` against the collaborator checkout; preserve the CSV. The existence of 600 PDBs is not proof that each is a fold of the correct final sequence.
3. If a baseline sequence changes, **do not reuse its old PDB or RMSD/pLDDT measurement**. Regenerate the affected ESMFold prediction and all dependent measurements. If no sequence changes, verify PDB correspondence before reuse.
4. Recompute both RMSD CSVs using the same pinned `compute_RSMD.py` revision; capture source hash and reference PDB hash.
5. Confirm whether the manuscript's target is experimental-reference CDR RMSD or RFdiffusion-backbone CDR RMSD; these are different estimands and must not be conflated.
6. Reconcile 9NFU domain boundaries and mapping, and report the 9NH7 AbMPNN outlier sensitivity.
7. Only after inputs and measurement provenance pass, rerun statistics and update manuscript text/figures. Preserve older outputs as historical, not current.

This handoff does **not** assert the corrected baseline has arrived or that GPU folding was rerun.
