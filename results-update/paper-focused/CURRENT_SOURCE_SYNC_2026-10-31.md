# Latest GitHub source alignment check — 2026-10-31

The checked-in paper-focused notebook **already includes the currently published corrected baseline inputs and RMSD results**. The source blobs in the collaborator's default branch match the files in this notebook branch byte-for-byte by Git blob SHA:

| Source | Notebook input | Git blob SHA | Match |
|---|---|---|---|
| `baseline_sequences/baseline_sequences_final_filled/9NFU.fa` | `inputs/9NFU_baseline.fa` | `6f75dc71f6cd43c042e8d6a5ef07311f215b8923` | Yes |
| `baseline_sequences/baseline_sequences_final_filled/9NH7_EBH.fa` | `inputs/9NH7_baseline.fa` | `fb76c39f574482dcb0c3c3ef8da2be8f9c33a681` | Yes |
| `RSMD_RESULTS/cdr_rmsd_results.csv` | `inputs/cdr_rmsd_results.csv` | `5aeced4700ceb76bfff1ea3d0ba43bd06521dadb` | Yes |

Source collaborator: https://github.com/sanithuheen/ProteinMPNN-Ablang-antibody-design

The framework QC CSV `RSMD_RESULTS/cdr_rmsd_framework.csv` currently has Git blob SHA `eaffe068aecd4c6949da9ba167eb49c821f67afb`. It is **not** among the notebook's current selected input files, so this sync check does not certify framework-outlier sensitivity.

**Conclusion:** The notebook is not stale with respect to these three current source files. Do **not** overwrite figures or rerun the notebook merely because these files changed relative to an older audit; the existing notebook package already carries the updated versions. See `inputs/provenance.json` and `HANDOFF_UPDATE.md` for the pinned input lineage and scope.

**Remaining checks before final paper:** verify source ESMFold PDB-to-FASTA correspondence and framework-outlier sensitivity, ensure any new Oasis workbooks correspond to the corrected sequences, and confirm the manuscript uses the updated figures and tables. This file is a Git blob-level synchronization receipt, **not** a claim of full end-to-end scientific validation.
