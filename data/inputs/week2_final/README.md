# Week 2 finalized sequence inputs

Source: https://github.com/sanithuheen/ProteinMPNN-Ablang-antibody-design (main branch, retrieved 2026-10-07)

These four FASTA files are exact copies of the teammate's **final_filled** baseline and ensemble FASTA files, not the similarly named non-filled versions.

| Local path | Original source path | Designs | Residues per design |
|---|---|---:|---:|
| baseline/9NFU.fa | baseline_sequences/baseline_sequences_final_filled/9NFU.fa | 100 | 247 |
| baseline/9NH7_EBH.fa | baseline_sequences/baseline_sequences_final_filled/9NH7_EBH.fa | 100 | 115 |
| ensemble/9NFU.fa | ensemble_sequences/ensemble_sequences_final_filled/9NFU.fa | 100 | 247 |
| ensemble/9NH7.fa | ensemble_sequences/ensemble_sequences_final_filled/9NH7.fa | 100 | 115 |

All 400 FASTA entries have no X placeholders. The finalized FASTA files have been audited, including a complete 0-based residue-mapping verification for 9NFU (45,000/45,000 matches). Downstream statistical calculations have been reproduced from saved measurements. The original ESMFold prediction files are not archived, so it has not been established that all saved structural measurements were generated from these exact finalized FASTA entries. Do not represent old contact-geometry preservation metrics as CDR loop RMSD. CDR numbering/position correspondence and the AbMPNN input must be verified before structural reruns.

Next: compare old prediction input sequences against these exact FASTA entries, rerun ESMFold where they differ, align structures appropriately, and compute per-loop RMSD relative to the RFdiffusion reference structure for all three methods.
