# Input format guide

The student-facing input and update instructions are in [README.md](README.md#updater-route--no-json-editing-required). Edit `input_selection.csv` in a spreadsheet editor, not analysis internals.

Keep a compatible template row for each format. FASTAs preserve sample/seed numbers; Jacob designs use `full_seq_filled`. BioPhi workbook reads only Overview, never peptide-level rows. Select actual measurement columns and units explicitly. Percentile is not Identity; promb is not BioPhi; framework RMSD is not a loop RMSD. Empty/NA values stay missing. No sequence is filled or ID shifted automatically.

For evidence repair, supply original scored/folded sequence files and manifests instead of copying desired FASTA hashes into a score table. The notebook deliberately does not turn these saved exports into authenticated sequence measurements.
