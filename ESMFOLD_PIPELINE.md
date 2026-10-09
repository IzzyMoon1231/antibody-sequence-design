# ESMFold Prediction Pipeline

## Purpose
Generate ESMFold structures for finalized antibody sequences.

## Dataset
- 9NFU baseline: 100 sequences, 247 residues
- 9NFU ensemble: 100 sequences, 247 residues
- 9NH7 baseline: 100 sequences, 115 residues
- 9NH7 ensemble: 100 sequences, 115 residues

Total: 400 finalized sequences.

## Scripts
- scripts/run_esmfold_predictions.py
  - Reads prediction manifest
  - Supports dry runs and limited batches
  - Skips validated existing predictions
  - Validates generated PDB sequences
  - Uses exclusive lock files during saving

- scripts/validate_esmfold_predictions.py
  - Extracts amino acid sequences from PDBs
  - Compares PDB sequences against the manifest
  - Produces CSV validation reports

## Example dry run

Run the prediction script with these arguments:

--manifest /path/to/prediction_manifest.csv
--target 9NH7
--method baseline
--limit 1
--dry-run

## Validation status
- Finalized FASTA validation: PASS
- PDB validator unit tests: PASS
- Restart and file-protection tests: PASS
- Synthetic PDB conversion test: PASS
- Real ESMFold GPU inference: NOT TESTED
- Finalized ESMFold predictions: NOT GENERATED

## Limitations
Original baseline/ensemble predicted PDBs were not recovered.
Existing AbMPNN PDBs cannot substitute for finalized sequences.
Real GPU inference and structural validation remain necessary.
Old lock files require manual review before removal.

## Environment
See results/summary/esmfold_development_environment.json.
