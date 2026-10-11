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
The original unqualified PASS entries were documentation-reported checks; their
executable tests were not present in the reviewed base revision. The following
results are newly observed on the bounded CPU serialization fixture (2026-10-10):

- Original conversion at `3be2136e0863b705aee42715a681600b64530ae4`: **FAIL**, pytest exit 1.
  The realistic `GAVW` fixture passes CA sequence validation but emits 24 atoms
  instead of the official export's 30. Atom identity/numbering assertion fails;
  this is an observed defect, not an expected-failure/xfail marker.
- Corrected conversion and temporary-directory storage regression suite:
  **41 passed**, pytest exit 0 (1.09 seconds in the recorded run).
- Real ESMFold inference on any device: **NOT TESTED** by this patch.
- Finalized ESMFold predictions: **NOT GENERATED** by this patch.
- Existing scientific measurements and biological claims: **NOT VALIDATED** by
  synthetic export or storage tests.

## Reproduce the CPU contract tests
Run from this repository with existing dependencies; no weights or model imports
are needed. Do not run the production prediction CLI for these tests.

```sh
# Intentionally RED: execute the actual old conversion block extracted by AST
# from the pinned Git revision. Requires that revision in local Git history.
PYTHONDONTWRITEBYTECODE=1 HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 CUDA_VISIBLE_DEVICES='' UIRJ_SERIALIZER_BASELINE=1 python3 -m pytest -p no:cacheprovider -q tests/test_esmfold_serialization.py::test_realistic_atom14_export_matches_official_oracle

# GREEN: corrected serializer plus unchanged storage and restart protections.
PYTHONDONTWRITEBYTECODE=1 HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 CUDA_VISIBLE_DEVICES='' python3 -m pytest -p no:cacheprovider -q tests/test_esmfold_serialization.py tests/test_esmfold_prediction_storage.py
```

Recorded environment: Python 3.14.7, torch 2.12.0, Transformers 5.8.1,
NumPy 2.4.6, pytest 9.0.3. This differs from the repository's development
metadata. Compatibility is established only for the exercised export contract:
[official Transformers v5.18.0 export](https://github.com/huggingface/transformers/blob/a906d3c4b65095f2308b6a6a193e934d03b8eb5d/src/transformers/models/esm/modeling_esmfold.py),
commit `a906d3c4b65095f2308b6a6a193e934d03b8eb5d`. Its export function is
AST-identical to the installed function; `feats.py`, `tensor_utils.py`,
`data_transforms.py` and `protein.py` are byte-identical. The tests enforce the
pinned export AST and helper hashes, failing on unreviewed dependency drift.
Residue constants were separately checked: the source differs only in comment
spelling. This is not whole-package or model equivalence.

The independent oracle extracts and executes the official export function,
without importing the modeling module. It never invokes the candidate seam.
By default it uses the hash-verified installed source; to use an already local
pinned source, set `UIRJ_ESMFOLD_ORACLE` to its `modeling_esmfold.py` path. The checked-in tests use the installed, hash-verified source by default; no
machine-specific oracle path is required. The commands above reproduce the
original failure and corrected behavior when the pinned Git revision and
verified dependency sources are available. For handoff, retain stdout, exit
codes and dependency versions alongside each test run; do not write test
receipts into scientific prediction directories.

Fixture positions have explicit `(structure stages, batch, residues, atom14,
xyz)` axes; the structure-stage dimension is distinct from outer recycling.
All fixture tensors are explicitly CPU tensors. `GAVW` covers real residue
mappings with unique coordinate sentinels, nonuniform atom-level pLDDT,
nonconsecutive residue numbering, a missing side-chain atom and masked padding
in a second batch member. Tests compare all atom identities, coordinates and
B-factors to the official oracle and check selected sentinels directly. A
sequence-correct but atom-coordinate-wrong control is rejected. These synthetic
coordinates are deliberately not physically plausible structures.

`serialize_esmfold_output` maps final atom14 coordinates through the official
`atom14_to_atom37` helper before PDB export. It rejects missing fields, shape,
index dtype, residue-range, residue-inconsistent mapping, invalid mask and
nonfinite-data errors. Missing atoms/padding are allowed as mask subsets. The
existing single-sequence `predict` path delegates only serialization to this
seam; model loading, CUDA settings, validator and storage behavior are unchanged.
Storage tests cover exclusive locking, a competing writer, no overwrite,
completion markers, interrupted/corrupt copies, temp cleanup, read-only lock
inspection and restart skip/block behavior before any model/device access.

## Limitations
The originating predictions for this repository's saved metric tables have not
been traceably recovered. Archived baseline/ensemble PDBs exist in the linked
ProteinMPNN-Ablang repository, but are not established as finalized-input
predictions. Their historical domain-level measurements require separate
sequence, numbering, reference and metric provenance checks; a full-length
mismatch alone does not establish that those historical measurements are invalid.
Existing AbMPNN PDBs cannot substitute for finalized sequences.
Real GPU inference and structural validation remain necessary.
Old lock files require manual review before removal.

## Environment
See results/summary/esmfold_development_environment.json.
