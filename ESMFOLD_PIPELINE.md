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

## P1b input and provenance contract (CPU validation only)

The v1 bundle extends this repository's existing prediction manifest and validator.
It does not approve biological claims or regenerate sequences, predictions, metrics,
or statistics. The four final input files contain **400 sampled design records** and
**363 distinct full sequences**. Repeated sequences retain separate design identities.
A sequence check, a matching run receipt, and scientific metric eligibility are
separate checks. Original generation-code revision remains `NOT_ESTABLISHED`.

Build a deterministic bundle into a new directory outside both source repositories:

```sh
python scripts/build_prediction_manifest.py \
  --repo-root /path/to/antibody-sequence-design \
  --upstream-root /path/to/ProteinMPNN-Ablang-antibody-design \
  --output-root /path/to/new-bundle
```

Inputs must match their named local Git commits and expected origin repositories.
The builder reads existing files; it does not execute upstream scripts, fill residues,
fetch sources, import model code, or calculate metrics. Its four outputs are:

- `prediction_manifest.csv`: existing `target,method,index,sequence,pdb_path`
  columns plus source/header/ordinal, file/sequence hashes, snapshot identities,
  map/reference hashes, design identity, and optional planned run/prediction identity.
- `canonical_residue_map.csv`: 362 reference coordinate positions, preserving
  missing atoms, SEQRES/ATOM differences and current/historical masks and domains.
- `input_provenance.json`: source hashes, unknown generation revision, run status
  and input-bundle digest; no invented scientific authority.
- `prediction_eligibility.csv`: input evidence with prediction state
  `NOT_INSPECTED` and final metric eligibility `NOT_ASSESSED`.

Existing `results/file_manifest.csv`, historical predictions and metric tables are
unchanged. Two builds from identical inputs produce identical identity artifact
bytes. Manifest row permutations do not change normalized records or input digest;
source FASTA order remains part of source-record provenance.

### Identity and path rules

`source_sequence_sha256` identifies exact uppercase ASCII sequence bytes without a
newline. `record_key` identifies a target/method/1-based FASTA ordinal within an input
cohort. `design_id` additionally binds that key to source repository, snapshot,
path and exact file hash; a numeric sample or seed is not a cross-method pairing.
`prediction_id` binds a design, planned run hash and `prediction_replicate_id`.
It is empty for an unresolved draft run. A prediction identity is a plan until a
matching output receipt is observed. Measurement IDs are separate from all of these.

One bundle contains one requested prediction replicate for each of 400 designs.
Use another bundle with `--prediction-replicate-id 2` for a separately identified
replicate; do not collapse independent folds or sampled designs by sequence hash.
No prediction cache/deduplication is implemented. Changing replicate ID changes the
input/run plan, while source design identity remains unchanged.

All v1 paths are repository-relative or relative to an explicit prediction output
root. Absolute paths, traversal, symlink escapes and destination collisions fail.
Predictions use `<run-hash>/<target>/<method>/<index>/replicate-<id>.pdb`.
Draft paths use `draft/` and cannot execute. Neither filenames nor repeated run-hash
references alone establish a measurement join. Repeated references to a single run
are legitimate; duplicate/conflicting run registry definitions are rejected.

### Coordinate and scientific authority

The map records physical chains (9NFU C, 9NH7 E), sequence/SEQRES indices,
PDB residue numbers and insertion codes, CA presence and source-qualified domains.
For the frozen inputs, 9NFU sequence index is PDB number+2 and SEQRES index is
sequence index+16; 9NH7 uses PDB number−2 and sequence index+4. Unresolved positions
remain explicit. A new insertion code requires an explicit coordinate mapping;
it is not silently compressed away. The parser records the actual first MODEL
serial and accepts a single CA with blank/A alternate-location code; duplicate CA
keys or other ambiguous alternate locations fail. Strict predicted sequences must
be one chain, contiguous numbering 1..length, no insertion codes, with CA per residue.

The 9NFU chain C residue −1/design index 1 discrepancy remains Q in ATOM versus
H in SEQRES. Current and historical masks/domains remain separate columns, with
`UNRESOLVED_FOR_FINAL_METRICS`. Current code is not proof of the final inputs'
generation version. The validator reconstructs map semantics from committed source
facts, not just map hashes. Full-length mismatch alone is not a verdict about the
validity of a partial-domain historical metric.

### Strict validation and shared API

```sh
python scripts/validate_esmfold_predictions.py \
  --manifest /path/to/new-bundle/prediction_manifest.csv \
  --input-root /path/to/antibody-sequence-design \
  --upstream-root /path/to/ProteinMPNN-Ablang-antibody-design \
  --output-root /path/to/predictions \
  --provenance-mode strict --report /path/to/new-report.csv
```

A report path must be new. This command validates source/map identity and reports
prediction availability; it creates no predictions. Explicit `--provenance-mode
legacy` retains historical validator path/sequence behavior without v1 assurances.
Existing Python sequence-validation functions remain available. The runner also
permits legacy restart inspection for compatibility, but legacy execution cannot
load a model under this contract.

P1c should consume `validate_prediction_bundle(...)`, which returns canonical
records, freshly checked eligibility, `generation_code_commit=None`,
`generation_code_status=NOT_ESTABLISHED` and unresolved mapping authority. Do not
construct eligibility from user-supplied PASS/status columns. The shared
`validate_measurement_joins` checks prediction-backed identity/cardinality only:
`record_key,design_id,prediction_id,source_sequence_sha256,output_pdb_sha256,
reference_sha256,residue_map_sha256,run_contract_sha256`, plus measurement,
metric-policy and extraction identifiers. It rejects duplicate/conflicting or
unverified prediction joins. It does not calculate metrics or approve an estimand,
replicate unit, pairing or multiplicity policy. Sequence-only score provenance is
not implemented by this prediction-backed helper; P1c must specify it separately.

### Future execution contract and receipts

P1b tests never perform inference. A future separately authorized pilot must supply
an optional `--run-contract` JSON with schema `uirj-run-v1` and exactly these fields:
`schema_version,model_id,model_revision,assets,serializer_sha256,runtime,device,
precision,inference_settings,seed,deterministic_algorithms,input_bundle_digest,
run_instance_id`. No `approved:true` field is interpreted as authority.

Supported contract values are model `facebook/esmfold_v1`, a 40-character declared
revision, an `assets` mapping of every local regular file to its SHA256,
`runtime` containing explicit Python/torch/transformers versions, `device=cuda`,
`precision={"esm":"float16","trunk":"float32"}`, and positive integer
`num_recycles`/`chunk_size` inference settings. Unsupported or incomplete settings
remain unresolved draft evidence when building; they cannot execute.
`serializer_sha256` binds the exact runner file including the CPU serializer.
`input_bundle_digest` excludes run-dependent output names to avoid circular IDs.
`run_instance_id` distinguishes separately commissioned runs with the same settings.

Local asset hashes pin bytes, not their authenticated upstream origin. A revision
string passed to a local loader does not establish where those bytes came from.
All local files must match the declared inventory; symlinks or extra files fail.
Both tokenizer/model calls use `local_files_only=True`; there is no download fallback.
Actual runtime/serializer identity is checked before execution and again before
receipt persistence. Floating model parameter precision and CUDA placement are
observed after model preparation. Runtime compatibility and inference correctness
remain untested until a real pilot; mocks establish only parameter binding.

The base seed is transformed per prediction as
`(base_seed + int(SHA256(prediction_id ASCII)[:8],16)) % 2**32`.
The runner applies that seed and the explicit deterministic-algorithms setting
before each prediction; row order does not change the requested per-prediction seed.
The exact runner digest binds this policy. It does not promise cross-device bitwise
reproducibility without pilot evidence.

The runner additionally requires `--execute-approved-run CONTRACT_SHA256` and an
explicit `--asset-root` before a new prediction. This argument is a workflow
acknowledgment, not cryptographic proof of owner authorization. It is not permission
from this document or from JSON. A dry run loads no model and writes no outputs:

```sh
python scripts/run_esmfold_predictions.py \
  --manifest /path/to/new-bundle/prediction_manifest.csv \
  --input-root /path/to/antibody-sequence-design \
  --upstream-root /path/to/ProteinMPNN-Ablang-antibody-design \
  --output-root /path/to/predictions --limit 1 --dry-run
```

Persistence remains exclusive-lock protected and never overwrites existing files:
validated PDB bytes are saved/fsynced, `<pdb>.provenance.json` is saved/fsynced, then
`.complete` is published. Interrupted files/receipts remain for inspection. Legacy
COMPLETE-only files can pass sequence checks but remain `RUN_UNVERIFIED`; they are
not silently accepted as v1 completed predictions. Receipts bind input/run/design/
replicate/output hashes and observed software identity, not proof of real inference,
correct atomic geometry, model accuracy or scientific validity. Final metric
eligibility remains `NOT_ASSESSED`, even for a provenance-verified prediction.

### CPU acceptance commands

Run within the externally monitored resource budget; do not run model inference:

```sh
python -m pytest -q tests/test_prediction_provenance.py \
  tests/test_esmfold_prediction_storage.py tests/test_esmfold_serialization.py
```

Fixtures cover input/hash/path mutations, duplicated sampled sequences, prediction
replicates, first-model/insertion/CA ambiguity, unresolved map facts, stale receipts,
interruptions, run/measurement join ambiguity, local asset inventory, observed
parameter settings, and model-free dry runs. Acceptance also includes two actual
bundle CLI runs, strict validator CLI, source-preservation checks and deliberate
rehashed map/source/run mutations. Report executed outcomes separately from this
checklist; its existence is not a passing result.
