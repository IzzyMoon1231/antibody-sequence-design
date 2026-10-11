# Reproducibility Guide

## Project Overview

This project compares two computational antibody sequence-design methods:

- **Baseline:** ProteinMPNN
- **Ensemble:** ProteinMPNN + AbLang

Two antibody-containing experimental structures were evaluated:

- **9NFU:** Single-chain variable fragment (scFv)
- **9NH7:** Single-domain antibody (VHH)

The evaluation examines sequence conservation, contact-position
preservation, predicted structural confidence (pLDDT), and
structural agreement with experimental references (RMSD).

## Reproducibility Status

The repository contains saved per-design measurements and verified
scripts for reproducing downstream statistical comparisons.

All four downstream reproducibility checks passed.

This does not mean the complete antibody generation and
structure-prediction pipeline can currently be reproduced from scratch.

## Requirements

Python packages:

- numpy
- pandas
- scipy
- statsmodels

Install with:

    pip install numpy pandas scipy statsmodels

## Running the Reproducibility Checks

From the repository root, run:

    python scripts/run_reproducibility_checks.py --mode legacy-arithmetic

This runs four scripts:

1. scripts/reproduce_contact_statistics.py
2. scripts/reproduce_sequence_conservation.py
3. scripts/reproduce_plddt_statistics.py
4. scripts/reproduce_rmsd_statistics.py

Successful execution ends with:

    SUCCESS: ALL DOWNSTREAM REPRODUCIBILITY CHECKS PASSED

## Dataset Sizes

| Structure | Baseline | Ensemble | Unique Baseline | Unique Ensemble |
|---|---:|---:|---:|---:|
| 9NFU | 100 | 100 | 100 | 100 |
| 9NH7 | 100 | 100 | 87 | 76 |

The 9NH7 pLDDT and RMSD comparisons use unique sequences.

Other analyses use the design counts specified in their
individual verification scripts.

## Statistical Methods

Baseline and ensemble measurements were compared using
two-sided Mann-Whitney U tests.

Where applicable, p-values were adjusted using the
Benjamini-Hochberg false discovery rate (FDR) procedure.

Mean difference was defined as:

    Ensemble mean - Baseline mean

The direction of a desirable difference depends on the metric.

Generated designs are computational observations, not
independent experimental biological replicates.

## Evaluation Metrics

### Sequence Conservation

Measures retention of reference amino acids at evaluated
sequence positions, according to the saved analysis data.

### Contact-Position Preservation

Uses saved mapped-residue C-alpha displacement measurements,
including a 3-angstrom threshold.

This is a geometric proxy and does not directly establish
preservation of antibody-antigen contacts or binding affinity.

### pLDDT

Measures predicted local structural confidence.

Higher pLDDT reflects greater model confidence, not
necessarily better antibody binding or biological function.

### RMSD

Measures structural deviation from experimental references
under the specified structural alignment.

Lower RMSD indicates closer agreement under that alignment,
not necessarily improved binding.

For 9NFU, the analysis includes domain-specific alignment
and a filtered heavy-framework comparison.

## Selected Results

### Contact-Position Preservation (3 Angstroms)

| Structure | Baseline | Ensemble |
|---|---:|---:|
| 9NFU | 75.21% | 68.86% |
| 9NH7 | 81.87% | 53.20% |

### Sequence Conservation

| Structure | Baseline | Ensemble |
|---|---:|---:|
| 9NFU | 65.50% | 43.86% |
| 9NH7 | 15.27% | 22.00% |

### 9NH7 H3 Structural Metrics

| Metric | Baseline | Ensemble |
|---|---:|---:|
| H3 pLDDT | 84.27 | 73.46 |
| H3 RMSD | 0.918 Angstroms | 2.995 Angstroms |

The ensemble showed mixed effects across structures
and CDR regions rather than consistent improvement.

## Limitations

1. Original predicted structure files used to calculate
   pLDDT and RMSD are not currently included in the repository.

2. The complete upstream structure-prediction and RMSD
   extraction workflows are not currently archived.

3. Statistical scripts reproduce calculations from saved
   measurements rather than regenerating the measurements.

4. Geometric preservation does not establish experimental
   antibody-antigen binding or binding affinity.

5. Predicted structural confidence does not demonstrate
   experimental stability or therapeutic effectiveness.

6. Only two structural targets were evaluated.

7. Full regeneration of ensemble sequences may require
   access to the original AbLang model weights.


## Sequence and Structural Provenance Audit

- Both local reference PDB files (9NFU and 9NH7) were verified
  byte-for-byte against the corresponding raw PDB files in
  Sanithu's GitHub repository.
- All 45,000 entries in the 9NFU IMGT residue-mapping CSV
  match the finalized FASTA sequences using 0-based indexing.
- The 9NFU heavy domain occupies indices 3-123 inclusive;
  the light domain occupies indices 143-246 inclusive.
- The finalized FASTA files are tracked in the repository.
  Their original generation-code revision has not been
  independently established.
- Older generation manifests and the current generation script
  specify different design-mask boundaries. This version
  discrepancy remains unresolved.
- Original ESMFold-predicted structures and complete
  prediction-input records are unavailable. Consequently,
  exact correspondence between every saved structural
  measurement and finalized FASTA entry cannot currently
  be independently verified.
- Passing downstream statistical checks does not establish
  end-to-end reproducibility or experimental binding activity.

## Overall Interpretation

Adding AbLang to ProteinMPNN did not consistently improve
the evaluated structural properties across both antibodies.

Changes differed across structures and CDR regions.

The repository enables verification of downstream statistical
calculations while identifying the remaining requirements
for complete end-to-end reproducibility.

## P1c contract-only eligibility mode

The earlier numerical results and reported arithmetic checks above remain historical
evidence. This addition does not recompute them or promote them to final findings.
**Technical provenance, operator-attested policy approval, and scientific validity
are separate.** A `VERIFIED` eligibility report means the supported contract checks
passed; it never establishes biological validity, correct causal isolation, a valid
sampling model, binding, or population generalization.

### Explicit mode selection and import safety

The former no-argument command executed four numerical scripts. Mode selection is
now required, including for legacy use; importing the runner performs no work.
This intentional CLI migration prevents an eligibility check from accidentally
launching scientific arithmetic. Historical compatibility mode is explicit:

```sh
python scripts/run_reproducibility_checks.py --mode legacy-arithmetic
```

That mode runs the existing four numerical scripts unchanged and labels their
assurance as historical arithmetic with unverified input lineage. Missing
numerical dependencies cause failure; nothing installs dependencies. Eligibility
mode never dispatches those scripts and imports no statistical or model packages.

Read-only inspection of a P1b bundle, without inventing an analysis plan:

```sh
python scripts/run_reproducibility_checks.py --mode eligibility \
  --manifest /path/to/bundle/prediction_manifest.csv \
  --input-root /path/to/antibody-sequence-design \
  --upstream-root /path/to/ProteinMPNN-Ablang-antibody-design \
  --output-root /path/to/predictions \
  --report /path/to/new-eligibility-report.json
```

The report path must be new. Exit code 0 means the supported contract is complete;
exit code 2 means missing, unresolved, historical, ineligible or conflicting evidence.
That nonzero result is expected for the current real bundle: generation-code/run
lineage and scientific mapping authority are unresolved, and measurements and an
approved scientific plan are unavailable. Per-context P1b prediction coverage is
reported even without a plan. Missing observations remain absent, never zero.

Multiple `--manifest` arguments may represent independently identified prediction
replicates. If run contracts are supplied, give one `--run-contract` per manifest,
in the same order. Shared input designs are retained once by their full P1b design
identity; distinct prediction/run/replicate identities remain separate. Conflicting
run-instance definitions cannot be merged.

### Existing shared provenance interface

`statistics_contract.py` calls P1b's `validate_prediction_bundle` for observed input,
map, run and prediction evidence and uses `validate_measurement_joins` for structural
measurement identity/cardinality. It does not trust a saved `PASS` column, recreate
input hashes, or implement another manifest registry. Source snapshot SHA is never
substituted for generation revision. Actual P1b `NOT_ESTABLISHED` generation state
blocks real scientific eligibility even if every file hash matches.

The public eligibility runner accepts no prevalidated-status JSON shortcut.
`assess_validated_context` is an internal policy seam used by synthetic fixtures;
injected fixture evidence is not a route to validating real study data.

### Scientific policy and approval binding

Optional inputs are `--plan`, `--approval-receipt`,
`--approved-plan-sha256`, `--measurement-receipts`, and `--evidence-root`.
There is no live study-plan default. A complete JSON body, a hash, or an
`approved:true` field alone supplies no approval. The separate receipt has schema
`uirj-plan-approval-v1` and fields `plan_sha256,decision_id,decision_source,
decision_authority,decision_timestamp`. The operator must additionally name the
same digest with `--approved-plan-sha256`.

This verifies binding to the operator-attested decision record, not authenticated
human identity, the truth of the decision, or scientific adequacy. P1c implementation
authorization is not approval of particular study endpoints or statistical choices.
The synthetic positive fixture uses explicitly fictional approval evidence.

The plan schema is `uirj-statistics-plan-v1`. Required top-level fields are:

- `analysis_scope="fixed_cases"`, the cases 9NFU/9NH7 and methods baseline/ensemble.
- Exact `design_ids`, `endpoints`, and `evaluation_members` freezing intended
  populations and prediction replicates before considering available measurements.
- `hypotheses` and `families`, including unavailable planned members.
- `generation_interpretation="implemented_protocol_contrast"` and explicit
  `analysis_timing` (retrospective/exploratory or a separately supported new experiment).

An endpoint specifies case, metric, region, kind, units/direction, explicit numeric
bounds, metric policy and extraction code identities, reference/map hashes, atom
selection, alignment, domain, mask source, score definition, estimands, effect,
uncertainty, replication, missingness and pairing. Omitted choices remain blockers.
Supported metric/unit labels are `rmsd` or `contact_displacement` in `angstrom`,
`plddt` in `confidence_0_100`, and `contact_identity`/`oasis` in `fraction` or
`percent`. These labels are explicit schema vocabulary, not aliases guessed from
historical column names. Plan bounds must respect their inherent nonnegative or
0–1/0–100 domains; the software does not select tighter scientific thresholds.

This bounded contract supports draw-weighted primary membership and a separately
named unique-sequence-support sensitivity, retaining all evaluation members. It
requires explicit `statistic=mean_difference` and records E−B effects in raw metric units without reversing sign for lower-is-better
metrics. It creates membership lists and coverage counts only: no mean, median,
rank statistic, confidence interval, adjusted p-value, plot or aggregate is computed.
Unsupported estimands remain unresolved rather than silently substituted.

Descriptive uncertainty requires explicit `method=descriptive_only`, `unit=none`,
`hierarchy=none`. An inferential bootstrap declaration additionally needs explicit
unit, hierarchy, resample count, seed, confidence, interval algorithm and quantile
rule. Supported combinations are unpaired draws with `independent_draws`, unpaired
per-draw folding repeats with `draws_then_prediction_replicates`, or genuine pairs
with `paired_draws` and one fold per draw. Shared-sequence bootstrap hierarchies
remain unsupported and unresolved; no blanket row-independent bootstrap is implied.
These are policy consistency checks, not validation of bootstrap coverage.

For Holm/BH families, alpha and each hypothesis's null, test, sidedness and sampling
assumptions must be supplied. The bounded vocabulary distinguishes Mann–Whitney's
`equal_distributions` null, paired signed-rank's `zero_symmetric_pair_difference`,
declarations. Weak-mean permutation tests remain unsupported pending a separate
exchangeability/studentization and calibration contract; raw permutation is not
automatically exact for a weak equal-means null. No tests execute.
BH also requires a recorded dependence justification. An approved family retains
all planned hypothesis IDs; missing outputs cannot shrink it to the successful ones.
A test name or assertion of assumptions is not empirical proof they hold.

### Measurements, repeated draws and prediction replicates

Metric input schema is `uirj-metric-receipts-v1`, with a `measurements` list. Each
structural receipt retains P1b's exact core identity fields plus design, prediction,
measurement, endpoint, metric/region/unit, extraction-run and extraction-code IDs,
value or null with reason, and an explicit current/historical evidence class.

A measurement additionally binds `source_path`, exact source file hash,
1-based source row index, canonical row hash, `value_column`, and a `source_columns`
mapping covering every required identity. Matching only the numeric value is
insufficient. Extraction source path bytes must match `extraction_code_id`.
All paths resolve under the explicit evidence root. These checks bind identities,
values and named code bytes; they do **not** prove that code executed or that its
metric measures the intended scientific construct.

Existing legacy CSVs lack many of these identities. They are not automatically
eligible, and no inferred filename/row/sequence-similarity join fills the gaps.
Recovering genuine extraction receipts is a separate evidence task. Sequence-only
score-run provenance, including OASis, is explicitly
`SEQUENCE_SCORE_PROVENANCE_UNSUPPORTED`; planned endpoints/hypotheses remain visible.
No synthetic PDB receipt is fabricated to work around that remaining gate.

Each `evaluation_members` entry names endpoint, target design, prediction owner
design, actual prediction ID, replicate ID and configuration identity. The latter
is derived from the verified run contract excluding run-instance and input-bundle
identity, while retaining model/settings/seed policy. Different configurations
cannot be averaged under a single approved replicate rule.

Supported declarations are one prediction per draw, `mean_within_draw` membership
for genuine per-draw repeats, and explicit `shared_deterministic` sequence-level
reuse. No averaging occurs. Reuse requires identical verified sequence identity
and preserves each target draw's weight while counting the shared prediction once.
It cannot be labeled as independently generated folds. Equal sequence strings
alone do not establish dependence, and repeated designs are never silently deduped.
Unique-support membership is a sensitivity description of observed sequence support,
not an unbiased uniform-library estimator or evidence about antibody populations.

Pairing on nominal numeric IDs or equal seed labels is rejected. Genuine pairing
requires a hash-bound coupling receipt naming the protocol and complete one-to-one
baseline/ensemble design pairs for the case. Different paired sequences are allowed.
An unpaired plan cannot carry pairs. More within-case draws do not establish
between-target or antibody-format generalization from these two cases.

### Evidence states and missingness

Every report preserves reason codes alongside `MISSING`, `UNRESOLVED`, `CONFLICTING`,
`HISTORICAL`, `INELIGIBLE`, or contract `VERIFIED` states. Technical measurement states
and scientific-policy states are separate. All blockers remain visible rather than
being replaced by a single optimistic status. Duplicate/conflicting identities are
quarantined without choosing a first/last value; normalized diagnostics and member
sets are invariant to receipt/manifest order when source identities remain fixed.
Raw file hashes may legitimately change when original source files change order.

A absent source row is not proof of technical failure. Null requires an explicit
missing reason and is never converted to zero. Complete-case descriptive policies
must declare `observed_valid_subset`; they do not establish effects for the full
planned generated distribution. Incomplete scientific evidence still blocks final
inferential eligibility. Structural preservation/noninferiority and clinical or
binding claims are not implemented by this contract.

### Bounded acceptance

Under the external CPU/RSS/time monitor:

```sh
python -m pytest -q tests/test_statistics_contract.py
```

Synthetic fixtures exercise the valid two-case contract and mutations of source
identity, masks, units/ranges, missingness, duplicates, ordering, coupling,
replicate/reuse rules, frozen families, and approval binding. Legacy dispatch is
mocked; no historical arithmetic runs in the test suite. A real-bundle eligibility
CLI check should retain all 400 source design records and 363 unique sequences
while returning explicit unresolved/missing evidence. Actual execution receipts,
including failures, belong to the project checkpoint; this checklist is not itself
proof that tests passed.


The report emits `analysis_plan_sha256` and the exact operator approval binding
(decision ID, source, authority and timestamp, or null). These are immutable evidence
identities, not semantic normalization: reordering arrays in a plan changes its digest
and requires a matching approval binding. Row-permutation acceptance compares the
semantic decisions/member sets separately from those explicit plan identities.

Unique-support sensitivity requires `unique_support_weight=
"equal_per_sequence_within_case_method"` and `within_sequence_rule=
"mean_over_distinct_prediction_ids"`. Each observed unique sequence receives equal
support weight within its case/method, with all distinct same-configuration evaluation
predictions retained inside that sequence; a reused deterministic prediction is not
counted again for each copied draw. This explicitly differs from equal draw weighting
and does not claim either computation has occurred. Missing aggregation policy is an
unresolved contract, not an implicit averaging choice.

Identity counts distinguish total known prediction identities from
`unambiguous_prediction_identity_count`. Conflicting variants are retained in named
quarantine diagnostics and cannot enter verified membership; they are not silently
removed from evidence counts. Equivalent source design references across independent
prediction-replicate bundles do not create new generated draws.
