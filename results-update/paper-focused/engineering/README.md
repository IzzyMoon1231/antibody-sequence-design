# Engineering verification — not the scientific reading path

The student-facing result is [READ_RESULTS.md](../READ_RESULTS.md). Engineering checks establish computation and artifact behavior, not independent biological or causal acceptance.

## Current four-panel revision

The main notebook now has four separately executed case/metric figures and local, data-derived captions. The combined manuscript export is 2×2, using the exact same `draw_panel` helper. `tests/test_panels.py` compares actual collection offsets, interval segments, colors, median lines, labels and limits; checks case/letter titles; mutates structural range and similarity; and checks notebook cell separation plus both reader/export surfaces. Existing numerical, provenance, missing-input, ablation, non-overwriting export and fresh-kernel synthetic controls remain.

Current receipts: `panels-notebook.json`, `panels-refresh.json`, `panels-suite-green.log`. Earlier receipts below remain historical and were not rewritten. The private maintainer-only `refresh_panels.py` command (not part of the public allowlist) refreshes the checked reading snapshot explicitly; it is never called by notebook Run All. It compares every retained table column at 1e-12 and preserves original input/table/expected-CSV bytes. Default notebook export remains off.

Local desktop/phone HTML previews were reviewed, not GitHub's remote renderer. The combined manuscript image is downloadable; on a phone its text needs zoom, so the reading surfaces use individual panels. Publication and independent scientific acceptance are separate from local execution.

## Previous transparent notebook revision

The sections below this notice are **historical receipts**, not descriptions of the
current updater. Current behavior is select → preview → validate → calculate → plot
from memory → reason → explicit opt-in export. A normal Run All writes no analysis
outputs. Failed input checks preserve earlier exports instead of deleting them.
The notebook never reads an existing PNG. `workflow.py` pins original qualification
to selected bytes and rejects mutated in-memory previews. Counts need not be 100;
main arms may differ, while complete loop coverage and explicit ablation pairs remain
mandatory. Ablation can be explicitly omitted; the main structural section cannot.

Current verification: [transparent-notebook.json](transparent-notebook.json) and
[transparent-suite-green.log](transparent-suite-green.log). Fresh kernels exercise
original real data, smaller unequal-count/reversed synthetic controls, explicit
ablation omission, the actual final export cell, and missing/invalid inputs.
Artist coordinates and image/prose changes are checked against selected rows;
no saved figures are copied into the cell-workflow test workspace. Input/output
hash snapshots verify no analysis writes before export and preservation afterwards.
The old deletion/path-Image assertions were deliberately replaced for the requested
non-overwriting policy, not weakened to admit missing measurements. All original
numeric expectations remain unchanged.

Semantic REDs for hidden publish-all setup, mutated-preview provenance inheritance,
wrong scoring threshold, and missing preview columns are retained locally (excluded
from publication because tracebacks contain private paths). The first synthetic
subset fixture failed because lexically ordered antibody IDs were sliced by row;
it was corrected to select by numeric identity. This is a fixture correction,
not a scientific result. No inference, upstream scripts, network or installs ran.

## Historical executed checks

- `test_analysis.py` compares every numeric Table 1 column to the independently checked retained table, and all numeric ablation summary columns to the checked chain-specific output (tolerance 1e-12). Expected tables are derived from public observations; they contain no manuscript text.
- The changed-input control edits actual copied workbooks, H3 CSV and ablation CSV. It reverses humanness/structural directions and ablation signs, then checks regenerated prose, HTML, plotted mean values and both main/ablation image bytes. These mutations are synthetic tests, never presented as observations.
- Missing-input tests invalidate the readers and remove stale figures/tables. Schema and scientific-gate controls cover unsupported atom semantics, unmatched keys, nonfinite values, confidence flags, duplicate/missing ablation pairs, incompatible peptide denominators and sequence mismatch.
- `build_notebook.py` executes from a new kernel, checks every code cell executed, requires three embedded PNGs and rejects error outputs. See [notebook.json](notebook.json).
- `run.json` records exact selected input hashes and whether they match the original public provenance. No environment dump or machine paths are included.

## Test-first chronology

Observed RED before the first implementation: missing portable analysis helper (API prerequisite). GREEN: real Table 1 agreement. Next RED: missing regeneration/chart interface; GREEN: retained ablation agreement and changed-input-to-artists/image/prose propagation. Later REDs exposed surviving stale table files and a nonspecific missing-column error; the implementation now removes stale tables and names the missing structural column. A separate narrative RED required the data-conditional no-structural-cost qualification. The final suite retains the resulting controls. These are sequential vertical implementation slices, not scientific experiments.

The scientific-review correction added a semantic RED for missing complete loop coverage, then a GREEN covering all nine case/loop summaries against original selected rows, corrected narrative on both readers, and changed L2 input propagation. A second RED demonstrated a stale loop table after a missing loop row; GREEN now invalidates all summaries. An initial test launch lacked PYTHONPATH and failed import before the semantic RED; it was a command setup error, not a scientific defect. No new model, scoring or coordinate experiment was performed.

## Final communication pass

Four further vertical slices were observed RED then GREEN: absent main comparison API; mean-only ablation bars rather than reported observations (semantic failure); absent second figure/experiment-map reading surfaces; and an opening/notebook hierarchy that did not meet the plain-language contract. Local raw RED logs are not published because tracebacks contain machine paths. `test_figures.py` compares main observed-point coordinates, all nine loop medians/IQR segments, ablation seed-pair dots and means directly with selected rows, including reversed-sign controls. It checks the experiment map, separated author checklist and invalidation of every figure after missing inputs. Main and ablation image bytes, prose and HTML change with mutated files. Thirteen test methods pass in the final suite.

No new biological or inferential experiments were performed. Main and ablation numeric summary tables retain the preceding checked values. Figure intervals describe observations; none estimates uncertainty. Table 1 and the full loop/ablation tables remain numerical references.

## Prior coordinate evidence

[prior_coordinate_replay.json](prior_coordinate_replay.json) is an unchanged, small receipt from the preceding coordinate replay. It is deliberately distinguished from this package's newly executed summary regeneration. It records the source/AST/reference hashes, CDR/framework row counts and zero numeric discrepancies. The large PDB collection is not shipped. [SOURCES.md](../SOURCES.md) describes the public paths and replay procedure.

## Reading-surface checks

The offline HTML was inspected in a local Chromium browser at phone width (390 CSS px) and desktop width. Phone body width equals the viewport; the main table has its own horizontal scroller. All three main images load; the spread/confidence disclosure has its own table scroller. Main Figure 1, seed-pair Figure 2 and the editable experiment map were inspected at both widths, with no clipping. Screenshot evidence is retained locally, excluded from the public allowlist. This is local-render verification, not a claim about a published GitHub renderer.

## Workflow for a successor

Keep science in the tested helper; regenerate every reading surface from the same result objects. Do not patch numerical prose by hand. Check both unchanged-input agreement and deliberately changed-input propagation, then clean-kernel execution and mobile rendering. Keep publication allowlists separate from generated caches and kernel connection files. Review the final allowlist before copying this directory into a public branch.
