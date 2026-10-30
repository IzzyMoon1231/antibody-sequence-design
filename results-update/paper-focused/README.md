# Run, inspect and explain your results

Open **[results.ipynb](results.ipynb)** in local JupyterLab or VS Code. Its saved outputs use the included **real public example**. For reading without Python: [GitHub reader](READ_RESULTS.md), or download `results.html` for a self-contained mobile/offline view (GitHub does not host that HTML).

## Reading the four main panels

The notebook and readers show **9NH7 similarity → 9NH7 structure → 9NFU similarity → 9NFU structure**, each followed by a caption calculated from the selected frames and a reflection question. The four plotting cells call `analysis.plot_panel(data, case, metric)` directly; no saved image is loaded. Colors are baseline teal and ensemble rust. Comparable panels share x limits; structural limits include the full range across both selected cases, so changed inputs cannot silently clip plotted intervals. Identity retains its meaningful 0–1 scale with padding at the endpoints.

Structural **dots are medians**, and thick lines are the **middle 50% (IQR), not confidence intervals**. Similarity dots retain each draw, with black median ticks. All loop changes remain visible, including the other heavy/light loops that H3 misses. Captions and counts update with the data; original audit qualifications do not transfer to replacement inputs.

Opt-in exports include `9NH7_similarity`, `9NH7_structure`, `9NFU_similarity`, and `9NFU_structure` as PNG/SVG, plus `figure1` as a compact **2×2 manuscript layout** (cases in columns, similarity/structure in rows, A–D). Both layouts use the same drawing helper and loaded frames. Use the individual panels on a phone; the combined download is intended for manuscript review and may need zoom. `results.html` embeds its reading images for offline use; keep the `figures/` folder alongside it for download links.

## First run — then your own results

1. Keep this folder together. Use Python 3.11+ with `requirements.txt`. Optional setup in your own environment: `python -m pip install -r requirements.txt`. No model/database downloads or credentials are needed. No dependency installation was performed for this revision.
2. Open the notebook with its folder as the working directory. Leave the three input selectors unchanged and **Restart Kernel and Run All**. Export defaults to **False**. Read the small actual input previews and validation counts before the figures.
3. Predict the direction of one comparison, check it against the visible group medians and differences, and explain it in your own words. Check the other loops and confidence counts rather than only the favorable endpoint. Dots, medians and IQRs are descriptive, not confidence intervals.
4. To use your own compatible files, create `inputs/my_run/` and **copy** your files using Finder/Explorer, VS Code's local folder explorer, or the Upload button in a **locally running** JupyterLab file browser. Do not use a remote kernel/cloud upload for private results. No widget or network transmission is required; selecting a path only reads that local file. Never overwrite the included example.
5. Copy `selection.csv` to `inputs/my_run/selection.csv`, retain all four case/method rows and update each workbook/FASTA path. In notebook step 1 choose that CSV, your structural CSV and your ablation CSV (or `None` to omit ablation). **Every path is relative to the notebook package folder**, even paths inside a nested selection CSV. Absolute paths, URLs, `..`, and symlink escapes are rejected. Example:

   ```python
   SELECTION = 'inputs/my_run/selection.csv'
   STRUCTURAL = 'inputs/my_run/loops.csv'
   ABLATION = 'inputs/my_run/ablation.csv'  # or None
   ```

   Example selection row (keep the other three rows too):

   ```csv
   structure,method,workbook,fasta
   9NH7,baseline,inputs/my_run/9NH7_baseline.xlsx,inputs/my_run/9NH7_baseline.fa
   ```

6. Restart and run again; inspect the rows and counts. If a check fails, stop, fix the named export/key/column, clear old notebook outputs, and restart. Do not fill missing measurements or discard adverse results. No failed validation changes earlier exported folders. A previous saved notebook/reader remains a **historical snapshot**, not a fallback result for new inputs.
7. Complete the notebook's five reasoning notes. With export still off, restart/run all and save. Then change the **last** cell to `EXPORT = True` and run it to create a new `exports/<UTC timestamp>-<unique suffix>/` folder. It contains figures (PNG and editable SVG), tables, same-source readers, selected input bytes, hashes/versions and your notes. Save a copy of your executed notebook into that folder as `results.ipynb`; `REPLAY.md` explains the selected-inputs root for replay. Reset export to False before your next restart. Check privacy and permission before sharing your own data.

### A short interpretation model

For the **included example only**: the 9NH7 ensemble heavy-domain similarity median is higher, but its H3 RMSD median is also higher. This supports a descriptive trade-off in these retained predictions, not preserved binding. For your own run, supply the actual medians, direction, counts and an unfavorable or heterogeneous result. A hypothesis about the cause needs separate evidence; do not ask AI to certify a scientific conclusion. Ask your mentor about a new metric, reference, case or statistical design.

## Supported schema — a compatible-result reader, not arbitrary-data analysis

The two fixed cases **9NH7 and 9NFU** and the named endpoint definitions remain mandatory. Sample counts are **any positive number**, including unequal main-arm counts; small counts have weak descriptive coverage and do not create biological replication. All loops within one arm must cover exactly its scored designs. No resampling, pairing, confidence filtering, deduplication or inferential test is introduced automatically.

| Input | Required content and checks |
|---|---|
| Selection CSV | Exactly `structure,method,workbook,fasta`; one row per 9NH7/9NFU × baseline/ensemble. Nonblank local relative paths. |
| BioPhi `.xlsx` | `Overview` sheet with unique nonblank `Antibody`, `Heavy OASis Identity`, `Threshold` = `relaxed`; `VH` with unique `Antibody`, `scheme` = `kabat`, `chain_type` = `H`, and numbered residue columns. Identity is finite 0–1. Heavy-only `_VH` rows are selected when present. Baseline numeric suffix is one-based; ensemble suffix zero-based. Reconstructed VH must occur in its indexed full FASTA. No pooled Identity or percentile substitution. |
| FASTA `.fa` / `.fasta` | Canonical nonempty amino-acid sequences, same indexed draw order and count as scored VH rows. Headers may repeat in the original source; record order is preserved. Scored indices must cover 0 through record-count minus one exactly. Repeated sequences remain separate draws. |
| Structural CSV (required) | `structure,method,design_id,cdr,domain,rmsd_angstrom,n_atoms_per_residue,cdr_mean_plddt,low_confidence,framework_rmsd_angstrom`. Original optional columns `n_residues,framework_outlier` are retained. `design_id` = `CASE_METHOD_N` with zero-based N matching the row's case/method (e.g. `9NH7_baseline_000`). Unique case/method/index/loop keys. All H1/H2/H3 for each case, plus L1/L2/L3 for 9NFU, with complete scored-design coverage. 9NH7 heavy domain VHH; 9NFU heavy VH, light VL; three N/CA/C atoms per residue. RMSDs finite/nonnegative; pLDDT finite 0–100; flag exactly equals mean loop pLDDT < 70. Only baseline/ensemble analyzed; known `abmpnn` comparator rows explicitly excluded. |
| Ablation CSV (optional as a whole) | `id,Heavy OASis Identity,Light OASis Identity`. IDs exactly `CASE_cdrh3_METHOD_sN` or `CASE_others_METHOD_sN`, unique canonical nonnegative seed N. Both methods for each seed in all four case/condition combinations; pair counts may differ across conditions. Heavy scores and 9NFU light scores finite 0–1; no missing partners. Existing qualifying-window denominators 107 (9NH7 H), 113 (9NFU H), 96 (9NFU L) enforced. 9NH7 light values unused. Set selector to `None`, not a missing filename, to omit this section explicitly. |

File-size limit is 50 MB per input. Use trusted data exports, not executable notebook/script uploads. Files with matching columns do **not** establish matching measurement protocols: structural values must use the documented domain-local N/CA/C deposited-reference fit, operational motif spans, and unchanged case/domain conventions. Different alignment references, full-scFv fits, new domains/cases or peptide denominators require a reviewed adapter; do not rename columns to bypass this boundary. See [SOURCES.md](SOURCES.md).

Changed bytes (including selection metadata) lose the original coordinate/execution qualification. The notebook validates only selected bytes and their computational correspondence. An uploaded provenance file cannot confer the pinned audit's status. Edited in-memory previews are rejected; change source exports and read again. Main arms are independent collections, not row-paired; ablation uses reported seed pairs with conditional causal linkage. Higher heavy-domain similarity is not reduced immunogenicity, and RMSD is not binding.

## Reproduction and implementation

- `analysis.py`: reviewed descriptive formulas and figure constructors, reused rather than replacing the science.
- `workflow.py`: local byte snapshots, input/identity/schema/range gates and opt-in new-folder export. It neither scores nor folds.
- `build_notebook.py`: generates the teaching cells and executes them in a clean IPC kernel. `python build_notebook.py` refreshes **only the notebook and an execution receipt**; it does not export figures/tables or replace a results bundle. On execution failure it leaves the previous saved snapshot intact and exits with an error; that snapshot is not the failed run's result.
- `python -m unittest discover -s tests -v`: actual original-input arithmetic, selected-row plot artists, changed/smaller **synthetic** fixtures, fresh-kernel cells, path/identity/range failures, and export preservation. Expected tables are frozen original-result checks; never edit them to make a new dataset pass.
- Included `tables/`, `figures/`, `READ_RESULTS.md`, `APPENDIX.md`, `results.html` remain the checked original example. Your new runs live in new export folders. No stale PNG is loaded by a notebook analysis cell.

The notebook and readers share calculation/plot helpers. No manuscript edits, publication or messages are part of this workflow. [Author checklist](AUTHOR_CHECKLIST.md) separates manuscript integration from computation; [engineering notes](engineering/README.md) distinguish historical receipts from current verification.
