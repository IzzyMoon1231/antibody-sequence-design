# Read or update these results

**Start with [READ_RESULTS.md](READ_RESULTS.md)** for the answer, one main table and one separate ablation figure. [results.ipynb](results.ipynb) is executed with embedded output. Download `results.html` for a self-contained phone/offline reading view; GitHub shows HTML source, not a hosted website.

## One update route

1. Keep this entire directory together. Open `selection.csv` in a spreadsheet: one row per case/method selects a public BioPhi workbook and its final-filled FASTA. To replace a workbook, copy the replacement into `inputs/` and change that row's relative filename; do the same for its matching FASTA. The structural export is `inputs/cdr_rmsd_results.csv`; the separate ablation export is `inputs/ablation_scores.csv`.
2. Use Python 3.11 or later with the packages in `requirements.txt`. In an environment where these are not already installed, the optional setup is `python -m pip install -r requirements.txt`. No models, weights or scoring databases are required. No installation was performed to prepare this package.
3. From this directory, run `python build_notebook.py`. This launches a fresh kernel, regenerates the table, figure, readout, supporting plots and mobile HTML, then saves the executed notebook. Alternatively open `results.ipynb` here and **Restart Kernel and Run All**, then save it. Do not edit generated prose or chart values by hand.

`python -m unittest discover -s tests -v` runs the portable tests from this directory. The expected tables are frozen checks for the supplied original inputs, not an oracle for your replacement dataset. If deliberate replacement data are used, do not rewrite expected tables merely to make tests pass: review the new cohort/endpoint and archive it separately.

### Supported data contract

This is a deliberately narrow reader for the selected cohorts, not a general antibody-analysis pipeline:

- Exactly 9NH7/9NFU × baseline/ensemble, 100 indexed draws per arm. Workbook `Overview` needs `Antibody` and `Heavy OASis Identity`; `VH` holds numbered residues plus `Antibody`, `scheme`, `chain_type`. The heavy-only row suffix is selected only where present. Baseline score IDs are one-based (subtract one); ensemble IDs are zero-based. Every reconstructed VH must occur in its indexed full FASTA. This operational correspondence does not validate biological domain boundaries.
- The structural CSV follows the public `compute_RSMD.py` export: `structure,method,design_id,cdr,domain,rmsd_angstrom,n_residues,n_atoms_per_residue,cdr_mean_plddt,low_confidence,framework_rmsd_angstrom,framework_outlier`. Main H3 rows must be VHH (9NH7) or VH (9NFU), three backbone atoms per residue, complete keyed coverage, finite nonnegative RMSD and pLDDT in 0–100. `low_confidence` must equal H3 mean pLDDT < 70. A file with the right columns is not evidence of the right measurement: preserve the domain-local N/CA/C protocol and deposited references.
- Ablation IDs must be `CASE_cdrh3_METHOD_sN` or `CASE_others_METHOD_sN`, N=0–99, with both methods per seed. Heavy and light scores stay separate; qualifying-window denominators are 107 (9NH7 H), 113 (9NFU H) and 96 (9NFU L). Missing pairs, duplicate IDs, invalid scores or incompatible denominators stop the update. The pooled `oasis` column is never substituted for Heavy Identity.

Changed bytes trigger a visible **Changed-input analysis** warning and do not inherit the original coordinate/execution qualification. A failed update invalidates the generated reader and removes charts/tables; the builder also clears notebook outputs. Plain schema errors identify the input to repair. Neither replacing data nor passing tests authorizes causal, binding-preservation, non-inferiority or clinical claims.

## Files that matter

- `analysis.py`: inspected, bounded refactor of retained-output calculations; no team code imports or network calls.
- `build_notebook.py`: notebook source and clean-kernel execution route.
- `inputs/`, `selection.csv`: selected actual public inputs and pinned provenance.
- `tables/`, `figures/`, `READ_RESULTS.md`, `APPENDIX.md`, `results.html`, `results.ipynb`: generated reading surfaces.
- [SOURCES.md](SOURCES.md): endpoint/cohort boundaries and coordinate-replay instructions.
- [engineering/README.md](engineering/README.md): verification, separate from scientific conclusions.

Publication is a separate parent review step. No manuscript edits, team messages or public pushes were made by this package.
