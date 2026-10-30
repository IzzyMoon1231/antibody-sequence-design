# UIRJ saved-results review

Open **Update_results.ipynb** on GitHub to read the saved, executed cells; no kernel is needed to view the embedded figures. On a phone, scroll through the OASis and structural sections. Download preview.html for a self-contained local output-only view (GitHub does not render HTML files as pages).

This is a review-branch package, not a final-paper update. Prepared with AI-assisted implementation and checks using Hermes (OpenAI model); human scientific review remains pending. The inputs were already public repository artifacts. All values remain **UNBOUND_SOURCE_EXPORT**: the exact scored/folded sequence and original run have not been authenticated. No inference or scoring is performed.

## Run or update

1. Keep the notebook, helpers, inputs.json and inputs folder together. Use a Python environment with requirements.txt dependencies.
2. Open this folder in Jupyter/VS Code. For the supplied snapshot, leave settings unchanged and Restart Kernel → Run All.
3. For new results, copy original files into inputs/new-delivery/. Copy inputs.json to my_inputs.json. Select that filename in the first notebook cell.
4. Edit paths, cohort labels and explicit value/domain columns. Use a new cohort for a genuinely different run. Remove unwanted entries. Do not relabel a percentile as Identity or a loop score as framework RMSD.
5. Read the selected-file preview, then run all cells. Every selected result group is plotted before export. Exports are new dated folders; old results are never replaced.

Supported inputs: sample/seed-bearing FASTAs; Jacob designs_filled CSV; BioPhi summary CSV or Overview XLSX; framework RMSD CSV; CDR RMSD v3 CSV; promb saved-score CSV. Fractions must already be 0–1; RMSD uses angstroms. Missing values remain missing. No automatic ID shifts, pairwise joins, scoring, statistical tests or model downloads occur.

## Read the figures

Curves are empirical cumulative distributions: the vertical coordinate is the fraction of nonmissing saved values at or below the horizontal value. Legends show present and missing counts. Color consistently identifies method. Each panel preserves source, cohort, case, domain and endpoint. Shared endpoint scales expose outliers; a narrow cluster can therefore look nearly vertical. Curves are descriptive saved draws, not independent biological replicates or a fair causal comparison.

BioPhi source-selected score, Overview Identity and promb peptide content are distinct definitions. The workbook repeats the 9NH7 ensemble CSV population and must not be counted as extra independent data. Framework and CDR-loop RMSDs are separate; no geometry or residue mapping was validated. Larger OASis is not demonstrated lower clinical immunogenicity; smaller RMSD is not better binding.

## Snapshot and limitations

Read-only GitHub HEAD checks still returned Sanithu bc3a1070a4e63489ef785763ae2193c4ffd5aa19 (Oct25) and Isabella f76c09fa3458da7c33daa76de150ab5f8c9a96c9. Selected inputs are unchanged. The latest Jacob Drive delivery was **not recovered**. source-check.json records the check; acquisition.json and xlsx-acquisition.json record selected public input hashes.

The selected files contain 1,000 sequences and 4,000 imported measurements (including the duplicate-view workbook), summarized into 40 source-specific groups and 22 figures. There are 400 same-ID sequence disagreements across collections and six possible structural ID-offset groups. Neither is repaired by guessing. No result is promoted into final paper evidence.
