# UIRJ results — start with the scientific analysis

**The current recommended notebook is [paper-focused/results.ipynb](paper-focused/results.ipynb). Start with [its scientific reading page](paper-focused/READ_RESULTS.md): one main table, one ablation figure, and the interpretation.**

The file-by-file explorer below is retained as historical diagnostic material. Its blanket sequence-linkage warnings are superseded for the selected, checked inputs in the new analysis.

## Earlier file explorer

Review-branch material prepared with AI-assisted implementation and checks using Hermes (OpenAI model). Human scientific review remains pending.

**Start with [the illustrated reading page](READ_RESULTS.md)** on GitHub. It has every plot, source label, computed readout and limitation. For a shorter phone view, download `mobile.html` and open it in a browser: its plots are embedded, so the reading view works offline. GitHub does not host HTML files as pages. Full-size image links need the adjacent `figures/` folder; embedded figures can also be zoomed offline.

[Executed notebook](Update_results.ipynb) · [Selection spreadsheet](input_selection.csv) · [Source check](source-check.json)

> **Scores not yet linked to selected sequences.** Obtain the actual scored/folded sequence file, run settings and structural mapping. Matching design IDs alone is not proof. The latest Jacob Drive delivery was not recovered; these are the Oct25 repository inputs pinned to Sanithu `bc3a1070`. No fresh source check is claimed by this presentation update.

## Reader route

1. Read the overview and file/sequence warnings.
2. Choose **9NH7** or **9NFU**, then humanness or structure.
3. Read the axes before each plot, the actual median/count readout afterwards, and its limitation.
4. Keep the warnings with any shared figure. These are descriptive saved draws, not a causal method ranking or final paper evidence.

## Updater route — no JSON editing required

1. Keep this folder together; use a Python environment with `requirements.txt`. Open the folder in Jupyter or VS Code.
2. Copy new original files into `inputs/new-delivery/`; do not overwrite old files.
3. Open **input_selection.csv** in a spreadsheet editor. One row selects one input file. Change the path and cohort for replacements, duplicate a compatible row for another file, or delete unneeded rows. A different run needs a different cohort. Save as CSV, retaining column names.
4. In the notebook's **Setup · the only section to edit**, choose that selection file. Leave export and previous-bundle settings alone unless needed.
5. **Restart Kernel → Run All.** Check the selected paths and mismatch warnings before interpreting plots.
6. Open the newly dated `exports/update-.../` folder. Its `mobile.html`, `READ_RESULTS.md`, figures, exact tables and copied inputs all come from that same execution. Old exports are never replaced. The top-level reading pages remain the supplied snapshot; a new run writes its own pages in its bundle.

Malformed sequences, duplicate IDs and invalid values stop with an error. Missing files remain explicit warnings. If an input changes after preview, restart the run. A run with missing inputs is not complete scientific evidence.

## Selection columns

- `type`: `fasta`, `designs` or `metrics`.
- `path`: relative to this folder; preserve original files.
- `cohort`: the collection/run label, not merely the method.
- `structure`, `method`: required for FASTAs and workbook inputs where not present in rows. Cases are 9NH7/9NFU; methods baseline/ensemble/AbMPNN.
- Metrics: preserve `id_column`, `value_column`, `domain_column`, `endpoint`, `unit` and workbook `sheet`. Empty domain means unspecified, not inferred. Fractions must already be 0–1; RMSD uses angstroms.

Supported files: sample/seed-bearing FASTAs; Jacob `designs_filled.csv`; BioPhi summary CSV or **Overview** XLSX; framework RMSD CSV; CDR RMSD v3 CSV; promb saved-score CSV. For a genuinely new schema, ask someone to review the parser rather than relabelling columns to force a match. `inputs.json` retains default metadata and limitations; it is not the routine editing interface.

## Interpretation boundaries

BioPhi source-selected score, Overview Identity and promb peptide content remain distinct. The 9NH7 ensemble workbook repeats a CSV population and is not another independent sample. Cases, cohorts, domains and sources are not pooled. Framework and loop RMSD remain separate; geometry and residue mappings are not validated. OASis is not clinical immunogenicity; RMSD is not binding affinity. No automatic ID shifts, paired tests, confidence intervals, correlations, model downloads or inference occur.

Counts are generated in [reading-data.json](reading-data.json), not manually maintained here. Full-precision data are retained in the notebook's dated export; display medians are rounded only for reading. [Design notes](DESIGN_NOTES.md) explain the presentation choices.
