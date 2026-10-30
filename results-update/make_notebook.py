"""Build the student notebook; all observations come from the selected run."""
from pathlib import Path
import nbformat as n
root=Path(__file__).resolve().parent
cells=[]
def md(text): cells.append(n.v4.new_markdown_cell(text))
def code(text): cells.append(n.v4.new_code_cell(text))
md('''# What do the saved antibody results tell us?

**Scores not yet linked to selected sequences.** Before using these as final paper results, obtain the exact sequence files used for scoring/folding, the run settings and structural mapping. Matching IDs are not enough.

**Reader route — no code needed:** [overview](#overview) → [file and sequence checks](#checks) → [9NH7](#9nh7) or [9NFU](#9nfu) → [next steps](#next). For a shorter GitHub reading view use [READ_RESULTS.md](READ_RESULTS.md); download `mobile.html` for a phone-friendly offline view. GitHub displays HTML source, not a hosted page.

**Updater route:** change only the marked setup section and the selection spreadsheet, then **Restart Kernel → Run All**. Read the checks before interpreting any graph. Exports never replace earlier bundles.

By the end you should be able to identify which files were used, explain a distribution plot, and name what still blocks a paper claim. Nothing here scores proteins, runs inference, uploads or sends messages.

Included data: the Sanithu repository snapshot `bc3a1070` (Oct25). **The latest Jacob Drive delivery was not recovered.** This notebook reviews saved exports; it does not establish their scientific lineage.''')
md('''## Setup · the only section to edit

First run? Leave everything unchanged. For a new delivery:
1. Copy the original files into a new folder under `inputs/`.
2. Open `input_selection.csv` in a spreadsheet editor. Each row selects one file. Change its path and cohort; a different run needs a distinct cohort. Delete rows you do not want, or copy a compatible row for another file.
3. Keep the correct ID, value and domain columns; do not rename an endpoint to make it comparable. [GUIDE.md](GUIDE.md) describes supported formats. JSON is an internal default, not a required editing interface.

Open this notebook's folder in Jupyter/VS Code. The next cell loads your selection; it does not run any scientific model.''')
code('''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Markdown, Image
import update_results as update
from reading_views import load_selection, write_reading

# EDIT HERE ONLY (paths relative to this notebook's folder)
PROJECT_FOLDER = Path.cwd()
SELECTION_FILE = "input_selection.csv"
EXPORT_FOLDER = PROJECT_FOLDER / "exports"
PREVIOUS_BUNDLE = None  # optional: "exports/update-..." for input-byte comparison

config = load_selection(PROJECT_FOLDER / SELECTION_FILE,
                        json.loads((PROJECT_FOLDER / "inputs.json").read_text()))
config["previous_bundle"] = PREVIOUS_BUNDLE
print("Selection loaded. Inspect the file list below before continuing.")''')
md('''<a id="checks"></a>
## 1. Are these the right files?

Check paths, cohort labels and value columns. A cohort is one named collection from one run, not a method name. Each listed file is actually read; missing files remain visible warnings rather than silently falling back to older data.''')
code('''for kind in ("fasta", "designs", "metrics"):
    display(Markdown("**" + kind + "**"))
    for item in config[kind]:
        display(Markdown(f'- `{item["path"]}` — **{item["cohort"]}**; value: `{item.get("value_column", "sequence")}`; domain: `{item.get("domain_column") or "not specified"}`'))
run = update.load_run(PROJECT_FOLDER, config)''')
md('''<a id="overview"></a>
## Overview · what is available, and what is not?

These counts describe selected files, not independent biological replicates. The workbook is a repeated view of the 9NH7 ensemble CSV. Keep it separate; never add it as another independent sample. The figures below preserve every source, cohort, case, domain and endpoint.''')
code('''mismatches = sum(r["status"] == "SEQUENCE_MISMATCH" for r in run["sequence_comparison"])
offsets = sum(r["possible_index_offset"] for r in run["coverage"])
display(Markdown(f"**{len(run['sequences']):,} sequences · {len(run['metrics']):,} imported measurements · {len(run['summaries'])} source-specific groups**"))
display(Markdown(f"**Check before interpreting:** {mismatches} same-ID sequence disagreements; {offsets} possible ID-offset groups. No offsets are applied. Scores remain unlinked to selected sequences."))
for item in run["missing"]:
    display(Markdown("- " + item))''')
md('''## 2. Do the sequences match?

The first table counts exact sequence text within each cohort. Repeated sequences stay as repeated draws. The second compares the same case/method/ID across collections, without joining scores or pairing methods.

**What to do:** if sequences disagree, ask for the sequence file actually scored and folded. Do not copy a hash from the desired FASTA or shift row numbers. An exact text match still does not prove a model used that input.''')
code('''seq = pd.DataFrame(run["sequences"])
if not seq.empty:
    display(seq.groupby(["cohort", "structure", "method"]).agg(draws=("record_id", "size"), unique_sequences=("sequence_sha256", "nunique")))
else:
    print("No selected sequence files were loaded. Sequence identity cannot be checked.")
comparison = pd.DataFrame(run["sequence_comparison"])
if not comparison.empty:
    display(comparison.groupby(["structure", "method", "status"]).size().rename("records").to_frame())
else:
    print("No cross-cohort comparison is available.")
coverage = pd.DataFrame(run["coverage"])
if not coverage.empty:
    display(coverage.groupby(["cohort", "structure", "method"]).agg(groups=("domain", "size"), possible_offsets=("possible_index_offset", "sum"), missing_values=("n_missing_values", "sum")))
print("Exact disagreements and unmatched IDs are preserved in sequence_comparison.csv and coverage.csv on export.")''')
md('''## 3. How should I read the plots?

Pick a case below, then humanness or structure. All panels are available; the [mobile reading view](mobile.html) groups them into expandable sections.

- **Horizontal:** saved score (0–1) or reported RMSD (Å), as labelled.
- **Vertical:** fraction of nonmissing saved values at or below the horizontal value. At height 0.5, read the median on the horizontal axis.
- **Legend:** present and missing counts; stable colors identify methods. Lines are descriptive, unpaired saved draws.
- **After each plot:** a computed median/count readout and a limitation. Display medians are rounded, not the exported values.

A rightward humanness curve means larger saved scores, not lower clinical immunogenicity. A leftward RMSD curve means smaller reported deviations, not stronger binding. BioPhi source-selected score, Overview Identity and promb peptide content are different definitions. Framework RMSD and CDR-loop RMSD are not interchangeable. Shared endpoint scales preserve tails; narrow groups may look almost vertical.''')
code('''from result_figures import build_figures

def render_sections(rows, case, unit):
    # Shared plot builder; only display selection changes, never source rows.
    for key, fig in build_figures(rows):
        source, cohort, structure, domain, endpoint, measure_unit = key
        if structure != case or measure_unit != unit:
            continue
        display(Markdown(f"### {domain} — {endpoint}\\n\\n**Question:** Where do these saved values fall? **Cohort:** {cohort}. Source: `{source}`.\\n\\nRead the horizontal measurement against the cumulative fraction on the vertical axis; inspect the median at height 0.5 and the legend counts."))
        import io
        buffer = io.BytesIO(); fig.savefig(buffer, format="png", dpi=140)
        display(Image(data=buffer.getvalue(), width=640))
        selected = [s for s in run["summaries"] if (s["source"],s["cohort"],s["structure"],s["domain"],s["endpoint"],s["unit"]) == key]
        for s in selected:
            median = "not available" if s["median"] is None else f"{s['median']:.4f}"
            display(Markdown(f"**Readout — {s['method']}:** median {median}; {s['n_present']} present, {s['n_missing']} missing."))
        display(Markdown("**Limit:** sequence/run linkage is unverified; saved draws are not independent biological replicates. This panel does not establish causality, clinical immunogenicity or binding affinity."))''')
for case in ('9NH7','9NFU'):
    md(f'<a id="{case.lower()}"></a>\n## {case} · read this case separately\n\nCompare only compatible measurements within their labelled source, cohort and domain. [Back to overview](#overview).')
    for unit,title in [('fraction','Humanness: where do the scores fall?'),('angstrom','Structure: how large are the reported deviations?')]:
        md(f'### {title}\n\nExpect source-specific distributions, then computed median/count readouts. No pooled ranking or paired test is performed.')
        code(f'render_sections(run["metrics"], "{case}", "{unit}")')
md('''<a id="next"></a>
## 4. What should we export or fix next?

The next cell saves a **new dated bundle**: full-precision tables, original input copies, hashes, all figures and the same-data reading pages. It stops if an input changed after loading. Old exports are never replaced.

Then obtain the scored/folded sequence manifest and structural provenance, resolve collection disagreements, and agree the cohort, endpoint and statistical design before a paper update. CSV arithmetic is not geometry validation. The latest Drive handoff and original paper-era loop-group ablation inputs remain unavailable here.''')
code('''bundle = update.export_run(run, EXPORT_FOLDER)
reading = write_reading(run, bundle)
# Include the reading surfaces in this new bundle's integrity manifest.
import hashlib
hashes = {str(p.relative_to(bundle)): hashlib.sha256(p.read_bytes()).hexdigest()
          for p in sorted(bundle.rglob("*")) if p.is_file() and p.name != "SHA256SUMS.json"}
(bundle / "SHA256SUMS.json").write_text(json.dumps(hashes, indent=2) + "\\n")
print("Saved:", bundle.relative_to(PROJECT_FOLDER))
print("Open READ_RESULTS.md or mobile.html in that folder.")
print("Figures:", reading["counts"]["figures"])
print("Previous input-byte comparison:", run["changed_inputs"])''')
md('''### Ready for review, not a final biological conclusion

Keep warnings with the plots when sharing. The machine-readable label `UNBOUND_SOURCE_EXPORT` means **the score has not been authenticated against its original sequence/run**; it is not a zero score or a failed arithmetic calculation. No publication occurs here.''')
nb=n.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}})
n.write(nb,root/'Update_results.ipynb')
