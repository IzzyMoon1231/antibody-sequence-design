"""Build the teaching notebook; execute with no server or analysis-file exports.
The CLI writes only results.ipynb and a new engineering execution receipt.
"""
from pathlib import Path
import json
import sys
from analysis import ORIENTATION
import tempfile
import nbformat as nbf
from nbclient import NotebookClient
from jupyter_client import KernelManager


def notebook():
    md=nbf.v4.new_markdown_cell; code=nbf.v4.new_code_cell
    cells=[
        md('# Repertoire similarity and structural cost\n\n**Question:** Do ensemble sequences score more like human antibodies, and what happens to predicted structure?\n\n**Reader:** outputs below are from the included real public results, not synthetic demonstrations. [Compact GitHub reader](READ_RESULTS.md). **Student:** run each step, inspect the rows, then explain the result in your own words. You will choose files → inspect → validate → calculate → plot → interpret → optionally export. Nothing is scored or folded here. [Student guide and schemas](README.md).' + '\n\n' + ORIENTATION),
        md('## 1 · Choose local inputs\n\nFirst run: leave these paths alone. For your own compatible results, copy files into `inputs/my_run/` using Finder/Explorer or a **local** JupyterLab file browser Upload button. In VS Code use the local folder explorer. Do not use a remote notebook server or cloud upload. This code never transmits files.\n\nEdit a copy of `selection.csv`: one workbook and matching FASTA for each case/method. Paths inside that CSV are relative to this notebook folder, **not** to the CSV location. Select structural and ablation CSVs below too. Set `ABLATION = None` to omit the separate ablation; the main analysis requires all measured loops. Only the two supported cases and existing endpoint conventions are supported. A new metric needs a mentor-reviewed design, not relabelled columns.'),
        code("from pathlib import Path\nimport sys\nsys.dont_write_bytecode = True\nimport pandas as pd\nimport matplotlib.pyplot as plt\nfrom IPython.display import display, Markdown\nimport analysis\nfrom workflow import read_inputs, validate_inputs, status, export_bundle\nROOT = Path.cwd()  # Open this notebook in its package folder.\nSELECTION = 'selection.csv'\nSTRUCTURAL = 'inputs/cdr_rmsd_results.csv'\nABLATION = 'inputs/ablation_scores.csv'  # or None; never silently substituted"),
        md('## 2 · What did we actually select?\n\nRun this cell to read the files, not compute results. Check the paths, antibody IDs, units and first rows. BioPhi’s heavy-domain Identity is a fraction from 0 to 1, including for scFv. The structural table has one row per design **and loop**. The original file also contains AbMPNN; the comparison explicitly excludes that arm. Do these files describe the same main designs?'),
        code("data = None  # invalidate the previous in-memory analysis before a new read\nloaded = None\nloaded = read_inputs(ROOT, SELECTION, STRUCTURAL, ABLATION)\ndisplay(loaded['selection'])\nfor key, sheets in loaded['books'].items():\n    print(key, 'Overview — first 3 actual rows')\n    display(sheets['Overview'][['Antibody', 'Heavy OASis Identity']].head(3))\ndisplay(loaded['structural'].head(5))\nif loaded['ablation'] is not None:\n    display(loaded['ablation'].head(4))"),
        md('## 3 · Do the identities and measurements pass the checks?\n\nThese checks ask whether the selected scores, sequences and structural measurements belong together and can be compared without missing or invalid values. Inspect the counts below before continuing. See [What is checked](README.md#what-is-checked) for identity rules, ranges and pairing mechanics.\n\n**If a check fails:** stop; read the named file/row/key, correct the source export or selected path, then restart and run from step 1. Do not delete difficult rows, fill missing measurements, or use old visible plots. Clear notebook outputs after a failed interactive run. Previously exported folders are historical and stay untouched. Passing these checks is computational compatibility, not certification of the experimental protocol.'),
        code("data = None\ndata = validate_inputs(loaded)\ndisplay(Markdown(status(data)))\ndisplay(data['main'].groupby(['structure','method']).agg(\n    draws=('index','size'), heavy_median=('oasis','median'),\n    H3_below70=('low_confidence','sum')))\nprint('Scored VH/FASTA and complete-loop checks passed. Main arms remain unpaired.')"),
        md('### Experiment map · two questions, separate cohorts\n\nCounts below come from the validated selection. Read the main comparison first; the ablation is not its next stage.'),
        code("map_figure = analysis.plot_experiment_map(data)\ndisplay(map_figure)\nplt.close(map_figure)"),
        md('## 4 · Predict, then calculate\n\n**Before running:** predict which arm has the higher heavy-domain score in each case. Does a higher similarity score force lower RMSD?\n\nThe code first groups actual observations to show medians and counts. A delta is **ensemble median − baseline median**, not the median of paired differences. All draws (including repeated sequences and low-confidence predictions) keep their weight. These are descriptive summaries, not confidence intervals or significance tests.'),
        code("main = data['main']\ngrouped = main.groupby(['structure','method'])[['oasis','rmsd_angstrom']].median()\ndisplay(grouped.rename(columns={'oasis':'Heavy Identity median','rmsd_angstrom':'H3 median (Å)'}))\ndeltas = grouped.xs('ensemble',level='method') - grouped.xs('baseline',level='method')\ndisplay(deltas.rename(columns={'oasis':'Δ marginal Identity median','rmsd_angstrom':'Δ marginal H3 median (Å)'}))\ntable, ablation_table = analysis.summarize(data)  # same formulas + separate domain denominators\nassert table.set_index('structure').biophi_heavy_difference_medians.sort_index().equals(deltas.oasis.sort_index())\nloops, loop_readout = analysis.loop_summary(data)\ndisplay(analysis.table_display(table))"),
        md('## 5 · Four views of the main comparison\n\nEach panel is drawn from the loaded frames, not a PNG. Colors and scales agree across cases. Structural dots are medians; lines show the middle 50%, not confidence intervals. Main arms remain unpaired. The optional export also creates a compact 2×2 manuscript figure: cases in columns, metrics in rows.'),
        md('### 9NH7 · Does ensemble increase heavy-domain repertoire similarity?'),
        code("panel_9NH7_similarity = analysis.plot_panel(data, '9NH7', 'similarity')\ndisplay(panel_9NH7_similarity)\nplt.close(panel_9NH7_similarity)"),
        code("display(Markdown(analysis.panel_readout(data, '9NH7', 'similarity')))"),
        md('### 9NH7 · What structural trade-off accompanies the similarity result?'),
        code("panel_9NH7_structure = analysis.plot_panel(data, '9NH7', 'structure')\ndisplay(panel_9NH7_structure)\nplt.close(panel_9NH7_structure)"),
        code("display(Markdown(analysis.panel_readout(data, '9NH7', 'structure')))"),
        md('### 9NFU · Does the heavy-only similarity pattern recur in the scFv case?'),
        code("panel_9NFU_similarity = analysis.plot_panel(data, '9NFU', 'similarity')\ndisplay(panel_9NFU_similarity)\nplt.close(panel_9NFU_similarity)"),
        code("display(Markdown(analysis.panel_readout(data, '9NFU', 'similarity')))"),
        md('### 9NFU · What does H3 miss about the other heavy and light loops?'),
        code("panel_9NFU_structure = analysis.plot_panel(data, '9NFU', 'structure')\ndisplay(panel_9NFU_structure)\nplt.close(panel_9NFU_structure)"),
        code("display(Markdown(analysis.panel_readout(data, '9NFU', 'structure')))"),
        code("readout, case_details, region_readout = analysis.describe(data, table, ablation_table)"),
        md('### Check your prediction\n\nWrite two sentences: (1) a numerical observation for one case, (2) a limitation that changes its interpretation. Compare H3 with the other loops before calling a structural cost small.\n\n**Worked example for the included dataset only:** “The 9NH7 ensemble heavy-domain median is higher, but its H3 median RMSD is also higher. This is a similarity/structure trade-off in retained predictions, not evidence of preserved binding.” Verify both directions against your current table; do **not** reuse this sentence if your values differ.\n\n**Hypothesis:** a protocol difference might contribute to the score shift; these summaries do not identify the cause. **Unestablished claims:** reduced immunogenicity, preserved binding, structural non-inferiority, a general antibody-format effect, or AbLang-only causation.'),
        md('## 6 · Supporting evidence: spread and confidence\n\nRead every loop, not only a favorable endpoint. Mean loop pLDDT < 70 flags low confidence; it does not establish why a structure changed and does not remove a draw. Domain-local fits omit relative VH/VL orientation. Mapping is operational/motif-based, not independently validated biological numbering.'),
        code("display(loops.round(4))\nfor detail in case_details:\n    display(Markdown(detail))"),
        md('## 7 · Separate ablation: calculate reported pair differences\n\nSee the experiment map after validation: this is another cohort, not a pairing of the main designs. Here the key is case + condition + reported seed. Look at the explicit join and subtraction below. `others` is H1/H2 for 9NH7, but H1/H2 **plus L1/L2/L3** for 9NFU; the latter is not an isolated heavy-loop intervention. Heavy and light endpoints are kept separate.\n\n**Predict:** should every pair have the same sign just because the mean is positive? Check rather than assume.'),
        code("ablation = data['ablation']\nif not ablation.empty:\n    pairs = ablation.pivot(index=['structure','condition','seed'], columns='method', values='Heavy OASis Identity')\n    pairs['difference_pp'] = 100 * (pairs['ensemble'] - pairs['baseline'])\n    display(pairs.head(6))\n    display(pairs.groupby(['structure','condition']).difference_pp.agg(['size','mean','min','max']))\n    display(ablation_table[['structure','condition','domain','n_pairs','mean_difference','positive','negative','zero']])\nelse:\n    print('Ablation not selected: no pair differences or causal interpretation computed.')"),
        md('## 8 · Figure 2: inspect the distribution, not just its mean\n\nDots are reported seed-pair differences in percentage points; diamonds are means; thick lines are IQRs. A positive value favors ensemble similarity. The source-defined zero-logit control’s executed helper/weights remain unauthenticated for the original dataset; new files need their own protocol evidence. No matched structural ablation measurements are supplied.'),
        code("if not ablation.empty:\n    ablation_figure = analysis.plot_ablation(ablation_table, ablation)\n    display(ablation_figure)\n    plt.close(ablation_figure)\n    display(Markdown('**Computed observations:** ' + region_readout))\nelse:\n    ablation_figure = None\n    print('Figure 2 omitted — no selected ablation.')"),
        md('## 9 · Explain the result yourself\n\nReplace the notes below before exporting. Use values from your selected tables, not an AI verdict.\n\n- **Observation:** which case, endpoint, contrast and draw count? Name an adverse or heterogeneous result too.\n- **Check:** do plotted dots/medians agree with the small tables? Are confidence counts plausible? Is any endpoint being confused with another?\n- **Hypothesis:** one possible explanation and what evidence would discriminate it.\n- **Limit:** what remains unestablished? What would you ask a mentor before designing a new metric or inferential test?\n- **Reflection:** did your prediction hold? What changed your mind?'),
        code("student_notes = {\n    'observation': '',\n    'sanity_check': '',\n    'hypothesis_and_test': '',\n    'limitation_and_mentor_question': '',\n    'reflection': '',\n}\ndisplay(Markdown(analysis.LIMIT))"),
        md('## 10 · Optional export — only when you choose\n\nFirst **Restart Kernel and Run All** with export off; compare your outputs, complete your notes and save the notebook. Then set `EXPORT = True` below and run this cell. It writes tables, figures, same-source GitHub/mobile readers, exact selected input bytes, hashes and versions into a **new UTC-dated folder** under `exports/`. No old export or included example is overwritten. Keep your saved executed notebook with it. Reset the switch to `False` before the next restart. A software check is not scientific endorsement; review privacy before sharing your own inputs.'),
        code("EXPORT = False\nif EXPORT:\n    if not all(value.strip() for value in student_notes.values()):\n        raise ValueError('Complete the five reasoning notes before exporting; use your own words.')\n    export_path = export_bundle(ROOT, data, confirm=True)\n    import json\n    (export_path / 'student-notes.json').write_text(json.dumps(student_notes, indent=2))\n    print('New local export:', export_path.relative_to(ROOT))\nelse:\n    print('Export is off. No analysis outputs were written.')"),
    ]
    nb=nbf.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}})
    for i,c in enumerate(nb.cells): c['id']=f'uirj-step-{i:02d}'
    return nb


def execute(nb, root):
    root=Path(root).resolve()
    # Short IPC names avoid macOS Unix socket path limits; temporary files stay task-owned.
    with tempfile.TemporaryDirectory(prefix='.k-',dir=Path(__file__).parent) as tmp:
        km=KernelManager(kernel_name='python3',transport='ipc',ip=str(Path(tmp)/'s'),connection_file=str(Path(tmp)/'k.json'))
        km.kernel_spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
        NotebookClient(nb,km=km,timeout=120,resources={'metadata':{'path':str(root)}}).execute(cleanup_kc=True)
    return nb


def build(root):
    root=Path(root).resolve()
    nb=execute(notebook(),root)
    # Write only after successful fresh-kernel execution. Existing snapshots survive a failure.
    nbf.write(nb,root/'results.ipynb')
    codes=[c for c in nb.cells if c.cell_type=='code']
    report={'clean_kernel':'PASS','executed_code_cells':len(codes),'embedded_png_outputs':sum('image/png' in o.get('data',{}) for c in codes for o in c.outputs),'error_outputs':0,'analysis_export_enabled':False}
    assert all(c.execution_count is not None for c in codes)
    assert not any(o.output_type=='error' for c in codes for o in c.outputs)
    (root/'engineering/panels-notebook.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__': build(Path(__file__).resolve().parent)
