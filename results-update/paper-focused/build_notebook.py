"""Regenerate and execute the reader from a clean kernel, with no server."""
from pathlib import Path
import json
import sys
import nbformat as nbf
from nbclient import NotebookClient
from jupyter_client import KernelManager


def build(root):
    root=Path(root).resolve()
    md=nbf.v4.new_markdown_cell; code=nbf.v4.new_code_cell
    nb=nbf.v4.new_notebook(cells=[
        md('# Repertoire similarity and structural cost\n\n## Does adding AbLang improve repertoire similarity without a structural cost?\n\nRead the computed answer, Table 1 and the separate ablation Figure 1 below. This notebook regenerates the same [GitHub reading page](READ_RESULTS.md) and offline mobile HTML. It recalculates retained public outputs, not models or scoring databases.'),
        code("from pathlib import Path\nfrom IPython.display import display, Markdown, Image\nfrom analysis import publish, SCOPE, LIMIT, ABLATION, LOOP_CAPTION, PAPER_CHANGES\nROOT = Path.cwd()  # Open this notebook from the package folder.\nresult = publish(ROOT)\ndisplay(Markdown('**Answer.** ' + result['readout'] + '\\n\\n' + result['status']))"),
        md('## Table 1 · Two objectives, kept distinct\n\nHigher Heavy OASis Identity means greater repertoire similarity; higher H3 RMSD means greater displacement from the deposited reference. B = baseline; E = ensemble. These are marginal medians, not paired changes.'),
        code("display(result['table'])\ndisplay(Markdown(result['caption'] + '\\n\\n' + result['loop_readout'] + '\\n\\n' + SCOPE))"),
        md('## What differs between cases?\n\nRead the spread and prediction confidence alongside the median. No designs are deleted for low confidence.'),
        code("display(Markdown('\\n\\n'.join(result['cases']) + '\\n\\n' + LIMIT))"),
        md('## Figure 1 · Where does the humanness difference occur?\n\nThis is the **separate loop-group ablation**, not a decomposition of Table 1. Bars show mean ensemble-minus-baseline differences; no uncertainty intervals are implied.'),
        code("display(Image(filename=str(ROOT/'figures/figure1.png')))\ndisplay(Markdown(result['regions'] + '\\n\\n' + ABLATION + '\\n\\n' + result['light_note']))"),
        md('## Supporting loop table'),
        code("display(result['loops'].round(4))\ndisplay(Markdown(LOOP_CAPTION + '\\n\\n## Paper changes\\n\\n' + PAPER_CHANGES))"),
        md('## Read further only as needed\n\n[Supporting distributions and alternative endpoint notes](APPENDIX.md) · [Sources and scientific limits](SOURCES.md) · [Engineering receipts](engineering/README.md)\n\nScored sequences are now linked for the pinned selected main inputs; the old all-unbound characterization does not apply. This notebook rechecks scored VH/FASTA correspondence and retained-summary arithmetic. Prior coordinate replay is not re-executed here. Biological mapping acceptance, historical execution attribution, binding preservation and immunogenicity are not established by software checks.')
    ],metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}})
    nb.cells.append(md('## Update this analysis\n\nKeep the package together. Edit `selection.csv` in a spreadsheet to choose workbooks/FASTAs; structural and ablation CSVs live under `inputs/`. Preserve the schemas in [README](README.md). Restart and Run All, or run `python build_notebook.py`. Existing outputs are a snapshot until this succeeds. Changed bytes lose the original coordinate-audit qualification. Unsupported or missing files stop the update.'))
    for i,c in enumerate(nb.cells): c['id']=f'uirj-{i:02d}'
    target=root/'results.ipynb'
    # Save cleared output before execution so a failed update cannot retain old notebook results.
    nbf.write(nb,target)
    km=KernelManager(kernel_name='python3',transport='ipc',ip=str(root/'.k'),connection_file=str(root/'.kernel.json'))
    km.kernel_spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
    try:
        NotebookClient(nb,km=km,timeout=120,resources={'metadata':{'path':str(root)}}).execute(cleanup_kc=True)
    except Exception:
        for cell in nb.cells:
            if cell.cell_type=='code': cell.outputs=[]; cell.execution_count=None
        nbf.write(nb,target)
        raise
    nbf.write(nb,target)
    codes=[c for c in nb.cells if c.cell_type=='code']
    assert all(c.execution_count is not None for c in codes)
    assert not any(o.output_type=='error' for c in codes for o in c.outputs)
    images=sum('image/png' in o.get('data',{}) for c in codes for o in c.outputs)
    assert images>=1
    report={'clean_kernel':'PASS','executed_code_cells':len(codes),'embedded_png_outputs':images,'error_outputs':0}
    (root/'engineering/notebook.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__': build(Path(__file__).resolve().parent)
