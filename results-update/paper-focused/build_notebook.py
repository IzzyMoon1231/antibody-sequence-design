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
        md('# Repertoire similarity and structural cost\n\nDo ensemble sequences look more like human antibodies, and what happens to predicted structure? Read the answer and main figure first; the ablation is a separate question. This notebook recalculates retained public outputs, not models or scoring databases. [GitHub reader](READ_RESULTS.md) · [Update guide](README.md).'),
        code("from pathlib import Path\nfrom IPython.display import display, Markdown, Image\nfrom analysis import publish, LIMIT, LOOP_CAPTION\nROOT = Path.cwd()  # Keep the package together.\nresult = publish(ROOT)\ndisplay(Markdown(result['readout'] + '\\n\\n' + result['status']))"),
        md('## Which experiment answers which question?'),
        code("display(Markdown(result['map_caption']))\ndisplay(Image(filename=str(ROOT/'figures/experiment_map.png')))"),
        md('## Figure 1 · Sequence similarity and every measured loop\n\nRead each case separately. Light sequence dots retain every draw. Thick lines show the middle half, not uncertainty in the median. No main designs are paired across methods.'),
        code("display(Image(filename=str(ROOT/'figures/figure1.png')))\ndisplay(Markdown(result['main_caption'] + '\\n\\n' + result['loop_readout']))"),
        md('## Table 1 · Numerical reference'),
        code("display(result['table'])\ndisplay(Markdown(result['caption']))"),
        md('## Figure 2 · Where are differences seen in the separate ablation?\n\nHere a dot is a reported seed-pair difference, not a main-cohort design. Diamonds are means; lines are IQRs. The zero line means no difference.'),
        code("display(Image(filename=str(ROOT/'figures/figure2.png')))\ndisplay(Markdown(result['ablation_caption'] + '\\n\\n' + result['regions'] + '\\n\\n' + result['light_note']))"),
        md('## What can we conclude?'),
        code("display(Markdown(LIMIT))"),
        md('## Details and next steps\n\n[Spread, prediction confidence and complete loop table](APPENDIX.md) · [Sources and measurement conventions](SOURCES.md) · [Author checklist](AUTHOR_CHECKLIST.md) · [Engineering receipts](engineering/README.md)\n\nSelected original sequences are linked; prior coordinate replay is not re-executed here. Software checks do not establish binding preservation, immunogenicity or authenticated causal execution.\n\nTo update, edit `selection.csv` for matching workbooks/FASTAs and preserve the CSV schemas in [README](README.md). Restart and Run All, then save; or run `python build_notebook.py`. Changed bytes lose the original coordinate-audit qualification. Missing inputs stop the update.')
    ],metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}})
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
    assert images==3
    report={'clean_kernel':'PASS','executed_code_cells':len(codes),'embedded_png_outputs':images,'error_outputs':0}
    (root/'engineering/notebook.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__': build(Path(__file__).resolve().parent)
