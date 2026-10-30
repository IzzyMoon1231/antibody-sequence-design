"""Portable reading surfaces, derived only from the loaded run."""
import csv
from copy import deepcopy


def write_selection(config, path):
    rows=[dict(type=kind, **spec) for kind in ('fasta','designs','metrics') for spec in config[kind]]
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields); writer.writeheader(); writer.writerows(rows)


def load_selection(path, defaults):
    config=deepcopy(defaults)
    for kind in ('fasta','designs','metrics'): config[kind]=[]
    with path.open(newline='') as f:
        for row in csv.DictReader(f):
            kind=row.pop('type')
            if kind not in ('fasta','designs','metrics'): raise ValueError('type must be fasta, designs or metrics')
            spec={k:v for k,v in row.items() if v!=''}
            if kind=='metrics': spec.setdefault('domain_column',None)
            config[kind].append(spec)
    return config


def write_reading(run, folder):
    """Generate all views and exact-data captions from the same in-memory run."""
    from pathlib import Path
    from html import escape as e
    import base64, json
    from result_figures import build_figures
    folder=Path(folder); folder.mkdir(parents=True,exist_ok=True)
    (folder/'figures').mkdir(exist_ok=True)
    panels=[]
    keys=('source','cohort','structure','domain','endpoint','unit')
    for i,(key,fig) in enumerate(build_figures(run['metrics'])):
        selected=[s for s in run['summaries'] if tuple(s[k] for k in keys)==key]
        source,cohort,case,domain,endpoint,unit=key
        image=f'figures/{i:02d}-{case}.png'; fig.savefig(folder/image,dpi=140)
        parts=[]
        for s in selected:
            median='not available' if s['median'] is None else f"{s['median']:.4f}"
            parts.append(f"{s['method']}: median {median}, {s['n_present']} present, {s['n_missing']} missing")
        panels.append(dict(id=f'plot-{i:02d}',case=case,domain=domain,endpoint=endpoint,cohort=cohort,
            source=source,unit=unit,image=image,readout='; '.join(parts)+'.'))
    counts=dict(sequences=len(run['sequences']),measurements=len(run['metrics']),groups=len(run['summaries']),figures=len(panels),
        mismatches=sum(r['status']=='SEQUENCE_MISMATCH' for r in run['sequence_comparison']),
        possible_offsets=sum(r['possible_index_offset'] for r in run['coverage']))
    warning='Scores not yet linked to selected sequences. Provide the sequence file used for each score export, its run settings and structural mapping before using these as final paper results.'
    overview=f"{counts['sequences']:,} sequence records · {counts['measurements']:,} imported measurements · {counts['groups']} source-specific groups · {counts['figures']} plots. Measurements include the repeated workbook view, not extra independent observations."
    identity=f"{counts['mismatches']} same-ID sequence disagreements; {counts['possible_offsets']} source/domain groups with possible ID offsets. Ask which sequence set was actually scored; do not shift IDs or attach scores by row number."
    limits='Read each source, cohort, case, domain and endpoint separately. BioPhi source-selected score, Overview Identity and promb peptide content are not interchangeable. The 9NH7 ensemble workbook repeats a CSV population. No pairing, causal effect, clinical immunogenicity or binding claim is established. Framework and CDR RMSD are different endpoints; geometry is unvalidated.'
    css='''*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:90px}body{margin:0;background:#f4f6f8;color:#192d3d;font:17px/1.65 system-ui,sans-serif}main{max-width:960px;margin:auto;padding:24px}header{padding:32px 0 20px}h1{font-size:clamp(2rem,6vw,3.4rem);line-height:1.1;letter-spacing:-.04em}h2{line-height:1.25;margin-top:38px}h3{line-height:1.35}a{color:#155a87}nav{position:sticky;top:0;background:#fff;border-bottom:1px solid #ccd6df;padding:12px;display:flex;gap:16px;flex-wrap:wrap;z-index:2;font-size:15px}nav a{min-height:30px}section,figure,.card{background:white;border:1px solid #dae1e7;border-radius:16px;padding:22px;margin:18px 0}figure{margin:18px 0}img{display:block;width:100%;height:auto}small,.source{color:#526576;overflow-wrap:anywhere;font-size:14px}.warning{border-left:5px solid #b47115;background:#fff6e7;padding:18px;border-radius:8px}.eyebrow{color:#155a87;font-weight:700;letter-spacing:.12em;font-size:13px}.readout{border-left:3px solid #2379ab;padding-left:12px}summary{cursor:pointer;font-weight:700;padding:12px 0}li{margin:8px 0}details>summary{min-height:48px}footer{margin:32px 0;color:#526576}@media(max-width:500px){main{padding:14px}section,figure,.card{padding:14px}nav{gap:10px;font-size:14px}h2{font-size:24px}}'''
    html=[f'<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>UIRJ · Read the saved results</title><style>{css}</style><nav aria-label="Reading sections"><a href="#overview">Overview</a><a href="#checks">Checks</a>']
    cases=sorted({p['case'] for p in panels})
    html += [f'<a href="#{e(c)}">{e(c)}</a>' for c in cases]
    html += ['<a href="#next">Next steps</a></nav><main><header><div class="eyebrow">UIRJ / SAVED-RESULTS REVIEW</div><h1>What do these results actually tell us?</h1><p>A guided look at the supplied files. Descriptive evidence, not a validated method ranking.</p></header>',f'<div class="warning"><strong>{e(warning)}</strong></div>',f'<section id="overview"><h2>Start with the evidence boundary</h2><p>{e(overview)}</p><p>{e(limits)}</p><p><strong>Reader:</strong> check the warnings, then choose a case below. No Python needed.</p><p><strong>Updater:</strong> open Update_results.ipynb, edit the marked setup section and input_selection.csv, then restart and run all. These reading views are regenerated from that run.</p></section>',f'<section id="checks"><h2>Are the files and sequences right?</h2><p>{e(identity)}</p><ul>']
    html += [f'<li>{e(x)}</li>' for x in run['missing']]
    html += ['</ul></section>']
    md=['# What do these results actually tell us?','', '> **'+warning+'**','',overview,'',limits,'','**Reader:** use the case links below. **Updater:** open [the executed notebook](Update_results.ipynb) and edit the marked setup section. No JSON editing required.','', '## Are the files and sequences right?',identity,'']+['- '+x for x in run['missing']]
    for case in cases:
        md += ['',f'- [{case} details](#{case.lower()})']
    for case in cases:
        html += [f'<section id="{e(case)}"><h2>{e(case)} · What changes across the saved distributions?</h2><p>Choose humanness or structure. Each plot keeps its original source and cohort; do not pool panels.</p>']
        md += ['',f'## {case}','']
        for unit,title in [('fraction','Humanness: where do the saved scores fall?'),('angstrom','Structure: how large are the reported deviations?')]:
            group=[p for p in panels if p['case']==case and p['unit']==unit]
            if not group: continue
            html += [f'<details><summary>{title} · {len(group)} plots</summary>']
            md += [f'### {title}','']
            for p in group:
                question=f"{p['domain']} — {p['endpoint']}"
                axes=('Horizontal: saved fraction (0–1). Rightward means larger scores, not lower clinical immunogenicity.' if unit=='fraction' else 'Horizontal: reported RMSD (Å). Leftward means smaller deviations under this export’s alignment, not better binding.')+' Vertical: fraction of nonmissing values at or below x. Read the median at height 0.5; legends give present and missing counts.'
                limitation='Unpaired saved draws, not independent biological replicates. Sequence/run linkage unverified. Shared endpoint scales can make narrow distributions look vertical; tails are not clipped.'
                data=base64.b64encode((folder/p['image']).read_bytes()).decode()
                html += [f'<figure id="{p["id"]}"><h3>{e(question)}</h3><p class="source">{e(p["cohort"])}<br>{e(p["source"])}</p><p>{e(axes)}</p><a href="{p["image"]}" aria-label="Open full size plot"><img alt="{e(question)}" src="data:image/png;base64,{data}"></a><figcaption><p class="readout"><strong>Readout:</strong> {e(p["readout"])}</p><p><strong>Limit:</strong> {e(limitation)}</p></figcaption><a href="#'+case+'">Back to case</a></figure>']
                md += [f'#### {question}',f'**Cohort:** {p["cohort"]}  ',f'**Source:** `{p["source"]}`','',axes,'',f'![Saved distribution]({p["image"]})','',f'**Readout:** {p["readout"]}','',f'**Limit:** {limitation}','']
            html += ['</details>']
        html += ['</section>']
    next_step='Obtain the actual scored/folded sequence set and run manifest; resolve conflicting IDs without guessing; agree the cohort, endpoint and statistical design. Re-run this notebook to review a new dated bundle. No model inference or publication occurs here.'
    html += [f'<section id="next"><h2>What should we fix or export?</h2><p>{next_step}</p><p>Full precision tables and source hashes are in the dated notebook export. Captions round display medians only; raw data are unchanged.</p></section><footer>Offline reading copy · same loaded data as the notebook · no tracking, external fonts or network requests.</footer></main></html>']
    md += ['## What should we fix or export?',next_step,'','[Executed notebook](Update_results.ipynb) · [Update guide](README.md) · Download mobile.html for an offline phone-friendly view. GitHub shows HTML source, not a hosted page.']
    (folder/'mobile.html').write_text(''.join(html))
    (folder/'READ_RESULTS.md').write_text('\n'.join(md)+'\n')
    report=dict(counts=counts,panels=panels,input_hashes=run['input_hashes'])
    (folder/'reading-data.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
