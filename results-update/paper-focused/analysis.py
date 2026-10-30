"""Offline descriptive analysis; bounded refactor of inspected analyze.py.
No upstream imports, generation, scoring, network access or inferential tests.
"""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd

CASES = ('9NH7', '9NFU')
METHODS = ('baseline', 'ensemble')


def fasta(path):
    records = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            records.append([line[1:].split()[0], ''])
        elif line.strip():
            if not records:
                raise ValueError(f'{path.name}: FASTA sequence precedes header')
            records[-1][1] += line.strip()
    return records


def load(root):
    root = Path(root)
    selection = pd.read_csv(root/'selection.csv')
    if set(selection.columns) != {'structure', 'method', 'workbook', 'fasta'}:
        raise ValueError('selection.csv: require structure, method, workbook, fasta columns')
    if set(zip(selection.structure, selection.method)) != {(s,m) for s in CASES for m in METHODS} or len(selection)!=4:
        raise ValueError('selection.csv: select exactly one baseline and ensemble for each supported case')
    paths = list(selection.workbook)+list(selection.fasta)+['inputs/cdr_rmsd_results.csv','inputs/ablation_scores.csv']
    for p in paths:
        if Path(p).is_absolute() or not (root/p).resolve().is_relative_to(root.resolve()):
            raise ValueError('Input paths must stay inside the package')
        if not (root/p).is_file():
            raise ValueError(f'Missing input: {p}. Restore it or correct selection.csv; no results regenerated.')
    rows = []
    for r in selection.itertuples():
        overview = pd.read_excel(root/r.workbook, sheet_name='Overview')
        vh = pd.read_excel(root/r.workbook, sheet_name='VH')
        seqs = {str(row['Antibody']): ''.join(str(row[c]) for c in vh.columns if c not in ['Antibody','scheme','chain_type'] and pd.notna(row[c])).replace('-','').replace(' ','') for _,row in vh.iterrows()}
        overview = overview[overview.Antibody.str.endswith('_VH')] if overview.Antibody.str.endswith('_VH').any() else overview
        sequences = fasta(root/r.fasta)
        for _, row in overview.iterrows():
            match = re.search(r'(\d+)(?:_VH)?$',str(row.Antibody))
            if not match:
                raise ValueError(f'{r.workbook}: unsupported Antibody ID {row.Antibody}')
            idx = int(match[1]) - (r.method == 'baseline')
            seq = seqs.get(row.Antibody)
            if idx < 0 or idx >= len(sequences) or not seq or seq not in sequences[idx][1]:
                raise ValueError(f'{r.workbook}: VH/FASTA mismatch at {row.Antibody}; correct sequence linkage first')
            rows.append(dict(structure=r.structure,method=r.method,index=idx,oasis=row['Heavy OASis Identity']))
    bio = pd.DataFrame(rows)
    if bio.duplicated(['structure','method','index']).any():
        raise ValueError('BioPhi: duplicate case/method/index')
    for _,g in bio.groupby(['structure','method']):
        if set(g['index']) != set(range(100)):
            raise ValueError('BioPhi: this selected-cohort reader requires exactly indices 0–99 per arm')
    raw = pd.read_csv(root/'inputs/cdr_rmsd_results.csv')
    required={'cdr','method','structure','design_id','domain','rmsd_angstrom','n_atoms_per_residue','cdr_mean_plddt','low_confidence','framework_rmsd_angstrom'}
    if required-set(raw.columns):
        raise ValueError('cdr_rmsd_results.csv: missing columns: '+', '.join(sorted(required-set(raw.columns))))
    h3 = raw[(raw.cdr=='H3') & raw.method.isin(METHODS) & raw.structure.isin(CASES)].copy()
    h3['index'] = h3.design_id.str.extract(r'(\d+)$').astype(int)
    main = bio.merge(h3,on=['structure','method','index'],how='outer',validate='one_to_one',indicator=True)
    if not main._merge.eq('both').all():
        raise ValueError('H3 metrics and scored VH must have identical case/method/index keys')
    if not main.n_atoms_per_residue.eq(3).all() or not main.domain.eq(main.structure.map({'9NH7':'VHH','9NFU':'VH'})).all():
        raise ValueError('Unsupported H3 endpoint: require per-domain N/CA/C backbone, VHH or VH')
    for col,lo,hi in [('oasis',0,1),('rmsd_angstrom',0,np.inf),('cdr_mean_plddt',0,100),('framework_rmsd_angstrom',0,np.inf)]:
        if not np.isfinite(main[col]).all() or not main[col].between(lo,hi).all():
            raise ValueError(f'Invalid {col}: require finite values in [{lo}, {hi}]')
    if not main.low_confidence.eq((main.cdr_mean_plddt<70).astype(int)).all():
        raise ValueError('low_confidence must match H3 mean pLDDT < 70')
    a = pd.read_csv(root/'inputs/ablation_scores.csv')
    meta = a.id.str.extract(r'^(?P<structure>9NH7|9NFU)_(?P<condition>cdrh3|others)_(?P<method>baseline|ensemble)_s(?P<seed>\d+)$')
    if not meta.notna().all().all() or a.id.duplicated().any():
        raise ValueError('Ablation: require unique case_condition_method_sN identifiers')
    a = pd.concat([a,meta],axis=1)
    for st in CASES:
        for condition in ('cdrh3','others'):
            g = a[(a.structure==st)&(a.condition==condition)]
            w = g.pivot(index='seed',columns='method',values='Heavy OASis Identity')
            if set(w.columns)!=set(METHODS) or set(w.index)!=set(map(str,range(100))) or not np.isfinite(w.to_numpy()).all() or not ((w>=0)&(w<=1)).all().all():
                raise ValueError(f'Ablation {st}/{condition}: require 100 complete finite reported seed pairs, scores 0–1')
    original = {x['path']:x['sha256'] for x in json.loads((root/'inputs/provenance.json').read_text())}
    hashes = {p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths}
    return dict(main=main.drop(columns='_merge'), loops=raw[raw.method.isin(METHODS) & raw.structure.isin(CASES)].copy(), ablation=a, hashes=hashes,
                original_inputs=all(original.get(p)==v for p,v in hashes.items()))


def summarize(data):
    tables=[]
    for st in CASES:
        g=data['main'].query('structure==@st')
        b=g[g.method=='baseline']; e=g[g.method=='ensemble']
        tables.append(dict(structure=st,humanness_n_per_arm=len(b),h3_n_per_arm=len(b),
            biophi_heavy_baseline_median=b.oasis.median(),biophi_heavy_ensemble_median=e.oasis.median(),
            biophi_heavy_difference_medians=e.oasis.median()-b.oasis.median(),
            h3_baseline_median_A=b.rmsd_angstrom.median(),h3_ensemble_median_A=e.rmsd_angstrom.median(),
            h3_difference_medians_A=e.rmsd_angstrom.median()-b.rmsd_angstrom.median()))
    abl=[]
    for (st,cond),g in data['ablation'].groupby(['structure','condition']):
        for domain,col,denom,npos in [('H','Heavy OASis Identity',107 if st=='9NH7' else 113,(15 if st=='9NH7' else 14) if cond=='cdrh3' else 16),('L','Light OASis Identity',96,0 if cond=='cdrh3' else 17)]:
            if st=='9NH7' and domain=='L': continue
            w=g.pivot(index='seed',columns='method',values=col)
            if not np.isfinite(w.to_numpy()).all() or not ((w>=0)&(w<=1)).all().all():
                raise ValueError(f'Ablation {st}/{cond}/{domain}: missing or invalid scores')
            if np.max(np.abs(w.to_numpy()*denom-np.rint(w.to_numpy()*denom)))>=1e-8:
                raise ValueError(f'Ablation {st}/{cond}/{domain}: scores do not match the declared peptide denominator {denom}')
            d=w.ensemble-w.baseline
            abl.append(dict(structure=st,condition=cond,domain=domain,n_pairs=len(w),designed_domain_residues=npos,peptide_denominator=denom,baseline_mean=w.baseline.mean(),ensemble_mean=w.ensemble.mean(),mean_difference=d.mean(),mean_gain_percentage_points=100*d.mean(),gain_per_designed_residue=d.mean()/npos if npos else np.nan,mean_additional_qualifying_9mers=d.mean()*denom,positive=int((d>0).sum()),negative=int((d<0).sum()),zero=int((d==0).sum())))
    return pd.DataFrame(tables),pd.DataFrame(abl)


def direction(value):
    return 'higher' if value > 1e-12 else 'lower' if value < -1e-12 else 'unchanged'


def describe(data, table, ablation):
    answers=[]; cases=[]; regions=[]
    for t in table.itertuples():
        answers.append(f'{t.structure}: {direction(t.biophi_heavy_difference_medians)} heavy-domain repertoire similarity and {direction(t.h3_difference_medians_A)} H3 backbone RMSD in the ensemble median.')
        g=data['main'].query('structure==@t.structure')
        b=g[g.method=='baseline']; e=g[g.method=='ensemble']
        bq=b.rmsd_angstrom.quantile([.25,.75]); eq=e.rmsd_angstrom.quantile([.25,.75])
        overlap=max(bq.iloc[0],eq.iloc[0]) <= min(bq.iloc[1],eq.iloc[1])
        cases.append(f"**{t.structure}.** H3 interquartile ranges are {bq.iloc[0]:.3f}–{bq.iloc[1]:.3f} Å (baseline) and {eq.iloc[0]:.3f}–{eq.iloc[1]:.3f} Å (ensemble); these ranges {'overlap' if overlap else 'do not overlap'}. The ensemble middle-half spread is {direction((eq.iloc[1]-eq.iloc[0])-(bq.iloc[1]-bq.iloc[0]))}. Mean H3 RMSD is {b.rmsd_angstrom.mean():.5f} → {e.rmsd_angstrom.mean():.5f} Å; framework medians are {b.framework_rmsd_angstrom.median():.3f} → {e.framework_rmsd_angstrom.median():.3f} Å. H3 mean pLDDT is below 70 in {int(b.low_confidence.sum())}/{len(b)} → {int(e.low_confidence.sum())}/{len(e)} designs. All remain in the primary summaries; low confidence weakens interpretation as a true fold, not the retained arithmetic.")
        h=ablation[(ablation.structure==t.structure)&(ablation.domain=='H')].set_index('condition')
        comparison=direction(h.loc['others','mean_difference']-h.loc['cdrh3','mean_difference'])
        regions.append(f"{t.structure}: the mean difference outside H3 is {comparison} than H3-only ({h.loc['others','mean_difference']:+.5f} versus {h.loc['cdrh3','mean_difference']:+.5f}); positive reported pairs: {int(h.loc['others','positive'])}/{int(h.loc['others','n_pairs'])} versus {int(h.loc['cdrh3','positive'])}/{int(h.loc['cdrh3','n_pairs'])}.")
    if (table.h3_difference_medians_A > 1e-12).any():
        answers.append('The retained median pattern does not support a blanket no-structural-cost claim; the structural shift must be read case by case, alongside spread and confidence.')
        if (table.h3_difference_medians_A > 1e-12).all():
            largest=table.loc[table.h3_difference_medians_A.idxmax()]
            smallest=table.loc[table.h3_difference_medians_A.idxmin()]
            if largest.structure != smallest.structure:
                answers.append(f'The H3 median increase is larger for {largest.structure} ({largest.h3_difference_medians_A:+.4f} Å) than {smallest.structure} ({smallest.h3_difference_medians_A:+.4f} Å); this is an observed scale difference, not a biological-significance threshold.')
    else:
        answers.append('H3 medians do not increase in these inputs; this alone does not establish structural preservation or non-inferiority.')
    return ' '.join(answers),cases,' '.join(regions)


def loop_summary(data):
    rows=[]
    for st in CASES:
        for cdr in ('H1','H2','H3') + (('L1','L2','L3') if st=='9NFU' else ()):
            g=data['loops'].query('structure==@st and cdr==@cdr')
            row=dict(structure=st,cdr=cdr)
            for method in METHODS:
                arm=g[g.method==method]
                if len(arm)!=100 or arm.design_id.nunique()!=100:
                    raise ValueError('Complete loops require 100 unique designs per arm')
                if not np.isfinite(arm.rmsd_angstrom).all() or not arm.rmsd_angstrom.ge(0).all():
                    raise ValueError('Invalid loop RMSD')
                row[method+'_n']=len(arm)
                row[method+'_median_A']=float(arm.rmsd_angstrom.median())
                row[method+'_below70']=int((arm.cdr_mean_plddt<70).sum())
            rows.append(row)
    loops=pd.DataFrame(rows)
    z=loops[loops.structure=='9NFU'].set_index('cdr')
    changes=z.ensemble_median_A-z.baseline_median_A
    # Highlight non-H3 changes larger than the focused H3 shift, not a significance margin.
    larger=', '.join(c for c in changes.index if c!='H3' and changes[c]>max(0,changes['H3'])+1e-12) or 'none'
    improved=', '.join(c for c in changes.index if changes[c]<-1e-12) or 'none'
    text=f'Across 9NFU loops, median RMSD increases larger than its H3 shift occur in {larger}; decreases (improvements relative to the reference) occur in {improved}. The H3 result is not overall scFv preservation. The complete supporting loop table retains all draws. Low confidence means mean loop pLDDT below 70, not any low-confidence residue, and does not explain away every worsening: 9NFU H2 below-70 counts are {int(z.loc["H2","baseline_below70"])} baseline and {int(z.loc["H2","ensemble_below70"])} ensemble.'
    return loops,text


LOOP_CAPTION='All selected baseline/ensemble loops, marginal medians in Å, with draw counts and mean loop pLDDT below-70 counts. No confidence filtering is applied. N/CA/C RMSD follows motif-defined spans after separate domain framework fits to deposited reference chains E (9NH7) and C (9NFU); CDR atoms are excluded from fitting. This does not measure relative VH/VL orientation or whole-scFv preservation.'
PAPER_CHANGES='Use Table 1 and the complete supporting loop table to replace stale humanness and loop values; retain the heterogeneous loop-level structural costs, not universal preservation, foldability or therapeutic claims. Describe deposited-reference, motif-defined N/CA/C domain fits and whole heavy-domain Identity consistently. Keep Figure 1 descriptive, with the joint heavy/light other-CDR condition and unchanged-light control; remove per-residue or cross-chain efficiency claims. Remove index-paired Wilcoxon and inherited significance claims, completed contact-preservation claims, and structural-ablation comparisons without matching outputs. H3 is a focused endpoint, not an authenticated prospective designation. Two fixed cases are not a format ablation. These are update instructions, not certification that a manuscript has been edited.'


def plot_ablation(ablation):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'svg.hashsalt':'uirj-paper-focused','svg.fonttype':'none'})
    fig,axes=plt.subplots(2,1,figsize=(6.2,7.2),layout='constrained')
    all_values=ablation[ablation.domain=='H'].mean_gain_percentage_points
    low=min(0,float(all_values.min())); high=max(0,float(all_values.max())); pad=max(high-low,1)*.22
    for ax,st in zip(axes,CASES):
        z=ablation[(ablation.structure==st)&(ablation.domain=='H')].set_index('condition').loc[['cdrh3','others']]
        values=z.mean_gain_percentage_points.to_numpy()
        bars=ax.bar(['H3 only','H1/H2' if st=='9NH7' else 'H1/H2 + L1/L2/L3'],values,color=['#287c8e','#bc643f'],width=.55)
        ax.axhline(0,color='#555555',linewidth=.7)
        ax.bar_label(bars,labels=[f'{v:+.2f} pp' for v in values],padding=5)
        ax.set_ylim(low-pad if low < 0 else 0,high+pad if high > 0 else 0)
        ax.set_title(st+' · heavy domain',loc='left',weight='bold')
        ax.set_ylabel('Heavy OASis Identity difference\n(percentage points)')
        ax.spines[['right','top']].set_visible(False)
    return fig


def table_display(table):
    return pd.DataFrame([{'Case':r.structure,'Heavy-domain Identity, B → E':f'{r.biophi_heavy_baseline_median:.4f} → {r.biophi_heavy_ensemble_median:.4f}', 'Δ Identity':f'{r.biophi_heavy_difference_medians:+.4f}', 'H3 RMSD, B → E (Å)':f'{r.h3_baseline_median_A:.4f} → {r.h3_ensemble_median_A:.4f}', 'Δ RMSD (Å)':f'{r.h3_difference_medians_A:+.4f}'} for r in table.itertuples()])


def markdown_table(df):
    return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in df.itertuples(index=False,name=None))


QUESTION='Does adding AbLang improve human-repertoire similarity without a structural cost?'
SCOPE='Table 1 compares independently sampled all-loop baseline and ensemble protocol cohorts, not paired seeds. Heavy OASis Identity is the whole numbered heavy-domain fraction of qualifying overlapping 9-mers at the reported relaxed threshold: not an OASis percentile, isolated CDR score, or immunogenicity measurement. H3 RMSD uses a per-domain framework N/CA/C backbone fit to deposited experimental coordinates, not the original RFdiffusion model. Loops are operational sequence-alignment/motif definitions, not independently validated biological numbering. Separate domain fits exclude relative VH/VL orientation and cannot establish whole-scFv orientation or interface preservation.'
LIMIT='Neither overlapping ranges nor a small mean shift establishes structural non-inferiority: a meaningful margin and justified uncertainty analysis are needed. Predicted RMSD does not measure binding. These are draw-weighted designs on two fixed cases, not independent biological experiments or evidence of an antibody-format effect. No p-values, confidence intervals, or causal acceptance are inferred here.'
ABLATION='Figure 1 is a separate loop-group generation cohort. Bars are mean ensemble-minus-baseline Heavy OASis Identity differences across reported seed-matched pairs; no error bars are implied. The source defines a same-decoder, same-seed zero-AbLang-logit control, but the executed helper revision/import resolution and weights are not authenticated. Thus these are descriptive retained outputs, not proof that AbLang alone caused the changes. Other CDRs means H1/H2 for 9NH7 and H1/H2 + L1/L2/L3 for 9NFU. The 9NFU heavy score is from joint heavy/light redesign, not an isolated heavy-loop intervention. No matched ablation structural results are supplied.'
LINKAGE='For the pinned original inputs, prior audit linked scored VH to final-filled FASTAs and full predicted sequences, and replayed the implemented coordinate measurements without numeric discrepancies. This package rechecks the scored-VH/FASTA link and recalculates retained metrics; it does not rerun coordinate replay, folding or BioPhi scoring. Computation verification is distinct from independent biological mapping or causal acceptance.'


def publish(root):
    import html
    import base64
    import matplotlib.pyplot as plt
    root=Path(root)
    try:
        data=load(root); table,ablation=summarize(data)
        loops,loop_readout=loop_summary(data)
    except (ValueError,KeyError,FileNotFoundError) as exc:
        message=f'Update failed: {exc}'
        (root/'READ_RESULTS.md').write_text('# Update failed\n\n'+message+'\n')
        (root/'results.html').write_text('<h1>Update failed</h1><p>'+html.escape(message)+'</p>')
        for name in ['figure1.png','figure1.svg','9NH7_distributions.png','9NFU_distributions.png']:
            (root/'figures'/name).unlink(missing_ok=True)
        for name in ['tables/table1.csv','tables/ablation.csv','tables/loops.csv','engineering/run.json','APPENDIX.md']:
            (root/name).unlink(missing_ok=True)
        raise ValueError(message) from exc
    (root/'figures').mkdir(exist_ok=True); (root/'tables').mkdir(exist_ok=True); (root/'engineering').mkdir(exist_ok=True)
    table.to_csv(root/'tables/table1.csv',index=False); ablation.to_csv(root/'tables/ablation.csv',index=False)
    answer,cases,regions=describe(data,table,ablation)
    loops.to_csv(root/'tables/loops.csv',index=False)
    light=ablation[ablation.domain=='L'][['structure','condition','n_pairs','mean_difference','positive','zero']]
    light_note='9NFU light-domain Identity: '+ '; '.join(f"{r.condition}: mean paired gain {r.mean_difference:+.7f}, {r.positive}/{r.n_pairs} positive, {r.zero}/{r.n_pairs} unchanged" for r in light.itertuples())+'. Heavy and light denominators are separate; no cross-chain efficiency inference follows.'
    fig=plot_ablation(ablation)
    bar_values=[float(p.get_height()) for ax in fig.axes for p in ax.patches]
    for ext in ['png','svg']: fig.savefig(root/f'figures/figure1.{ext}',dpi=170,metadata={'Date':None} if ext=='svg' else None)
    plt.close(fig)
    for st in CASES:
        g=data['main'].query('structure==@st')
        fig,axes=plt.subplots(2,1,figsize=(6.2,6.5),layout='constrained')
        for ax,col,label in zip(axes,['oasis','rmsd_angstrom'],['Heavy OASis Identity (fraction)','Domain-local H3 backbone RMSD (Å)']):
            for m,color in zip(METHODS,['#287c8e','#bc643f']):
                values=np.sort(g[g.method==m][col].to_numpy())
                ax.step(values,np.arange(1,len(values)+1)/len(values),where='post',label=f'{m} (n={len(values)})',color=color)
            ax.set_xlabel(label); ax.set_ylabel('Cumulative fraction'); ax.legend(); ax.grid(alpha=.15)
        fig.suptitle(st+' · all-loop cohort distributions')
        fig.savefig(root/f'figures/{st}_distributions.png',dpi=170); plt.close(fig)
    status='Pinned retained-input analysis' if data['original_inputs'] else 'Changed-input analysis — original coordinate and execution audit does not transfer to changed files.'
    link=LINKAGE if data['original_inputs'] else 'Scored VH/FASTA linkage was rechecked, but changed input bytes are not covered by the pinned coordinate audit. Interpretation remains provisional; updating a metric file is not new measurement validation.'
    display=table_display(table)
    n=int(table.humanness_n_per_arm.iloc[0]); pairs=int(ablation.n_pairs.iloc[0]); total=len(data['ablation'])
    caption=f'Each main cell is a median over {n} designs per arm per case; B = baseline, E = ensemble. Deltas are differences of marginal medians, not paired changes. Higher Identity indicates greater repertoire similarity; higher RMSD indicates greater displacement from the deposited reference.'
    abcaption=f'{pairs} reported pairs per condition; {total} ablation designs in total. '+ABLATION
    md='# Repertoire similarity and structural cost\n\n'+f'## {QUESTION}\n\n'+answer+'\n\n'+status+'\n\n## Table 1 · Two objectives, kept distinct\n\n'+markdown_table(display)+'\n\n'+caption+'\n\n'+loop_readout+'\n\n'+SCOPE+'\n\n## What differs between cases?\n\n'+'\n\n'.join(cases)+'\n\n'+LIMIT+'\n\n## Figure 1 · Where does the humanness difference occur?\n\n'+abcaption+'\n\n![Figure 1. Heavy-domain ablation differences](figures/figure1.png)\n\n'+regions+'\n\n'+light_note+'\n\n## Paper changes\n\n'+PAPER_CHANGES+'\n\n## What this evidence supports\n\n'+link+'\n\n[Supporting distributions and alternatives](APPENDIX.md) · [Sources and scientific limits](SOURCES.md) · [Update the analysis](README.md) · [Executed notebook](results.ipynb) · [Engineering checks](engineering/README.md)\n'
    (root/'READ_RESULTS.md').write_text(md)
    appendix='# Supporting analyses\n\nReturn to [Table 1 and Figure 1](READ_RESULTS.md). These plots show every retained main design; low-confidence designs are not removed. Read ECDFs as the fraction of designs at or below each x-axis value, not uncertainty intervals.\n\n'
    appendix+='## Complete selected loop table\n\n'+markdown_table(loops.round(4))+'\n\n'+LOOP_CAPTION+'\n\n[Machine-readable loop table](tables/loops.csv)\n\n'
    appendix+='## Light-domain ablation control\n\n'+markdown_table(light)+'\n\n'+light_note+'\n\n'
    for st in CASES: appendix+=f'## {st}\n\n![{st} all-loop distributions](figures/{st}_distributions.png)\n\n'
    appendix+='## Alternative endpoints are not interchangeable\n\nJacob’s `promb` human-oas peptide-content scores are not BioPhi Heavy OASis Identity. His v3 whole-framework CA fit uses different loop spans and predicted structures; in scFv it also includes inter-domain orientation. It is not an alignment-only sensitivity of Table 1. No alternate-cohort numbers are blended here.\n\nThe [chain-specific ablation table](tables/ablation.csv) keeps light-domain denominators separate. Gain per designed residue is only a heuristic: overlapping peptide windows are not additive residue effects. No structural-preservation claim follows from this ablation.\n'
    (root/'APPENDIX.md').write_text(appendix)
    def para(s): return '<p>'+html.escape(s).replace('**','')+'</p>'
    def image(name,alt):
        encoded=base64.b64encode((root/'figures'/name).read_bytes()).decode()
        return f'<img alt="{html.escape(alt)}" src="data:image/png;base64,{encoded}">'
    body='<header><p class="eyebrow">UIRJ · retained computational results</p><h1>Repertoire similarity and structural cost</h1></header><main><h2>'+html.escape(QUESTION)+'</h2>'+para(answer)+para(status)+'<h2>Table 1 · Two objectives</h2><div class="table-scroll" tabindex="0">'+display.to_html(index=False,border=0)+'</div>'+para(caption)+para(loop_readout)+para(SCOPE)+'<h2>What differs between cases?</h2>'+''.join(para(c) for c in cases)+para(LIMIT)+'<h2>Figure 1 · Where does the difference occur?</h2>'+para(abcaption)+image('figure1.png','Heavy-domain ablation differences')+para(regions)+para(light_note)+'<h2>Paper changes</h2>'+para(PAPER_CHANGES)+'<h2>What this evidence supports</h2>'+para(link)+'<details><summary>Supporting loops, distributions and alternatives</summary><div class="table-scroll">'+loops.round(4).to_html(index=False,border=0)+'</div>'+para(LOOP_CAPTION)+''.join(image(st+'_distributions.png',st+' distributions') for st in CASES)+para('ECDFs include every retained main design. Jacob promb and whole-framework CA results are different endpoints/cohorts, not an alignment-only sensitivity. See APPENDIX.md and SOURCES.md for interpretation.')+'</details><nav><a href="README.md">Update guide</a> · <a href="SOURCES.md">Sources and limits</a> · <a href="results.ipynb">Notebook</a></nav></main>'
    css='body{margin:0;background:#f7f5ef;color:#192a30;font:17px/1.65 system-ui,sans-serif}header,main{max-width:850px;margin:auto;padding:24px}header{border-bottom:3px solid #287c8e}h1{font-size:clamp(2rem,5vw,3rem);line-height:1.12}h2{font-size:1.35rem;margin-top:2em;line-height:1.3}.eyebrow{color:#287c8e;font-size:13px;text-transform:uppercase;letter-spacing:.08em}img{width:100%;height:auto;margin:20px 0}table{border-collapse:collapse;font-size:14px;min-width:650px}td,th{padding:12px;border-bottom:1px solid #bcc9c9;text-align:left}.table-scroll{overflow-x:auto}details{padding:16px;background:#e8efec}summary{cursor:pointer;font-weight:700}a{color:#116779}nav{margin:30px 0}@media(max-width:450px){header,main{padding:18px}body{font-size:16px}}'
    (root/'results.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Repertoire similarity and structural cost</title><style>'+css+'</style></head><body>'+body+'</body></html>')
    (root/'engineering/run.json').write_text(json.dumps({'input_sha256':data['hashes'],'pinned_inputs_unchanged':data['original_inputs'],'main_rows':len(data['main']),'ablation_rows':total,'bar_values_percentage_points':bar_values,'computation':'retained-output summaries; no folding, scoring or coordinate replay'},indent=2)+'\n')
    return dict(loops=loops,loop_readout=loop_readout,light_note=light_note,readout=answer,cases=cases,regions=regions,table=display,bar_values=bar_values,status=status,caption=caption,ablation_caption=abcaption)
