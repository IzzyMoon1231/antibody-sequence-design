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
    """Compatibility entry point; notebook exposes reading and validation separately."""
    from workflow import read_inputs, validate_inputs
    return validate_inputs(read_inputs(root))


def summarize(data):
    tables=[]
    for st in CASES:
        g=data['main'].query('structure==@st')
        b=g[g.method=='baseline']; e=g[g.method=='ensemble']
        tables.append(dict(structure=st,humanness_n_per_arm=len(b) if len(b)==len(e) else np.nan,h3_n_per_arm=len(b) if len(b)==len(e) else np.nan,baseline_n=len(b),ensemble_n=len(e),
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
    return pd.DataFrame(tables),pd.DataFrame(abl, columns=['structure','condition','domain','n_pairs','designed_domain_residues','peptide_denominator','baseline_mean','ensemble_mean','mean_difference','mean_gain_percentage_points','gain_per_designed_residue','mean_additional_qualifying_9mers','positive','negative','zero'])


def direction(value):
    return 'higher' if value > 1e-12 else 'lower' if value < -1e-12 else 'unchanged'


def describe(data, table, ablation):
    directions={direction(v) for v in table.biophi_heavy_difference_medians}
    opening=(f'Ensemble sequences score {next(iter(directions))} on human-antibody sequence similarity in both examples (whole heavy-domain medians).' if len(directions)==1 else 'The sequence-similarity comparison differs between the two examples.')
    answers=[opening]; cases=[]; regions=[]
    for t in table.itertuples():
        answers.append(f'{t.structure}: {direction(t.biophi_heavy_difference_medians)} heavy-domain repertoire similarity and {direction(t.h3_difference_medians_A)} H3 backbone RMSD in the ensemble median.')
        g=data['main'].query('structure==@t.structure')
        b=g[g.method=='baseline']; e=g[g.method=='ensemble']
        bq=b.rmsd_angstrom.quantile([.25,.75]); eq=e.rmsd_angstrom.quantile([.25,.75])
        overlap=max(bq.iloc[0],eq.iloc[0]) <= min(bq.iloc[1],eq.iloc[1])
        cases.append(f"**{t.structure}.** H3 interquartile ranges are {bq.iloc[0]:.3f}–{bq.iloc[1]:.3f} Å (baseline) and {eq.iloc[0]:.3f}–{eq.iloc[1]:.3f} Å (ensemble); these ranges {'overlap' if overlap else 'do not overlap'}. The ensemble middle-half spread is {direction((eq.iloc[1]-eq.iloc[0])-(bq.iloc[1]-bq.iloc[0]))}. Mean H3 RMSD is {b.rmsd_angstrom.mean():.5f} → {e.rmsd_angstrom.mean():.5f} Å; framework medians are {b.framework_rmsd_angstrom.median():.3f} → {e.framework_rmsd_angstrom.median():.3f} Å. H3 mean pLDDT is below 70 in {int(b.low_confidence.sum())}/{len(b)} → {int(e.low_confidence.sum())}/{len(e)} designs. All remain in the primary summaries; low confidence weakens interpretation as a true fold, not the retained arithmetic.")
        if ablation.empty:
            continue
        h=ablation[(ablation.structure==t.structure)&(ablation.domain=='H')].set_index('condition')
        comparison=direction(h.loc['others','mean_difference']-h.loc['cdrh3','mean_difference'])
        regions.append(f"{t.structure}: the mean heavy-domain Identity difference (ensemble minus baseline; pp = percentage points) outside H3 is {comparison} than H3-only ({100*h.loc['others','mean_difference']:+.2f} pp versus {100*h.loc['cdrh3','mean_difference']:+.2f} pp); positive reported pairs: {int(h.loc['others','positive'])}/{int(h.loc['others','n_pairs'])} versus {int(h.loc['cdrh3','positive'])}/{int(h.loc['cdrh3','n_pairs'])}.")
    if (table.h3_difference_medians_A > 1e-12).any():
        answers.append('The retained median pattern does not support a blanket no-structural-cost claim; the structural shift must be read case by case, alongside spread and confidence.')
        if (table.h3_difference_medians_A > 1e-12).all():
            largest=table.loc[table.h3_difference_medians_A.idxmax()]
            smallest=table.loc[table.h3_difference_medians_A.idxmin()]
            if largest.structure != smallest.structure:
                answers.append(f'The H3 median increase is larger for {largest.structure} ({largest.h3_difference_medians_A:+.4f} Å) than {smallest.structure} ({smallest.h3_difference_medians_A:+.4f} Å); this is an observed scale difference, not a biological-significance threshold.')
    else:
        answers.append('H3 medians do not increase in these inputs; this alone does not establish structural preservation or non-inferiority.')
    return ' '.join(answers),cases,'\n\n'.join(regions)


def loop_summary(data):
    rows=[]
    for st in CASES:
        for cdr in ('H1','H2','H3') + (('L1','L2','L3') if st=='9NFU' else ()):
            g=data['loops'].query('structure==@st and cdr==@cdr')
            row=dict(structure=st,cdr=cdr)
            for method in METHODS:
                arm=g[g.method==method]
                expected = set(data['main'].query('structure==@st and method==@method')['index'])
                observed = set(arm.design_id.str.extract(r'(\d+)$')[0].astype(int))
                if observed != expected or len(arm)!=len(expected):
                    raise ValueError('Complete loops require the same unique design indices as the scored arm')
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
PAPER_CHANGES='Use Table 1 and the complete supporting loop table to replace stale humanness and loop values; retain the heterogeneous loop-level structural costs, not universal preservation, foldability or therapeutic claims. Describe deposited-reference, motif-defined N/CA/C domain fits and whole heavy-domain Identity consistently. Keep Figure 2 descriptive, with the joint heavy/light other-CDR condition and unchanged-light control; remove per-residue or cross-chain efficiency claims. Remove index-paired Wilcoxon and inherited significance claims, completed contact-preservation claims, and structural-ablation comparisons without matching outputs. H3 is a focused endpoint, not an authenticated prospective designation. Two fixed cases are not a format ablation. These are update instructions, not certification that a manuscript has been edited.'


COLORS = ('#287c8e', '#bc643f')


def panel_limits(data):
    """Same scales across cases; include the full selected structural range."""
    return {'similarity': (-.025, 1.025),
            'structure': (-.02 * max(1., data['loops'].rmsd_angstrom.max()),
                          1.06 * max(1., data['loops'].rmsd_angstrom.max()))}


def draw_panel(ax, data, case, metric):
    """Draw from loaded frames, never a disk image. Shared by both layouts."""
    if case not in CASES or metric not in ('similarity', 'structure'):
        raise ValueError('Choose a supported case and similarity or structure')
    if metric == 'similarity':
        for j,(m,color) in enumerate(zip(METHODS,COLORS)):
            v=data['main'].query('structure==@case and method==@m').oasis.to_numpy()
            ax.scatter(v,j+np.linspace(-.16,.16,len(v)),s=12,alpha=.20,color=color)
            q=np.quantile(v,[.25,.5,.75])
            ax.hlines(j,q[0],q[2],color=color,lw=5)
            ax.plot(q[1],j,'|',color='black',ms=15,mew=2)
        ax.set(yticks=[0,1],yticklabels=METHODS,ylim=(1.5,-.5),
               xlabel='Heavy-domain BioPhi Identity (fraction)')
        title='sequence similarity'
    else:
        loops=['H1','H2','H3']+(['L1','L2','L3'] if case=='9NFU' else [])
        for j,(m,color) in enumerate(zip(METHODS,COLORS)):
            q=np.array([data['loops'].query('structure==@case and method==@m and cdr==@c').rmsd_angstrom.quantile([.25,.5,.75]).to_numpy() for c in loops])
            y=np.arange(len(loops))+(j-.5)*.25
            ax.hlines(y,q[:,0],q[:,2],color=color,lw=3)
            ax.scatter(q[:,1],y,color=color,s=35,label=m,zorder=3)
        ax.set(yticks=range(len(loops)),yticklabels=loops,ylim=(len(loops)-.5,-.5),
               xlabel='Domain-local backbone RMSD (Å)')
        ax.legend(loc='upper right',fontsize=10)
        title='every measured loop'
    ax.set_xlim(panel_limits(data)[metric])
    ax.set_title(case+' · '+title,loc='left',weight='bold',fontsize=13)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='x',alpha=.15)
    return ax


def plot_panel(data, case, metric):
    """One notebook output, explicit case and endpoint; no file writes."""
    import matplotlib.pyplot as plt
    with plt.rc_context({'font.size':12,'svg.hashsalt':'uirj-paper-focused','svg.fonttype':'none'}):
        fig,ax=plt.subplots(figsize=(6.2,3.2 if metric=='similarity' else 4.3),layout='constrained')
        draw_panel(ax,data,case,metric)
    return fig


def plot_main(data):
    """Compact manuscript layout: cases in columns, metrics in rows."""
    import matplotlib.pyplot as plt
    with plt.rc_context({'font.size':11,'svg.hashsalt':'uirj-paper-focused','svg.fonttype':'none'}):
        fig,axes=plt.subplots(2,2,figsize=(12,7.2),layout='constrained',gridspec_kw={'height_ratios':[1,1.5]})
        for row,metric in enumerate(('similarity','structure')):
            for col,case in enumerate(CASES):
                ax=draw_panel(axes[row,col],data,case,metric)
                ax.set_title(chr(65+2*row+col)+'  '+ax.get_title(loc='left'),loc='left',weight='bold',fontsize=13)
    return fig


def panel_readout(data, case, metric):
    """Case-local descriptive caption, recalculated for every selected input."""
    if case not in CASES or metric not in ('similarity','structure'):
        raise ValueError('Choose a supported case and similarity or structure')
    g=data['main'].query('structure==@case')
    counts=g.groupby('method').size()
    prefix=f"{case}: baseline n={counts['baseline']}, ensemble n={counts['ensemble']}. "
    if metric=='similarity':
        med=g.groupby('method').oasis.median()
        return prefix+f"Heavy-domain Identity medians: {med['baseline']:.4f} → {med['ensemble']:.4f}; ensemble is {direction(med['ensemble']-med['baseline'])} (Δ {med['ensemble']-med['baseline']:+.4f}). Light dots are all draws; black ticks are medians; thick lines are the middle 50%, not confidence intervals. This is repertoire similarity, not immune response. "+('The scFv score here is heavy-only, not a pooled VH/VL endpoint. ' if case=='9NFU' else '')+'Reflect: does the structural panel support the same favorable direction?'
    loops, all_loop_readout=loop_summary(data)
    z=loops[loops.structure==case]
    details='; '.join(f'{r.cdr} {r.baseline_median_A:.3f} → {r.ensemble_median_A:.3f} (Δ {r.ensemble_median_A-r.baseline_median_A:+.3f})' for r in z.itertuples())
    h=z[z.cdr=='H3'].iloc[0]
    result=prefix+'Loop median RMSD in Å, baseline → ensemble: '+details+'. Structural dots are medians; lines are the middle 50%, not confidence intervals. '
    if case=='9NFU':
        result+=all_loop_readout+' Reflect: which other loops change the conclusion you would draw from H3 alone?'
    else:
        delta=g.groupby('method').oasis.median()
        result+=f"The ensemble has {direction(delta['ensemble']-delta['baseline'])} similarity and {direction(h.ensemble_median_A-h.baseline_median_A)} H3 RMSD. Reflect: is this a similarity/structure trade-off in these inputs? It is not evidence of preserved binding."
    return result


def plot_ablation(ablation, observations):
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':14,'svg.hashsalt':'uirj-paper-focused','svg.fonttype':'none'})
    # Shared autoscaling unions every selected pair and the zero reference;
    # padding stays dynamic when replacement inputs change the range.
    fig,axes=plt.subplots(2,1,figsize=(6.2,6.5),layout='constrained',sharex=True)
    for ax in axes: ax.tick_params(axis='x',labelbottom=True)
    if observations.empty:
        for ax in axes: ax.axis('off')
        axes[0].text(.05,.5,'Ablation not selected — no comparison computed',fontsize=12)
        return fig
    for ax,st in zip(axes,CASES):
        ax.axvline(0,color='#555555',linewidth=.8)
        for j,cond in enumerate(('cdrh3','others')):
            w=observations.query('structure==@st and condition==@cond').pivot(index='seed',columns='method',values='Heavy OASis Identity')
            v=100*(w.ensemble-w.baseline).to_numpy()
            ax.scatter(v,j+np.linspace(-.18,.18,len(v)),s=16,alpha=.25,color='#287c8e')
            q=np.quantile(v,[.25,.75])
            ax.hlines(j,q[0],q[1],color='#287c8e',lw=4)
            ax.plot(v.mean(),j,'D',color='#172f38',ms=7)
            ax.annotate(f'mean {v.mean():+.2f} pp',(v.mean(),j),xytext=(0,22),textcoords='offset points',ha='center',fontsize=12)
        ax.set(yticks=[0,1],yticklabels=['H3 only','H1/H2' if st=='9NH7' else 'H1/H2 +\nL1/L2/L3'],ylim=(1.55,-.65),xlabel='Heavy-domain Identity difference\n(pp; ensemble minus baseline)')
        ax.set_title(st+' · separate ablation',loc='left',weight='bold')
        ax.spines[['top','right']].set_visible(False); ax.grid(axis='x',alpha=.15)
        ax.margins(x=.18)
    return fig


def plot_experiment_map(data):
    import matplotlib.pyplot as plt
    # Two independent cards: no pipeline or causal connector between cohorts.
    fig,axes=plt.subplots(2,1,figsize=(5.4,4.3))
    fig.subplots_adjust(left=.02,right=.98,bottom=.02,top=.98,hspace=.08)
    n=len(data['main']); a=len(data['ablation'])
    cards=[('Main comparison', f'{n} draws · all-loop redesign', [
        'Question: more similar, at what structural cost?',
        'Compare: baseline vs ensemble, in each case',
        'Measure: heavy Identity + every loop’s RMSD',
        'Contrast: unpaired group medians and spread']),
        ('Separate ablation', f'{a} draws · separate cohort' if a else 'not selected', [
        'Question: where do similarity differences appear?',
        'Compare: baseline vs ensemble within loop group',
        'Measure: heavy / light Identity separately',
        'Contrast: reported seed-pair differences',
        'H3 only vs other loops; no main-cohort pairing',
        'No structural ablation measurements'] if a else [
        'No ablation comparison computed.',
        'Select a compatible ablation file to include it.',
        'No structural ablation inference.'])]
    for ax,(title,count,lines),color,bg in zip(axes,cards,COLORS,('#eaf3f4','#faf0e9')):
        ax.set(xlim=(0,1),ylim=(0,1),xticks=[],yticks=[],facecolor=bg)
        for spine in ax.spines.values(): spine.set_visible(False)
        ax.text(.035,.94,title,fontsize=16,weight='bold',color=color,va='top')
        ax.text(.035,.79,count,fontsize=13,weight='bold',va='top')
        for i,line in enumerate(lines):
            ax.text(.035,.63-i*.105,line,fontsize=13,va='top')
    return fig


def table_display(table):
    return pd.DataFrame([{'Case':r.structure,'Heavy-domain Identity, B → E':f'{r.biophi_heavy_baseline_median:.4f} → {r.biophi_heavy_ensemble_median:.4f}', 'Δ Identity':f'{r.biophi_heavy_difference_medians:+.4f}', 'H3 RMSD, B → E (Å)':f'{r.h3_baseline_median_A:.4f} → {r.h3_ensemble_median_A:.4f}', 'Δ RMSD (Å)':f'{r.h3_difference_medians_A:+.4f}'} for r in table.itertuples()])


def markdown_table(df):
    return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in df.itertuples(index=False,name=None))


ORIENTATION = '''## Before you start

**Research question:** can adding an antibody language model make redesigned sequences more like human antibodies, and what happens to their predicted structures? Here **baseline** means ProteinMPNN alone; **ensemble** combines ProteinMPNN with AbLang. These labels describe protocols, not proof that AbLang alone caused a difference.

The two examples are **9NH7**, a heavy-only antibody domain, and **9NFU**, a linked heavy/light antibody fragment. They are two fixed cases, not biological replications or a representative antibody population. Many generated sequences (draws) from each case do not change that limit.

**Words used below**

- **Heavy / light (VH / VL):** the two kinds of antibody variable domain, the regions that help recognize a target.
- **VHH / scFv:** VHH is a heavy-only variable domain; an scFv joins VH and VL with a linker in one chain.
- **CDR / H3:** CDRs are variable loops involved in target recognition; H3 is the third heavy-domain loop. H1/H2 and L1/L2/L3 name the others.
- **BioPhi / OASis Identity:** the scoring tool and its human-repertoire similarity measure; here, the fraction of qualifying short sequence windows across the heavy domain. Higher is more similar, not proof of reduced immunogenicity.
- **RMSD (Å):** aligned backbone displacement from a reference, in ångströms. Lower is closer; predicted structure is not a binding test.
- **pLDDT:** prediction confidence (0–100), not experimental validation.
- **Ablation:** a separate controlled-change comparison, here redesigning H3 only or other loops. Its similarity scores do not establish structural preservation.

Read the main comparison first; the separate ablation asks where similarity differences appear, not whether the main designs bind their target.
'''

QUESTION='Are the sequences more human-like, and what happens to their predicted structures?'
SCOPE='Table 1 compares independently sampled all-loop baseline and ensemble protocol cohorts, not paired seeds. Heavy OASis Identity is the whole numbered heavy-domain fraction of qualifying overlapping 9-mers at the reported relaxed threshold: not an OASis percentile, isolated CDR score, or immunogenicity measurement. H3 RMSD uses a per-domain framework N/CA/C backbone fit to deposited experimental coordinates, not the original RFdiffusion model. Loops are operational sequence-alignment/motif definitions, not independently validated biological numbering. Separate domain fits exclude relative VH/VL orientation and cannot establish whole-scFv orientation or interface preservation.'
LIMIT='Neither overlapping ranges nor a small mean shift establishes structural non-inferiority: a meaningful margin and justified uncertainty analysis are needed. Predicted RMSD does not measure binding. These are draw-weighted designs on two fixed cases, not independent biological experiments or evidence of an antibody-format effect. No p-values, confidence intervals, or causal acceptance are inferred here.'
ABLATION='Figure 2 uses a separate loop-group generation cohort. Light dots show individual reported seed-pair differences; dark diamonds show the mean and thick lines the middle half (IQR), not confidence intervals. Positive differences favor ensemble sequence similarity. The source defines a same-decoder, same-seed zero-AbLang-logit control, but the executed helper revision/import resolution and weights are not authenticated. Thus these are descriptive retained outputs, not proof that AbLang alone caused the changes. Other CDRs means H1/H2 for 9NH7 and H1/H2 + L1/L2/L3 for 9NFU. The 9NFU heavy score is from joint heavy/light redesign, not an isolated heavy-loop intervention. No matched ablation structural results are supplied.'
LINKAGE='For the pinned original inputs, prior audit linked scored VH to final-filled FASTAs and full predicted sequences, and replayed the implemented coordinate measurements without numeric discrepancies. This package rechecks the scored-VH/FASTA link and recalculates retained metrics; it does not rerun coordinate replay, folding or BioPhi scoring. Computation verification is distinct from independent biological mapping or causal acceptance.'


def publish(root, data=None):
    import html
    import base64
    import matplotlib.pyplot as plt
    root=Path(root)
    if data is None:
        data=load(root)
    table,ablation=summarize(data)
    loops,loop_readout=loop_summary(data)
    (root/'figures').mkdir(exist_ok=True); (root/'tables').mkdir(exist_ok=True); (root/'engineering').mkdir(exist_ok=True)
    table.to_csv(root/'tables/table1.csv',index=False); ablation.to_csv(root/'tables/ablation.csv',index=False)
    answer,cases,regions=describe(data,table,ablation)
    loops.to_csv(root/'tables/loops.csv',index=False)
    light=ablation[ablation.domain=='L'][['structure','condition','n_pairs','mean_difference','positive','zero']]
    light_note='9NFU light-domain Identity: '+ '; '.join(f"{r.condition}: mean paired gain {r.mean_difference:+.7f}, {r.positive}/{r.n_pairs} positive, {r.zero}/{r.n_pairs} unchanged" for r in light.itertuples())+'. Heavy and light denominators are separate; no cross-chain efficiency inference follows.'
    fig=plot_ablation(ablation, data['ablation'])
    bar_values=[float(line.get_xdata()[0]) for ax in fig.axes for line in ax.lines[1:]]
    for ext in ['png','svg']: fig.savefig(root/f'figures/figure2.{ext}',dpi=170,metadata={'Date':None} if ext=='svg' else None)
    plt.close(fig)
    figures = [('figure1',plot_main(data)),('experiment_map',plot_experiment_map(data))]
    figures += [(f'{st}_{metric}', plot_panel(data,st,metric)) for st in CASES for metric in ('similarity','structure')]
    for name, figure in figures:
        for ext in ('png','svg'):
            figure.savefig(root/f'figures/{name}.{ext}',dpi=170,metadata={'Date':None} if ext=='svg' else None)
        plt.close(figure)
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
    n=' / '.join(f'{r.structure}: B={r.baseline_n}, E={r.ensemble_n}' for r in table.itertuples()); pairs=' / '.join(f'{r.structure}/{r.condition}/{r.domain}: {r.n_pairs}' for r in ablation.itertuples()) or 'none selected'; total=len(data['ablation'])
    caption=f'Each main cell is a marginal median; selected draw counts: {n}; B = baseline, E = ensemble. Deltas are differences of marginal medians, not paired changes. Higher Identity indicates greater repertoire similarity; higher RMSD indicates greater displacement from the deposited reference.'
    abcaption=f'Reported pair counts: {pairs}; {total} ablation designs in total. '+ABLATION
    main_caption=f'Light sequence dots retain all observed draws ({n}); black ticks mark medians. Structural dots are group medians for every measured loop. Thick lines show the middle half (25th–75th percentiles), not confidence intervals. All draws keep their original weight, including repeated sequences. Baseline and ensemble are separate collections: no row-pairing or connecting designs. Higher Identity means more human-repertoire sequence similarity; lower RMSD means closer to the deposited reference. These are whole heavy-domain scores and domain-local loop distances, not immune-response or binding measurements.'
    map_caption='Two questions, two separate collections. Main comparison: all-loop redesign, sequence similarity and predicted structure. Ablation: H3-only versus other-loop redesign, sequence similarity only. Cards show each question, measurement and comparison; neither is a stage of the other.'
    (root/'AUTHOR_CHECKLIST.md').write_text('# Author checklist — manuscript integration\n\n'+PAPER_CHANGES+'\n\n- [ ] Replace stale values using [Table 1](tables/table1.csv) and [all loops](tables/loops.csv).\n- [ ] Use [Figure 1](figures/figure1.svg) for the main question and [Figure 2](figures/figure2.svg) for the separate ablation.\n- [ ] Keep conventions in [Sources and methods](SOURCES.md); keep the experiment map in the notebook unless the paper needs it.\n- [ ] Review all claims against the qualifications above; notebook completion is not manuscript endorsement.\n')
    md='# Repertoire similarity and structural cost\n\n'+ORIENTATION+'\n\n'+f'## {QUESTION}\n\n'+answer+'\n\n'+status+'\n\n'+map_caption+'\n\n![Experiment map](figures/experiment_map.png)\n\n## Figure 1 · Sequence similarity and every measured loop\n\n'+main_caption+'\n\n'+''.join(f'### {st} · {metric}\n\n![{st} {metric}](figures/{st}_{metric}.png)\n\n'+panel_readout(data,st,metric)+'\n\n' for st in CASES for metric in ('similarity','structure'))+'[Download compact 2×2 Figure 1 (SVG)](figures/figure1.svg) · [PNG](figures/figure1.png)\n\n'+loop_readout+'\n\n## Table 1 · Numerical reference\n\n'+markdown_table(display)+'\n\n'+caption+'\n\n## Figure 2 · Where are differences seen in the separate ablation?\n\n'+abcaption+'\n\n![Figure 2. Reported seed-pair differences](figures/figure2.png)\n\n'+regions+'\n\n'+light_note+'\n\n## What this evidence supports\n\n'+LIMIT+'\n\n'+link+'\n\n[Spread and confidence details](APPENDIX.md) · [Sources and methods](SOURCES.md) · [Update guide](README.md) · [Executed notebook](results.ipynb) · [Author checklist — Paper changes](AUTHOR_CHECKLIST.md) · [Engineering checks](engineering/README.md)\n'
    (root/'READ_RESULTS.md').write_text(md)
    appendix='# Supporting analyses\n\nReturn to [Table 1 and Figure 1](READ_RESULTS.md). These plots show every retained main design; low-confidence designs are not removed. Read ECDFs as the fraction of designs at or below each x-axis value, not uncertainty intervals.\n\n'
    appendix+='## Spread and confidence by case\n\n'+'\n\n'.join(cases)+'\n\n'+SCOPE+'\n\n'
    appendix+='## Complete selected loop table\n\n'+markdown_table(loops.round(4))+'\n\n'+LOOP_CAPTION+'\n\n[Machine-readable loop table](tables/loops.csv)\n\n'
    appendix+='## Light-domain ablation control\n\n'+markdown_table(light)+'\n\n'+light_note+'\n\n'
    for st in CASES: appendix+=f'## {st}\n\n![{st} all-loop distributions](figures/{st}_distributions.png)\n\n'
    appendix+='## Alternative endpoints are not interchangeable\n\nJacob’s `promb` human-oas peptide-content scores are not BioPhi Heavy OASis Identity. His v3 whole-framework CA fit uses different loop spans and predicted structures; in scFv it also includes inter-domain orientation. It is not an alignment-only sensitivity of Table 1. No alternate-cohort numbers are blended here.\n\nThe [chain-specific ablation table](tables/ablation.csv) keeps light-domain denominators separate. Gain per designed residue is only a heuristic: overlapping peptide windows are not additive residue effects. No structural-preservation claim follows from this ablation.\n'
    (root/'APPENDIX.md').write_text(appendix)
    def para(s): return '<p>'+html.escape(s).replace('**','')+'</p>'
    def image(name,alt):
        encoded=base64.b64encode((root/'figures'/name).read_bytes()).decode()
        return f'<img alt="{html.escape(alt)}" src="data:image/png;base64,{encoded}">'
    orientation_html='<section aria-label="Before you start">'
    for block in ORIENTATION.strip().split('\n\n'):
        if block.startswith('## '):
            orientation_html+='<h2>'+html.escape(block[3:])+'</h2>'
        elif block.startswith('- '):
            orientation_html+='<ul>'+''.join('<li>'+html.escape(line[2:]).replace('**','')+'</li>' for line in block.splitlines())+'</ul>'
        else:
            orientation_html+=para(block)
    orientation_html+='</section>'
    body='<header><p class="eyebrow">UIRJ · retained computational results</p><h1>Repertoire similarity and structural cost</h1></header><main>'+orientation_html+'<h2>'+html.escape(QUESTION)+'</h2>'+para(answer)+para(status)+para(map_caption)+image('experiment_map.png','Two separate experiment collections')+'<h2>Figure 1 · Main comparison</h2>'+para(main_caption)+''.join('<h3>'+st+' · '+metric+'</h3>'+image(f'{st}_{metric}.png',f'{st}_{metric}')+para(panel_readout(data,st,metric)) for st in CASES for metric in ('similarity','structure'))+'<p><a href="figures/figure1.svg">Download compact 2×2 Figure 1 (SVG)</a> · <a href="figures/figure1.png">PNG</a></p>'+para(loop_readout)+'<h2>Table 1 · Numerical reference</h2><div class="table-scroll" tabindex="0">'+display.to_html(index=False,border=0)+'</div>'+para(caption)+'<h2>Figure 2 · Separate ablation</h2>'+para(abcaption)+image('figure2.png','Reported seed-pair differences, means and IQR')+''.join(para(region) for region in regions.split('\n\n'))+para(light_note)+'<h2>What this evidence supports</h2>'+para(LIMIT)+para(link)+'<details><summary>Spread, confidence and complete loop table</summary>'+''.join(para(c) for c in cases)+'<div class="table-scroll">'+loops.round(4).to_html(index=False,border=0)+'</div>'+para(LOOP_CAPTION)+'</details><nav><a href="README.md">Update guide</a> · <a href="SOURCES.md">Sources and methods</a> · <a href="APPENDIX.md">Supporting results</a> · <a href="results.ipynb">Notebook</a> · <a href="AUTHOR_CHECKLIST.md">Author checklist — Paper changes</a></nav></main>'
    css='body{margin:0;background:#f7f5ef;color:#192a30;font:17px/1.65 system-ui,sans-serif}header,main{max-width:850px;margin:auto;padding:24px}header{border-bottom:3px solid #287c8e}h1{font-size:clamp(2rem,5vw,3rem);line-height:1.12}h2{font-size:1.35rem;margin-top:2em;line-height:1.3}.eyebrow{color:#287c8e;font-size:13px;text-transform:uppercase;letter-spacing:.08em}img{width:100%;height:auto;margin:20px 0}img[alt="Two separate experiment collections"]{max-width:540px;display:block}table{border-collapse:collapse;font-size:14px;min-width:650px}td,th{padding:12px;border-bottom:1px solid #bcc9c9;text-align:left}.table-scroll{overflow-x:auto}details{padding:16px;background:#e8efec}summary{cursor:pointer;font-weight:700}a{color:#116779}nav{margin:30px 0}@media(max-width:450px){header,main{padding:18px}body{font-size:16px}}'
    (root/'results.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Repertoire similarity and structural cost</title><style>'+css+'</style></head><body>'+body+'</body></html>')
    (root/'engineering/run.json').write_text(json.dumps({'input_sha256':data['hashes'],'pinned_inputs_unchanged':data['original_inputs'],'main_rows':len(data['main']),'ablation_rows':total,'bar_values_percentage_points':bar_values,'computation':'retained-output summaries; no folding, scoring or coordinate replay'},indent=2)+'\n')
    return dict(main_caption=main_caption,map_caption=map_caption,loops=loops,loop_readout=loop_readout,light_note=light_note,readout=answer,cases=cases,regions=regions,table=display,bar_values=bar_values,status=status,caption=caption,ablation_caption=abcaption)
