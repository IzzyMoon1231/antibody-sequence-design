"""Saved-result figures: pure row selection; no cached globals, joins or scoring."""
from pathlib import Path
from collections import defaultdict
import io
import textwrap

COLORS={'baseline':'#2379AB','ensemble':'#D66B22','AbMPNN':'#7A5195'}


def build_figures(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    groups=defaultdict(list)
    keys=('source','cohort','structure','domain','endpoint','unit')
    for r in rows: groups[tuple(r[k] for k in keys)].append(r)
    figures=[]
    for key, group in sorted(groups.items()):
        source,cohort,case,domain,endpoint,unit=key
        fig,ax=plt.subplots(figsize=(6.4,5.0),layout=None)
        for method in sorted({r['method'] for r in group}):
            selected=[r for r in group if r['method']==method]
            values=sorted(r['value'] for r in selected if r['value'] is not None)
            label=f'{method} · n={len(values)}, missing={len(selected)-len(values)}'
            if values:
                ax.step(values,[(i+1)/len(values) for i in range(len(values))],where='post',
                        color=COLORS.get(method),linewidth=2.2,label=label)
            else: ax.plot([],[],label=label,color=COLORS.get(method))
        maximum=max((r['value'] for r in rows if r['endpoint']==endpoint and r['unit']==unit and r['value'] is not None),default=1)
        ax.set_xlim(0,1 if unit=='fraction' else max(1,maximum*1.05))
        ax.set_ylim(0,1.02)
        ax.set_xlabel('Saved score (fraction, 0–1)' if unit=='fraction' else 'Reported RMSD (Å)',fontsize=11)
        ax.set_ylabel('Fraction of saved values ≤ x',fontsize=11)
        ax.grid(alpha=.18)
        ax.spines[['top','right']].set_visible(False)
        ax.legend(loc='lower right',fontsize=9,framealpha=.95)
        shown_domain='domain not specified' if domain=='not-specified-in-export' else domain
        fig.suptitle(f'{case} · {shown_domain}\n'+textwrap.fill(endpoint,52),x=.12,y=.98,ha='left',fontsize=12,fontweight='bold')
        fig.text(.12,.12,cohort,fontsize=9)
        fig.text(.12,.08,Path(source).name,fontsize=8)
        fig.text(.12,.035,'Saved export · sequence/run linkage unverified',fontsize=9,color='#735125')
        fig.subplots_adjust(left=.12,right=.97,bottom=.25,top=.77)
        plt.close(fig)  # Keep the figure object without accumulating GUI managers.
        figures.append((key,fig))
    return figures


def render_sections(rows, unit):
    from IPython.display import display, Image, Markdown
    import matplotlib.pyplot as plt
    figures=build_figures([r for r in rows if r['unit']==unit])
    for key,fig in figures:
        source,cohort,case,domain,endpoint,_=key
        display(Markdown(f'### {case} / {domain} — {endpoint}\n\nCohort: **{cohort}**. Source: `{source}`. Unpaired saved-value distribution.'))
        buffer=io.BytesIO()
        fig.savefig(buffer,format='png',dpi=140)
        display(Image(data=buffer.getvalue(),width=640))
        plt.close(fig)
    if not figures: display(Markdown('No selected measurements in this section.'))
