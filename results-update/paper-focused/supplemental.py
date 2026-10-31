"""Opt-in exploratory appendix; no writes during calculation; no upstream runs."""
from pathlib import Path
from io import BytesIO
import hashlib
import json
import platform
import numpy as np
import pandas as pd
import scipy
import stats_helper as stats
import workflow

HERE = Path(__file__).resolve().parent
RECEIPT_SHA256 = '5f1caccd6ca48e4202c338d70a67f5cf34f4d2e4658b072ab9a8f390e0787c25'
ASSUMPTIONS = ('Exploratory conditional computational-draw inference, not biological replication. '
 'Independent draws within/across arms; exchangeable labels under the group null and exchangeable '
 'ranks under the association null. Duplicates retained with draw weight do not establish independence. '
 'Known shared trajectories/copies/clusters require a different design: do not interpret iid p or intervals. '
 '95% percentile independently resampled-arm HL intervals are nonsimultaneous. '
 'MW is not generally a median/HL null; location-shift interpretation requires comparable shapes. '
 'Positive H3 RMSD shift is more deviation, not improvement. No isolated AbLang causation.')
MISSING = ('Withheld: current scored VH -> full FASTA -> actual predicted-coordinate sequence -> '
 'content-bound H3 measurement lineage is missing. ID compatibility alone is insufficient. '
 'This package accepts only its byte-bound original join receipt; changed files need a new independently '
 'reviewed coordinate/measurement receipt and adapter, not an uploaded assertion or copied sequence column.')


def export(root, data, result, *, confirm=False, review_plan=''):
    """Separate opt-in dated bundle; recompute and compare before any write."""
    if confirm is not True or not review_plan.strip() or result is None:
        raise ValueError('Explicit export opt-in, current results and a review plan are required')
    settings=result['receipt']['settings']
    current=run(data,enabled=True,n_perm=settings['n_perm_requested'],n_boot=settings['n_boot'],seed=settings['seed'])
    if current['receipt']!=result['receipt'] or any(not current[k].equals(result[k]) for k in ('comparisons','associations','families')):
        raise ValueError('Statistics or inputs changed; restart and recompute before export')
    from datetime import datetime, timezone
    import uuid
    target=Path(root)/'exports'/('statistics-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')+'-'+uuid.uuid4().hex[:8])
    target.mkdir(parents=True,exist_ok=False)
    for key in ('comparisons','associations','families'):
        current[key].to_csv(target/(key+'.csv'),index=False)
    record=dict(current['receipt'],review_plan=review_plan,adopted_policy=False,selection=data['loaded']['paths'])
    record['code_sha256']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('stats_helper.py','supplemental.py','workflow.py','analysis.py')}
    (target/'receipt.json').write_text(json.dumps(record,indent=2)+'\n')
    for name,blob in data['loaded']['blobs'].items():
        p=workflow.local_path(target/'selected-inputs',name)
        p.parent.mkdir(parents=True,exist_ok=True); p.write_bytes(blob)
    for name in ('stats_helper.py','supplemental.py','workflow.py','analysis.py','STATISTICS.md'):
        (target/name).write_bytes((HERE/name).read_bytes())
    if record['joined_rows']:
        (target/'inputs').mkdir()
        for name in ('association_receipt.json','audited_join_rows.csv'):
            (target/'inputs'/name).write_bytes((HERE/'inputs'/name).read_bytes())
    (target/'REPLAY.md').write_text('Load exact files from selected-inputs with workflow.read_inputs using the selected paths in your saved notebook; validate, then supplemental.run(enabled=True, **your recorded settings). The receipt records the requested permutation budget, bootstrap count and seed. Keep both export switches off. Save your executed notebook alongside this bundle. Original association permission is byte-bound, not inferred from filenames.\n')
    return target


def verified_rows(data):
    """Reuse a historical audit only for precisely its input bytes and values."""
    checked = workflow.validate_inputs(data['loaded'])
    if not checked['main'].equals(data['main']):
        raise ValueError('Current main frame changed after validation; restart and read inputs again')
    path = HERE/'inputs/association_receipt.json'
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=RECEIPT_SHA256:
        return None, MISSING
    receipt=json.loads(path.read_bytes())
    blobs=data['loaded']['blobs']
    needed={k:v for k,v in receipt['input_sha256'].items() if not k.endswith('ablation_scores.csv')}
    if any(k not in blobs or hashlib.sha256(blobs[k]).hexdigest()!=v for k,v in needed.items()):
        return None, MISSING
    p=HERE/'inputs/audited_join_rows.csv'
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=receipt['join_sha256']:
        return None, MISSING
    rows=pd.read_csv(BytesIO(p.read_bytes()))
    rows['index']=rows.id.str.extract(r'_(\d+)$')[0].astype(int)
    keys=['structure','method','index']
    a=rows.sort_values(keys).reset_index(drop=True)
    b=checked['main'].sort_values(keys).reset_index(drop=True)
    if not a[keys].equals(b[keys]) or not np.allclose(a.score,b.oasis,rtol=0,atol=1e-15) or not np.array_equal(a.rmsd,b.rmsd_angstrom):
        return None, MISSING
    return rows, receipt['qualification']


def run(data, *, enabled=False, n_perm=100000, n_boot=10000, seed=42):
    if enabled is not True:
        return None
    joined, qualification=verified_rows(data)
    comparisons=[]; associations=[]
    for case,g in data['main'].groupby('structure'):
        for endpoint in ['oasis','rmsd_angstrom']:
            comparisons.append(dict(test_id=f'{case}:{endpoint}:MW',case=case,endpoint=endpoint,
                **stats.compare(g.loc[g.method=='baseline',endpoint],g.loc[g.method=='ensemble',endpoint],n_perm=n_perm,n_boot=n_boot,seed=seed)))
        for method in ['baseline','ensemble']:
            result=dict(test_id=f'{case}:{method}:rho',case=case,method=method)
            if joined is None:
                result.update(rho=None,p=None,status=qualification)
            else:
                a=joined[(joined.structure==case)&(joined.method==method)]
                result.update(stats.associate(a[['id','sequence','score']],a[['id','fold_sequence','rmsd']].rename(columns={'fold_sequence':'sequence'}),n_perm=n_perm,seed=seed))
            associations.append(result)
    allp={r['test_id']:r['p'] for r in comparisons+associations}
    six=[r['test_id'] for r in comparisons if r['endpoint']=='oasis']+[r['test_id'] for r in associations]
    families=[]
    for name,members,method in [('draft_six',six,'BH'),('expanded_eight',list(allp),'BH'),('expanded_eight',list(allp),'Holm')]:
        complete=all(allp[k] is not None and np.isfinite(allp[k]) for k in members)
        adjusted=stats.adjust({k:allp[k] for k in members},family=members,method=method) if complete else {}
        families.extend(dict(family=name,family_size=len(members),test_id=k,p=allp[k],candidate_adjusted_p=adjusted.get(k),adjustment=method,status='candidate only; not adopted' if complete else 'withheld: incomplete family') for k in members)
    return dict(comparisons=pd.DataFrame(comparisons),associations=pd.DataFrame(associations),families=pd.DataFrame(families),receipt=dict(
        input_sha256={k:hashlib.sha256(v).hexdigest() for k,v in data['loaded']['blobs'].items()},
        source_receipt_sha256=RECEIPT_SHA256 if joined is not None else None,
        joined_rows=0 if joined is None else len(joined),qualification=qualification,
        settings=dict(n_perm_requested=n_perm,n_boot=n_boot,seed=seed),assumptions=ASSUMPTIONS,
        method='Two-sided absolute rank/label permutation; ties average ranked; plus-one Monte Carlo; exact branches recorded per row; independent-arm percentile HL bootstrap',
        family_policy='Draft six: two Identity MW + four within-case/arm rho. Expanded eight adds two H3 MW. BH assumes independent/appropriate positively dependent tests; Holm controls FWER without dependence requirement. Neither adopted; never choose by favorable p.',
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)))
