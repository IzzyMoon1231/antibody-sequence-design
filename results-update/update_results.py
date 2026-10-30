"""Offline descriptive update helpers. No scoring, network, or upstream imports."""
import csv
import hashlib
import re
import math
from reviewed_stats import stats as reviewed_stats
from collections import defaultdict
from pathlib import Path
import json
import sys
import shutil
import uuid
from datetime import datetime, timezone


def load_run(root, config):
    root=Path(root)
    run={'config':config,'sequences':[],'metrics':[],'missing':list(config.get('missing_artifacts',[])),'input_hashes':{},'input_paths':{}}
    for category in ('fasta','designs','metrics'):
        for spec in config.get(category,[]):
            path=root/spec['path']
            if not path.is_file():
                run['missing'].append(f'MISSING FILE: {spec["path"]}; its results were not calculated.')
                continue
            run['input_hashes'][spec['path']]=hashlib.sha256(path.read_bytes()).hexdigest()
            run['input_paths'][spec['path']]=str(path.resolve())
            try:
                if category=='fasta':
                    run['sequences'].extend(read_fasta(path,spec['structure'],spec['method'],spec['cohort']))
                elif category=='designs':
                    run['sequences'].extend(read_designs(path,spec['cohort']))
                else:
                    run['metrics'].extend(read_metrics(root,spec))
            except (ValueError,KeyError) as exc:
                raise ValueError(f'{spec["path"]}: {exc}') from exc
    # Duplicate sequence identities are forbidden within each cohort.
    cohorts=sorted({r['cohort'] for r in run['sequences']})
    for c in cohorts: _index([r for r in run['sequences'] if r['cohort']==c])
    run['sequence_comparison']=[]
    for i,c in enumerate(cohorts):
        for d in cohorts[i+1:]:
            left=[r for r in run['sequences'] if r['cohort']==c]
            right=[r for r in run['sequences'] if r['cohort']==d]
            run['sequence_comparison'].extend(compare_sequences(left,right))
    run['coverage']=check_coverage(run['sequences'],run['metrics'])
    run['summaries']=summarize(run['metrics'])
    run['source_records']=[]
    for manifest_name in ('acquisition.json','xlsx-acquisition.json'):
        manifest=root/manifest_name
        if manifest.is_file():
            records=json.loads(manifest.read_text())
            for record in records if isinstance(records,list) else [records]:
                if run['input_hashes'].get(record['path'])==record['sha256']:
                    run['source_records'].append(record)
    run['code_hashes']={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() for name in ('update_results.py','reviewed_stats.py','reviewed-helper-provenance.json','Update_results.ipynb') if (Path(__file__).parent/name).is_file()}
    run['changed_inputs']={'status':'NO_PREVIOUS_BUNDLE_SELECTED','changed':[],'added':[],'removed':[]}
    if config.get('previous_bundle'):
        previous=json.loads((root/config['previous_bundle']/'provenance.json').read_text())['input_hashes']
        current=run['input_hashes']
        run['changed_inputs']=dict(status='COMPARED_INPUT_BYTES_ONLY',changed=sorted(k for k in previous.keys() & current.keys() if previous[k]!=current[k]),added=sorted(current.keys()-previous.keys()),removed=sorted(previous.keys()-current.keys()))
    return run


def export_run(run, output_root):
    """Append-only snapshot; never edit scientific inputs or historical bundles."""
    out=Path(output_root)/('update-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')+'-'+uuid.uuid4().hex[:8])
    out.mkdir(parents=True,exist_ok=False)
    for name,key in [('summary','summaries'),('measurements','metrics'),('sequence_comparison','sequence_comparison'),('coverage','coverage'),('sequences','sequences')]:
        rows=run[key]
        with (out/(name+'.csv')).open('x',newline='') as f:
            if rows:
                w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
    import importlib.metadata
    provenance=dict(scientific_status='SOURCE_EXPORT_REVIEW_ONLY',created_utc=datetime.now(timezone.utc).isoformat(),
                    source_records=run['source_records'],code_hashes=run['code_hashes'],
                    dependency_versions={name:importlib.metadata.version(name) for name in ('matplotlib','pandas','openpyxl')},
                    python=sys.version,helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    input_hashes=run['input_hashes'],config=run['config'],changed_inputs=run['changed_inputs'],missing=run['missing'])
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2,allow_nan=False)+'\n')
    copied={}
    for i,(label,path) in enumerate(run['input_paths'].items()):
        target=out/'inputs'/(f'{i:02d}-'+Path(path).name);target.parent.mkdir(exist_ok=True)
        data=Path(path).read_bytes()
        if hashlib.sha256(data).hexdigest()!=run['input_hashes'][label]:
            raise ValueError(f'Input changed after preview: {label}; restart and run all')
        target.write_bytes(data);copied[label]=str(target.relative_to(out))
    (out/'input-copy-map.json').write_text(json.dumps(copied,indent=2)+'\n')
    notes=['# Private results-update bundle','',run['config'].get('title',''),'',
           '**These are summaries of saved exports, not newly scored or structurally validated results.**',
           'No paired tests, correlations, confidence intervals, non-inferiority conclusions, or causal claims are calculated.',
           'Matching identifiers alone is not proof that a score belongs to a sequence. No index offsets were applied.',
           'OASis is a repertoire-similarity proxy, not clinical immunogenicity; RMSD is not binding affinity.',
           'Sources, cases, cohorts, domains and endpoints stay separate. Duplicate-sequence draws remain draws, not independent biological replicates.',
           '', '## What still needs attention']+['- '+x for x in run['missing']]
    mismatch=sum(r['status']=='SEQUENCE_MISMATCH' for r in run['sequence_comparison'])
    offsets=sum(r['possible_index_offset'] for r in run['coverage'])
    notes+=['',f'Same-ID sequence disagreements: {mismatch}. Coverage groups with possible ±1 ID offsets: {offsets}. See sequence_comparison.csv and coverage.csv for exact records.',
            '', '## What changed',json.dumps(run['changed_inputs'],indent=2),'', '## Figures','Saved-export distributions only; every panel is separate. Do not combine panels into an effect estimate.']
    figures=plot_exports(run['metrics'],out/'figures')
    notes += [f'- [{p.name}](figures/{p.name})' for p in figures]
    (out/'READ_ME_FIRST.md').write_text('\n'.join(notes)+'\n')
    hashes={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file()}
    (out/'SHA256SUMS.json').write_text(json.dumps(hashes,indent=2)+'\n')
    return out


def plot_exports(rows, folder):
    from result_figures import build_figures
    import matplotlib.pyplot as plt
    folder=Path(folder);folder.mkdir(exist_ok=True)
    paths=[]
    for i,(key,fig) in enumerate(build_figures(rows)):
        png=folder/f'{i:02d}-{key[2]}.png'
        fig.savefig(png,dpi=140);fig.savefig(png.with_suffix('.svg'))
        plt.close(fig);paths.append(png)
    return paths


def read_metrics(root, spec):
    """Explicit columns; raw scores are NOT authenticated sequence measurements."""
    path=Path(root)/spec['path']
    if path.suffix.lower()=='.xlsx':
        import pandas as pd
        if spec.get('sheet')!='Overview':
            raise ValueError('Use the BioPhi Overview sheet, not per-peptide rows')
        source=pd.read_excel(path,sheet_name='Overview',dtype=str,keep_default_na=False).to_dict('records')
    else:
        source=list(csv.DictReader(path.read_text().splitlines()))
    source=[dict(r,structure=r.get('structure',spec.get('structure')),method=r.get('method',spec.get('method'))) for r in source]
    out=[]
    aliases={'ProteinMPNN':'baseline','proteinmpnn':'baseline','baseline':'baseline','ensemble':'ensemble','AbMPNN':'AbMPNN','abmpnn':'AbMPNN'}
    for r in source:
        structure=r['structure']; method=aliases[r['method']]
        raw_id=r[spec['id_column']]
        patterns=[r'(\d+)',rf'{structure}_{r["method"]}_(\d+)']
        if method=='baseline': patterns.append(rf'{structure}_ProteinMPNN_sample(\d+)')
        if method=='ensemble': patterns.extend([rf'{structure}_seed(\d+)',rf'{structure}_ensemble_seed(\d+)'])
        matches=[m for pat in patterns if (m:=re.fullmatch(pat,raw_id))]
        if len(matches)!=1:
            raise ValueError(f'Unrecognized ID {raw_id}; specify real ID convention, never shift silently')
        value=r[spec['value_column']]
        value=None if value.strip() in ('','NA','NaN','nan') else float(value)
        if value is not None and (not math.isfinite(value) or value<0 or (spec['unit']=='fraction' and value>1)):
            raise ValueError(f'Invalid {spec["unit"]} value for {raw_id}: {value}')
        out.append(dict(structure=structure,method=method,cohort=spec['cohort'],record_id=str(int(matches[0][1])),
                        source_id=raw_id,source=spec['path'],domain=r[spec['domain_column']] if spec.get('domain_column') else 'not-specified-in-export',
                        endpoint=spec['endpoint'],unit=spec['unit'],value=value,status='UNBOUND_SOURCE_EXPORT',source_settings='Threshold='+r.get('Threshold','not retained'),source_row=dict(r)))
    seen=set()
    for r in out:
        key=tuple(r[k] for k in ('structure','method','record_id','domain','endpoint'))
        if key in seen: raise ValueError(f'Duplicate measurement: {key}')
        seen.add(key)
    if not out: raise ValueError(f'Empty metric table: {path}')
    return out


def check_coverage(sequences, metrics):
    """Compare literal IDs within cohort; a possible offset is a warning, not a map."""
    groups=defaultdict(list)
    for r in metrics:
        groups[tuple(r[k] for k in ('source','cohort','structure','method','domain','endpoint'))].append(r)
    out=[]
    for key,rows in sorted(groups.items()):
        source,cohort,structure,method,domain,endpoint=key
        expected={r['record_id'] for r in sequences if (r['cohort'],r['structure'],r['method'])==(cohort,structure,method)}
        observed={r['record_id'] for r in rows}
        offset=bool(expected and observed and expected!=observed and any({int(x)+d for x in observed}=={int(x) for x in expected} for d in (-1,1)))
        out.append(dict(source=source,cohort=cohort,structure=structure,method=method,domain=domain,endpoint=endpoint,
                        n_sequences=len(expected),n_metric_ids=len(observed),n_missing_values=sum(r['value'] is None for r in rows),
                        sequence_ids_without_metric=sorted(expected-observed),metric_ids_without_sequence=sorted(observed-expected),
                        possible_index_offset=offset,status='ID_CHECK_ONLY_NOT_SEQUENCE_PROOF'))
    return out


def bind_metrics(metrics, sequences):
    """Only explicit measured-input hashes permit sequence attachment; not scientific approval."""
    index={(r['cohort'],r['structure'],r['method'],r['record_id']):r for r in sequences}
    out=[]
    for r in metrics:
        key=tuple(r[k] for k in ('cohort','structure','method','record_id'))
        seq=index.get(key)
        if seq is None or r.get('sequence_sha256')!=seq['sequence_sha256']:
            raise ValueError(f'Missing or mismatched measured-input sequence hash: {key}. Obtain the actual scored/folded input; do not copy a hash from the desired FASTA.')
        out.append(dict(r,status='SEQUENCE_HASH_MATCH_ONLY'))
    return out


def summarize(rows):
    groups=defaultdict(list)
    keys=('source','cohort','structure','method','domain','endpoint','unit','status')
    for r in rows: groups[tuple(r[k] for k in keys)].append(r['value'])
    result=[]
    for key,values in sorted(groups.items()):
        present=[v for v in values if v is not None]
        values_summary=reviewed_stats(present) if present else dict(n=0,mean=None,median=None,min=None,max=None)
        result.append(dict(zip(keys,key),n_rows=len(values),n_present=values_summary.pop('n'),n_missing=len(values)-len(present),**values_summary))
    return result



def read_fasta(path, structure, method, cohort):
    """Preserve sample/seed numbers, never shift or pair them."""
    records=[]
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith('>'):
            records.append([line[1:], ''])
        elif records:
            records[-1][1] += line.strip()
        else:
            raise ValueError('Sequence before FASTA header')
    out=[]
    for header, seq in records:
        pattern={'baseline':r'\bsample=(\d+)\b','ensemble':r'\bseed=(\d+)\b','AbMPNN':rf'^{structure}_AbMPNN_(\d+)$'}[method]
        match=re.search(pattern,header)
        if not match:
            raise ValueError(f'Unrecognized {method} header: {header}; do not silently skip native records')
        if not seq or set(seq)-set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError(f'Empty or noncanonical sequence: {header}; supply resolved inputs, do not fill automatically')
        out.append(dict(structure=structure,method=method,cohort=cohort,record_id=str(int(match[1])),header=header,sequence=seq,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()))
    ids=[r['record_id'] for r in out]
    if not out or len(ids)!=len(set(ids)):
        raise ValueError('Empty FASTA or duplicate sample/seed IDs')
    return out


def read_designs(path, cohort):
    rows=[]
    for r in csv.DictReader(Path(path).read_text().splitlines()):
        method={'proteinmpnn':'baseline','ensemble':'ensemble'}[r['method']]
        seq=r['full_seq_filled']
        if not seq or set(seq)-set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('Noncanonical full_seq_filled')
        rows.append(dict(structure=r['structure'],method=method,record_id=r['sample'],cohort=cohort,sequence=seq,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()))
    _index(rows)
    return rows


def _index(rows):
    result={}
    for row in rows:
        key=(row['structure'],row['method'],row['record_id'])
        if key in result:
            raise ValueError(f'Duplicate sequence identity: {key}')
        result[key]=row
    return result


def compare_sequences(left, right):
    """Diagnostic only: no shifted joins, no cohort reassignment."""
    a,b=_index(left),_index(right)
    out=[]
    for key in sorted(a.keys() | b.keys()):
        x,y=a.get(key),b.get(key)
        status=('MISSING_LEFT' if x is None else 'MISSING_RIGHT' if y is None else
                'EXACT_SEQUENCE_MATCH' if x['sequence']==y['sequence'] else 'SEQUENCE_MISMATCH')
        out.append(dict(structure=key[0],method=key[1],record_id=key[2],status=status,
                        left_sha256=x['sequence_sha256'] if x else '',right_sha256=y['sequence_sha256'] if y else '',
                        left_cohort=x['cohort'] if x else '',right_cohort=y['cohort'] if y else ''))
    return out
