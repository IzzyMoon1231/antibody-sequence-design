"""Local, read-only intake and explicit new-folder export for results.ipynb.
No inference, scoring, network or automatic repairs. See README for the contract.
"""
from pathlib import Path
from io import BytesIO
from datetime import datetime, timezone
import hashlib
import json
import re
import uuid
import numpy as np
import pandas as pd
import analysis

# Trust anchor for the shipped public provenance list, not an upload's assertion.
PROVENANCE_SHA256 = '530fcef239e0f7c99f8ae6a5881114bbfc298720ac1a1f436c5a23e62704601f'


def local_path(root, name):
    if not isinstance(name, str) or not name or ':' in name or '\\' in name:
        raise ValueError('Choose a nonempty local relative filename, not a URL or drive path.')
    p = Path(name)
    if p.is_absolute() or '..' in p.parts or not (Path(root)/p).resolve().is_relative_to(Path(root).resolve()):
        raise ValueError('Input paths must stay inside the package; copy files into inputs/my_run/.')
    return Path(root)/p


def require(frame, columns, label):
    missing = set(columns)-set(frame.columns)
    if missing:
        raise ValueError(f'{label}: missing columns: '+', '.join(sorted(missing)))


def read_inputs(root, selection='selection.csv', structural='inputs/cdr_rmsd_results.csv', ablation='inputs/ablation_scores.csv'):
    """Read selected bytes once. Return actual tables for inspection; no scientific validation yet."""
    root = Path(root)
    blobs = {}
    def read(name, suffix):
        p = local_path(root, name)
        if p.suffix.lower() not in suffix:
            raise ValueError(f'{name}: expected file type {suffix}')
        if not p.is_file():
            raise ValueError(f'Missing input: {name}. Copy it locally or correct the selected path; previous exports are unchanged.')
        if p.stat().st_size > 50_000_000:
            raise ValueError(f'{name}: file exceeds the 50 MB local reader limit; ask your mentor before adapting the reader.')
        blobs[name] = p.read_bytes()
        return BytesIO(blobs[name])
    try:
        chosen = pd.read_csv(read(selection, ('.csv',)))
        require(chosen, ['structure','method','workbook','fasta'], selection)
        if set(chosen.columns) != {'structure','method','workbook','fasta'} or chosen.isna().any().any():
            raise ValueError('selection.csv: require only structure, method, workbook, fasta, without blanks')
        if len(chosen)!=4 or set(zip(chosen.structure,chosen.method)) != {(s,m) for s in analysis.CASES for m in analysis.METHODS}:
            raise ValueError('selection.csv: select exactly one baseline and ensemble for each supported case')
        books, sequences = {}, {}
        for r in chosen.itertuples():
            key = (r.structure,r.method)
            try:
                books[key] = pd.read_excel(read(r.workbook, ('.xlsx',)), sheet_name=['Overview','VH'])
            except (ValueError, OSError) as exc:
                raise ValueError(f'{r.workbook}: cannot read Overview/VH sheets: {exc}') from exc
            require(books[key]['Overview'], ['Antibody','Heavy OASis Identity'], r.workbook+' / Overview')
            text = read(r.fasta, ('.fa','.fasta')).getvalue().decode('utf-8')
            records = []
            for line in text.splitlines():
                if line.startswith('>'):
                    records.append([line[1:], ''])
                elif line.strip():
                    if not records:
                        raise ValueError(f'{r.fasta}: FASTA sequence precedes header')
                    records[-1][1] += line.strip()
            sequences[key] = records
        raw = pd.read_csv(read(structural, ('.csv',)))
        abl = None if ablation is None else pd.read_csv(read(ablation, ('.csv',)))
    except (OSError, KeyError, UnicodeError, pd.errors.ParserError) as exc:
        raise ValueError(f'Could not read selected local file: {exc}. Check format and required sheet names.') from exc
    # Original qualification is computed from a trusted package file, never an uploaded manifest.
    provenance = Path(__file__).with_name('inputs')/'provenance.json'
    pinned = {}
    if provenance.is_file() and hashlib.sha256(provenance.read_bytes()).hexdigest()==PROVENANCE_SHA256:
        pinned = {x['path']:x['sha256'] for x in json.loads(provenance.read_bytes())}
    hashes = {name:hashlib.sha256(b).hexdigest() for name,b in blobs.items()}
    original = hashes[selection] == '1a4f9010b08454bf83a3b4b5cfdb02da83e1bf5aa0569a91991311c4776e37f5' and all(pinned.get(name)==digest for name,digest in hashes.items() if name!=selection)
    loaded = dict(selection=chosen, books=books, sequences=sequences, structural=raw, ablation=abl,
                blobs=blobs, hashes=hashes, original_inputs=original, paths=dict(selection=selection,structural=structural,ablation=ablation))
    loaded['seal'] = frame_seal(loaded)
    return loaded


def frame_seal(loaded):
    """Bind previews and interpretation to the bytes parsed at intake."""
    h = hashlib.sha256()
    frames = [loaded['selection'], loaded['structural'], loaded['ablation']]
    frames += [frame for sheets in loaded['books'].values() for frame in sheets.values()]
    for frame in frames:
        if frame is not None:
            h.update(frame.to_json(orient='split', double_precision=15).encode())
    h.update(repr(loaded['sequences']).encode())
    h.update(repr(loaded['paths']).encode())
    h.update(repr(loaded['original_inputs']).encode())
    for name, blob in loaded['blobs'].items():
        h.update(name.encode()); h.update(blob)
    return h.hexdigest()


def numeric(frame, col, lo, hi, label):
    try:
        values = pd.to_numeric(frame[col], errors='raise')
    except (ValueError, TypeError):
        raise ValueError(f'{label}: Invalid {col}; expected numeric values, not text') from None
    bad = ~np.isfinite(values) | ~values.between(lo, hi)
    if bad.any():
        raise ValueError(f'{label}: Invalid {col} at rows {list(frame.index[bad][:5])}; require finite values in [{lo}, {hi}]')
    frame[col] = values


def validate_inputs(loaded):
    """Validate inspected in-memory frames, retaining independent arm sizes and all loops."""
    if loaded is None or loaded.get('seal') != frame_seal(loaded):
        raise ValueError('Inspected tables changed after reading. Save compatible input files, then restart and read again; original provenance cannot transfer.')
    rows = []
    for r in loaded['selection'].itertuples():
        key = (r.structure, r.method)
        overview = loaded['books'][key]['Overview'].copy()
        vh = loaded['books'][key]['VH'].copy()
        require(overview, ['Antibody','Heavy OASis Identity'], r.workbook+' / Overview')
        require(vh, ['Antibody','scheme','chain_type'], r.workbook+' / VH')
        for frame in (overview,vh):
            if frame.Antibody.isna().any() or frame.Antibody.duplicated().any():
                raise ValueError(f'{r.workbook}: Antibody identifiers must be present and unique')
        if not vh.chain_type.eq('H').all() or not vh.scheme.eq('kabat').all():
            raise ValueError(f'{r.workbook}: VH requires chain_type H and scheme kabat; do not mix metric definitions')
        require(overview, ['Threshold'], r.workbook+' / Overview')
        if not overview['Threshold'].eq('relaxed').all():
            raise ValueError(f'{r.workbook}: require reported relaxed OASis threshold')
        overview = overview[overview.Antibody.str.endswith('_VH')] if overview.Antibody.str.endswith('_VH').any() else overview
        seqs = {str(row['Antibody']): ''.join(str(row[c]) for c in vh.columns if c not in ['Antibody','scheme','chain_type'] and pd.notna(row[c])).replace('-','').replace(' ','') for _,row in vh.iterrows()}
        records = loaded['sequences'][key]
        if not records or any(not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', seq) for _,seq in records):
            raise ValueError(f'{r.fasta}: require nonempty canonical amino-acid sequences')
        indices = []
        for _,row in overview.iterrows():
            match = re.fullmatch(r'.*?(\d+)(?:_VH)?', str(row.Antibody))
            if not match:
                raise ValueError(f'{r.workbook}: unsupported Antibody ID {row.Antibody}')
            idx = int(match[1]) - (r.method=='baseline')
            seq = seqs.get(row.Antibody)
            if idx<0 or idx>=len(records) or not seq or seq not in records[idx][1]:
                raise ValueError(f'{r.workbook}: VH/FASTA mismatch at {row.Antibody}; correct sequence linkage first')
            indices.append(idx)
            rows.append(dict(structure=r.structure,method=r.method,index=idx,oasis=row['Heavy OASis Identity']))
        if len(indices)!=len(set(indices)) or set(indices)!=set(range(len(records))):
            raise ValueError(f'{r.workbook}: require one score per FASTA record, contiguous zero-based linked indices; no silent subset')
    bio = pd.DataFrame(rows)
    numeric(bio,'oasis',0,1,'BioPhi')
    raw = loaded['structural'].copy()
    require(raw, ['cdr','method','structure','design_id','domain','rmsd_angstrom','n_atoms_per_residue','cdr_mean_plddt','low_confidence','framework_rmsd_angstrom'], 'cdr_rmsd_results.csv')
    # Known original comparator rows are explicitly excluded, never relabelled.
    if not raw.structure.isin(analysis.CASES).all() or not raw.method.isin((*analysis.METHODS,'abmpnn')).all():
        raise ValueError('Structural CSV: unsupported structure or method; supported cases 9NH7/9NFU, methods baseline/ensemble (abmpnn excluded)')
    loops = raw[raw.method.isin(analysis.METHODS)].copy()
    ids = loops.design_id.astype(str).str.extract(r'^(9NH7|9NFU)_(baseline|ensemble)_(\d+)$')
    if ids.isna().any().any() or not ids[0].eq(loops.structure).all() or not ids[1].eq(loops.method).all():
        raise ValueError('Structural CSV: design_id must be CASE_METHOD_N and match structure/method columns; N is zero-based')
    loops['index'] = ids[2].astype(int)
    if loops.duplicated(['structure','method','index','cdr']).any():
        raise ValueError('Structural CSV: duplicate case/method/index/cdr')
    h3 = loops[loops.cdr=='H3']
    main = bio.merge(h3,on=['structure','method','index'],how='outer',validate='one_to_one',indicator=True)
    if not main._merge.eq('both').all():
        bad=main.loc[main._merge!='both',['structure','method','index','_merge']].head().to_dict('records')
        raise ValueError(f'H3 metrics and scored VH must have identical case/method/index keys: {bad}')
    for st in analysis.CASES:
        expected_loops = {'H1','H2','H3'} | ({'L1','L2','L3'} if st=='9NFU' else set())
        for method in analysis.METHODS:
            g=loops[(loops.structure==st)&(loops.method==method)]
            expected=set(bio[(bio.structure==st)&(bio.method==method)]['index'])
            if set(g.cdr)!=expected_loops or any(set(a['index'])!=expected for _,a in g.groupby('cdr')):
                raise ValueError(f'Complete loops {st}/{method}: each expected H/L loop must cover exactly the scored indices; restore missing rows, do not fill values')
    domains = np.where(loops.cdr.str.startswith('L'),'VL',loops.structure.map({'9NH7':'VHH','9NFU':'VH'}))
    if not loops.domain.eq(domains).all() or not loops.n_atoms_per_residue.eq(3).all():
        raise ValueError('Unsupported H3 endpoint or other loop endpoint: require per-domain N/CA/C backbone, VHH or VH / VL')
    for col,lo,hi in [('rmsd_angstrom',0,np.inf),('cdr_mean_plddt',0,100),('framework_rmsd_angstrom',0,np.inf)]:
        numeric(loops,col,lo,hi,'Structural CSV')
    if not loops.low_confidence.eq((loops.cdr_mean_plddt<70).astype(int)).all():
        raise ValueError('low_confidence must match mean loop pLDDT < 70; correct the flags, do not filter out the rows')
    main=bio.merge(loops[loops.cdr=='H3'],on=['structure','method','index'],validate='one_to_one')
    a=loaded['ablation']
    if a is None:
        a=pd.DataFrame(columns=['structure','condition','method','seed','Heavy OASis Identity','Light OASis Identity'])
    else:
        a=a.copy()
        require(a,['id','Heavy OASis Identity','Light OASis Identity'],'Ablation')
        meta=a.id.astype(str).str.extract(r'^(?P<structure>9NH7|9NFU)_(?P<condition>cdrh3|others)_(?P<method>baseline|ensemble)_s(?P<seed>0|[1-9]\d*)$')
        if not meta.notna().all().all() or a.id.duplicated().any() or set(meta.columns)&set(a.columns):
            raise ValueError('Ablation: require unique case_condition_method_sN identifiers; no conflicting metadata columns')
        a=pd.concat([a,meta],axis=1)
        for st in analysis.CASES:
            for cond in ('cdrh3','others'):
                g=a[(a.structure==st)&(a.condition==cond)]
                w=g.pivot(index='seed',columns='method',values='Heavy OASis Identity')
                if len(w)==0 or set(w.columns)!=set(analysis.METHODS) or w.isna().any().any():
                    raise ValueError(f'Ablation {st}/{cond}: require complete finite reported seed pairs (any positive count); restore missing partner')
        for domain,col in [('H','Heavy OASis Identity'),('L','Light OASis Identity')]:
            group=a if domain=='H' else a[a.structure=='9NFU'].copy()
            numeric(group,col,0,1,'Ablation '+domain)
        # Reuse the reviewed domain denominators and difference arithmetic.
    data=dict(main=main,loops=loops,ablation=a,hashes=loaded['hashes'],original_inputs=loaded['original_inputs'],loaded=loaded)
    analysis.summarize(data)  # also validates domain-specific ablation denominators
    return data


def status(data):
    return ('Pinned retained-input analysis. '+analysis.LINKAGE if data['original_inputs'] else
            'Changed-input analysis — original coordinate and execution audit does not transfer to changed files. Only the selected rows and scored-VH/FASTA correspondence were checked; metric lineage and biological mapping need review.')


def export_bundle(root, data, *, confirm=False):
    """Explicitly export to a unique UTC-dated folder; never overwrite earlier results."""
    if confirm is not True:
        raise ValueError('Export is opt-in: set EXPORT = True only after reviewing the current run.')
    # Revalidate the inspected raw frames; no file reread, stale PNG or previous result use.
    checked=validate_inputs(data['loaded'])
    for key in ('main','loops','ablation'):
        if not checked[key].equals(data[key]):
            raise ValueError('In-memory results changed since validation. Restart and Run All before export.')
    data = checked
    root=Path(root).resolve()
    destination=local_path(root,'exports')
    destination.mkdir(exist_ok=True)
    target=destination/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')+'-'+uuid.uuid4().hex[:8])
    target.mkdir(exist_ok=False)
    analysis.publish(target, data=data)
    for name,blob in data['loaded']['blobs'].items():
        p=local_path(target/'selected-inputs',name)
        p.parent.mkdir(parents=True,exist_ok=True); p.write_bytes(blob)
    import platform
    import importlib.metadata
    record=dict(input_sha256=data['hashes'],selection=data['loaded']['paths'],qualification=status(data),
                python=platform.python_version(),packages={n:importlib.metadata.version(n) for n in ('numpy','pandas','matplotlib','openpyxl')},
                computation='Retained-output descriptive analysis only; no inference, scoring or coordinate replay')
    (target/'input-record.json').write_text(json.dumps(record,indent=2)+'\n')
    for name in ('README.md','SOURCES.md','analysis.py','workflow.py','build_notebook.py'):
        (target/name).write_bytes((Path(__file__).parent/name).read_bytes())
    (target/'REPLAY.md').write_text('''# Replay this selected-data export

Save your executed notebook here as `results.ipynb`. In that copy set
`ROOT = Path.cwd() / "selected-inputs"`, retaining your three selected paths.
Keep EXPORT=False, then restart and run all. This recomputes from the exact
selected bytes without reading the exported figures. Do not run the default
example builder in this folder; its default paths describe the original package.

Alternatively, from this folder in Python:

```python
from pathlib import Path
import json, analysis
from workflow import read_inputs, validate_inputs
record = json.loads(Path('input-record.json').read_text())
data = validate_inputs(read_inputs(Path('selected-inputs'), **record['selection']))
table, ablation = analysis.summarize(data)
print(table)
fig = analysis.plot_main(data)
```

The input-record qualification describes the export-time check. If the trusted
original provenance file is absent on replay, qualification is conservatively
changed-input, even if numerical inputs are unchanged. That is not a numerical
failure. Never use an uploaded manifest to restore an original audit status.

`SOURCES.md` describes the original protocol and prior audit; its original-input
links refer to the original package, not proof about your replacements. Keep
that distinction when sharing this folder. The GitHub/mobile readers are static
snapshots of this export; student-notes.json records your reasoning separately.
''')
    # Reading surface links point to the executed notebook supplied by the student on save.
    (target/'WHAT_CHANGED.md').write_text('# Export of selected data\n\n'+status(data)+'\n\nInput hashes and versions: input-record.json. Exact selected bytes: selected-inputs/. Save the executed notebook alongside this export before sharing. No manuscript, causal or biological acceptance is implied.\n')
    return target
