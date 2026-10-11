"""Build the missing input/provenance interface; no inference or metric arithmetic."""
import argparse
import ast
import csv
import io
import json
import subprocess
from pathlib import Path

from validate_esmfold_predictions import (
    AA3, AUTHORITY, REQUIRED_V1, SCHEMA, canonical_bytes, design_identity,
    digest_bytes, file_sha256, input_bundle_digest, parse_pdb_residues,
    prediction_identity, read_fasta_records, rooted_path, validate_input_records,
    validate_run_contract,
)

REPO1 = 'https://github.com/sanithuheen/ProteinMPNN-Ablang-antibody-design'
REPO2 = 'https://github.com/IzzyMoon1231/antibody-sequence-design'
SOURCES = [('9NFU', 'baseline', '9NFU.fa'), ('9NFU', 'ensemble', '9NFU.fa'),
           ('9NH7', 'baseline', '9NH7_EBH.fa'), ('9NH7', 'ensemble', '9NH7.fa')]
MAP_FIELDS = ('reference_id,reference_repo,reference_commit,reference_path,reference_sha256,'
              'target,physical_chain,sequence_index0,seqres_index0,seqres_aa,'
              'reference_residue_number,reference_insertion_code,reference_atom_aa,'
              'reference_ca_present,reference_atom_seqres_agree,fill_required,'
              'current_design_mask,historical_rmsd_mask,current_domain,historical_domain,'
              'authority_status').split(',')
MANIFEST_FIELDS = ['target', 'method', 'index', 'sequence', 'pdb_path'] + sorted(
    REQUIRED_V1 - {'target', 'method', 'index', 'sequence', 'pdb_path'})


def git_commit(root):
    return subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()


def committed_source(root, relative, commit):
    path = rooted_path(root, relative)
    expected = subprocess.check_output(['git', '-C', str(root), 'show', f'{commit}:{relative}'])
    actual = path.read_bytes()
    if expected != actual:
        raise ValueError(f'Source differs from named commit: {relative}')
    return {'path': relative, 'commit': commit, 'sha256': digest_bytes(actual)}


def collect_source_evidence(repo_root, upstream_root):
    for root, expected in ((upstream_root, REPO1), (repo_root, REPO2)):
        remote = subprocess.check_output(['git', '-C', str(root), 'remote', 'get-url', 'origin'], text=True).strip().removesuffix('.git')
        if remote != expected:
            raise ValueError('Local repository origin does not match expected project')
    commits = {'repo1': git_commit(upstream_root), 'repo2': git_commit(repo_root)}
    files = {'repo1': set(), 'repo2': set()}
    for target, method, name in SOURCES:
        files['repo2'].add(f'data/inputs/week2_final/{method}/{name}')
        files['repo1'].update([f'{method}_sequences/{method}_sequences_final_filled/{name}',
                              f'{method}_sequences/{name}'])
    files['repo1'].update(['ensemble_pmpnn_ablang.py', '03_compute_rmsd.py',
                          'baseline_pmpnn_results/baseline_run_manifest.json'])
    for target in ('9NFU', '9NH7'):
        files['repo1'].update([f'RFDiffusion_structures/raw/{target}.pdb',
            f'ensemble_sequences/manifests/{target}_ensemble_manifest.json',
            f'Parameters_for_replication/{target}_run_manifest.json'])
        files['repo2'].add(f'data/pdb/{target}.pdb')
    return {'schema_version': 'uirj-input-evidence-v1', 'commits': commits,
            'generation_code_commit': None, 'generation_code_status': 'NOT_ESTABLISHED',
            'scientific_authority': AUTHORITY,
            'sources': {key: [committed_source(root, p, commits[key]) for p in sorted(files[key])]
                        for key, root in [('repo1', upstream_root), ('repo2', repo_root)]}}


def literal_assignment(path, name):
    found = [n.value for n in ast.parse(Path(path).read_text()).body
             if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in n.targets)]
    if len(found) != 1:
        raise ValueError(f'Ambiguous source assignment: {name}')
    return ast.literal_eval(found[0])


def domain_at(layout, index, length):
    hits = [domain or 'Unmapped' for domain, start, end in layout
            if start <= index < (length if end is None else end)]
    if len(hits) != 1:
        raise ValueError('Ambiguous or missing domain mapping')
    return hits[0]


def map_reference(path, target, chain, length, start_number, offset,
                  current, historical, current_layout, historical_layout, reference_commit):
    seqres = ''.join(AA3.get(a, 'X') for line in Path(path).read_text().splitlines()
                    if line.startswith('SEQRES') and line[11] == chain for a in line[19:].split())
    if len(seqres) < length + offset or 'X' in seqres:
        raise ValueError('Missing/ambiguous reference SEQRES')
    records = [r for r in parse_pdb_residues(path) if r['chain'] == chain]
    # The frozen affine mapping has no insertion codes. Preserve parser support,
    # but require a separately established map for any new inserted residue.
    if any(r['insertion_code'] != ' ' for r in records):
        raise ValueError('Insertion code needs an explicit approved coordinate map')
    if any(not start_number <= r['residue_number'] < start_number + length for r in records):
        raise ValueError('Reference residue outside frozen coordinate span')
    by_number = {r['residue_number']: r for r in records}
    sha = file_sha256(path)
    rows = []
    for index in range(length):
        residue = by_number.get(start_number + index)
        rows.append(dict(zip(MAP_FIELDS, [
            f'{target}:{sha}', REPO2, reference_commit, f'data/pdb/{target}.pdb', sha,
            target, chain, str(index), str(index + offset), seqres[index + offset],
            str(start_number + index), residue['insertion_code'] if residue else '',
            residue['aa'] if residue else '', str(bool(residue and residue['ca_present'])).lower(),
            str(residue['aa'] == seqres[index + offset]).lower() if residue else '',
            str(residue is None).lower(), str(index in current).lower(),
            str(index in historical).lower(), domain_at(current_layout, index, length),
            domain_at(historical_layout, index, length), AUTHORITY])))
    return rows


def build_residue_map(repo_root, upstream_root, evidence):
    repo_root, upstream_root = Path(repo_root), Path(upstream_root)
    configs = literal_assignment(upstream_root / 'ensemble_pmpnn_ablang.py', 'STRUCTURES')
    spans = literal_assignment(upstream_root / '03_compute_rmsd.py', 'CDR_SPANS')
    baseline = json.loads((upstream_root / 'baseline_pmpnn_results/baseline_run_manifest.json').read_text())
    result = []
    for config in sorted(configs, key=lambda c: c['name']):
        target = config['name']; chain, length, start, offset = {
            '9NFU': ('C', 247, -2, 16), '9NH7': ('E', 115, 2, 4)}[target]
        current = set(config['cdr_global_indices'])
        historical = {i for a, b in spans[target].values() for i in range(a, b)}
        old = json.loads((upstream_root / f'Parameters_for_replication/{target}_run_manifest.json').read_text())
        manifest = json.loads((upstream_root / f'ensemble_sequences/manifests/{target}_ensemble_manifest.json').read_text())
        if current != set(manifest['cdr_global_indices']) or canonical_bytes(config['domain_layout']) != canonical_bytes(manifest['domain_layout']):
            raise ValueError('Current source authorities conflict')
        if current != {p-1 for p in baseline['structures'][target]['cdr_positions_1indexed']}:
            raise ValueError('Baseline/current mask disagreement')
        if historical != set(old['cdr_global_indices']):
            raise ValueError('Historical mask authorities conflict')
        reference = repo_root / f'data/pdb/{target}.pdb'
        if reference.read_bytes() != (upstream_root / f'RFDiffusion_structures/raw/{target}.pdb').read_bytes():
            raise ValueError('Reference copies disagree')
        rows = map_reference(reference, target, chain, length, start, offset, current,
            historical, config['domain_layout'][chain], old['domain_layout'][chain], evidence['commits']['repo2'])
        expected_missing = {int(r['sequence_index0']) for r in rows if r['fill_required'] == 'true'}
        for _, method, name in [s for s in SOURCES if s[0] == target]:
            original = read_fasta_records(upstream_root / f'{method}_sequences/{name}', allow_x=True)
            final = read_fasta_records(repo_root / f'data/inputs/week2_final/{method}/{name}')
            if len(original) != len(final):
                raise ValueError('Original/final record count mismatch')
            for (header, sequence), (final_header, final_sequence) in zip(original, final):
                if header != final_header or len(sequence) != length:
                    raise ValueError('Original/final identity mismatch')
                if {i for i,a in enumerate(sequence) if a == 'X'} != expected_missing:
                    raise ValueError('Placeholder/resolved-coordinate mismatch')
                filled = ''.join(rows[i]['seqres_aa'] if a == 'X' else a for i,a in enumerate(sequence))
                if final_sequence != filled:
                    raise ValueError('Final fill lineage mismatch')
        result.extend(rows)
    return result


def csv_bytes(rows, fields):
    out = io.StringIO(newline='')
    writer = csv.DictWriter(out, fieldnames=fields, lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    return out.getvalue().encode('utf-8')


def build_manifest(repo_root, upstream_root, evidence, map_rows, map_digest, replicate_id=1):
    if type(replicate_id) is not int or replicate_id < 1:
        raise ValueError('Positive prediction replicate ID required')
    rows = []
    for target, method, name in SOURCES:
        source = f'data/inputs/week2_final/{method}/{name}'
        upstream = f'{method}_sequences/{method}_sequences_final_filled/{name}'
        records = read_fasta_records(Path(repo_root) / source)
        if len(records) != 100:
            raise ValueError('Frozen cohort requires100 source records per group')
        if (Path(repo_root) / source).read_bytes() != (Path(upstream_root) / upstream).read_bytes():
            raise ValueError('Upstream input differs')
        ref = next(r for r in map_rows if r['target'] == target)
        for ordinal, (header, sequence) in enumerate(records, 1):
            row = dict(schema_version=SCHEMA, target=target, method=method, index=str(ordinal),
                sequence=sequence, pdb_path=f'draft/{target}/{method}/{ordinal}/replicate-{replicate_id}.pdb',
                record_key=f'{target}/{method}/{ordinal}', source_header=header,
                source_design_id=str(ordinal - (method == 'ensemble')),
                source_record_ordinal=str(ordinal), source_sequence_sha256=digest_bytes(sequence.encode('ascii')),
                expected_length=str({'9NFU': 247, '9NH7': 115}[target]), source_repo=REPO2,
                source_commit=evidence['commits']['repo2'], source_path=source,
                source_file_sha256=file_sha256(Path(repo_root)/source), upstream_repo=REPO1,
                upstream_commit=evidence['commits']['repo1'], upstream_path=upstream,
                upstream_file_sha256=file_sha256(Path(upstream_root)/upstream),
                reference_id=ref['reference_id'], reference_sha256=ref['reference_sha256'],
                residue_map_sha256=map_digest, mapping_authority_status=AUTHORITY,
                planned_run_id='', prediction_replicate_id=str(replicate_id), prediction_id='')
            row['design_id'] = design_identity(row)
            rows.append(row)
    return sorted(rows, key=lambda r: (r['target'], r['method'], int(r['index'])))


def write_bundle(repo_root, upstream_root, output_root, run_contract=None, replicate_id=1):
    output = Path(output_root)
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError('Output directory must be new or empty')
    # Never place derived bundles in either source repository.
    for root in (Path(repo_root).resolve(), Path(upstream_root).resolve()):
        if output.resolve() == root or root in output.resolve().parents:
            raise ValueError('Bundle output must be outside source repositories')
    evidence = collect_source_evidence(repo_root, upstream_root)
    mapping = build_residue_map(repo_root, upstream_root, evidence)
    map_data = csv_bytes(mapping, MAP_FIELDS)
    rows = build_manifest(repo_root, upstream_root, evidence, mapping, digest_bytes(map_data), replicate_id)
    if run_contract is not None:
        # Incomplete policies remain explicit draft evidence, not invented run IDs.
        try:
            run_id = validate_run_contract(run_contract, rows)
        except (ValueError, KeyError, TypeError) as exc:
            evidence['run_contract_status'] = 'RUN_UNRESOLVED'
            evidence['run_contract_reason'] = str(exc)
        else:
            for row in rows:
                row['planned_run_id'] = run_id
                row['prediction_id'] = prediction_identity(row)
                row['pdb_path'] = row['pdb_path'].replace('draft/', f'{run_id}/', 1)
            evidence['run_contract_status'] = 'COMPLETE_NOT_EXECUTION_AUTHORIZATION'
    else:
        evidence['run_contract_status'] = 'RUN_UNRESOLVED'
    validate_input_records(rows, repo_root, upstream_root, output / 'predictions')
    evidence['input_bundle_digest'] = input_bundle_digest(rows)
    evidence['design_records'] = len(rows)
    evidence['unique_sequences'] = len({r['source_sequence_sha256'] for r in rows})
    evidence['prediction_replicate_id'] = str(replicate_id)
    eligibility = [{'record_key': r['record_key'], 'design_id': r['design_id'],
        'source_sequence_sha256': r['source_sequence_sha256'], 'prediction_id': r['prediction_id'],
        'prediction_replicate_id': r['prediction_replicate_id'],
        'input_integrity_status': 'VERIFIED', 'mapping_status': AUTHORITY,
        'prediction_presence_status': 'NOT_INSPECTED', 'sequence_validation_status': 'NOT_ASSESSED',
        'run_provenance_status': 'NOT_ASSESSED', 'serializer_contract_status': 'NOT_ASSESSED',
        'final_metric_eligibility': 'NOT_ASSESSED'} for r in rows]
    artifacts = {'prediction_manifest.csv': csv_bytes(rows, MANIFEST_FIELDS),
                 'canonical_residue_map.csv': map_data,
                 'input_provenance.json': canonical_bytes(evidence) + b'\n',
                 'prediction_eligibility.csv': csv_bytes(eligibility, list(eligibility[0]))}
    if run_contract is not None:
        artifacts['planned_run.json'] = canonical_bytes(run_contract) + b'\n'
    output.mkdir(parents=True, exist_ok=True)
    for name, data in artifacts.items():
        with (output / name).open('xb') as handle:
            handle.write(data)
    return evidence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', required=True)
    parser.add_argument('--upstream-root', required=True)
    parser.add_argument('--output-root', required=True)
    parser.add_argument('--run-contract')
    parser.add_argument('--prediction-replicate-id', type=int, default=1)
    args = parser.parse_args()
    contract = json.loads(Path(args.run_contract).read_text()) if args.run_contract else None
    evidence = write_bundle(args.repo_root, args.upstream_root, args.output_root,
                            contract, args.prediction_replicate_id)
    print(json.dumps({'design_records': evidence['design_records'],
                      'unique_sequences': evidence['unique_sequences'],
                      'run_contract_status': evidence['run_contract_status']}, sort_keys=True))


if __name__ == '__main__':
    main()
