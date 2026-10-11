
import argparse
import csv
from pathlib import Path

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

def extract_pdb_sequence(path):
    residues = []
    seen = set()
    models = 0

    with open(path, errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL "):
                models += 1
                if models > 1:
                    break
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM  "):
                continue
            if line[12:16].strip() != "CA":
                continue

            chain = line[21]
            residue_id = line[22:27]
            key = (chain, residue_id)

            if key in seen:
                continue

            seen.add(key)
            aa = AA3.get(line[17:20].strip(), "X")
            residues.append(aa)

    return "".join(residues)

def validate_sequence(expected, observed):
    if not observed:
        return "EMPTY_OR_INVALID_PDB", 0

    if len(expected) != len(observed):
        return "LENGTH_MISMATCH", sum(
            a == b for a, b in zip(expected, observed)
        )

    matching = sum(a == b for a, b in zip(expected, observed))

    if matching == len(expected):
        return "PASS", matching

    return "SEQUENCE_MISMATCH", matching

def validate_manifest(manifest_path, report_path=None):
    with open(manifest_path, newline="") as handle:
        records = list(csv.DictReader(handle))

    required = {"target", "method", "index", "sequence", "pdb_path"}
    if not records or not required.issubset(records[0]):
        raise ValueError("Manifest is empty or missing required columns")

    results = []

    for record in records:
        path = Path(record["pdb_path"])
        expected = record["sequence"].strip().upper()

        marker = path.with_name(path.name + ".complete")

        observed = ""
        matches = 0

        if not path.is_file():
            if marker.exists():
                status = "ORPHAN_COMPLETION_MARKER"
            else:
                status = "MISSING"

        elif not marker.is_file():
            status = "INCOMPLETE_NO_MARKER"

        else:
            try:
                if marker.read_text().strip() != "COMPLETE":
                    status = "INVALID_COMPLETION_MARKER"
                else:
                    observed = extract_pdb_sequence(path)
                    status, matches = validate_sequence(
                        expected, observed
                    )
            except Exception:
                status = "READ_ERROR"

        results.append({
            "target": record["target"],
            "method": record["method"],
            "index": record["index"],
            "pdb_path": str(path),
            "status": status,
            "expected_length": len(expected),
            "observed_length": len(observed),
            "matching_positions": matches,
        })

    if report_path:
        report_path = Path(report_path)
        report_path.parent.mkdir(parents=True, exist_ok=True)

        with open(report_path, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=results[0].keys()
            )
            writer.writeheader()
            writer.writerows(results)

    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--report")
    args = parser.parse_args()

    results = validate_manifest(args.manifest, args.report)

    from collections import Counter
    counts = Counter(r["status"] for r in results)

    print("ESMFOLD PREDICTION VALIDATION")
    print("=" * 50)
    print("Total records:", len(results))

    for status, count in sorted(counts.items()):
        print(f"{status}: {count}")

    if args.report:
        print("Report:", args.report)


# v1 provenance helpers are stdlib-only. PASS below always means sequence equality.
import hashlib
from collections import Counter
import json
import re
from pathlib import PurePosixPath

SCHEMA = 'uirj-prediction-input-v1'
AUTHORITY = 'UNRESOLVED_FOR_FINAL_METRICS'
AMINO_ACIDS = frozenset(AA3.values())
IDENTITY_FIELDS = ('record_key', 'source_sequence_sha256', 'output_pdb_sha256',
                   'reference_sha256', 'residue_map_sha256', 'run_contract_sha256')
REQUIRED_V1 = {
    'schema_version', 'target', 'method', 'index', 'sequence', 'pdb_path',
    'record_key', 'source_header', 'source_design_id', 'source_record_ordinal',
    'source_sequence_sha256', 'expected_length', 'source_repo', 'source_commit',
    'source_path', 'source_file_sha256', 'upstream_repo', 'upstream_commit',
    'upstream_path', 'upstream_file_sha256', 'reference_id', 'reference_sha256',
    'residue_map_sha256', 'mapping_authority_status', 'planned_run_id', 'design_id',
    'prediction_replicate_id', 'prediction_id',
}


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=True, allow_nan=False).encode('ascii')


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require_digest(value, name='digest'):
    if not isinstance(value, str) or not re.fullmatch('[0-9a-f]{64}', value):
        raise ValueError(f'Invalid {name}')
    return value


def rooted_path(root, relative):
    """Portable paths never acquire meaning from cwd or escape the bound root."""
    if root is None:
        raise ValueError('Explicit root required')
    if not isinstance(relative, str) or not relative or '\\' in relative:
        raise ValueError('Invalid relative path')
    p = PurePosixPath(relative)
    if p.is_absolute() or '..' in p.parts or ':' in relative or relative != p.as_posix():
        raise ValueError('Nonportable or escaping path')
    base = Path(root).resolve()
    result = (base / relative).resolve()
    if result == base or base not in result.parents:
        raise ValueError('Path escapes root')
    return result


def read_csv_records(path):
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Missing or duplicate CSV columns')
        rows = list(reader)
    if not rows or any(None in row or None in row.values() for row in rows):
        raise ValueError('Empty or malformed CSV')
    return rows


def read_fasta_records(path, allow_x=False):
    records, header, chunks = [], None, []
    def finish():
        if header is not None:
            seq = ''.join(chunks)
            allowed = AMINO_ACIDS | ({'X'} if allow_x else set())
            if not seq or set(seq) - allowed:
                raise ValueError(f'Invalid sequence in {path}')
            records.append((header, seq))
    for line in Path(path).read_text(encoding='ascii').splitlines():
        if line.startswith('>'):
            finish()
            header, chunks = line[1:], []
            if not header or header != header.strip():
                raise ValueError('Empty or padded FASTA header')
        elif line:
            if header is None:
                raise ValueError('Sequence before FASTA header')
            chunks.append(line)
    finish()
    if not records:
        raise ValueError('Empty FASTA')
    return records


def parse_pdb_residues(path):
    """First model, ATOM residues, explicit keys; ambiguous CA never collapses."""
    residues = {}
    model_id = 1
    model_seen = False
    for line in Path(path).read_text(encoding='ascii').splitlines():
        if line.startswith('MODEL '):
            if model_seen or residues:
                break
            model_id = int(line[10:14])
            model_seen = True
        if line.startswith('ENDMDL'):
            break
        if not line.startswith('ATOM  '):
            continue
        if len(line) < 54:
            raise ValueError('Truncated ATOM record')
        key = (line[21], int(line[22:26]), line[26])
        aa = AA3.get(line[17:20].strip(), 'X')
        item = residues.setdefault(key, {'chain': key[0], 'residue_number': key[1],
            'insertion_code': key[2], 'aa': aa, 'ca_present': False, 'model': model_id, 'ca_altloc': ''})
        if item['aa'] != aa:
            raise ValueError('Conflicting residue identities')
        if line[12:16].strip() == 'CA':
            if item['ca_present'] or line[16] not in (' ', 'A'):
                raise ValueError('Ambiguous CA assignment')
            item['ca_present'] = True
            item['ca_altloc'] = line[16]
    return list(residues.values())


def input_bundle_digest(rows):
    # Exclude run-dependent names to avoid circular input/run identity.
    return digest_bytes(canonical_bytes([
        {k: row[k] for k in sorted(REQUIRED_V1 - {'planned_run_id', 'pdb_path', 'prediction_id'})}
        for row in sorted(rows, key=lambda r: r['record_key'])]))


def validate_input_records(rows, input_root, upstream_root, output_root):
    if not rows:
        raise ValueError('Empty manifest')
    seen, sources, destinations, cached = set(), set(), set(), {}
    for row in rows:
        if not REQUIRED_V1 <= row.keys() or row['schema_version'] != SCHEMA:
            raise ValueError('Missing v1 columns or schema mismatch')
        target, method = row['target'], row['method']
        if target not in {'9NFU', '9NH7'} or method not in {'baseline', 'ensemble'}:
            raise ValueError('Unknown target/method')
        index = row['index']
        if not re.fullmatch('[1-9][0-9]*', index):
            raise ValueError('Index must be positive canonical integer')
        if row['source_record_ordinal'] != index:
            raise ValueError('Source ordinal differs from design index')
        key = f'{target}/{method}/{index}'
        if row['record_key'] != key or key in seen:
            raise ValueError('Duplicate or inconsistent record key')
        seen.add(key)
        if row['design_id'] != design_identity(row):
            raise ValueError('Design identity mismatch')
        if not re.fullmatch('[1-9][0-9]*', row['prediction_replicate_id']):
            raise ValueError('Invalid prediction replicate identity')
        seq = row['sequence']
        length = {'9NFU': 247, '9NH7': 115}[target]
        if set(seq) - AMINO_ACIDS or len(seq) != length or row['expected_length'] != str(length):
            raise ValueError('Invalid sequence or expected length')
        if digest_bytes(seq.encode('ascii')) != row['source_sequence_sha256']:
            raise ValueError('Sequence hash mismatch')
        for field in ('source_sequence_sha256', 'source_file_sha256', 'upstream_file_sha256',
                      'reference_sha256', 'residue_map_sha256'):
            require_digest(row[field], field)
        if row['reference_id'] != f"{target}:{row['reference_sha256']}":
            raise ValueError('Reference identity mismatch')
        for prefix, root in [('source', input_root), ('upstream', upstream_root)]:
            if not row[f'{prefix}_repo'] or not re.fullmatch('[0-9a-f]{40}', row[f'{prefix}_commit']):
                raise ValueError('Missing repository snapshot identity')
            path = rooted_path(root, row[f'{prefix}_path'])
            if path not in cached:
                cached[path] = (file_sha256(path), read_fasta_records(path))
            sha, records = cached[path]
            if sha != row[f'{prefix}_file_sha256']:
                raise ValueError('Source file hash mismatch')
            if int(index) > len(records) or records[int(index)-1] != (row['source_header'], seq):
                raise ValueError('Source record/header mismatch')
        if row['source_file_sha256'] != row['upstream_file_sha256']:
            raise ValueError('Upstream copy differs')
        pattern = r'(?:^|[, ])sample=(\d+)(?:[, ]|$)' if method == 'baseline' else r'(?:^| )seed=(\d+)(?: |$)'
        matches = re.findall(pattern, row['source_header'])
        expected_id = str(int(index) - (method == 'ensemble'))
        if matches != [expected_id] or row['source_design_id'] != expected_id:
            raise ValueError('Ambiguous or inconsistent source design ID')
        source_key = (row['source_repo'], row['source_commit'], row['source_path'], row['source_design_id'])
        if source_key in sources:
            raise ValueError('Duplicate source identity')
        sources.add(source_key)
        if row['mapping_authority_status'] != AUTHORITY:
            raise ValueError('Mapping scientific authority not established by v1')
        run = row['planned_run_id']
        if run:
            require_digest(run, 'planned_run_id')
        namespace = run or 'draft'
        if row['prediction_id'] != prediction_identity(row):
            raise ValueError('Prediction identity mismatch')
        replica = row['prediction_replicate_id']
        if row['pdb_path'] != f'{namespace}/{target}/{method}/{index}/replicate-{replica}.pdb':
            raise ValueError('Prediction namespace does not match identity')
        destination = rooted_path(output_root, row['pdb_path'])
        if destination in destinations:
            raise ValueError('Prediction destination collision')
        destinations.add(destination)
    if len({r['residue_map_sha256'] for r in rows}) != 1 or len({r['planned_run_id'] for r in rows}) != 1:
        raise ValueError('Mixed map or run identities')
    return rows


def validate_residue_map(path, rows, input_root, upstream_root):
    expected = {r['residue_map_sha256'] for r in rows}
    if expected != {file_sha256(path)}:
        raise ValueError('Residue map hash mismatch')
    mappings = read_csv_records(path)
    keys, counts = set(), Counter()
    for item in mappings:
        target = item['target']
        if target not in {'9NFU', '9NH7'}:
            raise ValueError('Unknown map target')
        index = int(item['sequence_index0'])
        key = (target, index)
        if key in keys or not 0 <= index < {'9NFU': 247, '9NH7': 115}[target]:
            raise ValueError('Duplicate/out-of-range map key')
        keys.add(key); counts[target] += 1
        refs = {r['reference_sha256'] for r in rows if r['target'] == target}
        if refs and refs != {item['reference_sha256']}:
            raise ValueError('Map/reference join mismatch')
        ref = rooted_path(input_root, item['reference_path'])
        if file_sha256(ref) != item['reference_sha256']:
            raise ValueError('Reference hash mismatch')
        if item['authority_status'] != AUTHORITY:
            raise ValueError('Unapproved map authority')
    if counts != {'9NFU': 247, '9NH7': 115}:
        raise ValueError('Incomplete canonical map')
    # Reconstruct this small coordinate table from pinned source facts. A map's
    # self-consistent hash alone is insufficient semantic evidence.
    from build_prediction_manifest import collect_source_evidence, build_residue_map, REPO1, REPO2
    evidence = collect_source_evidence(input_root, upstream_root)
    for row in rows:
        if (row['source_repo'], row['upstream_repo']) != (REPO2, REPO1):
            raise ValueError('Repository identity mismatch')
        if row['source_commit'] != evidence['commits']['repo2'] or row['upstream_commit'] != evidence['commits']['repo1']:
            raise ValueError('Repository snapshot mismatch')
    reconstructed = build_residue_map(input_root, upstream_root, evidence)
    order = lambda r: (r['target'], int(r['sequence_index0']))
    if sorted(mappings, key=order) != sorted(reconstructed, key=order):
        raise ValueError('Canonical map differs from source mapping evidence')
    return mappings


def load_input_manifest(path, input_root=None, upstream_root=None, output_root=None,
                        legacy=False):
    rows = read_csv_records(path)
    if legacy:
        if any(r.get('schema_version') for r in rows):
            raise ValueError('v1 cannot bypass validation in legacy mode')
        if not {'target', 'method', 'index', 'sequence', 'pdb_path'} <= rows[0].keys():
            raise ValueError('Manifest missing required columns')
        return rows
    validate_input_records(rows, input_root, upstream_root, output_root)
    from build_prediction_manifest import SOURCES
    expected_groups = {(target, method): name for target, method, name in SOURCES}
    for (target, method), name in expected_groups.items():
        group = [row for row in rows if (row['target'], row['method']) == (target, method)]
        if {row['index'] for row in group} != {str(i) for i in range(1, 101)}:
            raise ValueError('Incomplete frozen sampled-design cohort')
        for row in group:
            if row['source_path'] != f'data/inputs/week2_final/{method}/{name}' or row['upstream_path'] != f'{method}_sequences/{method}_sequences_final_filled/{name}':
                raise ValueError('Manifest source is not the canonical final input')
    validate_residue_map(Path(path).parent / 'canonical_residue_map.csv', rows,
                         input_root, upstream_root)
    return sorted(rows, key=lambda r: (r['target'], r['method'], int(r['index'])))


RUN_FIELDS = {'schema_version', 'model_id', 'model_revision', 'assets',
              'serializer_sha256', 'runtime', 'device', 'precision',
              'inference_settings', 'seed', 'deterministic_algorithms', 'input_bundle_digest', 'run_instance_id'}


def validate_run_contract(contract, rows):
    if set(contract) != RUN_FIELDS or contract['schema_version'] != 'uirj-run-v1':
        raise ValueError('Incomplete/unknown run contract fields')
    if contract['model_id'] != 'facebook/esmfold_v1' or not re.fullmatch('[0-9a-f]{40}', contract['model_revision']):
        raise ValueError('Model revision must be pinned')
    if not re.fullmatch('[A-Za-z0-9_-]{1,80}', contract['run_instance_id']):
        raise ValueError('Explicit run instance identity required')
    require_digest(contract['serializer_sha256'])
    if contract['input_bundle_digest'] != input_bundle_digest(rows):
        raise ValueError('Run input bundle mismatch')
    if contract['device'] != 'cuda' or contract['precision'] != {'esm': 'float16', 'trunk': 'float32'}:
        raise ValueError('Unsupported execution settings')
    settings = contract['inference_settings']
    if set(settings) != {'num_recycles', 'chunk_size'} or any(type(v) is not int or v < 1 for v in settings.values()):
        raise ValueError('Explicit positive recycle/chunk settings required')
    if type(contract['seed']) is not int or not 0 <= contract['seed'] < 2**32 or type(contract['deterministic_algorithms']) is not bool:
        raise ValueError('Invalid seed/determinism policy')
    if set(contract['runtime']) != {'python', 'torch', 'transformers'} or not all(isinstance(v, str) and v for v in contract['runtime'].values()):
        raise ValueError('Explicit runtime versions required')
    if not isinstance(contract['assets'], dict) or not contract['assets']:
        raise ValueError('Pinned model asset hashes required')
    for path, sha in contract['assets'].items():
        rooted_path('/tmp/uirj-assets-validation', path)
        require_digest(sha)
    return digest_bytes(canonical_bytes(contract))


def validate_receipt(receipt, row, contract, rows, path):
    run = validate_run_contract(contract, rows)
    expected = {'schema_version': 'uirj-prediction-receipt-v1',
                'record_key': row['record_key'], 'design_id': row['design_id'],
                'prediction_id': row['prediction_id'], 'prediction_replicate_id': row['prediction_replicate_id'],
                'source_sequence_sha256': row['source_sequence_sha256'],
                'input_bundle_digest': input_bundle_digest(rows),
                'run_contract_sha256': run,
                'output_pdb_sha256': file_sha256(path),
                'serializer_sha256': contract['serializer_sha256'],
                'runtime': contract['runtime']}
    if row['planned_run_id'] != run:
        raise ValueError('Manifest run identity mismatch')
    if set(receipt) != set(expected) | {'completed_at'} or not receipt['completed_at']:
        raise ValueError('Malformed prediction receipt')
    if any(receipt[k] != v for k, v in expected.items()):
        raise ValueError('Prediction receipt identity/hash mismatch')
    return receipt


def check_prediction_status(path, sequence, strict=False):
    path = Path(path); marker = path.with_name(path.name + '.complete')
    if not path.is_file():
        return 'ORPHAN_COMPLETION_MARKER' if marker.exists() else 'MISSING'
    if not marker.is_file():
        return 'INCOMPLETE_NO_MARKER'
    try:
        if marker.read_text().strip() != 'COMPLETE':
            return 'INVALID_COMPLETION_MARKER'
        if strict:
            residues = parse_pdb_residues(path)
            if len({r['chain'] for r in residues}) != 1 or any(not r['ca_present'] for r in residues):
                return 'AMBIGUOUS_PREDICTION_MAPPING'
            if [(r['residue_number'], r['insertion_code']) for r in residues] != [(i, ' ') for i in range(1, len(sequence)+1)]:
                return 'UNRESOLVED_PREDICTION_MAPPING'
            observed = ''.join(r['aa'] for r in residues)
        else:
            observed = extract_pdb_sequence(path)
        return validate_sequence(sequence, observed)[0]
    except Exception:
        return 'READ_ERROR'


def prediction_eligibility(row, output_root, rows, contract=None, *, input_root=None, upstream_root=None, manifest_path=None):
    if manifest_path is None:
        raise ValueError('Eligibility requires verified source/manifest context')
    checked = load_input_manifest(manifest_path, input_root, upstream_root, output_root)
    by_key = {r['record_key']: r for r in checked}
    if row != by_key.get(row['record_key']) or input_bundle_digest(rows) != input_bundle_digest(checked):
        raise ValueError('Eligibility input context mismatch')
    return _prediction_eligibility(row, output_root, rows, contract)


def _prediction_eligibility(row, output_root, rows, contract=None):
    # Internal only: caller has verified source identity/map semantics.
    path = rooted_path(output_root, row['pdb_path'])
    status = check_prediction_status(path, row['sequence'], strict=True)
    result = {'record_key': row['record_key'], 'design_id': row['design_id'],
              'prediction_id': row['prediction_id'], 'prediction_replicate_id': row['prediction_replicate_id'], 'source_sequence_sha256': row['source_sequence_sha256'],
              'input_integrity_status': 'VERIFIED', 'mapping_status': row['mapping_authority_status'],
              'prediction_presence_status': status, 'sequence_validation_status': status,
              'run_provenance_status': 'NOT_ASSESSED', 'serializer_contract_status': 'NOT_ASSESSED',
              'final_metric_eligibility': 'NOT_ASSESSED', 'output_pdb_sha256': '',
              'reference_sha256': row['reference_sha256'], 'residue_map_sha256': row['residue_map_sha256'],
              'run_contract_sha256': row['planned_run_id']}
    if status != 'PASS':
        return result
    result['output_pdb_sha256'] = file_sha256(path)
    sidecar = path.with_name(path.name + '.provenance.json')
    if contract is None or not sidecar.is_file():
        result['run_provenance_status'] = 'RUN_UNVERIFIED'
        return result
    try:
        validate_receipt(json.loads(sidecar.read_text()), row, contract, rows, path)
    except (ValueError, KeyError, TypeError, OSError):
        result['run_provenance_status'] = 'RECEIPT_MISMATCH'
    else:
        result['run_provenance_status'] = 'VERIFIED'
        result['serializer_contract_status'] = 'IDENTITY_VERIFIED_NOT_GEOMETRY_VALIDATED'
    return result


def validate_measurement_joins(measurements, identities):
    """Identity-only P1c seam. No arithmetic, metric approval or scientific inference."""
    index = {}
    for item in identities:
        if item.get('run_provenance_status') != 'VERIFIED' or item.get('sequence_validation_status') != 'PASS':
            raise ValueError('Unverified prediction identity cannot support measurement join')
        for field in ('design_id', 'prediction_id', 'source_sequence_sha256', 'output_pdb_sha256', 'reference_sha256', 'residue_map_sha256', 'run_contract_sha256'):
            require_digest(item.get(field), field)
        key = (item['design_id'], item['prediction_id'])
        if key in index:
            raise ValueError('Ambiguous design identity join')
        index[key] = item
    seen = set()
    for measurement in measurements:
        required = set(IDENTITY_FIELDS) | {'measurement_id', 'design_id', 'prediction_id', 'metric_policy_id', 'extraction_code_id', 'metric', 'region'}
        if not required <= measurement.keys() or not all(measurement[k] for k in required):
            raise ValueError('Missing measurement identity')
        key = (measurement['design_id'], measurement['prediction_id'])
        if key not in index:
            raise ValueError('Unmatched measurement design')
        identity = index[key]
        if any(measurement[k] != identity[k] for k in IDENTITY_FIELDS):
            raise ValueError('Measurement prediction/input/map/run mismatch')
        unique = (key, measurement['metric'], measurement['region'], measurement['metric_policy_id'], measurement['extraction_code_id'])
        if measurement['measurement_id'] in seen or unique in seen:
            raise ValueError('Duplicate or ambiguous measurement identity')
        seen.update((measurement['measurement_id'], unique))
    return sorted(measurements, key=lambda m: (m['design_id'], m['prediction_id'], m['measurement_id']))


def design_identity(row):
    fields = ('record_key', 'source_repo', 'source_commit', 'source_path', 'source_file_sha256')
    return digest_bytes(canonical_bytes({k: row[k] for k in fields}))


def prediction_identity(row):
    if not row['planned_run_id']:
        return ''  # no fictional prediction identity for an unresolved run
    return digest_bytes(canonical_bytes({k: row[k] for k in
        ('design_id', 'planned_run_id', 'prediction_replicate_id')}))


def validate_run_registry(contracts, rows):
    by_instance = {}
    for contract in contracts:
        digest = validate_run_contract(contract, rows)
        instance = contract['run_instance_id']
        if instance in by_instance:
            raise ValueError('Duplicate/conflicting run registry identity')
        by_instance[instance] = digest
    return by_instance


def validate_prediction_bundle(manifest_path, input_root, upstream_root, output_root, contract=None):
    """Public P1c seam: source/map-verified records plus freshly observed statuses."""
    rows = load_input_manifest(manifest_path, input_root, upstream_root, output_root)
    if contract is not None:
        validate_run_contract(contract, rows)
    return {'generation_code_commit': None, 'generation_code_status': 'NOT_ESTABLISHED',
            'mapping_authority_status': AUTHORITY, 'records': rows, 'eligibility': [
        _prediction_eligibility(row, output_root, rows, contract) for row in rows]}


def provenance_main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--report')
    parser.add_argument('--input-root')
    parser.add_argument('--upstream-root')
    parser.add_argument('--output-root')
    parser.add_argument('--run-contract')
    parser.add_argument('--provenance-mode', choices=['strict', 'legacy'], default='strict')
    args = parser.parse_args()
    if args.provenance_mode == 'legacy':
        load_input_manifest(args.manifest, legacy=True)
        results = validate_manifest(args.manifest, args.report)
    else:
        rows = load_input_manifest(args.manifest, args.input_root, args.upstream_root, args.output_root)
        contract = json.loads(Path(args.run_contract).read_text()) if args.run_contract else None
        if contract is not None:
            validate_run_contract(contract, rows)
        results = [_prediction_eligibility(r, args.output_root, rows, contract) for r in rows]
        if args.report:
            with open(args.report, 'x', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(results[0]))
                writer.writeheader(); writer.writerows(results)
    print(json.dumps(results, sort_keys=True))


if __name__ == '__main__':
    provenance_main()
