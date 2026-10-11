"""Eligibility contracts only: no scientific arithmetic or metric extraction.

P1b owns input/map/run/prediction identities. This module adds metric receipt and
analysis-policy relationships; VERIFIED means contract checks, not scientific truth.
"""
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import validate_esmfold_predictions as core

STATES = ('CONFLICTING', 'INELIGIBLE', 'HISTORICAL', 'MISSING', 'UNRESOLVED', 'VERIFIED')
CASES = {'9NFU', '9NH7'}
METHODS = {'baseline', 'ensemble'}
PLAN_SCHEMA = 'uirj-statistics-plan-v1'
RECEIPT_SCHEMA = 'uirj-metric-receipts-v1'
METRIC_UNITS = {'rmsd': {'angstrom'}, 'plddt': {'confidence_0_100'},
                'contact_displacement': {'angstrom'}, 'contact_identity': {'fraction', 'percent'},
                'oasis': {'fraction', 'percent'}}
UNIT_DOMAINS = {'angstrom': (0, None), 'confidence_0_100': (0, 100),
                'fraction': (0, 1), 'percent': (0, 100)}
SOURCE_ID_FIELDS = set(core.IDENTITY_FIELDS) | {
    'design_id', 'prediction_id', 'measurement_id', 'endpoint_id', 'metric', 'region',
    'unit', 'metric_policy_id', 'extraction_code_id', 'extraction_run_id'}


def issue(state, code, key='', detail=''):
    return {'state': state, 'reason_code': code, 'key': str(key), 'detail': str(detail)}


def normalized_issues(items):
    return [dict(zip(('state', 'reason_code', 'key', 'detail'), values))
            for values in sorted({tuple(i[k] for k in ('state', 'reason_code', 'key', 'detail')) for i in items})]


def evidence_state(items):
    return next((state for state in STATES if any(i['state'] == state for i in items)), 'VERIFIED')


def read_json(path):
    with open(path) as handle:
        return json.load(handle)


def diagnostic_key(value):
    # Ordering malformed evidence is not an identity hash; preserve nonfinite errors.
    return json.dumps(value, sort_keys=True, allow_nan=True)


def policy_digest(value):
    return core.digest_bytes(core.canonical_bytes(value))


def metric_configuration(contract):
    # Run instance and cohort/replicate membership change identity, not settings.
    return policy_digest({k: v for k, v in contract.items()
                          if k not in {'run_instance_id', 'input_bundle_digest'}})


def approval_issues(plan, approval, approved_plan_sha256):
    digest = policy_digest(plan)
    if not approval or not approved_plan_sha256:
        return [issue('UNRESOLVED', 'SCIENTIFIC_APPROVAL_NOT_SUPPLIED')]
    required = {'schema_version', 'plan_sha256', 'decision_id', 'decision_source',
                'decision_authority', 'decision_timestamp'}
    if set(approval) != required or approval.get('schema_version') != 'uirj-plan-approval-v1':
        return [issue('CONFLICTING', 'MALFORMED_APPROVAL_RECEIPT')]
    if digest != approval['plan_sha256'] or digest != approved_plan_sha256:
        return [issue('CONFLICTING', 'APPROVAL_PLAN_HASH_MISMATCH')]
    if any(not isinstance(approval[k], str) or not approval[k].strip() for k in required):
        return [issue('UNRESOLVED', 'APPROVAL_AUTHORITY_EVIDENCE_INCOMPLETE')]
    return []  # Operator attestation binding only; no human authentication claim.


def validate_plan(plan, designs, evidence_root):
    """Return all known policy blockers. Do not choose missing scientific defaults."""
    problems = []
    required = {'schema_version', 'analysis_scope', 'cases', 'methods', 'design_ids',
                'endpoints', 'evaluation_members', 'hypotheses', 'families',
                'generation_interpretation', 'analysis_timing'}
    if not isinstance(plan, dict) or not required <= plan.keys():
        return [issue('MISSING', 'ANALYSIS_PLAN_FIELDS_MISSING')]
    if plan['schema_version'] != PLAN_SCHEMA:
        problems.append(issue('INELIGIBLE', 'UNSUPPORTED_PLAN_SCHEMA'))
    if plan['analysis_scope'] != 'fixed_cases' or set(plan['cases']) != CASES or set(plan['methods']) != METHODS:
        problems.append(issue('INELIGIBLE', 'CASE_OR_COMPARISON_SCOPE_UNSUPPORTED'))
    if plan['analysis_timing'] not in {'retrospective_exploratory', 'prospective_new_experiment'}:
        problems.append(issue('UNRESOLVED', 'ANALYSIS_TIMING_UNRESOLVED'))
    if plan['generation_interpretation'] != 'implemented_protocol_contrast':
        problems.append(issue('UNRESOLVED', 'CONTROL_EQUIVALENCE_NOT_ESTABLISHED'))
    if len(plan['design_ids']) != len(set(plan['design_ids'])):
        problems.append(issue('CONFLICTING', 'DUPLICATE_PLANNED_DESIGN'))
    if set(plan['design_ids']) != set(designs):
        problems.append(issue('CONFLICTING', 'FROZEN_DESIGN_POPULATION_MISMATCH'))
    if any(not isinstance(plan[k], list) or any(not isinstance(x, dict) for x in plan[k])
           for k in ('endpoints', 'evaluation_members', 'hypotheses', 'families')):
        raise ValueError('Plan collections must be lists of objects')
    endpoints = {}
    required_endpoint = {'endpoint_id', 'case', 'metric', 'region', 'kind', 'unit', 'direction',
        'bounds', 'metric_policy_id', 'extraction_code_id', 'reference_sha256', 'residue_map_sha256',
        'definition', 'estimands', 'effect', 'uncertainty', 'replication', 'missingness', 'pairing'}
    endpoint_counts = Counter(e.get('endpoint_id') for e in plan['endpoints'])
    for endpoint in sorted(plan['endpoints'], key=diagnostic_key):
        key = endpoint.get('endpoint_id', '<missing>')
        if not required_endpoint <= endpoint.keys():
            problems.append(issue('MISSING', 'ENDPOINT_FIELDS_MISSING', key)); continue
        if endpoint_counts[key] > 1:
            problems.append(issue('CONFLICTING', 'DUPLICATE_ENDPOINT_ID', key)); continue
        endpoints[key] = endpoint
        if any(not isinstance(endpoint[k], dict) for k in
               ('bounds', 'definition', 'estimands', 'effect', 'uncertainty', 'replication', 'missingness', 'pairing')):
            raise ValueError('Nested endpoint policies must be objects')
        if endpoint['case'] not in CASES or endpoint['kind'] not in {'prediction', 'sequence_score'}:
            problems.append(issue('INELIGIBLE', 'ENDPOINT_SCOPE_UNSUPPORTED', key))
        if any(not isinstance(endpoint[k], str) or not endpoint[k].strip() for k in
               ('endpoint_id', 'metric', 'region', 'unit', 'metric_policy_id', 'extraction_code_id')):
            problems.append(issue('UNRESOLVED', 'ENDPOINT_DEFINITION_UNRESOLVED', key))
        if endpoint['metric'] not in METRIC_UNITS or endpoint['unit'] not in METRIC_UNITS[endpoint['metric']]:
            problems.append(issue('INELIGIBLE', 'UNSUPPORTED_METRIC_UNIT', key))
        if endpoint['direction'] not in {'higher', 'lower'}:
            problems.append(issue('UNRESOLVED', 'METRIC_DIRECTION_UNRESOLVED', key))
        bounds = endpoint['bounds']
        if set(bounds) != {'lower', 'upper'} or any(v is not None and
                (type(v) not in (int, float) or not math.isfinite(v)) for v in bounds.values()):
            problems.append(issue('CONFLICTING', 'INVALID_METRIC_BOUNDS', key))
        elif bounds['lower'] is not None and bounds['upper'] is not None and bounds['lower'] > bounds['upper']:
            problems.append(issue('CONFLICTING', 'INVALID_METRIC_BOUNDS', key))
        domain = UNIT_DOMAINS.get(endpoint['unit'])
        if domain and set(bounds) == {'lower', 'upper'}:
            if any(bound is not None and (bounds[name] is None or type(bounds[name]) not in (int, float))
                   for name, bound in zip(('lower', 'upper'), domain)):
                problems.append(issue('UNRESOLVED', 'INHERENT_UNIT_BOUNDS_UNRESOLVED', key))
            elif ((domain[0] is not None and bounds['lower'] < domain[0]) or
                  (domain[1] is not None and bounds['upper'] > domain[1])):
                problems.append(issue('INELIGIBLE', 'PLAN_BOUNDS_EXCEED_UNIT_DOMAIN', key))
        for field in ('reference_sha256', 'residue_map_sha256'):
            try: core.require_digest(endpoint[field], field)
            except ValueError: problems.append(issue('UNRESOLVED', 'METRIC_DEFINITION_HASH_UNRESOLVED', key, field))
        definition = endpoint['definition']
        for field in ('atom_selection', 'alignment', 'domain', 'mask_source', 'score_definition'):
            if not definition.get(field) or definition[field] == 'UNRESOLVED':
                problems.append(issue('UNRESOLVED', 'METRIC_DEFINITION_UNRESOLVED', key, field))
        estimands = endpoint['estimands']
        if (estimands.get('primary') != 'draw_weighted' or
            estimands.get('sensitivity') != 'unique_sequence_support' or
            estimands.get('unique_support_rule') != 'retain_all_evaluation_members' or
            estimands.get('unique_support_weight') != 'equal_per_sequence_within_case_method' or
            estimands.get('within_sequence_rule') != 'mean_over_distinct_prediction_ids'):
            problems.append(issue('UNRESOLVED', 'ESTIMAND_POLICY_UNSUPPORTED_OR_UNRESOLVED', key))
        if endpoint['effect'] != {'statistic': 'mean_difference', 'location': 'ensemble_minus_baseline', 'sign': 'raw_metric_units'}:
            problems.append(issue('UNRESOLVED', 'EFFECT_DEFINITION_UNSUPPORTED_OR_UNRESOLVED', key))
        uncertainty = endpoint['uncertainty']
        if uncertainty.get('method') == 'descriptive_only':
            if uncertainty.get('unit') != 'none' or uncertainty.get('hierarchy') != 'none':
                problems.append(issue('CONFLICTING', 'DESCRIPTIVE_UNCERTAINTY_CONTRADICTION', key))
        elif uncertainty.get('method') == 'bootstrap':
            if (uncertainty.get('unit') not in {'draw', 'pair'} or not uncertainty.get('hierarchy') or
                type(uncertainty.get('resamples')) is not int or uncertainty['resamples'] < 1 or
                type(uncertainty.get('seed')) is not int or
                type(uncertainty.get('confidence')) not in (int, float) or not 0 < uncertainty['confidence'] < 1 or
                not uncertainty.get('interval_algorithm') or not uncertainty.get('quantile_rule')):
                problems.append(issue('UNRESOLVED', 'UNCERTAINTY_CONTRACT_INCOMPLETE', key))
        else:
            problems.append(issue('UNRESOLVED', 'UNCERTAINTY_METHOD_UNSUPPORTED_OR_UNRESOLVED', key))
        replication = endpoint['replication']
        if replication.get('rule') not in {'single_per_draw', 'mean_within_draw', 'shared_deterministic'}:
            problems.append(issue('UNRESOLVED', 'REPLICATE_POLICY_UNRESOLVED', key))
        if replication.get('assignment') not in {'per_draw', 'shared_sequence'} or not replication.get('configuration_id'):
            problems.append(issue('UNRESOLVED', 'REPLICATE_CONFIGURATION_UNRESOLVED', key))
        if (replication.get('rule') == 'shared_deterministic') != (replication.get('assignment') == 'shared_sequence'):
            problems.append(issue('CONFLICTING', 'REUSE_REPLICATION_POLICY_CONTRADICTION', key))
        missingness = endpoint['missingness']
        if missingness.get('rule') not in {'block', 'descriptive_valid_subset'}:
            problems.append(issue('UNRESOLVED', 'MISSINGNESS_POLICY_UNRESOLVED', key))
        if missingness.get('rule') == 'descriptive_valid_subset' and (
                missingness.get('estimand') != 'observed_valid_subset' or uncertainty.get('method') != 'descriptive_only'):
            problems.append(issue('CONFLICTING', 'MISSINGNESS_CHANGES_ESTIMAND_WITHOUT_DECLARATION', key))
        if missingness.get('rule') == 'block' and missingness.get('estimand') != 'full_planned_population':
            problems.append(issue('CONFLICTING', 'MISSINGNESS_ESTIMAND_CONTRADICTION', key))
        pairing = endpoint['pairing']
        if pairing.get('mode') == 'unpaired':
            if set(pairing) != {'mode'} or uncertainty.get('unit') == 'pair':
                problems.append(issue('CONFLICTING', 'UNPAIRED_PLAN_HAS_PAIRING', key))
        elif pairing.get('mode') == 'paired':
            if pairing.get('basis') != 'prespecified_coupling':
                problems.append(issue('INELIGIBLE', 'NOMINAL_ID_OR_SEED_PAIRING_REJECTED', key))
            else:
                try:
                    path = core.rooted_path(evidence_root, pairing['evidence_path'])
                    if core.file_sha256(path) != pairing['evidence_sha256']:
                        raise ValueError('Coupling evidence hash mismatch')
                    receipt = read_json(path)
                    if receipt != {'coupling_protocol_id': pairing['coupling_protocol_id'], 'pairs': pairing['pairs']}:
                        raise ValueError('Coupling evidence content mismatch')
                    paired, coupling_ids = set(), set()
                    for pair in pairing['pairs']:
                        left, right = pair['baseline_design_id'], pair['ensemble_design_id']
                        if not pair['coupling_id'] or pair['coupling_id'] in coupling_ids or left in paired or right in paired:
                            raise ValueError('Duplicate/ambiguous coupling')
                        if designs[left]['method'] != 'baseline' or designs[right]['method'] != 'ensemble' or any(designs[d]['target'] != endpoint['case'] for d in (left, right)):
                            raise ValueError('Coupling group mismatch')
                        paired.update((left, right)); coupling_ids.add(pair['coupling_id'])
                    expected = {d for d, r in designs.items() if r['target'] == endpoint['case']}
                    if paired != expected:
                        raise ValueError('Incomplete coupling coverage')
                except (KeyError, OSError, TypeError, ValueError) as exc:
                    problems.append(issue('UNRESOLVED', 'PAIRING_EVIDENCE_UNVERIFIED', key, exc))
        else:
            problems.append(issue('UNRESOLVED', 'PAIRING_MODE_UNRESOLVED', key))
        if uncertainty.get('method') == 'bootstrap':
            if pairing.get('mode') == 'paired':
                if uncertainty.get('unit') != 'pair' or uncertainty.get('hierarchy') != 'paired_draws' or replication.get('rule') != 'single_per_draw':
                    problems.append(issue('UNRESOLVED', 'PAIRED_UNCERTAINTY_HIERARCHY_UNSUPPORTED', key))
            else:
                expected_hierarchy = {'single_per_draw':'independent_draws',
                    'mean_within_draw':'draws_then_prediction_replicates'}.get(replication.get('rule'))
                if uncertainty.get('unit') != 'draw' or uncertainty.get('hierarchy') != expected_hierarchy or expected_hierarchy is None:
                    problems.append(issue('UNRESOLVED', 'UNCERTAINTY_REPLICATION_HIERARCHY_UNSUPPORTED', key))
    if {e.get('case') for e in endpoints.values()} != CASES:
        problems.append(issue('MISSING', 'CASE_ENDPOINT_COVERAGE_MISSING'))
    hypothesis_ids, family_ids = set(), set()
    hypothesis_counts = Counter(h.get('hypothesis_id') for h in plan['hypotheses'])
    for hypothesis in sorted(plan['hypotheses'], key=diagnostic_key):
        hid = hypothesis.get('hypothesis_id')
        if not hid or hypothesis_counts[hid] > 1:
            problems.append(issue('CONFLICTING', 'DUPLICATE_OR_MISSING_HYPOTHESIS', hid)); continue
        hypothesis_ids.add(hid)
        endpoint = endpoints.get(hypothesis.get('endpoint_id'))
        if not endpoint or hypothesis.get('case') != endpoint.get('case') or hypothesis.get('comparison') != ['ensemble', 'baseline']:
            problems.append(issue('CONFLICTING', 'HYPOTHESIS_ENDPOINT_COMPARISON_MISMATCH', hid))
    family_counts = Counter(f.get('family_id') for f in plan['families'])
    for family in sorted(plan['families'], key=diagnostic_key):
        fid = family.get('family_id')
        if not fid or family_counts[fid] > 1:
            problems.append(issue('CONFLICTING', 'DUPLICATE_OR_MISSING_FAMILY', fid)); continue
        family_ids.add(fid)
        expected = {h['hypothesis_id'] for h in plan['hypotheses'] if h.get('family_id') == fid}
        declared = family.get('hypothesis_ids', [])
        if not expected or set(declared) != expected or len(declared) != len(set(declared)):
            problems.append(issue('CONFLICTING', 'HYPOTHESIS_FAMILY_MEMBERSHIP_MISMATCH', fid))
        if family.get('correction') not in {'descriptive_only', 'holm', 'bh'} or family.get('missingness') != 'block':
            problems.append(issue('UNRESOLVED', 'MULTIPLICITY_POLICY_UNRESOLVED', fid))
        if family.get('correction') == 'bh' and not family.get('dependence_justification'):
            problems.append(issue('UNRESOLVED', 'BH_DEPENDENCE_JUSTIFICATION_MISSING', fid))
        if family.get('correction') in {'holm', 'bh'}:
            if type(family.get('alpha')) not in (int, float) or not 0 < family['alpha'] < 1:
                problems.append(issue('UNRESOLVED', 'MULTIPLICITY_ALPHA_UNRESOLVED', fid))
            supported = {'mann_whitney': ('equal_distributions', 'unpaired'),
                'paired_signed_rank': ('zero_symmetric_pair_difference', 'paired')}
            for hypothesis in [h for h in plan['hypotheses'] if h.get('family_id') == fid]:
                hid = hypothesis.get('hypothesis_id')
                declaration = supported.get(hypothesis.get('test'))
                endpoint = endpoints.get(hypothesis.get('endpoint_id'), {})
                if (declaration is None or hypothesis.get('null') != declaration[0] or
                    endpoint.get('pairing', {}).get('mode') != declaration[1] or
                    hypothesis.get('sidedness') not in {'two_sided', 'greater', 'less'} or
                    not hypothesis.get('sampling_assumptions') or hypothesis['sampling_assumptions'] == 'UNRESOLVED'):
                    problems.append(issue('UNRESOLVED', 'INFERENTIAL_HYPOTHESIS_CONTRACT_UNRESOLVED', hid))

    if any(h.get('family_id') not in family_ids for h in plan['hypotheses']) or not family_ids:
        problems.append(issue('MISSING', 'HYPOTHESIS_FAMILY_MISSING'))
    return problems


def verify_measurement_sources(measurements, evidence_root):
    """Bind metric receipt values to immutable original source rows; no estimates."""
    problems, files = [], {}
    for m in measurements:
        key = m.get('measurement_id', '<missing>')
        required = {'source_path', 'source_file_sha256', 'source_row_index', 'source_row_sha256', 'value_column', 'value', 'source_columns', 'extraction_source_path'} | SOURCE_ID_FIELDS
        if not required <= m.keys():
            problems.append(issue('MISSING', 'MEASUREMENT_SOURCE_FIELDS_MISSING', key)); continue
        try:
            path = core.rooted_path(evidence_root, m['source_path'])
            if path not in files:
                files[path] = (core.file_sha256(path), core.read_csv_records(path))
            sha, rows = files[path]
            ordinal = m['source_row_index']
            if sha != m['source_file_sha256'] or type(ordinal) is not int or not 1 <= ordinal <= len(rows):
                raise ValueError('Source hash/ordinal mismatch')
            raw = rows[ordinal-1]
            if policy_digest(raw) != m['source_row_sha256']:
                raise ValueError('Source row hash mismatch')
            columns = m['source_columns']
            if set(columns) != SOURCE_ID_FIELDS or len(set(columns.values())) != len(columns):
                raise ValueError('Missing/ambiguous source identity columns')
            if any(raw[columns[k]] != str(m[k]) for k in SOURCE_ID_FIELDS):
                raise ValueError('Source row identity mismatch')
            extraction = core.rooted_path(evidence_root, m['extraction_source_path'])
            if core.file_sha256(extraction) != m['extraction_code_id']:
                raise ValueError('Extraction code identity mismatch')
            value = raw[m['value_column']]
            if m['value'] is None:
                if value != '' or not m.get('missing_reason'):
                    raise ValueError('Null value lacks matching missing source/reason')
            elif type(m['value']) not in (int, float) or not math.isfinite(m['value']) or float(value) != m['value']:
                raise ValueError('Nonfinite or mismatched source value')
        except FileNotFoundError as exc:
            problems.append(issue('MISSING', 'MEASUREMENT_SOURCE_MISSING', key, exc))
        except (OSError, KeyError, TypeError, ValueError, OverflowError) as exc:
            problems.append(issue('CONFLICTING', 'MEASUREMENT_SOURCE_MISMATCH', key, exc))
    return problems


def assess_validated_context(contexts, plan, approval, approved_plan_sha256,
                             measurements, evidence_root):
    """Internal policy seam for already P1b-validated contexts; fixtures may inject these."""
    global_issues, designs, predictions, configurations = [], {}, {}, {}
    design_variants, prediction_variants = defaultdict(list), defaultdict(list)
    ambiguous_designs, ambiguous_predictions, runs = set(), set(), {}
    generations, members, prediction_coverage = [], [], []
    for context in contexts:
        bundle = context['validated']
        generations.append(bundle.get('generation_code_status'))
        prediction_coverage.append({
            'design_ids': sorted({r['design_id'] for r in bundle['records']}),
            'run_contract_sha256': policy_digest(context['contract']) if context.get('contract') else None,
            'prediction_presence_counts': dict(sorted(Counter(r['prediction_presence_status'] for r in bundle['eligibility']).items())),
            'sequence_validation_counts': dict(sorted(Counter(r['sequence_validation_status'] for r in bundle['eligibility']).items())),
            'run_provenance_counts': dict(sorted(Counter(r['run_provenance_status'] for r in bundle['eligibility']).items()))})
        if bundle.get('generation_code_status') != 'VERIFIED' or not bundle.get('generation_code_commit'):
            global_issues.append(issue('UNRESOLVED', 'GENERATION_SOURCE_RUN_NOT_ESTABLISHED'))
        if bundle.get('mapping_authority_status') != 'VERIFIED':
            global_issues.append(issue('UNRESOLVED', 'SCIENTIFIC_MAPPING_AUTHORITY_UNRESOLVED'))
        for row in bundle['records']:
            stable = {k: row[k] for k in ('design_id', 'record_key', 'target', 'method', 'source_sequence_sha256')}
            design_variants[row['design_id']].append(stable)
            existing = designs.get(row['design_id'])
            if existing is not None and existing != stable:
                global_issues.append(issue('CONFLICTING', 'CONFLICTING_DESIGN_IDENTITY', row['design_id']))
                ambiguous_designs.add(row['design_id'])
            if existing is None:
                designs[row['design_id']] = stable
        contract = context.get('contract')
        if contract:
            instance, digest = contract.get('run_instance_id'), policy_digest(contract)
            if instance in runs and runs[instance] != digest:
                global_issues.append(issue('CONFLICTING', 'CONFLICTING_RUN_INSTANCE_DEFINITION', instance))
            runs[instance] = digest
        for observed in bundle['eligibility']:
            key = (observed['design_id'], observed['prediction_id'])
            if not observed['prediction_id']:
                continue
            prediction_variants[key].append(observed)
            if key in predictions:
                ambiguous_predictions.add(key)
                global_issues.append(issue('CONFLICTING', 'DUPLICATE_PREDICTION_CONTEXT', observed['prediction_id']))
            if key not in predictions:
                predictions[key] = observed
                configurations[key] = metric_configuration(contract) if contract else None
    for key in ambiguous_designs:
        designs.pop(key, None)
    for key in ambiguous_predictions:
        predictions.pop(key, None); configurations.pop(key, None)
    if {r['target'] for r in designs.values()} != CASES or any(
            {r['method'] for r in designs.values() if r['target'] == case} != METHODS for case in CASES):
        global_issues.append(issue('MISSING', 'CASE_METHOD_INPUT_COVERAGE_MISSING'))
    if plan is None:
        global_issues.append(issue('MISSING', 'ANALYSIS_PLAN_MISSING'))
        endpoints, hypotheses, families, planned = [], [], [], []
    else:
        try:
            global_issues.extend(validate_plan(plan, designs, evidence_root))
            global_issues.extend(approval_issues(plan, approval, approved_plan_sha256))
            endpoints, hypotheses, families = plan.get('endpoints', []), plan.get('hypotheses', []), plan.get('families', [])
            planned = plan.get('evaluation_members', [])
        except (KeyError, TypeError, ValueError) as exc:
            global_issues.append(issue('CONFLICTING', 'MALFORMED_ANALYSIS_PLAN', detail=exc))
            endpoints, hypotheses, families, planned = [], [], [], []
    metric_problems = verify_measurement_sources(measurements, evidence_root)
    measurement_issues = defaultdict(list)
    for problem in metric_problems:
        measurement_issues[problem['key']].append(problem)
    by_measurement, measurement_keys, ambiguous_measurements = {}, {}, set()
    measurement_variants, key_members = defaultdict(list), defaultdict(list)
    for m in measurements:
        mid = m.get('measurement_id', '<missing>')
        measurement_variants[mid].append(m)
        if mid in by_measurement:
            measurement_issues[mid].append(issue('CONFLICTING', 'DUPLICATE_OR_CONFLICTING_MEASUREMENT', mid))
            ambiguous_measurements.add(mid)
        if mid not in by_measurement:
            by_measurement[mid] = m
        key = (m.get('design_id'), m.get('prediction_id'), m.get('endpoint_id'))
        key_members[key].append(mid)
        if key in measurement_keys:
            for duplicate in (mid, measurement_keys[key]):
                measurement_issues[duplicate].append(issue('CONFLICTING', 'AMBIGUOUS_MEASUREMENT_JOIN', duplicate))
                ambiguous_measurements.add(duplicate)
        if key not in measurement_keys:
            measurement_keys[key] = mid
        identity = predictions.get((m.get('design_id'), m.get('prediction_id')))
        if identity is None:
            measurement_issues[mid].append(issue('MISSING', 'PREDICTION_IDENTITY_MISSING', mid))
        else:
            try: core.validate_measurement_joins([m], [identity])
            except (KeyError, TypeError, ValueError) as exc:
                measurement_issues[mid].append(issue('CONFLICTING', 'P1B_MEASUREMENT_JOIN_REJECTED', mid, exc))
        if m.get('evidence_class') == 'historical':
            measurement_issues[mid].append(issue('HISTORICAL', 'HISTORICAL_MEASUREMENT_NOT_PROMOTED', mid))
        elif m.get('evidence_class') != 'current':
            measurement_issues[mid].append(issue('UNRESOLVED', 'MEASUREMENT_EVIDENCE_CLASS_UNRESOLVED', mid))
    for mid in ambiguous_measurements:
        by_measurement.pop(mid, None)
    measurement_keys = {key: ids[0] for key, ids in key_members.items()
                        if len(ids) == 1 and ids[0] not in ambiguous_measurements}
    global_issues.extend(p for problems in measurement_issues.values() for p in problems)
    # Duplicate endpoints are quarantined; no last-row definition is selected.
    counts = Counter(e.get('endpoint_id') for e in endpoints)
    endpoints = [e for e in endpoints if counts[e.get('endpoint_id')] == 1]
    endpoint_index = {e.get('endpoint_id'): e for e in endpoints}
    slots, slots_seen, used_measurements = defaultdict(list), set(), set()
    for member in sorted(planned, key=diagnostic_key):
        key = (member.get('endpoint_id'), member.get('design_id'))
        slot_id = key + (member.get('replicate_id'),)
        if slot_id in slots_seen:
            global_issues.append(issue('CONFLICTING', 'DUPLICATE_PLANNED_REPLICATE_SLOT', '/'.join(map(str, slot_id))))
        slots_seen.add(slot_id)
        if key[0] not in endpoint_index or key[1] not in designs:
            global_issues.append(issue('CONFLICTING', 'UNMATCHED_PLANNED_EVALUATION_MEMBER', str(key)))
        slots[key].append(member)
    cells = []
    for endpoint in sorted(endpoints, key=lambda e: str(e.get('endpoint_id'))):
        eid = endpoint.get('endpoint_id')
        for did, design in sorted(designs.items()):
            if design['target'] != endpoint.get('case'):
                continue
            problems, local_members = [], []
            expected_slots = slots.get((eid, did), [])
            if endpoint.get('kind') == 'sequence_score':
                problems.append(issue('UNRESOLVED', 'SEQUENCE_SCORE_PROVENANCE_UNSUPPORTED', eid))
            if not expected_slots:
                problems.append(issue('MISSING', 'PLANNED_REPLICATE_MISSING', did))
            replication = endpoint.get('replication', {})
            if replication.get('rule') == 'single_per_draw' and len(expected_slots) > 1:
                problems.append(issue('CONFLICTING', 'SINGLE_PREDICTION_CARDINALITY_VIOLATION', did))
            for member in sorted(expected_slots, key=lambda s: str(s.get('replicate_id'))):
                owner, prediction_id = member.get('prediction_owner_design_id'), member.get('prediction_id')
                observed = predictions.get((owner, prediction_id))
                if observed is None:
                    problems.append(issue('MISSING', 'PLANNED_PREDICTION_MISSING', prediction_id)); continue
                if str(member.get('replicate_id')) != observed.get('prediction_replicate_id'):
                    problems.append(issue('CONFLICTING', 'PREDICTION_REPLICATE_IDENTITY_MISMATCH', prediction_id))
                if owner != did and (replication.get('rule') != 'shared_deterministic' or
                        designs.get(owner, {}).get('source_sequence_sha256') != design['source_sequence_sha256']):
                    problems.append(issue('CONFLICTING', 'REUSE_MISLABELED_AS_INDEPENDENT_FOLD', prediction_id))
                if member.get('configuration_id') != configurations.get((owner, prediction_id)) or member.get('configuration_id') != replication.get('configuration_id'):
                    problems.append(issue('CONFLICTING', 'REPLICATE_CONFIGURATION_MISMATCH', prediction_id))
                mid = measurement_keys.get((owner, prediction_id, eid))
                if mid is None:
                    problems.append(issue('MISSING', 'MEASUREMENT_MISSING', prediction_id)); continue
                used_measurements.add(mid)
                m = by_measurement[mid]
                problems.extend(measurement_issues[mid])
                if mid in ambiguous_measurements:
                    continue
                if any(m.get(k) != endpoint.get(k) for k in ('metric', 'region', 'unit', 'metric_policy_id', 'extraction_code_id', 'reference_sha256', 'residue_map_sha256')):
                    problems.append(issue('CONFLICTING', 'METRIC_DEFINITION_OR_UNIT_MISMATCH', mid))
                if m.get('value') is None:
                    problems.append(issue('MISSING', 'MEASUREMENT_NULL', mid, m.get('missing_reason', 'UNSPECIFIED')))
                else:
                    value = m['value']; bounds = endpoint.get('bounds', {})
                    if type(value) not in (int, float) or not math.isfinite(value):
                        problems.append(issue('INELIGIBLE', 'NONFINITE_MEASUREMENT', mid))
                    elif endpoint.get('unit') in UNIT_DOMAINS and (
                            value < UNIT_DOMAINS[endpoint['unit']][0] or
                            (UNIT_DOMAINS[endpoint['unit']][1] is not None and value > UNIT_DOMAINS[endpoint['unit']][1])):
                        problems.append(issue('INELIGIBLE', 'MEASUREMENT_OUTSIDE_UNIT_DOMAIN', mid))
                    elif ((bounds.get('lower') is not None and value < bounds['lower']) or
                          (bounds.get('upper') is not None and value > bounds['upper'])):
                        problems.append(issue('INELIGIBLE', 'METRIC_RANGE_VIOLATION', mid))
                local_members.append({'measurement_id': mid, 'prediction_id': prediction_id,
                    'prediction_owner_design_id': owner, 'replicate_id': str(member.get('replicate_id')),
                    'configuration_id': member.get('configuration_id')})
            technical = normalized_issues(problems)
            scientific = normalized_issues(problems + global_issues)
            cells.append({'endpoint_id': eid, 'design_id': did, 'record_key': design['record_key'],
                'case': design['target'], 'method': design['method'], 'source_sequence_sha256': design['source_sequence_sha256'],
                'planned_replicates': len(expected_slots), 'retained_members': local_members,
                'technical_state': evidence_state(technical), 'technical_reasons': technical,
                'scientific_state': evidence_state(scientific), 'scientific_reasons': scientific})
    for mid in set(by_measurement) - used_measurements:
        global_issues.append(issue('INELIGIBLE', 'UNPLANNED_MEASUREMENT_RETAINED', mid))
    # Recompute scientific states after all global conflicts/exclusions are known.
    for cell in cells:
        cell['scientific_reasons'] = normalized_issues(cell['technical_reasons'] + global_issues)
        cell['scientific_state'] = evidence_state(cell['scientific_reasons'])
    family_report = []
    for family in sorted(families, key=diagnostic_key):
        fid = family.get('family_id'); hs = [h for h in hypotheses if h.get('family_id') == fid]
        missing = []
        for h in hs:
            relevant = [c for c in cells if c['endpoint_id'] == h.get('endpoint_id')]
            if not relevant or any(c['scientific_state'] != 'VERIFIED' for c in relevant):
                missing.append(h.get('hypothesis_id'))
        family_report.append({'family_id': fid, 'planned_hypothesis_ids': sorted(family.get('hypothesis_ids', [])),
            'unavailable_hypothesis_ids': sorted(missing), 'state': 'MISSING' if missing else 'VERIFIED',
            'reason_codes': ['PLANNED_FAMILY_INCOMPLETE_NO_SHRINK'] if missing else [],
            'correction': family.get('correction'), 'adjusted_values_computed': False})
    all_problems = normalized_issues(global_issues + [p for c in cells for p in c['technical_reasons']] +
        [issue(f['state'], code, f['family_id']) for f in family_report for code in f['reason_codes']])
    support = defaultdict(list)
    support_predictions = defaultdict(set)
    for cell in cells:
        skey = (cell['endpoint_id'], cell['method'], cell['source_sequence_sha256'])
        support[skey].append(cell['design_id'])
        support_predictions[skey].update(m['prediction_id'] for m in cell['retained_members'])
    return {'schema_version': 'uirj-eligibility-report-v1', 'assurance': 'CONTRACT_ONLY_NOT_SCIENTIFIC_VALIDITY',
        'approval_assurance': 'OPERATOR_ATTESTATION_BINDING_NOT_AUTHENTICATED_AUTHORITY',
        'analysis_plan_sha256': policy_digest(plan) if plan is not None else None,
        'approval_binding': {'operator_plan_sha256': approved_plan_sha256, 'receipt': approval},
        'state': evidence_state(all_problems), 'reasons': all_problems,
        'generation_evidence_states': sorted(set(generations), key=str),
        'prediction_coverage': sorted(prediction_coverage, key=diagnostic_key),
        'designs': [designs[k] for k in sorted(designs)], 'design_count': len(design_variants), 'unambiguous_design_count': len(designs),
        'quarantined_design_variants': [{'design_id': k, 'variants': sorted(design_variants[k], key=diagnostic_key)} for k in sorted(ambiguous_designs)],
        'quarantined_prediction_variants': [{'design_id':k[0], 'prediction_id':k[1], 'variants':sorted(prediction_variants[k], key=diagnostic_key)} for k in sorted(ambiguous_predictions)],
        'quarantined_measurement_ids': sorted(ambiguous_measurements),
        'unique_sequence_count': len({r['source_sequence_sha256'] for r in designs.values()}),
        'prediction_identity_count': len(prediction_variants),
        'unambiguous_prediction_identity_count': len(predictions), 'measurement_identity_count': len(measurement_variants),
        'cells': cells, 'families': family_report,
        'estimand_memberships': {'draw_weighted': [{'endpoint_id':c['endpoint_id'], 'design_id':c['design_id'], 'member_prediction_ids':[m['prediction_id'] for m in c['retained_members']]} for c in cells],
            'unique_sequence_support': [{'endpoint_id':e,'method':m,'source_sequence_sha256':s,'design_ids':sorted(ids), 'prediction_ids':sorted(support_predictions[(e,m,s)]),
                'support_weight_policy':'equal_per_sequence_within_case_method',
                'within_sequence_policy':'mean_over_distinct_prediction_ids'} for (e,m,s),ids in sorted(support.items())]},
        'statistics_computed': False, 'scientific_validity_established': False}


def evaluate_eligibility(manifests, input_root, upstream_root, output_root,
                         contracts=None, plan=None, approval=None,
                         approved_plan_sha256=None, measurements=None, evidence_root=None):
    """Public entry: obtain observed P1b evidence rather than trusting saved PASS columns."""
    contracts = contracts if contracts is not None else [None] * len(manifests)
    if len(contracts) != len(manifests) or not manifests:
        raise ValueError('One optional run contract per input manifest required')
    if measurements is not None and (not isinstance(measurements, list) or any(not isinstance(m, dict) for m in measurements)):
        raise ValueError('Measurements must be a list of receipt objects')
    contexts = []
    for manifest, contract in zip(manifests, contracts):
        validated = core.validate_prediction_bundle(manifest, input_root, upstream_root, output_root, contract)
        contexts.append({'validated': validated, 'contract': contract})
    return assess_validated_context(contexts, plan, approval, approved_plan_sha256,
                                    measurements or [], evidence_root)
