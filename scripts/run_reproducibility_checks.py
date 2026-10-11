"""Explicit legacy arithmetic or evidence eligibility. Import has no side effects."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = (
    'scripts/reproduce_contact_statistics.py',
    'scripts/reproduce_sequence_conservation.py',
    'scripts/reproduce_plddt_statistics.py',
    'scripts/reproduce_rmsd_statistics.py',
)


def run_legacy_arithmetic():
    """Explicit compatibility mode, never called from eligibility mode."""
    print('HISTORICAL ARITHMETIC ONLY: INPUT LINEAGE AND SCIENTIFIC VALIDITY UNVERIFIED')
    failures = []
    for script in SCRIPTS:
        result = subprocess.run([sys.executable, str(ROOT / script)], cwd=ROOT,
                                capture_output=True, text=True)
        print(f'Running: {script}\n{result.stdout}')
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        if result.returncode:
            failures.append(script)
    if failures:
        print('LEGACY CHECK FAILED (no dependency installation attempted): ' + ', '.join(failures))
        return 1
    print('SUCCESS: ALL DOWNSTREAM REPRODUCIBILITY CHECKS PASSED (HISTORICAL ARITHMETIC ONLY)')
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', required=True, choices=['eligibility', 'legacy-arithmetic'])
    parser.add_argument('--manifest', action='append', default=[])
    parser.add_argument('--run-contract', action='append', default=[])
    parser.add_argument('--input-root')
    parser.add_argument('--upstream-root')
    parser.add_argument('--output-root')
    parser.add_argument('--evidence-root')
    parser.add_argument('--plan')
    parser.add_argument('--approval-receipt')
    parser.add_argument('--approved-plan-sha256')
    parser.add_argument('--measurement-receipts')
    parser.add_argument('--report')
    args = parser.parse_args(argv)
    if args.mode == 'legacy-arithmetic':
        if any([args.manifest, args.run_contract, args.plan, args.approval_receipt,
                args.approved_plan_sha256, args.measurement_receipts, args.report,
                args.input_root, args.upstream_root, args.output_root, args.evidence_root]):
            parser.error('Eligibility arguments cannot be combined with legacy arithmetic')
        return run_legacy_arithmetic()
    # Lazy import preserves a simple side-effect-free legacy entry point.
    import statistics_contract as contract
    if not args.manifest or not all((args.input_root, args.upstream_root, args.output_root)):
        parser.error('Eligibility requires manifest and explicit input/upstream/output roots')
    try:
        runs = [contract.read_json(p) for p in args.run_contract] if args.run_contract else None
        plan = contract.read_json(args.plan) if args.plan else None
        approval = contract.read_json(args.approval_receipt) if args.approval_receipt else None
        measurements = []
        if args.measurement_receipts:
            receipt = contract.read_json(args.measurement_receipts)
            if set(receipt) != {'schema_version', 'measurements'} or receipt['schema_version'] != contract.RECEIPT_SCHEMA or not isinstance(receipt['measurements'], list):
                raise ValueError('Unsupported metric receipt envelope')
            measurements = receipt['measurements']
        report = contract.evaluate_eligibility(args.manifest, args.input_root, args.upstream_root,
            args.output_root, runs, plan, approval, args.approved_plan_sha256,
            measurements, args.evidence_root)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        report = {'schema_version': 'uirj-eligibility-report-v1', 'state': 'CONFLICTING',
                  'assurance': 'CONTRACT_ONLY_NOT_SCIENTIFIC_VALIDITY',
                  'reasons': [contract.issue('CONFLICTING', 'EVIDENCE_VALIDATION_REJECTED', detail=exc)],
                  'statistics_computed': False, 'scientific_validity_established': False}
    text = json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + '\n'
    if args.report:
        # New task-local file only, no replacement of saved measurements/reports.
        with open(args.report, 'x') as handle:
            handle.write(text)
    print(text, end='')
    return 0 if report['state'] == 'VERIFIED' else 2


if __name__ == '__main__':
    raise SystemExit(main())
