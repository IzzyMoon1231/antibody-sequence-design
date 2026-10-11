
import argparse
import csv
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

import torch

from validate_esmfold_predictions import (
    extract_pdb_sequence,
    validate_sequence,
)

def load_manifest(path, **kwargs):
    from validate_esmfold_predictions import load_input_manifest
    # Existing API retains legacy inspection; v1 always uses strict context.
    with open(path, newline='') as handle:
        rows = list(csv.DictReader(handle))
    strict = bool(rows and rows[0].get('schema_version'))
    return load_input_manifest(path, legacy=not strict, **kwargs)


def check_prediction(path, sequence):
    from validate_esmfold_predictions import check_prediction_status
    return check_prediction_status(path, sequence)


def observed_runtime_versions():
    import platform
    import importlib.metadata
    return {'python': platform.python_version(), 'torch': importlib.metadata.version('torch'),
            'transformers': importlib.metadata.version('transformers')}


def verify_model_parameters(model):
    parameters = list(model.named_parameters())
    if not parameters:
        raise ValueError('Cannot verify an empty model')
    for name, parameter in parameters:
        expected = torch.float16 if name.startswith('esm.') else torch.float32
        if parameter.is_floating_point() and parameter.dtype != expected:
            raise ValueError(f'Model precision mismatch: {name}')
        if parameter.device.type != 'cuda':
            raise ValueError(f'Model device mismatch: {name}')


def verify_execution_contract(contract, rows, asset_root):
    """Verify locally observable settings/assets; does not load model or authorize work."""
    from validate_esmfold_predictions import validate_run_contract, rooted_path, file_sha256
    run_id = validate_run_contract(contract, rows)
    if any(r['planned_run_id'] != run_id for r in rows):
        raise ValueError('Manifest/run contract identity mismatch')
    actual = observed_runtime_versions()
    if actual != contract['runtime'] or file_sha256(__file__) != contract['serializer_sha256']:
        raise ValueError('Actual runtime/serializer differs from run contract')
    root = Path(asset_root).resolve()
    if not root.is_dir():
        raise ValueError('Explicit local model assets required')
    observed = {}
    for path in root.rglob('*'):
        if path.is_symlink():
            raise ValueError('Symlinked model assets are unsupported')
        if path.is_file():
            relative = path.relative_to(root).as_posix()
            observed[relative] = file_sha256(rooted_path(root, relative))
    if observed != contract['assets']:
        raise ValueError('Model asset inventory/hash mismatch')
    return run_id


def load_model(contract=None, asset_root=None):
    if contract is None or asset_root is None:
        raise ValueError('Pinned local run contract required before model loading')
    if not torch.cuda.is_available():
        raise RuntimeError('CUDA GPU unavailable. Cannot run ESMFold predictions.')
    from transformers import EsmForProteinFolding, AutoTokenizer
    # Approved assets were verified before this call; never use remote fallback.
    tokenizer = AutoTokenizer.from_pretrained(str(asset_root), local_files_only=True,
                                               revision=contract['model_revision'])
    model = EsmForProteinFolding.from_pretrained(
        str(asset_root), local_files_only=True, revision=contract['model_revision'],
        low_cpu_mem_usage=True)
    model = model.eval().cuda()
    model.esm = model.esm.half()
    model.set_chunk_size(contract['inference_settings']['chunk_size'])
    verify_model_parameters(model)
    return model, tokenizer

def serialize_esmfold_output(output, batch_index=0):
    """Export one batch member from ESMFold atom14 output, without model I/O.

    Positions have shape (structure stages, batch, residues, 14, 3).
    Mapping and masks follow the Hugging Face/OpenFold atom37 contract.
    Missing atoms and fully masked padding are allowed; mislabeled maps are not.
    """
    from transformers.models.esm.openfold_utils.data_transforms import make_atom14_masks
    from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
    from transformers.models.esm.openfold_utils.protein import Protein as OFProtein, to_pdb

    fields = (
        "positions", "aatype", "atom37_atom_exists",
        "residx_atom37_to_atom14", "residue_index", "plddt",
    )
    missing = [name for name in fields if name not in output]
    if missing:
        raise ValueError(f"Missing ESMFold output fields: {missing}")
    if any(not isinstance(output[name], torch.Tensor) for name in fields):
        raise ValueError("ESMFold output fields must be tensors")
    # Transfer only the final structure stage and fields needed for serialization.
    positions = output["positions"]
    if positions.ndim != 5 or positions.shape[0] < 1 or positions.shape[-2:] != (14, 3):
        raise ValueError("positions shape must be (stages, batch, residues, 14, 3)")
    batch, residues = positions.shape[1:3]
    if batch < 1 or residues < 1:
        raise ValueError("positions shape requires nonempty batch and residues")
    if not isinstance(batch_index, int) or not 0 <= batch_index < batch:
        raise ValueError("batch_index is outside the output batch")
    expected_shapes = {
        "aatype": (batch, residues),
        "atom37_atom_exists": (batch, residues, 37),
        "residx_atom37_to_atom14": (batch, residues, 37),
        "residue_index": (batch, residues),
        "plddt": (batch, residues, 37),
    }
    for name, shape in expected_shapes.items():
        if tuple(output[name].shape) != shape:
            raise ValueError(f"{name} shape must be {shape}")
    data = {name: output[name].detach().cpu() for name in expected_shapes}
    final_positions = positions[-1].detach().float().cpu()
    for name, value in {**data, "positions": final_positions}.items():
        if not torch.isfinite(value).all():
            raise ValueError(f"{name} must contain finite values")
    for name in ("aatype", "residue_index", "residx_atom37_to_atom14"):
        if data[name].dtype not in (torch.int32, torch.int64):
            raise ValueError(f"{name} must have integer dtype")
    if not ((data["aatype"] >= 0) & (data["aatype"] <= 20)).all():
        raise ValueError("aatype must be an OpenFold residue index in [0, 20]")
    canonical = make_atom14_masks({"aatype": data["aatype"]})
    if not torch.equal(data["residx_atom37_to_atom14"], canonical["residx_atom37_to_atom14"]):
        raise ValueError("atom37-to-atom14 mapping does not match aatype")
    mask = data["atom37_atom_exists"]
    if not ((mask == 0) | (mask == 1)).all() or (mask > canonical["atom37_atom_exists"]).any():
        raise ValueError("atom37_atom_exists must be a binary subset of the residue atom mask")

    atom37_positions = atom14_to_atom37(final_positions, data).numpy()
    i = batch_index
    protein = OFProtein(
        aatype=data["aatype"][i].numpy(),
        atom_positions=atom37_positions[i],
        atom_mask=mask[i].numpy(),
        residue_index=data["residue_index"][i].numpy() + 1,
        b_factors=data["plddt"][i].float().numpy(),
        chain_index=None,
    )
    return to_pdb(protein)


def predict(sequence, model, tokenizer, inference_settings=None):
    inputs = tokenizer(
        [sequence],
        return_tensors="pt",
        add_special_tokens=False,
    )

    inputs = {
        key: value.cuda()
        for key, value in inputs.items()
    }

    with torch.no_grad():
        kwargs = {} if inference_settings is None else {'num_recycles': inference_settings['num_recycles']}
        output = model(**inputs, **kwargs)

    return serialize_esmfold_output(output)


def save_prediction(pdb_text, destination, sequence, provenance=None):
    import os
    import tempfile
    from pathlib import Path

    destination = Path(destination)
    if provenance is not None:
        from validate_esmfold_predictions import validate_run_contract, file_sha256
        row, contract, rows = provenance['row'], provenance['contract'], provenance['rows']
        if row['sequence'] != sequence or row['planned_run_id'] != validate_run_contract(contract, rows):
            raise ValueError('Persistence input/run identity mismatch')
        if observed_runtime_versions() != contract['runtime'] or file_sha256(__file__) != contract['serializer_sha256']:
            raise ValueError('Persistence actual runtime/serializer mismatch')
    destination.parent.mkdir(parents=True, exist_ok=True)

    lock_path = destination.with_name(destination.name + ".lock")
    lock_fd = None
    temp_path = None

    try:
        # Exclusive lock prevents concurrent writers.
        lock_fd = os.open(
            lock_path,
            os.O_CREAT | os.O_EXCL | os.O_WRONLY,
            0o600,
        )

        if destination.exists():
            raise FileExistsError(destination)
        receipt_path = destination.with_name(destination.name + '.provenance.json')
        if receipt_path.exists():
            raise FileExistsError(receipt_path)

        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".pdb",
            prefix=".esmfold_tmp_",
            dir=destination.parent,
            delete=False,
        ) as temp:
            temp.write(pdb_text)
            temp_path = Path(temp.name)

        observed = extract_pdb_sequence(temp_path)
        status, _ = validate_sequence(sequence, observed)

        if status != "PASS":
            raise ValueError(
                f"Generated PDB failed validation: {status}. "
                f"Expected {len(sequence)} residues, "
                f"observed {len(observed)}."
            )

        # The lock serializes cooperating prediction runners.
        # Never replace a destination that already exists.
        if destination.exists():
            raise FileExistsError(destination)

        # Google Drive supports exclusive file creation, but not hard links.
        # Never overwrite an existing prediction.
        # A completion marker distinguishes finished files from interrupted copies.
        import shutil

        marker = destination.with_name(destination.name + ".complete")

        if marker.exists():
            raise FileExistsError(
                f"Completion marker already exists: {marker}"
            )

        fd = os.open(
            destination,
            os.O_CREAT | os.O_EXCL | os.O_WRONLY,
            0o600,
        )

        try:
            with os.fdopen(fd, "wb") as output:
                with open(temp_path, "rb") as source_file:
                    shutil.copyfileobj(source_file, output)
                output.flush()
                os.fsync(output.fileno())

            copied_sequence = extract_pdb_sequence(destination)
            copied_status, _ = validate_sequence(
                sequence, copied_sequence
            )

            if copied_status != "PASS":
                raise ValueError(
                    f"Copied PDB failed validation: {copied_status}"
                )

            if provenance is not None:
                import json
                from datetime import datetime, timezone
                from validate_esmfold_predictions import (file_sha256, input_bundle_digest,
                    validate_receipt, parse_pdb_residues)
                residues = parse_pdb_residues(destination)
                if (len({r['chain'] for r in residues}) != 1 or
                    any(not r['ca_present'] for r in residues) or
                    [(r['residue_number'], r['insertion_code']) for r in residues] !=
                        [(i, ' ') for i in range(1, len(sequence)+1)]):
                    raise ValueError('Ambiguous generated prediction mapping')
                receipt = dict(schema_version='uirj-prediction-receipt-v1',
                    record_key=row['record_key'], design_id=row['design_id'],
                    prediction_id=row['prediction_id'], prediction_replicate_id=row['prediction_replicate_id'],
                    source_sequence_sha256=row['source_sequence_sha256'],
                    input_bundle_digest=input_bundle_digest(rows),
                    run_contract_sha256=row['planned_run_id'], output_pdb_sha256=file_sha256(destination),
                    serializer_sha256=contract['serializer_sha256'], runtime=contract['runtime'],
                    completed_at=datetime.now(timezone.utc).isoformat())
                validate_receipt(receipt, row, contract, rows, destination)
                receipt_fd = os.open(receipt_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
                with os.fdopen(receipt_fd, 'w') as handle:
                    json.dump(receipt, handle, sort_keys=True)
                    handle.write('\n'); handle.flush(); os.fsync(handle.fileno())

            # Publish completion only after the copied PDB and optional receipt validate.
            marker_fd = os.open(
                marker,
                os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                0o600,
            )

            with os.fdopen(marker_fd, "w") as marker_file:
                marker_file.write("COMPLETE\n")
                marker_file.flush()
                os.fsync(marker_file.fileno())

        except Exception:
            # Leave any incomplete destination for manual inspection.
            # Never delete a potentially valuable prediction automatically.
            raise

        temp_path.unlink()
        temp_path = None

    finally:
        if temp_path is not None and temp_path.exists():
            temp_path.unlink()

        if lock_fd is not None:
            os.close(lock_fd)
            lock_path.unlink(missing_ok=True)


def inspect_prediction_locks(manifest_path, minimum_age_minutes=60, output_root=None):
    """Report prediction locks without deleting or changing them."""
    import csv
    import time
    from pathlib import Path

    if minimum_age_minutes < 0:
        raise ValueError("minimum_age_minutes must be nonnegative")

    with open(manifest_path, newline="") as handle:
        rows = list(csv.DictReader(handle))

    now = time.time()
    findings = []

    for row in rows:
        if row.get('schema_version'):
            from validate_esmfold_predictions import rooted_path
            destination = rooted_path(output_root, row['pdb_path'])
        else:
            destination = Path(row['pdb_path'])
        lock = destination.with_name(destination.name + ".lock")

        if not lock.exists():
            continue

        age_minutes = (now - lock.stat().st_mtime) / 60
        prediction_status = check_prediction(
            destination, row["sequence"]
        )

        findings.append({
            "lock_path": str(lock),
            "age_minutes": round(age_minutes, 2),
            "prediction_status": prediction_status,
            "review_recommended": age_minutes >= minimum_age_minutes,
        })

    return findings

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--manifest", required=True)
    parser.add_argument("--target", choices=["9NFU", "9NH7"])
    parser.add_argument(
        "--method", choices=["baseline", "ensemble"]
    )
    parser.add_argument("--limit", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument('--input-root')
    parser.add_argument('--upstream-root')
    parser.add_argument('--output-root')
    parser.add_argument('--run-contract')
    parser.add_argument('--asset-root')
    parser.add_argument('--execute-approved-run', help='Explicit approved contract SHA256; not inferred from JSON')

    args = parser.parse_args()

    if args.limit < 1:
        parser.error("--limit must be at least 1")

    rows = load_manifest(args.manifest, input_root=args.input_root,
                         upstream_root=args.upstream_root, output_root=args.output_root)
    strict = bool(rows[0].get('schema_version'))
    contract = None
    if args.run_contract:
        import json
        contract = json.loads(Path(args.run_contract).read_text())
    if strict:
        from validate_esmfold_predictions import _prediction_eligibility, rooted_path
        eligibility = {r['record_key']: _prediction_eligibility(r, args.output_root, rows, contract) for r in rows}

    selected = [
        row for row in rows
        if (args.target is None or row["target"] == args.target)
        and (args.method is None or row["method"] == args.method)
    ]

    statuses = Counter()
    pending = []

    for row in selected:
        if strict:
            observed = eligibility[row['record_key']]
            status = observed['prediction_presence_status']
            if status == 'PASS' and observed['run_provenance_status'] != 'VERIFIED':
                status = 'SEQUENCE_MATCH_RUN_UNVERIFIED'
        else:
            status = check_prediction(row["pdb_path"], row["sequence"])

        statuses[status] += 1

        if status == "MISSING":
            pending.append(row)

    print("ESMFOLD PREDICTION RUNNER")
    print("=" * 55)
    print("Selected sequences:", len(selected))
    print("Existing prediction statuses:", dict(statuses))
    print("Missing predictions:", len(pending))
    print("Prediction limit:", args.limit)

    if args.dry_run:
        print("\nDRY RUN: No predictions generated.")

        for row in pending[:args.limit]:
            print(
                row["target"],
                row["method"],
                row["index"],
                "->",
                row["pdb_path"],
            )

        return

    invalid_existing = sum(
        count for status, count in statuses.items()
        if status not in {"PASS", "MISSING"}
    )

    if invalid_existing:
        raise RuntimeError(
            f"{invalid_existing} existing PDB files failed "
            "validation. Investigate before continuing."
        )

    if not pending:
        print("All selected predictions already completed.")
        return

    if not strict or contract is None or not args.asset_root:
        raise ValueError('Execution requires v1 manifest, complete local run contract and asset root')
    run_id = verify_execution_contract(contract, rows, args.asset_root)
    if args.execute_approved_run != run_id:
        raise ValueError('Explicit execution authorization for this contract hash required')

    if not torch.cuda.is_available():
        print(
            "\nGPU UNAVAILABLE: No predictions attempted. "
            "Existing files are safe."
        )
        return

    model, tokenizer = load_model(contract, args.asset_root)

    completed = 0

    for row in pending[:args.limit]:
        destination = rooted_path(args.output_root, row["pdb_path"])
        sequence = row["sequence"]

        print(
            f"\nPredicting {row['target']} "
            f"{row['method']} #{row['index']}"
        )

        # Derive each prediction seed from its stable identity, not processing order.
        import hashlib
        seed = (contract['seed'] + int(hashlib.sha256(row['prediction_id'].encode('ascii')).hexdigest()[:8], 16)) % 2**32
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.use_deterministic_algorithms(contract['deterministic_algorithms'])
        pdb_text = predict(sequence, model, tokenizer, contract['inference_settings'])

        save_prediction(pdb_text, destination, sequence, {'row': row, 'contract': contract, 'rows': rows})

        print("SAVED AND VERIFIED:", destination)
        completed += 1

    print("\nNew predictions completed:", completed)

if __name__ == "__main__":
    main()
