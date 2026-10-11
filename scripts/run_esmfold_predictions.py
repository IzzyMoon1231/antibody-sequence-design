
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

def load_manifest(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError("Prediction manifest is empty")

    required = {"target", "method", "index", "sequence", "pdb_path"}
    if not required.issubset(rows[0]):
        raise ValueError("Manifest missing required columns")

    return rows

def check_prediction(path, sequence):
    path = Path(path)
    marker = path.with_name(path.name + ".complete")

    if not path.is_file():
        if marker.exists():
            return "ORPHAN_COMPLETION_MARKER"
        return "MISSING"

    if not marker.is_file():
        return "INCOMPLETE_NO_MARKER"

    try:
        if marker.read_text().strip() != "COMPLETE":
            return "INVALID_COMPLETION_MARKER"

        observed = extract_pdb_sequence(path)
        status, _ = validate_sequence(sequence, observed)
        return status
    except Exception:
        return "READ_ERROR"

def load_model():
    if not torch.cuda.is_available():
        raise RuntimeError(
            "CUDA GPU unavailable. Cannot run ESMFold predictions."
        )

    try:
        import transformers
    except ImportError as exc:
        raise RuntimeError(
            "Transformers is not installed. Install a compatible "
            "ESMFold environment before running predictions."
        ) from exc

    from transformers import EsmForProteinFolding, AutoTokenizer

    model_id = "facebook/esmfold_v1"

    print("Loading ESMFold model:", model_id)

    tokenizer = AutoTokenizer.from_pretrained(model_id)

    model = EsmForProteinFolding.from_pretrained(
        model_id,
        low_cpu_mem_usage=True,
    )

    model = model.eval().cuda()

    # Use half precision for the language model only.
    # The folding trunk remains in float32.
    model.esm = model.esm.half()

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


def predict(sequence, model, tokenizer):
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
        output = model(**inputs)

    return serialize_esmfold_output(output)


def save_prediction(pdb_text, destination, sequence):
    import os
    import tempfile
    from pathlib import Path

    destination = Path(destination)
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

            # Publish completion only after the copied PDB validates.
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


def inspect_prediction_locks(manifest_path, minimum_age_minutes=60):
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
        destination = Path(row["pdb_path"])
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

    args = parser.parse_args()

    if args.limit < 1:
        parser.error("--limit must be at least 1")

    rows = load_manifest(args.manifest)

    selected = [
        row for row in rows
        if (args.target is None or row["target"] == args.target)
        and (args.method is None or row["method"] == args.method)
    ]

    statuses = Counter()
    pending = []

    for row in selected:
        status = check_prediction(
            row["pdb_path"], row["sequence"]
        )

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

    if not torch.cuda.is_available():
        print(
            "\nGPU UNAVAILABLE: No predictions attempted. "
            "Existing files are safe."
        )
        return

    model, tokenizer = load_model()

    completed = 0

    for row in pending[:args.limit]:
        destination = Path(row["pdb_path"])
        sequence = row["sequence"]

        print(
            f"\nPredicting {row['target']} "
            f"{row['method']} #{row['index']}"
        )

        pdb_text = predict(sequence, model, tokenizer)

        save_prediction(pdb_text, destination, sequence)

        print("SAVED AND VERIFIED:", destination)
        completed += 1

    print("\nNew predictions completed:", completed)

if __name__ == "__main__":
    main()
