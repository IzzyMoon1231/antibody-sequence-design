
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

    if not path.is_file():
        return "MISSING"

    try:
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

def predict(sequence, model, tokenizer):
    from transformers.models.esm.openfold_utils.protein import (
        Protein as OFProtein,
        to_pdb,
    )

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

    atom37_positions = output.positions[-1][0].float().cpu().numpy()
    aatype = output.aatype[0].cpu().numpy()
    atom37_mask = output.atom37_atom_exists[0].cpu().numpy()
    residue_index = (
        output.residue_index[0].cpu().numpy() + 1
    )
    plddt = output.plddt[0].float().cpu().numpy()

    protein = OFProtein(
        aatype=aatype,
        atom_positions=atom37_positions,
        atom_mask=atom37_mask,
        residue_index=residue_index,
        b_factors=plddt,
        chain_index=None,
    )

    return to_pdb(protein)


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

        os.link(temp_path, destination)
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
