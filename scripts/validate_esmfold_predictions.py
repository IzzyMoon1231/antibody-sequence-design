
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

if __name__ == "__main__":
    main()
