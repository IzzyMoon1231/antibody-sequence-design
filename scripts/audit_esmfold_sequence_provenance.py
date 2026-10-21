#!/usr/bin/env python3
"""Compare archived ESMFold structures to the exact finalized FASTA inputs.

Run: python scripts/audit_esmfold_sequence_provenance.py /path/to/collaborator/repo
No third-party dependencies. Exit status 1 means at least one discrepancy.
"""
import argparse
import csv
from collections import Counter
from pathlib import Path

AA3 = dict(zip(
    "ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCEQGHILKMFPSTWYV",
))
INPUTS = {
    ("baseline", "9NFU"): "baseline_sequences/baseline_sequences_final_filled/9NFU.fa",
    ("baseline", "9NH7"): "baseline_sequences/baseline_sequences_final_filled/9NH7_EBH.fa",
    ("ensemble", "9NFU"): "ensemble_sequences/ensemble_sequences_final_filled/9NFU.fa",
    ("ensemble", "9NH7"): "ensemble_sequences/ensemble_sequences_final_filled/9NH7.fa",
    ("abmpnn", "9NFU"): "AbMPNN_sequences/9NFU_AbMPNN_100_FINAL_FILLED.fa",
    ("abmpnn", "9NH7"): "AbMPNN_sequences/9NH7_AbMPNN_100_FINAL_FILLED.fa",
}

def read_fasta(path):
    records, seq = [], []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq:
                records.append("".join(seq).upper().replace("-", ""))
            seq = []
        else:
            seq.append(line)
    if seq:
        records.append("".join(seq).upper().replace("-", ""))
    return records

def pdb_sequence(path):
    seen, sequence = set(), []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM  ") or line[12:16].strip() != "CA":
            continue
        key = (line[21], line[22:26], line[26])
        if key not in seen:
            seen.add(key)
            sequence.append(AA3.get(line[17:20].strip(), "X"))
    return "".join(sequence)

def audit(root):
    for (method, target), relative in INPUTS.items():
        fasta = root / relative
        if not fasta.is_file():
            yield method, target, "", "", "", "missing_fasta", str(fasta)
            continue
        records = read_fasta(fasta)
        if not records:
            yield method, target, "", 0, "", "empty_fasta", str(fasta)
        folder = root / "esmfold_results" / method / target
        for index, expected in enumerate(records):
            name = f"{target}_{method}_{index:03d}"
            path = folder / f"{name}.pdb"
            if not path.is_file():
                yield method, target, name, len(expected), "", "missing_pdb", str(path)
                continue
            observed = pdb_sequence(path)
            status = "match" if observed == expected else ("empty_pdb" if not observed else "sequence_mismatch")
            yield method, target, name, len(expected), len(observed), status, str(path)
        for path in sorted(folder.glob(f"{target}_{method}_*.pdb")):
            try:
                index = int(path.stem.rsplit("_", 1)[-1])
            except ValueError:
                continue
            if index >= len(records):
                yield method, target, path.stem, "", len(pdb_sequence(path)), "extra_pdb", str(path)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository", type=Path)
    parser.add_argument("--output", type=Path, default=Path("esmfold_sequence_audit.csv"))
    args = parser.parse_args()
    rows = list(audit(args.repository))
    with args.output.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(["method", "target", "design_id", "fasta_length", "pdb_length", "status", "pdb_path"])
        writer.writerows(rows)
    counts = Counter(row[5] for row in rows)
    print("Sequence audit:", dict(counts))
    print("Saved:", args.output)
    if not rows or any(row[5] != "match" for row in rows):
        raise SystemExit(1)

if __name__ == "__main__":
    main()
