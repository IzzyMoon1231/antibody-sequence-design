
"""
Run all downstream reproducibility checks for the project.

These checks verify statistical analyses from saved CSV outputs.
They do NOT reproduce upstream structure prediction or original
structural measurements.
"""

from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

scripts = [
    "scripts/reproduce_contact_statistics.py",
    "scripts/reproduce_sequence_conservation.py",
    "scripts/reproduce_plddt_statistics.py",
    "scripts/reproduce_rmsd_statistics.py",
]

failures = []

print("ANTIBODY DESIGN REPRODUCIBILITY CHECK")
print("=" * 70)

for script in scripts:
    print(f"\nRunning: {script}")
    print("-" * 70)

    result = subprocess.run(
        [sys.executable, str(ROOT / script)],
        cwd=ROOT,
        capture_output=True,
        text=True
    )

    print(result.stdout)

    if result.stderr.strip():
        print("STDERR:")
        print(result.stderr)

    if result.returncode != 0:
        failures.append(script)
        print(f"FAILED: {script}")
    else:
        print(f"PASSED: {script}")

print("\n" + "=" * 70)

if failures:
    print("REPRODUCIBILITY CHECK FAILED")
    print("Failed scripts:")
    for script in failures:
        print(" -", script)
    raise SystemExit(1)

print("SUCCESS: ALL DOWNSTREAM REPRODUCIBILITY CHECKS PASSED")
print("Verified analyses:")
print(" - Contact-position preservation")
print(" - Sequence conservation")
print(" - pLDDT statistics")
print(" - RMSD statistics")
