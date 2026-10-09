"""
Reproduce 3 Å C-alpha positional-preservation statistics.

Inputs:
    data/processed/week2_per_design_preservation.csv
    data/processed/week2_statistical_tests.csv

The preservation metric describes mapped C-alpha displacement
relative to the reference, not direct antibody-antigen binding.

Run:
    python scripts/reproduce_contact_statistics.py
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"

designs = pd.read_csv(
    DATA / "week2_per_design_preservation.csv"
)
saved = pd.read_csv(
    DATA / "week2_statistical_tests.csv"
)

for structure in sorted(designs["structure"].unique()):
    subset = designs[designs["structure"] == structure]

    baseline = subset.loc[
        subset["method"] == "baseline",
        "preservation_3A_pct"
    ].dropna()

    ensemble = subset.loc[
        subset["method"] == "ensemble",
        "preservation_3A_pct"
    ].dropna()

    assert len(baseline) == 100
    assert len(ensemble) == 100

    test = mannwhitneyu(
        baseline, ensemble, alternative="two-sided"
    )

    reference = saved.loc[
        saved["structure"] == structure
    ].iloc[0]

    difference = ensemble.mean() - baseline.mean()

    checks = {
        "baseline_mean": np.isclose(
            baseline.mean(),
            reference["baseline_mean_pct"],
            atol=0.0051
        ),
        "ensemble_mean": np.isclose(
            ensemble.mean(),
            reference["ensemble_mean_pct"],
            atol=0.0051
        ),
        "difference": np.isclose(
            difference,
            reference["difference_percentage_points"],
            atol=0.0051
        ),
        "U_statistic": np.isclose(
            test.statistic,
            reference["U_statistic"]
        ),
        "p_value": np.isclose(
            test.pvalue,
            reference["p_value"],
            rtol=0.0001,
            atol=0
        )
    }

    print(f"\\n{structure}")
    print(f"Baseline: {baseline.mean():.4f}%")
    print(f"Ensemble: {ensemble.mean():.4f}%")
    print(f"Difference: {difference:.4f} percentage points")
    print(f"U statistic: {test.statistic}")
    print(f"p-value: {test.pvalue:.10g}")
    print("Verification:", checks)

    assert all(checks.values()), (
        f"Verification failed for {structure}"
    )

print("\\nSUCCESS: Both structures verified.")
