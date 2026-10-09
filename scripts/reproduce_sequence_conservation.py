"""
Reproduce contact-position sequence-conservation statistics.

Uses saved per-design conservation measurements.
Does not measure binding affinity or reproduce original
contact-position identification.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed"

baseline_df = pd.read_csv(
    DATA / "week2_sequence_conservation_baseline.csv"
)
ensemble_df = pd.read_csv(
    DATA / "week2_sequence_conservation_ensemble.csv"
)
saved = pd.read_csv(
    DATA / "week2_sequence_conservation_stats.csv"
)

for structure in sorted(saved["structure"].unique()):
    baseline = baseline_df.loc[
        baseline_df["structure"] == structure,
        "contact_conservation_pct"
    ].dropna()

    ensemble = ensemble_df.loc[
        ensemble_df["structure"] == structure,
        "contact_conservation_pct"
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
        "baseline mean": np.isclose(
            round(baseline.mean(), 2),
            reference["baseline_mean_pct"]
        ),
        "ensemble mean": np.isclose(
            round(ensemble.mean(), 2),
            reference["ensemble_mean_pct"]
        ),
        "difference": np.isclose(
            round(difference, 2),
            reference["difference_percentage_points"]
        ),
        "U statistic": np.isclose(
            test.statistic,
            reference["U_statistic"]
        ),
        "p-value (4 significant figures)": (
            f"{test.pvalue:.4g}" ==
            f"{reference['p_value']:.4g}"
        )
    }

    print(f"\\n{structure}")
    print(f"Baseline: {baseline.mean():.6f}%")
    print(f"Ensemble: {ensemble.mean():.6f}%")
    print(f"Difference: {difference:.6f} percentage points")
    print(f"U statistic: {test.statistic}")
    print(f"p-value: {test.pvalue:.10g}")
    print("Verification:", checks)

    assert all(checks.values()), (
        f"Verification failed for {structure}"
    )

print("\\nSUCCESS: Both structures verified.")
