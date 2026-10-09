"""
Reproduce pLDDT statistical comparisons from saved per-design data.

9NFU: 100 designs per group.
9NH7: 87 unique baseline and 76 unique ensemble sequences.

This reproduces downstream statistics, not ESMFold predictions.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]

def verify_means_and_p(b, e, saved, label):
    test = mannwhitneyu(b, e, alternative="two-sided")
    difference = e.mean() - b.mean()

    checks = {
        "baseline_mean": np.isclose(
            b.mean(), saved["baseline_mean"], atol=1e-6
        ),
        "ensemble_mean": np.isclose(
            e.mean(), saved["ensemble_mean"], atol=1e-6
        ),
        "difference": np.isclose(
            difference, saved["difference"], atol=1e-6
        ),
        "p_value": np.isclose(
            test.pvalue, saved["p_value"],
            rtol=1e-5, atol=0
        )
    }

    print(
        f"{label:20s} "
        f"B={b.mean():.4f} "
        f"E={e.mean():.4f} "
        f"p={test.pvalue:.5g} "
        f"verified={all(checks.values())}"
    )

    if not all(checks.values()):
        print("Failed checks:", checks)

    assert all(checks.values()), f"Verification failed: {label}"

# 9NH7: unique sequences
nh = pd.read_csv(
    ROOT / "results/9NH7/9NH7_final_pLDDT_results.csv"
).drop_duplicates(["group", "sequence"])

nh_saved = pd.read_csv(
    ROOT / "results/9NH7/9NH7_unique_pLDDT_statistics.csv"
)

print("9NH7 — unique sequences")

for _, row in nh_saved.iterrows():
    metric = row["metric"]

    b = nh.loc[nh["group"] == "baseline", metric].dropna()
    e = nh.loc[nh["group"] == "ensemble", metric].dropna()

    assert len(b) == row["baseline_n"]
    assert len(e) == row["ensemble_n"]

    verify_means_and_p(b, e, row, metric)

# 9NFU: all designs
nf = pd.read_csv(
    ROOT / "results/9NFU/9NFU_CDR_pLDDT_results.csv"
)

nf_saved = pd.read_csv(
    ROOT / "results/9NFU/9NFU_pLDDT_statistical_comparison.csv"
)

print("\n9NFU — all designs")

for _, row in nf_saved.iterrows():
    region = row["region"]

    column = (
        "overall_pLDDT" if region == "overall"
        else "non_CDR_pLDDT" if region == "non_CDR"
        else f"{region}_pLDDT"
    )

    b = nf.loc[nf["group"] == "baseline", column].dropna()
    e = nf.loc[nf["group"] == "ensemble", column].dropna()

    assert len(b) == 100
    assert len(e) == 100

    reference = row.copy()
    reference["baseline_mean"] = row["baseline_mean"]
    reference["ensemble_mean"] = row["ensemble_mean"]
    reference["difference"] = row["mean_difference"]

    verify_means_and_p(b, e, reference, region)

print("\nSUCCESS: All pLDDT statistical comparisons verified.")
