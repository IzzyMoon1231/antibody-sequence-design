
"""
Reproduce downstream RMSD statistics from saved per-design CSVs.

9NFU: 100 designs per group, eight regions.
9NH7: 87 unique baseline and 76 unique ensemble designs,
       four regions, with Benjamini-Hochberg FDR correction.

Does not reproduce structural predictions, residue alignment,
or the original RMSD measurements.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]

# --------------------------------------------------
# 9NFU: verify against saved statistical comparisons
# --------------------------------------------------

raw = pd.read_csv(
    ROOT / "results/9NFU/9NFU_domain_aligned_CDR_RMSD.csv"
)

saved = pd.read_csv(
    ROOT / "results/9NFU/9NFU_domain_aligned_RMSD_statistics.csv"
)

regions = {
    "Heavy framework (filtered)": "H_framework_filtered_RMSD",
    "Light framework": "L_framework_original_RMSD",
    "H1": "H1_RMSD",
    "H2": "H2_RMSD",
    "H3": "H3_RMSD",
    "L1": "L1_RMSD",
    "L2": "L2_RMSD",
    "L3": "L3_RMSD",
}

print("9NFU RMSD verification")
print("=" * 65)

assert set(saved["region"]) == set(regions)

p_values_9nfu = []

for region, column in regions.items():
    b = raw.loc[
        raw["group"] == "baseline", column
    ].dropna()

    e = raw.loc[
        raw["group"] == "ensemble", column
    ].dropna()

    assert len(b) == 100 and len(e) == 100

    test = mannwhitneyu(b, e, alternative="two-sided")
    reference = saved.loc[saved["region"] == region].iloc[0]

    checks = {
        "baseline mean": np.isclose(
            b.mean(), reference["baseline_mean"], atol=1e-6
        ),
        "ensemble mean": np.isclose(
            e.mean(), reference["ensemble_mean"], atol=1e-6
        ),
        "difference": np.isclose(
            e.mean() - b.mean(), reference["difference"], atol=1e-6
        ),
        "p-value": np.isclose(
            test.pvalue, reference["p_value"],
            rtol=1e-5, atol=0
        )
    }

    assert all(checks.values()), f"{region}: {checks}"

    p_values_9nfu.append(test.pvalue)

    print(
        f"{region:28s} "
        f"B={b.mean():.4f} "
        f"E={e.mean():.4f} "
        f"p={test.pvalue:.5g} VERIFIED"
    )

# Verify FDR correction against saved results
_, adjusted_9nfu, _, _ = multipletests(
    p_values_9nfu, method="fdr_bh"
)

for region, adjusted_p in zip(regions, adjusted_9nfu):
    reference = saved.loc[saved["region"] == region].iloc[0]

    assert np.isclose(
        adjusted_p,
        reference["FDR_adjusted_p"],
        rtol=1e-5,
        atol=0
    ), f"FDR mismatch: {region}"

print("9NFU FDR correction verified.")

# --------------------------------------------------
# 9NH7: calculate unique-sequence RMSD comparisons
# --------------------------------------------------

baseline = pd.read_csv(
    ROOT / "results/9NH7/9NH7_baseline_unique_experimental_reference_RMSD.csv"
)

ensemble = pd.read_csv(
    ROOT / "results/9NH7/9NH7_ensemble_unique_RMSD_and_pLDDT.csv"
)

assert len(baseline) == 87
assert len(ensemble) == 76

columns = [
    "framework_rmsd",
    "H1_rmsd",
    "H2_rmsd",
    "H3_rmsd"
]

results = []

for column in columns:
    b = baseline[column].dropna()
    e = ensemble[column].dropna()

    test = mannwhitneyu(
        b, e, alternative="two-sided"
    )

    results.append({
        "region": column,
        "baseline_mean": b.mean(),
        "ensemble_mean": e.mean(),
        "difference": e.mean() - b.mean(),
        "U_statistic": test.statistic,
        "p_value": test.pvalue
    })

df = pd.DataFrame(results)

reject, adjusted, _, _ = multipletests(
    df["p_value"].to_numpy(),
    alpha=0.05,
    method="fdr_bh"
)

df["FDR_adjusted_p"] = adjusted
df["significant_FDR_0.05"] = reject

print("\n9NH7 RMSD statistics")
print("=" * 65)
print(df.to_string(index=False))

print("\nSUCCESS: RMSD statistical analysis completed.")
