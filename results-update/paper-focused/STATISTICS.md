# Optional exploratory statistics

The default `results.ipynb` runs the descriptive workflow with statistics and both exports off. First-time readers can skip section 9b. `results-exploratory.ipynb` is a separately executed example with statistics on and exports off, not a new default or approved protocol. Run either from this package folder. No private sibling directory, manuscript, scorer, model or network service is required.

## What is computed

Four unpaired comparisons: heavy relaxed BioPhi Identity and domain-local framework N/CA/C H3 RMSD, separately for 9NH7 and 9NFU. All retained draws keep their weight, including repeated sequences and low-confidence predictions. Ensemble minus baseline is the direction; larger RMSD means greater deviation. Main arms are never seed-paired.

Each row includes marginal-median difference, all-pair Hodges–Lehmann (HL) difference, 95% percentile independently resampled-arm bootstrap interval, U, rank-biserial effect and U/(nE*nB) = P(E>B) + half P(E=B), plus fresh two-sided raw permutation p. MW is not generally a median or HL test; a location-shift interpretation needs comparable distributional shapes. Narrow discrete intervals are possible; intervals are nonsimultaneous. An interval crossing zero and MW p below .05 need not contradict each other.

Default opt-in settings: 100,000 requested permutations, 10,000 bootstrap resamples, seed 42 per test. Exact label enumeration is used up to 50,000 allocations; equal-sized complete separation uses 2/C(2n,n). Otherwise label permutations retain ties and use plus-one correction. Association ranks are averaged; exhaustive rank permutations apply up to eight observations, otherwise Monte Carlo plus-one. Comparison `n_permutations` is the **requested** budget, not draws used on an exact branch. Each row names its p method. Monte Carlo floor is resolution-limited, never p=0.

These calculations assume independent computational draws and exchangeability under the relevant null. They do not establish that assumption. Repeated identical sequences do not by themselves prove dependence or independence. Known copies, shared trajectories, selection or cached repeated measurements require a reviewed sampling/cluster model before interpreting iid p/intervals. No biological replication, mechanism, binding preservation, reduced immunogenicity or AbLang-only causal claim follows.

## Associations: original only, fail closed for changed data

`inputs/audited_join_rows.csv` retains the 400 previously checked scored VH/full FASTA/actual PDB CA-sequence links and H3 measurements. `inputs/association_receipt.json` binds those rows to original workbook/FASTA/measurement hashes and the previously audited coordinate hashes. Its own hash is pinned in `supplemental.py`. This is reuse of a historical audit, not a fresh coordinate replay or a generation execution certificate. The original audit adapter's hash is recorded, not imported.

Current selected bytes are hashed and validated again, and the selected keys and numeric values must match the pinned join. Four Spearman associations are then calculated separately by case and arm. Missing/modified receipt, missing join, changed main files or differing selected rows withhold correlations. Changed ablation alone does not affect the main association cohort. A valid new workbook/CSV schema or an ID-only match does **not** grant association permission.

For new data the missing evidence is exactly: scored domain → selected full FASTA → sequence extracted from the actual predicted coordinates → content-bound measurement lineage specifying those coordinate bytes, reference, domain/framework fit, atom selection and loop convention. This public adapter does not accept a user-written assertion as proof; ask for a reviewed new receipt/adapter. Group effects remain computable conditional on the selected endpoint labels and sampling assumptions; intake compatibility does not authenticate new measurements.

## Review the whole family, not whichever p is favorable

- Draft six: 9NH7 Identity MW; 9NFU Identity MW; 9NH7 baseline rho; 9NH7 ensemble rho; 9NFU baseline rho; 9NFU ensemble rho.
- Expanded eight adds 9NH7 H3 MW and 9NFU H3 MW.
- Historical scripts corrected groups of two and four separately, not the draft joint six. The draft's six leaves the structural inferential scope unresolved.

The notebook displays every member of six/BH, eight/BH and eight/Holm as **candidate** adjustments. BH needs independence/appropriate positive dependence and targets FDR; Holm controls familywise error without a dependence requirement. Neither is adopted. Do not select the family/method by favorable results. An unavailable member withholds that entire family's correction. No stars, automatic significance conclusions, filtered correlations, other-loop tests or ablation p-values are supplied.

Fill `STATS_REVIEW_PLAN` with the intended claim, sampling unit, endpoints, proposed full family/method, reviewer and unresolved assumptions. This records a review proposal, not approval. Statistical-owner review and publication-policy review remain outstanding; this is not completion of all team statistics.

## What needs rerunning?

| Change | Earliest affected work | Then |
|---|---|---|
| Notes, captions or plotting choices only | Notebook presentation | Restart/run all; keep measurements fixed; new export if desired |
| Already measured scores/CSV rows selected differently | Intake and validation | Refresh descriptive tables/plots and optional effects; associations withheld without current lineage |
| Full sequences or scoring definition changes | Obtain/recompute affected compatible scores with the team's reviewed scorer | Establish score-to-sequence links; changed predictions need corresponding structural measurements; rerun intake onward |
| Predicted coordinates, reference, fit atoms, domains or loop mapping change | Recompute affected structural measurements using the reviewed measurement pipeline | New content-bound measurement/join audit, then downstream summaries and optional stats |
| Statistical estimand, sampling unit or family changes | Scientific review first | Recompute appropriate statistics; no need to generate/fold merely to change a correction |

This notebook does not generate, score, fold or measure coordinates. Recover original realised artifacts before proposing new experiments. No hardcoded upstream command is supplied because the correct runner depends on the accepted run configuration.

## Export and replay

`EXPORT_EXPLORATORY_STATS = True` is separate from the descriptive export. It requires current results and a nonblank review plan, recomputes before writing, and creates a new `exports/statistics-<UTC>-<unique>/` bundle. It includes exact tables, raw/candidate p, methods, settings, assumptions, all family members, input/code hashes, versions, selected bytes and (when qualified) the pinned source receipt/join rows. Old bundles are never overwritten. Save your executed notebook with it; reset both export switches before restarting. Nothing is published or sent.

## Verification boundaries

The helper is copied byte-for-byte from the tested extension. Public helper tests retain its three arithmetic/join/family tests; private original-adapter tests are not shipped as broken sibling dependencies. Public integration tests replace those adapter/bundle paths with the self-contained snapshot and changed/missing/stale controls. The original five-test suite is also rerun separately during maintainer verification. Automated execution is not a supervised beginner usability study or independent scientific endorsement.
