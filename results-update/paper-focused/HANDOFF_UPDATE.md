# Updating the October 12 handoff

**Use the current checked comparison below; keep the old handoff as a historical record.** The old “baseline rejected” wording described a mismatch with the uploaded sequence files used in that review. It did not establish that the team's scoring or folding was wrong. Those blanket labels do not apply to the selected current baseline/ensemble inputs checked for this notebook.

## Which outputs replace the old ones?

Paths in the left column refer to the October 12 handoff ZIP. Replacement links refer to this notebook package's included, checked snapshot—not automatically to newer files on someone's computer.

| Old handoff output | What to use now |
|---|---|
| `evidence/per-draw.csv`, `summary.json`, `TABLES.md`, `RESULTS.md` | For the selected **baseline-versus-ensemble, heavy-domain** comparison, use [Table 1](tables/table1.csv), the [complete loop table](tables/loops.csv), and the [results explanation](READ_RESULTS.md). The old files are not current results. This is a narrower replacement, not a refreshed three-method analysis. |
| `three-method-means-and-duplicates.svg` / `.pdf` | Use [main Figure 1](figures/figure1.svg) for the selected two-method comparison. It does not replace AbMPNN, domain-sensitivity or duplicate-weighting conclusions. |
| `00-ABLATION-DISCREPANCY.md` | For the current separate ablation summary, use the [chain-specific table](tables/ablation.csv) and [Figure 2](figures/figure2.svg). Keep the old discrepancy report as history. Current heavy-domain results must not be substituted for old pooled heavy/light endpoints without changing the description. |

The old rejection statements occur in the ZIP's `README.md`, `findings/COLLABORATOR_INTAKE.md`, and `engineering/ENGINEERING_HANDOFF.md`. Treat them as historical findings about the files compared then, not instructions to discard current measurements.

## What does not need repeating just because a FASTA upload was stale?

Do not automatically repeat OASis scoring, ESMFold prediction, RMSD calculation or ablation. First establish which sequence file actually produced each measurement. For the included selected baseline/ensemble inputs, scored heavy-domain sequences and retained predicted sequences were linked to the current FASTAs, and the implemented structural measurements were replayed successfully. See [source definitions and verification limits](SOURCES.md).

Keep completed reruns and their input records. A rerun on the same inputs and settings can check reproducibility; a run on different sequences is a corrected or different comparison. The current checked cohort is not identical to every older delivery, so an old verification result must not be transferred to new files.

## What remains outside this replacement?

Older **three-method comparisons, extended-domain sensitivity, unique-sequence weighting, peptide/context decomposition, and contact analyses** are not all refreshed by this notebook. Standalone unchanged ensemble or AbMPNN measurements are not invalidated merely because the baseline changed, but comparisons involving the old baseline need corresponding current results.

If these extra claims are no longer used in the paper, archive their outputs—no rerun is needed. If a claim is retained, use an existing corrected analysis where available or update only that specific calculation from matching inputs. The notebook does not establish contact preservation, binding retention or reduced immunogenicity.

## For future updates

Follow the [student guide and change-to-rerun instructions](README.md) and run [the notebook](results.ipynb) with the intended compatible files. Use the newly exported tables and figures together; the checked-in files remain the included snapshot. Optional [statistical calculations](STATISTICS.md) are exploratory: the team still needs to settle the test family and assumptions before adopting inferential claims.

This note supersedes the old reporting guidance for the named current comparison. It does not edit the manuscript, erase historical receipts, or certify other cohorts. Prepared with AI-assisted tracing and checks; author review remains necessary.
