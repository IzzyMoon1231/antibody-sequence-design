# Repertoire similarity and structural cost

## Does adding AbLang improve human-repertoire similarity without a structural cost?

9NH7: higher heavy-domain repertoire similarity and higher H3 backbone RMSD in the ensemble median. 9NFU: higher heavy-domain repertoire similarity and higher H3 backbone RMSD in the ensemble median. The retained median pattern does not support a blanket no-structural-cost claim; the structural shift must be read case by case, alongside spread and confidence. The H3 median increase is larger for 9NH7 (+1.3790 Å) than 9NFU (+0.1205 Å); this is an observed scale difference, not a biological-significance threshold.

Pinned retained-input analysis

## Table 1 · Two objectives, kept distinct

| Case | Heavy-domain Identity, B → E | Δ Identity | H3 RMSD, B → E (Å) | Δ RMSD (Å) |
| --- | --- | --- | --- | --- |
| 9NH7 | 0.5981 → 0.8131 | +0.2150 | 0.8960 → 2.2750 | +1.3790 |
| 9NFU | 0.6549 → 0.8584 | +0.2035 | 2.9390 → 3.0595 | +0.1205 |

Each main cell is a median over 100 designs per arm per case; B = baseline, E = ensemble. Deltas are differences of marginal medians, not paired changes. Higher Identity indicates greater repertoire similarity; higher RMSD indicates greater displacement from the deposited reference.

Across 9NFU loops, median RMSD increases larger than its H3 shift occur in H2, L2, L3; decreases (improvements relative to the reference) occur in H1. The H3 result is not overall scFv preservation. The complete supporting loop table retains all draws. Low confidence means mean loop pLDDT below 70, not any low-confidence residue, and does not explain away every worsening: 9NFU H2 below-70 counts are 0 baseline and 0 ensemble.

Table 1 compares independently sampled all-loop baseline and ensemble protocol cohorts, not paired seeds. Heavy OASis Identity is the whole numbered heavy-domain fraction of qualifying overlapping 9-mers at the reported relaxed threshold: not an OASis percentile, isolated CDR score, or immunogenicity measurement. H3 RMSD uses a per-domain framework N/CA/C backbone fit to deposited experimental coordinates, not the original RFdiffusion model. Loops are operational sequence-alignment/motif definitions, not independently validated biological numbering. Separate domain fits exclude relative VH/VL orientation and cannot establish whole-scFv orientation or interface preservation.

## What differs between cases?

**9NH7.** H3 interquartile ranges are 0.835–0.959 Å (baseline) and 1.121–5.047 Å (ensemble); these ranges do not overlap. The ensemble middle-half spread is higher. Mean H3 RMSD is 0.90177 → 2.91623 Å; framework medians are 1.453 → 1.433 Å. H3 mean pLDDT is below 70 in 0/100 → 32/100 designs. All remain in the primary summaries; low confidence weakens interpretation as a true fold, not the retained arithmetic.

**9NFU.** H3 interquartile ranges are 2.881–2.995 Å (baseline) and 2.722–3.384 Å (ensemble); these ranges overlap. The ensemble middle-half spread is higher. Mean H3 RMSD is 2.94351 → 2.98316 Å; framework medians are 1.008 → 1.020 Å. H3 mean pLDDT is below 70 in 0/100 → 18/100 designs. All remain in the primary summaries; low confidence weakens interpretation as a true fold, not the retained arithmetic.

Neither overlapping ranges nor a small mean shift establishes structural non-inferiority: a meaningful margin and justified uncertainty analysis are needed. Predicted RMSD does not measure binding. These are draw-weighted designs on two fixed cases, not independent biological experiments or evidence of an antibody-format effect. No p-values, confidence intervals, or causal acceptance are inferred here.

## Figure 1 · Where does the humanness difference occur?

100 reported pairs per condition; 800 ablation designs in total. Figure 1 is a separate loop-group generation cohort. Bars are mean ensemble-minus-baseline Heavy OASis Identity differences across reported seed-matched pairs; no error bars are implied. The source defines a same-decoder, same-seed zero-AbLang-logit control, but the executed helper revision/import resolution and weights are not authenticated. Thus these are descriptive retained outputs, not proof that AbLang alone caused the changes. Other CDRs means H1/H2 for 9NH7 and H1/H2 + L1/L2/L3 for 9NFU. The 9NFU heavy score is from joint heavy/light redesign, not an isolated heavy-loop intervention. No matched ablation structural results are supplied.

![Figure 1. Heavy-domain ablation differences](figures/figure1.png)

9NH7: the mean difference outside H3 is higher than H3-only (+0.15093 versus +0.05869); positive reported pairs: 100/100 versus 100/100. 9NFU: the mean difference outside H3 is higher than H3-only (+0.08504 versus +0.05549); positive reported pairs: 100/100 versus 100/100.

9NFU light-domain Identity: cdrh3: mean paired gain +0.0000000, 0/100 positive, 100/100 unchanged; others: mean paired gain +0.1403125, 100/100 positive, 0/100 unchanged. Heavy and light denominators are separate; no cross-chain efficiency inference follows.

## Paper changes

Use Table 1 and the complete supporting loop table to replace stale humanness and loop values; retain the heterogeneous loop-level structural costs, not universal preservation, foldability or therapeutic claims. Describe deposited-reference, motif-defined N/CA/C domain fits and whole heavy-domain Identity consistently. Keep Figure 1 descriptive, with the joint heavy/light other-CDR condition and unchanged-light control; remove per-residue or cross-chain efficiency claims. Remove index-paired Wilcoxon and inherited significance claims, completed contact-preservation claims, and structural-ablation comparisons without matching outputs. H3 is a focused endpoint, not an authenticated prospective designation. Two fixed cases are not a format ablation. These are update instructions, not certification that a manuscript has been edited.

## What this evidence supports

For the pinned original inputs, prior audit linked scored VH to final-filled FASTAs and full predicted sequences, and replayed the implemented coordinate measurements without numeric discrepancies. This package rechecks the scored-VH/FASTA link and recalculates retained metrics; it does not rerun coordinate replay, folding or BioPhi scoring. Computation verification is distinct from independent biological mapping or causal acceptance.

[Supporting distributions and alternatives](APPENDIX.md) · [Sources and scientific limits](SOURCES.md) · [Update the analysis](README.md) · [Executed notebook](results.ipynb) · [Engineering checks](engineering/README.md)
