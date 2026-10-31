# Repertoire similarity and structural cost

## Are the sequences more human-like, and what happens to their predicted structures?

Ensemble sequences score higher on human-antibody sequence similarity in both examples (whole heavy-domain medians). 9NH7: higher heavy-domain repertoire similarity and higher H3 backbone RMSD in the ensemble median. 9NFU: higher heavy-domain repertoire similarity and higher H3 backbone RMSD in the ensemble median. The retained median pattern does not support a blanket no-structural-cost claim; the structural shift must be read case by case, alongside spread and confidence. The H3 median increase is larger for 9NH7 (+1.3790 Å) than 9NFU (+0.1205 Å); this is an observed scale difference, not a biological-significance threshold.

Pinned retained-input analysis

Two questions, two separate collections. Main comparison: all-loop redesign, sequence similarity and predicted structure. Ablation: H3-only versus other-loop redesign, sequence similarity only. Arrows show which measurements belong to each collection, not a causal effect.

![Experiment map](figures/experiment_map.png)

## Figure 1 · Sequence similarity and every measured loop

Light sequence dots retain all observed draws (9NH7: B=100, E=100 / 9NFU: B=100, E=100); black ticks mark medians. Structural dots are group medians for every measured loop. Thick lines show the middle half (25th–75th percentiles), not confidence intervals. All draws keep their original weight, including repeated sequences. Baseline and ensemble are separate collections: no row-pairing or connecting designs. Higher Identity means more human-repertoire sequence similarity; lower RMSD means closer to the deposited reference. These are whole heavy-domain scores and domain-local loop distances, not immune-response or binding measurements.

### 9NH7 · similarity

![9NH7 similarity](figures/9NH7_similarity.png)

9NH7: baseline n=100, ensemble n=100. Heavy-domain Identity medians: 0.5981 → 0.8131; ensemble is higher (Δ +0.2150). Light dots are all draws; black ticks are medians; thick lines are the middle 50%, not confidence intervals. This is repertoire similarity, not immune response. Reflect: does the structural panel support the same favorable direction?

### 9NH7 · structure

![9NH7 structure](figures/9NH7_structure.png)

9NH7: baseline n=100, ensemble n=100. Loop median RMSD in Å, baseline → ensemble: H1 2.574 → 2.647 (Δ +0.074); H2 0.801 → 1.184 (Δ +0.383); H3 0.896 → 2.275 (Δ +1.379). Structural dots are medians; lines are the middle 50%, not confidence intervals. The ensemble has higher similarity and higher H3 RMSD. Reflect: is this a similarity/structure trade-off in these inputs? It is not evidence of preserved binding.

### 9NFU · similarity

![9NFU similarity](figures/9NFU_similarity.png)

9NFU: baseline n=100, ensemble n=100. Heavy-domain Identity medians: 0.6549 → 0.8584; ensemble is higher (Δ +0.2035). Light dots are all draws; black ticks are medians; thick lines are the middle 50%, not confidence intervals. This is repertoire similarity, not immune response. The scFv score here is heavy-only, not a pooled VH/VL endpoint. Reflect: does the structural panel support the same favorable direction?

### 9NFU · structure

![9NFU structure](figures/9NFU_structure.png)

9NFU: baseline n=100, ensemble n=100. Loop median RMSD in Å, baseline → ensemble: H1 1.266 → 1.083 (Δ -0.183); H2 1.407 → 2.087 (Δ +0.679); H3 2.939 → 3.059 (Δ +0.120); L1 1.286 → 1.307 (Δ +0.020); L2 1.882 → 2.333 (Δ +0.452); L3 1.032 → 1.639 (Δ +0.607). Structural dots are medians; lines are the middle 50%, not confidence intervals. Across 9NFU loops, median RMSD increases larger than its H3 shift occur in H2, L2, L3; decreases (improvements relative to the reference) occur in H1. The H3 result is not overall scFv preservation. The complete supporting loop table retains all draws. Low confidence means mean loop pLDDT below 70, not any low-confidence residue, and does not explain away every worsening: 9NFU H2 below-70 counts are 0 baseline and 0 ensemble. Reflect: which other loops change the conclusion you would draw from H3 alone?

[Download compact 2×2 Figure 1 (SVG)](figures/figure1.svg) · [PNG](figures/figure1.png)

Across 9NFU loops, median RMSD increases larger than its H3 shift occur in H2, L2, L3; decreases (improvements relative to the reference) occur in H1. The H3 result is not overall scFv preservation. The complete supporting loop table retains all draws. Low confidence means mean loop pLDDT below 70, not any low-confidence residue, and does not explain away every worsening: 9NFU H2 below-70 counts are 0 baseline and 0 ensemble.

## Table 1 · Numerical reference

| Case | Heavy-domain Identity, B → E | Δ Identity | H3 RMSD, B → E (Å) | Δ RMSD (Å) |
| --- | --- | --- | --- | --- |
| 9NH7 | 0.5981 → 0.8131 | +0.2150 | 0.8960 → 2.2750 | +1.3790 |
| 9NFU | 0.6549 → 0.8584 | +0.2035 | 2.9390 → 3.0595 | +0.1205 |

Each main cell is a marginal median; selected draw counts: 9NH7: B=100, E=100 / 9NFU: B=100, E=100; B = baseline, E = ensemble. Deltas are differences of marginal medians, not paired changes. Higher Identity indicates greater repertoire similarity; higher RMSD indicates greater displacement from the deposited reference.

## Figure 2 · Where are differences seen in the separate ablation?

Reported pair counts: 9NFU/cdrh3/H: 100 / 9NFU/cdrh3/L: 100 / 9NFU/others/H: 100 / 9NFU/others/L: 100 / 9NH7/cdrh3/H: 100 / 9NH7/others/H: 100; 800 ablation designs in total. Figure 2 uses a separate loop-group generation cohort. Light dots show individual reported seed-pair differences; dark diamonds show the mean and thick lines the middle half (IQR), not confidence intervals. Positive differences favor ensemble sequence similarity. The source defines a same-decoder, same-seed zero-AbLang-logit control, but the executed helper revision/import resolution and weights are not authenticated. Thus these are descriptive retained outputs, not proof that AbLang alone caused the changes. Other CDRs means H1/H2 for 9NH7 and H1/H2 + L1/L2/L3 for 9NFU. The 9NFU heavy score is from joint heavy/light redesign, not an isolated heavy-loop intervention. No matched ablation structural results are supplied.

![Figure 2. Reported seed-pair differences](figures/figure2.png)

9NH7: the mean heavy-domain Identity difference (ensemble minus baseline; pp = percentage points) outside H3 is higher than H3-only (+15.09 pp versus +5.87 pp); positive reported pairs: 100/100 versus 100/100.

9NFU: the mean heavy-domain Identity difference (ensemble minus baseline; pp = percentage points) outside H3 is higher than H3-only (+8.50 pp versus +5.55 pp); positive reported pairs: 100/100 versus 100/100.

9NFU light-domain Identity: cdrh3: mean paired gain +0.0000000, 0/100 positive, 100/100 unchanged; others: mean paired gain +0.1403125, 100/100 positive, 0/100 unchanged. Heavy and light denominators are separate; no cross-chain efficiency inference follows.

## What this evidence supports

Neither overlapping ranges nor a small mean shift establishes structural non-inferiority: a meaningful margin and justified uncertainty analysis are needed. Predicted RMSD does not measure binding. These are draw-weighted designs on two fixed cases, not independent biological experiments or evidence of an antibody-format effect. No p-values, confidence intervals, or causal acceptance are inferred here.

For the pinned original inputs, prior audit linked scored VH to final-filled FASTAs and full predicted sequences, and replayed the implemented coordinate measurements without numeric discrepancies. This package rechecks the scored-VH/FASTA link and recalculates retained metrics; it does not rerun coordinate replay, folding or BioPhi scoring. Computation verification is distinct from independent biological mapping or causal acceptance.

[Spread and confidence details](APPENDIX.md) · [Sources and methods](SOURCES.md) · [Update guide](README.md) · [Executed notebook](results.ipynb) · [Author checklist — Paper changes](AUTHOR_CHECKLIST.md) · [Engineering checks](engineering/README.md)
