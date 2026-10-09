# Antibody Sequence Design: Structural Evaluation

## Objective

Compare ProteinMPNN baseline sequences with sequences generated
using a ProteinMPNN + AbLang ensemble approach.

This analysis evaluates predicted structural confidence (ESMFold
pLDDT) and similarity to experimental reference structures (RMSD).

## Dataset

| Structure | Antibody type | Baseline | Ensemble |
|---|---|---:|---:|
| 9NH7 | VHH | 100 designs | 100 designs |
| 9NFU | scFv | 100 designs | 100 designs |

Unique sequences:
- 9NH7: 87 baseline, 76 ensemble
- 9NFU: 100 baseline, 100 ensemble

No identical sequences were shared between baseline and ensemble
within either structure.

## Methods

### Predicted structural confidence

ESMFold-predicted structures were sequence-verified against
their corresponding FASTA records.

Mean pLDDT was calculated for each design and evaluated by region.
Higher pLDDT indicates greater prediction confidence, not
experimentally confirmed accuracy or binding affinity.

### Experimental-reference RMSD

Predicted structures were compared with experimental reference
structures using matched residues and structural superposition.

9NH7 used framework-based alignment.

9NFU CDR comparisons used domain-specific framework alignment.
The heavy-domain alignment excluded H123 and H124 in a documented
sensitivity analysis because of anomalous reference mapping.
Original analyses are retained separately.

These values are experimental-reference RMSDs, NOT
RFdiffusion backbone RMSDs.

### Statistics

Exploratory comparisons used two-sided Mann-Whitney U tests
and Benjamini-Hochberg false discovery rate correction.

For 9NH7, primary comparisons used unique sequences to avoid
counting identical designs multiple times.

These computational designs are not independent biological
replicates.

## Main findings

### 9NH7

The ensemble increased mean H1 pLDDT by approximately 9.51
points but decreased H2 by 4.65 points and H3 by 10.80 points.

The baseline had lower mean experimental-reference RMSD
for H2 and H3. H1 RMSD was similar between groups.

### 9NFU

The ensemble had lower mean pLDDT across the evaluated regions.

For domain-aligned CDR RMSD, the ensemble was closer to the
experimental reference in H1, whereas the baseline was closer
in H2, L1, L2, and L3. H3 differences were small.

## Interpretation

AbLang integration did not consistently improve predicted
structural confidence or similarity to experimental references.
Effects varied by antibody structure and CDR region.

## Limitations

- Only two antibody structures were evaluated.
- No experimental binding or stability measurements were performed.
- Structural evaluations rely on computational predictions.
- RMSD depends on reference choice, residue mapping, and alignment.
- Sequence uniqueness is not a complete measure of diversity.
- Statistical tests are exploratory and do not represent
  independent biological replication.
- Improved humanness, antigen-contact recovery, and binding
  performance have not been established by these analyses.

## Files

- `9NH7/`: 9NH7 results and statistics
- `9NFU/`: 9NFU results and statistics
- `summary/`: sequence uniqueness results
- `figures/`: pLDDT and RMSD comparison figures
- `file_manifest.csv`: inventory of original consolidated files

## Reproducibility

Full reproduction additionally requires original FASTA sequences,
predicted PDB files, reference structures, model versions,
generation parameters, residue mapping, and analysis scripts.
