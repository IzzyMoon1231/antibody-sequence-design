# Reproducibility Guide

## Project Overview

This project compares two computational antibody sequence-design methods:

- **Baseline:** ProteinMPNN
- **Ensemble:** ProteinMPNN + AbLang

Two antibody-containing experimental structures were evaluated:

- **9NFU:** Single-chain variable fragment (scFv)
- **9NH7:** Single-domain antibody (VHH)

The evaluation examines sequence conservation, contact-position
preservation, predicted structural confidence (pLDDT), and
structural agreement with experimental references (RMSD).

## Reproducibility Status

The repository contains saved per-design measurements and verified
scripts for reproducing downstream statistical comparisons.

All four downstream reproducibility checks passed.

This does not mean the complete antibody generation and
structure-prediction pipeline can currently be reproduced from scratch.

## Requirements

Python packages:

- numpy
- pandas
- scipy
- statsmodels

Install with:

    pip install numpy pandas scipy statsmodels

## Running the Reproducibility Checks

From the repository root, run:

    python scripts/run_reproducibility_checks.py

This runs four scripts:

1. scripts/reproduce_contact_statistics.py
2. scripts/reproduce_sequence_conservation.py
3. scripts/reproduce_plddt_statistics.py
4. scripts/reproduce_rmsd_statistics.py

Successful execution ends with:

    SUCCESS: ALL DOWNSTREAM REPRODUCIBILITY CHECKS PASSED

## Dataset Sizes

| Structure | Baseline | Ensemble | Unique Baseline | Unique Ensemble |
|---|---:|---:|---:|---:|
| 9NFU | 100 | 100 | 100 | 100 |
| 9NH7 | 100 | 100 | 87 | 76 |

The 9NH7 pLDDT and RMSD comparisons use unique sequences.

Other analyses use the design counts specified in their
individual verification scripts.

## Statistical Methods

Baseline and ensemble measurements were compared using
two-sided Mann-Whitney U tests.

Where applicable, p-values were adjusted using the
Benjamini-Hochberg false discovery rate (FDR) procedure.

Mean difference was defined as:

    Ensemble mean - Baseline mean

The direction of a desirable difference depends on the metric.

Generated designs are computational observations, not
independent experimental biological replicates.

## Evaluation Metrics

### Sequence Conservation

Measures retention of reference amino acids at evaluated
sequence positions, according to the saved analysis data.

### Contact-Position Preservation

Uses saved mapped-residue C-alpha displacement measurements,
including a 3-angstrom threshold.

This is a geometric proxy and does not directly establish
preservation of antibody-antigen contacts or binding affinity.

### pLDDT

Measures predicted local structural confidence.

Higher pLDDT reflects greater model confidence, not
necessarily better antibody binding or biological function.

### RMSD

Measures structural deviation from experimental references
under the specified structural alignment.

Lower RMSD indicates closer agreement under that alignment,
not necessarily improved binding.

For 9NFU, the analysis includes domain-specific alignment
and a filtered heavy-framework comparison.

## Selected Results

### Contact-Position Preservation (3 Angstroms)

| Structure | Baseline | Ensemble |
|---|---:|---:|
| 9NFU | 75.21% | 68.86% |
| 9NH7 | 81.87% | 53.20% |

### Sequence Conservation

| Structure | Baseline | Ensemble |
|---|---:|---:|
| 9NFU | 65.50% | 43.86% |
| 9NH7 | 15.27% | 22.00% |

### 9NH7 H3 Structural Metrics

| Metric | Baseline | Ensemble |
|---|---:|---:|
| H3 pLDDT | 84.27 | 73.46 |
| H3 RMSD | 0.918 Angstroms | 2.995 Angstroms |

The ensemble showed mixed effects across structures
and CDR regions rather than consistent improvement.

## Limitations

1. Original predicted structure files used to calculate
   pLDDT and RMSD are not currently included in the repository.

2. The complete upstream structure-prediction and RMSD
   extraction workflows are not currently archived.

3. Statistical scripts reproduce calculations from saved
   measurements rather than regenerating the measurements.

4. Geometric preservation does not establish experimental
   antibody-antigen binding or binding affinity.

5. Predicted structural confidence does not demonstrate
   experimental stability or therapeutic effectiveness.

6. Only two structural targets were evaluated.

7. Full regeneration of ensemble sequences may require
   access to the original AbLang model weights.

## Overall Interpretation

Adding AbLang to ProteinMPNN did not consistently improve
the evaluated structural properties across both antibodies.

Changes differed across structures and CDR regions.

The repository enables verification of downstream statistical
calculations while identifying the remaining requirements
for complete end-to-end reproducibility.
