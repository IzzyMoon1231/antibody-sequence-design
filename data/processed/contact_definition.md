# Antigen Contact Residue Definition

## Purpose

Identify antibody CDR residues that interact with the antigen
in antibody-antigen complex structures.

## Chain mapping

H = antibody heavy chain
L = antibody light chain
T = antigen/target chain

## CDR loops

H1
H2
H3
L1
L2
L3

## Contact definition

Antigen-contact residues are antibody residues whose atoms fall within a
specified distance cutoff of any antigen atom.

For initial pipeline development and validation on 9NFU, a working cutoff
of 5.0 Å is used.

This cutoff is currently a pipeline parameter and should not be interpreted
as the definitive Bennett et al. contact definition until confirmed from
the published methodology or associated code.

Only contacts belonging to CDR residues are retained for the CDR-contact
analysis.

CDR assignments use IMGT numbering through AbNumber/ANARCII.

## 9NFU chain mapping

Antigen:
- Chain A — TcdB

Antibody:
- Chain C — designed scFv

## Current validation result

Using a 5.0 Å working cutoff:

- Total antibody-antigen contact residues: 16
- CDR contact residues: 14
- Non-CDR contact residues: 2

CDR contacts were detected in H1, H2, H3, L1, and L2.
No L3 contact was detected under the current 5.0 Å working definition.

## Output

For every contact residue, record:

- structure ID
- antibody chain
- residue number
- amino acid
- CDR loop
- antigen chain
- closest antigen residue
- minimum distance
- contact status
