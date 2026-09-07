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

To be finalized from the Bennett/RFantibody methodology.

For each antibody CDR residue, calculate its distance to
atoms in the antigen.

If the residue satisfies the selected distance criterion,
label it as an antigen-contact residue.

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
