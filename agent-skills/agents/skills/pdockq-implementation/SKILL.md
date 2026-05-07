---
name: pdockq-implementation
description: This skill calculates pDockQ and pDockQ2 scores based on interface contacts and PAE-derived probabilities, using CB-CA distances and logistic formulas for structural plausibility assessment.
---

# Goal

Implement pDockQ and pDockQ2.

# Contact Definition

Use:
- CB-CB distance <= 8 Å
- CA for glycine

# pDockQ

Compute:
- interface residues
- number of contacts
- interface mean pLDDT

Formula:
0.724 / (1 + exp(-0.052 * (x - 152.611))) + 0.018

# pDockQ2

Combine:
- interface contacts
- PAE-derived ptm values

Use:
ptm(PAE[i][j], 10.0)

Formula:
1.31 / (1 + exp(-0.075 * (x - 84.733))) + 0.005

# Constraints

- preserve interface residue lists
- validate contact counting
- use numpy implementations
- document formulas inline

# Dependencies

- numpy for distance calculations
