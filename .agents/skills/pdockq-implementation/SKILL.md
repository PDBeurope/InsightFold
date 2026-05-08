---
name: pdockq-implementation
description: This skill calculates pDockQ and pDockQ2 scores based on interface contacts and PAE-derived probabilities, using CB-CA distances and logistic formulas for structural plausibility assessment.
---

# Goal

Implement pDockQ and pDockQ2.

Use `$homodimer-confidence-scoring` as the project authority for contact definitions, edge cases, and current sigmoid constants. If a PRD, spec, or reference implementation disagrees with this skill, stop and flag the conflict before implementing.

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

For InsightFold homodimer notebooks, use the pDockQ2 constants from `$homodimer-confidence-scoring` or the checked reference implementation. Do not infer constants from memory.

# Constraints

- preserve interface residue lists
- validate contact counting
- use numpy implementations
- document formulas inline
- return floor scores for zero-contact cases

# Dependencies

- numpy for distance calculations
