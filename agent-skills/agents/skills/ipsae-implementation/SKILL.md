---
name: ipsae-implementation
description: This skill computes interface predicted TM-score (ipTM) and interface predicted aligned error (ipSAE) metrics, preserving PAE asymmetry and matching reference formulas for accurate confidence scoring.
---

# Goal

Implement ipTM and ipSAE metrics exactly as IPSAE version 4.

# Metrics

Implement:
- ipTM_d0chn
- ipSAE_d0res
- ipSAE_d0chn
- ipSAE_d0dom

# Core Formula

ptm(x, d0) = 1 / (1 + (x / d0)^2)

# d0 Rules

For L > 27:
1.24 * (L - 15)^(1/3) - 1.8

Otherwise:
1.0

# Important Constraints

- PAE matrix is asymmetric
- Preserve chain directionality
- Respect residue indexing
- Never symmetrize PAE

# ipTM Rules

- use all inter-chain PAE values
- no cutoff
- normalize by chain-pair length

# ipSAE Rules

- apply PAE cutoff
- compute residue-specific normalization
- preserve asymmetric scoring

# Validation

Target tolerance:
±0.001

Validate against ipsae.py.

# Dependencies

- numpy for matrix operations
