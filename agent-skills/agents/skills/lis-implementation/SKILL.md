---
name: lis-implementation
description: This skill computes the Local Interaction Score by averaging normalized low-PAE inter-chain interactions, providing a measure of interaction density without global docking confidence.
---

# Goal

Implement Local Interaction Score (LIS).

# Workflow

1. Extract inter-chain PAE blocks.
2. Filter PAE < 12.
3. Compute:
(12 - PAE) / 12
4. Average valid scores.

# Constraints

- preserve asymmetry
- do not symmetrize prematurely
- explain interpretation

# Dependencies

- numpy for array operations

# Educational Guidance

Explain LIS as:
"density of low-PAE inter-chain interactions"

Clarify:
- LIS measures interaction density
- not global docking confidence
