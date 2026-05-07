---
name: validation-against-reference
description: This skill compares computed metrics against ipsae.py and other references, ensuring numerical accuracy within tolerance and rejecting outputs that fail validation checks.
---

# Goal

Validate notebook outputs against reference implementations.

# Required Validations

- ipTM
- ipSAE variants
- pDockQ
- pDockQ2
- LIS

# Tolerance

±0.001

# Validation Workflow

1. Compare notebook outputs.
2. Compare intermediate calculations.
3. Validate residue indexing.
4. Validate asymmetric PAE handling.
5. Verify interface residue selection.

# Failure Conditions

Reject outputs if:
- formulas differ
- residue indexing mismatches
- cutoffs incorrect
- asymmetry lost
- tolerance exceeded

# Dependencies

- reference ipsae.py implementation
