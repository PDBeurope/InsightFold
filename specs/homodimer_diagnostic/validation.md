# Homodimer Diagnostic Validation Plan

## Validation Levels

Use four validation levels:

1. `static`: inspect notebook structure, unresolved TODOs, dependency policy, and spec traceability.
2. `smoke`: run the notebook with FX-001 happy-path candidate.
3. `full`: run happy-path, borderline, disagreement, negative, and unsupported fixtures.
4. `reference`: compare IPSAE-family scores against DunbrackLab/IPSAE reference calculations within +/- 0.001 where equivalent inputs are used.

Default target before beta: `full` plus `reference`.

Current fixture readiness:

- `smoke`: ready with FX-001
- `negative`: ready with FX-004 and FX-005
- `full`: blocked pending final FX-002/FX-003 curation and trusted numeric snapshots

## Required Checks

| Check ID | Check | Fixture(s) | Requirement |
|---|---|---|---|
| VAL-001 | Notebook restarts and runs top-to-bottom | FX-001 | REQ-020 |
| VAL-002 | Required AFDB fields are validated before use | FX-001, FX-002, FX-004, FX-005 | REQ-002, REQ-003 |
| VAL-003 | Non-homodimer or malformed inputs fail clearly | FX-004, FX-005 | REQ-001, REQ-007 |
| VAL-004 | PAE matrix shape matches chain metadata | FX-001, FX-002, FX-003 | REQ-005 |
| VAL-005 | mmCIF parser extracts required coordinates | FX-001, FX-002, FX-003 | REQ-004 |
| VAL-006 | pLDDT provenance is stated | FX-001, FX-002, FX-003 | REQ-006 |
| VAL-007 | Contact map and interface residues are produced | FX-001, FX-002, FX-003 | REQ-008, REQ-015 |
| VAL-008 | ipTM is computed and intermediate values are shown | FX-001, FX-002, FX-003 | REQ-009, REQ-014 |
| VAL-009 | ipSAE variants are computed and intermediate values are shown | FX-001, FX-002, FX-003 | REQ-010, REQ-014 |
| VAL-010 | pDockQ and pDockQ2 are computed | FX-001, FX-002, FX-003 | REQ-011, REQ-012 |
| VAL-011 | LIS is computed | FX-001, FX-002, FX-003 | REQ-013 |
| VAL-012 | PAE heatmaps and masks render | FX-001, FX-002, FX-003 | REQ-016 |
| VAL-013 | Per-residue profiles render | FX-001, FX-002, FX-003 | REQ-017 |
| VAL-014 | MolViewSpec minimum views or fallback exist | FX-001 | REQ-018 |
| VAL-015 | Diagnostic summary avoids clinical overclaiming | FX-001 to FX-003 | REQ-019, REQ-021 |
| VAL-016 | Runtime is less than 60 seconds on target runtime | FX-001 | REQ-020 |
| VAL-017 | Disallowed dependencies are absent | all | REQ-020 |

## Fixture-Specific Expected Results

### FX-001

Smoke checks:

- API response is non-empty and `isComplex` is true.
- Required URLs are present.
- PAE chain metadata has two chains of length 172.
- Dimer length is 344.
- Contact count is expected to be 116 using current notebook logic.
- Existing score snapshot from `notebooks/homodimer_diagnostic.ipynb`:
  - `ipSAE_d0res`: 0.9143
  - `ipSAE_d0chn`: 0.9529
  - `ipSAE_d0dom`: 0.9527
  - `ipTM` / `ipTM_d0chn`: 0.9529
  - `pDockQ`: 0.6913
  - `pDockQ2`: 0.9269
  - `LIS`: 0.7564

Use these as provisional smoke snapshots only. Reference validation requires independent regeneration after IPSAE commit/version pinning.

### FX-002

Candidate checks:

- API response is non-empty and `isComplex` is true.
- Required URLs are present.
- Notebook computes the full metric table.
- If `ipSAE_d0res` is near the intended 0.6 boundary, promote FX-002 to official borderline fixture.
- If not, retain FX-002 as a secondary successful complex fixture and request a better borderline accession.

### FX-003

Blocked until AFDB/PDBe domain reviewer supplies an accession with known or suspected metric disagreement.

### FX-004

Expected API response includes an identifier-format error. Notebook should surface this clearly and stop.

### FX-005

API response is valid, but `isComplex` is false. Notebook should reject it as unsupported by v1 before homodimer scoring.

## Reference Validation

Reference source:

- DunbrackLab/IPSAE formulas, version stated in PRD as IPSAE v4 / January 2026.

Tolerance:

- ipSAE-family values should match within +/- 0.001 where equivalent inputs and definitions are used.

Blocking gaps:

- Need frozen reference values for curated fixtures.
- Need decision on exact IPSAE reference version/commit.
- Need explicit handling of cases where local CB-CB parser produces a different contact set than reference due to parsing differences.
- Need final metric-disagreement fixture for summary validation.

## Validation Report Output

Expected path:

```text
specs/homodimer_diagnostic/validation-report.md
```

The report should include:

- command/tool used for execution
- notebook path
- fixture IDs run
- runtime
- dependency versions
- pass/fail table
- cell errors with context
- numeric comparisons
- unresolved risks
