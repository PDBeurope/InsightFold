# Synthetic Fixtures For NMR Restraints Notebook

This directory contains repo-authored synthetic fixtures used to validate edge-case behavior that is not reliably guaranteed by deposited entries.

## Files

- `wrapped_dihedral_minimal.tsv`
  - Hand-computed wrapped and unwrapped dihedral interval expectations.
  - Used to validate circular-geometry logic for intervals such as `170` to `-170`.

- `partial_mapping_warning.str`
  - Minimal NMR-STAR-like distance loop with 10 logical restraints.
  - Intended to produce warning-range logical mapping coverage (`0.8`, i.e., `8/10`) when mapped against `9L1V` model atoms.

- `partial_mapping_abort.str`
  - Minimal NMR-STAR-like distance loop with 10 logical restraints.
  - Intended to produce abort-range logical mapping coverage (`0.3`, i.e., `3/10`) when mapped against `9L1V` model atoms.

## Notes

- These fixtures are diagnostics and threshold-validation helpers, not biological interpretation examples.
- Downloaded external files remain under `specs/nmr_restraints/fixtures/cache/` (gitignored).
