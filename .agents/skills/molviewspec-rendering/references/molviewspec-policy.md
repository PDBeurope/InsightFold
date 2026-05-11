# MolViewSpec Policy

InsightFold notebooks must use MolViewSpec/Mol* for 3D protein structure visualization.

## Required Pattern

- Generate declarative MolViewSpec state with the `molviewspec` Python package.
- Prefer AFDB `bcifUrl` for remote AlphaFold Database structures when available.
- Render inline with `state.molstar_html()` and a base64 `IFrame` data URI when supported.
- Provide `.mvsj` export or Mol* viewer URL fallback when inline rendering is unavailable.
- Keep visualization failures isolated from data loading, scoring, validation, and interpretation.

## Forbidden Substitute Pattern

Do not add imperative notebook-only 3D protein viewers as an alternative backend. If encountered in legacy notebooks or helper modules, mark it for migration to MolViewSpec/Mol* rather than copying the pattern into new work.

## Migration Notes

- Move structure source, chain IDs, residue numbering, confidence arrays, and interface masks into explicit inputs.
- Replace imperative styling calls with MolViewSpec builder state generation.
- Preserve working inline rendering behavior from `notebooks/homodimer_diagnostic.ipynb`: `state.molstar_html()` -> base64 encode -> `IPython.display.IFrame`.
- Add a generated/skipped view summary so validation can inspect behavior.

