# NMR Restraints Visualization Documentation Plan

Source PRD: `specs/nmr_restraints/nmr_restraints_prd.md`

## Documentation Goals

Documentation must make the notebook usable by structural-biology users who understand residue-level structure interpretation but may not understand NMR restraint semantics or refinement workflows. It must also give implementers and reviewers enough reference material to validate data contracts, algorithms, fixture behavior, and scientific caveats.

The notebook should stand alone as an educational research-use workflow. External docs should support reuse, validation, and maintenance.

## Documentation Types

| Documentation Type | Artifact | Audience | Required Content |
|---|---|---|---|
| Tutorial | Notebook introduction and first-run section | New users, PDB users, educators | What the notebook does, RUO scope, how to run with default `9L1V`, expected outputs, how to read density and violation views. |
| How-to | Notebook parameter/input section and optional README snippet | Returning users, researchers | How to switch PDB ID, use local files, change model index, adjust violation thresholds, change local radius, and export tables. |
| Reference | Spec pack plus notebook appendix | Developers, maintainers, power users | Input parameters, endpoints, file formats, output table schemas, cache key, formulas, mapping thresholds, dependency versions. |
| Explanation | Notebook markdown sections | All users | NMR-derived structures, distance restraints, dihedral restraints, ambiguity, restraint density, violations, interpretation risks, and v1 limitations. |
| Validation notes | `validation.md` and notebook validation snapshot | Implementers and reviewers | Fixture expected values, runtime, dependency checks, algorithm checks, edge-case behavior, and known unresolved questions. |

## Notebook Markdown Requirements

The notebook must include concise explanatory markdown before each major computational output:

- Purpose and research-use scope.
- What NMR restraints are and why they matter.
- Difference between structure coordinates and experimental restraints.
- Distance restraint interpretation.
- Dihedral restraint interpretation.
- Ambiguous `OR` restraint behavior in v1.
- Atom mapping and why unmapped restraints matter.
- Restraint density as experimental coverage, not confidence.
- Violation magnitude as local inconsistency, not automatic error.
- No-violation caveat.
- Missing-restraint and unsupported-record caveats.
- Local evidence view interpretation.
- Export and provenance interpretation.

## Required User-Facing Copy

The implementation may adapt wording, but must preserve these meanings.

### Notebook Introduction

This notebook visualizes how deposited NMR restraints relate to a deposited structural model. It helps users explore where restraint support is dense, sparse, or locally inconsistent. It is an exploratory evidence viewer, not a formal validation engine.

### Density Caveat

Restraint density is the number of deposited restraints associated with a residue. It is a proxy for experimental coverage and should not be interpreted as structural confidence or correctness.

### Violation Caveat

A restraint violation occurs when measured model geometry falls outside the deposited allowed range. Violations may reflect ambiguity, conformational heterogeneity, local dynamics, refinement tradeoffs, or model-restraint inconsistency. They are not automatically errors.

### No-Violation Caveat

No thresholded violations were detected with the current settings. This does not mean the structure is fully validated, correct, or free of experimental limitations.

### Mapping Warning

Some restraint atoms could not be mapped to the selected coordinate model. Results may represent only a subset of deposited restraints.

### Unsupported Ambiguity Warning

Some ambiguous restraint semantics are not supported in v1 and were excluded. V1 supports explicit `OR` distance restraints by selecting the smallest measured candidate distance.

### Empty Restraint State

No compatible NMR restraint data were available for this structure or input file. The notebook cannot compute restraint density or violations without compatible restraints.

## External Documentation Updates

| Artifact | Update Needed | Timing |
|---|---|---|
| Notebook README or index | Add notebook purpose, path, default fixture, dependency notes, and RUO scope | After notebook implementation path is final |
| Fixture manifest | Add observed expected snapshots and selected local-view residue | Before validation sign-off |
| Data contracts | Update if real fixtures expose schema changes or parser-specific column names | During implementation |
| Validation report | Add executed fixture results, runtime, warnings, and known limitations | After execution validation |
| Lifecycle runbook, if present | Add notebook stage, review status, and graduation evidence | Before final review |

## Export Documentation

Exported tables or summaries should include:

- source PDB ID or local label
- model index
- source model file URL/path
- source restraint file URL/path
- retrieval/load timestamp
- cache status
- configuration thresholds
- mapping coverage
- warning count
- RUO/non-validation caveat

## Maintainer Notes

- Keep examples tied to pinned fixtures, not arbitrary current API behavior.
- Do not describe density values as confidence scores.
- Do not describe violation counts as pass/fail validation.
- Keep ambiguous restraint limitations visible wherever ambiguous results are shown.
- Treat parser changes as scientific changes unless fixture snapshots prove equivalent behavior.

## Documentation Acceptance Criteria

- A new user can run the default fixture without reading the PRD.
- A returning user can identify where to change PDB ID, local file paths, thresholds, radius, and model index.
- A developer can find every output table schema in `data-contracts.md` or notebook reference cells.
- A reviewer can trace formulas and assumptions from notebook text to `requirements.md`, `notebook-design.md`, and `validation.md`.
- The notebook includes the required caveats for density, violations, no-violation states, missing restraints, mapping warnings, and unsupported ambiguity.
