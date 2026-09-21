# Requirement Traceability Matrix

Source PRD: `specs/nmr_restraints/nmr_restraints_prd.md`

PRD goal labels used here:

- `G-001` Make deposited NMR restraint evidence interpretable.
- `G-002` Support reproducible exploratory analysis.
- `G-003` Reduce the specialist barrier to NMR restraint interpretation.
- `G-004` Maintain scientific caution.
- `G-005` Prepare an implementation-ready notebook specification.

| Requirement | PRD Goal | Notebook Section / Cell | Fixture | Validation Check | Review Evidence |
|---|---|---|---|---|---|
| REQ-001 | G-002 | Sections 2-3 / C003-C004 | `happy-path-9l1v` | Input validation | `requirements.md` remote-mode contract |
| REQ-002 | G-002 | Sections 2-4 / C003-C005 | `local-copy-9l1v` | Local-file execution | `notebook-ux-contract.md` first runnable cell |
| REQ-003 | G-002 | Section 4 / C005 | `happy-path-9l1v`, `negative-missing-restraints-1crn` | Retrieval contract | `data-contracts.md` external source table |
| REQ-004 | G-002 | Section 5 / C006 | `happy-path-9l1v` | Structure parsing | `fixture-manifest.md` model counts |
| REQ-005 | G-002 | Section 5 / C006 | `happy-path-9l1v` | Structure parsing | `notebook-design.md` selected-model design |
| REQ-006 | G-001, G-002 | Section 6 / C007 | `happy-path-9l1v`, `edge-dihedral-9l1v` | Restraint parsing | `data-contracts.md` observed source tags |
| REQ-007 | G-002, G-004 | Section 6 / C007 | `negative-missing-restraints-1crn`, `edge-low-mapping-synthetic` | Restraint parsing, empty restraint behavior | `docs-plan.md` unsupported-record messaging |
| REQ-008 | G-001, G-002 | Section 7 / C008 | `happy-path-9l1v`, `edge-or-9l1v` | Atom mapping | `spec-pack-overview.md` preserved logical-restraint decision |
| REQ-009 | G-002, G-004 | Section 7 / C008 | `edge-low-mapping-synthetic` | Mapping thresholds | `requirements.md` threshold table |
| REQ-010 | G-001 | Section 8 / C009 | `edge-or-9l1v` | Ambiguous `OR` semantics | `notebook-ux-contract.md` trust limits |
| REQ-011 | G-001, G-004 | Section 8 / C009 | `happy-path-9l1v` | Distance violation formula | `spec-pack-overview.md` preserved bound requirement |
| REQ-012 | G-001 | Section 8 / C009 | `edge-dihedral-9l1v`, `edge-wrapped-dihedral-synthetic` | Dihedral formula | `fixture-manifest.md` wrapped-interval fixture |
| REQ-013 | G-001, G-003 | Section 9 / C009 | `happy-path-9l1v` | Residue density | `docs-plan.md` density caveat |
| REQ-014 | G-001, G-002 | Sections 8-9 / C009-C010 | `happy-path-9l1v` | Structure parsing, restraint parsing, residue density | `data-contracts.md` computed output contracts |
| REQ-015 | G-001, G-003 | Section 8 / C009 | `happy-path-9l1v` | Violation prioritization | `requirements.md` acceptance criteria |
| REQ-016 | G-001, G-003 | Section 10 / C011 | `happy-path-9l1v` | Visualization state build | `spec-pack-overview.md` inline Mol* decision |
| REQ-017 | G-001, G-003 | Section 10 / C011 | `happy-path-9l1v` | Visualization state build | `notebook-design.md` clutter-control design |
| REQ-018 | G-001, G-003 | Sections 11-12 / C012-C014 | `happy-path-9l1v` | Local residue snapshot | `notebook-ux-contract.md` user flow |
| REQ-019 | G-001, G-003 | Section 12 / C013-C014 | `edge-dihedral-9l1v` | Visualization state build | `notebook-design.md` separate local modes |
| REQ-020 | G-002 | Sections 11-12 / C012-C014 | `happy-path-9l1v` | Hidden-state hazard | `cell-blueprint.md` authoritative `selection_state` |
| REQ-021 | G-003 | Section 11 / C012 | `happy-path-9l1v` | Visualization fallback | `cell-blueprint.md` widget fallback rules |
| REQ-022 | G-002, G-003 | Sections 10-12 / C011-C014 | `happy-path-9l1v`, `local-copy-9l1v` | Visualization fallback | `validation.md` render-isolation checks |
| REQ-023 | G-004 | Sections 1, 13 / C001, C010 | all user-facing fixtures | Scientific language review | `docs-plan.md` required wording |
| REQ-024 | G-004 | Sections 1, 13 / C001, C010 | `negative-missing-restraints-1crn` | Empty restraint behavior, no-violation behavior | `requirements.md` edge-case table |
| REQ-025 | G-002 | Sections 4, 14 / C005, C015 | `local-copy-9l1v` | Retrieval contract, local-file execution | `data-contracts.md` provenance and cache contract |
| REQ-026 | G-002, G-003 | Section 14 / C015 | `happy-path-9l1v` | Export artifacts | `docs-plan.md` export documentation |

## Coverage Notes

- No major notebook output is intentionally orphaned from a requirement.
- `spec-review.md` is retained as review history but is not the primary traceability artifact anymore.
- Mutation and AlphaFold work remain intentionally absent from this matrix because they are out of scope for this notebook.
