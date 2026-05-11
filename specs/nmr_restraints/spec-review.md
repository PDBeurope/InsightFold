# Notebook Spec Review Report

| Field | Value |
|---|---|
| Spec | `specs/nmr_restraints/` |
| Reviewer | Codex using `notebook-spec-review` |
| Review Mode | readiness |
| Date | 2026-05-08 |
| Verdict | pass-with-assumptions |

## Summary

The spec pack is complete enough to start prototype/scaffold implementation. It includes requirements, design, tasks, data contracts, fixture planning, validation, documentation, assumptions, and non-goals. Most requirements trace cleanly to notebook sections, tasks, and validation checks.

It is not ready for full scientific implementation sign-off or validation claims. The main remaining gaps are fixture readiness and exact NMR-STAR source-tag mapping. The mapping denominator, distance-bound, and embedded-rendering decisions have been resolved in the spec after review.

## Blocking Findings

| ID | Finding | Evidence | Required Change | Owner |
|---|---|---|---|---|
| B-001 | Implementation-derived fixture snapshots are still incomplete | `fixture-manifest.md` now pins `9L1V`, `1CRN`, and synthetic edge fixture roles with retrieval, checksum, structure, and parser-count expectations; mapping coverage, top violations, density summary, selected local-view residue, and runtime remain `TBD after implementation parser run` | During implementation validation, record mapping coverage, top violations, density summary, selected local-view residue, and runtime for `happy-path-9l1v`; create the synthetic wrapped-dihedral and low-mapping fixtures | Fixture curator / implementation lead |

## Major Findings

| ID | Finding | Evidence | Required Change | Owner |
|---|---|---|---|---|
| M-001 | NMR-STAR source-column mapping is underspecified for implementation | `data-contracts.md` defines normalized output columns but does not list accepted source tag names or fallback aliases for bounds, IDs, atom identifiers, saveframes, or loop categories; `tasks.md` asks to implement parsers directly | Add a parser source-tag table for distance and dihedral restraints, including required tags, accepted aliases, unit assumptions, missing-value handling, and unsupported-record behavior | Implementation lead / NMR reviewer |
| M-002 | Fixture readiness is treated as both implementation scaffolding and validation blocker | `notebook-design.md` says implementation may begin with provisional choices; `requirements.md` Q-001 is marked blocking; `fixture-manifest.md` says provisional fixtures are sufficient to start scaffolding but not final validation | Clarify implementation gate language: what may start before fixture-selection, and which tasks must wait for frozen fixture snapshots | Product owner / implementation lead |

## Resolved Review Notes

| ID | Decision | Spec Update |
|---|---|---|
| R-001 | Mapping coverage is counted at logical-restraint level. A non-ambiguous restraint maps as `1/1` when all required atoms map; an `OR` ambiguous logical restraint maps as `1/1` when at least one candidate member maps and can be evaluated. | Updated `requirements.md`, `data-contracts.md`, `notebook-design.md`, `tasks.md`, `validation.md`, and `spec-pack-overview.md`. |
| R-002 | V1 distance restraints are expected to have both upper and lower bounds. Rows missing either bound are unsupported and counted in diagnostics. | Updated `requirements.md`, `data-contracts.md`, `tasks.md`, `validation.md`, and `spec-pack-overview.md`. |
| R-003 | MolViewSpec/Mol* views must be embedded inline into notebook cells using the `state.molstar_html()` plus base64 `IFrame` pattern from `notebooks/homodimer_diagnostic.ipynb`. | Updated `requirements.md`, `notebook-design.md`, `data-contracts.md`, `tasks.md`, `validation.md`, and `spec-pack-overview.md`. |

## Minor Findings

| ID | Finding | Evidence | Suggested Change |
|---|---|---|---|
| m-001 | Violation color classes are configurable but not concretely defined | `requirements.md` Q-006 leaves low/medium/high thresholds open; `data-contracts.md` requires `display_color_class` | Add provisional color-class rules, even if marked reviewable, so implementers do not invent inconsistent thresholds. |
| m-002 | Residue identity and insertion-code handling should be made more explicit | `data-contracts.md` allows `auth_seq_id` string/int and mentions preserving insertion semantics where needed | Define a canonical string `residue_key` format including chain, sequence ID, insertion code if present, component ID, and model index. |
| m-003 | Install cell behavior is implied but not explicitly specified | `notebook-design.md` lists dependencies and install source; `tasks.md` has setup/dependency checks | Add whether the notebook should install missing packages automatically, fail with instructions, or provide a commented install cell. |
| m-004 | Export contract covers violation table but optional satisfied-restraint export remains unresolved | `requirements.md` REQ-026 requires violation CSV; Q-008 leaves satisfied restraints optional | Decide whether `evaluated_restraints` export is required for scientific reproducibility or stays a stretch task. |
| m-005 | Non-protein and ligand behavior is open but local context design includes centroid fallback | `requirements.md` Q-005 says ligand/nucleic-acid/non-protein handling is unresolved; `notebook-design.md` allows centroid fallback if enabled | Keep v1 protein-only unless fixture-selection deliberately includes non-protein restraints; otherwise mark centroid fallback as future work. |

## Human Review Questions

| ID | Question | Why Human Judgment Is Needed | Recommended Owner |
|---|---|---|---|
| H-001 | Is `9L1V` representative enough for the happy-path story, or should another NMR entry be selected? | Fixture suitability depends on deposited restraint quality, parseability, mapping behavior, and educational value | NMR-aware structural biologist / fixture curator |
| H-002 | Are the provisional thresholds `0.5 A` for distance and `5 degrees` for dihedral highlighting appropriate for an educational exploratory viewer? | Thresholds affect user interpretation and visual salience | NMR domain reviewer |
| H-003 | Are synthetic edge fixtures acceptable for ambiguity, wrapped dihedrals, and partial mapping, or must every edge fixture come from deposited entries? | Synthetic fixtures improve algorithm validation but may not represent real deposition complexity | Product owner / domain reviewer |

## Traceability Check

| Requirement ID | Design Coverage | Task Coverage | Validation Coverage | Status |
|---|---|---|---|---|
| REQ-001 | Sections 2-3 | T022, T023, T030 | Input validation | pass |
| REQ-002 | Section 3 | T023, T032 | Local-file execution, input validation | pass |
| REQ-003 | Section 4 | T031 | Retrieval contract | pass |
| REQ-004 | Section 5 | T033 | Structure parsing | pass |
| REQ-005 | Section 5 | T034 | Structure parsing | pass |
| REQ-006 | Section 6 | T035, T036 | Restraint parsing | pass-with-assumption: source tags need detail |
| REQ-007 | Section 6 | T037 | Restraint parsing, scientific language review | pass-with-assumption |
| REQ-008 | Section 7 | T040, T041 | Atom mapping | pass |
| REQ-009 | Section 7 | T042, T016 | Mapping thresholds | pass |
| REQ-010 | Section 8 | T043, T014 | Ambiguous `OR` semantics | pass-with-assumption: `edge-or-9l1v` selected; selected-member snapshot pending implementation run |
| REQ-011 | Section 8 | T044 | Distance violation formula | pass |
| REQ-012 | Section 8 | T045, T046, T015 | Dihedral formula | pass-with-assumption: `edge-dihedral-9l1v` selected for real parsing; wrapped synthetic fixture still to be created |
| REQ-013 | Section 9 | T048 | Residue density | pass |
| REQ-014 | Sections 8-9 | T047, T048 | Output snapshot checks | pass |
| REQ-015 | Section 8 | T049 | Violation prioritization | pass |
| REQ-016 | Section 10 | T060 | Visualization state build | pass |
| REQ-017 | Section 10 | T061 | Visualization state build | pass |
| REQ-018 | Sections 11-12 | T063, T065, T066, T067 | Visualization state build, selected local residue snapshot | pass-with-assumption: local residue fixture missing |
| REQ-019 | Section 12 | T068 | Visualization state build | pass |
| REQ-020 | Sections 11-12 | T063, T105 | Hidden-state hazard | pass |
| REQ-021 | Section 11 | T064, T069 | Visualization fallback | pass |
| REQ-022 | Sections 10, 12 | T062 | Visualization fallback | pass |
| REQ-023 | Sections 1, 13 | T021, T080, T107 | Scientific language review | pass |
| REQ-024 | Sections 1, 13 | T080, T103, T107 | Empty restraint and no-violation behavior | pass |
| REQ-025 | Sections 4, 14 | T031 | Retrieval contract, local-file parity | pass |
| REQ-026 | Section 14 | T082, T085 | Export artifacts | pass |

## Readiness Checklist

- [x] Requirements are observable and testable
- [x] Notebook design maps to requirements
- [x] Tasks are executable and ordered
- [ ] Fixtures are present and justified
- [x] Data contracts are explicit
- [x] Scientific assumptions are documented
- [x] Validation plan is concrete
- [x] Documentation plan covers tutorial, how-to, reference, and explanation
- [ ] Blocking questions are resolved or assigned

## Implementation Readiness

Prototype implementation may start if the team accepts these assumptions:

- Begin with scaffolding, input handling, retrieval/local loading, dependency checks, and parser spike work.
- Treat `9L1V` as the selected primary fixture; implementation-derived snapshots still need to be recorded after parser/geometry code exists.
- Do not finalize parser source-tag support or fixture-derived expected values until edge fixtures are pinned.
- Do not claim scientific validation until fixture snapshots and NMR domain review are complete.

Before full implementation sign-off, resolve B-001 and M-001. Clarify M-002 gate language before assigning implementation work broadly.
