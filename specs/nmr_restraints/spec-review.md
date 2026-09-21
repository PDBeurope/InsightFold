# Notebook Spec Review Report

| Field | Value |
|---|---|
| Spec | `specs/nmr_restraints/` |
| Reviewer | Codex using `notebook-spec-review` |
| Review Mode | readiness |
| Date | 2026-05-11 |
| Verdict | pass-with-assumptions |

## Summary

The upgraded pack now meets the current README lifecycle artifact contract for a notebook spec pack. It includes the required UX contract, cell blueprint, traceability matrix, fixture manifest, data contracts, validation plan, and tasks, and those artifacts are materially connected rather than placeholder-only.

The first runnable user input cell is explicit in `notebook-ux-contract.md` and is carried into `cell-blueprint.md` as C003. The pack also now cleanly separates the real user workflow from validation fixtures, makes the MolViewSpec/Mol* route explicit as inline embedded rendering rather than external-link handoff, and provides requirement-to-cell-to-validation traceability across the full requirement set.

Implementation may start for notebook scaffolding and core build work. It should start with the documented assumptions still in place: fixture-derived expected snapshots are not frozen yet, `pynmrstar` remains provisional until confirmed on the pinned fixtures, and NMR domain review is still required before scientific/public readiness or validation sign-off.

## Blocking Findings

No blocking findings for implementation-start readiness under the current lifecycle contract.

## Major Findings

| ID | Finding | Evidence | Required Change | Owner |
|---|---|---|---|---|
| M-001 | Validation readiness still depends on frozen fixture-derived snapshots and synthetic edge fixtures, so the pack is ready for implementation start but not execution validation sign-off | `fixture-manifest.md` leaves `happy-path-9l1v` mapping coverage, top violations, density summary, selected local-view residue, and runtime as implementation-derived snapshots; `validation.md` still lists snapshot artifacts as required evidence; `tasks.md` keeps this work in T001, T012, T015, and T016 | Freeze `happy-path-9l1v` expected outputs and create the wrapped-dihedral and low-mapping synthetic fixtures before execution validation is treated as complete | Fixture curator / implementation lead |
| M-002 | The parser contract is much stronger than before, but `pynmrstar` is still a provisional dependency rather than a confirmed fixture-backed parser choice | `spec-pack-overview.md` lists `pynmrstar` suitability as a blocking item before scientific sign-off; `notebook-design.md` and `requirements.md` both call it provisional; `tasks.md` T003 still asks for confirmation or replacement/wrapper documentation | Confirm `pynmrstar` on the pinned `9L1V` fixture set and document any wrapper or replacement decision before validation claims are made | Implementation lead / domain reviewer |
| M-003 | NMR domain review remains an explicit prerequisite for scientific wording, threshold interpretation, and public educational release | `spec-pack-overview.md` lists NMR domain review as required; `requirements.md` keeps Q-006 and Q-009 open at the review/owner level; `docs-plan.md` and `validation.md` both require cautious interpretation language review | Define and complete the NMR domain-review checklist before release or scientific sign-off; keep implementation scoped to RUO exploratory behavior until then | Domain reviewer / product owner |

## Minor Findings

| ID | Finding | Evidence | Suggested Change |
|---|---|---|---|
| m-001 | The review language around “implementation may start” and “validation may not start” is mostly clear, but it appears in several files with slightly different emphasis | `spec-pack-overview.md`, `requirements.md`, `fixture-manifest.md`, and `tasks.md` all describe the implementation-versus-validation boundary | Add one short canonical readiness note in `spec-pack-overview.md` and mirror that wording elsewhere to reduce future gate ambiguity |
| m-002 | The export surface is still slightly open-ended between required and optional advanced artifacts | `requirements.md` makes violation CSV required while `tasks.md` and `docs-plan.md` leave some extra exports optional | Keep violation CSV as the only required export in v1 and explicitly label other exports as optional stretch outputs in one place |

## Human Review Questions

| ID | Question | Why Human Judgment Is Needed | Recommended Owner |
|---|---|---|---|
| H-001 | Are the provisional distance and dihedral display thresholds appropriate for an exploratory educational notebook? | Thresholds affect user interpretation and visual salience even when they are not validation rules | NMR domain reviewer |
| H-002 | Is `9L1V` still the best primary fixture after parser confirmation, or should the happy-path fixture change before fixture freezing? | The right happy-path fixture depends on deposited restraint quality, parseability, mapping behavior, and teaching value | Fixture curator / NMR-aware structural biologist |
| H-003 | Is the v1 ambiguity explanation sufficiently cautious for public-facing educational use? | Smallest-distance `OR` selection is a pragmatic simplification and may need domain-approved wording | Domain reviewer / product owner |

## Traceability Check

| Requirement ID | Design Coverage | Task Coverage | Validation Coverage | Status |
|---|---|---|---|---|
| REQ-001 | Sections 2-3 / C003-C004 | T022, T023, T030 | Input validation | pass |
| REQ-002 | Sections 2-4 / C003-C005 | T023, T032 | Local-file execution | pass |
| REQ-003 | Section 4 / C005 | T031 | Retrieval contract | pass |
| REQ-004 | Section 5 / C006 | T033 | Structure parsing | pass |
| REQ-005 | Section 5 / C006 | T034 | Structure parsing | pass |
| REQ-006 | Section 6 / C007 | T035, T036 | Restraint parsing | pass-with-assumption: parser choice still provisional |
| REQ-007 | Section 6 / C007 | T037 | Restraint parsing, empty restraint behavior | pass |
| REQ-008 | Section 7 / C008 | T040, T041 | Atom mapping | pass |
| REQ-009 | Section 7 / C008 | T041, T042 | Mapping thresholds | pass |
| REQ-010 | Section 8 / C009 | T043 | Ambiguous `OR` semantics | pass |
| REQ-011 | Section 8 / C009 | T044 | Distance violation formula | pass |
| REQ-012 | Section 8 / C009 | T045, T046 | Dihedral formula | pass-with-assumption: wrapped synthetic fixture still needs creation |
| REQ-013 | Section 9 / C009 | T048 | Residue density | pass |
| REQ-014 | Sections 8-9 / C009-C010 | T047, T048 | Computed output tables | pass |
| REQ-015 | Section 8 / C009 | T047, T049 | Violation prioritization | pass |
| REQ-016 | Section 10 / C011 | T060, T062 | Visualization state build | pass |
| REQ-017 | Section 10 / C011 | T061, T062 | Visualization state build | pass |
| REQ-018 | Sections 11-12 / C012-C014 | T063, T065, T066, T067 | Local residue snapshot | pass-with-assumption: frozen local snapshot still pending |
| REQ-019 | Section 12 / C013-C014 | T068 | Visualization state build | pass |
| REQ-020 | Sections 11-12 / C012-C014 | T063, T105 | Hidden-state hazard | pass |
| REQ-021 | Section 11 / C012 | T064, T069 | Visualization fallback | pass |
| REQ-022 | Sections 10-12 / C011-C014 | T062, T104 | Visualization fallback | pass |
| REQ-023 | Sections 1, 13 / C001, C010 | T021, T080, T107 | Scientific language review | pass-with-assumption: final wording still needs domain review |
| REQ-024 | Sections 1, 13 / C001, C010 | T080, T103, T107 | Empty restraint behavior, no-violation behavior | pass |
| REQ-025 | Sections 4, 14 / C005, C015 | T031, T082-T084 | Retrieval contract, cache/provenance behavior | pass |
| REQ-026 | Section 14 / C015 | T082, T085 | Export artifacts | pass |

## Readiness Checklist

- [x] Requirements are observable and testable
- [x] Notebook design maps to requirements
- [x] Tasks are executable and ordered
- [x] Fixtures are present and justified
- [x] Data contracts are explicit
- [x] Scientific assumptions are documented
- [x] Validation plan is concrete
- [x] Documentation plan covers tutorial, how-to, reference, and explanation
- [x] Blocking questions are resolved or assigned

## Implementation Readiness

Implementation may start: yes, for notebook scaffolding and core implementation work.

Implementation may not yet claim: fixture-complete execution validation, parser finalization, or scientific/public readiness.

The next lifecycle expectations are:

1. use fixture-curation work to freeze expected snapshots and synthetic edge fixtures
2. confirm or replace the provisional `pynmrstar` parser decision on pinned fixtures
3. complete NMR domain review before release-facing interpretation or scientific sign-off

Notebook build was not started as part of this review.
