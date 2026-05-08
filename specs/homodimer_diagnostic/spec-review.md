# Notebook Spec Review Report

| Field | Value |
|---|---|
| Spec | `specs/homodimer_diagnostic/` |
| Reviewer | `agents/lifecycle/spec-reviewer.md` using `skills/notebook-spec-review/SKILL.md` |
| Review Mode | readiness |
| Date | 2026-05-07 |
| Verdict | pass-with-assumptions |

## Summary

The spec pack is good enough to start implementation of a prototype notebook because requirements, notebook design, data contracts, validation checks, documentation intent, and task breakdown are present and traceable.

It is not ready to claim full validation, beta readiness, or scientific sign-off. The main blockers are fixture curation, frozen numeric reference outputs, explicit score interpretation bands, and an IPSAE reference version/commit.

## Blocking Findings

| ID | Finding | Evidence | Required Change | Owner |
|---|---|---|---|---|
| B-001 | Full validation fixture set is incomplete | `fixture-manifest.md` marks FX-002 borderline, FX-003 disagreement, and FX-005 unsupported valid AFDB input as TBD/needs curation | Select accessions and record provenance, expected outputs, and validation purpose | Human/domain reviewer with fixture curator |
| B-002 | Reference numeric outputs are not frozen | `fixture-manifest.md` says FX-001 numeric expected values are TBD; `validation.md` requires reference comparison | Generate trusted score snapshots and record tolerances, especially for ipSAE-family values | Notebook validator / domain reviewer |
| B-003 | Diagnostic threshold bands are unresolved | `requirements.md` OQ-003 is blocking for final user-facing interpretation thresholds | Define high/moderate/low bands per metric or mark summary language as qualitative only | Human/domain reviewer |

## Major Findings

| ID | Finding | Evidence | Required Change | Owner |
|---|---|---|---|---|
| M-001 | IPSAE reference version is not pinned | `validation.md` says IPSAE v4 / January 2026 but also requires exact version or commit | Record the exact DunbrackLab/IPSAE commit/tag or archived source used for validation | Notebook validator |
| M-002 | Happy-path accession is only a candidate | `fixture-manifest.md` marks `AF-0000000065889468` as candidate | Confirm it is publicly retrievable, two-chain homodimeric, and suitable for the intended story | Fixture curator |
| M-003 | MolViewSpec fallback is specified but not acceptance-tested | `notebook-design.md` allows fallback; `validation.md` checks minimum views or fallback | Define what counts as acceptable fallback output in Colab and local notebook contexts | Notebook builder |

## Minor Findings

| ID | Finding | Evidence | Suggested Change |
|---|---|---|---|
| m-001 | Requirements are comprehensive but dense | `requirements.md` contains 21 requirements | During implementation, keep markdown explanations short and user-facing rather than mirroring requirement text |
| m-002 | Disallowed dependencies are clear but allowed dependency versions are not | `notebook-design.md` lists package names only | Record versions during validation and consider a lightweight requirements cell or environment note |
| m-003 | pLDDT fallback policy is an assumption | `requirements.md` OQ-002 | Confirm with domain owner before final review |

## Human Review Questions

| ID | Question | Why Human Judgment Is Needed | Recommended Owner |
|---|---|---|---|
| H-001 | Which accessions best represent borderline and disagreement cases? | Fixture representativeness depends on biological and AFDB production knowledge | AFDB/PDBe domain reviewer |
| H-002 | What score bands should be shown to non-specialists? | Thresholds influence interpretation and risk overclaiming | AFDB/PDBe domain reviewer |
| H-003 | Should clinical/RUO wording be stronger for public-facing use? | Audience includes clinical researchers, but notebook is RUO only | Product/domain reviewer |

## Traceability Check

| Requirement ID | Design Coverage | Task Coverage | Validation Coverage | Status |
|---|---|---|---|---|
| REQ-001 | Sections 2, 4 | T013, T014 | VAL-003 | pass |
| REQ-002 | Sections 3, 4 | T006 | VAL-002 | pass |
| REQ-003 | Sections 3, 4 | T007 | VAL-002 | pass |
| REQ-004 | Section 5 | T010 | VAL-005 | pass |
| REQ-005 | Sections 4, 5 | T008 | VAL-004 | pass |
| REQ-006 | Section 5 | T009 | VAL-006 | pass |
| REQ-007 | Section 5 | T011 | VAL-003 | pass |
| REQ-008 | Section 6 | T016, T017 | VAL-007 | pass |
| REQ-009 | Section 8 | T018, T019 | VAL-008 | pass |
| REQ-010 | Section 8 | T018, T020 | VAL-009 | pass |
| REQ-011 | Section 8 | T021 | VAL-010 | pass |
| REQ-012 | Section 8 | T022 | VAL-010 | pass |
| REQ-013 | Section 8 | T023 | VAL-011 | pass |
| REQ-014 | Section 8 | T024 | VAL-008 to VAL-011 | pass |
| REQ-015 | Section 9 | T025, T026 | VAL-007 | pass |
| REQ-016 | Section 7 | T027, T028 | VAL-012 | pass |
| REQ-017 | Section 9 | T029, T030 | VAL-013 | pass |
| REQ-018 | Section 10 | T031 | VAL-014 | pass-with-assumption |
| REQ-019 | Section 11 | T032 | VAL-015 | pass-with-assumption |
| REQ-020 | Sections 2, 12 | T013, T033 to T038 | VAL-001, VAL-016, VAL-017 | pass |
| REQ-021 | Sections 1, 11 | T012, T032 | VAL-015 | pass |

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

- Use FX-001 as a smoke-test accession while fixture curation continues.
- Treat score-band language as qualitative until thresholds are approved.
- Do not claim reference validation until numeric snapshots and IPSAE commit/version are pinned.

Before beta or graduation, B-001 through B-003 must be resolved.

