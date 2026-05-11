# Notebook Spec Review Report

| Field | Value |
|---|---|
| Spec | `specs/protein_model_chem/` |
| Reviewer | Codex |
| Review Mode | readiness |
| Date | 2026-05-08 |
| Verdict | pass-with-assumptions |

## Summary

The spec pack is complete as a document set and has good workflow coverage, data contracts, task decomposition, RUO framing, validation intent, and pinned prototype fixtures.

Post-fixture update: the manifest now pins representative PDBe/PDB, AFDB, interface, unsupported-context, invalid-identifier, mutation-mismatch, and optional missing-residue fixtures. Prototype implementation may start with explicit assumptions.

Remaining assumptions: exact interaction counts, perturbation counts, graph deltas, and severity labels must be frozen after the first trusted notebook run and reviewed by a structural biology domain reviewer.

## Blocking Findings

| ID | Finding | Evidence | Required Change | Owner |
|---|---|---|---|---|
| B-001 | None active for prototype implementation. | Minimum fixtures are now pinned in `fixture-manifest.md`. | Freeze exact interaction and perturbation snapshots after first trusted implementation run. | Implementation lead / domain reviewer |

## Major Findings

| ID | Finding | Evidence | Required Change | Owner |
|---|---|---|---|---|
| M-001 | Dependency policy now follows the revised visualization decision but still needs dependency-risk review. | Spec now requires BioPython, SciPy, NetworkX, Plotly, Requests, and Python 3.11+; MolViewSpec is optional visualization-only code with required graceful fallback. | Confirm this dependency budget is acceptable before notebook implementation. | Implementation lead |
| M-002 | Interaction weighting is specified but still scientifically unreviewed. | Spec now includes baseline weights and relative-weight factors. | Domain-review the weights and ensure outputs label them heuristic/non-thermodynamic. | Domain reviewer |
| M-003 | Severity ranking now has a provisional rubric but remains uncalibrated. | Spec defines low/moderate/high descriptive rules from clash, buried mismatch, interaction loss, cavity, centrality, fragmentation, and secondary-structure incompatibility. | Validate against fixtures and domain-review before scientific sign-off. | Domain reviewer |
| M-004 | PDBe source handling has a provisional policy, but exact suitability cutoffs remain undefined. | Spec now prefers PDBe when resolution and residue coverage are acceptable, and builds biological assembly if available. | Domain-review resolution and residue-coverage cutoffs during fixture selection. | Implementation lead / domain reviewer |
| M-005 | Validation plan has pinned fixture metadata but not exact computed interaction snapshots. | Fixture manifest includes stable accessions and behavioral expectations; exact interaction counts and perturbation counts remain implementation-derived. | After first trusted implementation run, add exact/tolerance snapshots for local WT interaction count, interaction types, perturbation signal counts, graph deltas, and severity rationale. | Implementation lead / domain reviewer |
| M-006 | Export contract now names required formats but still needs fixture validation. | Spec now requires JSON, CSV, PNG/SVG, and graph edge/node or graph object export. | Validate generated artifacts after implementation. | Implementation lead |

## Minor Findings

| ID | Finding | Evidence | Suggested Change |
|---|---|---|---|
| m-001 | Runtime target now includes full-notebook and subruntime targets. | Spec has under 2 minutes full run, under 30 seconds structure analysis, and under 10 seconds mutation recomputation. | Validate with fixtures. |
| m-002 | Mutation input formats are now enumerated in data contracts. | Spec includes `R273H`, `A:W128A`, and `CHAIN_A:W128A`. | Add examples to notebook markdown during implementation. |
| m-003 | Visualization color conventions are now specified. | Spec defines colors for hydrogen bond, salt bridge, hydrophobic, aromatic, disulfide, and steric clash. | Validate rendered colors where practical. |
| m-004 | Logging/observability now has a notebook-local log artifact. | Spec includes `notebook_log`. | Ensure all modules append relevant warnings/failures during implementation. |
| m-005 | Recommended dependencies include `mdtraj` and `freesasa`, but V1 scope is meant to stay lightweight. | Supporting doc recommends these, while spec disallows heavyweight simulation and treats packing/SASA as optional. | Keep `freesasa` optional for SASA/burial only and reject `mdtraj` unless a concrete lightweight use is justified. |

## Human Review Questions

| ID | Question | Why Human Judgment Is Needed | Recommended Owner |
|---|---|---|---|
| H-001 | Are the supporting-doc interaction thresholds scientifically acceptable as provisional v1 defaults? | H-bond, salt bridge, hydrophobic, aromatic, disulfide, and clash cutoffs determine all downstream results. | Structural biology domain reviewer |
| H-002 | Are baseline interaction weights acceptable as heuristic relative constraint weights? | The weights look quantitative and may be overinterpreted unless reviewed and clearly labeled. | Structural biology domain reviewer |
| H-003 | Should V1 use BioPython as the parser, or is another parser preferred for mmCIF/PDB robustness? | Parser choice affects residue numbering, alternate conformations, biological assemblies, and missing residue handling. | Implementation lead with domain input |
| H-004 | Should local sidechain substitution be represented with atom templates, residue-class rules, or interaction removal/addition heuristics? | This determines whether mutation outputs are geometry-derived or rule-derived and affects scientific defensibility. | Domain reviewer / implementation lead |
| H-005 | What minimum PDBe resolution and residue coverage qualify as "suitable experimental structure" for PDBe-first behavior? | Experimental structures are not uniformly reliable; the source preference needs explicit criteria. | Structural biology domain reviewer |
| H-006 | Can low/moderate/high structural disruption likelihood use provisional heuristic bands, or should severity remain descriptive until validated fixtures exist? | Severity labels are useful but high-risk for overinterpretation. | Product owner / domain reviewer |
| H-007 | Which fixtures should represent happy path, low-confidence AFDB, interface perturbation, missing biological context, invalid ID, mutation mismatch, and missing residue? | Fixture choices define the scientific and implementation test surface. | Fixture selection / domain reviewer |

## Traceability Check

| Requirement ID | Design Coverage | Task Coverage | Validation Coverage | Status |
|---|---|---|---|---|
| REQ-001 | Sections 2-3, `structure_request` | T010 | Retrieval and top-to-bottom checks | pass |
| REQ-002 | Section 4, AFDB data contract | T011, T015 | Structure retrieval contract | pass |
| REQ-003 | Section 4, PDBe data contract | T012, T015 | Structure retrieval contract | pass |
| REQ-004 | Section 4, source preference | T013 | Needs PDBe suitability fixture/criteria | pass-with-human-review |
| REQ-005 | Section 5, residue/atom tables | T020-T024 | Structure parsing contract | pass-with-fixtures |
| REQ-006 | Section 6, mutation validation | T025 | Mutation validation, FX-006 | pass-with-fixtures |
| REQ-007 | Section 7, interaction config, function candidates | T030-T038 | Interaction detection and configuration audit | pass-with-human-review |
| REQ-008 | Section 7, graph contract | T040 | Graph construction | pass |
| REQ-009 | Section 8, network metrics and weighted graph design | T041 | Graph construction and interaction weighting | pass-with-human-review |
| REQ-010 | Sections 7 and 11 | T061, T063 | Visualization existence | pass |
| REQ-011 | Section 9, sidechain-only 8 Angstrom mutation model | T050-T056 | Mutation perturbation output | pass-with-human-review |
| REQ-012 | Section 10 | T054-T056, T064 | Mutation perturbation output | pass-with-fixtures |
| REQ-013 | Section 12 | T067, T072 | Explanation guardrails | pass |
| REQ-014 | Section 12, provisional severity rubric | T068-T069 | Severity validation | pass-with-human-review |
| REQ-015 | Sections 5 and 8 | T042-T043 | Confidence propagation | pass-with-fixtures |
| REQ-016 | Error handling, confidence caveats | T043 | Unsupported-context warnings | pass-with-fixtures |
| REQ-017 | Section 11 visualization plan | T060-T066 | Visualization existence | pass |
| REQ-018 | Section 13, export contract | T070 | Export completeness and export-format validation | pass |
| REQ-019 | Section 1 and guardrails | T072, T087 | Explanation guardrails | pass |

## Remaining Required Edits

Resolved in this pass:

- `requirements.md` now converts parser, threshold, mutation model, weighting, and visualization questions into provisional decisions or human-review questions.
- `notebook-design.md` now includes Python 3.11+, BioPython, MolViewSpec/Mol*, Plotly, NetworkX, configuration defaults, thresholds, weights, severity rubric, colors, and export formats.
- `data-contracts.md` now includes mutation notation parsing, weighting fields, perturbation schema fields, interpretation record fields, and export formats.
- `validation.md` now includes configuration, weighting, export-format, and subruntime validation checks.
- `fixture-manifest.md` now pins prototype fixtures FX-001 through FX-007.

Still required before scientific sign-off:

1. During implementation and domain review:
   - Freeze exact local interaction counts and perturbation counts.
   - Confirm chain/residue numbering for `6M0J` under BioPython.
   - Add cache paths and checksums if fixture caching is introduced.

2. During domain review:
   - Confirm interaction thresholds and baseline weights.
   - Confirm BioPython parser behavior on selected fixtures.
   - Confirm PDBe suitability criteria.
   - Confirm provisional severity labels or downgrade severity to descriptive-only output.

## Readiness Checklist

- [x] Requirements are observable and testable
- [x] Notebook design maps to requirements
- [x] Tasks are executable and ordered
- [x] Fixtures are present and justified
- [x] Data contracts are explicit
- [x] Scientific assumptions are documented with concrete thresholds and weights in the spec pack
- [x] Validation plan is concrete enough for prototype run against pinned fixtures
- [x] Documentation plan covers tutorial, how-to, reference, and explanation
- [x] Non-fixture blocking questions are resolved or assigned

## Implementation Readiness Decision

Prototype implementation may start with assumptions. The next action is `$notebook-from-spec` after acknowledging that exact computed snapshots and domain-reviewed severity labels are still pending.

Full pass still requires domain review of thresholds, weights, mutation interpretation, PDBe suitability, severity labels, and implementation-derived fixture snapshots.
