# Homodimer Diagnostic Notebook Tasks

## Phase 1: Spec and Fixture Readiness

- [x] T001 Convert PRD into directory-style notebook spec pack.
- [ ] T002 Resolve fixture blockers for borderline and disagreement accessions. FX-002 is a provisional lower-confidence candidate; FX-003 still needs AFDB/PDBe domain input.
- [ ] T003 Record numeric expected outputs for FX-001 after trusted reference run. Provisional smoke snapshots are recorded in `fixture-manifest.md` and `validation.md`.
- [ ] T004 Decide high/moderate/low threshold bands for diagnostic summary.
- [ ] T005 Confirm IPSAE reference version or commit for validation.

## Phase 2: Data Contracts

- [x] T006 Implement AFDB metadata fetcher and response validation.
- [x] T007 Implement download helpers for text and JSON resources.
- [x] T008 Implement PAE JSON parser and schema checks.
- [x] T009 Implement confidence JSON parser and pLDDT validation.
- [x] T010 Implement lightweight mmCIF atom-site parser.
- [x] T011 Implement homodimer compatibility validation.

## Phase 3: Notebook Architecture

- [x] T012 Create notebook introduction, RUO framing, and audience guidance.
- [x] T013 Create setup/import section with allowed dependencies only.
- [x] T014 Create accession input and provenance display section.
- [x] T015 Create validation summary section for fetched data.

## Phase 4: Scientific Computation

- [x] T016 Implement residue table construction with CA and CB/CA fallback coordinates.
- [x] T017 Implement CB-CB inter-chain contact detection with 8.0 Angstrom cutoff.
- [x] T018 Implement `d0` and `ptm` shared functions.
- [x] T019 Implement ipTM_d0chn calculations and intermediate table.
- [x] T020 Implement ipSAE_d0res, ipSAE_d0chn, and ipSAE_d0dom calculations.
- [x] T021 Implement pDockQ calculation.
- [x] T022 Implement pDockQ2 calculation.
- [x] T023 Implement LIS calculation.
- [x] T024 Build final score table.

## Phase 5: Visualizations and Interpretation

- [x] T025 Implement contact map visualization.
- [x] T026 Implement interface coverage visualization.
- [x] T027 Implement full PAE heatmap with quadrant labels.
- [x] T028 Implement score-specific PAE masks.
- [x] T029 Implement per-residue contribution profiles.
- [x] T030 Implement interface versus non-interface pLDDT comparison.
- [x] T031 Implement minimum MolViewSpec chain overview and pLDDT mapping views.
- [x] T032 Implement plain-language diagnostic summary.

## Phase 6: Notebook Validation

- [ ] T033 Run static notebook validation.
- [ ] T034 Run smoke execution with FX-001.
- [ ] T035 Run negative fixture validation with malformed accession.
- [ ] T036 Run full validation once FX-002, FX-003, and FX-005 are curated.
- [ ] T037 Compare IPSAE-family values with reference implementation within +/- 0.001.
- [ ] T038 Record validation report in `specs/homodimer_diagnostic/validation-report.md`.

## Phase 7: Final Review

- [ ] T039 Run notebook review after execution validation.
- [ ] T040 Resolve blocking or major review findings.
- [ ] T041 Decide whether notebook is ready for beta, continued prototype iteration, or archive.
