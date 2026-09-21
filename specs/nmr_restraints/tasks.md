# NMR Restraints Visualization Notebook Tasks

Source requirements: `specs/nmr_restraints/requirements.md`

Current-contract note: this task list assumes the upgraded spec-pack artifacts exist before notebook implementation starts. It does not authorize notebook construction in this lifecycle step.

## Phase A: Spec Preparation

- [x] T001 Resolve blocking fixture snapshot questions in `fixture-manifest.md` and record expected snapshots for the preserved fixture set.
- [ ] T002 Confirm whether PDBe entry-file endpoints remain canonical for v1 remote retrieval and document any fallback policy changes.
- [x] T003 Confirm `pynmrstar` parser suitability on the selected fixture set; document replacement/wrapper decision if needed.
- [ ] T004 Define domain-review checklist for NMR restraint interpretation language before public educational release.

## Phase B: Data Contracts And Fixtures

- [x] T010 Create or identify fixture data directory for NMR restraint notebook artifacts.
- [x] T011 Fetch/cache `9L1V` mmCIF and NMR-STAR files or document non-committed retrieval instructions, keeping raw downloaded files under `specs/nmr_restraints/fixtures/cache/` or another repo-local gitignored cache path.
- [x] T012 Generate implementation-derived expected snapshots for `happy-path-9l1v`: mapping coverage, density summary, top violations, selected local-view residue, and runtime.
- [ ] T013 [P] Validate negative missing-restraints fixture `negative-missing-restraints-1crn` using model endpoint `200` and restraint endpoint `404`.
- [ ] T014 [P] Validate ambiguous `OR` fixture `edge-or-9l1v` using 451 logical `OR` groups from `9L1V`.
- [x] T015 [P] Create `edge-wrapped-dihedral-synthetic` with hand-computed circular interval expectations.
- [x] T016 [P] Create `edge-low-mapping-synthetic` fixtures for warning and abort thresholds.
- [ ] T017 Update `data-contracts.md` if fixture parsing exposes required field or schema changes beyond the observed `9L1V` source tags.

## Phase C: Notebook Architecture

- [x] T020 Create notebook skeleton at `notebooks/nmr_restraints_visualization.ipynb`.
- [x] T021 Add purpose, RUO framing, and interpretation caveats markdown.
- [x] T022 Add setup cell with imports, dependency checks, version reporting, and `CONFIG`.
- [x] T023 Add input-mode section supporting remote PDB ID and local file paths.
- [x] T024 Add central `notebook_log` and structured status objects for failures/warnings.
- [x] T025 Add validation snapshot section scaffold with placeholder checks linked to fixtures.

## Phase D: Retrieval And Parsing Implementation

- [x] T030 Implement `validate_pdb_id(value)` and input validation tests in notebook cells.
- [x] T031 Implement `fetch_with_retries(url, cache_key, config)` with timeout, retry, cache, and provenance.
- [ ] T032 Implement `load_local_file(path, expected_suffixes)` for local model/restraint mode.
- [x] T033 Implement Gemmi structure parser producing `atom_table`, `residue_table`, `atom_lookup`, and `model_summary`.
- [x] T034 Implement model selection and v1 altloc filtering with ignored-altloc diagnostics.
- [x] T035 Implement NMR-STAR parser for compatible distance restraint loops.
- [x] T036 Implement NMR-STAR parser for compatible dihedral restraint loops.
- [x] T037 Implement unsupported-record collection with source saveframe/loop provenance.

## Phase E: Mapping And Scientific Computation

- [x] T040 Implement canonical atom-key normalization and exact `auth_*` atom mapping.
- [x] T041 Implement `mapping_report` with parsed, mapped, unmapped, excluded, coverage, and threshold status at logical-restraint level; report raw ambiguous member counts separately.
- [x] T042 Implement abort/warn/normal branch behavior for mapping coverage thresholds.
- [x] T043 Implement ambiguous `OR` distance grouping by logical restraint ID and smallest-distance selection.
- [x] T044 Implement Euclidean distance measurement and distance violation calculation requiring both lower and upper bounds for evaluated v1 distance restraints; exclude missing-bound rows with diagnostics.
- [x] T045 Implement torsion-angle calculation normalized to `[-180, +180]`.
- [x] T046 Implement circular interval and wrapped dihedral violation calculation.
- [x] T047 Implement `evaluated_restraints` and `violation_table` output tables.
- [x] T048 Implement residue absolute and normalized density calculation.
- [x] T049 Add configurable distance and dihedral display thresholds and color-class assignment.

## Phase F: Notebook UI And Visualization

- [x] T060 Implement global density-colored MolViewSpec state builder.
- [x] T061 Implement global violated-restraint overlay with `max_visible_restraints` cap.
- [x] T062 Implement `render_mvs_state(state, label)` with isolated rendering failures and fallback messages, using `state.molstar_html()` embedded inline as a base64 data-URI `IFrame` as in `notebooks/homodimer_diagnostic.ipynb`; remote/cache MolViewSpec states must use the public PDBe model URL rather than embedding the full cached mmCIF as a nested data URI.
- [x] T063 Implement residue selector and authoritative `selection_state`.
- [ ] T064 [P] Implement sequence/density navigation view with non-color encoding where feasible.
- [x] T065 Implement local neighborhood computation using CA coordinates for proteins and documented centroid fallback if enabled.
- [x] T066 Implement local evidence table for selected residue, including satisfied, violated, and ambiguous restraints.
- [x] T067 Implement local MolViewSpec state builder.
- [x] T068 Implement separate distance and dihedral local view modes.
- [ ] T069 Implement sortable/searchable/filterable violation table where notebook UI support is available, with static fallback.

## Phase G: Interpretation And Export

- [x] T080 Add explanatory markdown/copy for density, normalized density, violations, no violations, mapping warnings, unsupported ambiguity warnings, and empty restraint states.
- [ ] T081 Implement interpretation summary from computed diagnostics and outputs.
- [x] T082 Implement CSV export for violation table with provenance columns.
- [x] T083 Implement optional export for residue density table, mapping report, and evaluated restraints.
- [ ] T084 Implement optional MVS JSON state export.
- [ ] T085 Add notebook HTML export instructions or helper cell if repository conventions allow it.

## Phase H: Validation

- [ ] T100 Run restart-and-run-all on `happy-path-9l1v` remote mode.
- [ ] T101 Run restart-and-run-all on local-copy fixture mode.
- [ ] T102 Validate distance, dihedral, ambiguity, mapping, and density outputs against expected snapshots.
- [ ] T103 Validate missing-restraints, partial-mapping warning, and low-mapping abort behavior.
- [ ] T104 Validate visualization state existence and fallback behavior.
- [ ] T105 Validate no hidden state by running notebook from a clean kernel twice and comparing key outputs.
- [x] T106 Record dependency versions and runtime budget in validation snapshot.
- [ ] T107 Review all notebook and exported text for RUO/non-validation language compliance.

## Phase I: Documentation And Review

- [ ] T120 Update `docs-plan.md` checklist after notebook implementation.
- [ ] T121 Add tutorial run instructions inside the notebook.
- [ ] T122 Add how-to notes for switching PDB IDs, local files, thresholds, and selected model.
- [ ] T123 Add reference section for inputs, outputs, parameters, endpoints, and table schemas.
- [ ] T124 Add explanation section for NMR restraints, density, violations, ambiguity, and scientific limits.
- [ ] T125 Perform NMR domain review of formulas, terminology, caveats, and fixture behavior.
- [ ] T126 Perform final notebook review against `validation.md`.
