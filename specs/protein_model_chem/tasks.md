# Protein Model Chemistry Notebook Tasks

## Phase 1: Spec And Fixture Readiness

- [x] T001 Convert PRD into directory-style notebook spec pack.
- [x] T002 Run `$notebook-spec-review` and resolve non-fixture blocking spec findings.
- [x] T003 Run `$fixture-selection` to pin FX-001 happy-path, FX-002 AFDB confidence caveat, FX-005 invalid identifier, and FX-006 mutation mismatch fixtures.
- [ ] T004 Obtain domain review for interaction definitions, thresholds, mutation model, graph weighting, and severity language.
- [x] T005 Adopt provisional v1 structure parser and supported file formats: BioPython wrapper, mmCIF preferred, PDB supported where parser allows.
- [x] T006 Adopt provisional PDBe policy: prefer PDBe when resolution and residue coverage are acceptable, build biological assembly if available, otherwise show caveat.
- [ ] T007 Decide whether local user-uploaded structures are in v1 or deferred.
- [x] T008 Adopt provisional thresholds, weights, mutation radius, sidechain-only mutation model, MolViewSpec/Mol*/Plotly/NetworkX visualization, and export formats.

## Phase 2: Data Contracts And Retrieval

- [x] T010 Implement structure request validation for AFDB, PDBe, and optional local inputs.
- [x] T011 [P] Implement AFDB metadata and coordinate retrieval with provenance capture.
- [x] T012 [P] Implement PDBe coordinate and metadata retrieval with provenance capture.
- [ ] T013 Implement source preference logic, including PDBe-first behavior and manual override.
- [x] T014 Implement retrieval failure handling and fixture/cache fallback messages.
- [x] T015 Implement raw payload validation for coordinate text, metadata, and confidence fields.

## Phase 3: Structure Parsing And Normalization

- [x] T020 Implement selected mmCIF/PDB parser wrapper.
- [x] T021 Build atom table with chain IDs, residue IDs, atom names, coordinates, alternate conformations, and source provenance.
- [x] T022 Build residue table with normalized residue identifiers, residue names, representative coordinates, and quality flags.
- [x] T023 Map AFDB pLDDT or confidence metadata to residues.
- [ ] T024 Map PDBe experimental metadata, missing residues, alternate conformations, and quality caveats to residues.
- [x] T025 Validate mutation input against residue table and stop clearly on mismatches.

## Phase 4: Wild-Type Interaction Detection

- [x] T030 Implement documented interaction configuration table.
- [x] T031 Implement hydrogen-bond detection with geometric evidence and confidence.
- [x] T032 Implement salt-bridge detection with geometric evidence and confidence.
- [x] T033 Implement hydrophobic-contact detection with geometric evidence and confidence.
- [x] T034 Implement aromatic-interaction detection with geometric evidence and confidence.
- [x] T035 Implement disulfide detection with geometric evidence and confidence.
- [x] T036 Implement steric-clash detection with geometric evidence and confidence.
- [x] T037 Implement local-packing contact detection and optional burial proxy.
- [x] T038 Combine interaction records into a validated WT interaction table.

## Phase 5: Graph And Confidence Analysis

- [x] T040 Build residue-level interaction graph or equivalent edge/node tables.
- [x] T041 Compute local degree, interaction density, centrality-like values, fragmentation, and constraint redundancy.
- [x] T042 Implement structure confidence, interaction confidence, and mechanistic confidence classification.
- [x] T043 Add low-confidence and unsupported-context warnings.

## Phase 6: Mutation Perturbation

- [x] T050 Implement provisional sidechain-substitution-only mutation model with no global minimization and no conformational sampling.
- [x] T051 Estimate lost or weakened compatible interactions near mutation site.
- [x] T052 Estimate new steric incompatibility or clash signals.
- [x] T053 Estimate packing, cavity, or local constraint-density changes where supported.
- [x] T054 Compare WT and mutant local interaction environments.
- [x] T055 Compute local network delta around mutation site.
- [x] T056 Produce compatibility delta table with evidence and confidence.

## Phase 7: Visualizations And Interpretation

- [x] T060 Create structure overview and provenance/confidence display.
- [x] T061 Create mutation-site local neighborhood visualization.
- [x] T062 Create interaction-type summary plot.
- [x] T063 Create residue interaction graph or ego-network visualization.
- [ ] T064 Create WT versus mutant perturbation visualization.
- [x] T065 Create severity/confidence summary panel.
- [x] T066 Implement optional MolViewSpec/Mol* views with API verification, generated/skipped summary, `.mvsj` export where requested, and Mol* viewer URL fallback.
- [x] T067 Generate templated observation, interpretation, and hypothesis summary.
- [x] T068 Implement structural disruption likelihood assignment with visible threshold/caveat policy.
- [x] T069 Implement provisional severity rubric and label it uncalibrated until domain review.

## Phase 8: Export And Notebook UX

- [x] T070 Add export bundle for tables, summaries, figures, provenance, and caveats.
- [x] T071 Add markdown explanations before every major code section.
- [x] T072 Add interpretation vocabulary and forbidden-claim guardrails.
- [x] T073 Add final "what to trust / what not to trust" section.
- [x] T074 Add runtime and version snapshot cell.

## Phase 9: Validation

- [x] T080 Run static notebook validation for hidden-state and execution-order hazards.
- [x] T081 Run restart-and-run-all validation on FX-001.
- [ ] T082 Run confidence downgrade validation on FX-002.
- [ ] T083 Run invalid identifier validation on FX-005.
- [ ] T084 Run mutation mismatch validation on FX-006.
- [ ] T085 Validate interaction table and graph contracts.
- [ ] T086 Validate explanation outputs contain observation/interpretation/hypothesis labels.
- [ ] T087 Validate outputs avoid clinical, pathogenicity, exact-energy, and deterministic folding claims.
- [x] T088 Record validation report under `specs/protein_model_chem/validation-report.md`.
- [ ] T089 Validate config thresholds, weights, mutation radius, export formats, and subruntime targets against the spec.
- [ ] T089a Freeze exact local interaction counts and perturbation signal counts for FX-001 after first trusted implementation run.
- [ ] T089b Confirm FX-003/FX-004 chain and residue numbering under the selected BioPython parser.

## Phase 10: Review And Graduation

- [x] T090 Run `$notebook-execution-validation` after implementation.
- [ ] T091 Run `$notebook-review` after execution validation.
- [ ] T092 Resolve domain review findings on thresholds, mutation interpretation, and confidence language.
- [ ] T093 Decide readiness state: continue prototype, ready for beta, or archive/re-scope.
