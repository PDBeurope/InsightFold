# NMR Restraints Visualization Notebook Requirements

Source PRD: `specs/nmr_restraints/nmr_restraints_prd.md`

## Goal

Build a self-contained InsightFold notebook that retrieves or loads an NMR structure model and its deposited restraint file, maps distance and dihedral restraints onto deposited coordinates, computes residue-level restraint density and restraint violations, and provides synchronized tabular, sequence, and 3D MolViewSpec/Mol* views for exploratory evidence inspection.

The notebook must remain an exploratory research-use evidence viewer. It must not become a formal NMR validation engine, automated quality scorer, clinical tool, or structure-refinement workflow.

## Users

- Structural biologists inspecting where deposited restraints support or conflict with an NMR model.
- Protein biochemists and functional researchers who need residue-level evidence context without specialist NMR restraint tooling.
- Computational biology researchers who need reproducible parsing, mapping, geometry, and exportable tables.
- PDB users exploring public NMR entries.
- InsightFold maintainers and reviewers who need explicit data contracts, fixtures, and scientific caveats.

## User Stories

| ID | User Story |
|---|---|
| US-001 | As a structural biologist, I want to enter a PDB ID and see restraint density mapped onto the structure, so that I can quickly identify experimentally dense and sparse regions. |
| US-002 | As a protein biochemist, I want to select a residue and see its local restraint network, so that I can understand whether a functional region is directly supported by deposited evidence. |
| US-003 | As a PDB user, I want violated restraints listed and highlighted spatially, so that I can inspect local inconsistencies without reading raw restraint files. |
| US-004 | As a notebook developer, I want unsupported records and mapping failures reported explicitly, so that partial analyses do not appear complete. |
| US-005 | As an educator, I want no-violation outputs to include caveats, so that students do not mistake absence of highlighted violations for structural correctness. |

## Functional Requirements

| ID | Requirement | Acceptance Criteria |
|---|---|---|
| REQ-001 | Accept remote PDB ID input | Notebook exposes a four-character PDB ID input, validates it with `^[A-Za-z0-9]{4}$`, normalizes it to uppercase for labels and lowercase for PDBe URLs where needed, and stops before retrieval when invalid. |
| REQ-002 | Accept local-file input | Notebook supports local mmCIF model files and local NMR-STAR `.str` restraint files; missing or unsupported local paths fail before parsing with a readable message. |
| REQ-003 | Retrieve remote files | Remote mode fetches a model file and NMR restraint file from documented endpoints, records URLs and retrieval timestamps, retries transient failures up to three times, and reports unavailable files clearly. |
| REQ-004 | Parse structure coordinates | Gemmi-based parsing produces `atom_table`, `residue_table`, `model_summary`, and an atom lookup keyed by deposited `auth_asym_id`, `auth_seq_id`, `auth_comp_id`, and `auth_atom_id`. |
| REQ-005 | Filter model and alternate conformations | Notebook analyzes a selected `model_index` and ignores unsupported alternate conformations except primary/default, blank, or `.` altloc records; ignored altloc counts are reported. |
| REQ-006 | Parse compatible restraints | Notebook extracts compatible distance restraints and dihedral restraints from NMR-STAR loops/saveframes into unified `distance_restraints` and `dihedral_restraints` tables with source-loop provenance. |
| REQ-007 | Identify unsupported restraint records | Unsupported categories, incomplete rows, unsupported ambiguity semantics, and unparseable bounds are counted and shown in diagnostics without being silently dropped. |
| REQ-008 | Map restraint atoms | Restraint atoms map to structure atoms by exact `auth_*` identity; mapping coverage is calculated over logical restraints, not raw ambiguous member rows. A non-ambiguous logical restraint maps as `1/1` when all required atoms map. An `OR` ambiguous logical distance restraint maps as `1/1` when at least one candidate member maps and can be evaluated; member-level mapping counts are reported separately as diagnostics. |
| REQ-009 | Enforce mapping coverage thresholds | Mapping coverage above 95% proceeds normally, 70% to 95% proceeds with warning, and below 70% aborts analysis before geometry calculations. |
| REQ-010 | Evaluate ambiguous distance restraints | For restraints with shared `_Gen_dist_constraint.ID` and `Member_logic_code = OR`, all candidate atom-pair distances are computed, the smallest measured distance is selected, and only that member is used for violation, reporting, and visualization in v1. |
| REQ-011 | Compute distance violations | Euclidean distances are measured in Angstroms. All v1 distance restraints are expected to provide both lower and upper bounds; rows missing either bound are unsupported and excluded with diagnostics. Violation magnitude is zero when inside bounds, `lower - measured` when below lower bound, and `measured - upper` when above upper bound. |
| REQ-012 | Compute dihedral violations | Torsion angles are computed from four mapped atoms, normalized to `[-180, +180]`, and evaluated against allowed intervals with circular geometry including wrapped intervals. |
| REQ-013 | Compute residue density | Absolute density counts unique restraints associated with each residue; normalized density is `absolute_density / structure_mean_density`, where `structure_mean_density = N_restraints / N_residues`. |
| REQ-014 | Produce documented output tables | Notebook produces residue density, violation, mapping diagnostics, unsupported-record, and provenance tables with documented columns and stable names. |
| REQ-015 | Prioritize violations | Violation table defaults to descending `violation_magnitude`, then `restraint_type`, then residue number; thresholds for displayed distance and dihedral violations are configurable. |
| REQ-016 | Render global density view | Initial MolViewSpec/Mol* view displays the complete structure as cartoon representation, colored by normalized density by default, with a toggle for absolute density. |
| REQ-017 | Control global restraint clutter | Global view renders only violated restraints by default and caps visible restraints with `max_visible_restraints`, showing highest-magnitude violations first. |
| REQ-018 | Render local evidence view | Selecting a residue computes a spatial neighborhood within `local_context_radius` Angstroms and displays selected residue, neighboring residues, local violations, and all restraints involving the selected residue. |
| REQ-019 | Separate distance and dihedral views | Distance and dihedral overlays are shown in separate panels or modes by default to reduce clutter. |
| REQ-020 | Synchronize navigation | Residue selector, sequence view, violation table selection, and local 3D view share one authoritative Python notebook state and regenerate deterministic visualization states. |
| REQ-021 | Provide notebook-friendly interactivity | Where dependencies permit, tables are sortable/searchable/filterable and selections update local context; fallback static tables remain available when widgets fail. |
| REQ-022 | Isolate visualization failures | Failure in one visualization panel produces a visible warning and does not prevent other tables or visualizations from rendering. |
| REQ-023 | Provide scientifically cautious explanations | Notebook text states that density is experimental coverage, not confidence, and violations are local inconsistencies, not automatic structural errors. |
| REQ-024 | Handle no-violation and missing-restraint states | No-violation outputs state that absence of thresholded violations is not a validation verdict; missing compatible restraints exit gracefully with diagnostics. |
| REQ-025 | Cache optional artifacts | When enabled, downloaded files, parsed structures, parsed restraints, computed tables, and visualization state JSON are cached with key `{pdb_id}_{model_index}_{config_hash}`. |
| REQ-026 | Export user artifacts | Notebook can export violation tables as CSV and the executed notebook as HTML; exported tables include provenance fields. |

## Non-Functional Requirements

| ID | Requirement | Acceptance Criteria |
|---|---|---|
| NFR-001 | Runtime target | Happy-path fixture runs top-to-bottom in under 30 seconds after dependencies are available on local Jupyter CPU for structures up to about 500 residues and 10,000 restraints. |
| NFR-002 | Dependency budget | Uses Python 3.11+, Gemmi, NumPy, pandas, pynmrstar or equivalent, MolViewSpec 1.8.1, and lightweight notebook visualization dependencies; no GPU-only, MD, refinement, or ML dependency is required. |
| NFR-003 | Reproducibility | Notebook is restart-and-run-all compatible and does not depend on manual execution order or hidden frontend state. |
| NFR-004 | Provenance | Every retrieved or loaded source records file path or URL, retrieval/load timestamp, parser version where available, and whether cached data was used. |
| NFR-005 | Graceful degradation | Network, parser, mapping, and visualization failures produce explicit diagnostics and leave downstream variables in documented empty or aborted states. |
| NFR-006 | Accessibility | Visual encodings do not rely exclusively on color; tables and text summaries expose the same core evidence as 3D views. |
| NFR-007 | Security | Remote and local files are parsed as data only; notebook never executes downloaded or uploaded file content. |
| NFR-008 | RUO scope | Notebook text and exports retain research-use, non-clinical, non-validation framing. |

## Non-Goals

- Automated quality scores, global confidence metrics, pass/fail labels, or validation verdicts.
- NMR refinement, re-refinement, structure repair, or energy minimization.
- Ensemble-wide statistics or conformer comparison in v1.
- Weighted ambiguity semantics, pseudoatom expansion, dynamic averaging, or full assignment-resolution semantics in v1.
- Machine-learning interpretation.
- Cross-structure comparison, batch workflows, or production platform features.
- Clinical, diagnostic, regulated, or patient-specific interpretation.

## Edge Cases And Failure Behavior

| Case | Required Behavior |
|---|---|
| Invalid PDB ID | Stop before retrieval; display expected four-character format. |
| Remote model missing | Stop retrieval-dependent workflow; show endpoint and exception/status. |
| Remote restraint file missing | Show missing-restraint state and do not attempt geometry analysis. |
| Local file path missing or unsupported extension | Stop before parsing; show accepted formats. |
| Malformed mmCIF | Stop before restraint mapping; show parser diagnostic. |
| Malformed NMR-STAR or no compatible restraints | Exit gracefully with empty-restraint explanation. |
| Mapping coverage below 70% | Abort geometry analysis and show mapping table with example unmapped restraints. |
| Mapping coverage 70% to 95% | Continue with warning and mark outputs as partial. |
| Unsupported ambiguity semantics | Exclude unsupported records, count them, and warn that v1 only supports explicit `OR` smallest-distance behavior. |
| Wrapped dihedral interval | Use circular interval logic and include a validation check. |
| No violations above thresholds | Show no-violation message with caveat that this is not a validation verdict. |
| Visualization backend unavailable | Still produce tables, diagnostics, exports, and text explanations. |
| Large structure or dense restraint set | Apply visibility limits and grouped residue coloring; warn if runtime or rendering may exceed target. |

## Open Questions

| ID | Question | Type | Owner | Resolution |
|---|---|---|---|---|
| Q-001 | Which exact fixture entries should cover happy path, partial mapping, ambiguous restraints, dihedral restraints, no violations, and missing/unsupported restraints? | blocking | Fixture selection / domain reviewer | Provisional happy path is `9L1V`; fixture-selection must record expected counts before implementation validation. |
| Q-002 | Should PDBe entry-file endpoints be canonical over RCSB endpoints for all remote retrieval? | assumption | Implementation lead | Proceed with PDBe canonical for v1 model/restraint retrieval and retain RCSB model endpoint as documented fallback. |
| Q-003 | Which NMR-STAR parser is canonical if `pynmrstar` fails on selected fixtures? | assumption | Implementation lead / domain reviewer | Start with `pynmrstar`; switch or wrap only if fixture parsing fails. |
| Q-004 | Should v1 analyze only one selected NMR model, or expose multi-model controls while avoiding ensemble analysis? | assumption | Product owner / domain reviewer | Use `model_index = 0` default with a visible selector; no ensemble-wide summaries in v1. |
| Q-005 | How should ligand, nucleic-acid, and non-protein restraints be handled? | non-blocking | Domain reviewer | Treat as unsupported or warn in v1 unless fixture-selection explicitly includes support. |
| Q-006 | What thresholds define low, medium, and high violation colors? | assumption | Domain reviewer | Use configurable provisional thresholds in notebook design; require review before release. |
| Q-007 | Should local-file mode use file path inputs, upload widgets, or both? | non-blocking | Implementation lead | Prefer file path input for local Jupyter and upload widgets for Colab/JupyterHub; visual outputs must be embedded into notebook cells using the inline MolViewSpec HTML/IFrame pattern from `notebooks/homodimer_diagnostic.ipynb`. |
| Q-008 | Should satisfied restraints be exportable even though primary violation table emphasizes violations? | non-blocking | Product owner | Export mapped restraint table as optional advanced artifact if implementation time permits. |
| Q-009 | What domain-review checklist is required before public educational release? | blocking for release | Domain reviewer | Must be defined before public release, but not before initial notebook implementation. |
