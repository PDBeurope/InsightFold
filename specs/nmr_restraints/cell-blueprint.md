# Jupyter Cell Blueprint

Source design: `specs/nmr_restraints/notebook-design.md`

Kernel assumption: Python 3.11+ local Jupyter/JupyterLab. The notebook must run top-to-bottom without relying on hidden widget or Mol* frontend state.

| Cell ID | Type | Purpose | Inputs | Outputs | User Editable? | Hidden-State Risk | Validation Hook | Requirements |
|---|---|---|---|---|---|---|---|---|
| C001 | markdown | Introduce scope, RUO framing, and interpretation guardrails | none | intro text | no | Low; static text only | Review required caveat language | REQ-023, REQ-024, NFR-008 |
| C002 | code | Import packages, detect versions, and report dependency availability | environment | `dependency_report` | no | Medium if later cells re-import conditionally; keep imports centralized here | Dependency audit and runtime baseline | NFR-001, NFR-002, NFR-003 |
| C003 | code | Define `CONFIG`, mode selector, and input validation helpers | user-edited config | `CONFIG`, `input_request`, `runtime_parameters` | yes | Medium if later cells mutate `CONFIG`; downstream cells must treat it as read-only | Input validation checks for remote/local modes | REQ-001, REQ-002 |
| C004 | markdown | Explain how remote mode, local mode, and fixture defaults relate | `CONFIG` | user guidance | no | Low | Manual UX review | REQ-001, REQ-002, REQ-023 |
| C005 | code | Retrieve or load model/restraint files and record provenance | `CONFIG`, `input_request` | `model_source`, `restraint_source`, `provenance`, `notebook_log` | no | Medium if cache policy or paths are mutated ad hoc; keep cache helper pure; retain browser visualization URL separately from parsed cache path | Retrieval contract and local-file parity | REQ-003, NFR-004, NFR-005 |
| C006 | code | Parse structure and normalize residue/atom tables | `model_source`, `runtime_parameters` | `atom_table`, `residue_table`, `atom_lookup`, `model_summary` | no | Medium if altloc filtering depends on previous cell state; compute from sources each run | Structure count snapshot and altloc diagnostics | REQ-004, REQ-005 |
| C007 | code | Parse NMR-STAR distance and dihedral restraints plus unsupported rows | `restraint_source` | `distance_restraints_raw`, `dihedral_restraints_raw`, `unsupported_records` | no | Medium if parser holds implicit saveframe state; keep output tables explicit | Parser count snapshot and unsupported-record review | REQ-006, REQ-007 |
| C008 | code | Map restraint atoms and compute logical-restraint coverage report | parsed tables, `atom_lookup` | `mapped_distance_restraints`, `mapped_dihedral_restraints`, `mapping_report` | no | High if logical coverage is inferred from previously filtered tables; recompute deterministically here | Mapping threshold and ambiguous-group diagnostics | REQ-008, REQ-009 |
| C009 | code | Evaluate ambiguity, geometry, violations, and density | mapped tables, `residue_table`, thresholds | `evaluated_restraints`, `violation_table`, `residue_density_table`, `geometry_diagnostics` | no | High if selected ambiguous member or density counts are cached separately; derive from mapped tables each run | Distance, dihedral, ambiguity, density snapshots | REQ-010, REQ-011, REQ-012, REQ-013, REQ-015 |
| C010 | markdown | Explain mapping, density, and violation outputs before visualization | `mapping_report`, computed tables | interpretation text | no | Low | Scientific language review | REQ-023, REQ-024 |
| C011 | code | Build and render global MolViewSpec/Mol* state inline | `model_source`, `residue_density_table`, `violation_table`, `CONFIG` | `global_mvs_state`, global render output | no | High if render helper stores frontend-only selection, browser-local paths, or nested full-file data URIs; use `state.molstar_html()` plus base64 `IFrame` each run and a public PDBe model URL for remote/cache mode | Visualization state existence and fallback | REQ-016, REQ-017, REQ-022 |
| C012 | code | Create authoritative `selection_state` from residue selector or violation selection | `residue_table`, `violation_table`, widgets or fallback controls | `selection_state` | yes | High if widget callbacks mutate hidden globals; keep one explicit state object shown in output | Hidden-state and deterministic selection checks | REQ-018, REQ-020, REQ-021 |
| C013 | code | Compute local context and local evidence tables from `selection_state` | `selection_state`, structure tables, `evaluated_restraints`, `CONFIG` | `local_context_table`, filtered tables | no | Medium if local context incrementally mutates prior outputs; recompute from authoritative state | Local residue snapshot and threshold behavior | REQ-018, REQ-019, REQ-020 |
| C014 | code | Build and render local MolViewSpec/Mol* state inline | local context inputs, `model_source` | `local_mvs_state`, local render output | no | High for same reasons as global render; no filesystem-only browser assumptions and no nested full cached mmCIF data URI for remote/cache mode | Local view presence and render fallback | REQ-018, REQ-019, REQ-022 |
| C015 | code | Export CSV and optional JSON/state artifacts with provenance | computed tables, states, `provenance`, `CONFIG` | `export_manifest` | yes | Medium if export cell depends on unstated previous selections; document whether export is global or selection-specific | Export artifact check | REQ-026, NFR-004 |
| C016 | markdown / validation | Record validation snapshot, known caveats, and review status | all major outputs | `validation_snapshot` summary | no | Low if snapshot is derived directly from visible tables; avoid manual edits in notebook body | Restart-and-run-all evidence and docs completeness | NFR-001, NFR-003, NFR-005 |

## Cell Order Rules

- Cells C002 through C016 must execute in order with a fresh kernel.
- No later cell may silently redefine parser, geometry, or rendering helpers created earlier.
- Widget state is optional; when widgets are unavailable, the notebook must expose equivalent visible fallback controls without changing downstream variable names.
- Cache behavior must be controlled only by `CONFIG`, not by ad hoc local variables introduced later in the notebook.
