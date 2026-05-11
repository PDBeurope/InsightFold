# NMR Restraints Visualization Notebook Design

Source requirements: `specs/nmr_restraints/requirements.md`

## Notebook Target

Path: `notebooks/nmr_restraints_visualization.ipynb`

Runtime target: Python 3.11+ in local Jupyter or JupyterLab. Colab support is desirable but secondary. MolViewSpec/Mol* views must render inline inside notebook cells using the embedded HTML/IFrame pattern already used in `notebooks/homodimer_diagnostic.ipynb`; the notebook must not depend on opening external Mol* links for the primary experience.

Lifecycle status: draft spec. Implementation may begin with the provisional choices below, but fixture expected-value snapshots and NMR domain review are required before scientific sign-off.

## Configuration Defaults

The notebook should expose one `CONFIG` dictionary and display a readable configuration table near the top of the notebook.

```python
CONFIG = {
    "input_mode": "remote",
    "pdb_id": "9L1V",
    "model_index": 0,
    "violation_threshold_distance": 0.5,
    "violation_threshold_dihedral": 5.0,
    "local_context_radius": 6.0,
    "max_visible_restraints": 250,
    "cache_enabled": True,
    "remote_model_endpoint": "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif",
    "remote_restraint_endpoint": "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_nmr-data.str",
    "fallback_model_endpoint": "https://files.rcsb.org/download/{pdb_id}.cif",
    "mapping_warn_threshold": 0.70,
    "mapping_normal_threshold": 0.95,
    "density_mode": "normalized",
}
```

## Dependency Policy

| Dependency | Required? | Purpose | Install Source | Constraint |
|---|---|---|---|---|
| Python 3.11+ | yes | Runtime baseline | environment | Notebook target from PRD |
| `gemmi` | yes | Canonical mmCIF coordinate parsing and atom iteration | pip/conda | Parser source for geometry coordinates |
| `numpy` | yes | Distances, torsions, vector math | environment / pip | Lightweight scientific core |
| `pandas` | yes | Output tables, diagnostics, exports | environment / pip | DataFrame outputs |
| `pynmrstar` | yes, provisional | NMR-STAR parsing | pip | Replace or wrap only if fixtures fail |
| `requests` | yes | Remote retrieval | environment / pip | All calls wrapped with retry and timeout |
| `molviewspec` | yes | Declarative visualization state | pip | Target version 1.8.1 |
| `ipywidgets` | optional but preferred | Notebook controls and synchronized selection | pip | Static fallback required |
| `plotly` | optional | Interactive plots or compact sequence/density views | pip | Static/table fallback required |
| `matplotlib` | optional | Static density and violation plots for validation | environment / pip | Useful fallback |
| `scipy` | optional | Spatial neighborhood queries | pip | Use NumPy fallback if absent |
| Disallowed v1 | no | Refinement, MD, ML, clinical scoring | n/a | No OpenMM/Rosetta/GROMACS/PyTorch requirement |

## Notebook Sections

| Section | Purpose | Key Outputs | Requirement IDs |
|---|---|---|---|
| 1. Purpose and RUO framing | Explain exploratory scope, users, and non-validation limits | Intro markdown, interpretation vocabulary | REQ-023, REQ-024, NFR-008 |
| 2. Setup and configuration | Imports, dependency checks, `CONFIG`, display options | `CONFIG`, `dependency_report`, `notebook_log` | REQ-001, REQ-002, NFR-002 |
| 3. Input mode and parameters | Capture remote PDB ID or local file paths plus thresholds | `input_request`, `runtime_parameters` | REQ-001, REQ-002 |
| 4. Retrieval and provenance | Fetch or load model/restraint files and record provenance | `model_source`, `restraint_source`, `provenance` | REQ-003, NFR-004 |
| 5. Structure parsing | Parse mmCIF with Gemmi and normalize atom/residue records | `atom_table`, `residue_table`, `atom_lookup`, `model_summary` | REQ-004, REQ-005 |
| 6. Restraint parsing | Parse NMR-STAR compatible distance/dihedral restraints and unsupported records | `distance_restraints_raw`, `dihedral_restraints_raw`, `unsupported_records` | REQ-006, REQ-007 |
| 7. Atom mapping and coverage | Resolve restraint atoms against `atom_lookup` and enforce thresholds | `mapped_distance_restraints`, `mapped_dihedral_restraints`, `mapping_report` | REQ-008, REQ-009 |
| 8. Geometry and violation calculation | Evaluate distance, dihedral, ambiguous restraints, and thresholded violations | `evaluated_restraints`, `violation_table`, `geometry_diagnostics` | REQ-010, REQ-011, REQ-012, REQ-015 |
| 9. Residue density calculation | Compute absolute and normalized density | `residue_density_table`, `density_summary` | REQ-013, REQ-014 |
| 10. Global visualization | Build global density-colored MolViewSpec state and violation overlay | `global_mvs_state`, rendered global view | REQ-016, REQ-017, REQ-022 |
| 11. Sequence and residue navigation | Provide residue selector, sequence/density view, and authoritative selection state | `selection_state`, filtered tables | REQ-018, REQ-020, REQ-021 |
| 12. Local evidence visualization | Build local neighborhood view and local restraint overlays | `local_context_table`, `local_mvs_state`, rendered local view | REQ-018, REQ-019, REQ-022 |
| 13. Interpretation summary | Display density, violation, no-violation, mapping, and unsupported-record explanations | `interpretation_summary` | REQ-023, REQ-024 |
| 14. Export artifacts | Export CSV tables, optional MVS state JSON, notebook HTML instructions | `export_manifest` | REQ-026 |
| 15. Validation snapshot | Show fixture checks, runtime, dependency versions, and known caveats | `validation_snapshot` | NFR-001, NFR-003, NFR-005 |

Each major code section should be preceded by markdown that explains the biological interpretation and limitations in plain structural-biology language.

## Data Flow

```text
input mode + PDB ID or local paths + runtime parameters
  -> retrieve/load mmCIF and NMR-STAR files
  -> record provenance and cache status
  -> parse structure with Gemmi
  -> parse NMR-STAR distance and dihedral restraints
  -> normalize atom identifiers and build atom lookup
  -> map restraints to atoms and enforce mapping coverage
  -> evaluate ambiguous distance restraints
  -> compute distance and dihedral geometry
  -> compute violation magnitudes and priority ordering
  -> compute residue absolute and normalized density
  -> build global MolViewSpec state
  -> select residue or violation
  -> compute local neighborhood and local evidence tables
  -> build local MolViewSpec state
  -> explain, export, and validate
```

## Variable And Artifact Handoff

| Name | Produced By | Consumed By | Type / Shape | Notes |
|---|---|---|---|---|
| `CONFIG` | Setup | All sections | dict | Runtime defaults, thresholds, endpoints, cache policy |
| `dependency_report` | Setup | Validation | DataFrame/dict | Package presence and versions |
| `input_request` | Input mode | Retrieval | dict | `input_mode`, `pdb_id`, `model_path`, `restraint_path` |
| `runtime_parameters` | Input mode | Geometry, visualization | dict | `model_index`, thresholds, radius, density mode |
| `model_source` | Retrieval/local load | Structure parsing, visualization | dict | Raw path/URL, bytes/text path, format, cache status |
| `restraint_source` | Retrieval/local load | Restraint parsing | dict | Raw path/URL, text path, cache status |
| `provenance` | Retrieval/local load | All outputs, export | DataFrame/dict | Source URLs, timestamps, parser versions |
| `atom_table` | Structure parsing | Mapping, geometry, visualization | DataFrame | One row per atom in selected model |
| `residue_table` | Structure parsing | Density, navigation, visualization | DataFrame | One row per residue with chain and coordinate summary |
| `atom_lookup` | Structure parsing | Mapping | dict | Keyed by `(auth_asym_id, auth_seq_id, auth_comp_id, auth_atom_id)` |
| `distance_restraints_raw` | Restraint parsing | Mapping, ambiguity grouping | DataFrame | Parsed distance rows with bounds and source loop |
| `dihedral_restraints_raw` | Restraint parsing | Mapping | DataFrame | Parsed dihedral rows with bounds and atom identifiers |
| `unsupported_records` | Restraint parsing | Diagnostics, interpretation | DataFrame | Unsupported categories, missing fields, parse failures |
| `mapped_distance_restraints` | Mapping | Geometry, density, visualization | DataFrame | Distance restraints with mapped atom coordinates or unmapped reason |
| `mapped_dihedral_restraints` | Mapping | Geometry, density, visualization | DataFrame | Dihedral restraints with mapped atom coordinates or unmapped reason |
| `mapping_report` | Mapping | Validation, interpretation | DataFrame/dict | Coverage, counts, examples, abort/warn status |
| `evaluated_restraints` | Geometry | Density, export | DataFrame | All evaluated mapped restraints including satisfied restraints |
| `violation_table` | Geometry | Tables, visualization, export | DataFrame | Thresholded and unthresholded violation fields |
| `residue_density_table` | Density | Global visualization, export | DataFrame | Absolute and normalized density by residue |
| `selection_state` | Navigation | Local context and views | dict | Selected residue or selected violation ID |
| `local_context_table` | Local evidence | Local visualization, export | DataFrame | Neighboring residues and local restraints |
| `global_mvs_state` | Global visualization | Rendering/export | MVS builder/state or JSON | Deterministic state generated from notebook data |
| `local_mvs_state` | Local visualization | Rendering/export | MVS builder/state or JSON | Deterministic state regenerated from `selection_state` |
| `interpretation_summary` | Interpretation | Export/review | Markdown/DataFrame | Caveats and evidence summaries |
| `export_manifest` | Export | Validation/review | DataFrame/dict | Exported artifact paths and provenance |
| `notebook_log` | All sections | Validation/export | list/DataFrame | Warnings, failures, elapsed times, excluded counts |

## Function Candidates

| Function | Responsibility |
|---|---|
| `validate_pdb_id(value)` | Validate and normalize four-character PDB IDs. |
| `fetch_with_retries(url, cache_key, config)` | Retrieve remote files with timeout, retries, and optional cache. |
| `load_local_file(path, expected_suffixes)` | Load local mmCIF or NMR-STAR file as data with validation. |
| `parse_structure_with_gemmi(model_source, model_index)` | Produce atom and residue tables plus lookup from selected model. |
| `parse_nmrstar_restraints(restraint_source)` | Extract distance, dihedral, and unsupported restraint records. |
| `normalize_restraint_atoms(restraint_table)` | Normalize atom identity columns to canonical deposited identifiers. |
| `map_restraints(restraints, atom_lookup)` | Map restraint atoms to atom coordinates and record unmapped reasons. |
| `compute_mapping_report(mapped_tables)` | Produce coverage counts and threshold decision. |
| `group_ambiguous_distance_restraints(distance_table)` | Identify shared-ID `OR` groups and candidate members. |
| `evaluate_distance_restraints(mapped_distance_table)` | Compute measured distances and distance violation magnitudes. |
| `compute_dihedral_angle(coords4)` | Compute torsion angle normalized to `[-180, +180]`. |
| `evaluate_dihedral_restraints(mapped_dihedral_table)` | Compute torsions and circular violation magnitudes. |
| `compute_residue_density(evaluated_restraints, residue_table)` | Produce absolute and normalized residue density table. |
| `build_global_mvs_state(model_source, residue_density_table, violation_table, config)` | Create global density and violation visualization state. |
| `compute_local_context(selection_state, atom_table, residue_table, evaluated_restraints, config)` | Produce local neighborhood and local restraint subset. |
| `build_local_mvs_state(local_context_table, evaluated_restraints, config)` | Create local evidence visualization state. |
| `render_mvs_state(state, label)` | Render MolViewSpec/Mol* HTML inline with fallback and error isolation. Use `state.molstar_html()`, base64-encode the HTML, and display it through an `IPython.display.IFrame` data URI as in `notebooks/homodimer_diagnostic.ipynb`. |
| `export_outputs(tables, states, provenance, config)` | Write CSV/JSON artifacts and return export manifest. |

## Parser And Algorithm Choices

| Decision | Rationale | Alternatives Considered | Risk |
|---|---|---|---|
| Use PDBe entry-file endpoints as canonical remote source | PRD includes PDBe model and NMR restraint endpoints and they provide paired entry files for NMR workflows | RCSB model endpoints, manual download only | Endpoint availability can change; keep fallback and local-file mode. |
| Use Gemmi for mmCIF parsing | PRD explicitly identifies Gemmi as canonical coordinate parser | BioPython, ad hoc text parsing | Gemmi object model must be wrapped carefully to preserve deposited auth identifiers. |
| Use `pynmrstar` provisionally for NMR-STAR parsing | Purpose-built parser for STAR files and compatible with notebook use | Custom parser, BMRB tooling | Fixtures may expose category/loop variants; parser choice remains reviewable. |
| Use exact deposited `auth_*` mapping | Matches PRD and avoids introducing external atom-remapping ambiguity | Label identifiers, fuzzy residue matching, external atom mapping | Some entries may map poorly; coverage diagnostics must be clear. Mapping coverage is measured over logical restraints; an evaluable ambiguous `OR` group counts as `1/1` mapped even if some candidate members are unmapped. |
| Analyze selected model only | Keeps v1 from becoming ensemble-wide analysis | Ensemble aggregation | NMR users may expect ensemble behavior; UI must label selected model clearly. |
| Evaluate `OR` ambiguous distances by smallest measured candidate | Matches PRD and common v1 simplification | Weighted ambiguity, average over members, pseudoatom expansion | Oversimplifies some NMR semantics; warn that unsupported semantics are excluded. |
| Treat Python state as authoritative | Prevents Mol* frontend state drift and supports reproducibility | Rely on frontend selections | Widget/Mol* state can desynchronize; regenerate views from notebook state. |
| Use global view for density and thresholded violations only | Reduces visual clutter | Render all restraints globally | Dense restraint networks become unreadable and slow. |
| Embed MolViewSpec views inline in notebook cells | Matches the direction set by `notebooks/homodimer_diagnostic.ipynb` and keeps the notebook self-contained | External Mol* links only, local filesystem paths in browser | Inline rendering still needs a browser-accessible model payload; implementation must use URL-backed remote files or embedded/served local model content rather than arbitrary local paths. |

## Visualization Plan

| View | Contents | Interaction | Fallback |
|---|---|---|---|
| Global density view | Complete structure as cartoon colored by normalized or absolute residue density | Rotate/zoom; density mode toggle; thresholded violation overlay | Static density table and optional 2D residue plot |
| Global violation overlay | Highest-magnitude violated restraints up to `max_visible_restraints` | Display distance/dihedral type, magnitude, involved residues | Violation table only |
| Sequence/density view | Linear residue navigation with density and local violation indicators | Click or select residue to update `selection_state` | Dropdown/numeric residue selector |
| Local evidence view | Selected residue, spatial neighbors, all restraints involving selected residue | Radius and threshold controls update local context | Local context and restraint tables |
| Distance restraint overlay | Pairwise local distance restraints colored by violation magnitude | Separate distance mode | Distance rows filtered for selected residue |
| Dihedral restraint view | Local torsion-associated residues/atoms and dihedral violation indicators | Separate dihedral mode | Dihedral rows filtered for selected residue |

## Error Handling And Graceful Degradation

- Retrieval functions return structured success/failure objects; downstream parsing does not run on failed retrieval.
- Parser failures append to `notebook_log` and set explicit empty outputs rather than leaving undefined variables.
- Mapping abort below 70% sets `analysis_status = "aborted_low_mapping"` and skips geometry/visualization sections that depend on mapped restraints.
- Visualization functions use isolated `try/except` blocks and show fallback tables when rendering fails.
- Rendering should follow the homodimer diagnostic helper pattern: call `state.molstar_html()`, base64-encode the result, and display it as an inline `IFrame` so the view appears in the executed notebook cell. Local-file mode must embed or serve the mmCIF content in a notebook-accessible way before building the MVS state; it must not pass an arbitrary local filesystem path to browser-side Mol*.
- Empty results are valid artifacts when explicitly labeled, especially no compatible restraints and no thresholded violations.

## Mandatory And Optional Sections

Mandatory for v1:

- Setup and configuration.
- Remote or local input.
- Retrieval/loading and provenance.
- Structure parsing.
- Restraint parsing.
- Mapping and coverage diagnostics.
- Geometry and density calculations.
- Residue density and violation tables.
- Global density view or fallback.
- Scientific interpretation caveats.
- Validation snapshot.

Optional/stretch for v1:

- Full sequence click interaction if widget support is unstable.
- Plotly sequence/density view.
- CSV export for all satisfied restraints.
- MVS state JSON export.
- Colab file-upload widgets.
