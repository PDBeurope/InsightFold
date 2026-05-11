# Protein Model Chemistry Notebook Design

## Notebook Target

Path: `notebooks/protein_model_chem.ipynb`

Runtime target: local Jupyter or Google Colab CPU. The happy-path fixture should run in under 2 minutes after dependencies are installed. The notebook should avoid GPU-only, MD, and heavyweight force-field dependencies.

Subtargets from `protein-bond-type_energy.md`:

- Structure analysis should complete in under 30 seconds for an average monomer fixture.
- Local mutation perturbation recomputation should complete in under 10 seconds after structure normalization.

Runtime platform: Python 3.11+ in Jupyter Notebook, JupyterLab, or Colab CPU.

Lifecycle status: draft spec with provisional v1 implementation choices. Implementation should not begin beyond scaffolding until fixtures are pinned, but parser, interaction-threshold, weighting, visualization, and mutation-model choices are now specified for v1 and require domain review before scientific sign-off.

## Dependency Policy

| Dependency | Required? | Purpose | Install Source | Constraint |
|---|---|---|---|---|
| Python 3.11+ | yes | Runtime baseline | environment | Required target from supporting spec |
| `biopython` | yes | mmCIF/PDB parsing and residue/atom access | pip | Parser wrapper must normalize numbering and alternate conformations |
| `numpy` | yes | Geometry, distances, arrays | environment / pip | Lightweight scientific core |
| `pandas` | yes | Residue, interaction, mutation, and export tables | environment / pip | DataFrame outputs |
| `scipy` | yes | Distance calculations and spatial queries | environment / pip | Avoid heavyweight simulation use |
| `plotly` | yes | Interactive network and perturbation plots | pip | Must have static fallback where needed |
| `matplotlib` | yes | Static plots and fallback visualizations | environment / pip | Required for validation-friendly plots |
| `requests` | yes | AFDB and PDBe retrieval | environment / pip | Network calls must be wrapped |
| `networkx` | yes | Multi-edge interaction graph and network metrics | pip | Metrics must be documented and interpretable |
| `molviewspec` | optional | MolViewSpec state generation and inline Mol* rendering or viewer links | pip | Visualization-only; failure must not block data loading, scoring, validation, or interpretation |
| `seaborn` | optional | Statistical/static figure styling | pip | Optional convenience only |
| `freesasa` | optional | SASA/burial factor exploration | pip | Optional; not required for v1 pass |
| `mdtraj` | no for v1 | Trajectory-oriented analysis | n/a | Excluded unless a later spec adds a concrete lightweight use |
| Disallowed v1 | no | MD, exact energy, heavy ML | n/a | No PyTorch, OpenMM, Rosetta, GROMACS, GROMACS-like simulation, training pipeline, or GPU-only workflow in v1 |

MolViewSpec is optional visualization code. Imports must remain local to the visualization section, and no required downstream cells may depend on MolViewSpec-only variables.

## Configuration Defaults

The notebook must expose one `CONFIG` dictionary and display it in a human-readable threshold table.

```python
CONFIG = {
    "structure_source_priority": ["pdbe", "alphafold"],
    "mutation_radius_angstrom": 8.0,
    "hbond_distance_cutoff_angstrom": 3.5,
    "hbond_angle_min_degrees": 120,
    "salt_bridge_cutoff_angstrom": 4.0,
    "hydrophobic_cutoff_angstrom": 5.0,
    "aromatic_cutoff_angstrom": 6.0,
    "disulfide_cutoff_angstrom": 2.2,
    "af_low_confidence_threshold": 50,
    "af_high_confidence_threshold": 80,
    "interaction_weights": {
        "disulfide": 1.00,
        "salt_bridge": 0.75,
        "aromatic": 0.65,
        "hydrogen_bond": 0.60,
        "hydrophobic": 0.50,
        "vdw_packing": 0.40,
        "steric_clash": -1.00,
    },
}
```

## Interaction Configuration

These thresholds are provisional v1 defaults from `protein-bond-type_energy.md`. They are reviewable heuristics, not physical energy definitions.

| Interaction | Detection Rule | Required Evidence | Baseline Weight |
|---|---|---|---:|
| `hydrogen_bond` | donor-acceptor distance <= 3.5 Angstrom and angle >= 120 degrees | donor/acceptor atoms, distance, angle | 0.60 |
| `salt_bridge` | opposite charge groups within <= 4.0 Angstrom | charged groups, distance | 0.75 |
| `hydrophobic` | hydrophobic sidechain centroid distance <= 5.0 Angstrom | residues A, V, I, L, M, F, W, Y; centroid distance | 0.50 |
| `aromatic` | aromatic ring centroid distance <= 6.0 Angstrom | residues F, W, Y, H; centroid distance | 0.65 |
| `disulfide` | cysteine sulfur-sulfur distance <= 2.2 Angstrom | SG atoms and distance | 1.00 |
| `steric_clash` | van der Waals overlap threshold exceeded | atom pair, radii, overlap | -1.00 |
| `vdw_packing` | local neighbor density / packing contact heuristic | neighbor count, local buriedness proxy | 0.40 |

## Interaction Weighting

Interaction weights are relative structural-constraint heuristics. They must not be described as kcal/mol, free energies, thermodynamic measurements, or validated ddG.

```text
relative_weight =
base_interaction_weight
* geometry_quality
* burial_factor
* confidence_factor
* cooperativity_factor
```

Factor ranges:

- `geometry_quality`: 0.0 to 1.0 from distance, angle, and steric validity.
- `burial_factor`: 0.0 to 1.0 approximate solvent exclusion or packing depth.
- `confidence_factor`: 0.0 to 1.0 from AF pLDDT or experimental quality/atom completeness.
- `cooperativity_factor`: 0.0 to 1.0 from local interaction density and network reinforcement.

## Provisional Severity Rubric

Severity is an uncalibrated structural disruption likelihood label. It is useful for triage only and remains non-clinical.

| Contributor | Default Effect |
|---|---|
| New steric clash | strong concern |
| Buried polarity mismatch | strong concern |
| Disulfide or salt-bridge loss | moderate to high concern depending on confidence and burial |
| Aromatic or hydrophobic packing loss | moderate concern, higher if buried or central |
| Cavity or packing reduction | moderate concern |
| Weighted degree or local centrality loss | concern proportional to local role |
| Local fragmentation or redundancy loss | concern proportional to graph impact |
| Secondary-structure incompatibility | moderate to high concern |

Provisional label rule:

- `LOW`: no clash, no buried polarity mismatch, limited low-weight interaction changes, and low network impact.
- `MODERATE`: one or more meaningful interaction/packing/network signals without strong clash or strong buried incompatibility.
- `HIGH`: steric clash, buried polarity mismatch, major central interaction loss, or multiple reinforcing disruption signals.

Domain review and fixture calibration may revise these labels.

## Notebook Sections

| Section | Purpose | Key Outputs | Requirements |
|---|---|---|---|
| 1. Purpose and RUO framing | State audience, research-use scope, and non-clinical limits | Intro markdown, interpretation vocabulary | REQ-013, REQ-019 |
| 2. Setup and configuration | Imports, dependency checks, display settings, source preference, thresholds | `CONFIG`, imports, threshold table | REQ-001, REQ-016, NFR-002 |
| 3. Structure input | Accept AFDB accession, PDBe identifier, or optional local path | `structure_request` | REQ-001 |
| 4. Retrieve structure and metadata | Fetch AFDB/PDBe files and provenance | `structure_payload`, `metadata`, `provenance` | REQ-002, REQ-003, REQ-004 |
| 5. Validate and normalize structure | Parse coordinates, chains, residues, missing atoms, confidence | `residue_table`, `atom_table`, `quality_report` | REQ-005, REQ-015 |
| 6. Mutation input and validation | Capture chain/residue/wild-type/mutant residue and validate mapping | `mutation_spec`, `mutation_validation` | REQ-006 |
| 7. Wild-type interaction detection | Detect interactions and build graph around whole structure and mutation site | `wt_interactions`, `wt_graph`, `wt_local_context` | REQ-007, REQ-008, REQ-010 |
| 8. Network and confidence summaries | Compute local density, centrality, redundancy, and confidence tiers | `wt_network_metrics`, `confidence_summary` | REQ-009, REQ-015 |
| 9. Mutant compatibility estimate | Apply conservative mutation model and estimate local perturbations | `mutant_interactions`, `compatibility_delta` | REQ-011 |
| 10. WT versus mutant comparison | Compare lost/gained interactions, clashes, packing changes, and graph shifts | `comparison_tables`, `network_delta` | REQ-012 |
| 11. Visualizations | Render structure overview, local residue neighborhood, network, perturbation map, and confidence/severity panels | MolViewSpec/Mol* states, Plotly/NetworkX views, static fallback figures | REQ-017 |
| 12. Mechanistic interpretation | Generate observation/interpretation/hypothesis layers and disruption likelihood | `mechanistic_summary`, `severity_summary` | REQ-013, REQ-014, REQ-019 |
| 13. Export artifacts | Write or display tables, figure paths, markdown summary, and provenance bundle | export bundle | REQ-018 |
| 14. Validation snapshot | Show fixture checks, runtime, dependency versions, and open caveats | validation summary | NFR-001, NFR-003 |

Each code section should be preceded by explanatory markdown for structural-biology users. Implementation code should isolate reusable logic in notebook-local functions so later refactoring into `src/insightfold/` is straightforward.

## Data Flow

```text
source preference + accession/path + mutation spec
  -> fetch AFDB/PDBe/local structure and metadata
  -> parse atoms and residues
  -> normalize chains, residue IDs, confidence, missing coordinates
  -> validate mutation mapping
  -> detect wild-type interactions
  -> build residue interaction graph
  -> compute local network and confidence signals
  -> estimate mutant compatibility perturbation
  -> compare WT vs mutant interactions and graph context
  -> visualize evidence
  -> generate observation / interpretation / hypothesis summary
  -> export tables, figures, caveats, provenance
```

## Variable And Artifact Handoff

| Name | Produced By | Consumed By | Type / Shape | Notes |
|---|---|---|---|---|
| `CONFIG` | Setup | All sections | dict | Source preference, thresholds, dependency flags, display options |
| `structure_request` | Structure input | Retrieval | dict | `source`, `identifier`, optional `chain_filter`, optional `local_path` |
| `structure_payload` | Retrieval | Parsing | dict | Raw mmCIF/PDB text, metadata JSON, source URLs |
| `provenance` | Retrieval | All outputs | dict/DataFrame | Source, accession, retrieval date, version, confidence provenance |
| `atom_table` | Parsing | Interaction detection | DataFrame | One row per parsed atom with coordinates and identifiers |
| `residue_table` | Normalization | Mutation validation, interaction detection, visualization | DataFrame | One row per residue with chain, residue ID, residue name, coordinates, confidence |
| `quality_report` | Normalization | Confidence summaries, warnings | DataFrame/dict | Missing residues, missing atoms, low confidence, assembly caveats |
| `mutation_spec` | Mutation input | Mutation model | dict | Chain, residue identifier, expected WT, mutant residue |
| `mutation_validation` | Mutation validation | All mutation sections | dict | Pass/fail, observed residue, confidence around site |
| `wt_interactions` | WT interaction detection | Graph, comparison, export | DataFrame | Residue pair, type, geometry, confidence, evidence text |
| `wt_graph` | Graph construction | Network metrics, visualization | NetworkX graph or adjacency table | Use transparent edge attributes |
| `wt_network_metrics` | Network section | Comparison, summary | DataFrame | Degree, local density, centrality, redundancy, component info |
| `mutant_interactions` | Mutation model | Comparison, summary | DataFrame | Estimated compatible/incompatible local interactions |
| `compatibility_delta` | Mutation model | Summary, severity | DataFrame | Lost, gained, weakened, clash, cavity/packing signals |
| `network_delta` | Comparison | Summary, visualization | DataFrame | Local connectivity and centrality changes |
| `mechanistic_summary` | Interpretation | Export | DataFrame/Markdown | Observation, interpretation, hypothesis records |
| `severity_summary` | Interpretation | Export | DataFrame/dict | Low/moderate/high structural disruption likelihood and rationale |
| `notebook_log` | Setup and all modules | Validation/export | list/DataFrame | Retrieval failures, parsing failures, validation failures, warnings, confidence downgrades, export failures |
| `mvs_states` | Visualization | Export/validation | dict/list | Optional MolViewSpec states for generated 3D views; skipped when `molviewspec` is unavailable |
| `export_bundle` | Export | Validation/review | files or in-memory dict | JSON, CSV, PNG/SVG, graph objects, optional `.mvsj`, markdown, provenance |

## Function Candidates

| Function | Responsibility |
|---|---|
| `fetch_afdb_structure(identifier)` | Return AFDB structure text, confidence metadata, and provenance. |
| `fetch_pdbe_structure(identifier, assembly_policy)` | Return PDBe structure text, experimental metadata, assembly notes, and provenance. |
| `parse_structure(structure_text, format_hint)` | Produce atom-level records from mmCIF/PDB input. |
| `build_residue_table(atom_table, metadata)` | Normalize chain/residue identifiers, residue names, coordinates, and confidence fields. |
| `validate_mutation(residue_table, mutation_spec)` | Confirm chain/residue/wild-type mapping before mutation analysis. |
| `detect_hydrogen_bonds(atom_table, residue_table, config)` | Return hydrogen-bond interaction records and confidence. |
| `detect_salt_bridges(atom_table, residue_table, config)` | Return salt-bridge interaction records and confidence. |
| `detect_hydrophobic_contacts(atom_table, residue_table, config)` | Return hydrophobic contact records and confidence. |
| `detect_aromatic_interactions(atom_table, residue_table, config)` | Return aromatic stacking or edge interaction records and confidence. |
| `detect_disulfides(atom_table, residue_table, config)` | Return disulfide records and confidence. |
| `detect_steric_clashes(atom_table, residue_table, config)` | Return clash records and confidence. |
| `detect_local_packing(atom_table, residue_table, config)` | Return packing-contact records and optional burial proxy. |
| `build_interaction_graph(interactions)` | Build residue-level graph or adjacency tables with edge metadata. |
| `compute_network_metrics(graph, mutation_site)` | Compute local and global graph summaries using documented formulas. |
| `estimate_mutant_compatibility(residue_table, wt_interactions, mutation_spec, config)` | Apply conservative mutation model and return perturbation records. |
| `compare_wt_mutant(wt_interactions, mutant_interactions, wt_metrics)` | Produce interaction and network delta tables. |
| `assign_confidence(record, quality_report, config)` | Assign structure, interaction, and mechanistic confidence. |
| `rank_structural_disruption(comparison_tables, confidence_summary, config)` | Produce provisional structural disruption likelihood. |
| `build_mechanistic_summary(evidence_tables)` | Generate templated observation, interpretation, and hypothesis statements. |

## Parser And Algorithm Choices

| Decision | Rationale | Alternatives Considered | Risk |
|---|---|---|---|
| Use BioPython parser wrappers with mmCIF preferred | AFDB and PDBe support mmCIF and BioPython is the provisional v1 parser from the supporting spec | gemmi, biotite, PDB-only parsing | BioPython may expose residue numbering or assembly edge cases; fixtures must test them. |
| Keep mutation model to sidechain substitution only within 8 Angstrom local radius | PRD forbids exact mutant structures and exact energetics; supporting spec excludes minimization and sampling | Repacking, minimization, MD, full atomistic mutant prediction | Too-simple model may miss effects; output must be framed as compatibility perturbation. |
| Report qualitative compatibility signals | Maintains scientific defensibility without ddG claims | Numeric stability score | Users may want one score; severe overclaiming risk if exact-looking. |
| Use transparent weighted MultiGraph metrics | Supports multiple interaction classes between the same residue pair and preserves interpretability | SimpleGraph, learned graph model | Graph metrics can look authoritative without biological validation. |
| Use template-based explanations | Keeps text grounded in computed observations | Free-form generated narrative | Free-form text risks hallucinated mechanisms. |
| Prefer PDBe when resolution and residue coverage are acceptable | Experimental structures may provide stronger coordinate evidence than AFDB | AFDB-first, manual-only | Exact cutoffs still need domain review; caveats must be shown. |
| Use MolViewSpec/Mol* for 3D molecular views | Matches InsightFold visualization practice and supports portable `.mvsj` states plus Mol* viewer fallbacks | Direct viewer-specific imperative rendering | MolViewSpec is version-sensitive; API verification and graceful fallback are required. |

## Visualization Plan

- Structure overview with chain coloring and source/provenance caption.
- Confidence overlay showing pLDDT for AFDB or experimental quality/missing regions for PDBe.
- Local mutation-site neighborhood with residue labels and interaction edges.
- Interaction-type summary table and stacked bar chart.
- Residue-level network view focused on mutation-site ego network.
- WT versus mutant interaction delta table and perturbation map.
- Compatibility signal panel showing lost interactions, clashes, packing/cavity proxy, and confidence.
- Severity/confidence summary with explicit non-clinical wording.
- MolViewSpec/Mol* views for structure overview, mutation-site focus, confidence overlay, and interaction perturbation.
- Plotly/NetworkX views for weighted network graphs, node centrality scaling, interaction filtering, local neighborhood highlighting, and weighted edge rendering.

MolViewSpec view set:

| View | Generate When | Required Inputs | Fallback |
|---|---|---|---|
| Chain overview | Structure URL or local structure data is available | `structure_source`, `structure_url`, `structure_format`, `chain_ids` | Static structure summary plus skipped-view message |
| Confidence mapping | AFDB pLDDT or PDBe/PDB B-factor data is available | `plddt_per_chain` or `bfactor_per_chain`, residue IDs per chain | Skip confidence view and show reason |
| Mutation-site focus | Mutation validates and residue IDs align | chain ID, residue ID, local residues within configured radius | Static local interaction table |
| Interaction perturbation | WT/mutant comparison has local interaction deltas | residue IDs, interaction deltas, chain IDs | Static perturbation table/plot |

MolViewSpec implementation requirements:

- Verify installed API before generating views: `mvs.create_builder()`, `mvs.ComponentExpression`, `.download(...).parse(...).model_structure()`, `.representation(...)`, `.color(...)`, `builder.get_state().molstar_html()`.
- Confirm accepted parse format strings for the installed version, such as `bcif`, `mmcif`, or `cif`.
- Prefer AFDB `bcifUrl` when available; otherwise use mmCIF.
- Do not assume chain IDs; use parsed chain labels from `residue_table`.
- Guard per-residue component generation: per-residue components are acceptable for <= 500 rendered residues; larger structures should collapse to contiguous residue runs or score bands.
- Wrap the whole MolViewSpec section in `try/except ImportError`, and wrap each view in its own `try/except Exception`.
- Emit a generated/skipped view summary and append skipped views to `notebook_log`.
- Display inline with Mol* when possible; otherwise emit a Mol* viewer URL fallback.
- Export `.mvsj` states only when requested by the export section and when state generation succeeds.

Interaction colors:

| Interaction | Color |
|---|---|
| hydrogen bond | blue |
| salt bridge | red |
| hydrophobic | gold |
| aromatic | purple |
| disulfide | green |
| steric clash | orange |

## Error Handling And Graceful Degradation

- Retrieval failures stop retrieval-dependent sections and show source URLs plus fixture/cache guidance.
- Parser failures show file format, source, and first actionable validation issue.
- Missing optional confidence metadata allows analysis only with lower confidence and visible provenance note.
- Mutation validation failures stop mutation sections but still allow structure and WT interaction inspection where possible.
- Unsupported biological contexts emit warnings and suppress high-confidence mechanistic claims.
- Optional MolViewSpec/Mol* visualization failures do not block tables, static plots, scoring, validation, or summaries.
- All warnings and failures are appended to `notebook_log`.

## Mandatory Versus Optional Sections

Mandatory for MVP:

- RUO framing, setup, source input, retrieval, parsing, mutation validation, WT interaction detection, mutation perturbation comparison, confidence summary, interpretation, validation snapshot.

Optional or stretch:

- Local file upload, multi-structure comparison, batch mutations, additional MolViewSpec diagnostic views beyond the core view set, and calibrated disruption thresholds.

## Implementation Assumptions

- Proceed with AFDB and PDBe as primary sources; local files are optional.
- Use one mutation at a time for v1.
- Use transparent rule-based interaction detection and sidechain-substitution-only compatibility perturbation.
- Use provisional qualitative severity only if calibration remains unresolved, and label it as uncalibrated.
- Do not use conservation, clinical evidence, or genomics evidence in v1 scoring.
- Use JSON, CSV, PNG/SVG, graph-object, optional `.mvsj`, markdown, and provenance exports where the implementation can generate them without hidden state.
