# Product Requirements Document (PRD)
## NMR Restraints Visualization Notebook

---

# 1. Document Overview

| Field | Details |
|---|---|
| Product Name | NMR Restraints Visualization Notebook |
| Version | v1 Draft |
| Date | 2026-05-08 |
| Author(s) | Codex, based on `specs/nmr_restraints/nmr_restraints_visualization_complete_specification_v2.md` |
| Stakeholders | Structural biologists, protein biochemists, computational biology researchers, PDB users, NMR structure interpreters, InsightFold maintainers |
| Status | Draft |
| Related Documents | `specs/nmr_restraints/nmr_restraints_visualization_complete_specification_v2.md` |
| Regulatory Scope | RUO; non-clinical; not intended for diagnostic, regulatory, clinical, or automated structure-validation decisions |

---

# 2. Executive Summary

## Problem Statement

Researchers inspecting NMR-derived macromolecular structures need to understand how deposited experimental restraints support, constrain, or locally conflict with the deposited coordinate model. Existing workflows often require specialist NMR tooling, manual restraint parsing, and separate structure viewers, making it difficult for non-specialist structural biology users to answer:

> Which regions of an NMR structure are experimentally supported, weakly constrained, or locally inconsistent with the deposited NMR restraints?

The notebook should transform deposited NMR restraint and coordinate data into an interpretable spatial evidence model. It must help users explore local experimental support while avoiding misleading validation claims.

---

## Product Vision

Create an interactive, reproducible InsightFold notebook that retrieves or accepts NMR structure and restraint files, maps restraints onto deposited coordinates, computes residue-level restraint density and restraint violations, and renders the evidence in synchronized tables, sequence controls, and MolViewSpec/Mol* 3D views.

The notebook should occupy a clear niche between raw restraint files, general molecular viewers, and formal validation engines: an exploratory evidence viewer that makes deposited NMR support visible and understandable.

---

## Success Criteria

- Users can provide a PDB ID or local model/restraint files and generate a reproducible exploratory report.
- The notebook retrieves or loads mmCIF structure coordinates and NMR-STAR restraint files.
- Distance and dihedral restraints are parsed, mapped, and summarized with clear mapping coverage diagnostics.
- Residue restraint density is computed and visualized on the structure.
- Distance and dihedral violations are computed and displayed in sortable/filterable outputs.
- Users can move from a global evidence view into local residue-level evidence inspection.
- Ambiguous `OR` distance restraints are handled using the specified smallest-distance evaluation semantics.
- The notebook fails gracefully for missing files, malformed files, unsupported records, or low mapping coverage.
- User-facing text consistently frames density as experimental coverage and violations as local inconsistency, not correctness or failure.
- Typical structures up to about 500 residues and 10,000 restraints run in under 30 seconds on commodity hardware.

---

# 3. Background & Context

## Scientific Context

- Biological domain: NMR-derived macromolecular structure interpretation.
- Core object: deposited structure coordinates plus deposited NMR restraints.
- Structural model type: NMR structures are inferred by optimizing models against experimental constraints rather than directly observing a single final coordinate set.
- Restraint types in v1: distance restraints, usually NOE-derived, and dihedral restraints.
- Core interpretation need: relate restraints, violations, density, and local geometry in a single spatial context.
- Scientific limit: restraint density is a proxy for experimental coverage, not a confidence score; violations are evidence of inconsistency, not automatic errors.

---

## Business / Institutional Context

- Aligns with InsightFold's notebook-driven lifecycle for reproducible structural biology analysis.
- Supports PDB and PDBe users who need transparent interpretation of deposited NMR evidence.
- Provides a foundation for later NMR evidence visualization in web applications, automated reports, or batch workflows.
- Complements, but does not replace, formal validation systems, refinement software, or expert NMR analysis tools.

---

## Existing Solutions

| Solution | Strengths | Weaknesses |
|---|---|---|
| PDB/PDBe structure pages | Provide authoritative deposited structures and metadata | Restraint evidence is not always exposed in an integrated local 3D evidence workflow |
| General molecular viewers such as Mol*, PyMOL, and ChimeraX | Strong 3D inspection and structure navigation | Do not automatically parse NMR-STAR restraints into evidence density and violation views |
| NMR refinement and validation software | Domain-specific restraint handling and refinement support | Often specialist-oriented, workflow-heavy, or not optimized for exploratory notebook teaching |
| Raw NMR-STAR restraint files | Preserve deposited experimental information | Difficult for non-specialists to interpret spatially without parsing and mapping |
| Aggregate validation reports | Useful summary diagnostics | Can hide residue-level evidence networks and local restraint context |

---

# 4. Goals & Objectives

## Product Goals

| Goal | Metric |
|---|---|
| Make deposited NMR restraint evidence interpretable | Users can identify dense, sparse, violated, and locally inconsistent regions from notebook outputs |
| Support reproducible exploratory analysis | A PDB ID or local file pair produces deterministic tables, diagnostics, and visualization states |
| Reduce specialist barrier to NMR restraint interpretation | Non-specialist structure users can inspect restraint density and violations without raw NMR-STAR expertise |
| Maintain scientific caution | Outputs avoid validation verdicts, confidence scores, and correctness labels |
| Prepare for implementation-ready notebook specification | PRD captures inputs, outputs, non-goals, data needs, risks, and open questions for downstream spec work |

---

## Non-Goals

Clearly define what this product will **not** address.

- It will not produce automated quality scores, global confidence metrics, pass/fail labels, or validation verdicts.
- It will not perform NMR refinement, restraint recalculation, model correction, or energy minimization.
- It will not perform ensemble-wide analysis in v1.
- It will not implement weighted ambiguity semantics, pseudoatom handling, dynamic averaging semantics, or full NMR assignment logic in v1.
- It will not compare multiple structures or run cross-entry benchmarking in v1.
- It will not provide clinical, diagnostic, regulatory, or production validation output.
- It will not implement collaborative sharing, authentication, RBAC, or platform-scale user management.

---

# 5. Target Users

## Primary Users

| User Type | Goals | Pain Points |
|---|---|---|
| Structural biologists | Inspect where deposited restraints support or conflict with an NMR structure | Restraint evidence is often separated from 3D structural context |
| Protein biochemists | Interpret local experimental support around residues, regions, or functional sites | May understand structures but not detailed NMR restraint semantics |
| Computational biology researchers | Generate reproducible restraint-density and violation analyses | Need parsable outputs, diagnostics, and deterministic notebooks |
| PDB users | Explore deposited NMR evidence for a public structure | Raw restraint files and specialist tooling create friction |
| InsightFold maintainers | Build and validate a scoped notebook workflow | Need explicit data contracts, failure modes, fixtures, and scientific caveats |

---

## Secondary Users

- Educators teaching NMR structure interpretation.
- Research collaborators reviewing local structural evidence.
- PDBe/PDB support or documentation teams.
- Developers building future structural evidence viewers.

---

# 6. User Personas

## Persona 1 - Structural Biology Researcher

### Background
- Understands macromolecular structure, residue-level interpretation, and experimental uncertainty.
- Uses PDB structures and molecular viewers, but may not routinely parse NMR-STAR restraint files.

### Needs
- Fast visibility into restraint-supported and weakly restrained regions.
- Local inspection of violated restraints in 3D context.
- Exportable tables and figures for research discussions.

### Frustrations
- Raw restraint files are hard to interpret spatially.
- Aggregate validation summaries can obscure local evidence patterns.
- Manual mapping between restraints and atoms is slow and error-prone.

---

## Persona 2 - Protein Biochemist or Functional Researcher

### Background
- Understands residues, domains, and protein function.
- Uses deposited structures to reason about mechanisms or experimental design.

### Needs
- Plain-language explanations of restraint density and violations.
- Clear local views around residues of interest.
- Guardrails that prevent overinterpreting lack of violations as correctness.

### Frustrations
- NMR-specific terminology can be opaque.
- Existing viewers rarely explain what deposited restraints imply locally.
- Visual clutter can make dense restraint networks unreadable.

---

## Persona 3 - InsightFold Notebook Developer

### Background
- Builds reproducible notebooks from PRDs, specs, fixtures, and validation criteria.
- Needs explicit functional boundaries and data assumptions.

### Needs
- Clear input modes, data sources, parsing rules, and failure conditions.
- Defined output tables, visualization states, and warning messages.
- Testable success criteria for downstream implementation and review.

### Frustrations
- NMR restraint semantics contain edge cases that can lead to hidden implementation assumptions.
- Visualization notebooks can become fragile if state ownership is unclear.

---

# 7. Product Scope

## In Scope

- Remote retrieval by PDB ID.
- Local-file mode for user-provided mmCIF model files and NMR-STAR `.str` restraint files.
- Structure parsing from mmCIF using Gemmi.
- NMR-STAR parsing for compatible distance and dihedral restraint categories.
- Exact atom mapping using deposited `auth_*` identifiers.
- Mapping coverage diagnostics and threshold-based continuation or abort behavior.
- Geometry calculations for distances and dihedral angles.
- Violation magnitude computation for distance and dihedral restraints.
- `OR` ambiguous distance restraint evaluation using smallest candidate distance.
- Residue-level absolute and normalized restraint density computation.
- MolViewSpec/Mol* global structure visualization colored by restraint density.
- Local evidence view around a selected residue.
- Distance and dihedral restraint overlays with clutter controls.
- Sequence view, residue navigation widget, and violation table interaction.
- CSV export for violation tables and notebook HTML export.
- Safe parsing practices and clear user-facing warnings.

---

## Out of Scope

- Clinical diagnosis, clinical reporting, or regulated structural validation.
- Automated structure correctness scoring.
- NMR refinement, re-refinement, or structure repair.
- Ensemble-wide statistics or conformer comparison in v1.
- Pseudoatom expansion, weighted ambiguous restraints, dynamic averaging, and full assignment-resolution semantics in v1.
- Machine-learning interpretation.
- Cross-structure comparison and batch workflows.
- Formal MolViewSpec state interchange schema.
- Custom Mol* camera choreography beyond standard viewer behavior.

---

# 8. Use Cases

| ID | Use Case | Priority |
|---|---|---|
| UC-01 | Enter a PDB ID and retrieve structure plus NMR restraint files | High |
| UC-02 | Load local mmCIF and NMR-STAR files for offline or custom analysis | High |
| UC-03 | Parse and map distance and dihedral restraints onto deposited coordinates | High |
| UC-04 | Review mapping coverage, ignored records, and unsupported restraint warnings | High |
| UC-05 | Compute residue-level restraint density and inspect global coverage patterns | High |
| UC-06 | Compute and prioritize distance and dihedral restraint violations | High |
| UC-07 | Select a residue and inspect all local restraints within spatial context | High |
| UC-08 | Select a violation row and center local visualization on involved atoms/residues | High |
| UC-09 | Export violation tables or notebook HTML for communication and review | Medium |
| UC-10 | Use explanatory copy to avoid interpreting density or violations as validation verdicts | High |

---

# 9. User Stories

## User Story Template

> As a [user type], I want to [action], so that [benefit].

---

## Example

> As a structural biologist, I want to enter a PDB ID and see restraint density mapped onto the structure, so that I can quickly identify experimentally dense and sparse regions.

### Acceptance Criteria

- The notebook accepts and normalizes a four-character PDB ID.
- Structure and restraint files are retrieved or a clear diagnostic is shown.
- Residue density is computed from mapped restraints.
- The initial 3D view shows the full structure colored by normalized density.
- The notebook explains that density is experimental coverage, not confidence.

---

## Additional User Stories

> As a protein biochemist, I want to select a residue and see its local restraint network, so that I can understand whether a functional region is directly supported by deposited evidence.

> As a PDB user, I want violated restraints listed and highlighted spatially, so that I can inspect local inconsistencies without reading raw restraint files.

> As a notebook developer, I want unsupported records and mapping failures reported explicitly, so that partial analyses do not appear complete.

> As an educator, I want no-violation outputs to include caveats, so that students do not mistake the absence of highlighted violations for structural correctness.

---

# 10. Functional Requirements

| ID | Requirement | Priority |
|---|---|---|
| FR-01 | Notebook shall support remote retrieval using a validated four-character PDB ID | High |
| FR-02 | Notebook shall support local mmCIF model files and local NMR-STAR restraint files | High |
| FR-03 | PDB IDs shall be normalized to uppercase before retrieval | High |
| FR-04 | Notebook shall retrieve or load structure coordinates in mmCIF-compatible form | High |
| FR-05 | Notebook shall parse structure coordinates with Gemmi and build canonical atom lookup tables | High |
| FR-06 | Notebook shall parse compatible distance restraint and dihedral restraint categories from NMR-STAR files | High |
| FR-07 | Notebook shall merge compatible restraint loops, saveframes, and tables into unified collections | High |
| FR-08 | Notebook shall map restraint atoms to coordinate atoms using exact deposited `auth_*` identifiers | High |
| FR-09 | Notebook shall continue normally when mapping coverage is greater than 95%, warn between 70% and 95%, and abort below 70% | High |
| FR-10 | Notebook shall ignore alternate conformations outside primary/default, blank, or `.` altloc records in v1 | Medium |
| FR-11 | Notebook shall compute Euclidean distances for distance restraints | High |
| FR-12 | Notebook shall compute dihedral angles normalized to `[-180, +180]` degrees | High |
| FR-13 | Notebook shall compute distance and dihedral violation magnitudes from allowed bounds | High |
| FR-14 | Notebook shall handle wrapped dihedral intervals using circular angular geometry | High |
| FR-15 | Notebook shall identify `OR` ambiguous distance restraints using shared IDs and `Member_logic_code = OR` | High |
| FR-16 | Notebook shall evaluate ambiguous distance restraints by selecting the smallest measured candidate distance | High |
| FR-17 | Notebook shall compute absolute and normalized residue restraint density | High |
| FR-18 | Notebook shall produce residue density and violation tables with documented fields | High |
| FR-19 | Notebook shall sort violations by descending magnitude, restraint type, and residue number by default | Medium |
| FR-20 | Notebook shall render an initial global MolViewSpec/Mol* structure view colored by normalized restraint density | High |
| FR-21 | Notebook shall allow users to toggle normalized and absolute density display | Medium |
| FR-22 | Notebook shall render only violated restraints in the global view to reduce clutter | High |
| FR-23 | Notebook shall cap globally visible restraints using a configurable `max_visible_restraints` value | Medium |
| FR-24 | Notebook shall provide local evidence inspection for a selected residue using a configurable spatial radius | High |
| FR-25 | Local evidence view shall display satisfied, violated, and ambiguous restraints involving the selected residue | High |
| FR-26 | Notebook shall separate distance and dihedral views by default to reduce visual clutter | Medium |
| FR-27 | Notebook shall provide a sequence or residue navigation control synchronized with 3D and tabular outputs | High |
| FR-28 | Notebook shall provide sortable, searchable, filterable violation inspection where feasible in notebook UI | Medium |
| FR-29 | Notebook shall isolate visualization panel failures so one failed view does not prevent other views from rendering | High |
| FR-30 | Notebook shall display scientifically cautious explanatory text for density, violations, no-violation states, and mapping warnings | High |
| FR-31 | Notebook shall optionally cache downloaded files, parsed structures, parsed restraints, and computed tables | Medium |
| FR-32 | Notebook shall export violation tables as CSV and support notebook HTML export | Medium |

---

# 11. Non-Functional Requirements

## Performance

- Target runtime: under 30 seconds for typical structures on commodity hardware.
- Target structure size: up to about 500 residues.
- Target restraint count: up to about 10,000 restraints.
- Visualization state generation should avoid excessive per-residue component creation for structures above 500 residues by grouping contiguous ranges where possible.

---

## Scalability

- V1 is scoped to one structure at a time.
- Batch processing, multi-entry comparison, and service-scale workloads are deferred.
- Internal module boundaries should permit later reuse in batch pipelines or web applications.

---

## Reliability

- Network retrieval should retry up to three times with exponential backoff.
- Hard failures should include structure unavailable, restraint file unavailable, malformed coordinate file, no parseable restraints, and mapping coverage below threshold.
- Warning conditions should include partially unmapped restraints, unsupported ambiguity semantics, missing restraint categories, ignored alternate conformations, and unsupported records.
- Partial visualization failures should be isolated and reported without collapsing the whole notebook.

---

## Security

- The notebook should avoid arbitrary code execution from downloaded files.
- Parsers should treat remote structure and restraint files as data only.
- No authentication is required for public PDB/PDBe retrieval.
- Local-file mode should not execute file contents.

---

## Compliance

- RUO only.
- No clinical, diagnostic, regulated validation, GxP, HIPAA, MDR, or IVDR use is claimed.
- Public structure and restraint retrieval should preserve source provenance.

---

# 12. Data Requirements

## Data Types

| Data Type | Details |
|---|---|
| Structure coordinates | Atom records, residues, chains, coordinates, altloc information, model/conformer selection |
| Distance restraints | Atom pair identifiers, lower bound, upper bound, ambiguity metadata |
| Dihedral restraints | Four atom identifiers, lower bound, upper bound |
| Computed geometry | Distances, torsion angles, selected ambiguous members |
| Computed evidence summaries | Residue density, violation tables, mapping statistics |
| Visualization states | Deterministic MolViewSpec state definitions generated from notebook state |

---

## Supported Formats

| Data Type | Formats |
|---|---|
| Structure model | mmCIF; Binary CIF may be used for remote MolViewSpec visualization when available |
| Restraints | NMR-STAR `.str` files |
| Tables | pandas DataFrames; CSV export for violation tables |
| Report | Jupyter notebook and HTML export |

---

## Metadata Standards

- PDB/mmCIF atom identifiers.
- NMR-STAR restraint categories.
- Deposited `auth_asym_id`, `auth_seq_id`, `auth_comp_id`, and `auth_atom_id` fields for canonical mapping.
- PDB or PDBe entry-file provenance for retrieved inputs.

---

# 13. Architecture Overview

## System Components

- Retrieval module for remote model and restraint downloads.
- Local input module for user-provided model and restraint files.
- Gemmi-based structure parser.
- NMR-STAR restraint parser.
- Atom mapping module.
- Geometry engine for distances, dihedrals, and violations.
- Density engine for residue-level restraint coverage.
- Visualization state builder using MolViewSpec.
- Notebook UI layer for widgets, tables, sequence navigation, filtering, and rendering.
- Diagnostics and logging layer for parsing, mapping, warnings, and runtime reporting.

---

## Integrations

| System | Purpose |
|---|---|
| PDBe entry-file endpoints | Remote retrieval of mmCIF model files and NMR restraint files |
| RCSB PDB model endpoints | Optional/fallback structure coordinate retrieval |
| BMRB/PDB NMR-STAR data | Deposited NMR restraint source |
| Gemmi | Canonical coordinate parsing |
| MolViewSpec 1.8.1 | Declarative reproducible molecular visualization state |
| Mol* | Interactive molecular rendering |
| JupyterLab / notebook environment | User execution, widgets, report generation |

---

# 14. Workflow & UX

## Key Workflows

1. User selects remote retrieval or local-file mode.
2. User provides a PDB ID or local model/restraint file paths.
3. Notebook validates inputs and retrieves or loads data.
4. Notebook parses coordinates and restraints.
5. Notebook maps restraint atoms and reports mapping coverage.
6. Notebook computes geometry, violations, and residue density.
7. Notebook displays a global density-colored structure and violation summaries.
8. User selects residues, sequence regions, or violation rows.
9. Notebook updates the local evidence view and synchronized tables.
10. User exports tables or notebook output for review.

---

## UX Requirements

- Initial view should immediately show global restraint-support distribution with low visual clutter.
- Density and violation explanations should appear near the outputs they qualify.
- Users should be able to move from global density overview to local residue inspection without editing code.
- Tables should remain sortable/searchable/filterable where notebook tooling permits.
- Visual encodings must not rely exclusively on color; saturation, borders, labels, tooltips, or table values should supplement color.
- No-violation and missing-restraint states must be explicit, non-alarming, and scientifically cautious.
- Controls should include PDB ID input, local model file input, local restraint file input, residue selector, local radius, violation thresholds, density mode, and restraint visibility limits.

---

# 15. AI / ML Components (Optional)

## AI Use Cases

- None in v1.

---

## Model Requirements

| Requirement | Details |
|---|---|
| Explainability | Not applicable; no AI/ML model is planned |
| Validation datasets | Not applicable for AI/ML; notebook validation should use pinned structural fixtures |
| Monitoring | Not applicable |
| Bias mitigation | Not applicable |

---

# 16. Risks & Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Users interpret density as structural correctness | High | Use repeated explanatory copy that density is experimental coverage, not confidence |
| Users interpret violations as errors | High | Frame violations as local inconsistency that may reflect ambiguity, heterogeneity, dynamics, or refinement tradeoffs |
| Missing or incompatible restraint files prevent analysis | High | Provide graceful empty states and local-file mode |
| Atom identifier mismatches reduce mapping coverage | High | Report mapping statistics, warn for partial mapping, abort below threshold |
| Ambiguous restraint semantics are oversimplified | Medium | Restrict v1 to explicit `OR` smallest-distance behavior and warn/exclude unsupported semantics |
| Visualization clutter hides important evidence | Medium | Render only violated restraints globally, provide local view for full restraint context, cap visible restraints |
| Wrapped dihedral intervals are mishandled | Medium | Require circular angular geometry and fixture validation |
| Remote endpoint availability changes | Medium | Support local-file mode and document retrieval provenance |
| Notebook state diverges from Mol* frontend state | Medium | Treat Python notebook state as authoritative and regenerate deterministic MolViewSpec states |
| Large structures or dense restraint sets exceed runtime targets | Medium | Use caching, visibility limits, and grouped visualization ranges |

---

# 17. Dependencies

- Public PDB/PDBe structure and restraint availability.
- Gemmi for coordinate parsing.
- `pynmrstar` or equivalent NMR-STAR parsing support.
- `numpy` and `pandas` for geometry and tabular outputs.
- MolViewSpec 1.8.1 and Mol* for interactive visualization.
- JupyterLab or compatible notebook runtime.
- Optional: `ipywidgets`, `plotly`, and local cache storage.
- Domain review by an NMR-aware structural biologist before release beyond internal draft use.

---

# 18. Milestones & Timeline

| Milestone | Owner | Target Date |
|---|---|---|
| PRD draft complete | Codex / InsightFold maintainers | 2026-05-08 |
| Notebook spec pack generated | TBD | TBD |
| Fixture selection complete | TBD | TBD |
| MVP notebook implemented | TBD | TBD |
| Execution validation complete | TBD | TBD |
| Domain/scientific review complete | TBD | TBD |

---

# 19. KPIs & Success Metrics

| Metric | Target |
|---|---|
| Input success for supported fixtures | 100% for selected remote and local-file fixtures |
| Mapping diagnostics present | 100% of runs report parsed counts, mapped counts, excluded counts, and mapping coverage |
| Runtime | Under 30 seconds for typical supported structures |
| Visualization availability | Global density view renders for all successful fixture runs |
| Local inspection availability | Local view updates from residue or violation selection for all supported fixtures |
| Scientific copy compliance | No output describes density as confidence or violations as automatic errors |
| Export availability | CSV violation export and notebook HTML export available for successful runs |
| Graceful failure coverage | Missing structure, missing restraints, malformed inputs, no parseable restraints, and low mapping coverage produce clear diagnostics |

---

# 20. Open Questions

- Which exact PDBe/RCSB retrieval source should be canonical when both endpoint families are available?
- Which NMR-STAR parser should be canonical if `pynmrstar` cannot robustly handle all required deposited files?
- What fixture PDB entries should cover normal execution, partial mapping, ambiguous restraints, dihedral restraints, no violations, and missing/unsupported restraint cases?
- Should v1 support a single selected model/conformer only, or provide clearer UI for multiple NMR models while still avoiding ensemble-wide analysis?
- How should ligand, nucleic-acid, or non-protein restraints be handled in v1?
- What exact thresholds define low, medium, and high distance or dihedral violation colors?
- Should local-file mode use file path inputs, upload widgets, or both across target notebook environments?
- What level of accessibility is feasible for sequence and 3D synchronized views inside Jupyter?
- Should satisfied restraints be exportable even though the primary violation table emphasizes violated restraints?
- What human/domain review checklist is required before considering the notebook ready for public educational use?

---

# Appendix

## Glossary

| Term | Definition |
|---|---|
| NMR | Nuclear magnetic resonance, an experimental method used to infer structural models from constraints and measurements |
| Restraint | A deposited experimental constraint used during NMR structure determination or refinement |
| Distance restraint | A bound on expected spatial distance between two atoms |
| Dihedral restraint | A bound on a torsion angle defined by four atoms |
| Violation | A measured geometry outside the deposited allowed range |
| Restraint density | Number of restraints associated with a residue, used here as a proxy for experimental coverage |
| NMR-STAR | Text format used for NMR data and restraint deposition |
| mmCIF | Macromolecular Crystallographic Information File format used for deposited structural models |
| MolViewSpec | Declarative specification format for reproducible molecular visualization states |
| Mol* | Interactive molecular visualization engine |
| RUO | Research Use Only |

---

## References

- `specs/nmr_restraints/nmr_restraints_visualization_complete_specification_v2.md`
- PDBe entry-file model endpoint: `https://www.ebi.ac.uk/pdbe/entry-files/download/{PDB_ID}.cif`
- PDBe NMR restraint endpoint: `https://www.ebi.ac.uk/pdbe/entry-files/download/{PDB_ID}_nmr-data.str`
- RCSB model endpoints: `https://models.rcsb.org/{PDB_ID}.bcif`; `https://files.rcsb.org/download/{PDB_ID}.cif`
- Gemmi documentation: `https://gemmi.readthedocs.io/en/latest/mol.html`
- MolViewSpec 1.8.1 documentation and Mol* viewer ecosystem
