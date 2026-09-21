# Notebook UX Contract

Source PRD: `specs/nmr_restraints/nmr_restraints_prd.md`

## Primary User Question

Which regions of this deposited NMR structure are experimentally supported, weakly constrained, or locally inconsistent with the deposited NMR restraints?

## Target User And Expertise

- Primary users: structural biologists, protein biochemists, computational biology researchers, and PDB users inspecting deposited NMR evidence.
- Assumed expertise: users understand residues, chains, and 3D structure inspection, but may not understand NMR-STAR files or ambiguity semantics.
- Notebook stance: exploratory research-use evidence viewer only. It is not a refinement workflow, not a formal validation engine, and not a mutation or AlphaFold analysis notebook.

## First Runnable User Input Cell

The first runnable cell must define one visible configuration object plus a concise mode selector. Defaults stay pinned to the happy-path fixture but the cell represents the real user workflow, not fixture-only execution.

| Parameter | Default | User Editable? | Validation | Notes |
|---|---|---|---|---|
| `input_mode` | `remote` | yes | one of `remote`, `local` | Remote is the primary v1 story; local supports parity and offline use. |
| `pdb_id` | `9L1V` | yes | `^[A-Za-z0-9]{4}$` | Remote mode default only. |
| `model_index` | `0` | yes | integer `>= 0`; model must exist | Single-model analysis only in v1. |
| `model_path` | empty | yes in local mode | existing `.cif` or `.mmcif` | Local mode only. |
| `restraint_path` | empty | yes in local mode | existing `.str` | Local mode only. |
| `violation_threshold_distance` | `0.5` | yes | float `>= 0` | Display threshold only; not a pass/fail rule. |
| `violation_threshold_dihedral` | `5.0` | yes | float `>= 0`, `<= 180` | Display threshold only. |
| `local_context_radius` | `6.0` | yes | float `> 0` | Defines neighborhood for local evidence view. |
| `max_visible_restraints` | `250` | yes | integer `>= 0` | Caps global overlay clutter. |
| `cache_enabled` | `true` | yes | boolean | Raw downloaded fixture files remain under repo-local gitignored cache paths. |

## Happy-Path Default Example

Default runnable example:

- `input_mode = "remote"`
- `pdb_id = "9L1V"`
- `model_index = 0`
- thresholds and radius left at defaults

Expected first visible confirmation that the notebook is working:

1. a dependency/config summary renders without errors
2. model and restraint provenance appears with recorded source URLs or local paths
3. parser counts for the selected model and the preserved `9L1V` restraint loops appear
4. an inline MolViewSpec/Mol* global view renders through the `state.molstar_html()` plus base64 `IFrame` embedding pattern, or the notebook shows the documented fallback warning and still continues with tables

## User Flow

| Stage | User Action | System Action | User Sees | Interpretation Need |
|---|---|---|---|---|
| 1 | Set remote PDB ID or local file paths | Validate input mode and parameters | Readable config summary and any early validation errors | Explain accepted formats and selected model scope. |
| 2 | Run retrieval/load step | Fetch or load mmCIF and NMR-STAR files, record provenance, consult repo-local cache if enabled | Provenance table and source status | Make clear whether data came from remote endpoints, local files, or cache. |
| 3 | Run parsing step | Parse structure, distance restraints, dihedral restraints, and unsupported records | Counts, warnings, parser diagnostics | State that unsupported records are excluded explicitly, not silently ignored. |
| 4 | Run mapping step | Map restraint atoms by exact `auth_*` identity and compute logical-restraint coverage | Mapping report with threshold decision | Explain that mapping coverage is counted at logical-restraint level. |
| 5 | Run geometry and density step | Evaluate distances, dihedrals, ambiguity, violations, and residue density | Violation table and residue density table | State that density is coverage, not confidence; violations are local inconsistency, not automatic error. |
| 6 | Inspect global evidence | Render density-colored structure and thresholded violated restraints | Inline Mol* global view | Global view is overview only; dense overlays are capped for readability. |
| 7 | Select residue or violation | Update one authoritative Python `selection_state` and regenerate local tables/views | Local evidence table and inline local Mol* view | Local context is the main place for residue-level interpretation. |
| 8 | Export results | Write CSV and optional JSON/state outputs with provenance | Export manifest | Remind user that outputs preserve RUO framing and provenance. |

## Fixtures Versus User Flow

- Fixtures validate the notebook; they do not replace the user workflow.
- `happy-path-9l1v` is the pinned default example because it exercises real distance, dihedral, and `OR` ambiguous restraints.
- `local-copy-9l1v` validates parity for local-file mode.
- `edge-or-9l1v`, `edge-dihedral-9l1v`, `negative-missing-restraints-1crn`, `edge-wrapped-dihedral-synthetic`, and `edge-low-mapping-synthetic` exist to validate preserved v1 decisions and failure behavior.
- Mutation fixtures and AlphaFold examples remain out of scope and must not be introduced as alternate default stories in this notebook.

## Trust, Caveats, And Limitations To Show In Notebook

- State that restraint density reflects deposited experimental coverage, not structure confidence or correctness.
- State that violations reflect local inconsistency with the selected model and bounds, not automatic errors.
- State that v1 distance restraints are expected to have both lower and upper bounds; missing-bound rows are unsupported and excluded with diagnostics.
- State that ambiguous `OR` distance restraints use smallest measured candidate distance in v1.
- State that mapping coverage is counted at logical-restraint level and low coverage can make results partial or abort analysis.
- State that MolViewSpec/Mol* views are a visualization layer over authoritative notebook state.
- State that local cached raw fixture files live only under repo-local paths such as `specs/nmr_restraints/fixtures/cache/` and that cache stays gitignored.
