# NMR Restraints Visualization Notebook Spec Pack

Source PRD: `specs/nmr_restraints/nmr_restraints_prd.md`

## Metadata

| Field | Value |
|---|---|
| Spec ID | `nmr-restraints-visualization-notebook` |
| Product | NMR Restraints Visualization Notebook |
| Lifecycle Stage | Scope / notebook specification |
| Status | Draft |
| Owner | InsightFold maintainers / TBD |
| Last Updated | 2026-05-08 |
| Target Runtime | Local Jupyter/JupyterLab Python 3.11+; Colab desirable but secondary |
| Source PRD | `specs/nmr_restraints/nmr_restraints_prd.md` |
| Related Source Spec | `specs/nmr_restraints/nmr_restraints_visualization_complete_specification_v2.md` |

## Pack Contents

| File | Purpose |
|---|---|
| `requirements.md` | Observable notebook requirements, user stories, non-goals, edge cases, and open questions. |
| `notebook-design.md` | Section-by-section notebook architecture, variable handoffs, dependencies, algorithms, visualization design, and degradation behavior. |
| `data-contracts.md` | External source, input, parsed table, computed output, visualization state, provenance, and cache contracts. |
| `fixture-manifest.md` | Provisional fixtures, required expected snapshots, tolerances, and fixture-selection blockers. |
| `tasks.md` | Atomic implementation, validation, documentation, and review tasks for an implementation agent. |
| `validation.md` | Execution, fixture, algorithm, visualization, hidden-state, runtime, and scientific-language validation plan. |
| `docs-plan.md` | Tutorial, how-to, reference, explanation, export, and maintainer documentation plan. |

## Implementation Entry Point

Start with:

1. `requirements.md` for acceptance criteria and open questions.
2. `fixture-manifest.md` to pin or confirm fixtures and expected snapshots.
3. `notebook-design.md` for architecture and function boundaries.
4. `tasks.md` for implementation sequencing.

Do not begin final validation until fixture expected snapshots are recorded.

## Blocking Items Before Scientific Sign-Off

- Confirm `9L1V` or replace it as the happy-path fixture after parsing and mapping snapshots are observed.
- Add fixture coverage for ambiguous `OR` restraints, wrapped dihedrals, missing restraints, and mapping threshold behavior.
- Confirm `pynmrstar` suitability or document parser replacement.
- Define and complete NMR domain review of formulas, terminology, ambiguity behavior, and interpretation caveats.

## Core Assumptions

- PDBe entry-file endpoints are the canonical v1 remote retrieval source.
- RCSB mmCIF endpoint remains a model-file fallback, not the primary restraint source.
- Gemmi is the canonical coordinate parser.
- `pynmrstar` is the provisional NMR-STAR parser.
- V1 analyzes a selected model index, defaulting to `0`, and does not perform ensemble-wide analysis.
- V1 supports explicit `OR` ambiguous distance restraints using smallest-distance selected-member behavior.
- Mapping coverage is computed over logical restraints; an evaluable `OR` ambiguous group counts as `1/1` mapped even if only one candidate member maps.
- V1 distance restraints are expected to include both upper and lower bounds; missing-bound rows are unsupported and reported in diagnostics.
- Python notebook state is authoritative; Mol* frontend state is regenerated from notebook variables.
- MolViewSpec/Mol* views are embedded inline into notebook cells using the `state.molstar_html()` plus base64 `IFrame` pattern from `notebooks/homodimer_diagnostic.ipynb`.
