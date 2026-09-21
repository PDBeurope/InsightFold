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
| Last Updated | 2026-05-11 |
| Target Runtime | Local Jupyter/JupyterLab Python 3.11+; Colab desirable but secondary |
| Source PRD | `specs/nmr_restraints/nmr_restraints_prd.md` |
| Related Source Spec | `specs/nmr_restraints/nmr_restraints_visualization_complete_specification_v2.md` |

## Pack Contents

| File | Purpose |
|---|---|
| `requirements.md` | Observable notebook requirements, user stories, non-goals, edge cases, and open questions. |
| `notebook-ux-contract.md` | First runnable user flow, editable parameters, trust messaging, and the line between user workflow and fixtures. |
| `notebook-design.md` | Section-by-section notebook architecture, variable handoffs, dependencies, algorithms, visualization design, and degradation behavior. |
| `cell-blueprint.md` | Cell-level execution plan, hidden-state controls, and validation hooks for restart-and-run-all notebook execution. |
| `traceability-matrix.md` | Mapping from PRD goals and requirements to notebook sections, cells, fixtures, and validation evidence. |
| `data-contracts.md` | External source, input, parsed table, computed output, visualization state, provenance, and cache contracts. |
| `fixture-manifest.md` | Provisional fixtures, required expected snapshots, tolerances, and fixture-selection blockers. |
| `tasks.md` | Atomic implementation, validation, documentation, and review tasks for an implementation agent. |
| `validation.md` | Execution, fixture, algorithm, visualization, hidden-state, runtime, and scientific-language validation plan. |
| `docs-plan.md` | Tutorial, how-to, reference, explanation, export, and maintainer documentation plan. |
| `spec-review.md` | Prior review notes retained as advisory history; not part of the required output contract. |

## Summary

This spec pack upgrades the NMR restraints notebook scope to the current lifecycle contract without reopening the product intent in the PRD. The preserved v1 product remains an exploratory NMR evidence viewer for deposited structures and deposited restraints. Mutation workflows and AlphaFold examples remain out of scope.

## Implementation Entry Point

Start with:

1. `requirements.md` for acceptance criteria and open questions.
2. `notebook-ux-contract.md` for the first runnable user flow and trust messaging.
3. `fixture-manifest.md` to use the preserved fixture set and record missing expected snapshots.
4. `notebook-design.md` and `cell-blueprint.md` for architecture and execution order.
5. `tasks.md` for implementation sequencing.

Do not begin final validation until fixture expected snapshots are recorded.

## Blocking Items Before Scientific Sign-Off

- Record implementation-derived expected snapshots for the preserved fixture set, especially `happy-path-9l1v` mapping coverage, density summary, top violations, and a selected local-view residue.
- Confirm `pynmrstar` suitability or document parser replacement.
- Define and complete NMR domain review of formulas, terminology, ambiguity behavior, and interpretation caveats.

## Blocking Questions

- Which implementation-derived snapshots should become the frozen expected values for `happy-path-9l1v` before execution validation begins?
- Does `pynmrstar` preserve the required `9L1V` source tags and counts without wrapper logic on the retained fixture set?
- What domain-review checklist is required before public educational sharing of the notebook?

## Assumptions Accepted For Build

- The retained fixture set is sufficient for implementation scaffolding: `happy-path-9l1v`, `local-copy-9l1v`, `edge-or-9l1v`, `edge-dihedral-9l1v`, `negative-missing-restraints-1crn`, `edge-wrapped-dihedral-synthetic`, and `edge-low-mapping-synthetic`.
- Raw downloaded fixture files live only under the repo, for example `specs/nmr_restraints/fixtures/cache/`, and that cache stays gitignored.
- Local cached copies support parity validation but are not the user workflow themselves.
- Prototype notebook work may start before frozen expected snapshots exist, but restart-and-run-all validation and scientific sign-off may not.

## Required Advisory Reviews

- NMR-aware domain review for restraint semantics, interpretation wording, and threshold choices.
- Lifecycle spec review after fixture snapshots are frozen and before notebook implementation is treated as scientifically ready.

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
- Mutation analysis and AlphaFold-specific examples are excluded from this notebook scope.
