---
name: prd-to-notebook-spec
description: Convert an approved InsightFold PRD into an implementation-ready notebook spec pack. Use after `$concept-to-prd` when Codex needs requirements, notebook UX contract, section design, cell blueprints, traceability matrix, data contracts, fixture manifest, validation plan, tasks, documentation plan, and handoff artifacts for spec-driven Jupyter notebook development.
---

# PRD to Notebook Spec

## Overview

Turn a PRD into a spec pack that can drive notebook implementation without relying on undocumented assumptions. This skill owns the implementation-spec material intentionally kept out of `$concept-to-prd`.

Use `assets/notebook-spec-pack-template.md` as the output structure. Use `scripts/scaffold_spec_pack.py` to create a directory skeleton when starting new work. Use `scripts/check_spec_pack.py` before handoff or review.

## Inputs

Use:

- an approved PRD, usually `prd/<feature>.md`
- optional scoping decision capture or advisory notes
- relevant README or lifecycle context
- existing specs, notebooks, fixtures, source modules, or validation reports if updating prior work
- explicit user constraints on runtime, dependencies, audience, visualization, or target environment

## Output Contract

Create a spec pack for each notebook feature:

```text
specs/<feature-slug>/
  spec-pack-overview.md
  requirements.md
  notebook-ux-contract.md
  notebook-design.md
  cell-blueprint.md
  traceability-matrix.md
  data-contracts.md
  fixture-manifest.md
  validation.md
  tasks.md
  docs-plan.md
```

Flat `specs/<feature>_notebook_spec.md` files are acceptable only for small prototypes and must contain equivalent sections.

## Workflow

### 1. Read Product Intent

Extract the product intent without turning it into implementation prematurely:

- target users and audience expertise
- core problem, goals, non-goals, and success criteria
- in-scope and out-of-scope workflows
- expected inputs, outputs, visualizations, and integrations
- scientific formulas, biological assumptions, and interpretation risks
- runtime, dependency, privacy, provenance, and review constraints

### 2. Resolve Ambiguity

List unresolved questions before writing implementation tasks. Mark each as:

- `blocking`: implementation should not start
- `assumption`: implementation can proceed only if recorded
- `non-blocking`: can be resolved during implementation

Do not hide scientific uncertainty. Missing formula sources, data schemas, accession examples, validation references, visualization requirements, or target environments must be visible.

### 3. Write The Notebook UX Contract

Create `notebook-ux-contract.md` before fixture selection. Include:

- primary biological or user question
- target user and assumed expertise
- first runnable user input cell
- user-editable parameters versus internal variables
- happy-path default example
- how fixtures relate to the user flow without becoming the user flow
- expected first visible confirmation that the notebook is working
- interpretation, caveat, and trust messaging requirements

### 4. Write Requirements

Write observable requirements with IDs such as `REQ-001`.

Include:

- user stories or workflow stories
- functional requirements
- non-functional requirements
- acceptance criteria examples
- non-goals
- edge cases and failure behavior
- educational and documentation requirements
- visualization requirements, including MolViewSpec/Mol* for 3D structure viewers

### 5. Write Notebook Design

Translate requirements into notebook architecture:

- section-by-section notebook outline
- section purpose, inputs, outputs, and requirement links
- variable and artifact handoff table
- data-flow diagram in text or Mermaid if useful
- dependency policy and install cell requirements
- parser, API, and algorithm choices with rationale
- visualization plan and fallback behavior
- provenance and reproducibility plan
- optional versus mandatory sections

For 3D molecular visualization, require `$molviewspec-rendering`. Do not specify alternative 3D protein viewer libraries.

### 6. Write Cell Blueprints

Create `cell-blueprint.md` with one row per planned cell or compact cell group:

- cell ID and type: markdown, code, validation, visualization, or export
- user-facing purpose
- inputs consumed
- outputs produced
- hidden-state risks
- restart-run-all expectation
- requirement IDs
- validation hook or expected output

Explicitly cover kernel assumptions, dependency install cells, cell order, cache behavior, and whether the cell is user-editable.

### 7. Define Data Contracts And Fixtures

Every runnable notebook needs at least one happy-path fixture. Prefer:

- one representative happy path
- one edge or negative fixture when feasible
- pinned accessions, versions, files, URLs, or cached snapshots
- expected output snapshots for key values
- network/cache behavior
- provenance shown to the user

For each external API, uploaded file, cached artifact, or generated export, document:

- source endpoint or file type
- required and optional fields
- validation checks
- failure behavior
- provenance and licensing

### 8. Write Traceability Matrix

Create `traceability-matrix.md` mapping:

- PRD goals to requirements
- requirements to notebook sections and cells
- requirements to fixtures
- requirements to validation checks
- requirements to review evidence

No major notebook output should be orphaned from a requirement.

### 9. Write Tasks

Tasks must be executable by an agent or developer:

```markdown
- [ ] T001 Create setup and dependency cell
- [ ] T002 Implement AFDB metadata fetcher producing `metadata`
- [ ] T003 Add validation for missing `paeDocUrl`
```

Mark independent tasks with `[P]` only when they can run in parallel without conflicting files or shared assumptions. Group tasks by phase:

- spec preparation
- data contracts and fixtures
- notebook architecture
- implementation
- validation
- documentation
- review

### 10. Write Validation Plan

Validation must cover:

- restart kernel and run all cells in order
- hidden state and cell-order hazards
- fixture agreement and expected output snapshots
- data contract checks
- scientific formula or algorithm checks
- dependency install/import budget
- output artifact checks
- MolViewSpec/Mol* rendering or documented fallback for 3D views
- visualization interpretability and human review needs
- markdown pedagogy and documentation completeness

If a reference implementation exists, state tolerance and comparison method.

### 11. Write Documentation Plan

Use documentation types deliberately:

- tutorial: how a new user runs the notebook
- how-to: how to adapt it to a new accession, dataset, or mode
- reference: inputs, outputs, parameters, APIs, formulas, schemas
- explanation: scientific background, limitations, interpretation

The notebook should stand alone, but the spec pack should state which external docs, README sections, or integration notes are expected.

## Quality Gate

Reject or revise the spec if:

- implementation tasks are not traceable to requirements
- the notebook UX contract is missing
- fixture data is missing for a runnable notebook
- external API/file contracts are implicit
- scientific assumptions are undocumented
- cell order, hidden state, install cells, or kernel assumptions are unspecified
- acceptance criteria cannot be tested
- visualization outputs lack validation or interpretation checks
- 3D structure visualization does not route through MolViewSpec/Mol*
- documentation is planned only as an afterthought
- dependency/runtime constraints are omitted

## Recommended Skill Chain

Use this skill after:

1. `$idea-scoping-interview`
2. `$scoping-decision-capture`
3. `$concept-to-prd`

Then hand off to:

1. `$notebook-spec-review`
2. `$fixture-selection`
3. `$notebook-from-spec`
4. `$notebook-execution-validation`
5. `$notebook-review`
6. domain-specific scientific implementation skills

