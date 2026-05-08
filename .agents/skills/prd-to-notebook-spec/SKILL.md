---
name: prd-to-notebook-spec
description: Convert an approved InsightFold PRD into spec-driven notebook development documentation. Use after `$convert-to-prd` or whenever Codex needs to create implementation-ready requirements, notebook design, tasks, fixtures, validation criteria, and documentation plans for a notebook-driven development lifecycle stage.
---

# PRD to Notebook Spec

## Overview

Turn a PRD into a spec pack that can drive notebook implementation without relying on undocumented assumptions. Use a hybrid of spec-driven development and notebook-driven development: requirements first, then technical/notebook design, then executable tasks, with InsightFold-specific additions for fixtures, data contracts, reproducibility, and documentation.

Load `references/research-basis.md` only when you need the source rationale. Use `assets/notebook-spec-pack-template.md` as the output structure.

## Inputs

Use:

- an approved PRD, usually `prd/<feature>.md`
- relevant README or lifecycle context
- existing specs, notebooks, fixtures, source modules, or validation reports if updating prior work
- explicit user constraints on runtime, dependencies, audience, or target environment

## Output Contract

Create a spec pack for each notebook feature. Prefer a directory for new work:

```text
specs/<feature-slug>/
  requirements.md
  notebook-design.md
  tasks.md
  validation.md
  docs-plan.md
  fixture-manifest.md
  data-contracts.md
```

If the repo or user expects the older flat pattern, create a consolidated file instead:

```text
specs/<feature-slug>_notebook_spec.md
```

The spec pack must be clear enough that a separate implementation agent can build the notebook from it without rereading the PRD.

## Workflow

### 1. Read the PRD

Extract:

- product/use-case title
- target users and audience expertise
- problem statement and success criteria
- in-scope and out-of-scope items
- required inputs, outputs, and external APIs
- scientific formulas, algorithms, or biological assumptions
- non-functional constraints such as runtime, dependencies, privacy, and execution target
- risks, unknowns, and open questions

Keep the PRD as product intent. Do not invent implementation details until the design section.

### 2. Resolve Ambiguity

Before writing implementation tasks, list unresolved questions. Mark each as:

- `blocking`: implementation should not start
- `assumption`: proceed only if explicitly recorded
- `non-blocking`: can be resolved during implementation

Do not hide scientific uncertainty. If formula sources, data schemas, accession examples, validation references, or target environments are missing, record them.

### 3. Write Requirements

Write requirements as observable behavior. Include:

- user stories or workflow stories
- acceptance criteria
- non-goals
- edge cases and failure behavior
- audience and educational requirements
- runtime and dependency requirements

Use requirement IDs such as `REQ-001`. Acceptance criteria should be testable.

### 4. Write Notebook Design

Translate requirements into a notebook architecture:

- section-by-section notebook outline
- cell intent for each section
- variable and artifact handoff table
- data-flow diagram in text or Mermaid if useful
- dependency policy and install cell requirements
- parser/algorithm choices with rationale
- visualization plan
- error handling and graceful degradation
- optional versus mandatory sections

For notebook-driven development, every major output should trace back to a requirement and every reusable computation should be isolated as a function or module candidate.

### 5. Define Fixtures and Data Contracts

Every notebook spec needs at least one runnable fixture. Prefer:

- one happy-path fixture
- one edge or negative fixture when feasible
- expected output snapshots for key values
- clear network/cache behavior

For each external data source or uploaded file, document:

- source endpoint or file type
- required fields
- optional fields
- validation checks
- failure behavior
- provenance shown to the user

### 6. Write Tasks

Tasks must be executable by an agent or developer. Use atomic tasks with clear outputs:

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

### 7. Write Validation Plan

Validation must cover:

- top-to-bottom notebook execution
- fixture agreement
- scientific formula or algorithm checks
- data contract checks
- hidden-state and execution-order hazards
- dependency/runtime budget
- visualization existence and interpretability
- documentation completeness

If a reference implementation exists, state tolerance and comparison method.

### 8. Write Documentation Plan

Use documentation types deliberately:

- tutorial: how a new user runs the notebook
- how-to: how to adapt it to a new accession, dataset, or mode
- reference: inputs, outputs, parameters, APIs, formulas, schemas
- explanation: scientific background, limitations, interpretation

The notebook itself should include enough explanation to stand alone, but the spec pack should also state which external docs or README updates are expected.

## Quality Gate

Reject or revise the spec if:

- implementation tasks are not traceable to requirements
- fixture data is missing for a notebook that claims to be runnable
- external API/file contracts are implicit
- scientific assumptions are undocumented
- the notebook design depends on hidden state
- acceptance criteria cannot be tested
- documentation is only planned as an afterthought
- dependency/runtime constraints are omitted

## Recommended First Skill Chain

Use this skill after:

1. `$convert-to-prd`

Then hand off to later skills such as:

1. `$fixture-selection`
2. `data-contract-validation`
3. `$notebook-execution-validation`
4. domain-specific scientific implementation skills
