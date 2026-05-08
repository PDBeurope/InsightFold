---
name: notebook-spec-review
description: Review InsightFold notebook spec packs before implementation. Use after `$prd-to-notebook-spec` and before notebook building to find missing requirements, vague acceptance criteria, weak fixtures, implicit data contracts, scientific ambiguity, dependency risks, validation gaps, and documentation omissions.
---

# Notebook Spec Review

## Overview

Review a notebook spec pack as an implementation gate. The goal is to stop weak specs before they become confident but wrong notebooks.

Use this skill on either a spec-pack directory:

```text
specs/<feature>/
  requirements.md
  notebook-design.md
  tasks.md
  validation.md
  docs-plan.md
  fixture-manifest.md
  data-contracts.md
```

or a consolidated `*_notebook_spec.md`.

Use `assets/spec-review-report-template.md` for the review output. Read `references/research-basis.md` only when the rationale is needed.

## Review Modes

Use one of three modes:

- `readiness`: decide whether implementation may start
- `scientific`: focus on formulas, biological assumptions, fixtures, and interpretation
- `implementation`: focus on tasks, dependencies, data contracts, notebook architecture, and validation

Default to `readiness` if no mode is specified.

## Review Workflow

### 1. Inventory the Spec Pack

Confirm the spec includes:

- requirements or user stories
- notebook design or section plan
- implementation tasks
- fixture manifest
- data contracts
- validation plan
- documentation plan
- open questions or assumptions

Flag missing artifacts as blocking unless the consolidated spec clearly covers the same content.

### 2. Trace Requirements to Design and Tasks

For every requirement:

- find the notebook section that satisfies it
- find at least one task that implements it
- find validation that proves it

Flag requirements that are unimplemented, untested, or only vaguely represented.

### 3. Check Acceptance Criteria

Acceptance criteria must be observable and testable. Reject criteria such as:

- "works well"
- "is user friendly"
- "handles data"
- "produces useful plots"

Rewrite suggestions should be concrete:

- "Given fixture X, Section 3 displays a PAE heatmap with chain boundaries and a color bar."
- "The notebook run-all completes without errors in the target runtime."
- "The computed score matches reference value Y within tolerance Z."

### 4. Check Scientific and Domain Assumptions

Flag:

- missing formula/source citations
- unreviewed thresholds or cutoffs
- ambiguous biological interpretation
- mismatch between fixture and intended use case
- missing limitations or RUO framing where relevant
- hidden assumptions about organisms, chains, accessions, sequence numbering, variants, structures, or APIs

If a domain expert must decide, mark the issue `human-review`.

### 5. Check Fixtures and Data Contracts

A runnable notebook needs fixtures. Confirm:

- at least one happy-path fixture exists
- edge or negative fixture is included when failure behavior matters
- each fixture has source/provenance
- expected outputs are defined for key values
- network/cache assumptions are explicit

For each external API or file:

- required fields are listed
- optional fields and fallbacks are documented
- schema or shape checks are specified
- error behavior is specified

### 6. Check Notebook Architecture

Confirm:

- notebook sections match user workflow
- variable handoff table is present or easy to infer
- hidden-state risks are addressed
- optional sections cannot break required downstream sections
- install/setup cells are explicit
- dependency budget matches the target runtime
- reusable computations are function/module candidates

### 7. Check Tasks

Tasks should be executable, ordered, and sized appropriately.

Flag tasks that:

- combine unrelated work
- omit outputs
- depend on unclear prior work
- lack validation
- are impossible to parallelize but marked parallel
- write over unrelated files

### 8. Check Documentation Plan

Confirm the plan covers:

- tutorial/run instructions
- how-to adaptation guidance
- reference for inputs, outputs, parameters, APIs, formulas
- explanation of scientific context and limitations

Notebook markdown can satisfy some documentation requirements, but the spec must say where that documentation will live.

## Severity

Use these severity labels:

- `blocking`: implementation should not start
- `major`: implementation can start only with an explicit assumption or owner
- `minor`: fix before final notebook review
- `suggestion`: improves clarity but does not block

## Output

Produce:

- verdict: `pass`, `pass-with-assumptions`, or `fail`
- findings ordered by severity
- required edits
- human-review questions
- updated implementation readiness checklist

Do not rewrite the full spec unless asked. Provide targeted changes.
