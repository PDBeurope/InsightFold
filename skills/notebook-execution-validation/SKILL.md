---
name: notebook-execution-validation
description: Validate InsightFold notebooks against their spec, fixtures, data contracts, and runtime constraints. Use after notebook implementation to check restart-and-run-all behavior, hidden state, fixture outputs, dependency budget, visualization rendering, documentation presence, and reproducibility readiness.
---

# Notebook Execution Validation

## Overview

Validate that a notebook is runnable and faithful to its spec. This is a mechanical and evidence-producing gate, distinct from qualitative scientific review.

Use this skill on a notebook plus its spec pack. Use `assets/validation-report-template.md` for the output. Read `references/research-basis.md` only when rationale is needed.

## Validation Levels

Use the highest feasible level:

- `static`: inspect notebook JSON/cells/spec without execution
- `smoke`: execute setup and happy-path fixture only
- `full`: restart kernel and run all cells against required fixtures
- `regression`: compare outputs against stored expected snapshots/tolerances

Default to `full` when execution is possible.

## Workflow

### 1. Identify Inputs

Collect:

- notebook path
- source spec pack
- fixture manifest
- validation plan
- target runtime
- dependency policy
- expected outputs/tolerances

If any are missing, record a validation limitation.

### 2. Static Notebook Checks

Inspect:

- first cell explains purpose/input
- setup/import cells are near top
- no unresolved TODOs in critical cells
- no absolute local paths unless explicitly allowed
- no hidden dependence on prior manual state
- optional sections are guarded
- markdown explains major code sections
- plots are labelled by code or surrounding markdown
- notebook metadata does not contain unexpected secrets

### 3. Dependency and Runtime Checks

Check:

- imports match allowed dependencies
- install cells are explicit and minimal
- target runtime is stated
- heavyweight dependencies are justified
- network requirements are documented
- credentials are not required unless stated

When possible, record package versions from execution.

### 4. Execute Notebook

Prefer:

```text
restart kernel -> run all cells -> capture outputs/errors/runtime
```

Acceptable tools include Jupyter execution, papermill, nbclient, or project-local runners. Use the tool already present in the repo/environment when possible.

Do not modify notebooks just to make validation pass unless the user asked for fixes.

### 5. Validate Fixtures

For each required fixture:

- set input/parameters
- run the specified cells or full notebook
- compare expected outputs
- record tolerance mode
- record failures and warnings

If only one fixture can be run, explain why others were skipped.

### 6. Validate Data Contracts

Confirm required schema/field checks exist and are exercised. Negative fixtures should produce clear errors or warnings, not stack traces, unless the spec accepts hard failure.

### 7. Validate Outputs and Documentation

Confirm:

- required tables/figures/summaries appear
- visualization cells produce outputs
- interpretation/limitations are present
- validation summary is present when required
- docs plan artifacts are updated or explicitly deferred

### 8. Produce Evidence

Return a validation report with:

- command/tool used
- fixture(s) run
- runtime
- pass/fail checks
- errors with cell context
- dependency notes
- unresolved risks

## Pass Criteria

A notebook passes execution validation only if:

- required cells execute top-to-bottom
- required fixture checks pass
- required outputs are present
- dependency/runtime constraints are satisfied or explicitly waived
- validation limitations are non-blocking

## Failure Handling

Classify failures:

- `execution-error`: cell failed
- `fixture-mismatch`: output differs from expected
- `contract-gap`: missing schema/field validation
- `hidden-state`: depends on prior state or manual action
- `dependency-drift`: imports/install mismatch spec
- `documentation-gap`: notebook runs but cannot be understood/reused

Do not hide flaky behavior. If rerun passes after a failure, record both.
