---
name: notebook-from-spec
description: Build an InsightFold Jupyter notebook from a reviewed spec pack. Use after `$prd-to-notebook-spec`, `$notebook-spec-review`, and `$fixture-selection` when Codex needs to implement a notebook while preserving requirements, data contracts, fixtures, validation criteria, documentation intent, and target-runtime constraints.
---

# Notebook From Spec

## Overview

Implement notebooks from the reviewed spec, not from memory of the PRD. The spec pack is the implementation contract.

Use this skill when the user asks to build, scaffold, or update a notebook from:

- `requirements.md`
- `notebook-design.md`
- `tasks.md`
- `fixture-manifest.md`
- `data-contracts.md`
- `validation.md`
- `docs-plan.md`

Use `assets/notebook-outline-template.md` when creating a new notebook outline. Read `references/research-basis.md` only when rationale is needed.

## Implementation Rules

- Follow the reviewed spec exactly unless the user approves a scope change.
- Preserve requirement IDs in comments or markdown where useful.
- Build the smallest notebook that satisfies the spec.
- Keep mandatory sections independent of optional sections.
- Make cells executable top-to-bottom after kernel restart.
- Use explicit variable handoffs; avoid hidden state.
- Keep setup, data loading, computation, visualization, interpretation, and validation separable.
- Include markdown that explains what each section does and how to interpret outputs.
- Treat reusable computations as pure functions or module candidates.

## Workflow

### 1. Confirm Readiness

Before editing notebooks, check:

- spec review verdict is `pass` or `pass-with-assumptions`
- blocking questions are resolved
- fixture manifest exists
- data contracts exist
- validation plan exists

If any are missing, stop and report what must be fixed.

### 2. Build a Notebook Blueprint

Create or confirm:

- notebook filename
- target runtime
- dependency list
- section order
- cell purposes
- variables produced and consumed
- outputs expected from fixtures
- optional sections and fallbacks

Do not start implementation until the blueprint matches the spec.

### 3. Implement Section by Section

For each section:

1. Add explanatory markdown.
2. Add code cell with one logical responsibility.
3. Produce named outputs expected by downstream sections.
4. Add lightweight checks or assertions where the spec requires them.
5. Keep visualization cells self-contained and labelled.

Prefer this section pattern:

```text
Markdown: goal and interpretation
Code: inputs and computation
Code/Markdown: output display and explanation
```

### 4. Respect Data Contracts

For every API/file input:

- validate required fields before use
- handle optional fields explicitly
- show clear user-facing errors
- record provenance in the notebook output when useful

Never silently assume a field exists if the spec says it is optional.

### 5. Integrate Fixtures

Use the happy-path fixture as the default notebook input unless the spec says otherwise.

Add fixture-aware checks:

- expected shapes
- expected metadata fields
- expected numeric values or tolerances
- expected warnings/failures for negative fixtures

Do not hard-code fixture outputs into computation logic.

### 6. Add Validation Section

The notebook should include an internal validation or sanity-check section when the spec requires it. This section should summarize:

- fixture used
- key checks passed/failed
- runtime or dependency notes
- known limitations

### 7. Update Tasks

As implementation progresses, update `tasks.md` only for tasks actually completed. Do not mark review/validation tasks complete until they are performed.

## Notebook Quality Requirements

The produced notebook must:

- run top-to-bottom after kernel restart
- have a clear first cell explaining purpose and inputs
- have explicit setup/import cells
- avoid hidden local file assumptions
- show provenance for fetched or uploaded data
- label plots and include legends/color bars where appropriate
- include interpretation and limitations
- avoid unexplained magic constants
- avoid broad dependencies not listed in the spec

## Output

Return:

- notebook path
- changed files
- completed tasks
- assumptions used
- validation still required

Do not claim the notebook is validated unless `$notebook-execution-validation` has been run.
