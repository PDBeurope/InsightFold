---
name: fixture-selection
description: Select and document pinned fixtures for InsightFold notebooks. Use after notebook spec review or whenever a notebook needs representative accessions, files, APIs, expected outputs, edge cases, provenance, and validation snapshots before implementation or regression testing.
---

# Fixture Selection

## Overview

Choose fixtures that make a notebook runnable, reviewable, and regression-testable. A fixture is not just an example input; it is a pinned input plus provenance, why it matters, expected outputs, and validation notes.

Use `assets/fixture-manifest-template.md` for the output. Read `references/research-basis.md` only when source rationale is needed.

## Fixture Set Policy

Prefer a minimal but meaningful set:

- `happy-path`: demonstrates the intended core workflow
- `edge-case`: stresses an important boundary, ambiguity, or known failure mode
- `negative`: confirms graceful failure for invalid/missing/unsupported data
- `reference`: matches a published or trusted reference implementation when available

For early prototypes, one happy-path fixture is acceptable only if the spec explicitly states why edge fixtures are deferred.

## Selection Workflow

### 1. Read the Spec

Extract:

- target biological object or dataset type
- required inputs and APIs
- success criteria
- expected computations and visualizations
- target runtime and dependency constraints
- risks and open questions

The fixture set must test the spec, not merely demonstrate a convenient example.

### 2. Define Fixture Roles

For each fixture role, state the property being tested:

- typical successful input
- small input for fast runtime
- large input near performance boundary
- missing optional field
- malformed or unsupported accession/file
- conflicting scientific signals
- known benchmark value

Avoid choosing examples only because they are familiar.

### 3. Record Provenance

For every fixture, record:

- source database or file path
- accession or stable identifier
- retrieval endpoint or command
- retrieval date if network data is used
- cached artifact path if present
- license or use constraint if relevant
- expected schema version or data format

Do not rely on an unpinned web query if the notebook must be reproducible.

### 4. Define Expected Outputs

Expected outputs can be numeric, structural, visual, or behavioral.

Examples:

- metadata fields present
- matrix shape
- chain count
- residue count
- score values and tolerances
- number of plotted panels
- expected warning for missing optional field
- successful top-to-bottom execution under runtime budget

Mark values as:

- `exact`
- `tolerance`
- `range`
- `presence`
- `manual-review`

### 5. Check Runtime and Practicality

A good fixture should:

- run quickly enough for routine validation
- be small enough for repo or cache policy
- exercise real code paths
- avoid unnecessary external dependencies
- not require credentials unless the whole notebook requires credentials

If a scientifically ideal fixture is too expensive, add a smaller smoke-test fixture and document the larger fixture as an optional benchmark.

### 6. Update Spec and Validation

Add fixture IDs to:

- `fixture-manifest.md`
- `validation.md`
- `tasks.md`
- notebook parameter/default input section

Every fixture should have at least one validation check.

## Output Rules

Return a fixture manifest and a short rationale. If fixtures cannot be finalized because source data must be fetched or a domain expert must choose examples, return candidate fixtures and mark blockers.

## Quality Gate

Reject a fixture set if:

- no happy-path fixture exists
- fixture has no stable identifier or source path
- expected outputs are absent
- selected fixture does not exercise the notebook's main claim
- external data contract is implicit
- runtime is incompatible with the notebook target
