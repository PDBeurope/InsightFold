# Notebook Spec Pack: <Notebook / Feature Name>

## Files

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

## spec-pack-overview.md

```markdown
# <Feature> Spec Pack Overview

| Field | Value |
|---|---|
| Spec ID | |
| Source PRD | |
| Lifecycle Stage | Scope |
| Status | Draft / Review / Approved |
| Owner | |
| Last Updated | |
| Target Runtime | Colab / local Jupyter / hosted beta / other |
| Target Notebook | notebooks/<feature>.ipynb |

## Summary

## Blocking Questions

## Assumptions Accepted For Build

## Required Advisory Reviews
```

## requirements.md

```markdown
# Requirements

## User Story

As a `<user type>`, I want to `<action>`, so that `<benefit>`.

## Functional Requirements

| ID | Requirement | Acceptance Criteria | Priority |
|---|---|---|---|
| REQ-001 | | | Must |

## Non-Functional Requirements

| ID | Requirement | Acceptance Criteria |
|---|---|---|
| NFR-001 | Runtime and dependency budget | |

## Non-Goals

## Open Questions

| ID | Question | Type | Owner | Resolution |
|---|---|---|---|---|
| Q-001 | | blocking / assumption / non-blocking | | |
```

## notebook-ux-contract.md

```markdown
# Notebook UX Contract

## Primary User Question

## Target User And Expertise

## First Runnable User Input Cell

| Parameter | Default | User Editable? | Validation | Notes |
|---|---|---|---|---|
| | | yes / no | | |

## Happy-Path Default Example

## User Flow

| Stage | User Action | System Action | User Sees | Interpretation Need |
|---|---|---|---|---|
| 1 | | | | |

## Fixtures Versus User Flow

## Trust, Caveats, And Limitations To Show In Notebook
```

## notebook-design.md

```markdown
# Notebook Design

## Section Outline

| Section | Purpose | Inputs | Outputs | Requirement IDs |
|---|---|---|---|---|
| 1. Setup | | | | |

## Data Flow

    input -> fetch/parse -> compute -> visualize -> interpret -> validate

## Variable And Artifact Handoff

| Name | Produced By | Consumed By | Type / Shape | Notes |
|---|---|---|---|---|
| | | | | |

## Dependency Policy

| Dependency | Required? | Purpose | Install Source | Constraint |
|---|---|---|---|---|
| | | | | |

## Visualization Plan

For 3D molecular views, use MolViewSpec/Mol* through `$molviewspec-rendering`.

## Design Decisions

| Decision | Rationale | Alternatives Considered | Risk |
|---|---|---|---|
| | | | |
```

## cell-blueprint.md

```markdown
# Jupyter Cell Blueprint

| Cell ID | Type | Purpose | Inputs | Outputs | User Editable? | Hidden-State Risk | Validation Hook | Requirements |
|---|---|---|---|---|---|---|---|---|
| C001 | code | | | | yes / no | | | REQ-001 |
```

## traceability-matrix.md

```markdown
# Requirement Traceability Matrix

| Requirement | PRD Goal | Notebook Section / Cell | Fixture | Validation Check | Review Evidence |
|---|---|---|---|---|---|
| REQ-001 | G-001 | | | | |
```

## data-contracts.md

```markdown
# Data Contracts

| Source / Artifact | Required Fields | Optional Fields | Validation | Failure Behavior | Provenance |
|---|---|---|---|---|---|
| | | | | | |
```

## fixture-manifest.md

```markdown
# Fixture Manifest

| Fixture ID | Source | Why Chosen | Expected Outputs | Cache / Version | Notes |
|---|---|---|---|---|---|
| happy-path | | | | | |
| edge-case | | | | | |
```

## validation.md

```markdown
# Validation Plan

| Check | Method | Pass Criteria | Fixture / Artifact | Requirement |
|---|---|---|---|---|
| Restart-and-run-all | Restart kernel and run all | No errors | Notebook | NFR-001 |
| Fixture agreement | Compare expected outputs | | | |
| Dependency audit | Inspect imports/install cells | | | |
| Visualization rendering | Verify expected plots/views or fallback | | | |
| Documentation completeness | Review against docs plan | | | |
```

## tasks.md

```markdown
# Tasks

## Phase A: Spec Preparation

- [ ] T001 Resolve blocking open questions

## Phase B: Data And Fixtures

- [ ] T010 Create fixture manifest
- [ ] T011 Document data contracts

## Phase C: Notebook Architecture

- [ ] T020 Create section skeleton
- [ ] T021 Define variable handoff table

## Phase D: Implementation

- [ ] T030 Implement first executable section

## Phase E: Validation

- [ ] T040 Run notebook top-to-bottom on pinned fixture

## Phase F: Documentation

- [ ] T050 Add user-facing run instructions
- [ ] T051 Add limitation and interpretation notes
```

## docs-plan.md

```markdown
# Documentation Plan

| Documentation Type | Artifact | Audience | Notes |
|---|---|---|---|
| Tutorial | Notebook intro / README | New users | |
| How-to | Adaptation section | Returning users | |
| Reference | Inputs, outputs, parameters, APIs | Developers / power users | |
| Explanation | Scientific background and limitations | All users | |
```
