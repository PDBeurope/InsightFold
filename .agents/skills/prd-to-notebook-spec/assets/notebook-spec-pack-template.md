# Notebook Spec Pack: <Notebook / Feature Name>

## 0. Metadata

| Field | Value |
|---|---|
| Spec ID | |
| Source PRD | |
| Lifecycle Stage | Scope |
| Status | Draft / Review / Approved |
| Owner | |
| Last Updated | |
| Target Runtime | Colab / local Jupyter / hosted beta / other |

## 1. Requirements

### 1.1 User Story

As a `<user type>`, I want to `<action>`, so that `<benefit>`.

### 1.2 Functional Requirements

| ID | Requirement | Acceptance Criteria |
|---|---|---|
| REQ-001 | | |

### 1.3 Non-Functional Requirements

| ID | Requirement | Acceptance Criteria |
|---|---|---|
| NFR-001 | Runtime and dependency budget | |

### 1.4 Non-Goals

- 

### 1.5 Open Questions

| ID | Question | Type | Owner | Resolution |
|---|---|---|---|---|
| Q-001 | | blocking / assumption / non-blocking | | |

## 2. Notebook Design

### 2.1 Notebook Sections

| Section | Purpose | Inputs | Outputs | Requirement IDs |
|---|---|---|---|---|
| 1. Setup | | | | |

### 2.2 Data Flow

```text
input -> fetch/parse -> compute -> visualize -> interpret -> validate
```

### 2.3 Variable and Artifact Handoff

| Name | Produced By | Consumed By | Type / Shape | Notes |
|---|---|---|---|---|
| | | | | |

### 2.4 Dependency Policy

| Dependency | Required? | Purpose | Install Source | Constraint |
|---|---|---|---|---|
| | | | | |

### 2.5 Design Decisions

| Decision | Rationale | Alternatives Considered | Risk |
|---|---|---|---|
| | | | |

## 3. Data Contracts

| Source | Required Fields | Optional Fields | Validation | Failure Behavior |
|---|---|---|---|---|
| | | | | |

## 4. Fixture Manifest

| Fixture ID | Source | Why Chosen | Expected Outputs | Notes |
|---|---|---|---|---|
| happy-path | | | | |
| edge-case | | | | |

## 5. Tasks

### Phase A: Spec Preparation

- [ ] T001 Resolve blocking open questions

### Phase B: Data and Fixtures

- [ ] T010 Create fixture manifest
- [ ] T011 Document data contracts

### Phase C: Notebook Architecture

- [ ] T020 Create section skeleton
- [ ] T021 Define variable handoff table in notebook comments or spec docs

### Phase D: Implementation

- [ ] T030 Implement first executable section

### Phase E: Validation

- [ ] T040 Run notebook top-to-bottom on pinned fixture

### Phase F: Documentation

- [ ] T050 Add user-facing run instructions
- [ ] T051 Add limitation and interpretation notes

## 6. Validation Plan

| Check | Method | Pass Criteria |
|---|---|---|
| Top-to-bottom execution | Restart kernel and run all | No errors |
| Fixture agreement | Compare expected outputs | |
| Dependency audit | Inspect imports/install cells | |
| Documentation completeness | Review against docs plan | |

## 7. Documentation Plan

| Documentation Type | Artifact | Audience | Notes |
|---|---|---|---|
| Tutorial | Notebook intro / README | New users | |
| How-to | Adaptation section | Returning users | |
| Reference | Inputs, outputs, parameters, APIs | Developers / power users | |
| Explanation | Scientific background and limitations | All users | |

## 8. Handoff Notes

- Implementation agent entry point:
- Required domain skills:
- Review agent checklist:
- Graduation evidence to collect:
