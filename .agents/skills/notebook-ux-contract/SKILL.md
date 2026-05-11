---
name: notebook-ux-contract
description: Define the first runnable user workflow for an InsightFold notebook before fixture selection and implementation. Use when a notebook needs a clear user-input entrypoint, a section storyboard, and a separation between user-editable parameters and validation fixtures.
---

# Notebook UX Contract

## Overview

Capture the user-facing contract of a notebook before fixtures and implementation details take over. The output should make it obvious what the user enters first, what the notebook does by default, and how fixtures/examples relate to the main workflow.

Use this skill after a PRD or notebook intent document and before fixture selection, spec review, or implementation.

## Output

Produce a concise UX contract containing:

- primary biological question or diagnostic goal
- target user and assumed expertise
- target explanation level
- knowledge gaps to support
- first runnable input cell
- user-editable parameters and safe defaults
- notebook section storyboard
- fixture/example relationship
- what users should edit versus not edit
- UX acceptance criteria

## Workflow

### 1. Define the User Entry Point

Specify the first runnable input cell as a general pattern, not a notebook-specific hardcode.

Required qualities:

- the user can run the notebook without editing internal dicts or function definitions
- defaults are valid and safe
- parameter names match the notebook domain
- optional inputs are clearly marked

### 2. Define the Main Question

State the one biological or diagnostic question the notebook is answering. If the notebook has more than one major question, identify the primary one and mark others as secondary.

### 3. Define the Explanation Level

Calibrate the notebook prose for the intended user. Pick a target explanation level such as:

- `concise expert`
- `guided practitioner`
- `beginner-aware lightweight`

Record:

- what the user is already expected to know
- what the notebook should explain briefly
- what should be deferred to references or external documentation

The goal is educational enough to support use and interpretation without turning the notebook into a full tutorial.

### 4. Define the Storyboard

List the notebook sections in user order:

- input
- retrieval or loading
- parsing and normalization
- analysis
- comparison or interpretation
- visualization
- validation or export

Keep the storyboard aligned to the user workflow, not the internal code layout.

### 5. Define the Markdown Rule

For each major section, define a short markdown contract:

- what this step does
- why it matters
- how to read the output
- caveats only when relevant

Prefer brief, readable prose over long tutorial blocks unless the target explanation level explicitly calls for more detail.

### 6. Separate Fixtures From UX

Describe fixtures as:

- validation inputs
- examples for learning
- regression checks

Do not let fixtures become the primary user interface unless the notebook is explicitly a test harness.

### 7. Write UX Acceptance Criteria

Include concrete checks such as:

- a non-developer can start the analysis from the first input cell
- the explanation level matches the target user without overwhelming them
- the notebook explains what each editable field means
- the notebook shows what changed when the user edits the input
- validation fixtures do not hide the real user entrypoint

## Suggested Template

```markdown
## Notebook UX Contract

Primary question:
- ...

Target user:
- ...

Explanation level:
- ...

Knowledge gaps to support:
- ...

First runnable input cell:
- ...

User-editable parameters:
- ...

Storyboard:
1. ...
2. ...

Per-section markdown rule:
- ...

Fixture relationship:
- ...

UX acceptance criteria:
- ...
```

## Quality Gate

Reject the contract if:

- the first runnable input cell is still a fixture selector
- the user must edit internal variables to run the main workflow
- fixtures are not clearly separated from the user flow
- the notebook question is vague enough to support multiple incompatible entrypoints
- the explanation level is undefined or obviously mismatched to the target user
