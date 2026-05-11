---
name: notebook-review
description: Perform final qualitative review of InsightFold notebooks. Use after implementation and execution validation to assess scientific correctness, educational clarity, reproducibility, visualization quality, maintainability, documentation completeness, lifecycle readiness, and whether human/domain approval is needed.
---

# Notebook Review

## Overview

Review the notebook as a scientific and user-facing artifact. This complements `$notebook-execution-validation`: execution validation proves it runs; notebook review decides whether it is good enough to share, beta test, or graduate.

Use `assets/notebook-review-report-template.md` for the output. Read `references/research-basis.md` only when rationale is needed.

## Review Inputs

Use:

- notebook
- spec pack
- notebook UX contract, if present
- execution validation report
- fixture manifest
- data contracts
- PRD if needed for intent
- source code modules touched by the notebook
- cell blueprint and traceability matrix when present

If execution validation failed, perform review only if the user explicitly asks for partial review.

## Review Dimensions

### 1. Scientific Correctness

Check:

- formulas/algorithms match cited sources
- constants, thresholds, and units are documented
- biological assumptions are explicit
- uncertainty and limitations are honest
- outputs do not overclaim
- interpretation matches evidence

Mark judgment-heavy issues for human/domain review.

### 2. Reproducibility

Check:

- notebook can be rerun with fixtures
- dependencies and runtime are clear
- data provenance is displayed or documented
- random seeds/cache/network assumptions are stated
- outputs are not manually edited
- reusable functions are deterministic where expected

### 3. Pedagogy and Documentation

Check:

- target audience is clear
- jargon is defined before use
- each major section explains purpose before code
- results are interpreted, not just displayed
- tutorial/how-to/reference/explanation needs are covered or planned
- users know what to change and what not to change
- the first runnable input cell is obvious and matches the intended user workflow

### 4. Notebook Structure and Maintainability

Check:

- cells are modular and ordered
- functions are small and named clearly
- hidden state is avoided
- optional sections are isolated
- repeated logic should be factored only when useful
- code likely belongs in a module if reused across notebooks

### 5. Visualization Quality

Check:

- plots answer a user question
- axes, legends, color bars, and labels are present
- color choices are consistent and accessible
- figure size and text are readable
- visual encodings are explained
- 3D protein structure views use MolViewSpec/Mol* and degrade gracefully

### 6. Lifecycle Readiness

Decide recommended next lifecycle state:

- keep iterating in notebook prototype
- move to instrumented beta
- prepare graduation review
- archive

Consider user value, scientific risk, maintenance burden, feedback readiness, and handoff clarity.

Also consider whether the notebook still feels like a user-facing analysis tool rather than a fixture harness.

## Severity

Use:

- `blocking`: do not share or beta
- `major`: fix before beta/graduation
- `minor`: fix before broad publication
- `suggestion`: polish or future improvement
- `human-review`: requires domain or product judgment

## Output

Lead with findings, ordered by severity. Include:

- verdict
- lifecycle recommendation
- scientific risks
- documentation gaps
- required fixes
- human-review questions
- short summary of strengths only after findings

Do not rewrite the notebook unless asked.
