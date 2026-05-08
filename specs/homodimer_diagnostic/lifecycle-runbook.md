# Homodimer Diagnostic Lifecycle Runbook

This runbook shows how to test the InsightFold lifecycle skills and agents using:

```text
prd/homodimer_diagnostic_notebook_prd.md
specs/homodimer_diagnostic/
notebooks/homodimer_diagnostic.ipynb
```

## Current Status

Completed:

- PRD exists.
- Directory-style spec pack exists at `specs/homodimer_diagnostic/`.
- Spec review exists at `specs/homodimer_diagnostic/spec-review.md`.

Current gate:

- `pass-with-assumptions`.

Main blockers before full validation or beta:

- curate borderline fixture accession
- curate disagreement fixture accession
- curate valid non-homodimer/unsupported fixture
- freeze expected numeric outputs
- pin IPSAE reference version/commit
- define score interpretation bands

## 1. Fixture Curation

Use this when you are ready to resolve the fixture blockers.

```text
Act as agents/lifecycle/fixture-curator.md.
Use skills/fixture-selection/SKILL.md.

Input spec pack:
specs/homodimer_diagnostic/

Review:
specs/homodimer_diagnostic/spec-review.md

Task:
Finish fixture selection for the homodimer diagnostic notebook.
Keep FX-001 as the happy-path candidate unless there is a better accession.
Choose or request candidate accessions for:
- borderline confidence
- metric disagreement
- valid AFDB input that is unsupported by v1
- malformed/negative input

Update:
specs/homodimer_diagnostic/fixture-manifest.md
specs/homodimer_diagnostic/validation.md
specs/homodimer_diagnostic/tasks.md
```

Expected result:

- fixture manifest moves from partial to ready for smoke or full validation
- unresolved fixture choices are clearly marked for human/domain review

## 2. Notebook Build

Use this when the fixture manifest is good enough for at least smoke implementation.

```text
Act as agents/lifecycle/notebook-builder.md.
Use skills/notebook-from-spec/SKILL.md.

Input spec pack:
specs/homodimer_diagnostic/

Spec review:
specs/homodimer_diagnostic/spec-review.md

Target notebook:
notebooks/homodimer_diagnostic.ipynb

Task:
Implement or update the notebook from the reviewed spec.
Preserve the existing notebook if useful, but make the final notebook satisfy:
- requirements.md
- notebook-design.md
- data-contracts.md
- fixture-manifest.md
- validation.md

Do not mark validation or review tasks complete.
```

Expected result:

- notebook updated or rebuilt
- implementation tasks in `tasks.md` marked complete only where work was actually done
- validation remains pending

## 3. Execution Validation

Use this after notebook implementation.

```text
Act as agents/lifecycle/notebook-validator.md.
Use skills/notebook-execution-validation/SKILL.md.

Notebook:
notebooks/homodimer_diagnostic.ipynb

Spec pack:
specs/homodimer_diagnostic/

Task:
Run the highest feasible validation level.
Start with static validation if notebook execution is blocked.
Prefer restart-and-run-all validation with FX-001.
Record all skipped fixture checks as limitations.

Write the report to:
specs/homodimer_diagnostic/validation-report.md
```

Expected result:

- validation report with pass/fail/limitations
- runtime and dependency evidence
- failures classified as execution, fixture, contract, hidden-state, dependency, or documentation gaps

## 4. Final Notebook Review

Use this only after execution validation exists.

```text
Act as agents/lifecycle/notebook-reviewer.md.
Use skills/notebook-review/SKILL.md.

Notebook:
notebooks/homodimer_diagnostic.ipynb

Spec pack:
specs/homodimer_diagnostic/

Validation report:
specs/homodimer_diagnostic/validation-report.md

Task:
Review the notebook as a scientific and user-facing artifact.
Assess scientific correctness, reproducibility, pedagogy, visualization quality, maintainability, and lifecycle readiness.

Write the report to:
specs/homodimer_diagnostic/notebook-review.md
```

Expected result:

- recommendation: continue prototype iteration, move to beta, prepare graduation review, or archive
- human/domain review questions
- required fixes before beta or graduation

## 5. Placeholder Future Stages

The following lifecycle stages are not implemented yet:

```text
beta-feedback-instrumentation
usage-signal-summary
beta-release-packaging
graduation-decision-brief
engineering-handoff-package
standing-notebook-maintenance-plan
archive-retrospective
```

Until these exist, use manual documents:

```text
specs/homodimer_diagnostic/beta-plan.md
specs/homodimer_diagnostic/feedback-summary.md
specs/homodimer_diagnostic/graduation-decision.md
specs/homodimer_diagnostic/maintenance-plan.md
specs/homodimer_diagnostic/archive-retrospective.md
```

