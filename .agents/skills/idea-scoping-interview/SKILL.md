---
name: idea-scoping-interview
description: Guide early InsightFold concept discovery before PRD writing. Use when Codex needs to interview the user about a rough notebook, biological use case, scientific question, product opportunity, target user, constraints, risks, visualizations, validation needs, or scope boundaries before creating a PRD or spec.
---

# Idea Scoping Interview

## Purpose

Run a structured discovery conversation for an early InsightFold notebook idea. The goal is to surface the scientific, user, product, technical, validation, and lifecycle dimensions before committing to a PRD.

Use `assets/opening-scoping-prompt-template.md` when the user wants a reusable prompt. Use `references/scoping-question-bank.md` when the idea needs a deeper panel-style interview.

## Inputs

Use any combination of:

- a rough notebook or feature idea
- the intended biological object, system, dataset, or workflow
- target users and their expected expertise
- known constraints, dependencies, APIs, files, or runtime environments
- candidate expert personas or advisory panel roles
- existing PRDs, specs, notebooks, or decision captures

## Outputs

Produce one of:

- a guided interview with sequenced questions
- a scoped concept brief ready for `$scoping-decision-capture`
- a list of blocking uncertainties before PRD writing
- a recommendation for which advisory agents should join the conversation

Do not produce a full PRD unless the user explicitly asks to continue with `$concept-to-prd`.

## Workflow

1. Identify the single most important question the notebook or feature is trying to answer.
2. Confirm the target user, their expertise, and what they should learn or decide.
3. Probe the biological and scientific basis: data sources, assays, structures, confidence signals, interpretation risks, and failure modes.
4. Probe the product workflow: inputs, outputs, visualizations, controls, ranking/filtering logic, exports, and what is out of scope.
5. Probe implementation constraints: runtime, dependencies, APIs, caching, fixtures, privacy, reproducibility, and validation.
6. Probe lifecycle fit: whether this is a standing notebook, a candidate for AFDB/PDBe integration, or an experiment likely to be archived.
7. Keep conclusions tentative until decisions are explicitly made by the user.

## Advisory Panel Routing

When the idea is biology-heavy, recommend advisory agents before or during discovery:

- `scientific-product-manager` for user problem, scope, and PRD readiness
- `computational-structural-biologist` for structural interpretation and overclaim risk
- `bioinformatics-data-engineer` for APIs, identifiers, files, provenance, and reproducibility
- `wet-lab-liaison`, `protein-biochemist`, `enzymologist`, or `protein-purification-scientist` when the notebook claims experimental relevance
- `molecular-visualization-specialist` for MolViewSpec/Mol* 3D viewer requirements
- `evaluation-benchmarking-specialist` for metrics, fixtures, baselines, leakage, and validation

Avoid invoking business/startup personas unless the user is making investment, hiring, commercialization, or company-strategy decisions.

## Quality Gate

Stop and ask targeted questions if the concept lacks:

- a specific user and use context
- a clear biological or scientific question
- a reason the work should be a notebook first
- expected input and output types
- a validation or fixture strategy
- explicit scope boundaries
- known assumptions and risks

