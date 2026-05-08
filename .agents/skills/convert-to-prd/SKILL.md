---
name: convert-to-prd
description: Convert rough InsightFold notebook ideas, summary notes, decisions, or planning captures into a structured life-sciences PRD. Use when Codex needs to scope an idea before notebook specification, define users, goals, non-goals, risks, success criteria, dependencies, data needs, and lifecycle readiness.
---

# Convert To PRD

## Purpose

Turn unstructured product or scientific notebook ideas into a PRD that can drive the rest of the InsightFold lifecycle.

Use `assets/PRD-TEMPLATE.md` as the output structure when a full PRD is needed. Keep the PRD specific enough that `$prd-to-notebook-spec` can create implementation-ready notebook requirements without guessing.

## Inputs

Use:

- a raw notebook idea, summary, decisions capture, or planning note
- any known biological object, target user, data source, or expected notebook output
- existing PRDs or specs if the task is an update
- README or lifecycle context when the idea needs InsightFold framing

## Outputs

Create or update:

- a PRD, usually `prd/<feature-slug>.md`
- a short list of blocking questions if the idea is not ready for PRD conversion
- explicit assumptions when proceeding with incomplete source material

## Instructions

1. Read the source idea, summary, or decisions capture.
2. Extract the user problem, target audience, biological object, expected notebook output, and why the work belongs in InsightFold.
3. Separate goals from non-goals before adding implementation detail.
4. Record assumptions, risks, unknowns, dependencies, and human/domain review needs.
5. Use `assets/PRD-TEMPLATE.md` for the final shape when producing a full PRD.
6. Leave fields blank or explicitly mark `TBD` when the source material does not support an inference.

## Output Requirements
- Output must be valid markdown.
- Prefer `prd/<feature-slug>.md` for new PRDs.
- Keep tables and section headers from the template unless the user requests a lighter PRD.
- Maintain life-sciences, structural-biology, AFDB, PDBe, or bioinformatics context where relevant.
- Do not overstate clinical, production, or public-facing readiness.

## Mapping Guidelines

- `Document Overview`: product name, status, stakeholders, related artifacts, and regulatory or RUO notes.
- `Executive Summary`: problem, proposed notebook, target users, and success criteria.
- `Background & Context`: scientific context, existing workflows, and why now.
- `Goals & Non-goals`: scope boundaries before implementation begins.
- `Use Cases` and `User Stories`: observable user behavior, not vague value statements.
- `Functional Requirements`: notebook inputs, outputs, calculations, visualizations, and user-facing behaviors.
- `Data Requirements`: APIs, files, fixtures, provenance, schemas, and caching assumptions.
- `Risks & Mitigations`: scientific, reproducibility, dependency, and adoption risks.
- `Open Questions`: unresolved items that `$prd-to-notebook-spec` must not silently invent.

## Quality Gate

Reject or revise the PRD if it lacks:

- a specific target user
- a clear notebook outcome
- explicit non-goals
- success criteria that can be tested later
- data-source or fixture expectations
- risk and assumption tracking
