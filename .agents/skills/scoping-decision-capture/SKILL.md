---
name: scoping-decision-capture
description: Convert an InsightFold scoping conversation, advisory panel discussion, meeting notes, or ideation transcript into a structured decision capture. Use before `$concept-to-prd` when Codex needs to preserve the core problem, system behavior, constraints, scope decisions, user journey, open questions, and graduation signals without prematurely writing a PRD.
---

# Scoping Decision Capture

## Purpose

Summarize a discovery conversation into explicit decisions, rationale, open questions, and lifecycle signals. This skill is the bridge between free-form ideation and `$concept-to-prd`.

Use `assets/decision-capture-template.md` as the output shape.

## Inputs

Use:

- conversation transcript or notes
- advisory panel feedback
- rough feature or notebook description
- existing concept brief
- user corrections or decisions from a planning session

## Outputs

Create a structured markdown decision capture containing:

- core problem or question
- how the proposed system works
- constraints, caveats, and reliability considerations
- explicit scope decisions and tradeoffs
- user journey or narrative flow
- open questions and deferred decisions
- success criteria and graduation or stop signals

Prefer `prompts/templates/<feature>-decision-capture.md` or a user-specified path when saving is requested. Otherwise return the markdown directly.

## Workflow

1. Separate resolved decisions from open questions.
2. Preserve user language where it captures intent, but rewrite vague statements into implementation-relevant prose.
3. Record why each major decision was made and what tradeoff was accepted.
4. Capture reliability, trust, caveats, and failure modes as first-class content.
5. Describe the intended user journey step by step.
6. Identify which advisory agents should review the capture before PRD writing.
7. End with a readiness judgment for `$concept-to-prd`.

## Quality Gate

Reject or revise the capture if it:

- presents unresolved items as decisions
- omits user journey or interpretation flow
- lacks explicit scope boundaries
- fails to name validation or fixture needs
- ignores scientific caveats or trust limits
- is too generic for another agent to write a PRD from it

