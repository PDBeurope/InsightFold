# Homodimer Diagnostic Documentation Plan

## Documentation Audiences

- First-time non-specialist users.
- AFDB/PDBe internal reviewers.
- Platform engineers maintaining the notebook.
- Future InsightFold agents reusing the pattern.

## Notebook-Inline Documentation

The notebook must include:

- purpose and RUO framing
- explanation of each input
- explanation of PAE, pLDDT, interface residues, contacts, ipTM, ipSAE, pDockQ, pDockQ2, and LIS
- interpretation notes below each visualization
- limitations section
- provenance summary
- validation snapshot

## External Documentation

Expected future docs:

| Document | Purpose | Status |
|---|---|---|
| `notebooks/README.md` | Explain available notebooks and target users | placeholder |
| `specs/homodimer_diagnostic/validation-report.md` | Evidence from notebook execution | generated after validation |
| `specs/homodimer_diagnostic/notebook-review.md` | Qualitative scientific/user-facing review | generated after final review |
| `specs/homodimer_diagnostic/beta-plan.md` | Plan for instrumented beta | placeholder |
| `specs/homodimer_diagnostic/graduation-decision.md` | Integrate, maintain, or archive decision | placeholder |

## Documentation Types

- Tutorial: how to run the notebook with the default accession.
- How-to: how to change accession and rerun.
- Reference: APIs, formulas, fields, metrics, thresholds.
- Explanation: why metrics can disagree and what each metric captures.

## Quality Requirements

- Define jargon before use.
- Avoid clinical claims.
- Link formulas to sources.
- State where values came from.
- Explain visual encodings.
- Make known limitations easy to find.

