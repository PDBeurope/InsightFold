Based on the summary below, produce a complete markdown specification for 
the EnzyMM AlphaFold notebook. This spec will be used for spec-driven 
development — it should be detailed enough that an engineer can implement 
the notebook from it without ambiguity.

[PASTE SUMMARY HERE]

Structure the spec as follows:

---

# Notebook Spec: EnzyMM on AlphaFold Structures

## Overview
One paragraph. The biological question, the tool, the audience, and what 
the notebook delivers.

## Scientific Context
The biology a developer needs to understand to implement this correctly. 
Include: what catalytic motifs are and why they matter, how EnzyMM finds 
them, what it means to find one in a predicted structure, and the key 
caveats the notebook must communicate to the user.

## Input Specification
- Accepted input format(s)
- Validation logic (what to check, what to reject, what to warn)
- Data fetching steps (AFDB REST API endpoint, expected response format, 
  error handling)

## EnzyMM Execution
- Installation / import approach in Colab
- CLI invocation or Python API call with exact parameters
- Parameter rationale (especially pLDDT cutoff if used)
- Expected output files

## Output Parsing and Interpretation
- Which fields from the TSV to use and how to interpret each
- Hit filtering and prioritisation logic
- How to communicate match quality to the user (RMSD thresholds, 
  predicted_correct flag, completeness)
- How to handle the zero-hit case

## Visualisations
For each visualisation, specify:
- What it shows
- Library to use (matplotlib/seaborn/MolViewSpec)
- What the user should understand from it
- Any pLDDT overlay logic

## Notebook Sections (with cell-level detail)
List each notebook section with:
- Section title and purpose
- What the user sees / reads
- What code is executed and what it produces

## User-Facing Explanations
Key explanatory text the notebook must include. Write the actual text 
(or a close draft) for:
- What EnzyMM is and how it works (1 paragraph, non-expert audience)
- What a catalytic motif match means (and doesn't mean) in a predicted 
  structure
- How to interpret the confidence overlay
- The caveat section on AlphaFold limitations

## Flywheel Integration
If applicable: what results should be submittable back to the InsightFold 
pre-computed layer, in what format, and under what conditions.

## Out of Scope for v1
Explicit list of features deferred to later versions, with one-line 
rationale for each.

## Open Questions Requiring Resolution Before Build
Anything still blocking implementation.

---

Produce complete, publication-quality spec text. Do not use placeholder 
language like "[describe visualisation here]" — write the actual content 
based on the decisions in the summary.