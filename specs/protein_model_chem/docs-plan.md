# Protein Model Chemistry Documentation Plan

Source PRD: `specs/protein_model_chem/protein_model_chem_prd.md`

## Documentation Goals

- Make the notebook runnable by structural biology users without requiring code expertise.
- Explain structural chemistry assumptions without overstating thermodynamic or clinical certainty.
- Preserve enough reference detail for implementation, validation, and domain review.
- Keep limitations near the outputs they qualify, not only in a final disclaimer.

## Documentation Matrix

| Documentation Type | Artifact | Audience | Notes |
|---|---|---|---|
| Tutorial | Notebook introduction and first runnable fixture walkthrough | New users, experimentalists | Show how to run the happy-path fixture and read the outputs. |
| How-to | Notebook section for changing accession/source/mutation | Returning users | Include examples for AFDB, PDBe, and optional local files if supported. |
| Reference | Input schema, output tables, interaction thresholds, graph metrics, confidence rules | Developers, reviewers, power users | Must align with `data-contracts.md`. |
| Explanation | Structural constraint networks, compatibility perturbation, confidence layers, limitations | All users | Explain observation vs interpretation vs hypothesis. |
| Review checklist | Domain review notes in notebook or companion markdown | Domain reviewers | Focus on thresholds, mutation model, severity language, and caveats. |

## Notebook Documentation Requirements

The notebook should include:

- Purpose and RUO framing.
- Audience assumptions: users understand protein chemistry and structural biology, but do not need to code or interpret raw force-field outputs.
- Input examples for AFDB, PDBe, and mutation notation.
- Source provenance explanation for AlphaFold and PDBe.
- Structure confidence explanation, including pLDDT and experimental caveats.
- Interaction class definitions with thresholds and limitations.
- Graph metric definitions and warnings against overinterpreting centrality.
- Mutation model assumptions and what the notebook does not model.
- Explanation of output layers: observation, interpretation, hypothesis.
- Structural disruption likelihood wording and non-clinical limitation.
- Export instructions and artifact descriptions.
- Fixture and validation notes.

## External Docs Or README Updates

Recommended after implementation:

- Add notebook entry to the project notebook index or README if one exists.
- Link to `specs/protein_model_chem/requirements.md` for scope and acceptance criteria.
- Link to `specs/protein_model_chem/data-contracts.md` for input/output reference.
- Link to final validation report after execution validation.

## Interpretation Language Guide

Preferred language:

- "suggests"
- "is consistent with"
- "may reduce local compatibility"
- "may require structural accommodation"
- "possible destabilization signal"
- "lower confidence because"
- "structural hypothesis"

Avoided language:

- "causes"
- "abolishes"
- "guarantees"
- "pathogenic"
- "benign"
- "disease severity"
- "diagnostic"
- "ddG"
- "unfolding risk" as a deterministic claim
- exact free-energy or total-stability claims

## Documentation Acceptance Criteria

- A user can run the happy-path fixture from the notebook instructions alone.
- A user can change accession/source/mutation using documented fields.
- A reviewer can find every interaction threshold and confidence rule.
- A reviewer can trace every summary statement to observation, interpretation, or hypothesis.
- Limitations are visible near relevant outputs.
- Documentation does not imply clinical, pathogenicity, exact thermodynamic, or deterministic predictive readiness.

## Deferred Documentation

- Production deployment documentation.
- Clinical interpretation guidance.
- Batch processing instructions.
- MD, force-field, or repacking engine documentation.
- Full API reference for a reusable package, unless notebook functions are later extracted into `src/insightfold/`.
