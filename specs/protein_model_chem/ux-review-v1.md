# UX Review — Protein Model Chemistry Notebook v1

| Field | Value |
|---|---|
| Notebook | `notebooks/protein_model_chem.ipynb` |
| Reviewer | Claude Code |
| Date | 2026-05-08 |
| Status | findings-and-recommendation |

---

## Summary

The v1 notebook was built as a fixture-driven test harness rather than as the educational,
user-facing tool described in the PRD. A user cannot currently enter their own protein and
mutation without editing Python source. This report identifies the three compounding causes and
recommends the minimal change set to reach the intended UX.

---

## Root Cause Analysis

### Cause 1 — `$fixture-selection` scope was validation, not UX

The skill's purpose is to pick representative proteins for correctness testing: a happy-path
case, edge cases for confidence downgrading, interface perturbation, unsupported context, and
negatives for invalid IDs and mutation mismatches. It did that correctly, producing FX-001
through FX-007.

The problem is what happened next. The `$notebook-from-spec` skill wired those fixtures as the
**primary interface** of the notebook. The setup cell exposes:

```python
FIXTURE_ID = "FX-001"   # change this to run another fixture
```

This makes the notebook behave like a test runner. A user who wants to analyse their own
protein has no obvious entry point.

### Cause 2 — `$notebook-spec-review` did not challenge the educational UX gap

The spec review checked traceability, scientific defensibility, and implementation readiness
thoroughly. None of the seven human-review questions (H-001 through H-007) asked how a
non-fixture user would actually input a protein or mutation. The relevant design spec line
(Section 3 of `notebook-design.md`) reads:

> Accept AFDB accession, PDBe identifier, or optional local path → `structure_request`

This was never challenged as an implementation requirement during review. It became a dict the
user must manually edit, not a surfaced input cell.

### Cause 3 — PRD intent diluted through the spec pipeline

The PRD is explicit about the user experience:

- UC-01: _"User can provide a structure accession or supported structure file."_
- UC-04: _"User can specify at least one mutation with chain, residue, original residue, and
  substituted residue."_
- Section 14 (UX Requirements): _"Inputs should be explicit and validated before analysis."_
- Story 1–3: Structural biologists and experimentalists entering their own proteins and mutations.

That intent survived into `notebook-design.md` but was not implemented as input widgets. Cell 24
of the built notebook acknowledges the workaround in a footnote:

> To analyze a new structure, update `structure_request` with `source`, `identifier`, and
> `mutation_notation`, then restart-run-all.

This buries the user's entry point in reuse documentation rather than surfacing it as the first
thing they see.

---

## What Is Working Well

The analysis cells are well-built and reusable. The following are correct and should not be
changed:

- Structure retrieval and provenance recording.
- BioPython parsing and residue/atom normalization.
- Interaction detection functions (hydrogen bonds, salt bridges, hydrophobic, aromatic,
  disulfide, steric clashes, VdW packing).
- Interaction graph construction and network metrics.
- Mutation perturbation and WT-vs-mutant comparison tables.
- Mechanistic summary with observation / interpretation / hypothesis separation.
- Confidence propagation and severity labelling.
- MolViewSpec visualization with graceful fallback.
- Export bundle and validation snapshot sections.

The only layer that needs reworking is the **input layer**.

---

## Recommended Fix — v2 Input Layer

### 1. Add a User Input section as the first code cell

Use `ipywidgets` (already permitted by `CLAUDE.md` dependency policy) to expose:

| Widget | Type | Default | Notes |
|---|---|---|---|
| Protein source | Dropdown | `PDBe` | Options: `PDBe`, `AlphaFold` |
| Protein identifier | Text | `2LZM` | PDB ID or UniProt accession; show format hint |
| Mutation | Text | `A:L99A` | Format hint: `CHAIN:WTresNUMMUTres` e.g. `A:L99A` or `R273H` |
| Run button | Button | — | Triggers or documents restart-and-run-all |

Render examples from the fixture registry below the widgets so users can click-to-fill a known
case when learning.

### 2. Feed widget values into the existing pipeline variables

No analysis code should change. The widgets write into:

```python
structure_request = {
    "source": source_widget.value.lower(),
    "identifier": identifier_widget.value.strip(),
}
mutation_notation = mutation_widget.value.strip()
```

These are already the variables consumed by retrieval and mutation validation cells.

### 3. Demote fixtures to an Examples appendix

Move `FIXTURES` and `FIXTURE_ID` out of the setup cell and into a clearly labelled
"Examples and Validation Fixtures" section near the end of the notebook. Fixtures become
pre-filled suggestions a learner can copy into the input widgets, not the mechanism that drives
analysis.

### 4. Update Section 1 framing markdown

The opening markdown should lead with:
_"Enter your protein identifier and mutation below to begin."_
Not with fixture-selection instructions.

---

## Scope of Change

| Area | Change Needed | Effort |
|---|---|---|
| Cell 0 (intro markdown) | Reframe opening to lead with user input | Small |
| Cell 2 (setup) | Remove `FIXTURE_ID`; add `ipywidgets` input section | Medium |
| Cells 3–19 (analysis pipeline) | None | None |
| Cell 22–23 (validation snapshot) | Decouple fixture-specific checks from general checks | Small |
| Cell 24 (reuse notes) | Replace fixture-change instructions with input widget guidance | Small |
| New appendix cell | Move `FIXTURES` dict here as examples for learners | Small |

---

## Process Recommendation

To prevent this pattern in future notebooks, the `$notebook-spec-review` skill should include
an explicit UX readiness check:

> **UX-01**: Is there a user-facing input cell that does not require the user to edit function
> definitions or internal pipeline variables?

This check would have caught the gap at review time, before implementation.

---

## References

- `specs/protein_model_chem/protein_model_chem_prd.md` — original PRD, Sections 8, 9, 14
- `specs/protein_model_chem/notebook-design.md` — Section 3 (Structure input), Section 6 (Mutation input)
- `specs/protein_model_chem/spec-review.md` — human review questions H-001 through H-007
- `specs/protein_model_chem/fixture-manifest.md` — FX-001 through FX-007
- `notebooks/protein_model_chem.ipynb` — built notebook, cells 2 and 24
