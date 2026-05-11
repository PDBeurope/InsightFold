# InsightFold Advisory Roles

Use these advisory agents for ideation, scoping, PRD review, spec review, notebook review, and graduation planning. They are role instructions, not auto-invoked skills.

Prefer small panels. Invoke only the expertise needed for the current decision. Startup roles such as Founder, CTO, or investor are intentionally excluded from the default panel; use them only if the user is making commercialization or organization-design decisions.

## Core Product And Delivery Roles

| Agent | Use For | Skills To Pair |
|---|---|---|
| `scientific-product-manager.md` | user problem, scope, PRD readiness, lifecycle fit | `$idea-scoping-interview`, `$scoping-decision-capture`, `$concept-to-prd` |
| `project-manager.md` | delivery plan, integration planning, graduation readiness, database handoff | `$notebook-review` |
| `machine-learning-engineer.md` | ML method choice, baselines, embeddings, model limits | `$idea-scoping-interview`, `$notebook-spec-review` |
| `scientific-software-engineer.md` | repository structure, notebook maintainability, tests, code quality | `$notebook-from-spec`, `$notebook-execution-validation` |
| `bioinformatics-data-engineer.md` | APIs, identifiers, data contracts, provenance, fixtures | `$fixture-selection`, `$prd-to-notebook-spec` |
| `evaluation-benchmarking-specialist.md` | metrics, baselines, leakage, tolerances, validation evidence | `$notebook-spec-review`, `$notebook-execution-validation` |
| `technical-writer.md` | notebook pedagogy, README, figure captions, limitations | `$notebook-review` |

## Biology And Scientific Roles

| Agent | Use For | Skills To Pair |
|---|---|---|
| `computational-structural-biologist.md` | structural interpretation, chain/residue assumptions, overclaims | `$notebook-spec-review`, `$notebook-review` |
| `protein-design-scientist.md` | design-adjacent framing, mutation effects, developability relevance | `$idea-scoping-interview`, `$notebook-review` |
| `wet-lab-liaison.md` | experimental realism, assay feasibility, validation constraints | `$idea-scoping-interview`, `$concept-to-prd` |
| `protein-biochemist.md` | protein mechanism, stability, binding, cofactors, oligomerization | `$idea-scoping-interview`, `$notebook-spec-review` |
| `enzymologist.md` | enzyme kinetics, active sites, substrate/product interpretation | `$idea-scoping-interview`, `$notebook-spec-review` |
| `protein-purification-scientist.md` | construct design, expression, purification, aggregation, QC | `$idea-scoping-interview`, `$concept-to-prd` |
| `biophysical-assay-scientist.md` | SPR/ITC/DSF/MST/SEC-MALS/CD evidence and assay caveats | `$idea-scoping-interview`, `$notebook-review` |
| `molecular-visualization-specialist.md` | MolViewSpec/Mol* views, 3D interpretation, viewer UX | `$molviewspec-rendering` |

## Invocation Pattern

Use one to four agents for most work:

```text
Act as .agents/agents/advisory/scientific-product-manager.md and use $concept-to-prd to review this concept capture before PRD writing.
```

For panel review, ask each agent for:

- highest-risk assumption
- missing evidence
- concrete recommendation
- whether the lifecycle should proceed, pause, or change scope
