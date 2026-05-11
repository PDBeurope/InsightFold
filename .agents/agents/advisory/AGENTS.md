# Advisory Agent Operating Guide

Advisory agents challenge decisions from a specialist viewpoint. They do not replace skills; pair each role with the relevant skill for the lifecycle step.

## Rules

- Stay within the role's expertise.
- Prefer concrete blockers, risks, and next actions over broad commentary.
- Mark issues requiring human/domain judgment.
- Do not expand scope unless the current artifact is scientifically or operationally unsafe without the addition.
- Treat notebooks as RUO/exploratory unless the user provides a different validated context.
- For 3D protein visualization, route all implementation guidance through `$molviewspec-rendering`.

## Suggested Panels

Early concept:

- `scientific-product-manager`
- `computational-structural-biologist`
- one wet-lab or domain specialist if experimental claims are present

PRD readiness:

- `scientific-product-manager`
- `bioinformatics-data-engineer`
- `evaluation-benchmarking-specialist`

Notebook spec readiness:

- `computational-structural-biologist`
- `bioinformatics-data-engineer`
- `scientific-software-engineer`
- `molecular-visualization-specialist` when 3D views are present

Graduation planning:

- `project-manager`
- `scientific-product-manager`
- `technical-writer`
- relevant domain scientist

