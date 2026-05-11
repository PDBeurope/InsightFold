# Research Basis

This skill adapts current spec-driven development practices to InsightFold notebook-driven development.

## Sources Consulted

- GitHub Spec Kit documentation: emphasizes intent-driven development, "what before how", multi-step refinement, and specifications as implementation-driving artifacts. Source: https://github.github.io/spec-kit/
- GitHub Spec Kit quick start: uses a sequence of constitution, specify, clarify, plan, tasks, analyze, implement. Source: https://github.github.io/spec-kit/quickstart.html
- Kiro Specs documentation: defines the core spec structure as `requirements.md`, `design.md`, and `tasks.md`, with requirements/design/tasks as a three-phase workflow. Source: https://kiro.dev/docs/specs/
- Kiro Feature Specs documentation: supports requirements-first and design-first variants depending on whether product behavior or technical feasibility is the starting point. Source: https://kiro.dev/docs/specs/feature-specs/
- Diataxis documentation framework: separates documentation into tutorials, how-to guides, reference, and explanation based on user needs. Source: https://diataxis.fr/
- Jupyter Book / MyST executable documentation: supports executing notebook content, labeling outputs, and embedding/reusing computational outputs. Source: https://jupyterbook.org/stable/execution/execution/
- GenePattern notebook best practices and reproducible notebook guidance: emphasize short cells, modular code, explicit dependencies, version control, run-all validation, input data explanation, and enabling notebooks to be read, run, and explored. Source: https://notebooks.genepattern.org/best-practices/

## Adaptation for InsightFold

The best fit for InsightFold is not a plain software feature spec. Notebooks need additional spec artifacts:

- fixture manifest: the notebook must run on pinned examples
- data contracts: external APIs and uploaded files must be explicit
- notebook design: sections, cells, and variable handoffs matter because notebooks have hidden-state risks
- validation plan: top-to-bottom execution is a first-class acceptance criterion
- documentation plan: narrative, scientific interpretation, inputs/outputs, and user adaptation guidance must be planned before implementation

## Recommended Spec Spine

Use the spec-driven spine:

```text
requirements -> notebook design -> tasks
```

Extend it with notebook-specific documents:

```text
data contracts -> fixtures -> validation -> documentation plan
```

## Practical Rules

- Keep PRD and spec separate: PRD explains product intent; spec explains executable notebook behavior and implementation boundaries.
- Do not start tasks until blocking ambiguities are resolved.
- Treat fixture selection as part of specification, not an implementation detail.
- Require every major notebook section to trace back to a requirement.
- Require every external input to have a data contract.
- Require every spec pack to include validation and documentation plans.
- Prefer domain-specific reference files over making the skill itself domain-specific.
