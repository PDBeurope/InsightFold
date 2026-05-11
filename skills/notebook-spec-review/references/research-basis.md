# Research Basis

This review skill adapts specification review practices from:

- Kiro Specs and Kiro best-practice guidance: requirements, design, and tasks should be version-controlled, traceable, reviewable, and task-sized; source: https://kiro.dev/docs/specs/ and https://kiro.dev/docs/specs/best-practices/
- GitHub Spec Kit: specs should be refined before planning and implementation; source: https://github.github.io/spec-kit/quickstart.html
- Diataxis: documentation should deliberately cover tutorials, how-to guides, reference, and explanation; source: https://diataxis.fr/
- Reproducible notebook guidance: notebooks must be readable, runnable, dependency-aware, and checked by restart/run-all; source: https://notebooks.genepattern.org/best-practices/

InsightFold adds notebook-specific gates:

- fixture manifest
- explicit data contracts
- variable handoff and hidden-state checks
- scientific assumption review
- target-runtime dependency audit
