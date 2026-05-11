# Research Basis

Notebook implementation from spec borrows from:

- GitHub Spec Kit: implementation should follow clarified specs and generated plans/tasks rather than one-shot prompts; source: https://github.github.io/spec-kit/quickstart.html
- Kiro Specs: requirements, design, and tasks provide implementation context and task tracking; source: https://kiro.dev/docs/specs/
- nbdev: notebooks can combine code, tests, and documentation, but this requires discipline around cells and docs; source: https://pypi.org/project/nbdev/
- Jupytext: text notebooks can improve version control and diffs; source: https://jupytext.readthedocs.io/
- papermill: notebooks can be parameterized and executed programmatically; source: https://papermill.readthedocs.io/en/latest/

InsightFold adaptation:

- The spec pack is the implementation contract.
- Notebook sections should map to requirements.
- Fixtures and data contracts are first-class implementation inputs.
- Validation is not complete until a separate execution-validation pass runs.
