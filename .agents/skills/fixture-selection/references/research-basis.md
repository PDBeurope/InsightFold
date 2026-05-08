# Research Basis

Fixture selection borrows from:

- Kiro's recommendation to keep specs version-controlled and scoped per feature, with tasks and validation tied to requirements; source: https://kiro.dev/docs/specs/best-practices/
- nbval's model of validating notebook execution against stored expected outputs; source: https://nbval.readthedocs.io/en/latest/index.html
- papermill's parameterized execution model for running notebooks against different inputs; source: https://papermill.readthedocs.io/en/latest/
- reproducible notebook guidance that data must be documented, available, or tiered so readers can reproduce at least downstream parts of the workflow; source: https://notebooks.genepattern.org/best-practices/

InsightFold fixture policy:

- Fixtures belong in the spec, not as late implementation choices.
- Each fixture needs provenance and expected outputs.
- At least one fixture must make the notebook runnable end-to-end.
- Edge and negative fixtures are required when the notebook claims robust behavior.
