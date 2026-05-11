# Using agents and skills to build upon a PRD

This directory is retained as historical/source material. Canonical reusable skills now live under the repository-level `skills/` directory.

The MolViewSpec rendering skill has moved to `skills/molviewspec-rendering/SKILL.md`.

## Step 1

Given a specific PRD (see [homodimer_diagnostic_notebook_spec](https://github.com/PDBeurope/InsightFold/blob/main/specs/homodimer_diagnostic_notebook_spec.md)), asked an AI assistant (**Codex**, Claude, Copilot): *which agents and skills do you need to make it happen?*

**Codex reply:**

```sh
.agents/
├── AGENTS.md
├── agents/
│   ├── orchestrator.md
│   ├── structural-biologist.md
│   ├── visualization-engineer.md
│   ├── scientific-python-engineer.md
│   └── notebook-reviewer.md
└── skills/
    ├── afdb-api-fetching/
    │   └── SKILL.md
    ├── ipsae-implementation/
    │   └── SKILL.md
    ├── pdockq-implementation/
    │   └── SKILL.md
    ├── lis-implementation/
    │   └── SKILL.md
    ├── pae-visualization/
    │   └── SKILL.md
    ├── notebook-assembly/
    │   └── SKILL.md
    └── validation-against-reference/
        └── SKILL.md
```

## Step 2

Asked AI assistant (Codex, Claude, **Copilot**) to populate the *.md files with what it think actually fits best to the agents/skills. 

NOTE: this is an important step, should be reviewed by a domain expert to get a better performance of the AI assistant.

## Step 3 

Asked AI assistant (Codex, Claude, **Copilot**) to create a Jupyter notebook for building the PRD using the agents and skills previsously created.

Jupyter notebook can be found in: `specs/homodimer_diagnostic_notebook.ipynb`
