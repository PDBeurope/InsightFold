---
name: prd-notebook-generation
description: Generate the PRD homodimer confidence diagnostic notebook from the spec using agent workflows and educational visualizations.
---

# Goal
Create a reproducible Jupyter notebook for the PRD homodimer confidence metric diagnostic use case.

# Requirements
- follow the PRD spec in `specs/homodimer_diagnostic_notebook_spec.md`
- implement AFDB API fetching, mmCIF parsing, PAE parsing, and confidence score formulas
- include ipTM, ipSAE variants, pDockQ, pDockQ2, and LIS
- use beginner-friendly explanations and annotated visual diagnostics
- keep dependencies minimal and notebook executable

# Workflow
1. load AFDB prediction data
2. parse atomic coordinates and PAE/confidence files
3. detect the inter-chain interface
4. recompute metrics and score variants
5. visualise PAE, contacts, and interfaces
6. summarise the diagnostic findings for a non-expert audience
