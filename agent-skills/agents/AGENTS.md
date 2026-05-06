# Homodimer Diagnostic Notebook Project

## Project Goal

Build a scientifically rigorous, educational, reproducible diagnostic notebook for AlphaFold Database homodimers.

The notebook must:
- recompute confidence metrics from scratch
- explain score disagreements
- generate educational visualizations
- remain understandable to non-specialists
- match the IPSAE reference implementation

## Scientific Rules

- Never alter published formulas
- Preserve asymmetric PAE logic
- Validate against ipsae.py formulas
- Maintain floating-point reproducibility
- Prefer clarity over optimization
- Avoid hidden heuristics unless documented

## Engineering Rules

- Prefer numpy-first implementations
- Avoid GPU dependencies
- No BioPython unless explicitly required
- Keep notebook self-contained
- All notebook cells must run sequentially
- Prefer deterministic code paths

## Visualization Rules

- All figures require legends and labels
- Use consistent AFDB colour conventions
- All heatmaps require colour bars
- Prioritize interpretability over aesthetics
- Minimum font size: 12pt

## Educational Rules

The target audience is unfamiliar with structural bioinformatics.

Every notebook section must:
- explain concepts before code
- define jargon
- interpret visualizations
- connect formulas to intuition

## Validation Rules

All scores must match reference implementations within ±0.001.

Required validations:
- ipTM
- ipSAE_d0res
- ipSAE_d0chn
- ipSAE_d0dom
- pDockQ
- pDockQ2
- LIS

## Coding Style

- Use descriptive variable names
- Prefer small pure functions
- Document formulas inline
- Keep plotting logic modular
- Avoid excessive abstraction

## Project Structure

- agents/ contains specialized reasoning agents
- skills/ contains reusable workflows
- notebook generation should be modular
- reusable calculations belong in skills
