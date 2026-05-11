# InsightFold

InsightFold is a notebook-driven development project for turning ideas about AlphaFold Database (AFDB), PDBe, protein structures, variants, interfaces, confidence metrics, and related biological questions into reviewed, runnable, and eventually reusable notebooks.

The goal is not just to make individual notebooks. The goal is to create a repeatable lifecycle for scientific notebook development:

1. capture an idea
2. decide whether it is worth pursuing
3. scope it into a PRD
4. convert the PRD into a notebook specification
5. review the specification
6. select pinned fixtures
7. build the notebook
8. execute and validate it
9. review it as a scientific artifact
10. beta test it
11. decide whether to integrate, maintain, or archive it

This matters because notebooks can easily become one-off experiments with hidden assumptions, fragile dependencies, unclear scientific claims, and no obvious path to maintenance. InsightFold aims to make notebook work faster while also making it easier to review, reproduce, explain, and hand off.

## What InsightFold Hopes To Achieve

InsightFold should help the project team:

- move from rough biological/product ideas to structured notebook prototypes
- keep notebook scope explicit before implementation starts
- make scientific assumptions, formulas, thresholds, and data sources traceable
- use pinned AFDB, PDBe, PDB, UniProt, or local fixtures instead of ad hoc examples
- validate notebooks by restart-and-run-all execution, not only visual inspection
- distinguish "the notebook runs" from "the notebook is scientifically and user-facingly good"
- decide whether a notebook should graduate into AFDB/PDBe, remain a standing notebook, or be archived

## Repository Layout

Key project areas:

```text
.agents/skills/
  idea-scoping-interview/
  scoping-decision-capture/
  concept-to-prd/
  prd-to-notebook-spec/
  notebook-spec-review/
  fixture-selection/
  notebook-from-spec/
  notebook-execution-validation/
  notebook-review/

.agents/agents/
  lifecycle/
    AGENTS.md
    spec-reviewer.md
    fixture-curator.md
    notebook-builder.md
    notebook-validator.md
    notebook-reviewer.md
  homodimer/
  advisory/

prd/
specs/
notebooks/
src/
skills/          # legacy/source copy retained for comparison
agents/          # legacy/source copy retained for comparison
agent-skills/
```

`.agents/skills/` contains reusable Codex skill instructions. Each skill has a `SKILL.md` file and, where useful, templates under `assets/` and rationale under `references/`.

`.agents/agents/lifecycle/` contains role definitions for agents that apply those skills at specific lifecycle stages. These are not standalone programs; they are operating instructions for how an AI agent should behave when assigned that role.

`prd/`, `specs/`, and `notebooks/` are the expected artifact path from product idea to notebook implementation.

`skills/`, `agents/`, and `agent-skills/` are retained as original/source material. New Codex-facing work should use `.agents/`.

`agent-skills/insightfold-lifecycle-skill-roadmap.md` is the broader lifecycle roadmap, including future skills that have not yet been implemented.

## Current Lifecycle Chain

The current implemented chain is:

```text
idea-scoping-interview
  -> scoping-decision-capture
  -> concept-to-prd
  -> prd-to-notebook-spec
  -> spec-reviewer
  -> fixture-curator
  -> notebook-builder
  -> notebook-validator
  -> notebook-reviewer
```

The corresponding skills are:

| Stage | Skill | Purpose |
|---|---|---|
| Concept discovery | `.agents/skills/idea-scoping-interview/SKILL.md` | Interview and scope a rough notebook idea |
| Decision capture | `.agents/skills/scoping-decision-capture/SKILL.md` | Summarize decisions, tradeoffs, and open questions |
| PRD creation | `.agents/skills/concept-to-prd/SKILL.md` | Convert a scoped concept into a true PRD |
| Notebook spec | `.agents/skills/prd-to-notebook-spec/SKILL.md` | Convert the PRD into an implementation-ready notebook spec pack |
| Spec review | `.agents/skills/notebook-spec-review/SKILL.md` | Review the spec before implementation starts |
| Fixture selection | `.agents/skills/fixture-selection/SKILL.md` | Choose pinned examples, edge cases, expected outputs, and provenance |
| Notebook build | `.agents/skills/notebook-from-spec/SKILL.md` | Build the notebook from the reviewed spec and fixture manifest |
| Execution validation | `.agents/skills/notebook-execution-validation/SKILL.md` | Validate restart-and-run-all behavior, fixture outputs, dependencies, and reproducibility |
| Final notebook review | `.agents/skills/notebook-review/SKILL.md` | Review scientific quality, interpretation, pedagogy, maintainability, and lifecycle readiness |

The corresponding agent roles are:

| Agent role | File | Responsibility |
|---|---|---|
| Shared lifecycle guidance | `.agents/agents/lifecycle/AGENTS.md` | Common rules for all lifecycle agents |
| Spec reviewer | `.agents/agents/lifecycle/spec-reviewer.md` | Blocks vague, incomplete, or untestable specs |
| Fixture curator | `.agents/agents/lifecycle/fixture-curator.md` | Selects fixtures and expected outputs |
| Notebook builder | `.agents/agents/lifecycle/notebook-builder.md` | Implements notebooks from reviewed specs |
| Notebook validator | `.agents/agents/lifecycle/notebook-validator.md` | Runs mechanical validation and reports evidence |
| Notebook reviewer | `.agents/agents/lifecycle/notebook-reviewer.md` | Reviews whether a notebook is ready to share, beta, graduate, or archive |

## Recommended Way Of Working

Use one lead thread first, then parallel subagents only after the core question is stable.

Mental model:

- `skill` = reusable method
- `agent` = role lens, priorities, and guardrails
- `spawned subagent` = a separate worker that can use both

Recommended sequence in a new conversation:

1. Start with one lead role in a single thread.
   Use `.agents/agents/advisory/scientific-product-manager.md` with `$idea-scoping-interview`.
   Goal: define the user, the biological question, the expected evidence, and the rough notebook outcome.

2. Do not parallelize yet if the core sentence is still unclear.
   The idea should be stable enough to say:
   `this notebook is for X user to answer Y question using Z evidence`

3. Once that sentence is stable, spawn 3 to 5 advisory subagents in parallel.
   Good default panel:
   - `scientific-product-manager`
   - `computational-structural-biologist`
   - `bioinformatics-data-engineer`
   - `evaluation-benchmarking-specialist`
   - optional domain specialist such as `enzymologist`, `wet-lab-liaison`, `protein-biochemist`, or `molecular-visualization-specialist`

4. Give every subagent the same idea statement, but a different role.
   Each subagent should return:
   - key opportunities
   - highest-risk assumptions
   - blocking questions
   - concrete recommendations
   - whether the idea is ready for `$scoping-decision-capture`

5. Synthesize the panel with `$scoping-decision-capture`.
   This is the bridge from exploratory discussion to a structured concept record.

6. Only then use `$concept-to-prd`.
   The PRD should stay a true PRD. Do not push notebook section order, cell plans, traceability tables, or fixture manifests into the PRD.

7. Use `$prd-to-notebook-spec` after the PRD is approved.
   This skill owns the implementation-ready spec pack.

8. Use lifecycle agents after the spec exists.
   The normal order is:
   - `spec-reviewer`
   - `fixture-curator`
   - `notebook-builder`
   - `notebook-validator`
   - `notebook-reviewer`

When to stay in one thread:

- very early ideation
- clarifying the user and problem
- deciding whether an idea is worth pursuing at all

When to use spawned subagents:

- once the idea has enough shape for role-specific critique
- when you want parallel expert views
- when you want role separation instead of one blended answer

When not to use too many subagents:

- if the idea is still vague
- if multiple roles would ask the same basic questions
- if you do not yet know which biological or technical domain matters most

## How To Use The Skills From Start To Finish

Use the lifecycle as a sequence of gates. Do not jump straight from idea to notebook unless the work is intentionally exploratory.

### 1. Start With An Idea

Write a short idea in plain language:

```text
Idea: Build a notebook that helps users inspect whether AFDB structure confidence changes around clinically interesting variant positions.
```

Recommended first prompt:

```text
Use $idea-scoping-interview.

Act as .agents/agents/advisory/scientific-product-manager.md.

I want to explore a new InsightFold notebook idea:
<idea>

Do not write a PRD yet.
Interview me until the user, the biological question, the evidence, and the rough notebook outcome are clear.
```

Expected outcome:

- stable problem statement
- target user
- expected notebook output
- obvious risks and unknowns
- judgment about whether parallel advisory review should start

### 2. Run Parallel Advisory Review

After the core question is stable, spawn advisory subagents in parallel.

Recommended default panel:

- `.agents/agents/advisory/scientific-product-manager.md`
- `.agents/agents/advisory/computational-structural-biologist.md`
- `.agents/agents/advisory/bioinformatics-data-engineer.md`
- `.agents/agents/advisory/evaluation-benchmarking-specialist.md`
- optional: `.agents/agents/advisory/molecular-visualization-specialist.md` when 3D structure views matter

Recommended prompt:

```text
Spawn parallel subagents for advisory review.

Shared idea statement:
<stable idea statement>

Subagent 1:
Act as .agents/agents/advisory/scientific-product-manager.md and use $idea-scoping-interview.
Focus on scope, user value, PRD readiness, and lifecycle fit.

Subagent 2:
Act as .agents/agents/advisory/computational-structural-biologist.md.
Focus on biological assumptions, interpretation risk, and overclaims.

Subagent 3:
Act as .agents/agents/advisory/bioinformatics-data-engineer.md.
Focus on APIs, identifiers, provenance, fixtures, and reproducibility.

Subagent 4:
Act as .agents/agents/advisory/evaluation-benchmarking-specialist.md.
Focus on validation, expected outputs, fixture quality, and failure cases.

Have each subagent return:
- key opportunities
- highest-risk assumptions
- blocking questions
- concrete recommendations
- whether this is ready for $scoping-decision-capture
```

Expected outcome:

- role-specific critique
- surfaced blockers
- better-defined scope boundaries
- clear recommendation on whether to proceed

### 3. Capture The Decisions

Use:

```text
.agents/skills/scoping-decision-capture/SKILL.md
```

Ask Codex:

```text
Use $scoping-decision-capture to synthesize the ideation discussion and advisory-panel outputs into a structured decision capture.
Do not write the PRD yet.
```

Expected output:

- core problem statement
- scope decisions
- user journey summary
- open questions
- success and stop signals

### 4. Convert The Concept To A PRD

Use:

```text
.agents/skills/concept-to-prd/SKILL.md
```

Ask Codex something like:

```text
Use $concept-to-prd to turn this scoped concept into an InsightFold PRD:
<decision capture or approved scoped concept>
```

Expected output:

```text
prd/<feature>.md
```

The PRD should define the user problem, target audience, scope, non-goals, success criteria, assumptions, risks, and expected notebook artifact.

### 5. Convert The PRD To A Notebook Spec Pack

Use:

```text
.agents/skills/prd-to-notebook-spec/SKILL.md
```

Ask Codex:

```text
Use $prd-to-notebook-spec to convert prd/<feature>.md into a notebook spec pack.
```

Expected output:

```text
specs/<feature>/
  spec-pack-overview.md
  requirements.md
  notebook-ux-contract.md
  notebook-design.md
  cell-blueprint.md
  traceability-matrix.md
  tasks.md
  validation.md
  docs-plan.md
  fixture-manifest.md
  data-contracts.md
```

For very small prototypes, a single consolidated spec is acceptable if it contains equivalent sections.

### 6. Review The Spec Before Building

Use:

```text
.agents/skills/notebook-spec-review/SKILL.md
.agents/agents/lifecycle/spec-reviewer.md
```

Ask Codex:

```text
Act as .agents/agents/lifecycle/spec-reviewer.md and use $notebook-spec-review to review specs/<feature>/ for implementation readiness.
```

The review should decide whether implementation can start. It should flag:

- missing acceptance criteria
- vague tasks
- weak fixtures
- missing data contracts
- scientific ambiguity
- dependency risks
- unclear validation
- documentation gaps

Do not build the notebook until blocking spec review findings are resolved.

### 7. Select Fixtures

Use:

```text
.agents/skills/fixture-selection/SKILL.md
.agents/agents/lifecycle/fixture-curator.md
```

Ask Codex:

```text
Act as .agents/agents/lifecycle/fixture-curator.md and use $fixture-selection to create a fixture manifest for specs/<feature>/.
```

Expected output:

```text
specs/<feature>/fixture-manifest.md
```

Good fixture manifests include:

- happy-path fixture
- edge-case or negative fixture where relevant
- stable identifiers or local paths
- source/provenance
- retrieval date for network data
- expected outputs
- tolerances
- validation checks

### 8. Build The Notebook From The Spec

Use:

```text
.agents/skills/notebook-from-spec/SKILL.md
.agents/agents/lifecycle/notebook-builder.md
```

Ask Codex:

```text
Act as .agents/agents/lifecycle/notebook-builder.md and use $notebook-from-spec to implement the notebook described in specs/<feature>/.
```

Expected output:

```text
notebooks/<feature>.ipynb
```

The notebook should:

- run top-to-bottom after kernel restart
- have explicit setup/import cells
- use clear markdown sections
- keep code cells focused
- validate required fields before use
- show provenance for fetched or uploaded data
- include interpretation and limitations
- avoid hidden local paths and hidden state

### 9. Validate Notebook Execution

Use:

```text
.agents/skills/notebook-execution-validation/SKILL.md
.agents/agents/lifecycle/notebook-validator.md
```

Ask Codex:

```text
Act as .agents/agents/lifecycle/notebook-validator.md and use $notebook-execution-validation to validate notebooks/<feature>.ipynb against specs/<feature>/.
```

Expected output:

```text
specs/<feature>/validation-report.md
```

Validation should check:

- restart-and-run-all execution
- fixture outputs
- expected warnings or failures
- data contracts
- dependency/runtime constraints
- hidden state
- unresolved TODOs in critical cells
- visualization outputs
- documentation presence

Execution validation only proves the notebook runs and matches its mechanical checks. It does not prove the notebook is scientifically ready.

### 10. Review The Notebook As A Scientific Artifact

Use:

```text
.agents/skills/notebook-review/SKILL.md
.agents/agents/lifecycle/notebook-reviewer.md
```

Ask Codex:

```text
Act as .agents/agents/lifecycle/notebook-reviewer.md and use $notebook-review to review notebooks/<feature>.ipynb after execution validation.
```

Expected output:

```text
specs/<feature>/notebook-review.md
```

The review should assess:

- scientific correctness
- formulas and cited sources
- biological assumptions
- uncertainty and limitations
- reproducibility
- pedagogy and documentation
- visualization quality
- maintainability
- readiness for beta, graduation, or archive

Use a human/domain review gate when interpretation, thresholds, clinical/RUO framing, or AFDB/PDBe publication risk requires judgment.

### 11. Move To Instrumented Beta

Current status: placeholder stage.

Placeholder skills to develop:

- `beta-feedback-instrumentation`
- `usage-signal-summary`
- `beta-release-packaging`

Expected future outputs:

```text
specs/<feature>/beta-plan.md
specs/<feature>/feedback-summary.md
specs/<feature>/beta-release-notes.md
```

The beta stage should answer:

- who tried the notebook
- whether they understood it
- whether they returned to it
- where they failed
- what they preferred over existing workflows
- what must change before graduation

### 12. Graduation Review

Current status: placeholder stage.

Placeholder skills to develop:

- `graduation-decision-brief`
- `engineering-handoff-package`
- `standing-notebook-maintenance-plan`
- `archive-retrospective`

Expected future outputs:

```text
specs/<feature>/graduation-decision.md
specs/<feature>/engineering-handoff.md
specs/<feature>/maintenance-plan.md
specs/<feature>/archive-retrospective.md
```

Possible outcomes:

- integrate into AFDB or PDBe
- keep as a standing InsightFold notebook
- archive with rationale and preserved state

## Recommended Prompts

### Full Ideation

```text
Use $idea-scoping-interview.

Act as .agents/agents/advisory/scientific-product-manager.md.

I want to explore a new InsightFold notebook idea:
<idea>

Do not write a PRD yet.
Interview me until the user, the question, the evidence, and the rough notebook outcome are clear.
```

### Parallel Advisory Review

```text
Spawn parallel subagents for advisory review of this stable idea:
<stable idea statement>

Use:
- .agents/agents/advisory/scientific-product-manager.md
- .agents/agents/advisory/computational-structural-biologist.md
- .agents/agents/advisory/bioinformatics-data-engineer.md
- .agents/agents/advisory/evaluation-benchmarking-specialist.md

Have each subagent return:
- key opportunities
- highest-risk assumptions
- blocking questions
- concrete recommendations
- whether this is ready for $scoping-decision-capture
```

### Ideation To PRD

```text
Use $scoping-decision-capture to synthesize the ideation discussion and advisory-panel outputs.
Then use $concept-to-prd only if the concept is ready.
```

### Full Lifecycle

```text
Take this idea through the InsightFold lifecycle.

1. Use $idea-scoping-interview in one lead thread first.
2. When the idea is stable, spawn parallel advisory subagents.
3. Use $scoping-decision-capture to synthesize the discussion.
4. Use $concept-to-prd only if the concept is ready.
5. Use $prd-to-notebook-spec to build the spec pack.
6. Then use the lifecycle agents in order.
7. Stop at each gate if there are blocking findings.

Idea:
<idea>
```

### From Existing Spec

```text
Act as .agents/agents/lifecycle/fixture-curator.md and complete fixture selection for specs/<feature>/.
Then act as .agents/agents/lifecycle/notebook-builder.md and build the notebook only if fixtures and validation criteria are ready.
```

### Validate Existing Notebook

```text
Act as .agents/agents/lifecycle/notebook-validator.md and use $notebook-execution-validation to validate notebooks/<feature>.ipynb against specs/<feature>/.
```

### Final Review

```text
Act as .agents/agents/lifecycle/notebook-reviewer.md and use $notebook-review to decide whether notebooks/<feature>.ipynb is ready for beta, graduation review, continued iteration, or archive.
```

## Development Status

Implemented lifecycle skills:

- `idea-scoping-interview`
- `scoping-decision-capture`
- `concept-to-prd`
- `prd-to-notebook-spec`
- `notebook-spec-review`
- `fixture-selection`
- `notebook-from-spec`
- `notebook-execution-validation`
- `notebook-review`

Implemented lifecycle agents:

- `spec-reviewer`
- `fixture-curator`
- `notebook-builder`
- `notebook-validator`
- `notebook-reviewer`

Important placeholders:

- `idea-capture`
- `idea-quality-check`
- `idea-to-triage-brief`
- `evidence-scan`
- `risk-and-assumption-log`
- `acceptance-criteria-authoring`
- `dependency-policy`
- `structure-parsing-strategy`
- `data-contract-validation`
- `scientific-computation-patterns`
- `notebook-result-interpretation`
- `scientific-figure-style`
- `fixture-regression-harness`
- `dependency-runtime-audit`
- `beta-feedback-instrumentation`
- `usage-signal-summary`
- `beta-release-packaging`
- `graduation-decision-brief`
- `engineering-handoff-package`
- `standing-notebook-maintenance-plan`
- `archive-retrospective`

## Guiding Principle

InsightFold should make notebook development faster, but speed is not the only target. The lifecycle exists so that each notebook carries enough context, evidence, validation, and review to be useful beyond the first experiment.
