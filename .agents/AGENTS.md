# InsightFold Repo Agents And Skills

This repository keeps Codex-discoverable skills and repo-local agent role documents under `.agents/`.

## Directory Contract

```text
.agents/
  skills/
    <skill-name>/
      SKILL.md
      agents/openai.yaml
      assets/
      references/
      scripts/
  agents/
    lifecycle/
    homodimer/
    advisory/
```

`SKILL.md` frontmatter must contain only `name` and `description`. Put inputs, outputs, workflows, quality gates, and examples in the Markdown body.

Agent role documents are Markdown operating instructions. They are not auto-invoked like skills; invoke them explicitly by path and pair them with the relevant `$skill-name`.

## Invocation

Use skills directly:

```text
Use $prd-to-notebook-spec to convert prd/<feature>.md into a spec pack.
```

Use an agent role plus a skill:

```text
Act as .agents/agents/lifecycle/spec-reviewer.md and use $notebook-spec-review to review specs/<feature>/.
```

Use the homodimer role suite when building the AFDB homodimer diagnostic notebook:

```text
Act as .agents/agents/homodimer/orchestrator.md and use $homodimer-confidence-scoring with the domain implementation skills.
```

## InsightFold Lifecycle

The project follows notebook-driven development:

```text
idea -> triage -> scope -> notebook prototype -> instrumented beta -> graduation review
```

Graduation outcomes are:

- integrate into AFDB or PDBe
- keep as a standing notebook
- archive with rationale

Supporting workstreams are:

- development harness
- feedback
- documentation

The current implemented repo chain starts at PRD creation and continues through notebook review:

```text
$convert-to-prd
  -> $prd-to-notebook-spec
  -> .agents/agents/lifecycle/spec-reviewer.md + $notebook-spec-review
  -> .agents/agents/lifecycle/fixture-curator.md + $fixture-selection
  -> .agents/agents/lifecycle/notebook-builder.md + $notebook-from-spec
  -> .agents/agents/lifecycle/notebook-validator.md + $notebook-execution-validation
  -> .agents/agents/lifecycle/notebook-reviewer.md + $notebook-review
```

Do not treat a runnable notebook as scientifically approved. Execution validation and qualitative scientific review are separate gates.
