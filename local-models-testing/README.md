# LLM Spec Evaluator

This workspace contains a small notebook-style evaluator for comparing Markdown
project specifications written by smaller LLMs.

## Candidate File Convention

Put candidate specs in a directory such as `candidate_mds/` and name each file:

```text
<project_name>--<model_name>.md
```

Examples:

```text
spec_nmr_restraints--gemma_31b.md
spec_nmr_restraints--gemma_4eb.md
```

The explicit `--` separator lets project and model names safely contain hyphens
or underscores.

## Optional References

Reference specs are optional. When used, place them in a reference directory with
one Markdown file per project:

```text
reference_mds/<project_name>.md
```

Example:

```text
reference_mds/spec_nmr_restraints.md
```

If no matching reference exists for a project, the evaluator runs intrinsic
evaluation only.

## Running

Set your API key:

```sh
export OPENAI_API_KEY="..."
```

Open `evaluator_notebook.ipynb` in Jupyter, update the configuration cell, then
run the cells top to bottom.

`evaluator_notebook.py` is also available as a lightweight percent-cell script
for editors that prefer `# %%` notebooks.

The default configuration writes:

```text
results/evaluations.json
results/summary.md
results/plots/<project>_category_bars.svg
results/plots/<project>_score_heatmap.svg
```

The Markdown summary embeds the SVG plots so you can compare individual rubric
categories across models, not just the overall score.

## Testing

Run the dependency-free test suite with:

```sh
python3 -m unittest discover -s tests
```
