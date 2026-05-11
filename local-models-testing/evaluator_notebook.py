# %%
"""Notebook-style runner for evaluating Markdown project specs.

Open this file in VS Code, Jupyter, or any editor that supports `# %%` cells.
Edit the configuration cell, then run the cells top to bottom.
"""

# %%
from pathlib import Path

from IPython.display import SVG, display

from spec_evaluator import evaluate_directory


# %%
# Configuration
INPUT_DIR = Path("candidate_mds")
REFERENCE_DIR = Path("reference_mds")
OUTPUT_DIR = Path("results")

JUDGE_MODEL = "gpt-5.5"
FILENAME_SEPARATOR = "--"

# Set to None to run without references.
ACTIVE_REFERENCE_DIR = REFERENCE_DIR


# %%
# Run evaluation.
#
# Required candidate filename format:
#   <project_name>--<model_name>.md
#
# Example:
#   spec_nmr_restraints--gemma_31b.md
#
# The OpenAI judge uses OPENAI_API_KEY from the environment.
result = evaluate_directory(
    input_dir=INPUT_DIR,
    output_dir=OUTPUT_DIR,
    reference_dir=ACTIVE_REFERENCE_DIR,
    model=JUDGE_MODEL,
    separator=FILENAME_SEPARATOR,
)


# %%
# Inspect project rankings in the notebook.
result["rankings"]


# %%
# Display generated category comparison plots.
for project, project_plots in result.get("plots", {}).items():
    print(project)
    for plot_path in project_plots.values():
        display(SVG(filename=str(OUTPUT_DIR / plot_path)))


# %%
# Inspect full normalized evaluations.
result["evaluations"]
