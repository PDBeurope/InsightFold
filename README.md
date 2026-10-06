# InsightFold

Source code for the prototype notebooks published on the **[PDBe Workbench](https://www.ebi.ac.uk/pdbe/workbench/)**.

The Workbench is a public dashboard for exploring computational notebooks for structural analysis. You can browse and filter the notebooks there and open any of them in Google Colab, with nothing to install. This repository holds the notebooks themselves and the libraries behind them, for anyone who wants to read the implementation, reuse it, or contribute.

These are **prototypes**: exploratory analyses published so the community can use them and build on them, not production services.

## The notebooks

| Notebook | The question it answers | |
|---|---|---|
| **Dimer Confidence Metric Diagnostic** | How confident is an AlphaFold prediction about a dimer's *interface*, as opposed to the folds either side of it? Computes ipTM, three ipSAE variants, pDockQ, pDockQ2 and LIS, with published thresholds and the provenance of each. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/PDBeurope/InsightFold/blob/main/notebooks/homodimer_diagnostic.ipynb) |
| **Cluster Quality Diagnostic** | How is pLDDT distributed across my protein's AFDB clusters, and is a relative better predicted than my protein? Includes TM-align superposition against the best and worst-predicted members. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/PDBeurope/InsightFold/blob/main/notebooks/cluster_quality_diagnostic.ipynb) |
| **Interface Interaction Analysis** | When several deposited assemblies map to the same PDBe-KB complex, how similar are their protein-protein interfaces? | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/PDBeurope/InsightFold/blob/main/notebooks/interface_analysis.ipynb) |

## Running them

**Through the Workbench** is the easiest route: filter for the analysis you want and open it, and it launches in Colab against the current version here.

**Straight from a badge above** does the same thing. Set the input in the form at the top of the notebook, then `Runtime` → `Run all`. The first cell fetches the code; nothing else needs editing.

**Locally**, clone and open the notebook from inside the checkout. The same first cell finds your checkout and uses it, with no clone and no network fetch of the code.

```
git clone https://github.com/PDBeurope/InsightFold.git
cd InsightFold
jupyter lab notebooks/homodimer_diagnostic.ipynb
```

The notebooks are committed without outputs, so a fresh clone shows the code and the prose, and you produce the figures by running it.

## How the code is organised

Each notebook is deliberately thin. The cells are a readable sequence of single calls; the analysis lives in an importable library, so it can be reused and tested without running a notebook.

```
notebooks/     the three notebooks above
src/insightfold/
  complex_interface_utils.py    dimer confidence metrics, figures, 3D views
  cluster_quality_utils.py      AFDB cluster analysis, TM-align, figures, 3D views
  notebook_setup.py             the shared bootstrap every notebook's first cell calls
specs/         requirements, data contracts, pinned fixtures and expected values
prd/           what each notebook is for, and what it deliberately does not do
CLAUDE.md      conventions, API contracts and the traps worth knowing before editing
```

## Contributing

We are working out how to open this up, most likely by inviting individual contributors to propose changes to the notebooks through pull requests. That process is still being developed, so there is nothing to join yet. Please check back as the Workbench grows.
