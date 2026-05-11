# Notebook Execution Validation Report

| Field | Value |
|---|---|
| Notebook | `notebooks/protein_model_chem.ipynb` |
| Spec Pack | `specs/protein_model_chem/` |
| Validation Level | full |
| Validator | Codex / `$notebook-execution-validation` |
| Date | 2026-05-08 |
| Verdict | fail |

## Execution Evidence

| Item | Value |
|---|---|
| Execution Tool | `nbclient` |
| Command / Method | Restarted a local Jupyter kernel and executed all cells with `NotebookClient(...).execute()` using `.venv/bin/python` |
| Kernel / Python Version | Python 3.11.9, kernel `python3` |
| Runtime | 21.462 seconds |
| Fixtures Run | FX-001 default fixture (`2LZM`, `A:L99A`) |

Package versions in the executed environment:

- `numpy 2.3.3`
- `pandas 2.3.3`
- `requests 2.32.5`
- `matplotlib 3.10.7`
- `networkx 3.6.1`
- `scipy 1.17.1`
- `biopython 1.87`
- `molviewspec 1.8.1`
- `plotly 6.7.0`
- `nbclient 0.10.4`
- `nbformat 5.10.4`
- `ipykernel 7.2.0`

## Check Results

| Check | Status | Evidence | Notes |
|---|---|---|---|
| Restart and run all | pass | All 25 cells executed successfully; 12 code cells have execution counts `1..12`. | Reproducible after dependency fix and `.venv` sync. |
| Fixture agreement | fail | Final validation cell reports `FX-001 local interactions present = False`. | Spec expects at least one hydrophobic or packing contact near residue 99. |
| Data contracts | pass | Notebook produced structured runtime artifacts: `residue_rows=164`, `interaction_rows=123`, `delta_rows=1`, `summary_rows=3`, `mvs_states=1`. | Negative fixtures were not exercised in this run. |
| Dependency policy | pass | Project kernel now imports `networkx`, `scipy`, `Bio`, `molviewspec`, and `plotly`; the forbidden legacy 3D viewer package is absent. | `networkx` was added to `pyproject.toml`; lockfile and `.venv` were refreshed. |
| Required outputs | pass | Mutation delta table, mechanistic summary table, notebook validation summary, and runtime summary all executed. | Export writing remained disabled in this fixture run (`write_exports=False`). |
| Visualization rendering | pass | Visualization cells executed without error and produced one MolViewSpec state (`mvs_states=1`). | This run confirms runtime generation, not human visual quality review. |
| Documentation presence | pass | First cell states notebook purpose, default fixture, and limitations; markdown sections still precede major code blocks. | Qualitative scientific review is still separate. |

## Failures

| ID | Type | Cell / Section | Error or Mismatch | Required Fix |
|---|---|---|---|---|
| F-001 | fixture-mismatch | Final validation cell / FX-001 acceptance check | Internal validation reports `FX-001 local interactions present = False`, while the fixture manifest expects at least one local hydrophobic or packing contact near `A:99`. | Review the mutation-neighborhood interaction filtering and fixture expectation together. Either the interaction detector is missing a local contact or the fixture expectation needs to be revised based on the implemented heuristics. |

## Static Checks Completed

- Notebook JSON parses successfully with `python -m json.tool`.
- All 12 code cells parse with `ast.parse`; no syntax issues found.
- No unresolved `TODO`, `FIXME`, absolute local workspace paths, obvious secret strings, or forbidden legacy 3D viewer references were found in `notebooks/protein_model_chem.ipynb`, `specs/protein_model_chem/`, or `pyproject.toml`.
- `git diff --check` passed for `pyproject.toml` and `uv.lock`.

## Rerun History

- Initial validation attempt failed before execution because the active environment was missing `networkx`, `scipy`, and `biopython`.
- The blocker was fixed by adding `networkx` to `pyproject.toml`, refreshing `uv.lock`, and syncing `.venv`.
- The rerun completed successfully at the kernel level, which closes the reproducibility blocker but exposes the remaining FX-001 fixture mismatch above.

## Limitations

- Only FX-001 was executed in this run. FX-002, FX-005, FX-006, and the optional FX-007 fixture still need separate validation runs.
- The notebook emitted an `nbformat` `MissingIDFieldWarning`; this did not block execution, but notebook cell IDs should eventually be normalized.
- This report validates runtime behavior and required artifacts, not domain correctness of thresholds or mechanistic claims.

## Recommended Next Action

- [ ] Investigate the FX-001 local-interaction mismatch and decide whether to change the detector or the fixture expectation.
- [ ] Rerun `$notebook-execution-validation` after the FX-001 mismatch is resolved.
- [ ] Run FX-002, FX-005, and FX-006 validation passes.
- [ ] Send to `$notebook-review` only after fixture validation is clean.
