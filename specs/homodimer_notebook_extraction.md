# Homodimer Diagnostic Notebook — Structured Extraction

**Source spec:** `specs/homodimer_diagnostic_notebook_spec.md` (2026-03-25)
**Ground truth:** `notebooks/homodimer_diagnostic.ipynb`
**Extracted:** 2026-05-06

---

## 1. Architecture Decisions

### 1.1 Data Fetching Strategy

**Chosen:** Two-tier online/offline mode with AFDB REST API as the default.

- Primary: `GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}` returns metadata, then three separate downloads (mmCIF, PAE JSON, pLDDT JSON) via direct URL fields.
- Fallback: `ipywidgets.FileUpload` widgets allow local file upload (mmCIF required; PAE and pLDDT optional). Controlled by `USE_LOCAL_FILE = False` flag.
- The fallback was **not in the spec** — it was added to support offline/air-gapped execution and PyCharm notebook usage without Colab runtime.

**Why REST over direct download:** The AFDB API response contains all file URLs plus metadata (organism, gene name, model version, sequence) in one call. Parsing that first avoids hardcoding URL patterns.

**Why not BinaryCIF for parsing:** The notebook downloads BinaryCIF (`bcifUrl`) only for MolViewSpec 3D views (performance). Coordinate parsing uses the text mmCIF (`cifUrl`) to avoid a binary decoding dependency.

### 1.2 Score Formula Sources

All five scores reimplemented from scratch in Python/NumPy, sourced from the [DunbrackLab/IPSAE repository](https://github.com/DunbrackLab/IPSAE) (`ipsae.py`, v4, January 2026). The notebook does **not** subprocess or import the script — it inlines the formulas.

**Why inline:** Subprocess calls would create an external file dependency, violating the self-contained notebook requirement. Inlining also allows per-residue decomposition (the script only outputs aggregate values).

### 1.3 Interface Detection Method

**Chosen:** CB-CB Euclidean distance ≤ 8.0 Å (CA for GLY), vectorised via NumPy broadcasting.

```python
diff = coords_A[:, np.newaxis, :] - coords_B[np.newaxis, :, :]
dist_matrix = np.sqrt((diff ** 2).sum(axis=-1))
```

**Known divergence from AFDB production:** The production `interface.py` (Majewski, Apache 2.0) uses CA-CA distances at 8.0 Å via PyTorch `radius_graph` on GPU. The notebook uses CB-CB to match the IPSAE scoring code's contact definition. The spec documents this explicitly; the notebook does not repeat the warning in a markdown cell (minor omission).

**Why no torch-geometric:** Eliminated as a dependency per the non-functional requirement (no GPU-dependent libraries, Colab free tier).

### 1.4 Visualisation Stack

| Layer | Library | Rationale |
|-------|---------|-----------|
| 2D plots (heatmaps, profiles, histograms) | `matplotlib` + `seaborn` | Zero extra install on Colab; full control of layout |
| 3D molecular viewer | `molviewspec` → Mol* | Same viewer as PDBe/RCSB; generates portable MVSJ JSON; pip-installable |
| 3D display method | Base64-encoded `IFrame` | Works in PyCharm notebooks, not just Colab; avoids external URL dependency |
| Statistical style | `seaborn.set_style('white')` | Matches spec §7.2 `whitegrid`/`white` — clean background for publication |

### 1.5 mmCIF Parsing

Manual line-by-line parser with no BioPython dependency. The `parse_mmcif_atoms` function locates the `_atom_site` loop by scanning for the `loop_` keyword followed by `_atom_site.` prefixed lines. Field indices are built dynamically from column headers, making the parser robust to column order changes across AFDB mmCIF versions.

---

## 2. Key Code Patterns

### 2.1 Primitive Scoring Functions

Two functions underpin all PAE-based scores and should be extracted for any future scoring pipeline:

```python
def d0_func(L: int) -> float:
    """Length-dependent TM-score normalisation."""
    if L <= 27:
        return 1.0
    return max(1.0, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)

def ptm_func(pae_vals, d0: float):
    """Vectorised TM-score transform. Works on scalars or NumPy arrays."""
    return 1.0 / (1.0 + (pae_vals / d0) ** 2)
```

Both are pure functions with no side effects. `ptm_func` is broadcast-safe and handles the full PAE matrix or a scalar.

### 2.2 Score Functions Return Dicts with Intermediates

All five `compute_*` functions return a flat dict containing both the final score and all intermediate values. This pattern is key — it avoids recomputation when building visualisations and the summary table.

Example shape:
```python
{
    'score': float,
    'n_contacts': int,
    'mean_plddt': float,
    'x': float,
    'if_A': np.ndarray(bool),  # interface mask, shape (nA,)
    'if_B': np.ndarray(bool),
}
```

The interface masks `if_A` / `if_B` are computed inside the scoring functions and reused downstream for pLDDT analysis and 3D colouring — no recomputation of contacts.

### 2.3 mmCIF Parser Structure

Two-function split:
1. `parse_mmcif_atoms(cif_text) → (records, col_idx)` — tokenises the atom site loop into a list of field lists and a column-name → index mapping.
2. `extract_cb_coords(records, col_idx) → (ch_coords, ch_resids, ch_resnames, ch_plddt)` — applies atom-type filtering (CB preferred, CA fallback for GLY), builds per-chain dicts sorted by residue number.

This split lets you reuse `parse_mmcif_atoms` for any downstream mmCIF task (e.g., extracting HETATM records) without coupling to the CB-selection logic.

### 2.4 PAE Quadrant Extraction Pattern

```python
nA, nB = len_A_pae, len_B_pae
pae_AA = pae_matrix[:nA,      :nA]
pae_AB = pae_matrix[:nA,      nA:nA+nB]
pae_BA = pae_matrix[nA:nA+nB, :nA]
pae_BB = pae_matrix[nA:nA+nB, nA:nA+nB]
```

Chain lengths are derived from the PAE JSON `chains` field, with a fallback to structure-derived lengths. This is the canonical decomposition and should be the template for any PAE-based analysis.

### 2.5 MolViewSpec Builder Pattern

Each 3D view follows the same four-step chain:

```python
builder = mvs.create_builder()
structure = (
    builder
    .download(url=struct_url)
    .parse(format=struct_fmt)      # 'bcif' or 'mmcif'
    .model_structure()
)
(structure
 .component(selector=mvs.ComponentExpression(label_asym_id='A'))
 .representation(type='cartoon')
 .color(color='#009688'))
```

Per-residue colouring loops over residue indices and builds one `ComponentExpression` per residue with `beg_label_seq_id` / `end_label_seq_id`. This is verbose but avoids needing annotation files.

### 2.6 Traffic-Light Classifier

```python
def traffic_light(val, green_thresh, amber_thresh):
    if val >= green_thresh: return 'green', 'HIGH'
    if val >= amber_thresh: return 'amber', 'MODERATE'
    return 'red', 'LOW'
```

Paired with a `thresholds` dict keyed by score name, this pattern generalises cleanly to any set of scores with different scales.

### 2.7 Inline HTML Display in Notebooks

```python
def show_mol_view(state, label, width=950, height=600):
    html = state.molstar_html()
    encoded = base64.b64encode(html.encode()).decode()
    display(HTML(f'<div ...>{label}</div>'))
    display(IFrame(src=f'data:text/html;base64,{encoded}', width=width, height=height))
```

The base64 `data:` URI pattern works in any Jupyter-compatible environment (Colab, JupyterLab, PyCharm) without serving a local HTTP server.

---

## 3. Data Contracts

### 3.1 AFDB Metadata API

**Endpoint:** `GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}`

**Returns:** JSON array; `entries[0]` is consumed.

| Field consumed | Type | Used for |
|----------------|------|----------|
| `cifUrl` | string URL | mmCIF download (coordinate parsing) |
| `bcifUrl` | string URL | BinaryCIF for MolViewSpec (3D viewer) |
| `paeDocUrl` | string URL | PAE matrix JSON download |
| `plddtDocUrl` | string URL | pLDDT JSON download |
| `uniprotAccession` | string | Metadata display |
| `sequence` | string | Monomer length display |
| `organismScientificName` | string | Metadata (falls back to `organism`) |
| `proteinFullName` | string | Metadata (falls back to `uniprotDescription`) |
| `geneNames` | string | Metadata display |
| `modelVersion` | string | Metadata display |

**Accession format:** `AF-{numeric_id}` (e.g., `AF-0000000065889468`).

### 3.2 PAE JSON

```json
[{
  "predicted_aligned_error": [[float, ...], ...],
  "max_predicted_aligned_error": 27.47,
  "chains": [
    {"name": "...", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 172},
    {"name": "...", "label_asym_id": "B", "sequenceStart": 1, "sequenceEnd": 172}
  ]
}]
```

- Matrix is `(nA+nB) × (nA+nB)`, dtype float32 after `np.array()`.
- `chains` field used to derive `nA`, `nB`; sorted by `label_asym_id` for determinism.
- `max_predicted_aligned_error` used as `vmax` for PAE heatmap colour scale (defaults to 31.75 if absent).
- Matrix is asymmetric in values: `PAE[i][j] != PAE[j][i]`.

### 3.3 pLDDT JSON

```json
{
  "residueNumber": [1, 2, ..., 344],
  "confidenceScore": [float, ...],
  "confidenceCategory": ["Low", "Confident", ...],
  "chains": [{"label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 172}, ...]
}
```

- `confidenceScore` sliced as `[:nA]` for chain A, `[nA:nA+nB]` for chain B.
- pLDDT from JSON and from `B_iso_or_equiv` in mmCIF are consistent; the notebook uses JSON values for scoring and mmCIF B-factors for display (both in cells, JSON takes precedence in score functions).

### 3.4 mmCIF Atom Site Fields Required

`group_PDB`, `pdbx_PDB_model_num` (optional), `label_atom_id`, `label_asym_id`, `label_seq_id`, `label_comp_id`, `Cartn_x`, `Cartn_y`, `Cartn_z`, `B_iso_or_equiv`.

- `label_seq_id` values of `'.'` or `'?'` are skipped.
- Only `group_PDB == 'ATOM'` records are consumed (HETATM skipped).
- Only `label_atom_id` in `('CA', 'CB')` are parsed.
- Model filtering: if `pdbx_PDB_model_num` exists, only model `'1'` is kept.

---

## 4. Deviations from Spec

| # | Spec requirement | Actual implementation | Reason / impact |
|---|------------------|-----------------------|-----------------|
| 1 | No local file mode mentioned | `USE_LOCAL_FILE` flag + ipywidgets FileUpload | Enables offline use (PyCharm, air-gapped); adds `ipywidgets` as an implicit dependency |
| 2 | MolViewSpec: save `.mvsj` and provide URL or inline widget | Base64-encoded HTML `IFrame` only | URL approach requires serving infrastructure; `.mvsj` file write not included. The inline display works universally but does not produce a portable file |
| 3 | All 4 MolViewSpec views must be generated | View 4 (disagreement) uses only chain A masks | Chain B disagreement not visualised; functionally partial. Spec listed 3D views 3 and 4 as "stretch goals" |
| 4 | ipSAE: `distance_cutoff ≤ 15` AFDB parameter documented | Not used in compute_ipsae | The 15 Å cutoff appears in the AFDB production ipSAE pipeline (domain-range filtering) but is not part of the per-residue TM-score formula as implemented from the IPSAE repo. No functional gap for the score values |
| 5 | Markdown cell documents CB-CB vs CA-CA distinction explicitly | Not added as a standalone markdown cell | The `parse_mmcif_atoms` docstring omits this; the distinction is only in the spec. Low priority for correctness, matters for user education |
| 6 | Verify against three test accessions | Only one accession (`AF-0000000065889468`) present in notebook | TBD cases not included; the notebook runs a single user-supplied accession |
| 7 | Score intermediates displayed as a table | Printed as formatted text, not an HTML table | Functional equivalent; easier to maintain, harder to scan visually |

---

## 5. Failure Modes and Edge Cases

### 5.1 Explicitly Handled

| Condition | Handling |
|-----------|----------|
| Zero inter-chain contacts | `compute_pdockq` returns `score=0.018` (sigmoid minimum); `compute_pdockq2` returns `score=0.005`. No division-by-zero. |
| Residue with no valid PAE pairs for ipSAE | `_per_res` inner loop: `if n0 == 0: continue`. Score stays at 0.0 for that residue. |
| GLY (no CB atom) | `extract_cb_coords` prefers CB, falls back to CA for any residue. Works for GLY automatically. |
| mmCIF residue with `label_seq_id = '.'` or `'?'` | Skipped via `if res_str in ('.', '?') or not res_str.lstrip('-').isdigit()` |
| Multiple NMR models in mmCIF | Filtered to `pdbx_PDB_model_num == '1'` when the field exists |
| PAE chain info absent | Falls back to `len(ch_resids[chain_ids[0]])` / `[chain_ids[1]]` for nA/nB |
| MolViewSpec not installed | Entire Section 6 wrapped in `try/except`; prints failure message and continues |
| Local file mode with no PAE or pLDDT uploaded | Sets `pae_raw = None` / `plddt_raw = None` with warning — **downstream cells will break** (no null-guard in later cells). Known limitation. |

### 5.2 Known Limitations (Not Handled)

- **Single monomer accessions:** If the AFDB entry has only one chain, the chain-split logic (`nA`, `nB`) will produce a zero-length chain B and scores will be 0.0 without an informative error.
- **Very small interfaces (< 2 contacts):** pDockQ sigmoid behaviour is tested but `log10(1) = 0` gives `x = 0`, producing the minimum score. Not an error, but may be surprising.
- **Symmetric homodimer assumption:** The notebook assumes chains are labelled `A` and `B` and that `chain_ids = sorted(ch_coords.keys())` gives `['A', 'B']`. A dimer with different labelling (e.g., `A` and `C`) will parse correctly but the PAE quadrant slicing may misalign if the sort order doesn't match the PAE JSON chain order.
- **Memory for very large dimers:** The distance matrix is `(nA, nB, 3)` floats in memory. For a 1000-residue monomer this is ~24 MB — fine for Colab. Larger assemblies would need chunking.

---

## 6. Generalisable Prompt Patterns

These patterns from the spec-to-notebook workflow transferred effectively and should be reused.

### 6.1 Spec-First with Exact Formula Transcription

The spec included complete mathematical definitions with variable names matching the reference implementation:

```
d0(L) = max(1.0, 1.24 * (L - 15)^(1/3) - 1.8)   for L > 27
ptm(x, d0) = 1 / (1 + (x / d0)^2)
```

This level of precision eliminated ambiguity and allowed the notebook to match the reference output within ±0.001. **Pattern:** For any numerical reimplementation, include the formula, the variable semantics, and a reference implementation URL in the spec.

### 6.2 Separate "What it measures" from Formula from "Key insight for users"

Each score section in the spec had three sub-sections: definition, formula, and a plain-English key insight. The notebook's markdown cells directly consumed the key insight text. **Pattern:** Write the user-facing explanation in the spec; it flows directly into the notebook's educational markdown cells.

### 6.3 JSON Schema Inclusion in Spec

Spec §4.3 and §4.4 provided exact JSON schemas with field names, types, and nesting. This allowed the data parsing code to be written without exploratory API calls during generation. **Pattern:** For any API-driven notebook, include a concrete JSON sample in the spec, not just a prose description of the response.

### 6.4 Dependency List with Rationale

Spec §8.2 listed exactly: `numpy, matplotlib, seaborn, requests, molviewspec`. Prohibition on BioPython, PyTorch, torch-geometric was explicit. This produced a notebook that runs on a free Colab instance in < 60 seconds. **Pattern:** State allowed and prohibited dependencies separately; prohibited is more valuable than allowed.

### 6.5 Acceptance Criteria with Numerical Tolerances

```
Score accuracy: all recomputed scores must match ipsae.py output. Tolerance: ±0.001
```

This defined a testable contract. **Pattern:** Always include a tolerance when requiring numerical accuracy; "must match" is ambiguous for floating-point outputs.

### 6.6 Heuristics Written as Conditional Prose

The diagnostic logic in spec §6 Section 7 was written as if-then statements in plain English with score variable names embedded. This mapped directly to the `if tl_ipsae == 'green' and tl_pdockq in ('amber', 'red'):` pattern in the notebook. **Pattern:** Write diagnostic logic in pseudo-code prose, not flowcharts; Claude translates prose conditionals into code more reliably than it interprets diagram descriptions.

### 6.7 Known Implementation Discrepancy Documentation

The spec explicitly documented the CB-CB vs CA-CA discrepancy (§5.2 and §9.5) without requiring the notebook to resolve it — just to document it. This avoided the spec-writer needing to make a technical choice they were uncertain about. **Pattern:** "Document and continue" is a valid spec instruction for known discrepancies between reference implementations.

### 6.8 Colour Scheme as a Table

The spec's colour table (§7.1) with hex codes was reproduced verbatim as Python constants at the top of the notebook (`COLOUR_A = '#009688'`, etc.). **Pattern:** Provide hex values, not colour names, in specs for UI work; they transfer directly to code without lookup.

---

## Appendix: Score Reference Card

| Score | Input | PAE cutoff | Distance cutoff | d0 based on | Aggregation |
|-------|-------|-----------|-----------------|-------------|-------------|
| ipTM | Full inter-chain PAE | None | None | `nA + nB` (fixed) | max per-residue mean |
| ipSAE_d0res | Inter-chain PAE | < 10 | None | valid pairs per residue | max per-residue mean |
| ipSAE_d0chn | Inter-chain PAE | < 10 | None | `nA + nB` (fixed) | max per-residue mean |
| ipSAE_d0dom | Inter-chain PAE | < 10 | None | residues with any inter-chain PAE < 10 | max per-residue mean |
| pDockQ | CB-CB distances + pLDDT | None | ≤ 8 Å | — | sigmoid(mean\_pLDDT × log10(n\_contacts)) |
| pDockQ2 | CB-CB + PAE at contacts | Fixed d0=10 | ≤ 8 Å | fixed 10.0 | sigmoid(mean\_pLDDT × mean\_ptm) |
| LIS | Inter-chain PAE | < 12 | None | — | mean((12 − PAE) / 12) |

**Primary AFDB classifier:** `ipSAE_d0res ≥ 0.6` (combined with backbone clashes ≤ 10, pLDDT ≥ 70).
