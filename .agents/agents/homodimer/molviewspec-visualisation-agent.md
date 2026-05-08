# MolViewSpec Visualisation Agent

## Purpose

This is a functional task-executing agent, not an advisory persona. Its sole responsibility is to produce correct MolViewSpec visualisation code for Section 6 of the homodimer diagnostic notebook. It is the only agent that touches Section 6. It does not touch score computation, interface detection, or pLDDT parsing — it consumes those results as inputs.

The agent exists because the MolViewSpec Python API is version-sensitive, the per-residue coloring pattern is non-obvious, and the display method (`molstar_html()`) may not be present in all installed versions. It must verify the API surface before writing any code.

---

## Inputs

The agent requires the following before starting. All are produced by upstream notebook sections.

| Input | Type | Source | Notes |
|---|---|---|---|
| `structure_source` | `'afdb'` or `'pdb'` | Detected from metadata | Determines which views are available |
| `structure_url` | string | `meta['bcifUrl']` or `meta['cifUrl']` (AFDB); RCSB URL (PDB) | bcif preferred for performance |
| `structure_format` | `'bcif'` or `'mmcif'` | Derived from URL | Must match what `.parse(format=...)` accepts |
| `chain_ids` | list of str | Parsed from mmCIF, sorted | e.g. `['A', 'B']` |
| `coloring_mode` | `'plddt'`, `'bfactor'`, or `'neutral'` | Set at notebook top based on source | Controls which views are generated |
| `plddt_per_chain` | dict of `{chain_id: np.ndarray}` | Section 1 | Required for `'plddt'` mode |
| `bfactor_per_chain` | dict of `{chain_id: np.ndarray}` | Section 1 | Required for `'bfactor'` mode |
| `ipsae_d0res_per_chain` | dict of `{chain_id: np.ndarray}` | Section 4 | Required for View 3 (AFDB only) |
| `interface_mask_per_chain` | dict of `{chain_id: np.ndarray[bool]}` | Section 2 | Required for Views 1 and 3 |
| `resids_per_chain` | dict of `{chain_id: list[int]}` | Section 1 | `label_seq_id` values for `ComponentExpression` |

---

## API Verification Steps

The agent performs these checks via `mcp__plugin_context7_context7__resolve-library-id` and `mcp__plugin_context7_context7__query-docs` **before writing any view code**.

1. **Resolve version** — confirm the installed version of `molviewspec` and identify the matching documentation version.

2. **Verify builder API surface** — confirm that these calls exist and have the expected signatures:
   - `mvs.create_builder()`
   - `.download(url=...).parse(format=...).model_structure()`
   - `.component(selector=mvs.ComponentExpression(...))`
   - `.representation(type=...)`
   - `.color(color=...)`
   - `builder.get_state().molstar_html()`

3. **Verify `ComponentExpression` field names** — confirm that `label_asym_id`, `beg_label_seq_id`, and `end_label_seq_id` are the correct fields for per-residue selection. These names have changed between versions.

4. **Check annotation-file support** — determine whether the installed version supports annotation-file-based coloring (a single JSON payload covering all residues) as an alternative to per-residue `ComponentExpression` loops. Annotation files are preferred if available.

5. **Verify `parse(format=...)` strings** — confirm accepted values for both `'bcif'` and `'mmcif'`. Some versions use `'cif'` rather than `'mmcif'`.

---

## Views Produced

| View | Description | AFDB (`plddt`) | PDB crystal (`bfactor`) | PDB NMR / `neutral` |
|---|---|---|---|---|
| 1 | Chain overview — teal (A), coral (B), amber (interface) | ✓ | ✓ | ✓ |
| 2 | Confidence mapping | pLDDT 4-band discrete | B-factor continuous, percentile-normalised | skip |
| 3 | Interface residues coloured by ipSAE d0res score | ✓ | skip | skip |
| 4 | PAE/contact disagreement regions | ✓ | skip | skip |

For skipped views, the agent emits `print('View N skipped: requires PAE / not available for this source.')` rather than raising an error.

### View 1 — Chain Overview
Cartoon representation. Chain A: `#009688` (teal). Chain B: `#FF7043` (coral). Interface residues overlaid as ball-and-stick in `#FFC107` (amber). Always generated.

### View 2 — Confidence Mapping

**pLDDT mode (AFDB):** Four discrete bands applied per residue.

| pLDDT range | Colour |
|---|---|
| > 90 | `#1565C0` (dark blue) |
| 70–90 | `#42A5F5` (light blue) |
| 50–70 | `#FFCA28` (yellow) |
| < 50 | `#EF6C00` (orange) |

**B-factor mode (PDB crystal):** Continuous colormap (`viridis` or `coolwarm`). Normalise using the 5th–95th percentile of B-factor values across all ATOM records (HETATM excluded). Low B-factor = ordered = cool end of colormap. Invert relative to pLDDT semantics.

```python
b_min, b_max = np.percentile(all_bfactors_atom_only, [5, 95])
norm = (b_val - b_min) / max(b_max - b_min, 1e-6)
norm = np.clip(norm, 0, 1)
```

**Neutral mode:** Skip.

### View 3 — Interface ipSAE Gradient (AFDB only)
Non-interface residues: cartoon in `#BDBDBD` (grey). Interface residues: ball-and-stick coloured by per-residue `ipsae_d0res` score using `RdYlGn` colormap (0=red, 1=green).

### View 4 — PAE/Contact Disagreement (AFDB only)
Three classes for chain A interface residues, mapped to ball-and-stick:

| Condition | Colour |
|---|---|
| ipSAE d0res > 0.5 **and** in contact | `#4CAF50` (green) — both signals agree |
| ipSAE d0res > 0.5, **not** in contact | `#2196F3` (blue) — PAE confident, no contact |
| ipSAE d0res ≤ 0.5, in contact | `#F44336` (red) — contact exists, low PAE confidence |

---

## Per-Residue Coloring Strategy

Two approaches, ranked by preference:

**1. Annotation file** (use if supported by installed version) — express all per-residue colors as a single JSON annotation attached to the structure node. No per-residue loop. One `mcp__plugin_context7_context7__query-docs` call confirms availability: query `"annotation file per-residue color ComponentAnnotation custom color scheme"`.

**2. Per-residue `ComponentExpression` loop** (fallback) — one `.component()` call per residue. Include a size guard:

- If total residues across all chains ≤ 500: loop per residue.
- If total residues > 500: collapse to band-level coloring — group residues into score bands (e.g. quartiles), one `ComponentExpression` per band using `beg_label_seq_id`/`end_label_seq_id` spanning contiguous runs. This keeps the MolViewSpec state tree manageable.

---

## Display Function

The agent always emits a `show_mol_view()` function. Before writing it, the agent confirms whether `molstar_html()` exists in the installed version.

- **Primary:** `state.molstar_html()` → base64-encode → `IFrame` with `data:text/html;base64,...` URI. Works inline in Colab and JupyterLab without network access.
- **Fallback (if `molstar_html()` absent):** Serialise the state to `.mvsj` JSON and construct a `https://molstar.org/viewer/?...` URL for the user to open.

```python
def show_mol_view(state, label, width=950, height=600):
    try:
        html = state.molstar_html()
        encoded = base64.b64encode(html.encode()).decode()
        display(HTML(f'<div style="font-weight:bold;margin:8px 0 4px">{label}</div>'))
        display(IFrame(src=f'data:text/html;base64,{encoded}', width=width, height=height))
    except AttributeError:
        # molstar_html() not available in this version — fall back to viewer URL
        import json
        mvsj = json.dumps(state.as_dict())
        encoded = base64.b64encode(mvsj.encode()).decode()
        url = f'https://molstar.org/viewer/?mvs-data={encoded}&mvs-data-format=mvsj'
        display(HTML(f'<b>{label}</b><br><a href="{url}" target="_blank">Open in Mol* viewer</a>'))
```

---

## Output Structure

The agent produces exactly two notebook cells for Section 6:

**Cell 6a — Markdown explanation cell:** Plain-language description of MolViewSpec, what views are generated, and how to interpret each. Written for a non-structural-biology audience (see spec §3).

**Cell 6b — Code cell:** A single self-contained Python cell with:
- `try/except ImportError` wrapping all molviewspec imports
- `show_mol_view()` definition
- `detect_structure_source()` helper
- `COLORING_MODE` branch logic
- Each view in its own `try/except Exception` block so one failure does not prevent others from rendering
- A closing `print(f'Views generated: {generated_views}')` summary line

```python
try:
    import molviewspec as mvs
    MVS_AVAILABLE = True
except ImportError:
    MVS_AVAILABLE = False
    print('molviewspec not installed — skipping Section 6.')

if MVS_AVAILABLE:
    # show_mol_view(), detect_structure_source(), view generation ...
    generated_views = []

    try:
        # View 1: chain overview
        ...
        generated_views.append('View 1 (chain overview)')
    except Exception as e:
        print(f'View 1 failed: {e}')

    # ... views 2–4 with same pattern ...

    print(f'Views generated: {", ".join(generated_views) or "none"}')
```

---

## What This Agent Does Not Touch

- **Section 1** — mmCIF parsing, AFDB API fetch, pLDDT/B-factor array construction
- **Section 2** — interface detection, contact mask, distance matrix
- **Section 4** — score computation, per-residue ipSAE arrays

The agent consumes `interface_mask_per_chain`, `ipsae_d0res_per_chain`, `plddt_per_chain`, and `resids_per_chain` as read-only inputs. If any are missing, it downgrades gracefully to the nearest available view rather than raising.

---

## Invocation

When invoking this agent in a new session, provide the inputs table above and specify the structure source and coloring mode explicitly.

Example:
```
MolViewSpec Visualisation Agent — produce Section 6 code for:
  structure_source: 'afdb'
  structure_url: meta['bcifUrl']  (bcif format)
  chain_ids: ['A', 'B']
  coloring_mode: 'plddt'
  plddt_per_chain: {plddt_A, plddt_B} from Section 1
  ipsae_d0res_per_chain: {per_res_AB_d0res, per_res_BA_d0res} from Section 4
  interface_mask_per_chain: {interface_A, interface_B} from Section 2
  resids_per_chain: {ch_resids['A'], ch_resids['B']} from Section 1

Verify the MolViewSpec API against the installed version before writing any code.
All four views are required.
```

For a PDB crystal structure:
```
MolViewSpec Visualisation Agent — produce Section 6 code for:
  structure_source: 'pdb'
  structure_url: 'https://files.rcsb.org/download/1ABC.cif'
  structure_format: 'mmcif'
  chain_ids: ['A', 'B']
  coloring_mode: 'bfactor'
  bfactor_per_chain: {bfactor_A, bfactor_B} from Section 1
  interface_mask_per_chain: {interface_A, interface_B} from Section 2
  resids_per_chain: {ch_resids['A'], ch_resids['B']} from Section 1

Views 3 and 4 are not available (no PAE). View 2 uses B-factor coloring.
```
