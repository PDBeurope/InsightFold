# InsightFold — CLAUDE.md

## Project Overview

InsightFold analyses AlphaFold DB structure predictions for protein–protein interactions.
The primary use case is homodimer interface quality assessment using AFDB REST APIs and
inline NumPy implementations of published scoring functions.

## Directory Layout

```
src/insightfold/
  interface.py        ← CB/CA-distance interface detection (no torch dependency)
  variants/           ← variant mapping submodule
notebooks/
  homodimer_diagnostic.ipynb    ← reference implementation
  analysis_template.ipynb       ← parametrised scaffold for new analysis types
specs/                ← extraction documents and notebook specs
agents/llm_advisory_roles       ← LLM advisory role definitions
```

## Dependencies

Allowed: `numpy`, `matplotlib`, `seaborn`, `requests`, `molviewspec`, `ipywidgets`.
Prohibited: `biopython`, `torch`, `torch-geometric` (GPU dependency, breaks Colab free tier).
The notebook must run on Colab free tier in < 60 s of wall-clock install time.

---

## Skill: Homodimer Confidence Scoring

Use this skill whenever computing interface quality scores from an AlphaFold DB homodimer
prediction. All formulas are sourced from DunbrackLab/IPSAE (`ipsae.py` v4, Jan 2026).
Tolerance for numerical agreement with the reference: ±0.001.

### Primitive Functions (inline in every scoring cell)

```python
def d0_func(L: int) -> float:
    """TM-score length-dependent normalisation. L = total residues in scoring domain."""
    if L <= 27:
        return 1.0
    return max(1.0, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)

def ptm_func(pae_vals, d0: float):
    """Vectorised TM-score transform. Broadcast-safe (scalar or ndarray)."""
    return 1.0 / (1.0 + (pae_vals / d0) ** 2)
```

### PAE Quadrant Extraction

Always extract quadrants this way — chain lengths from the PAE JSON `chains` field:

```python
nA = pae_json[0]['chains'][0]['sequenceEnd'] - pae_json[0]['chains'][0]['sequenceStart'] + 1
nB = pae_json[0]['chains'][1]['sequenceEnd'] - pae_json[0]['chains'][1]['sequenceStart'] + 1
pae_matrix = np.array(pae_json[0]['predicted_aligned_error'], dtype=np.float32)
pae_AB = pae_matrix[:nA, nA:nA+nB]   # rows=A, cols=B
pae_BA = pae_matrix[nA:nA+nB, :nA]   # rows=B, cols=A
pae_AA = pae_matrix[:nA, :nA]
pae_BB = pae_matrix[nA:nA+nB, nA:nA+nB]
```

### Score Formulas

All `compute_*` functions must return a flat dict with both the final score and
all intermediate values (masks, counts) so callers can reuse them without recomputation.

#### ipTM

```python
def compute_iptm(pae_AB, pae_BA, nA, nB):
    L = nA + nB
    d0 = d0_func(L)
    score_A = ptm_func(pae_AB, d0).mean(axis=1).max()  # per-row max over mean
    score_B = ptm_func(pae_BA, d0).mean(axis=1).max()
    score = max(score_A, score_B)
    return {'score': float(score), 'd0': d0}
```

#### ipSAE (three d0 variants)

The key difference is **what L is passed to `d0_func`**:

| Variant | d0 based on | L passed to d0_func |
|---------|-------------|---------------------|
| `d0res` | per-residue valid pair count | `n_valid` for that residue |
| `d0chn` | full chain pair total | `nA + nB` (fixed) |
| `d0dom` | residues with any inter-chain PAE < 10 | count of such residues |

PAE cutoff for all three: only pairs where `pae < 10` contribute.

```python
def compute_ipsae(pae_AB, pae_BA, nA, nB):
    PAE_CUTOFF = 10.0

    def _per_res(pae_block, n_self, n_other, d0_mode, n_dom=None):
        scores = []
        for i in range(n_self):
            row = pae_block[i]
            mask = row < PAE_CUTOFF
            n0 = int(mask.sum())
            if n0 == 0:
                scores.append(0.0)
                continue
            if d0_mode == 'res':
                d0 = d0_func(n0)
            elif d0_mode == 'chn':
                d0 = d0_func(n_self + n_other)
            else:  # dom
                d0 = d0_func(n_dom)
            scores.append(float(ptm_func(row[mask], d0).mean()))
        return np.array(scores)

    # domain size = residues with at least one inter-chain PAE < cutoff
    n_dom_A = int((pae_AB < PAE_CUTOFF).any(axis=1).sum())
    n_dom_B = int((pae_BA < PAE_CUTOFF).any(axis=1).sum())
    n_dom = max(n_dom_A + n_dom_B, 1)

    for mode in ('res', 'chn', 'dom'):
        sA = _per_res(pae_AB, nA, nB, mode, n_dom)
        sB = _per_res(pae_BA, nB, nA, mode, n_dom)
        yield mode, float(max(sA.max(), sB.max()))
```

#### pDockQ

Contact definition: CB–CB ≤ 8 Å (CA fallback for GLY — see interface.py).

```python
def compute_pdockq(if_A, if_B, plddt_arr_A, plddt_arr_B):
    n = int(if_A.sum()) + int(if_B.sum())
    if n == 0:
        return {'score': 0.018, 'n_contacts': 0}
    mean_plddt = np.concatenate([plddt_arr_A[if_A], plddt_arr_B[if_B]]).mean()
    x = mean_plddt * np.log10(max(n, 1))
    score = 0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018
    return {'score': float(score), 'n_contacts': n, 'mean_plddt': float(mean_plddt), 'x': float(x),
            'if_A': if_A, 'if_B': if_B}
```

#### pDockQ2

Uses PAE at contact pairs (fixed d0 = 10.0).

```python
def compute_pdockq2(coords_A, coords_B, pae_AB, plddt_A, plddt_B,
                    dist_cutoff=8.0, d0_fixed=10.0):
    diff = coords_A[:, np.newaxis, :] - coords_B[np.newaxis, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    in_contact = dist <= dist_cutoff
    if not in_contact.any():
        return {'score': 0.005, 'n_contacts': 0}
    pae_contacts = pae_AB[in_contact]
    plddt_contacts = np.concatenate([
        np.repeat(plddt_A, in_contact.sum(axis=1)),
        np.tile(plddt_B, (len(plddt_A), 1))[in_contact]
    ])
    mean_plddt = plddt_contacts.mean()
    mean_ptm = ptm_func(pae_contacts, d0_fixed).mean()
    x = mean_plddt * mean_ptm
    score = 0.715 / (1 + np.exp(-12.3 * (x - 0.605))) + 0.005
    return {'score': float(score), 'n_contacts': int(in_contact.sum()),
            'mean_plddt': float(mean_plddt), 'mean_ptm': float(mean_ptm), 'x': float(x)}
```

#### LIS (Local Interaction Score)

```python
def compute_lis(pae_AB, pae_BA):
    PAE_CUTOFF = 12.0
    def _lis_half(block):
        mask = block < PAE_CUTOFF
        if not mask.any():
            return 0.0
        return float(((PAE_CUTOFF - block[mask]) / PAE_CUTOFF).mean())
    score = max(_lis_half(pae_AB), _lis_half(pae_BA))
    return {'score': score}
```

### Traffic-Light Thresholds

```python
THRESHOLDS = {
    'ipsae_d0res': (0.6, 0.4),   # (green_min, amber_min)
    'iptm':        (0.7, 0.5),
    'pdockq':      (0.23, 0.09),
    'pdockq2':     (0.5, 0.23),
    'lis':         (0.15, 0.09),
}

def traffic_light(val, score_name):
    g, a = THRESHOLDS[score_name]
    if val >= g: return 'green',  'HIGH'
    if val >= a: return 'amber',  'MODERATE'
    return        'red',   'LOW'
```

Primary AFDB classifier: `ipSAE_d0res ≥ 0.6` (combined with backbone clashes ≤ 10, pLDDT ≥ 70).

---

## AFDB REST API Contracts

### Metadata Endpoint

```
GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}
```

Returns a JSON **array**; consume `result[0]`. Accession format: `AF-{id}` (e.g. `AF-0000000065889468`).

Fields consumed:

| Field | Type | Used for |
|-------|------|----------|
| `cifUrl` | URL string | mmCIF download (coordinate parsing) |
| `bcifUrl` | URL string | BinaryCIF for MolViewSpec (3D viewer only) |
| `paeDocUrl` | URL string | PAE JSON download |
| `plddtDocUrl` | URL string | pLDDT JSON download |
| `uniprotAccession` | string | Display |
| `sequence` | string | Monomer length |
| `organismScientificName` | string | Display (fallback: `organism`) |
| `proteinFullName` | string | Display (fallback: `uniprotDescription`) |
| `geneNames` | string | Display |
| `modelVersion` | string | Display |

### PAE JSON Schema

```json
[{
  "predicted_aligned_error": [[float, ...]],
  "max_predicted_aligned_error": 27.47,
  "chains": [
    {"label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 172},
    {"label_asym_id": "B", "sequenceStart": 1, "sequenceEnd": 172}
  ]
}]
```

- Matrix is `(nA+nB) × (nA+nB)`, asymmetric (`PAE[i][j] ≠ PAE[j][i]`).
- Sort `chains` by `label_asym_id` for determinism.
- Use `max_predicted_aligned_error` as `vmax` for heatmap (default 31.75 if absent).

### pLDDT JSON Schema

```json
{
  "residueNumber": [1, 2, ...],
  "confidenceScore": [float, ...],
  "confidenceCategory": ["Low", "Confident", ...],
  "chains": [{"label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 172}, ...]
}
```

Slice: chain A = `[:nA]`, chain B = `[nA:nA+nB]`.

### mmCIF Atom Site Required Fields

`group_PDB`, `label_atom_id`, `label_asym_id`, `label_seq_id`, `label_comp_id`,
`Cartn_x`, `Cartn_y`, `Cartn_z`, `B_iso_or_equiv`, `pdbx_PDB_model_num` (optional).

Filtering rules:
- Only `group_PDB == 'ATOM'` (skip HETATM).
- Only `label_atom_id in ('CA', 'CB')`.
- Skip `label_seq_id` of `'.'` or `'?'`.
- If `pdbx_PDB_model_num` present, keep only model `'1'`.

---

## Interface Detection Contract

See `src/insightfold/interface.py` for the canonical implementation.

Key semantics:
- CB–CB Euclidean distance ≤ 8.0 Å defines a contact.
- GLY has no CB; use CA instead (handled transparently by `extract_cb_coords`).
- **This uses CB–CB, not CA–CA.** AFDB production `interface.py` uses CA–CA with
  PyTorch `radius_graph`. The CB–CB choice matches the IPSAE scoring code's contact
  definition. Do not "fix" this without updating the scoring functions too.
- Distance matrix is `(nA, nB, 3)` in memory; fine up to ~1000-residue monomers (~24 MB).

---

## MolViewSpec Viewer Pattern

```python
import molviewspec as mvs, base64
from IPython.display import display, HTML, IFrame

def show_mol_view(state, label, width=950, height=600):
    html = state.molstar_html()
    encoded = base64.b64encode(html.encode()).decode()
    display(HTML(f'<h4>{label}</h4>'))
    display(IFrame(src=f'data:text/html;base64,{encoded}', width=width, height=height))

builder = mvs.create_builder()
structure = (
    builder
    .download(url=bcif_url)
    .parse(format='bcif')
    .model_structure()
)
structure.component(
    selector=mvs.ComponentExpression(label_asym_id='A')
).representation(type='cartoon').color(color='#009688')
show_mol_view(builder.get_state(), 'Chain A')
```

Per-residue colouring: build one `ComponentExpression` per residue with
`beg_label_seq_id` / `end_label_seq_id`. Wrap the entire Section 6 in
`try/except ImportError` so the notebook degrades gracefully if `molviewspec`
is not installed.

---

## Known Edge Cases

| Condition | Handling |
|-----------|----------|
| Zero inter-chain contacts | pDockQ → 0.018, pDockQ2 → 0.005 (sigmoid minimum) |
| Residue with no valid PAE pairs for ipSAE | Skip that residue (score stays 0.0) |
| GLY (no CB atom) | `extract_cb_coords` prefers CB, falls back to CA |
| `label_seq_id` = `'.'` or `'?'` | Skip row |
| Multiple NMR models | Filter to `pdbx_PDB_model_num == '1'` |
| PAE `chains` field absent | Fall back to structure-derived chain lengths |
| Single monomer accession | Chain-split produces zero-length chain B; scores will be 0.0 (no informative error) |
| Asymmetric chain labels (e.g., A+C) | PAE quadrant slicing may misalign — sort `chain_ids` and verify against PAE `chains` field |
