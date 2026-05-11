---
name: molviewspec-rendering
description: Generate and review MolViewSpec/Mol* 3D molecular visualization code for InsightFold notebooks. Use when Codex needs to build notebook sections with chain coloring, pLDDT or B-factor overlays, interface highlighting, ipSAE/PAE/contact diagnostic views, .mvsj state generation, inline Mol* rendering, or Mol* viewer URL fallbacks.
---

# MolViewSpec Rendering

## Purpose

Produce correct, reproducible MolViewSpec visualizations for InsightFold notebooks.

MolViewSpec/Mol* is the required InsightFold path for 3D protein structure viewers. Do not use imperative notebook 3D viewers as substitutes. If an existing notebook or spec uses another 3D protein viewer, treat that as legacy code to migrate or a policy violation to flag.

Use this skill for notebook sections that render protein structures with MolViewSpec and Mol*, especially when the notebook needs:

- chain overview views
- pLDDT or B-factor coloring
- interface residue highlighting
- ipSAE or score-driven residue coloring
- PAE/contact disagreement views
- `.mvsj` state export
- inline notebook rendering or Mol* viewer links

MolViewSpec code is optional visualization code. It must not be required for upstream data loading, scoring, validation, or interpretation cells to run.

Load `references/molviewspec-policy.md` when reviewing or migrating existing viewer code. Use `assets/molviewspec-section-template.md` when creating a new notebook visualization section.

## Required Inputs

Confirm these inputs before writing view code:

| Input | Type | Notes |
|---|---|---|
| `structure_source` | `'afdb'`, `'pdb'`, or `'local'` | Determines available views |
| `structure_url` | string or `None` | Prefer `bcifUrl` for AFDB, otherwise use mmCIF |
| `structure_format` | `'bcif'`, `'mmcif'`, or version-supported equivalent | Must match installed MolViewSpec API |
| `chain_ids` | list of str | Use parsed chain labels, do not assume `A` and `B` unless contracted |
| `resids_per_chain` | dict[str, list[int]] | Residue IDs for `ComponentExpression` |
| `interface_mask_per_chain` | dict[str, bool array] | Required for interface views |
| `plddt_per_chain` | dict[str, numeric array] | Required for AFDB pLDDT coloring |
| `bfactor_per_chain` | dict[str, numeric array] | Required for PDB B-factor coloring |
| `ipsae_d0res_per_chain` | dict[str, numeric array] | Required for AFDB ipSAE interface coloring |

If a required input is absent, skip only the affected view and emit a clear message. Do not raise unless the spec says that view is mandatory.

## API Verification

MolViewSpec is version-sensitive. Before writing or changing code, inspect the installed API and, when needed, verify against current documentation.

Minimum checks:

```python
import molviewspec as mvs

print(getattr(mvs, "__version__", "unknown"))
print(hasattr(mvs, "create_builder"))
print(hasattr(mvs, "ComponentExpression"))

builder = mvs.create_builder()
state = builder.get_state()
print(hasattr(state, "molstar_html"))
```

Confirm these calls and names for the installed version:

- `mvs.create_builder()`
- `.download(url=...).parse(format=...).model_structure()`
- `.component(selector=mvs.ComponentExpression(...))`
- `.representation(type='cartoon')`
- `.color(color='#RRGGBB')`
- `builder.get_state().molstar_html()`
- `mvs.ComponentExpression(label_asym_id=..., beg_label_seq_id=..., end_label_seq_id=...)`

Also confirm accepted parse format strings. Some versions may prefer `cif` over `mmcif`.

## View Set

### View 1: Chain Overview

Generate whenever structure data is available.

- Cartoon representation.
- Chain A default color: `#009688`.
- Chain B default color: `#FF7043`.
- Additional chains: use distinct, muted colors and label them in the notebook output.
- Interface residues: overlay as ball-and-stick in `#FFC107`.

Do not assume exactly two chains unless the data contract requires it.

### View 2: Confidence Mapping

Generate for AFDB structures with pLDDT or PDB crystal structures with B-factors.

AFDB pLDDT bands:

| Range | Color |
|---|---|
| `>= 90` | `#1565C0` |
| `70 <= x < 90` | `#42A5F5` |
| `50 <= x < 70` | `#FFCA28` |
| `< 50` | `#EF6C00` |

PDB B-factor mode:

- Use ATOM records only.
- Normalize B-factors by 5th and 95th percentiles.
- Treat low B-factor as more ordered.
- State explicitly that B-factor is not pLDDT.

```python
b_min, b_max = np.percentile(all_bfactors_atom_only, [5, 95])
norm = (b_val - b_min) / max(b_max - b_min, 1e-6)
norm = np.clip(norm, 0, 1)
```

Skip confidence mapping for NMR or neutral structures unless the spec defines a substitute confidence metric.

### View 3: Interface Score Coloring

Generate for AFDB structures when per-residue ipSAE or another contracted score is available.

- Non-interface residues: neutral grey `#BDBDBD`.
- Interface residues: ball-and-stick colored by the score.
- Use a red-yellow-green scale unless the spec defines another scale.
- Include the scale meaning in adjacent markdown.

### View 4: PAE/Contact Disagreement

Generate only when both PAE-derived confidence and contact/interface masks are available.

Default classes:

| Condition | Color | Meaning |
|---|---|---|
| score high and in contact | `#4CAF50` | PAE and contact agree |
| score high and not in contact | `#2196F3` | PAE confident, no close contact |
| score low and in contact | `#F44336` | Contact exists, low PAE confidence |

Do not generate this view for PDB-only structures unless the spec defines a PAE-like signal.

## Coloring Strategy

Prefer annotation-file based coloring if the installed MolViewSpec version supports it. This avoids one component node per residue.

Fallback to per-residue `ComponentExpression` loops when annotation coloring is unavailable.

Use this size guard for per-residue loops:

- `<= 500` residues across rendered chains: per-residue components are acceptable.
- `> 500` residues: collapse to contiguous residue runs or score bands.
- If banding cannot preserve the intended view, emit a warning and produce a simpler chain-level or interface-only view.

Never generate thousands of per-residue component nodes without a guard.

## Display Pattern

Each notebook should define a single display helper and use it for every view.

```python
import base64
import json
from IPython.display import HTML, IFrame, display


def show_mol_view(state, label, width=950, height=600):
    try:
        html = state.molstar_html()
        encoded = base64.b64encode(html.encode()).decode()
        display(HTML(f'<div style="font-weight:bold;margin:8px 0 4px">{label}</div>'))
        display(IFrame(src=f'data:text/html;base64,{encoded}', width=width, height=height))
    except AttributeError:
        mvsj = state.as_dict() if hasattr(state, "as_dict") else state.to_dict()
        encoded = base64.b64encode(json.dumps(mvsj).encode()).decode()
        url = f'https://molstar.org/viewer/?mvs-data={encoded}&mvs-data-format=mvsj'
        display(HTML(f'<b>{label}</b><br><a href="{url}" target="_blank">Open in Mol* viewer</a>'))
```

If the notebook runtime blocks data URI iframes, use the Mol* viewer URL fallback.

## Notebook Integration Rules

- Wrap the full MolViewSpec section in `try/except ImportError`.
- Keep MolViewSpec imports local to the visualization section.
- Do not let later required cells depend on MolViewSpec-only variables.
- Do not introduce non-MolViewSpec 3D protein viewers.
- Put each view in its own `try/except Exception` block.
- Print a summary of generated and skipped views.
- Keep structure URLs and generated state provenance visible.
- Export `.mvsj` only when the spec requests portable visualization states.
- Avoid clinical or biological conclusions in visualization labels unless the reviewed spec supports them.

Recommended section shape:

```text
Markdown: what the 3D views show and their limits
Code: import MolViewSpec, define helpers, verify API assumptions
Code: generate View 1
Code: generate View 2
Code: generate optional Views 3 and 4
Markdown or code output: generated/skipped view summary
```

## Validation Checklist

Before calling the view implementation complete:

- The section degrades gracefully when `molviewspec` is unavailable.
- The chosen parse format is accepted by the installed version.
- Chain selectors match parsed chain IDs.
- Residue selectors use the same residue numbering as the structure file.
- pLDDT, B-factor, and score arrays align with `resids_per_chain`.
- PAE-dependent views are skipped when PAE is absent.
- Large structures do not create unbounded per-residue component trees.
- Inline rendering works or the Mol* viewer fallback is emitted.
- No scoring or data-loading logic depends on MolViewSpec state.

## Boundaries

This skill owns visualization state generation only. It does not own:

- AFDB or PDB fetching
- mmCIF parsing
- PAE parsing
- score formulas
- interface detection
- scientific interpretation beyond visualization labels

If any of those upstream artifacts are missing or inconsistent, report the missing contract instead of patching unrelated notebook sections.
