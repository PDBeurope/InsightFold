# Notebook Spec: EnzyMM on AlphaFold Structures

**Version:** 1.0  
**Status:** Ready for engineering  
**Last updated:** April 2026  
**Project:** InsightFold — Notebook-Driven Development Layer, AFDB/PDBe, EMBL-EBI

---

## Overview

This notebook enables a researcher to take any protein in the AlphaFold
Database (AFDB) — or supply their own AlphaFold-predicted mmCIF file — and
determine whether its three-dimensional structure contains a known arrangement
of catalytic residues, using the Enzyme Motif Miner (EnzyMM). EnzyMM performs
geometric template matching against a curated library of 6780 catalytic site
templates derived from the Mechanism and Catalytic Site Atlas (M-CSA), returning
hits described by RMSD, orientation score, log E-value, matched residue
identities, EC numbers, CATH domain identifiers, and M-CSA entry IDs. The
notebook is aimed at researchers who are familiar with their protein of
interest but are not structural bioinformaticians: it translates raw EnzyMM
output into a layered biological interpretation, overlays AlphaFold confidence
(pLDDT) on every matched residue, generates an adaptive plain-language
narrative that changes depending on whether the protein is a characterised
enzyme, a reviewed protein with unknown function, or a completely uncharacterised
entry, and renders interactive 3D visualisations using MolViewSpec and Mol*.
The notebook is a validation instrument for the InsightFold framework: it is
designed to be cheap to run, easy to interpret, and to surface the biological
signal needed to decide whether population-scale pre-computation across the
full AFDB is justified.

---

## Scientific Context

### What catalytic motifs are and why they matter

Enzymes catalyse reactions by positioning a small set of amino acid side
chains — the catalytic residues — in a precise three-dimensional arrangement
that lowers the activation energy of a specific chemical transformation. This
arrangement, the catalytic motif or active site geometry, is among the most
conserved features in biology. Because the chemistry is constrained by the
twenty amino acids available to all life, the same geometric solution to a
catalytic problem emerges repeatedly across evolutionary history, sometimes in
proteins with no detectable sequence similarity and entirely different global
folds. A serine protease catalytic triad (Ser, His, Asp) appears in chymotrypsin-
like and subtilisin-like families separated by billions of years of evolution,
with no shared ancestry. This conservation of geometry in the absence of
sequence similarity is the biological premise that makes geometric template
matching both possible and scientifically meaningful.

The Mechanism and Catalytic Site Atlas (M-CSA) is the primary curated resource
for enzyme catalytic mechanisms and active sites. As of April 2026 it contains
1003 hand-curated entries covering 895 EC numbers, each describing the catalytic
residues, their roles in the reaction mechanism, and the experimental structures
from which the annotations were derived. Templates in EnzyMM are built from
clusters of homologous experimental structures so that each template represents
a consensus active site conformation rather than a single crystal snapshot.
A given template residue may permit a small set of chemically equivalent amino
acids — for example Asp or Glu, or Ser, Thr, or Tyr — where conservative
substitutions are observed across an enzyme family. Larger active sites are
subdivided via k-means clustering into partial-site templates describing
sub-arrangements of fewer residues. This gives a total library of 6780 templates
from 762 M-CSA enzyme families, each template annotated with the EC number,
CATH domain identifiers, UniProt accession, and M-CSA entry ID of the enzyme
from which it was derived.

### How EnzyMM finds catalytic motifs

EnzyMM uses PyJess, a Cython wrapper of the JESS geometric hashing algorithm,
to search a query protein structure against the template library. Matching is
purely geometric: it compares the three-dimensional positions and orientations
of functional atoms from amino acid side chains against template constraints,
with no sequence alignment or global fold comparison involved. This makes it
capable of detecting functional similarity across evolutionary distances that
are completely invisible to BLAST, HMM-based methods, or global fold comparison
tools such as Foldseek.

Each residue in a template is represented by three functional atoms selected
according to the residue's annotated catalytic role. Templates define match
mode codes that control the strictness of atom-type and residue-type matching,
and a dynamic per-atom distance threshold that reflects known conformational
variation within the cluster. EnzyMM applies RMSD and residue orientation
filters to eliminate spurious matches before returning results.

### What it means to find a catalytic motif in a predicted structure

A match between an AF2-predicted structure and an M-CSA template means that
the three-dimensional arrangement of side chain atoms in the prediction
geometrically resembles a known catalytic geometry from experimental crystal
structures. This is a structural hypothesis, not a confirmed activity. Its
strength depends on: the geometric quality of the match (RMSD, orientation
score, log E-value); the size and selectivity of the matched template; whether
all members of the same conformational cluster also matched (completeness);
the per-residue pLDDT confidence at the matched residues; and whether the
protein has prior annotation that is consistent or inconsistent with the hit.

A strong match — low RMSD, low orientation score, log E-value below −3,
completeness true, predicted_correct true, matched residues all above pLDDT 70
— is meaningful positive evidence for a specific enzymatic function and
mechanism. A weak or partial match should be presented as a lead for
investigation, not a conclusion.

### Key caveats the notebook must communicate

1. **M-CSA coverage is partial.** The atlas covers 1003 enzyme families.
   No EnzyMM hit does not mean the protein is not an enzyme; it may mean the
   mechanism is not yet curated.

2. **EnzyMM detects geometry, not chemistry.** Geometric similarity to a known
   active site is a structural hypothesis. Experimental validation is required
   to confirm catalytic activity.

3. **AF2 predicts a single rotamer state.** Catalytic residues often adopt
   multiple functionally relevant conformations in crystal structures. A
   mismatch between the predicted rotamer and the template's preferred
   orientation inflates the orientation score even when the backbone is correct.

4. **Cofactors are absent.** AF2 does not model metal ions or small molecules.
   Metalloenzyme geometry may match on the protein side while the catalytic
   mechanism requires a metal that is not present in the model.

5. **Induced fit and flexible active sites.** Some active sites are disordered
   in the apo state and only become catalytically competent upon substrate
   binding. AF2 predicts the apo ground state. A matched residue in the 50–70
   pLDDT range may represent a real active site in its open conformation.

6. **Multimeric active sites require multimeric models.** If a template spans
   multiple chains (template_multimeric = True) but the query is a monomer
   model, only part of the active site can be detected.

7. **Validated performance on the human proteome:** 42% recall of annotated
   enzymes, 6.4% false positive rate on non-enzymes, 83.7% precision at the
   residue level for matched enzymes with UniProt active site annotations.

---

## Input Specification

### Accepted input modes

**Mode A — UniProt accession (default):**  
The user provides a raw UniProt accession (e.g. `P04063`). The notebook
constructs all required URLs internally. Only accessions present in the AFDB
are valid; the notebook validates this via the AFDB prediction API before
proceeding.

**Mode B — Local mmCIF upload:**  
The user uploads a locally stored mmCIF file. Only `.cif` and `.mmcif`
extensions are accepted. PDB format is rejected. The file must contain pLDDT
values in the B-factor field (standard for all AFDB and ColabFold outputs) for
the pLDDT filtering and annotation logic to function. If B-factors are not in
the pLDDT range (0–100), a warning is shown.

### Validation logic

| Check | Condition | Action |
|---|---|---|
| Accession format | Does not match `^[A-Z][0-9][A-Z0-9]{3}[0-9]$` or known AF multi-fragment pattern | Raise error with message |
| AFDB API response | HTTP 4xx or empty list | Raise error: accession not found in AFDB |
| UniProt API response | HTTP 4xx | Warn: annotation tier cannot be determined; proceed as `uncharacterised` |
| mmCIF file (local) | Extension not `.cif` or `.mmcif` | Reject with message |
| mmCIF B-factors | All B-factors outside 0–100 range | Warn: file may not be an AF model; pLDDT filtering will be unreliable |
| COMPLEX_TYPE = multimer, single chain detected | Only one chain in the mmCIF | Warn: multimer mode selected but only one chain found |
| EnzyMM TSV empty | Zero rows after header | Handle as zero-hit case (see Output Parsing section) |

### Data fetching steps

**Step 1 — AFDB prediction metadata**

```
GET https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT_ACCESSION}
```

Expected response: a JSON array containing one object per AF model fragment.
Use `entries[0]` for single-fragment proteins. Extract:

- `cifUrl` — AFDB mmCIF download URL (always present)
- `bcifUrl` — binary CIF URL (preferred for MolViewSpec; may be absent in
  older entries, fall back to `cifUrl` with format `mmcif`)
- `uniprotAccession` — for confirmation
- `uniprotDescription` — protein name for display
- `organismScientificName` — for display

Download the mmCIF from `cifUrl` using a streaming GET request with a 60-second
timeout. Write to a temporary file at `/tmp/{accession}.cif`. This file is
passed directly to EnzyMM.

**Step 2 — UniProt entry**

```
GET https://rest.uniprot.org/uniprotkb/{UNIPROT_ACCESSION}.json
```

Expected response: a JSON object. Extract:

- `proteinDescription.recommendedName.fullName.value` — protein name
- `organism.scientificName` — organism
- `keywords` — check for `Enzyme` keyword and EC-related keywords
- `features` — filter for `type == "Active site"` to extract
  `location.start.value` (residue positions, 1-indexed) into
  `uniprot_active_site_residues: list[int]`
- `proteinExistence` — for information
- Determine `PROTEIN_MODE`:
  - Entry is in Swiss-Prot (`entryType == "UniProtKB reviewed (Swiss-Prot)"`)
    AND `recommendedName` contains an EC number or `ecNumbers` field is
    non-empty → `annotated_enzyme`
  - Entry is Swiss-Prot AND no EC number → `annotated_no_function`
  - Entry is TrEMBL or request fails → `uncharacterised`

**Step 3 — Error handling**

All HTTP requests must be wrapped in try/except. On connection error or timeout:
print a warning, set sensible defaults (`PROTEIN_MODE = 'uncharacterised'`,
`uniprot_active_site_residues = []`), and continue. The notebook must not
abort on a failed UniProt fetch — the EnzyMM run is independent.

---

## EnzyMM Execution

### Installation

```python
import subprocess, sys
subprocess.run(
    [sys.executable, '-m', 'pip', 'install', '-q', 'enzymm==0.3.1', 'molviewspec'],
    check=False
)
```

Pin to version 0.3.1, which is the validated version used during notebook
development. Include a version check after import:

```python
import enzymm
assert enzymm.__version__ == '0.3.1', \
    f"Unexpected enzymm version: {enzymm.__version__}. Results may differ."
```

### CLI invocation

EnzyMM is invoked via subprocess using its command-line interface. Do not use
the Python API — the CLI provides stable, documented parameter semantics and
produces the TSV output format this spec depends on.

**Primary run (with pLDDT cutoff):**

```python
import subprocess, os, tempfile

tmp_dir    = tempfile.mkdtemp()
out_tsv    = os.path.join(tmp_dir, f'{accession}_filtered.tsv')
pdb_dir    = os.path.join(tmp_dir, 'pdbs_filtered')
os.makedirs(pdb_dir, exist_ok=True)

cmd_filtered = [
    'enzymm',
    '-i', cif_path,
    '-o', out_tsv,
    '--pdbs', pdb_dir,
    '--conservation-cutoff', str(PLDDT_CUTOFF),
    '--skip-smaller-hits',
    '--jobs', '1',           # single thread in Colab to avoid OOM
]
result = subprocess.run(cmd_filtered, capture_output=True, text=True)
if result.returncode != 0:
    raise RuntimeError(f'EnzyMM failed:\n{result.stderr}')
```

**Comparison run (without pLDDT cutoff), executed only when RUN_COMPARISON = True:**

```python
out_tsv_unfiltered = os.path.join(tmp_dir, f'{accession}_unfiltered.tsv')
pdb_dir_unfiltered = os.path.join(tmp_dir, 'pdbs_unfiltered')
os.makedirs(pdb_dir_unfiltered, exist_ok=True)

cmd_unfiltered = [
    'enzymm',
    '-i', cif_path,
    '-o', out_tsv_unfiltered,
    '--pdbs', pdb_dir_unfiltered,
    '--skip-smaller-hits',
    '--jobs', '1',
]
result2 = subprocess.run(cmd_unfiltered, capture_output=True, text=True)
if result2.returncode != 0:
    raise RuntimeError(f'EnzyMM (unfiltered) failed:\n{result2.stderr}')
```

### Parameter rationale

| Parameter | Value | Rationale |
|---|---|---|
| `--conservation-cutoff` | 70 (default, user-overridable) | The boundary between the light blue (70–90) and yellow (50–70) pLDDT bands. Below 70, backbone geometry becomes unreliable and side chain positions cannot be trusted for geometric matching. |
| `--skip-smaller-hits` | True (default) | Suppresses smaller partial templates when a match to a larger template from the same region has already been found. Reduces noise in the output. |
| `--jobs` | 1 | Colab instances have limited CPU and memory. Single-threaded execution prevents resource exhaustion. |
| `--pdbs` | temp directory | Writes one PDB file per match in the query coordinate frame. Used for residue position extraction; not directly shown to the user. |

### Expected output files

- `{accession}_filtered.tsv` — full results table with pLDDT cutoff applied.
  May be empty (header only) if no hits pass the filter.
- `{accession}_unfiltered.tsv` — full results table without cutoff.
  Produced only when `RUN_COMPARISON = True`.
- `pdbs_filtered/` — one PDB per match from the filtered run. Each PDB
  contains MODEL 0 (the full query structure) and one additional MODEL per
  match (matched residues only in the query coordinate frame).

---

## Output Parsing and Interpretation

### TSV parsing

Parse the TSV using pandas, skipping comment lines (those beginning with `#`):

```python
import pandas as pd

def parse_enzymm_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', comment='#')
    if df.empty:
        return df
    # Parse list-valued fields from string representation
    df['template_ec']   = df['template_ec'].apply(
        lambda x: [e.strip().strip("'") for e in x.strip('[]').split(',')]
        if isinstance(x, str) and x.strip('[]') else []
    )
    df['template_cath'] = df['template_cath'].apply(
        lambda x: [e.strip().strip("'") for e in x.strip('[]').split(',')]
        if isinstance(x, str) and x.strip('[]') else []
    )
    # Parse matched_residues into a list of dicts
    def parse_residues(s):
        parts = [r.strip().strip("'") for r in s.strip('[]').split(',')]
        result = []
        for p in parts:
            tokens = p.split('_')
            if len(tokens) == 3:
                result.append({'name': tokens[0], 'chain': tokens[1],
                               'resnum': int(tokens[2])})
        return result
    df['matched_residues_parsed'] = df['matched_residues'].apply(parse_residues)
    return df
```

### Fields used and their interpretation

The following fields are consumed by the notebook. All others remain in the
DataFrame but are not surfaced in the default user-facing output.

| Field | Type | Interpretation |
|---|---|---|
| `template_ec` | list[str] | EC numbers of the matched template. May contain multiple entries. |
| `template_mcsa_id` | int | M-CSA entry ID. Used to generate hyperlinks and group hits from the same enzyme family. |
| `template_cath` | list[str] | CATH domain identifiers. Used to determine whether multiple hits share a structural superfamily. |
| `rmsd` | float | Atom-wise superposition distance in Å. Hard ceiling 2 Å. Lower is better. Colour-graded in the display table: < 0.5 green, 0.5–1.0 yellow, > 1.0 orange. |
| `log_evalue` | float | Statistical match quality. More negative is better. Below −3 is acceptable; below −4 is strong. |
| `orientation` | float | Mean pairwise side chain orientation angle in radians. Lower is better. Primary sort key within tiers. The most sensitive metric for distinguishing genuine catalytic geometry from coincidental superposition. |
| `completeness` | bool | Whether all cluster members matched. False indicates a partial active site detection. |
| `predicted_correct` | bool | EnzyMM's internal binary quality verdict. Tier 1 requires this to be True. |
| `template_effective_size` | int | Number of specific side-chain-interacting residues. Larger = more selective template. |
| `template_multimeric` | bool | Whether the template spans multiple chains. |
| `query_multimeric` | bool | Whether matched residues span multiple chains in the query. |
| `number_of_metal_ligands` | tuple(int, int) | Metal-coordinating residues in template and reference. High values indicate a metal-binding template, which is less predictive of full catalytic function. |
| `matched_residues` | str | Matched residues in format `[3-letter-code]_[chain]_[resnum]`. Parsed by `parse_residues()`. |
| `total_reference_residues` | int | Total catalytic residues in the reference structure. If `template_effective_size` < `total_reference_residues`, the template represents a partial site. |

### Hit filtering and prioritisation

**Step 1 — Tier assignment**

```python
df['tier'] = df.apply(
    lambda r: 1 if (r['completeness'] and r['predicted_correct']) else 2,
    axis=1
)
tier1 = df[df['tier'] == 1].copy()
tier2 = df[df['tier'] == 2].copy()
```

**Step 2 — Sort within Tier 1**

```python
tier1 = tier1.sort_values(
    by=['orientation', 'rmsd', 'log_evalue'],
    ascending=[True, True, False]
).reset_index(drop=True)
```

**Step 3 — Group by M-CSA entry**

Within each tier, group rows sharing the same `template_mcsa_id`. When the
same M-CSA entry appears multiple times (different conformational cluster
members), collapse into a single representative row using the member with the
lowest orientation score. Retain `cluster_size` in the display to indicate
how many members matched.

**Step 4 — Comparison delta (when RUN_COMPARISON = True)**

```python
filtered_ids   = set(df_filtered['template_mcsa_id'].astype(str)
                     + '_' + df_filtered['match_index'].astype(str))
unfiltered_ids = set(df_unfiltered['template_mcsa_id'].astype(str)
                     + '_' + df_unfiltered['match_index'].astype(str))
low_conf_hits  = df_unfiltered[
    (df_unfiltered['template_mcsa_id'].astype(str)
     + '_' + df_unfiltered['match_index'].astype(str))
    .isin(unfiltered_ids - filtered_ids)
].copy()
low_conf_hits['tier'] = 'low_confidence'
```

These hits are appended to the Tier 2 display with the label:
"Detected only without pLDDT filter — matched residues are in a low-confidence
region (pLDDT < {PLDDT_CUTOFF})."

**Step 5 — Multi-EC summary logic**

Computed over all Tier 1 hits and used in the narrative section:

```python
all_ec     = [ec for row in tier1['template_ec'] for ec in row]
all_cath   = [c  for row in tier1['template_cath'] for c in row]
ec_classes = list({ec.split('.')[0] for ec in all_ec})  # first digit only

shared_cath = len(set(all_cath)) < len(all_cath)         # any CATH repeated

# Shared residue core across hits
from collections import Counter
all_resnums = [r['resnum'] for row in tier1['matched_residues_parsed']
               for r in row]
shared_resnums = {k for k, v in Counter(all_resnums).items() if v > 1}

ec_related = len(ec_classes) == 1       # all same first digit
ec_divergent = len(ec_classes) > 1
```

### Interpreting match quality

| Criterion | Strong | Acceptable | Weak |
|---|---|---|---|
| RMSD (Å) | < 0.5 | 0.5–1.0 | 1.0–2.0 |
| log E-value | < −4 | −4 to −3 | −3 to 0 |
| Orientation (rad) | < 0.3 | 0.3–0.6 | > 0.6 |
| template_effective_size | ≥ 5 | 3–4 | ≤ 2 |
| completeness | True | — | False |
| predicted_correct | True | — | False |

A hit is considered reliable for biological interpretation when it meets
"strong" or "acceptable" on RMSD and orientation, has predicted_correct = True,
and has matched residues with pLDDT ≥ 70.

A hit is considered unreliable and must carry an explicit warning when:
- Any matched residue has pLDDT < 70
- template_multimeric = True and query_multimeric = False
- number_of_metal_ligands (template value) ≥ 3 (template likely describes a
  metal binding site rather than a full active site)
- The hit is present only in the unfiltered run

### Zero-hit case

When the filtered TSV is empty, the notebook displays:

> **No catalytic motif matches found**
>
> EnzyMM did not find any known catalytic residue arrangements in this
> structure that passed the quality filters. This does not mean the protein
> is not an enzyme. There are two common reasons for a no-hit result:
>
> 1. The protein's catalytic mechanism is not yet described in the Mechanism
>    and Catalytic Site Atlas (M-CSA), which currently covers 1003 enzyme
>    families. Novel or poorly characterised enzyme families are often absent.
>
> 2. The matched residues fall in a low-confidence region (pLDDT < 70) and
>    were excluded by the pLDDT filter. Try setting PLDDT_CUTOFF to 0 and
>    enabling RUN_COMPARISON to see if any matches appear in the unfiltered run.
>
> If RUN_COMPARISON = True and the unfiltered run also returned no hits, the
> structure is genuinely outside the current M-CSA coverage.

---

## Visualisations

### Visualisation 1 — pLDDT sequence strip with matched residue annotations

**What it shows:** The per-residue pLDDT score across the full protein
sequence, coloured by the four standard AlphaFold confidence bands. Tier 1
matched residues are annotated as vertical markers with residue labels and
their pLDDT values. The user can see at a glance whether the matched residues
are in a high-confidence region and how they relate to the rest of the
sequence.

**Library:** matplotlib

**Data source:** pLDDT values extracted from the B-factor field of the mmCIF
file. Use the first model and chain A (or the only chain for monomers).

```python
PLDDT_BANDS = [
    (90, 100, '#1565C0'),  # dark blue
    (70,  90, '#42A5F5'),  # light blue
    (50,  70, '#FFCA28'),  # yellow
    ( 0,  50, '#EF6C00'),  # orange
]
```

**Implementation:**

The plot is a filled step chart or scatter plot of residue index (x-axis) vs
pLDDT (y-axis, range 0–100). Each point is coloured according to which band
it falls in. A horizontal dashed grey line marks pLDDT = 70 (the filter
threshold). A horizontal dashed grey line at pLDDT = 90 marks the high-
confidence boundary.

For each Tier 1 matched residue, draw a vertical dashed line from y=0 to the
pLDDT value at that residue. At the top of the line, place a text label:
`{3-letter-code}{resnum}\n{pLDDT:.0f}`. Labels for adjacent residues should
be staggered vertically to avoid overlap.

Below the plot, render a small table (as printed text or a styled DataFrame):

| Residue | Position | pLDDT | Band |
|---|---|---|---|
| ARG | 201 | 98.9 | Very high (>90) |
| ASP | 203 | 98.8 | Very high (>90) |

**What the user should understand:** Whether the residues that EnzyMM matched
are in a region AlphaFold predicted confidently. A matched residue with
pLDDT > 90 means the backbone and likely the rotamer are reliable. A matched
residue with pLDDT 70–90 means the backbone is probably correct but the
rotamer may differ from the true position. A matched residue below 70 means
the geometry should not be trusted; such hits appear only in the unfiltered
run and are flagged accordingly.

---

### Visualisation 2 — 3D structure coloured by pLDDT with matched residues

**What it shows:** The full protein structure rendered as a cartoon, coloured
by pLDDT using the four standard confidence bands. Tier 1 hit 1 matched
residues are shown as ball-and-stick in a neutral highlight colour (`#FFFFFF`
white, or the user may prefer `#E040FB` magenta for visibility — to be decided
during engineering) overlaid on the confidence-coloured cartoon.

**Library:** MolViewSpec / Mol* (molviewspec Python package)

**Implementation pattern** (adapted from the homodimer notebook):

```python
import molviewspec as mvs
import base64
from IPython.display import display, HTML, IFrame

def show_mol_view(state, label, width=950, height=600):
    html = state.molstar_html()
    encoded = base64.b64encode(html.encode()).decode()
    display(HTML(f'<div style="margin:10px 0 4px;font-weight:bold;">{label}</div>'))
    display(IFrame(src=f'data:text/html;base64,{encoded}', width=width, height=height))

struct_url = meta.get('bcifUrl', meta.get('cifUrl'))
struct_fmt  = 'bcif' if 'bcifUrl' in meta else 'mmcif'

builder = mvs.create_builder()
structure = (
    builder
    .download(url=struct_url)
    .parse(format=struct_fmt)
    .model_structure()
)

# Apply pLDDT colouring band by band
for chain_id in chain_ids:
    for lo, hi, colour in PLDDT_BANDS:
        resnums_in_band = [
            r for r in all_resnums_for_chain
            if lo <= plddt_by_resnum[r] < hi
        ]
        for r in resnums_in_band:
            (structure
             .component(selector=mvs.ComponentExpression(
                 label_asym_id=chain_id,
                 beg_label_seq_id=r, end_label_seq_id=r))
             .representation(type='cartoon')
             .color(color=colour))

# Highlight matched residues as ball-and-stick
top_hit_residues = tier1.iloc[0]['matched_residues_parsed']
for res in top_hit_residues:
    (structure
     .component(selector=mvs.ComponentExpression(
         label_asym_id=res['chain'],
         beg_label_seq_id=res['resnum'],
         end_label_seq_id=res['resnum']))
     .representation(type='ball_and_stick')
     .color(color='#E040FB'))   # magenta; to be confirmed during engineering

show_mol_view(
    builder.get_state(),
    f'View 1: pLDDT confidence map — matched catalytic residues shown as sticks '
    f'(top hit: M-CSA {tier1.iloc[0]["template_mcsa_id"]}, EC {tier1.iloc[0]["template_ec"]})'
)
```

**Note on label_seq_id vs author_seq_id:** In AFDB mmCIF files, label_seq_id
and author_seq_id are typically identical. The matched_residues field from
EnzyMM uses author residue numbers from the PDB format output. The engineering
session must verify this alignment for the specific AF model version in use,
or add an explicit mapping step. See Open Questions.

**What the user should understand:** Where in the 3D structure the matched
catalytic residues sit, and how confident AlphaFold was in that region. If the
matched residues are coloured dark blue and sit in a well-defined pocket, the
geometry is trustworthy. If they are in an orange or yellow region, caution
is needed.

---

### Visualisation 3 — Multi-hit residue overlap (conditional, when Tier 1 has > 1 hit)

**What it shows:** The structure rendered as a grey cartoon. Matched residues
from all Tier 1 hits are shown as ball-and-stick coloured by hit identity.
Residues shared across two or more hits are coloured green (`#43A047`). Hit-
specific residues for hit 1 are coloured blue (`#1565C0`) and for hit 2 orange
(`#EF6C00`). This directly visualises the shared catalytic core versus hit-
specific peripheral residues, making convergent evidence immediately legible.

**Library:** MolViewSpec / Mol*

**Implementation:**

```python
if len(tier1) > 1:
    from collections import Counter

    all_matched = [
        (res['chain'], res['resnum'])
        for _, row in tier1.iterrows()
        for res in row['matched_residues_parsed']
    ]
    counts = Counter(all_matched)
    shared = {k for k, v in counts.items() if v > 1}

    builder3 = mvs.create_builder()
    structure3 = (builder3.download(url=struct_url)
                           .parse(format=struct_fmt)
                           .model_structure())

    # Full structure — grey cartoon
    (structure3
     .component(selector=mvs.ComponentExpression())
     .representation(type='cartoon')
     .color(color='#BDBDBD'))

    # Hit-specific residues
    hit_colours = ['#1565C0', '#EF6C00', '#AB47BC', '#00838F']
    for hit_idx, (_, row) in enumerate(tier1.iterrows()):
        colour = hit_colours[min(hit_idx, len(hit_colours) - 1)]
        for res in row['matched_residues_parsed']:
            key = (res['chain'], res['resnum'])
            use_colour = '#43A047' if key in shared else colour
            (structure3
             .component(selector=mvs.ComponentExpression(
                 label_asym_id=res['chain'],
                 beg_label_seq_id=res['resnum'],
                 end_label_seq_id=res['resnum']))
             .representation(type='ball_and_stick')
             .color(color=use_colour))

    show_mol_view(
        builder3.get_state(),
        'View 2: Matched residues by hit — shared residues in green, '
        'hit-specific residues in blue/orange'
    )
```

**What the user should understand:** When two hits match largely the same
residues (as in the P04063 example, where 4 of 5 residues were shared), this
is convergent evidence for a conserved catalytic strategy, not a conflicting
interpretation. The green residues are the consensus core of the detected
active site geometry.

---

## Notebook Sections (with cell-level detail)

---

### Section 0 — Title and user guide

**Cell 0.1 — Markdown: Title block**

```markdown
# EnzyMM — Catalytic Motif Search on AlphaFold Structures

**Purpose:** This notebook searches your AlphaFold protein structure for known
arrangements of catalytic residues — the specific amino acids that carry out
enzymatic reactions. It tells you whether your protein's predicted 3D structure
matches any of the 6780 catalytic templates curated in the Mechanism and
Catalytic Site Atlas (M-CSA), and what those matches imply about the protein's
possible function and mechanism.

**Audience:** Researchers familiar with their protein of interest but not
necessarily expert in structural bioinformatics. Every section includes a
plain-language explanation before any code.

**Usage:** Enter your UniProt accession in Section 2, adjust any parameters,
then run all cells (Runtime → Run all).
```

This cell contains no code. It is the entry point for the user.

---

### Section 1 — Installation and imports

**Cell 1.1 — Code: Install packages**

Installs `enzymm==0.3.1` and `molviewspec` via pip in quiet mode. Runs as a
subprocess to suppress output. Prints a single confirmation line on completion.
This cell must be the first code cell so that imports in subsequent cells
succeed. Takes approximately 30–60 seconds in a fresh Colab runtime.

```python
import subprocess, sys
subprocess.run(
    [sys.executable, '-m', 'pip', 'install', '-q',
     'enzymm==0.3.1', 'molviewspec'],
    check=False
)
print('Packages installed.')
```

**Cell 1.2 — Code: Imports and global constants**

Imports all libraries used across the notebook. Defines the pLDDT band colour
scheme and the `show_mol_view()` helper function. Defines the `parse_residues()`
helper. Sets matplotlib display parameters. Prints "Imports OK" on success.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import requests, json, io, os, re, tempfile, subprocess, sys
import textwrap, warnings, base64
from collections import Counter
from pathlib import Path
from IPython.display import display, HTML, IFrame
import gemmi
import molviewspec as mvs
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 12, 'axes.labelsize': 12, 'axes.titlesize': 13,
    'xtick.labelsize': 10, 'ytick.labelsize': 10, 'figure.dpi': 150,
})

PLDDT_BANDS = [
    (90, 100, '#1565C0'),
    (70,  90, '#42A5F5'),
    (50,  70, '#FFCA28'),
    ( 0,  50, '#EF6C00'),
]

def plddt_colour(value: float) -> str:
    for lo, hi, colour in PLDDT_BANDS:
        if lo <= value < hi or (value == 100 and hi == 100):
            return colour
    return '#EF6C00'

def plddt_band_label(value: float) -> str:
    if value >= 90: return 'Very high (>90)'
    if value >= 70: return 'Confident (70–90)'
    if value >= 50: return 'Low (50–70)'
    return 'Very low (<50)'

def show_mol_view(state, label: str, width: int = 950, height: int = 600):
    html = state.molstar_html()
    encoded = base64.b64encode(html.encode()).decode()
    display(HTML(
        f'<div style="margin:10px 0 4px;font-weight:bold;'
        f'font-family:sans-serif;">{label}</div>'
    ))
    display(IFrame(
        src=f'data:text/html;base64,{encoded}', width=width, height=height
    ))

def parse_enzymm_tsv(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep='\t', comment='#')
    except Exception:
        return pd.DataFrame()
    if df.empty or len(df.columns) < 5:
        return pd.DataFrame()
    def parse_list(s):
        if not isinstance(s, str): return []
        return [e.strip().strip("'\"") for e in s.strip('[]').split(',')
                if e.strip().strip("'\"")] 
    def parse_residues(s):
        if not isinstance(s, str): return []
        result = []
        for p in [x.strip().strip("'\"") for x in s.strip('[]').split(',')]:
            tokens = p.split('_')
            if len(tokens) == 3:
                try:
                    result.append({'name': tokens[0], 'chain': tokens[1],
                                   'resnum': int(tokens[2])})
                except ValueError:
                    pass
        return result
    df['template_ec']             = df['template_ec'].apply(parse_list)
    df['template_cath']           = df['template_cath'].apply(parse_list)
    df['matched_residues_parsed'] = df['matched_residues'].apply(parse_residues)
    return df

print('Imports OK.')
```

**Cell 1.3 — Code: Build M-CSA enzyme name lookup table from EnzyMM template directory**

EnzyMM ships ~16 MB of template files, each a PDB-like text file containing
`REMARK MCSA_ID` and `REMARK ENZYME` lines. This cell reads those files once
at startup and builds a `dict[int, str]` mapping M-CSA entry ID to enzyme
name. This replaces any runtime API call to the M-CSA or UniProt. The lookup
is used in Section 10 to display enzyme names alongside M-CSA hyperlinks.

The template directory location is resolved from the installed `enzymm`
package path. If the directory cannot be found (e.g. unexpected install
layout), the dict remains empty and Section 10 falls back to displaying only
the M-CSA ID and hyperlink with no enzyme name.

```python
import enzymm as _enzymm_pkg

def _build_mcsa_lookup() -> dict:
    lookup = {}
    try:
        tmpl_dir = Path(_enzymm_pkg.__file__).parent / 'data' / 'templates'
        if not tmpl_dir.exists():
            # Try alternate install path
            tmpl_dir = Path(_enzymm_pkg.__file__).parent / 'templates'
        for tmpl_file in tmpl_dir.glob('*.pdb'):
            mcsa_id = None
            enzyme_name = None
            for line in tmpl_file.read_text(errors='replace').splitlines():
                if line.startswith('REMARK MCSA_ID'):
                    try:
                        mcsa_id = int(line.split()[-1])
                    except (ValueError, IndexError):
                        pass
                elif line.startswith('REMARK ENZYME'):
                    enzyme_name = line.replace('REMARK ENZYME', '', 1).strip()
                if mcsa_id is not None and enzyme_name is not None:
                    lookup[mcsa_id] = enzyme_name
                    break
    except Exception as e:
        print(f'Warning: could not build M-CSA lookup table: {e}')
    return lookup

MCSA_NAME_LOOKUP = _build_mcsa_lookup()
print(f'M-CSA lookup table built: {len(MCSA_NAME_LOOKUP)} entries.')
```

**Cell 1.4 — Code: EnzyMM version check**

After import, verifies the installed EnzyMM version matches the expected
version. Prints a warning (not an error) if the version differs, so the
notebook continues but the user is informed.

```python
try:
    _ver = getattr(_enzymm_pkg, '__version__', 'unknown')
    if _ver != '0.3.1':
        print(f'WARNING: enzymm version {_ver} found; this notebook was '
              f'developed against 0.3.1. Results may differ.')
    else:
        print(f'enzymm {_ver} OK.')
except ImportError:
    raise ImportError('enzymm could not be imported. Re-run Cell 1.1.')
```

---

### Section 2 — User parameters

**Cell 2.1 — Markdown: Parameter guide**

```markdown
## Section 2 — Set your parameters

Enter your UniProt accession below (e.g. `P04063`) and adjust any settings.
Then run all cells (Runtime → Run all).

- **UNIPROT_ACCESSION** — the UniProt ID of your protein. Must be present in
  the AlphaFold Database.
- **USE_LOCAL_FILE** — set to `True` if your protein is not in AFDB and you
  want to upload your own mmCIF file.
- **PLDDT_CUTOFF** — residues with pLDDT below this value are excluded from
  EnzyMM matching. Default 70. Set to 0 to disable the filter entirely.
- **COMPLEX_TYPE** — `'monomer'` for single-chain AF models; `'multimer'`
  if you are uploading a multi-chain AF3 or ColabFold model.
- **RUN_COMPARISON** — if True, EnzyMM runs twice: once with the pLDDT
  filter, once without. Shows you which hits disappear under the filter.
```

**Cell 2.2 — Code: Parameter block**

This cell holds all user-facing parameters as simple Python assignments with
Colab `@param` widget annotations. No computation happens here.

```python
# ── USER PARAMETERS ─────────────────────────────────────────────────────────

UNIPROT_ACCESSION = 'P04063'     # @param {type:"string"}
USE_LOCAL_FILE    = False         # @param {type:"boolean"}
PLDDT_CUTOFF      = 70           # @param {type:"integer"}
COMPLEX_TYPE      = 'monomer'    # @param ["monomer", "multimer"]
RUN_COMPARISON    = True         # @param {type:"boolean"}
SKIP_SMALLER_HITS = True         # @param {type:"boolean"}

# ────────────────────────────────────────────────────────────────────────────
print(f'Accession  : {UNIPROT_ACCESSION}')
print(f'pLDDT cut  : {PLDDT_CUTOFF}')
print(f'Mode       : {COMPLEX_TYPE}')
print(f'Comparison : {RUN_COMPARISON}')
```

**Cell 2.3 — Code: File upload widget (conditional)**

This cell renders an ipywidgets `FileUpload` widget when `USE_LOCAL_FILE=True`,
accepting only `.cif` and `.mmcif` files. When `USE_LOCAL_FILE=False`, it
prints a confirmation that the notebook will download from AFDB. The widget
value is read in Section 3.

```python
import ipywidgets as _widgets

if USE_LOCAL_FILE:
    _upload_cif = _widgets.FileUpload(
        accept='.cif,.mmcif', multiple=False,
        description='Upload mmCIF'
    )
    display(_widgets.VBox([
        _widgets.HTML(
            '<b>Upload your mmCIF file, then run the next section.</b><br>'
            'Only AlphaFold-format mmCIF files are supported '
            '(pLDDT stored in B-factor field).'
        ),
        _upload_cif,
    ]))
else:
    _upload_cif = None
    print(f'Online mode — structure will be downloaded from AFDB '
          f'for UniProt accession {UNIPROT_ACCESSION}.')
```

---

### Section 3 — Data loading

**Cell 3.1 — Markdown: What happens here**

```markdown
## Section 3 — Loading your protein

We fetch the predicted structure from the AlphaFold Database and the protein
annotation from UniProt. This determines what prior information exists about
your protein: is it a characterised enzyme, a reviewed protein with unknown
function, or completely uncharacterised? The answer shapes how we interpret
the EnzyMM results.
```

**Cell 3.2 — Code: AFDB fetch**

Constructs the AFDB API URL from `UNIPROT_ACCESSION`. Issues a GET request
with a 30-second timeout. Validates that the response is a non-empty JSON
array. Extracts `cifUrl`, `bcifUrl` (if present), `uniprotDescription`, and
`organismScientificName` from `entries[0]`. Downloads the mmCIF file from
`cifUrl` with a 60-second timeout. Writes it to a temporary file at
`/tmp/{UNIPROT_ACCESSION}.cif`. Sets `cif_path` to this path. Sets
`struct_url` to `bcifUrl` if present (for MolViewSpec), else `cifUrl`. Sets
`struct_fmt` accordingly.

When `USE_LOCAL_FILE=True`, reads the mmCIF bytes from the upload widget and
writes them to `/tmp/uploaded_structure.cif`. Sets `cif_path` to this path.
Sets `struct_url = None` (MolViewSpec views will be skipped or degraded —
to be handled in Section 8 with a warning). Sets `UNIPROT_ACCESSION` to
`'uploaded_structure'` for downstream labelling.

On any HTTP error, raises a `RuntimeError` with a user-readable message
naming the accession and the HTTP status code.

After a successful download, prints:
```
Protein    : {uniprotDescription}
Organism   : {organismScientificName}
Structure  : {cifUrl}
mmCIF size : {len(cif_text)} characters
```

**Cell 3.3 — Code: UniProt annotation fetch**

Constructs the UniProt REST URL as
`https://rest.uniprot.org/uniprotkb/{UNIPROT_ACCESSION}.json`.
Issues a GET request with a 15-second timeout. On success, parses the JSON
and determines `PROTEIN_MODE`:

- `annotated_enzyme`: `entryType` contains "reviewed" AND at least one of:
  (a) `proteinDescription.recommendedName.ecNumbers` is non-empty, or (b)
  any item in `keywords` has `name == "Enzyme"`.
- `annotated_no_function`: `entryType` contains "reviewed" AND neither
  condition above is met.
- `uncharacterised`: `entryType` does not contain "reviewed", or the request
  fails.

Extracts `uniprot_active_site_residues` from `features` where
`type == "Active site"`, taking `location.start.value` (1-indexed integers).

On any exception, sets `PROTEIN_MODE = 'uncharacterised'` and
`uniprot_active_site_residues = []`. Prints a warning but does not raise.

Prints:
```
UniProt mode : {PROTEIN_MODE}
Active sites : {uniprot_active_site_residues}
```

**Cell 3.4 — Code: Structure parsing with gemmi**

Uses `gemmi.read_structure()` to parse the mmCIF file. gemmi is installed as
a dependency of EnzyMM and is therefore guaranteed to be available. This single
cell produces three outputs used throughout the notebook:

1. **`plddt_by_chain_resnum: dict[str, dict[int, float]]`** — per-residue
   pLDDT values keyed by chain ID and author residue number.

2. **`auth_to_label: dict[str, dict[int, int]]`** — mapping from
   (chain_id, auth_seq_id) to label_seq_id for every residue. Used in all
   MolViewSpec `ComponentExpression` calls to ensure the residue numbers from
   the EnzyMM `matched_residues` field (author numbering) are correctly
   translated to label_seq_id for the Mol* viewer.

3. **`chain_ids: list[str]`** — all unique chain IDs found in model 0, for
   use in MolViewSpec structure colouring loops.

```python
import gemmi

st = gemmi.read_structure(cif_path)
model = st[0]

plddt_by_chain_resnum: dict[str, dict[int, float]] = {}
auth_to_label:         dict[str, dict[int, int]]   = {}
chain_ids: list[str] = []

for chain in model:
    cid = chain.name
    chain_ids.append(cid)
    plddt_by_chain_resnum[cid] = {}
    auth_to_label[cid]         = {}
    for res in chain:
        auth_num  = res.seqid.num
        label_num = res.label_seq           # gemmi exposes label_seq directly
        # B-factor (pLDDT) — take from the first atom of the residue
        if len(res) > 0:
            bfactor = res[0].b_iso
            plddt_by_chain_resnum[cid][auth_num] = bfactor
        auth_to_label[cid][auth_num] = label_num

# Validate pLDDT range
all_bfactors = [v for chain_dict in plddt_by_chain_resnum.values()
                  for v in chain_dict.values()]
if all_bfactors:
    if max(all_bfactors) <= 1.0:
        print('WARNING: B-factors appear to be on the 0–1 scale rather than '
              '0–100. This file may not be an AlphaFold model. '
              'pLDDT filtering and confidence colouring will be unreliable.')
    elif max(all_bfactors) > 100.0:
        print('WARNING: B-factors exceed 100. This file may not be an '
              'AlphaFold model. pLDDT-dependent features may be unreliable.')
    else:
        print(f'pLDDT values OK — range: '
              f'{min(all_bfactors):.1f}–{max(all_bfactors):.1f}')

print(f'Chains found   : {chain_ids}')
print(f'Total residues : {len(all_bfactors)}')
```

**Cell 3.5 — Code: Inter-chain confidence metric extraction from mmCIF
(executes only when COMPLEX_TYPE = "multimer")**

Parses the `_ma_qa_metric` and `_ma_qa_metric_global` loops directly from the
mmCIF text to extract ipSAE, ipTM, pDockQ2, and LIS. This approach reads the
raw mmCIF text rather than using gemmi's QA metric API, because gemmi does not
expose a stable high-level interface for `_ma_qa_metric_global` in the version
shipped with EnzyMM 0.3.1.

The mmCIF structure of these sections is:

```
loop_
_ma_qa_metric.id
_ma_qa_metric.name
_ma_qa_metric.type
_ma_qa_metric.mode
_ma_qa_metric.software_group_id
1 pLDDT pLDDT global 1
2 pLDDT pLDDT local 1
3 ipTM other global 1
5 'ipsae_AB' other global 1
6 'ipsae_BA' other global 1
13 'pDockQ2_AB' other global 1
14 'pDockQ2_BA' other global 1
15 'LIS_AB' other global 1
16 'LIS_BA' other global 1
...

loop_
_ma_qa_metric_global.ordinal_id
_ma_qa_metric_global.model_id
_ma_qa_metric_global.metric_id
_ma_qa_metric_global.metric_value
1 1 1 89.52
2 1 3 0.81
4 1 5 0.870943
5 1 6 0.867331
...
```

Derived values (as specified in the notebook design decisions):
- **ipSAE** = `max(ipsae_AB, ipsae_BA)` — the more optimistic direction gives
  the tighter lower bound on inter-chain confidence.
- **ipTM** = value of metric named `ipTM`.
- **pDockQ2** = `min(pDockQ2_AB, pDockQ2_BA)` — conservative estimate.
- **LIS** = `min(LIS_AB, LIS_BA)` — conservative estimate.

```python
def parse_ma_qa_metrics(cif_text: str) -> dict[str, float]:
    """
    Parse _ma_qa_metric and _ma_qa_metric_global from mmCIF text.
    Returns a dict of metric_name -> float value for all global metrics.
    """
    import re

    # ── Step 1: build metric_id → name map ──────────────────────────────────
    id_to_name: dict[int, str] = {}
    in_metric_def = False
    metric_def_lines = []

    for line in cif_text.splitlines():
        stripped = line.strip()
        if stripped == 'loop_':
            in_metric_def = False
            metric_def_lines = []
        if '_ma_qa_metric.id' in stripped:
            in_metric_def = True
            continue
        if in_metric_def:
            if stripped.startswith('_'):
                continue          # still in column headers
            if not stripped or stripped.startswith('#'):
                in_metric_def = False
                continue
            if stripped.startswith('loop_'):
                in_metric_def = False
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    metric_id   = int(parts[0])
                    metric_name = parts[1].strip("'")
                    id_to_name[metric_id] = metric_name
                except ValueError:
                    pass

    # ── Step 2: build metric_id → value map ─────────────────────────────────
    id_to_value: dict[int, float] = {}
    in_global = False

    for line in cif_text.splitlines():
        stripped = line.strip()
        if '_ma_qa_metric_global.ordinal_id' in stripped:
            in_global = True
            continue
        if in_global:
            if stripped.startswith('_'):
                continue
            if not stripped or stripped.startswith('#') \
                    or stripped.startswith('loop_'):
                in_global = False
                continue
            parts = stripped.split()
            if len(parts) == 4:
                try:
                    metric_id = int(parts[2])
                    value     = float(parts[3])
                    id_to_value[metric_id] = value
                except ValueError:
                    pass

    # ── Step 3: name → value ─────────────────────────────────────────────────
    return {id_to_name[mid]: val
            for mid, val in id_to_value.items()
            if mid in id_to_name}


if COMPLEX_TYPE == 'multimer':
    with open(cif_path, 'r', errors='replace') as _f:
        _cif_text = _f.read()
    _raw_metrics = parse_ma_qa_metrics(_cif_text)

    # Derive the four target metrics
    ipsae_AB   = _raw_metrics.get('ipsae_AB',   None)
    ipsae_BA   = _raw_metrics.get('ipsae_BA',   None)
    pdockq2_AB = _raw_metrics.get('pDockQ2_AB', None)
    pdockq2_BA = _raw_metrics.get('pDockQ2_BA', None)
    lis_AB     = _raw_metrics.get('LIS_AB',     None)
    lis_BA     = _raw_metrics.get('LIS_BA',     None)
    iptm_val   = _raw_metrics.get('ipTM',       None)

    multimer_scores = {}
    if ipsae_AB is not None and ipsae_BA is not None:
        multimer_scores['ipSAE']   = max(ipsae_AB, ipsae_BA)
    if iptm_val is not None:
        multimer_scores['ipTM']    = iptm_val
    if pdockq2_AB is not None and pdockq2_BA is not None:
        multimer_scores['pDockQ2'] = min(pdockq2_AB, pdockq2_BA)
    if lis_AB is not None and lis_BA is not None:
        multimer_scores['LIS']     = min(lis_AB, lis_BA)

    if multimer_scores:
        print('Inter-chain confidence metrics extracted from mmCIF:')
        for name, val in multimer_scores.items():
            print(f'  {name:<10}: {val:.4f}')
    else:
        print('WARNING: No inter-chain confidence metrics found in this mmCIF. '
              'The file may be a monomer model or an older AF format that does '
              'not embed QA metrics.')
        multimer_scores = {}
else:
    multimer_scores = {}
    print('Monomer mode — inter-chain confidence metrics not applicable.')
```

The `multimer_scores` dict is consumed in Section 5 (Cell 5.3) to display an
inter-chain confidence summary panel whenever `COMPLEX_TYPE == 'multimer'` and
`query_multimeric == True` for any hit. The panel displays the four derived
metric values with plain-language thresholds drawn from the homodimer
diagnostic notebook's reference thresholds (ipSAE ≥ 0.5 = confident, ipTM ≥
0.5 = confident, pDockQ2 ≥ 0.23 = acceptable, LIS ≥ 0.3 = acceptable).



---

### Section 4 — Run EnzyMM

**Cell 4.1 — Markdown: What EnzyMM is and how it works**

```markdown
## Section 4 — Searching for catalytic motifs

EnzyMM (Enzyme Motif Miner) searches your protein's 3D structure for known
arrangements of catalytic residues — the specific atoms that carry out enzyme
reactions. It does this by comparing the positions and orientations of amino
acid side chains against a library of 6780 templates derived from the
Mechanism and Catalytic Site Atlas, a curated database of experimentally
characterised enzyme mechanisms.

Unlike sequence search tools (BLAST, HMMs), EnzyMM does not look at the
protein sequence at all. It looks only at the three-dimensional arrangement
of atoms. This means it can detect functional similarity between proteins
that share no detectable evolutionary history — a property called convergent
evolution. The same catalytic geometry evolved independently many times across
the tree of life.

A "hit" means that the arrangement of side chain atoms in your protein
geometrically resembles a known active site from a characterised enzyme. This
is a structural hypothesis about your protein's possible function, not a
confirmed activity. The quality of the match — measured by RMSD, orientation
score, and statistical significance — determines how much weight to give it.
```

**Cell 4.2 — Code: Build output paths and run EnzyMM (filtered)**

Creates a temporary working directory. Defines output paths for the filtered
TSV and the PDB match directory. Builds the CLI command as a Python list (see
EnzyMM Execution section). Runs EnzyMM via `subprocess.run()` with
`capture_output=True`. On non-zero return code, prints `stderr` and raises a
`RuntimeError`. On success, prints the number of lines in the output TSV
(excluding the comment header) as a proxy for hit count.

**Cell 4.3 — Code: Run EnzyMM (unfiltered, conditional)**

Executes only when `RUN_COMPARISON=True`. Identical invocation to Cell 4.2
but without `--conservation-cutoff`. Stores result in a separate TSV path.
Prints the unfiltered hit count.

**Cell 4.4 — Code: Parse TSV outputs**

Calls `parse_enzymm_tsv()` on both TSV paths (or just the filtered one if
`RUN_COMPARISON=False`). Stores results as `df_filtered` and `df_unfiltered`.
Handles the zero-hit case: if `df_filtered` is empty, sets `tier1 = None`
and `tier2 = None`, sets a flag `NO_HITS = True`, and prints the zero-hit
message (see Output Parsing section). All subsequent cells check `NO_HITS`
before executing.

**Cell 4.5 — Code: Apply tiering and comparison delta**

Applies the tier assignment logic to `df_filtered`. Sorts Tier 1 by
orientation, RMSD, log E-value. Groups by `template_mcsa_id`. When
`RUN_COMPARISON=True`, computes `low_conf_hits` as the set difference between
unfiltered and filtered hits. Prints a summary:

```
Tier 1 hits (complete + predicted correct) : N
Tier 2 hits (partial or unverified)        : M
Hits found only without pLDDT filter       : K
```

---

### Section 5 — Hit summary table

**Cell 5.1 — Markdown: Reading the results**

```markdown
## Section 5 — Top hits

The table below shows the highest-quality matches between your protein's
structure and known catalytic geometries from the M-CSA. Each row is a match
to one template. The columns most important for interpretation are:

- **EC class** — the enzyme class this template belongs to. This tells you
  what type of reaction the matched geometry is associated with.
- **M-CSA entry** — click the link to read the full curated mechanism.
- **RMSD** — how closely the atom positions match (lower = better, max 2 Å).
- **Orientation** — how closely the side chain directions match (lower =
  better). This is the most sensitive quality indicator.
- **Complete** — whether all sub-patterns of this template matched. A green
  tick means the full active site geometry was found; an orange cross means
  only part of it was detected.

Results in this section passed all quality filters. Partial or lower-
confidence results appear in Section 9.
```

**Cell 5.2 — Code: Render Tier 1 summary table**

Skips if `NO_HITS`. Constructs a display DataFrame from `tier1` with the
following columns: `match_index`, `template_ec` (joined as a string),
`template_mcsa_id` (with an HTML hyperlink to
`https://www.ebi.ac.uk/thornton-srv/m-csa/entry/{id}/`), `rmsd`,
`log_evalue`, `orientation`, `template_effective_size`, `completeness` (as
tick/cross emoji: ✓ / ✗), `predicted_correct` (as tick/cross), and
`matched_residues` (as a clean comma-separated string of
`{3-letter-code}{resnum}` format).

Applies conditional styling via `pandas.io.formats.style`:
- `rmsd` column: background green for < 0.5, yellow for 0.5–1.0, orange
  for > 1.0.
- `orientation` column: background green for < 0.3, yellow for 0.3–0.6,
  no colour for > 0.6.
- `completeness` and `predicted_correct`: green for True, red for False.

Renders via `display(df_display.to_html(escape=False))` to allow the
hyperlinks to render.

Also prints a plain-text fallback summary for non-HTML environments.

**Cell 5.3 — Code: Matched residue pLDDT table**

For each Tier 1 hit, extracts the pLDDT value for each matched residue from
`plddt_by_chain_resnum`. Builds and displays a small table:

```
Hit 1 (M-CSA 45, EC 2.4.1.19):
  ARG 201  — pLDDT 98.9  [Very high >90]
  ASP 203  — pLDDT 98.8  [Very high >90]
  GLU 228  — pLDDT 98.9  [Very high >90]
  HIS 312  — pLDDT 98.7  [Confident 70–90]
  ASP 313  — pLDDT 98.6  [Confident 70–90]
```

For any residue with pLDDT < 70 (only possible in an unfiltered run or if
the cutoff was set below 70), prints a warning box in HTML:

```html
<div style="background:#FFF3E0;border-left:4px solid #EF6C00;padding:8px;
margin:8px 0;">
⚠️ <b>Low-confidence residue detected:</b> {name}{resnum} has pLDDT {value:.1f},
below the reliability threshold of 70. The geometry at this residue is
uncertain; treat this match with caution.
</div>
```

Also detects `template_multimeric=True` + `query_multimeric=False` and prints:

```html
<div style="background:#E8F4FD;border-left:4px solid #42A5F5;padding:8px;
margin:8px 0;">
ℹ️ <b>Multimeric template on monomer model:</b> This template represents an
active site that spans multiple protein subunits. The match was found on a
monomer model; the complete site requires the assembled complex.
</div>
```

Detects metal-dominant templates (`number_of_metal_ligands` template value ≥ 3)
and prints:

```html
<div style="background:#F3E5F5;border-left:4px solid #AB47BC;padding:8px;
margin:8px 0;">
ℹ️ <b>Metal-ligand template:</b> This template consists primarily of
metal-coordinating residues. It may indicate a metal binding site rather than
a complete catalytic active site. AlphaFold models do not include metal ions,
so the geometry shown is the protein-side coordination scaffold only.
</div>
```

---

### Section 6 — Biological interpretation

**Cell 6.1 — Markdown: Interpretation framework**

```markdown
## Section 6 — Biological interpretation

The following analysis interprets the EnzyMM results in the context of what
is already known about your protein. The interpretation adapts based on
whether your protein has a confirmed enzyme annotation, a reviewed entry with
unknown function, or no prior characterisation.
```

**Cell 6.2 — Code: Auto-generate narrative**

This cell constructs and displays the auto-generated narrative block. It
consists of three paragraphs, each generated from conditional logic.

**Paragraph 1 — Protein context (adapts by PROTEIN_MODE):**

```python
if PROTEIN_MODE == 'annotated_enzyme':
    p1 = (
        f"{protein_name} ({UNIPROT_ACCESSION}) is annotated in UniProt as a "
        f"characterised enzyme with EC classification. EnzyMM found "
        f"{len(tier1)} Tier 1 catalytic geometry match(es). "
        f"We compare the matched residues against the known active site "
        f"positions from UniProt below."
    )
elif PROTEIN_MODE == 'annotated_no_function':
    p1 = (
        f"{protein_name} ({UNIPROT_ACCESSION}) has a reviewed UniProt entry "
        f"but no enzyme classification. EnzyMM found {len(tier1)} potential "
        f"catalytic arrangement(s). These are functional hypotheses, not "
        f"confirmed activities, and should be treated as leads for experimental "
        f"follow-up."
    )
else:
    p1 = (
        f"{UNIPROT_ACCESSION} has no reviewed annotation in UniProt. "
        f"EnzyMM found {len(tier1)} geometric match(es) to known catalytic "
        f"templates. All results should be treated as candidate hypotheses "
        f"requiring experimental validation before any functional conclusion "
        f"is drawn."
    )
```

**Paragraph 2 — Hit interpretation (EC and residue analysis):**

Logic, in order:

1. If all Tier 1 hits share at least one CATH domain identifier:
   "All matches belong to the same structural superfamily (CATH {cath_id}),
   indicating that the detected active site geometry is conserved across
   this enzyme fold."

2. If EC classes are all the same first digit:
   "The matched enzymes all belong to EC class {digit} ({class_name}),
   strengthening the functional inference."
   EC class names: 1=Oxidoreductases, 2=Transferases, 3=Hydrolases,
   4=Lyases, 5=Isomerases, 6=Ligases, 7=Translocases.

3. If EC classes diverge at the first digit:
   "The matched templates span different enzyme classes ({classes}).
   Evaluate the RMSD and orientation scores for each hit independently
   and consult the M-CSA entries before drawing functional conclusions."

4. If shared_resnums is non-empty (residues matched by > 1 hit):
   "The following residues are matched by multiple hits, providing
   convergent geometric evidence for this catalytic core: {residue_list}.
   This convergence strengthens the inference that these residues form a
   genuine active site arrangement."

5. If `PROTEIN_MODE == 'annotated_enzyme'` and
   `uniprot_active_site_residues` is non-empty:
   Compute the intersection of matched resnums with UniProt annotated
   positions. If overlap ≥ 1:
   "EnzyMM matched residues {matched} include position(s) {overlap} which
   are annotated as active site residues in UniProt. This is consistent
   with the known catalytic mechanism."
   If overlap == 0:
   "None of the EnzyMM matched residues correspond to UniProt-annotated
   active site positions ({uniprot_positions}). This may indicate a
   secondary binding site, a cryptic active site, or a geometric
   coincidence. Inspect the 3D visualisation carefully."

**Paragraph 3 — Confidence statement:**

```python
all_plddt = [plddt_by_chain_resnum[res['chain']][res['resnum']]
             for _, row in tier1.iterrows()
             for res in row['matched_residues_parsed']
             if res['resnum'] in plddt_by_chain_resnum.get(res['chain'], {})]

if all_plddt:
    min_plddt = min(all_plddt)
    if min_plddt >= 90:
        p3 = ("All matched residues are in the highest confidence band "
              "(pLDDT > 90). The predicted geometry at these positions is "
              "reliable and the match quality can be trusted.")
    elif min_plddt >= 70:
        p3 = (f"Most matched residues are in the confident band (pLDDT 70–90). "
              f"The lowest pLDDT among matched residues is {min_plddt:.1f}. "
              f"Backbone geometry is likely correct at these positions; side "
              f"chain rotamers may deviate from the true conformation.")
    else:
        p3 = (f"One or more matched residues fall below pLDDT 70 "
              f"(minimum: {min_plddt:.1f}). The geometry at these positions "
              f"is uncertain. This match should be treated with caution until "
              f"validated experimentally or by a higher-confidence model.")
```

Display all three paragraphs as styled HTML in an output box. Each paragraph
is separated by a horizontal rule.

**Cell 6.3 — Code: Multi-EC summary table**

When Tier 1 contains more than one unique EC class at the first digit level,
renders a small summary table:

| EC class | Reaction type | Hits | Top RMSD | Shared with other hits? |
|---|---|---|---|---|
| 2 (Transferase) | Transfer of functional groups | 1 | 0.55 | 4/5 residues shared |
| 3 (Hydrolase) | Hydrolysis reactions | 1 | 0.47 | 4/5 residues shared |

Printed beneath a heading: "Multiple EC classes detected — summary."

---

### Section 7 — pLDDT sequence strip

**Cell 7.1 — Markdown: Reading the confidence plot**

```markdown
## Section 7 — AlphaFold confidence at matched residues

The plot below shows AlphaFold's per-residue confidence score (pLDDT) across
your full protein sequence. Confidence is shown using the standard AlphaFold
colour scheme:

- **Dark blue** (pLDDT > 90): very high confidence — backbone and side chains
  are predicted reliably.
- **Light blue** (70–90): confident — backbone is reliable; side chains may
  deviate slightly.
- **Yellow** (50–70): low confidence — treat structural details with caution.
- **Orange** (< 50): very low confidence — the structure in this region is
  largely unreliable.

Matched catalytic residues are marked with vertical lines. A match in the dark
or light blue region is geometrically trustworthy. A match in yellow or orange
should be treated as provisional.
```

**Cell 7.2 — Code: Render the pLDDT sequence strip**

Skips if `NO_HITS`. Extracts per-residue pLDDT for all residues in chain A
(or the primary chain) from `plddt_by_chain_resnum`. Builds a numpy array of
(resnum, plddt) pairs sorted by resnum.

Creates a matplotlib figure (width=14, height=4). Plots a filled step chart
or scatter plot with each point coloured by its band. Draws horizontal dashed
grey lines at y=70 and y=90. Sets y-axis limits to (0, 100), labels axes
"Residue position" and "pLDDT".

For each unique matched residue in Tier 1 (all hits combined), draws a
vertical dashed line at x=resnum from y=0 to the pLDDT at that residue.
Places a text label above the line: `{name}{resnum}` with the pLDDT value
below it. Staggers label heights for adjacent residues (alternating 75 and 85
in y-axis coordinates) to prevent overlap.

Adds a legend with four coloured patches for the four bands, plus a dashed
line marker for "Matched catalytic residue."

Saves the figure at 150 DPI and calls `plt.show()`.

Beneath the figure, renders the residue pLDDT table (Cell 5.3 content) again
for easy reference without scrolling.

---

### Section 8 — 3D visualisation

**Cell 8.1 — Markdown: Viewing the structure**

```markdown
## Section 8 — 3D structure visualisation

The interactive views below show your protein structure coloured by AlphaFold
confidence, with the matched catalytic residues highlighted. You can rotate,
zoom, and pan the structure using your mouse or trackpad.

**View 1** colours the entire structure by pLDDT (same colour scheme as the
sequence plot above) and highlights the matched catalytic residues as spheres
and sticks in magenta.

**View 2** (shown when multiple hits were found) uses a grey structure with
hit-specific residue colours. Residues matched by more than one hit are shown
in green — these are the most robustly detected part of the active site.
```

**Cell 8.2 — Code: View 1 — pLDDT confidence map with matched residues**

Skips if `NO_HITS` or `struct_url is None` (local file mode without a URL).
Prints a warning if struct_url is None: "3D views require a URL and are not
available for local file uploads in this version."

Implements the MolViewSpec View 1 as described in the Visualisations section.
Applies pLDDT colouring residue by residue for chain(s) in the structure,
then adds ball-and-stick representation for Tier 1 hit 1 matched residues in
magenta (`#E040FB`). Calls `show_mol_view()`.

**Cell 8.3 — Code: View 2 — multi-hit residue overlap (conditional)**

Executes only when `len(tier1) > 1`. Implements MolViewSpec View 2 as
described in the Visualisations section. Grey cartoon for the full structure,
ball-and-stick coloured by hit membership and shared/unique status. Calls
`show_mol_view()`.

---

### Section 9 — Secondary hits

**Cell 9.1 — Markdown: Partial and unverified matches**

```markdown
## Section 9 — Secondary hits

The following matches did not pass all quality filters. They may represent:

- **Partial active site detections:** only a subset of a known template's
  residues were found (completeness = False). This is not always noise — it
  can indicate that part of an active site is geometrically present even if
  the full arrangement is not.
- **Matches not predicted correct:** EnzyMM's internal classifier judged these
  as lower confidence. The raw geometric metrics are still shown.
- **Low-confidence region matches:** hits found only in the run without the
  pLDDT filter. The matched residues are in a region where AlphaFold's
  prediction is uncertain.

Consult the RMSD, orientation, and template size for each before drawing
any conclusions.
```

**Cell 9.2 — Code: Render Tier 2 table**

Skips if `tier2` is empty and `low_conf_hits` is empty. Renders the same
styled table format as Section 5 but without pLDDT annotation and with the
caveat markdown above. `low_conf_hits` rows are appended to the table with
a column `flag` set to `'Low-confidence region'`. Tier 2 rows have
`flag = ''`. The `flag` column is displayed in orange text.

---

### Section 10 — M-CSA cross-links

**Cell 10.1 — Markdown + Code: M-CSA entry links**

For each unique `template_mcsa_id` in Tier 1, displays:

```
M-CSA Entry {id}: {enzyme_name}
EC: {ec_numbers}
Organism: {organism}
→ https://www.ebi.ac.uk/thornton-srv/m-csa/entry/{id}/
```

Where `enzyme_name` and organism are taken from the `template_uniprot_id`
column (or looked up from a pre-fetched UniProt call for the template UniProt
accession — to be decided during engineering depending on API availability).
If the lookup is too slow or unreliable, display only the ID and link.

Renders as HTML with each entry in a styled card-like div.

---

### Section 11 — Caveats and limitations

**Cell 11.1 — Markdown: Fixed caveats block**

```markdown
## Section 11 — Caveats and limitations

Please read these before drawing conclusions from EnzyMM results on an
AlphaFold model.

**1. A match is a hypothesis, not a confirmed activity.**
EnzyMM identifies geometric similarity between your protein's structure and
a known catalytic arrangement. This is a structural hypothesis. Only
experimental characterisation — mutagenesis, activity assays, structural
studies of the protein with substrate — can confirm enzymatic function.

**2. The M-CSA does not cover all enzymes.**
The Mechanism and Catalytic Site Atlas currently describes 1003 enzyme
families. Thousands of known enzyme families are not yet curated. A no-hit
result does not mean the protein is not an enzyme; it may mean its mechanism
is outside current M-CSA coverage.

**3. AlphaFold predicts one conformation.**
AlphaFold2 predicts a single structural model, typically corresponding to
the apo (substrate-free) state. Catalytic residues often move upon substrate
binding (induced fit). Some active sites are disordered in the apo state and
are only organised in the presence of substrate. EnzyMM may miss these sites
or match them only partially.

**4. Side chain rotamers may be inaccurate.**
Even in high-confidence regions, AlphaFold predicts one rotamer per residue.
Catalytic residues frequently adopt multiple rotamer states. A geometric
mismatch between the predicted rotamer and the template orientation inflates
the orientation score. A match with moderate orientation score in a
high-pLDDT region may still be biologically meaningful.

**5. Metal ions and cofactors are absent.**
AlphaFold does not model metal ions, heme groups, flavins, or other
cofactors. For metalloenzymes, the protein-side coordination geometry may
match a template while the full catalytic mechanism requires a metal that
is not present in the model. Hits to metal-ligand templates (flagged in
Section 5) should be interpreted with this in mind.

**6. Multimeric active sites require multimeric models.**
Some enzymes form their active sites at the interface between two or more
subunits. EnzyMM detects these with multimeric templates. If your query is
a monomer model but the template is multimeric (flagged in Section 5), the
detected geometry is incomplete. A ColabFold or AF3 multimer prediction would
be needed to assess the full site.

**7. Performance benchmarks.**
EnzyMM was validated on the AlphaFold2-predicted human proteome: 42% of
annotated enzymes were detected, 6.4% of non-enzymes received at least one
hit, and 83.7% of matched enzymes with UniProt active site annotations were
matched at the correct residues. False positives and missed enzymes both
occur; results should inform, not replace, expert judgement.
```

---

## User-Facing Explanations

### What EnzyMM is and how it works

> EnzyMM — the Enzyme Motif Miner — is a tool that searches a protein's
> three-dimensional structure for known arrangements of catalytic residues:
> the specific atoms responsible for carrying out enzymatic reactions. It does
> this by comparing the positions and orientations of amino acid side chains
> in your protein against a library of 6780 templates derived from the
> Mechanism and Catalytic Site Atlas, a curated database of experimentally
> characterised enzyme mechanisms covering 895 enzyme classes. The search is
> purely geometric — it does not look at the protein sequence at all — which
> means it can detect functional similarities between proteins that share no
> detectable evolutionary history. The same catalytic geometry has evolved
> independently many times across the tree of life; EnzyMM is designed to
> find these patterns wherever they occur, including in proteins predicted
> by AlphaFold.

### What a catalytic motif match means (and doesn't mean) in a predicted structure

> Finding a catalytic motif match means that the three-dimensional arrangement
> of side chain atoms in your predicted protein structure geometrically
> resembles a known active site from a characterised enzyme, as recorded in
> the M-CSA. This is a structural hypothesis about your protein's possible
> function, not a confirmed catalytic activity. The strength of the hypothesis
> depends on the quality of the match (RMSD, orientation score), the size and
> specificity of the matched template, and the AlphaFold confidence at the
> matched residues. A strong match — tight geometry, high pLDDT at the matched
> residues, the full template pattern found rather than just a subset — is
> meaningful positive evidence that deserves experimental follow-up. A weak or
> partial match is a lead for investigation. Neither outcome is definitive
> without experimental validation.

### How to interpret the confidence overlay

> The colour scheme on the sequence strip and 3D viewer follows the standard
> AlphaFold confidence scale. Dark blue (pLDDT > 90) means AlphaFold is very
> confident about both the backbone and the likely side chain position at that
> residue — geometry in this region can be trusted. Light blue (70–90) means
> the backbone is reliable but the exact side chain rotamer may deviate from
> the true structure. Yellow (50–70) means the local structure is uncertain;
> treat positional details with caution. Orange (below 50) means the region
> is largely unstructured or the model is highly uncertain — structural
> conclusions in this region are not reliable. When catalytic residues matched
> by EnzyMM fall in the dark or light blue region, the geometric match is
> trustworthy. When they fall in yellow or orange, the match may reflect
> uncertain geometry rather than a genuine active site, and should be
> interpreted accordingly.

### The caveat section on AlphaFold limitations

> AlphaFold2 is a remarkable tool but it has specific limitations that matter
> when interpreting EnzyMM results. First, it predicts a single conformation
> — typically the apo state without substrate — and some active sites only
> adopt their functional geometry when substrate is bound. Second, it does not
> model metal ions or cofactors, which are essential to the catalytic mechanism
> of many enzyme classes. Third, its side chain predictions become less
> accurate in flexible loops and near protein surfaces, even when the backbone
> is correct. Finally, the pLDDT score measures AlphaFold's confidence in its
> own prediction, not the accuracy of the prediction relative to experiment —
> a high pLDDT means AlphaFold is certain about what it predicted, not
> necessarily that the prediction is correct. In practice, for well-folded
> globular proteins with pLDDT > 70 throughout, AlphaFold backbone accuracy
> is high and catalytic residue positions are usually reliable. For disordered
> regions, termini, and surface loops, caution is warranted.

---

## Flywheel Integration

Flywheel integration (opt-in submission of results to the InsightFold
pre-computed layer) is deferred to a post-v1 sprint. No code for result
submission is included in this version. The notebook does not collect user
data. See Out of Scope.

---

## Out of Scope for v1

| Feature | Rationale for deferral |
|---|---|
| Flywheel opt-in result submission | No backend infrastructure in place; usage signal collection deferred to a later sprint once notebook adoption is established |
| M-CSA mechanism text (full reaction description) | M-CSA does not expose a public per-entry text API; hyperlinks to the M-CSA web interface are used instead |
| Custom user-supplied template libraries | EnzyMM supports `--template-dir`; not exposed in v1 to keep the interface simple |
| Batch mode (multiple accessions in one run) | Out of scope for this use case; batch users should use the EnzyMM CLI directly |
| PDB format input | AFDB native format is mmCIF; PDB format is incomplete for modern structures and omits important metadata |
| Foldseek or global fold comparison integration | Different question from catalytic geometry; a separate notebook |
| Population-scale pre-computation | The v1 notebook answers "what templates does my protein match?" not "what proteins in AFDB match this template?" Pre-computation deferred until usage signal justifies it |
| Rosetta or MD validation of matched geometries | Out of scope for a lightweight notebook; relevant for downstream experimental design |

---

## Open Questions Requiring Resolution Before Build

### 1. `--skip-smaller-hits` interaction with partial site detection

The `--skip-smaller-hits` flag suppresses smaller template matches when a
match to a larger template from the same structural region has already been
found. The exact interaction with the cluster/completeness system when a large
hit exists alongside a small partial hit from a **different** M-CSA entry was
not tested during spec development. It is possible the flag suppresses
genuinely independent hits from different enzyme families if their matched
residues spatially overlap in the query structure.

The engineering session should run a test case on a protein that produces both
large and small template hits — with and without `--skip-smaller-hits` — to
confirm the flag does not silently discard Tier 1 hits from distinct M-CSA
entries. If suppression of cross-entry hits is confirmed, expose a second
boolean parameter `SKIP_SMALLER_HITS` (default True) with an explicit warning
in the user parameter section explaining the trade-off: True reduces noise but
may miss independent partial sites; False returns the full hit set at the cost
of more manual filtering.
