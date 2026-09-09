# InsightFold — CLAUDE.md

## Project Overview

InsightFold analyses AlphaFold DB (AFDB) structure predictions for protein–protein
interactions. The primary use case is **dimer** interface quality assessment (homodimer *and*
heterodimer, automatically, with no user-facing switch) using AFDB REST APIs and NumPy
implementations of published scoring functions.

**This file is a map, not a reference.** The numbers and formulas live in the spec pack and in
`src/insightfold/complex_interface_utils.py`. Earlier revisions of this file carried its own
copies of both, they drifted, and every copy was wrong by 2026-09. Read the authorities below
before writing a formula or a threshold; do not transcribe them back into this file.

| Authority | Owns |
|-----------|------|
| `specs/homodimer_diagnostic/formula-reference.md` | **All seven formulas**, line-cited to `ipsae.py` v4 |
| `specs/homodimer_diagnostic/threshold-reference.md` | **All thresholds, bands and provenance**, plus the AFDB joint criterion |
| `specs/homodimer_diagnostic/rework-plan.md` | Locked decisions D1–D10, task history, every measured finding |
| `specs/homodimer_diagnostic/validation-report.md` | The 2026-09-08 end-to-end validation, 231 numerical comparisons |
| `specs/homodimer_diagnostic/fixture-manifest.md` | FX-001 … FX-010, what each one is for |
| `src/insightfold/complex_interface_utils.py` | The executable transcription of all of the above |

`threshold-reference.md` is transcribed *into* the module. Edit the document first and
re-transcribe; never the other way round.

---

## Directory Layout

```
src/insightfold/
  complex_interface_utils.py   ← THE module. AFDB access, parsing, interface
                                 detection, scoring, thresholds, plots, 3D views.
  variants/                    ← variant mapping submodule (unrelated pipeline)
notebooks/
  homodimer_diagnostic.ipynb   ← reference implementation; orchestration only
  analysis_template.ipynb      ← parametrised scaffold for new analysis types
specs/homodimer_diagnostic/    ← the spec pack listed above
  references/                  ← GITIGNORED. Publisher PDFs + ipsae_v4.py. Never commit.
agents/, agent-skills/         ← LLM advisory role definitions
```

`src/insightfold/interface.py` was **deleted 2026-09-08** (R095). Everything it did is in
`complex_interface_utils.py`, which additionally fixes a latent parser bug it had: it treated a
blank line inside the `_atom_site` loop as a terminator and silently dropped every atom after
it. If you find a reference to `interface.py`, it is stale.

### Module layout — `complex_interface_utils.py`

One flat module (D1), named by **domain** rather than by notebook (D8), living in
`src/insightfold/` (D7). Sections, in file order:

| Section | Contents |
|---------|----------|
| Constants | cutoffs, palettes, `AFDB_NEWS_URL`, `AFDB_RELEASE_SCALE` |
| AFDB access | `fetch_afdb_metadata`, `AFDBPrediction`, `download_*`, local-file mode (R025), assembly detection (R023), `format_metadata_report` (R020) |
| Structure parsing | `parse_mmcif_atoms`, `extract_chain_coords`, `parse_structure`, `ChainCoords` |
| PAE / pLDDT parsing | `parse_pae`, `parse_plddt`, `PAEMatrix`, `PLDDTScores`, `ChainPairPAE`, `verify_chain_identity` / `verify_document_agreement` / `verify_chain_lengths`, chain labelling (R021) |
| Interface detection | `detect_interface`, `InterfaceContacts` |
| Scoring primitives | `d0_scalar`, `d0_array`, `ptm_func` |
| Score functions | `compute_iptm_d0chn`, `compute_ipsae`, `compute_pdockq`, `compute_pdockq2`, `compute_lis`, directional transparency (R008) |
| Thresholds | `THRESHOLDS`, `Threshold`, `traffic_light`, `afdb_high_confidence`, `summarise_scores` |
| Plotting | one `plot_*` per figure, each returning a bare `Figure` |
| MolViewSpec views | one builder per view, plus `show_mol_view` |

**D8 — naming and promotion.** Modules are named by domain so several notebooks can share one.
A function is promoted out of a domain module into a shared one when its **second** consumer
appears, not in advance.

**D4 — dimer-scoped, chain-pair-generic.** Every API takes an *explicit ordered chain pair* (a
directional PAE block plus the two chain lengths, or a chain-pair object) rather than assuming
"A and B". PAE is asymmetric, so "the chain pair" always means an ordered pair. An N-chain
generalisation is then an extension over all ordered pairs, not a rewrite.

**D10 — deletions land last.** Redundant code is identified when found but deleted only after
all dependent work is verified, so an early deletion cannot masquerade as a later bug.

**Two names that must never be confused** (they were, before R015):
`n_interface_residues` (residues with at least one contact) vs `n_contact_pairs` (CB–CB pairs
within cutoff). On FX-001 they are 86 and 116. **pDockQ needs the pair count.**

---

## Dependencies

Allowed in the notebook and the module: `numpy`, `matplotlib`, `seaborn`, `requests`,
`molviewspec`, `ipywidgets`. `molviewspec` is imported **lazily**, inside the function that
needs it.

Prohibited in `complex_interface_utils.py` and in the dimer notebooks: `biopython` / `Bio`,
`torch`, `torch-geometric`, `gemmi`, `pandas`, `scipy`, `plotly`. The notebook must run on
Colab free tier in < 60 s of wall-clock install time.

### The `pyproject.toml` contradiction — resolved, read this before "fixing" either side

`pyproject.toml` declares `biopython`, `gemmi`, `scipy`, `networkx`, `plotly`, `pandas` and
`jupyter` as hard dependencies of the `insightfold` distribution. That is **not** a
contradiction of the rule above once the scope of each is stated, and both statements are
kept deliberately:

- The **prohibition is per-consumer**, not repo-wide. It binds `complex_interface_utils.py`,
  `homodimer_diagnostic.ipynb` and `analysis_template.ipynb`. Those three import none of the
  prohibited packages, verified by grep.
- The **`pyproject.toml` dependency set belongs to the other pipelines** in this repo:
  `src/insightfold/variants/` and `src/pdbe_interfaces/` import `pandas` and `scipy`, and
  `notebooks/protein_model_chem.ipynb` is the only file in the repo that imports `Bio`.
  `networkx` is used by that same notebook; `plotly` is imported by nothing in the repo and is
  a candidate for removal.
- **This is precisely why D9 forbids `pip install`** for the Colab bootstrap: installing the
  package would resolve `pyproject.toml` and drag the whole set in, breaking both the
  dependency rule and the 60 s budget. The bootstrap clones and extends `sys.path` instead, so
  the heavy dependencies are never resolved on the notebook path.

If the dimer work is ever split into its own distribution, move these to an optional extra.
Until then: **do not add a prohibited package to the module or the dimer notebooks, and do not
delete one from `pyproject.toml` without checking `src/insightfold/variants/` and
`notebooks/protein_model_chem.ipynb` first.**

### Colab bootstrap (D9)

`git clone --depth 1` plus `sys.path.insert(0, root / 'src')`. **Never `pip install`.** The
one bootstrap cell works unchanged locally and on Colab: `find_repo_root()` walks up for a
directory holding both `pyproject.toml` and `src/`; on Colab a missing or stale clone is
fetched and hard-reset into `/content/InsightFold`. `git reset --hard` is reachable only for
that disposable scratch directory, never for a checkout the user is working in. The repo is
public, so there is no token, no auth header and no `getpass` (which would block forever in a
Run-all notebook).

> **⚠️ TODO(merge) — ACTION REQUIRED WHEN THIS BRANCH MERGES.**
> The bootstrap cell of `notebooks/homodimer_diagnostic.ipynb` pins
> `REPO_BRANCH = 'homodimer-notebook-rework'`, because `complex_interface_utils.py` exists
> only on that branch. **Once it merges to `main`, flip `REPO_BRANCH` to `'main'` or, better,
> to a release tag** so a Colab run a year later reproduces rather than picking up drift.
> Grep for `TODO(merge)` before releasing the notebook. Until this is done, **every Colab run
> clones a feature branch**, and no Colab result predicts post-merge behaviour.

**Colab is currently unverified.** Every validation run had `IN_COLAB = False`, so the clone,
stale-clone-refresh and `pip install molviewspec` paths have never executed and the 60 s
budget is unmeasured. Any deliverable that runs in two environments must be exercised in both
or the untested one must be declared untested.

---

## Skill: Complex Confidence Scoring

**Do not reimplement these. Import them.**

```python
import sys; sys.path.insert(0, 'src')
from insightfold.complex_interface_utils import (
    compute_iptm_d0chn, compute_ipsae, compute_pdockq, compute_pdockq2, compute_lis,
    THRESHOLDS, traffic_light, afdb_high_confidence,
)
```

**Ground truth (D2):** `DunbrackLab/IPSAE` `ipsae.py` **v4** (3 Jan 2026, md5
`a48df7adc64afa36b4b475429b0a8aff`, sha256 `10cf9b08…`), stored gitignored at
`specs/homodimer_diagnostic/references/ipsae_v4.py`. Tolerance ±0.001. Publications supply the
explanatory prose and threshold rationale only. **Every formula is transcribed and line-cited
in `formula-reference.md`; read that file, not this section, before touching a score.**

### The seven values, and the trap in each

`ipsae.py` emits three rows per chain pair: `A→B`, `B→A`, and a `max` row. **The number a user
quotes is the `max` row.** PAE is asymmetric (`pae_AB != pae_BA.T`), so direction matters
everywhere.

| Value | One line | The trap |
|-------|----------|----------|
| `ipTM_d0chn` | mean `ptm(PAE, d0chn)` over **all** partner residues, no PAE cutoff; max over directions | It is **not** AlphaFold's own ipTM. Label it `ipTM_d0chn` |
| `ipSAE_d0chn` | same, restricted to `PAE < 10`; `d0 = d0_scalar(nA + nB)` | provably the **most permissive** of the three ipSAE variants |
| `ipSAE_d0res` | per row, `d0` from that row's own valid-pair count, using **`d0_array`** | uses the *array* `d0` helper, which differs from the scalar one at `L == 27` |
| `ipSAE_d0dom` | `d0` from the residues participating in any sub-cutoff pair, **per direction** | rows-with-any **plus cols-with-any of that same direction's block**; recompute for `B→A` |
| `pDockQ` | `x = mean_plddt * log10(npairs)` through Bryant's sigmoid | `npairs` is the **contact-pair count**, not the interface residue count. Symmetric: no max |
| `pDockQ2` | `x = mean_plddt * mean_ptm(PAE_contacts, d0=10)` through Zhu's sigmoid | directional (`mean_ptm` reads one block); `mean_plddt` is an **unweighted mean over unique interface residues** |
| `LIS` | per direction `mean((12 - PAE)/12)` over `PAE < 12` | combined as the **mean** of the two directions, not the max. Cutoff 12 is hardcoded and independent of the PAE cutoff |

`ipSAE_d0chn >= ipSAE_d0dom >= ipSAE_d0res` is a **theorem**, not an observation (`d0` is
monotone in `L`, `ptm` increases in `d0`, and `n0res_i <= n0dom <= nA + nB`). Consequence: the
three variants are one measurement seen three ways and must not be presented as three
independent confirmations, and a threshold set that gives `d0chn` a *lower* bar has the bias
backwards.

### The two `d0` helpers are not interchangeable

```python
d0_scalar(L)   # used for d0chn (ipTM_d0chn, ipSAE_d0chn) and d0dom
d0_array(L)    # used for d0res ONLY, per residue
```

They are bit-exact for every integer `L` except **`L == 27`**, where the scalar returns `1.0`
and the array form returns `1.038891`. On a residue with exactly 27 valid inter-chain pairs
this moves the score by 0.019 — 19× the tolerance. `ptm_func(x, d0) = 1 / (1 + (x/d0)**2)`.

### PAE quadrant extraction

Chain lengths come from the PAE JSON `chains` field, but **never index the `chains` array
positionally**: resolve chains by `label_asym_id` against the structure and the pLDDT document
(`verify_chain_identity`), and fail loudly on a mismatch. Neither the documents' `chains`
arrays nor the prediction endpoint's entry list carries a trustworthy order, and chain length
is no tie-breaker because a homodimer's chains are equal by definition. A positional guess
mis-slices every quadrant into plausible-looking wrong scores rather than an error.

```python
pae_AB = pae_matrix[:nA, nA:nA+nB]   # rows = A, cols = B  (the UPPER-RIGHT quadrant)
pae_BA = pae_matrix[nA:nA+nB, :nA]
pae_AA = pae_matrix[:nA, :nA]
pae_BB = pae_matrix[nA:nA+nB, nA:nA+nB]
```

### Thresholds and the classifier

**`specs/homodimer_diagnostic/threshold-reference.md` is the single source of truth.** The
module's `THRESHOLDS` is its transcription. Reproduced here only so a stale copy elsewhere is
recognisable as stale:

| Key | Green | Amber | Green provenance | Bands |
|-----|-------|-------|------------------|-------|
| `ipsae_d0res` | 0.70 | 0.60 | PUBLISHED — AFDB 2026, p. 12 | **4, published** |
| `ipsae_d0chn` | 0.70 | 0.60 | DERIVED (transferred from `d0res`) | 3 |
| `ipsae_d0dom` | 0.70 | 0.60 | DERIVED (transferred from `d0res`) | 3 |
| `iptm_d0chn` | 0.70 | 0.30 | DERIVED — Dunbrack 2025, p. 14 | 3 |
| `pdockq` | 0.23 | 0.12 | DERIVED — Bryant 2022, pp. 1, 3, 10 | 3 |
| `pdockq2` | 0.23 | 0.10 | PUBLISHED — Zhu 2023 p. 6 **and** AFDB 2026 p. 5 | 3 |
| `lis` | 0.21 | 0.10 | PUBLISHED — Kim 2024, p. 30 | 3 |

Semantics: `value >= green` → green; `green > value >= amber` → amber; else red. **Amber is a
stored number and must never be recomputed as `green / 2`.** The two amber values with no
literature basis at all (`pdockq2` 0.10, `lis` 0.10) are judgement calls and must be presented
as such; `ipsae_d0res` carries AFDB's four *named* bands, so its `traffic_light()` labels are
`VERY HIGH-CONFIDENCE` / `CONFIDENT` / `LOW-CONFIDENCE` / `BELOW AFDB THRESHOLD` rather than
HIGH / MODERATE / LOW.

**The AFDB classifier is a conjunction, not a single threshold:**

> `ipSAE_d0res >= 0.6` **AND** `pDockQ2 >= 0.23` → PASS / FAIL — `afdb_high_confidence()`

Precision 0.924 (FPR 0.043) on AFDB's homodimer benchmark, 0.958 (FPR 0.004) on the
heterodimer one (AFDB 2026, p. 5). Three things about it that prose must not smooth over:

1. **0.6 is AFDB's threshold, not Dunbrack's.** Dunbrack 2025 states no ipSAE cutoff anywhere,
   for any variant, at any PAE cutoff — verified by full-text search. Cite the *metric* to
   Dunbrack and the *threshold* to Han, Tsenkov, Venanzi et al. 2026. The metric belongs to
   whoever published it whatever cutoff is later applied; the same split applies to pDockQ
   (Bryant), pDockQ2 (Zhu) and LIS (Kim).
2. **0.6 is conservative, not optimal.** AFDB adopts it by citation as "community-established"
   and derives only the validation. Its own MCC-optimal cutoffs are **0.104** (homodimer) and
   **0.520** (heterodimer). A FAIL means "not selected for high-confidence release", **not**
   "not an interaction"; below-threshold dimers are published with their scores.
3. **On an AFDB accession the traffic light barely discriminates.** The entry exists because it
   already passed the 0.6 filter, so red essentially never appears. The full range is only
   reachable through the local-file path. Say so.

The homodimer caveat belongs beside the traffic light: AlphaFold's false-positive rate is
markedly worse for homodimers than heterodimers (true-positive rate at 1% FPR drops 63% → 18%),
because proteins that do not homodimerise in vivo often have homologs that do.

### Numbers that are frequently misquoted

`pDockQ2`'s sigmoid is bounded on `[0.005, 1.315]` and **can exceed 1.0** — never present it as
a probability or a percentage. Kim 2024's `0.911` / `0.891` are **AUC values, not thresholds**.
`~31 million` is the count of *candidate complexes compiled*, not dimers analysed
(19,148,379 homo + 7,561,477 hetero); 1,735,475 homo and 79,392 hetero are high-confidence,
1.81 M assemblies in total. `AFDB_RELEASE_SCALE` in the module carries all of these with their
page citations so prose never has to invent one.

---

## AFDB REST API Contracts

### Metadata endpoint

```
GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}
```

Accession format `AF-{id}` (e.g. `AF-0000000065889468`), or a bare UniProt accession for a
monomer.

> **The response array holds ONE ENTRY PER CHAIN, in non-deterministic order.**
> `result[0]` is **not safe**. The same accession returns `['A','B']` on one call and
> `['B','A']` on the next — measured over eight consecutive calls per accession. For a
> heterodimer that makes the reported UniProt ID, organism, protein name and monomer length
> flip between chains on consecutive runs of an unchanged notebook; on a homodimer it is
> invisible. Key by `chainId` (`AFDBPrediction` sorts `entries` by `chainId` in
> `__post_init__`, so arrival order is discarded before any caller can observe it), and use
> `shared_field()` — which *raises* when entries disagree — rather than trusting entry 0.

| Field | Type | Used for | Note |
|-------|------|----------|------|
| `chainId` | `'A'` / `'B'` | **the key for every per-chain field** | |
| `cifUrl` | URL | mmCIF download (coordinate parsing) | identical across entries; agreement is asserted |
| `bcifUrl` | URL | BinaryCIF for MolViewSpec only | may be present but **empty**; treat empty as absent |
| `paeDocUrl` | URL | PAE JSON download | served gzipped — use `requests`, not raw `urllib` |
| `plddtDocUrl` | URL | pLDDT JSON download | served gzipped |
| `isComplex` | bool | complex vs monomer | present on monomers too, as `false` |
| `assemblyType` | `'Homo'` / `'Hetero'` | homodimer vs heterodimer statement | absent on monomers |
| `oligomericState` | `'dimer'` | early refusal gate | absent on monomers |
| `oligomericStateDescription` | string | display, e.g. `'Heterodimer'` | **optional; absent on both fixtures.** Read defensively |
| `complexComposition` | list | stoichiometry, e.g. `[{identifierType, identifier, stoichiometry: 2}]` | absent on monomers |
| `uniprotAccession`, `uniprotId`, `uniprotDescription` | string | display | per chain |
| `sequence`, `sequenceStart`, `sequenceEnd` | | monomer length | per chain |
| `organismScientificName` | string | display | per chain |
| `gene` | string | display | **not `geneNames`** |
| `latestVersion` | int | display | **not `modelVersion`** |
| `globalMetricValue` | float | global pLDDT | can differ slightly between chain records |

**`geneNames` and `modelVersion` do not exist.** Earlier revisions of this file documented
those names; both printed `N/A` on every run for months as a result. The live service sends
`gene` and `latestVersion`.

**Measured 2026-09-08: every AFDB complex is a dimer.** On the search endpoint `isComplex:true`
and `oligomericState:dimer` both return **2,010,763** (Homo 1,930,523 + Hetero 80,240, summing
exactly); `trimer` and `tetramer` return nothing. So more than two chains is reachable only via
a local mmCIF. D4's future N-chain work is a local-file feature, not an AFDB one.

### Search endpoint — undocumented elsewhere, and a validation reference

```
GET https://alphafold.ebi.ac.uk/api/search?q=<solr query>&type=main&rows=<n>
    q=modelEntityId:AF-0000000065889468                     → one entry
    q=assemblyType:Hetero AND oligomericState:dimer         → numFound 80240
```

Used to source fixtures (R024), and — more importantly — **it publishes AFDB's own production
IPSAE-derived scores**, which are a second independent reference alongside running `ipsae.py`:

```
complexPredictionAccuracy_ipsae_{AB,BA}          complexPredictionAccuracy_pDockQ2_{AB,BA}
complexPredictionAccuracy_ipsae_d0chn_{AB,BA}    complexPredictionAccuracy_pDockQ
complexPredictionAccuracy_ipsae_d0dom_{AB,BA}    complexPredictionAccuracy_LIS{,_AB,_BA}
complexPredictionAccuracy_ipsae_d0res_{AB,BA}    complexPredictionAccuracy_ipsae_n0dom_{AB,BA}
complexPredictionAccuracy_ipsae_pae_cutoff       complexPredictionAccuracy_N_clash_backbone
complexPredictionAccuracy_iptm_af                complexPredictionAccuracy_ipTM
```

This module agrees with those to a worst |Δ| of **2.4e-6** over 231 comparisons. Four things to
know before comparing against them:

1. > **`complexPredictionAccuracy_ipTM` is AlphaFold's own ipTM, NOT `ipTM_d0chn`.**
   > It tracks `complexPredictionAccuracy_iptm_af` and is up to **0.054** away from
   > `ipTM_d0chn` (measured on FX-010). Substituting it produces a false failure on 6 of 7
   > fixtures. The field that corresponds to this project's value is
   > `complexPredictionAccuracy_ipsae_iptm_d0chn_{AB,BA}`.
2. **AFDB's stored precision has two regimes.** Some entries (FX-001, FX-002) store *every*
   field at 2 dp and **truncate**; others (FX-006…FX-010) store per-direction fields at full
   precision and **round** the rolled-up 2 dp fields. Applying either rule uniformly produces
   11 false failures. Method: use the full-precision per-direction fields wherever they exist;
   fall back to a 2 dp field only when there is nothing better, and then truncate.
3. **AFDB rolls `LIS` up as the max; `ipsae.py` and this module use the mean.** A convention
   difference, not a discrepancy. Compare `_LIS_AB` / `_LIS_BA` directionally.
4. `complexPredictionAccuracy_ipsae_dist_cutoff` is 15.0 and governs different counts than this
   project's contact `DIST_CUTOFF` of 8.0. `complexPredictionAccuracy_ipsae_pae_cutoff` is
   **10.0** — which settles the PAE-cutoff question: AFDB's production pipeline used 10, so the
   0.6 threshold transfers to these numbers directly and no caveat is needed.

### PAE JSON schema

```json
[{
  "predicted_aligned_error": [[float, ...]],
  "max_predicted_aligned_error": 27.47,
  "chains": [
    {"name": "Interferon-stimulated gene 20 kDa protein",
     "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 181}
  ]
}]
```

- Matrix is `(nA+nB) × (nA+nB)` and **asymmetric** (`PAE[i][j] != PAE[j][i]`).
- `chains[].name` is the chain's protein name. Use it for axis labels and view titles; a
  heterodimer labelled "Chain A / Chain B" tells the reader nothing.
- Use `max_predicted_aligned_error` as heatmap `vmax` (default 31.75 if absent).
- Assert that structure-derived chain lengths equal the PAE-derived `nA`/`nB`, and fail loudly.

### pLDDT JSON schema

```json
{"residueNumber": [...], "confidenceScore": [...], "confidenceCategory": [...],
 "chains": [{"name": "...", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 181}]}
```

Slice chain A `[:nA]`, chain B `[nA:nA+nB]` — **after** resolving chain identity, not before.

### mmCIF atom-site required fields

`group_PDB`, `label_atom_id`, `label_asym_id`, `label_seq_id`, `label_comp_id`,
`Cartn_x`, `Cartn_y`, `Cartn_z`, `B_iso_or_equiv`, `pdbx_PDB_model_num` (optional).

Column-key convention, chosen once in the module: the `_atom_site.` prefix is **stripped** and
keys are **lower-cased**, because mmCIF tags are case-insensitive by specification and AFDB
files mix cases within one loop (`group_PDB`, `Cartn_x`, `label_seq_id`).

Filtering: `group_PDB == 'ATOM'` only; `label_atom_id in ('CA','CB')` only; skip `label_seq_id`
of `'.'` or `'?'`; keep only model `'1'` when `pdbx_PDB_model_num` is present. **A blank line
inside the loop is not a terminator** — the loop ends at the next data name, `loop_`, or `#`.

---

## Interface Detection Contract

`detect_interface()` in `complex_interface_utils.py` is the canonical implementation.

- CB–CB Euclidean distance **≤ 8.0 Å** defines a contact (note `<=`, matching `ipsae_v4.py:652`).
- GLY has no CB; CA is substituted transparently.
- **This is CB–CB, not CA–CA.** AFDB's production pipeline (`interface.py`, Majewski,
  Apache-2.0, not in this repo) uses CA–CA with PyTorch `radius_graph`. The CB–CB choice
  matches the IPSAE scoring code's contact definition. Do not "fix" it without updating the
  scoring functions too.
- The distance matrix is `(nA, nB)`, built from an `(nA, nB, 3)` difference; fine up to
  ~1000-residue monomers (~24 MB).

---

## MolViewSpec Viewer Pattern

Use the module's builders (`show_mol_view` plus one builder per view); the pattern below is the
shape they follow.

```python
import molviewspec as mvs, base64
from IPython.display import display, HTML, IFrame

def show_mol_view(state, label, width=950, height=600):
    html = state.molstar_html()
    encoded = base64.b64encode(html.encode()).decode()
    display(HTML(f'<h4>{label}</h4>'))
    display(IFrame(src=f'data:text/html;base64,{encoded}', width=width, height=height))
```

- **Collapse contiguous same-colour residue runs into one `ComponentExpression`** using
  `beg_label_seq_id` / `end_label_seq_id`. One component per residue is what the notebook used
  to do; collapsing cut the pLDDT view from 344 components to 12.
- Resolve `label_asym_id` from the verified chain identities. Never hard-code `'A'` / `'B'`.
- Decide URL and format together (`resolve_structure_source`): an empty `bcifUrl` labelled
  `bcif` is a silent blank viewer.
- **Each view is a live Mol\* instance holding a WebGL context.** Browsers cap concurrent
  contexts at 16 and silently blank the oldest past the cap. The notebook sits at six; adding
  a seventh is a real cost, not a free addition.
- Drop non-finite values before colouring — a colormap paints "bad" values **opaque black**,
  the most emphatic colour in the scene, on exactly the residues with no measurement.
- Wrap the section in `try/except ImportError` so the notebook degrades gracefully.

---

## Known Edge Cases

| Condition | Handling |
|-----------|----------|
| Zero inter-chain contacts | **pDockQ → `0.0`, pDockQ2 → `0.0`** (`ipsae_v4.py:669, 700`). `0.018` / `0.005` are the sigmoid minima that Bryant's and Zhu's own code returns; `ipsae.py` deliberately short-circuits to `0.0` instead. Since ±0.001 agreement with `ipsae.py` is the acceptance criterion, `0.0` is required |
| No inter-chain `PAE < 12` | LIS → `0.0` |
| Residue with no valid PAE pair for ipSAE | that residue's per-residue entry stays `0.0` |
| GLY (no CB atom) | CA substituted transparently |
| `label_seq_id` = `'.'` or `'?'` | skip row |
| Multiple models | filter to `pdbx_PDB_model_num == '1'` |
| Blank line inside the `_atom_site` loop | **not** a terminator; keep parsing |
| PAE `chains` field absent | fall back to structure-derived chain lengths, then verify |
| Monomer or non-dimer accession | `UnsupportedAssemblyError` at the metadata gate, **before any download**. (Pre-2026-09 this was a bare `IndexError` at `pae.chain_ids[1]`; this file previously claimed it produced a zero-length chain B and scores of 0.0, which was never true) |
| Malformed / absent accession | `AccessionLookupError`, distinguishing malformed (HTTP 400, quoting the service's own message) from absent (404) |
| Chain ids disagree across structure / PAE / pLDDT | `verify_chain_identity` **raises**. Never map positionally |
| Prediction entries disagree on a shared field | `shared_field()` raises rather than trusting entry 0 |
| Local-file mode without PAE/pLDDT | fails **at upload**, naming every missing file. Six of seven values need the PAE matrix and the seventh needs pLDDT, so a "skip" would be a different notebook with a banner on it |
| Highly asymmetric chains (up to 8.16:1 measured) | supported; use `aspect='equal'` on every PAE panel — `aspect='auto'` rendered an 8.16:1 block and a 1:7.11 block as the same near-square |
| Short chain first | ~41% of AFDB heterodimers; FX-010 covers it |
| Raw `<` in interpolated HTML | escape it. `PAE < 10` in a description silently ate the rest of the table in a browser |

---

## Conventions

- **British spelling** in prose and identifiers where it is already used (`colour`, `normalise`).
- **No em dashes** in markdown cells, plot titles/labels or printed output. En dashes inside
  numeric ranges (`70–90`, `CB–CB`, `pp. 5–6`) are fine. Check for `&mdash;` and `&#8212;` in
  generated HTML too, which no Unicode grep will find.
- Plot functions return a bare `matplotlib.figure.Figure` built without pyplot, so nothing
  leaks into pyplot's registry and nothing is double-displayed. The notebook therefore needs
  `%matplotlib inline` to register the figure formatters.
- Importing the module never mutates global matplotlib state; `apply_plot_style()` is explicit.
- Default PAE colormap is `Greens_r` (dark green = low PAE = confident, legible on white);
  `viridis` is the colourblind-safe alternative and `RdBu_r` is kept for continuity but is
  **not** recommended — it is the only one of the three that loses monotone lightness under
  deuteranopia, protanopia and tritanopia.
- Notebook outputs are stripped before commit (D5).
