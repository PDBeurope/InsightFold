# Homodimer Diagnostic Notebook — Validation Report

Deliverable of **R090 + R091 + R092** in `specs/homodimer_diagnostic/rework-plan.md`, and of
**T038** in the original `specs/homodimer_diagnostic/tasks.md`, which called for this file and
which never had one.

**Date of run:** 2026-09-08.
**Subject:** `notebooks/homodimer_diagnostic.ipynb` at `homodimer-notebook-rework @ e162059`
("R040: em dash sweep. Closes M7"), against `src/insightfold/complex_interface_utils.py` at the
same revision.

**Result:** the notebook executes end to end, with zero errors and zero warnings, on all six
online positive fixtures and on the local-file path. All four negative paths refuse cleanly
before any score is produced. **All seven values agree with both independent references on all
seven scored fixtures. Nothing is outside the +-0.001 tolerance.** Worst deviation from
`ipsae.py` v4 is **4.7e-5**; worst deviation from AlphaFold DB's own published per-direction
values is **2.4e-6**.

Two things a reader should not mistake for failures, both documented and measured below:
`complexPredictionAccuracy_ipTM` is *not* ipSAE's `ipTM_d0chn` (Section 5.7), and AFDB's
rolled-up `LIS` is the max where this notebook reports the mean (Section 5.5).

One methodology correction is recorded against the R060 note (Section 5.4): AFDB's rounding
convention is **not** uniform. Truncation holds for entries stored entirely at 2 dp, but the
rolled-up 2 dp fields on full-precision entries are **rounded**. Applying R060's rule uniformly
produces 11 false failures across this fixture set; applying the regime-aware rule produces none.

---

## 1. Environment

Nothing here was pinned by the notebook; it is recorded so a later re-run can tell a genuine
regression from an environment change.

| Component | Version |
|---|---|
| OS | macOS 26.6.2, arm64 (Darwin 25.6.0) |
| Python | 3.11.9 (`main`, Aug 14 2024, Clang 18.1.8), from `.venv` |
| numpy | 2.3.3 |
| matplotlib | 3.10.7 |
| seaborn | 0.13.2 |
| requests | 2.32.5 |
| molviewspec | 1.8.1 |
| ipywidgets | 8.1.7 |
| nbclient / nbformat / ipykernel | 0.10.4 / 5.10.4 / 7.2.0 |
| Reference implementation | `specs/homodimer_diagnostic/references/ipsae_v4.py`, v4 (3 Jan 2026), sha256 `10cf9b08c68c91e06cb28526cf2026f47a3980c9048fd3226d13e3304eaf1c27` |
| AFDB endpoints | `/api/prediction/{acc}`, `/api/search?q=modelEntityId:{acc}&type=main&rows=1`, live on 2026-09-08 |

The bootstrap cell reported `Environment: local`, `Revision: homodimer-notebook-rework @ e162059`,
`molviewspec: available` on every run, so the Colab clone branch of that cell was never taken
(see Section 8, "What was not verified").

## 2. How the runs were driven

Each run reads `notebooks/homodimer_diagnostic.ipynb` from disk, changes **only** the
`ACCESSION_ID` string and, for the local-file cases, the `USE_LOCAL_FILE` flag in the user-input
cell, and executes the result under `nbclient` with the working directory set to `notebooks/`.
**The notebook on disk was not modified**; every patch is applied to an in-memory copy.

Two execution modes were used deliberately:

- **`allow_errors=False`** for the negative fixtures. This is what Jupyter and Colab "Run All"
  actually do: execution halts at the first failing cell. It is the mode that tells you what a
  user sees.
- **`allow_errors=True`** for the positive fixtures, so that a failure late in the notebook could
  not be masked by an early halt.

Local-file mode was driven by replacing `ipywidgets.FileUpload` with a subclass that arrives
pre-populated from files on disk. The notebook's own upload cells, its own
`require_local_documents` gate and its own parsing all run unchanged; only the click is
simulated.

Full-precision values were read back by appending one extra cell to the in-memory copy that
prints `repr`-level floats off the objects the notebook had already built (`res_ipsae`,
`res_pdockq`, `res_pdockq2`, `res_lis`, `scores`). That cell computes nothing. It exists because
the notebook prints to 4 dp, and a +-0.001 tolerance checked against a 4 dp print is a tolerance
checked against its own rounding.

## 3. R090 + R091 — execution matrix

All ten registered fixtures, plus the three local-file cases the manifest lists.

| # | Fixture | Accession | Mode | Cells run | Errors | Warnings (stderr) | Figures | Mol* views | Wall | Outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | FX-001 homodimer 172+172 | `AF-0000000065889468` | online | 33 / 33 | 0 | 0 | 6 | 6 | 4.66 s | clean |
| 2 | FX-002 homodimer 123+123 | `AF-0000000066503175` | online | 33 / 33 | 0 | 0 | 6 | 6 | 4.16 s | clean |
| 3 | FX-006 heterodimer 181+101 | `AF-0000000211034637` | online | 33 / 33 | 0 | 0 | 6 | 6 | 3.94 s | clean |
| 4 | FX-007 heterodimer 563+69 (8.16:1) | `AF-0000000204661110` | online | 33 / 33 | 0 | 0 | 6 | 6 | 4.39 s | clean |
| 5 | FX-008 borderline heterodimer 296+88 | `AF-0000000211619209` | online | 33 / 33 | 0 | 0 | 6 | 6 | 3.98 s | clean |
| 6 | FX-009 heterodimer 699+450 | `AF-0000000211157965` | online | 33 / 33 | 0 | 0 | 6 | 6 | 5.28 s | clean |
| 7 | FX-010 heterodimer 108+768 (short-first) | `AF-0000000211026350` | online | 33 / 33 | 0 | 0 | 6 | 6 | 4.75 s | clean |
| 8 | FX-004 malformed accession | `AF-NOT_A_REAL_ACCESSION` | online | 5 / 33 | halts at cell 6 | 0 | 0 | 0 | 2.49 s | refused: `AccessionLookupError` |
| 9 | FX-005 monomer | `AF-O15552-F1` | online | 5 / 33 | halts at cell 6 | 0 | 0 | 0 | 2.49 s | refused: `UnsupportedAssemblyError` |
| 10a | local-file, all three documents | FX-001 documents | local | 34 / 34 | 0 | 0 | 6 | 0 (by design) | 3.49 s | clean; scores bit-identical to run 1 |
| 10b | local-file, mmCIF only | FX-001 mmCIF | local | 7 / 33 | halts at cell 8 | 0 | 0 | 0 | 2.34 s | refused: `MissingLocalDocumentError` naming both missing files |
| 10c | local-file, mmCIF + PAE, no pLDDT | FX-001 mmCIF + PAE | local | 7 / 33 | halts at cell 8 | 0 | 0 | 0 | 2.33 s | refused: `MissingLocalDocumentError` naming the pLDDT file |
| 10d | local-file, monomer (extra) | FX-005 documents | local | 9 / 33 | halts at cell 10 | 0 | 0 | 0 | 2.35 s | refused: `UnsupportedAssemblyError` from the structural gate |

Notes on the matrix.

- **"Cells run"** counts code cells with an execution count. The positive runs carry 33 code
  cells; the local-file run in row 10a shows 34 because the harness prelude cell is counted.
- **Figures** are the six matplotlib PNGs: interface contact map (cell 13), PAE matrix (17),
  PAE score masks (18), per-residue score profiles (35), interface pLDDT (37), threshold
  margins (55). Cell 53 additionally renders the summary table as HTML.
- **Mol\* views** are the six `<iframe>` viewers: View 1 in Section 2 (cell 15) and Views 2-6 in
  Section 6 (cells 41, 43, 45, 47, 49). All six rendered on every online fixture.
- **Zero warnings** means zero bytes were written to stderr on any run: no `DeprecationWarning`,
  no `RuntimeWarning` from numpy, no matplotlib layout warning.
- **Local-file mode renders no Mol\* views**, and says so: *"Skipping the 3D views: Mol\*
  downloads the structure itself, so local file mode has no URL to hand it. The contact map
  above, and every score in this notebook, are unaffected."* That is the designed behaviour, not
  a failure.
- **Chain identity and assembly detection were correct on every fixture**: `Homodimer -- 2 chains
  (A, B)` on FX-001/FX-002, `Heterodimer -- 2 chains (A, B)` on the five heterodimers, with real
  protein labels (`ISG20 (A)` / `Sumo1 (B)`, `Sptlc3 (A)` / `Gm6993 (B)`, `PACRG (A)` /
  `MEIG1 (B)`, `UVRAG (A)` / `BECN1 (B)`, `Rbx1 (A)` / `CUL3 (B)`). Chain lengths matched the
  manifest exactly in every case, including the short-first orientation of FX-010 (108 then 768).

### 3.1 Scores and bands, as the notebook reported them

| Fixture | ipSAE_d0res | ipSAE_d0chn | ipSAE_d0dom | ipTM_d0chn | pDockQ | pDockQ2 | LIS | AFDB rule |
|---|---|---|---|---|---|---|---|---|
| FX-001 | 0.9143 VERY HIGH-CONFIDENCE | 0.9529 HIGH | 0.9527 HIGH | 0.9529 HIGH | 0.6913 HIGH | 0.9269 HIGH | 0.7564 HIGH | PASS |
| FX-002 | 0.7699 CONFIDENT | 0.8660 HIGH | 0.8553 HIGH | 0.8131 HIGH | 0.6877 HIGH | 0.7282 HIGH | 0.6592 HIGH | PASS |
| FX-006 | 0.7057 CONFIDENT | 0.8055 HIGH | 0.7891 HIGH | 0.7730 HIGH | 0.1452 MODERATE | 0.7054 HIGH | 0.6004 HIGH | PASS |
| FX-007 | 0.7120 CONFIDENT | 0.8181 HIGH | 0.7856 HIGH | 0.6687 MODERATE | 0.2737 HIGH | 0.5914 HIGH | 0.3980 HIGH | PASS |
| FX-008 | 0.6366 LOW-CONFIDENCE | 0.8313 HIGH | 0.7735 HIGH | 0.8183 HIGH | 0.2418 HIGH | 0.3254 HIGH | 0.4783 HIGH | PASS |
| FX-009 | 0.7119 CONFIDENT | 0.8947 HIGH | 0.8691 HIGH | 0.5769 MODERATE | 0.7191 HIGH | 0.2899 HIGH | 0.4278 HIGH | PASS |
| FX-010 | 0.6888 LOW-CONFIDENCE | 0.8597 HIGH | 0.7992 HIGH | 0.5958 MODERATE | 0.5805 HIGH | 0.2721 HIGH | 0.3841 HIGH | PASS |

Every fixture passes the AFDB joint release rule (`ipSAE_d0res >= 0.60 AND pDockQ2 >= 0.23`),
which is expected: all seven are AFDB-released complexes, so all seven already cleared that
filter upstream. Contact counts and interface sizes:

| Fixture | Contact pairs | Interface residues (chain A / chain B) | Mean interface pLDDT |
|---|---:|---|---:|
| FX-001 | 116 | 43 / 172 (25.0%) and 43 / 172 (25.0%) | 98.01 |
| FX-002 | 152 | 52 / 123 (42.3%) and 52 / 123 (42.3%) | 92.08 |
| FX-006 | 22 | 11 / 181 (6.1%) and 11 / 101 (10.9%) | 91.54 |
| FX-007 | 42 | 21 / 563 (3.7%) and 16 / 69 (23.2%) | 86.85 |
| FX-008 | 43 | 20 / 296 (6.8%) and 16 / 88 (18.2%) | 83.96 |
| FX-009 | 501 | 181 / 699 (25.9%) and 184 / 450 (40.9%) | 80.90 |
| FX-010 | 155 | 32 / 108 (29.6%) and 60 / 768 (7.8%) | 80.63 |

### 3.2 FX-002 is scored, and is not the borderline fixture

`fixture-manifest.md` left FX-002 at `candidate-needs-scoring` and asked R090 to "score it or
retire it". Scored: `ipSAE_d0res` **0.7699**, which lands in the `CONFIDENT` band, well clear of
the 0.6 release cutoff and of the 0.7 green threshold. It is **not** borderline. FX-008
(`ipSAE_d0res` 0.6366, `LOW-CONFIDENCE`) remains the only fixture in the set that sits in AFDB's
lowest released band, and it should keep the borderline role. **Recommendation: re-label FX-002
in the manifest as a second high-confidence homodimer** — it is still useful as the smallest and
fastest positive fixture, and as a second entry in the all-2 dp AFDB storage regime — **and drop
the "needs scoring" status.** (Manifest edit deferred: this task may not modify spec files.)

### 3.3 The R080 defect recorded in the manifest no longer reproduces

`fixture-manifest.md` records, under "What is still open", that Section 7 prints "consistently
HIGH confidence across all metrics" on FX-008 and FX-009 while a traffic light is amber. **That
does not happen at this revision.** Section 7 now reads, verbatim:

- FX-008: *"OVERALL: mixed. 4 of 5 independent scores are green; ipSAE_d0res 0.637
  (LOW-CONFIDENCE) is not. The scores disagree, and the paragraphs below say where."*
- FX-009: *"OVERALL: mixed. 4 of 5 independent scores are green; ipTM_d0chn 0.577 (MODERATE) is
  not."*
- FX-010: *"OVERALL: mixed. 3 of 5 independent scores are green; ipTM_d0chn 0.596 (MODERATE) and
  ipSAE_d0res 0.689 (LOW-CONFIDENCE) are not."*

The count is taken from the same bands the summary table paints. The open item in the manifest
can be closed by whoever next edits it.

---

## 4. R092 — method

Seven values, on seven fixtures, against two references that share no code with this notebook.

### 4.1 Reference (a): `ipsae.py` v4 as a subprocess

`references/ipsae_v4.py` is a CLI. It was invoked exactly as specified, once per fixture:

```
python ipsae_v4.py <af3_pae.json> <model.cif> 10 8
```

AFDB does not serve the shape this CLI expects. It serves the PAE matrix under
`predicted_aligned_error` in one document and per-residue pLDDT in another, while the CLI's
AF3 path wants a single `{"pae": ..., "atom_plddts": ...}` object. The repackaging was
mechanical and lossless:

- `pae` is `pae.json[0]["predicted_aligned_error"]`, copied verbatim.
- `atom_plddts` is the mmCIF `_atom_site.B_iso_or_equiv` column, one entry per atom, ordered by
  `_atom_site.id` — because `ipsae.py` indexes it as `atom_num - 1` (`ipsae_v4.py:531-534`) to
  pull the CA and CB pLDDT it needs. The atom ids were asserted contiguous from 1 on every
  fixture before writing.
- **No number was altered, rounded or re-derived.** The confidence JSON was not used for this
  reference at all, precisely so that the two paths stay independent: the notebook reads
  per-residue pLDDT from `confidence.json`, and `ipsae.py` reads it from the structure's
  B-factors.

Two properties of the CLI's own output limit how tightly it can be compared, and both are
recorded rather than worked around:

- It prints `pDockQ`, `pDockQ2` and `LIS` with `%8.4f` and the ipSAE variants with `%8.6f`
  (`ipsae_v4.py:948-967`). **Every non-zero delta below is at or under the reference's own print
  precision.** The largest, 4.7e-5, is just under half a unit in the last printed place of a
  4 dp field, which is the most a correctly rounded 4 dp print can be wrong by.
- Its `ipTM_af` column reads `0.000` on every fixture, because that column is copied from an
  AlphaFold 3 `summary_confidences.json` that AFDB does not publish. It is not compared. The
  quantity this notebook calls `ipTM_d0chn` is the CLI's `ipTM_d0chn` column, which is populated.

Rows consumed: the CLI emits one row per direction plus a combined `max` row. The `max` row is
what the notebook's reported values are compared against, because it is the row that applies
`ipsae.py`'s own combination rules — `max` for the ipSAE variants, `ipTM_d0chn` and `pDockQ2`,
the plain symmetric value for `pDockQ`, and the **mean** for LIS (`ipsae_v4.py:982`).

### 4.2 Reference (b): AlphaFold DB's own published production values

```
GET https://alphafold.ebi.ac.uk/api/search?q=modelEntityId:<accession>&type=main&rows=1
```

returned `numFound: 1` for all seven fixtures. The `complexPredictionAccuracy_*` fields are
produced by AFDB's own IPSAE run over its own pipeline, so they are a genuinely independent
check rather than a self-comparison.

The parameters AFDB published alongside the scores were checked first, because a mismatch there
would invalidate every value comparison that follows:

| Parameter | AFDB, all seven fixtures | This notebook |
|---|---|---|
| `complexPredictionAccuracy_ipsae_pae_cutoff` | 10.0 | `PAE_CUTOFF` 10.0 |
| `complexPredictionAccuracy_ipsae_dist_cutoff` | 15.0 | `DIST_CUTOFF` 8.0 — **different quantity, not a mismatch** |

AFDB's `ipsae_dist_cutoff` governs its `ipsae_dist_nres*` residue counts. This notebook's
`DIST_CUTOFF` governs the CB-CB contact set that feeds pDockQ and pDockQ2, and is the same 8 that
was passed to `ipsae.py` as the fourth argument. The two never meet. That pDockQ agrees with
AFDB at every published decimal place on all seven fixtures is the evidence that the contact
definitions have not diverged.

---

## 5. R092 — results

### 5.1 Reference (a): the seven reported values vs `ipsae.py` v4

| Fixture | Value | Notebook | `ipsae.py` v4 | \|delta\| | Verdict |
|---|---|---:|---:|---:|---|
| FX-001 | ipSAE_d0res | 0.914309 | 0.914309 | 0.0000001 | PASS |
| FX-001 | ipSAE_d0chn | 0.952902 | 0.952902 | 0.0000002 | PASS |
| FX-001 | ipSAE_d0dom | 0.952672 | 0.952672 | 0.0000004 | PASS |
| FX-001 | ipTM_d0chn | 0.952902 | 0.952902 | 0.0000002 | PASS |
| FX-001 | pDockQ | 0.691274 | 0.691300 | 0.0000264 | PASS |
| FX-001 | pDockQ2 | 0.926932 | 0.926900 | 0.0000318 | PASS |
| FX-001 | LIS | 0.756447 | 0.756400 | 0.0000470 | PASS |
| FX-002 | ipSAE_d0res | 0.769888 | 0.769888 | 0.0000001 | PASS |
| FX-002 | ipSAE_d0chn | 0.865987 | 0.865987 | 0.0000004 | PASS |
| FX-002 | ipSAE_d0dom | 0.855267 | 0.855267 | 0.0000004 | PASS |
| FX-002 | ipTM_d0chn | 0.813102 | 0.813102 | 0.0000002 | PASS |
| FX-002 | pDockQ | 0.687671 | 0.687700 | 0.0000287 | PASS |
| FX-002 | pDockQ2 | 0.728214 | 0.728200 | 0.0000137 | PASS |
| FX-002 | LIS | 0.659187 | 0.659200 | 0.0000131 | PASS |
| FX-006 | ipSAE_d0res | 0.705718 | 0.705718 | 0.0000001 | PASS |
| FX-006 | ipSAE_d0chn | 0.805532 | 0.805532 | 0.0000003 | PASS |
| FX-006 | ipSAE_d0dom | 0.789086 | 0.789086 | 0.0000005 | PASS |
| FX-006 | ipTM_d0chn | 0.772954 | 0.772954 | 0.0000004 | PASS |
| FX-006 | pDockQ | 0.145233 | 0.145200 | 0.0000326 | PASS |
| FX-006 | pDockQ2 | 0.705404 | 0.705400 | 0.0000035 | PASS |
| FX-006 | LIS | 0.600423 | 0.600400 | 0.0000233 | PASS |
| FX-007 | ipSAE_d0res | 0.711961 | 0.711961 | 0.0000001 | PASS |
| FX-007 | ipSAE_d0chn | 0.818091 | 0.818091 | 0.0000000 | PASS |
| FX-007 | ipSAE_d0dom | 0.785619 | 0.785619 | 0.0000004 | PASS |
| FX-007 | ipTM_d0chn | 0.668721 | 0.668721 | 0.0000005 | PASS |
| FX-007 | pDockQ | 0.273730 | 0.273700 | 0.0000301 | PASS |
| FX-007 | pDockQ2 | 0.591426 | 0.591400 | 0.0000259 | PASS |
| FX-007 | LIS | 0.397974 | 0.398000 | 0.0000258 | PASS |
| FX-008 | ipSAE_d0res | 0.636590 | 0.636590 | 0.0000004 | PASS |
| FX-008 | ipSAE_d0chn | 0.831314 | 0.831314 | 0.0000004 | PASS |
| FX-008 | ipSAE_d0dom | 0.773544 | 0.773544 | 0.0000001 | PASS |
| FX-008 | ipTM_d0chn | 0.818349 | 0.818349 | 0.0000003 | PASS |
| FX-008 | pDockQ | 0.241786 | 0.241800 | 0.0000138 | PASS |
| FX-008 | pDockQ2 | 0.325409 | 0.325400 | 0.0000089 | PASS |
| FX-008 | LIS | 0.478336 | 0.478300 | 0.0000364 | PASS |
| FX-009 | ipSAE_d0res | 0.711865 | 0.711865 | 0.0000001 | PASS |
| FX-009 | ipSAE_d0chn | 0.894745 | 0.894745 | 0.0000000 | PASS |
| FX-009 | ipSAE_d0dom | 0.869095 | 0.869095 | 0.0000001 | PASS |
| FX-009 | ipTM_d0chn | 0.576863 | 0.576863 | 0.0000002 | PASS |
| FX-009 | pDockQ | 0.719113 | 0.719100 | 0.0000129 | PASS |
| FX-009 | pDockQ2 | 0.289935 | 0.289900 | 0.0000347 | PASS |
| FX-009 | LIS | 0.427799 | 0.427800 | 0.0000008 | PASS |
| FX-010 | ipSAE_d0res | 0.688783 | 0.688783 | 0.0000002 | PASS |
| FX-010 | ipSAE_d0chn | 0.859679 | 0.859679 | 0.0000003 | PASS |
| FX-010 | ipSAE_d0dom | 0.799233 | 0.799233 | 0.0000001 | PASS |
| FX-010 | ipTM_d0chn | 0.595829 | 0.595829 | 0.0000004 | PASS |
| FX-010 | pDockQ | 0.580514 | 0.580500 | 0.0000141 | PASS |
| FX-010 | pDockQ2 | 0.272080 | 0.272100 | 0.0000199 | PASS |
| FX-010 | LIS | 0.384140 | 0.384100 | 0.0000403 | PASS |

**49 of 49 comparisons pass. Worst |delta| = 4.7e-5, 21x inside the +-0.001 tolerance.**

Every delta above 1e-6 is on `pDockQ`, `pDockQ2` or `LIS` — the three the reference prints to
4 dp. Every ipSAE variant and every `ipTM_d0chn`, printed to 6 dp, agrees to <= 6e-7. The
pattern is the reference's print width, not a numerical disagreement.

The per-direction values were compared too (12 per fixture, 84 in total, not tabulated here for
length): **worst |delta| 4.8e-5**, same distribution, all inside tolerance.

### 5.2 Reference (b): per-direction values vs AFDB

| Fixture | Regime | Value | Notebook | AFDB field | AFDB | Comparison | \|delta\| | Verdict |
|---|---|---|---:|---|---:|---|---:|---|
| FX-001 | TWO-DP | ipSAE_d0res A->B | 0.914309 | `ipsae_AB` | 0.91 | truncate to 0.91 | 0.0000000 | PASS |
| FX-001 | TWO-DP | ipSAE_d0chn A->B | 0.952902 | `ipsae_d0chn_AB` | 0.95 | truncate to 0.95 | 0.0000000 | PASS |
| FX-001 | TWO-DP | ipSAE_d0dom A->B | 0.952672 | (not published) | -- | -- | -- | n/a |
| FX-001 | TWO-DP | ipTM_d0chn A->B | 0.952902 | `ipsae_iptm_d0chn_AB` | 0.95 | truncate to 0.95 | 0.0000000 | PASS |
| FX-001 | TWO-DP | pDockQ2 A->B | 0.926932 | `pDockQ2_AB` | 0.92 | truncate to 0.92 | 0.0000000 | PASS |
| FX-001 | TWO-DP | LIS A->B | 0.756442 | `LIS_AB` | 0.75 | truncate to 0.75 | 0.0000000 | PASS |
| FX-001 | TWO-DP | ipSAE_d0res B->A | 0.913848 | `ipsae_BA` | 0.91 | truncate to 0.91 | 0.0000000 | PASS |
| FX-001 | TWO-DP | ipSAE_d0chn B->A | 0.952635 | `ipsae_d0chn_BA` | 0.95 | truncate to 0.95 | 0.0000000 | PASS |
| FX-001 | TWO-DP | ipSAE_d0dom B->A | 0.952404 | (not published) | -- | -- | -- | n/a |
| FX-001 | TWO-DP | ipTM_d0chn B->A | 0.952635 | `ipsae_iptm_d0chn_BA` | 0.95 | truncate to 0.95 | 0.0000000 | PASS |
| FX-001 | TWO-DP | pDockQ2 B->A | 0.926889 | `pDockQ2_BA` | 0.92 | truncate to 0.92 | 0.0000000 | PASS |
| FX-001 | TWO-DP | LIS B->A | 0.756452 | `LIS_BA` | 0.75 | truncate to 0.75 | 0.0000000 | PASS |
| FX-002 | TWO-DP | ipSAE_d0res A->B | 0.769888 | `ipsae_AB` | 0.76 | truncate to 0.76 | 0.0000000 | PASS |
| FX-002 | TWO-DP | ipSAE_d0chn A->B | 0.865987 | `ipsae_d0chn_AB` | 0.86 | truncate to 0.86 | 0.0000000 | PASS |
| FX-002 | TWO-DP | ipSAE_d0dom A->B | 0.855267 | (not published) | -- | -- | -- | n/a |
| FX-002 | TWO-DP | ipTM_d0chn A->B | 0.813102 | `ipsae_iptm_d0chn_AB` | 0.81 | truncate to 0.81 | 0.0000000 | PASS |
| FX-002 | TWO-DP | pDockQ2 A->B | 0.727737 | `pDockQ2_AB` | 0.72 | truncate to 0.72 | 0.0000000 | PASS |
| FX-002 | TWO-DP | LIS A->B | 0.659194 | `LIS_AB` | 0.65 | truncate to 0.65 | 0.0000000 | PASS |
| FX-002 | TWO-DP | ipSAE_d0res B->A | 0.768161 | `ipsae_BA` | 0.76 | truncate to 0.76 | 0.0000000 | PASS |
| FX-002 | TWO-DP | ipSAE_d0chn B->A | 0.864932 | `ipsae_d0chn_BA` | 0.86 | truncate to 0.86 | 0.0000000 | PASS |
| FX-002 | TWO-DP | ipSAE_d0dom B->A | 0.854133 | (not published) | -- | -- | -- | n/a |
| FX-002 | TWO-DP | ipTM_d0chn B->A | 0.812029 | `ipsae_iptm_d0chn_BA` | 0.81 | truncate to 0.81 | 0.0000000 | PASS |
| FX-002 | TWO-DP | pDockQ2 B->A | 0.728214 | `pDockQ2_BA` | 0.72 | truncate to 0.72 | 0.0000000 | PASS |
| FX-002 | TWO-DP | LIS B->A | 0.659179 | `LIS_BA` | 0.65 | truncate to 0.65 | 0.0000000 | PASS |
| FX-006 | FULL | ipSAE_d0res A->B | 0.555450 | `ipsae_AB` | 0.555450 | exact | 0.0000005 | PASS |
| FX-006 | FULL | ipSAE_d0chn A->B | 0.805532 | `ipsae_d0chn_AB` | 0.805532 | exact | 0.0000003 | PASS |
| FX-006 | FULL | ipSAE_d0dom A->B | 0.789086 | (not published) | -- | -- | -- | n/a |
| FX-006 | FULL | ipTM_d0chn A->B | 0.668714 | `ipsae_iptm_d0chn_AB` | 0.668714 | exact | 0.0000002 | PASS |
| FX-006 | FULL | pDockQ2 A->B | 0.705404 | `pDockQ2_AB` | 0.705403 | exact | 0.0000005 | PASS |
| FX-006 | FULL | LIS A->B | 0.608776 | `LIS_AB` | 0.608774 | exact | 0.0000022 | PASS |
| FX-006 | FULL | ipSAE_d0res B->A | 0.705718 | `ipsae_BA` | 0.705718 | exact | 0.0000001 | PASS |
| FX-006 | FULL | ipSAE_d0chn B->A | 0.785167 | `ipsae_d0chn_BA` | 0.785167 | exact | 0.0000001 | PASS |
| FX-006 | FULL | ipSAE_d0dom B->A | 0.768065 | (not published) | -- | -- | -- | n/a |
| FX-006 | FULL | ipTM_d0chn B->A | 0.772954 | `ipsae_iptm_d0chn_BA` | 0.772954 | exact | 0.0000004 | PASS |
| FX-006 | FULL | pDockQ2 B->A | 0.685271 | `pDockQ2_BA` | 0.685271 | exact | 0.0000000 | PASS |
| FX-006 | FULL | LIS B->A | 0.592070 | `LIS_BA` | 0.592071 | exact | 0.0000006 | PASS |
| FX-007 | FULL | ipSAE_d0res A->B | 0.312761 | `ipsae_AB` | 0.312761 | exact | 0.0000002 | PASS |
| FX-007 | FULL | ipSAE_d0chn A->B | 0.818091 | `ipsae_d0chn_AB` | 0.818091 | exact | 0.0000000 | PASS |
| FX-007 | FULL | ipSAE_d0dom A->B | 0.785619 | (not published) | -- | -- | -- | n/a |
| FX-007 | FULL | ipTM_d0chn A->B | 0.646488 | `ipsae_iptm_d0chn_AB` | 0.646488 | exact | 0.0000003 | PASS |
| FX-007 | FULL | pDockQ2 A->B | 0.591426 | `pDockQ2_AB` | 0.591426 | exact | 0.0000001 | PASS |
| FX-007 | FULL | LIS A->B | 0.417852 | `LIS_AB` | 0.417853 | exact | 0.0000008 | PASS |
| FX-007 | FULL | ipSAE_d0res B->A | 0.711961 | `ipsae_BA` | 0.711961 | exact | 0.0000001 | PASS |
| FX-007 | FULL | ipSAE_d0chn B->A | 0.752181 | `ipsae_d0chn_BA` | 0.752181 | exact | 0.0000004 | PASS |
| FX-007 | FULL | ipSAE_d0dom B->A | 0.726331 | (not published) | -- | -- | -- | n/a |
| FX-007 | FULL | ipTM_d0chn B->A | 0.668721 | `ipsae_iptm_d0chn_BA` | 0.668721 | exact | 0.0000005 | PASS |
| FX-007 | FULL | pDockQ2 B->A | 0.556196 | `pDockQ2_BA` | 0.556197 | exact | 0.0000006 | PASS |
| FX-007 | FULL | LIS B->A | 0.378096 | `LIS_BA` | 0.378097 | exact | 0.0000009 | PASS |
| FX-008 | FULL | ipSAE_d0res A->B | 0.544477 | `ipsae_AB` | 0.544477 | exact | 0.0000002 | PASS |
| FX-008 | FULL | ipSAE_d0chn A->B | 0.831314 | `ipsae_d0chn_AB` | 0.831315 | exact | 0.0000006 | PASS |
| FX-008 | FULL | ipSAE_d0dom A->B | 0.773544 | (not published) | -- | -- | -- | n/a |
| FX-008 | FULL | ipTM_d0chn A->B | 0.818349 | `ipsae_iptm_d0chn_AB` | 0.818349 | exact | 0.0000003 | PASS |
| FX-008 | FULL | pDockQ2 A->B | 0.250344 | `pDockQ2_AB` | 0.250344 | exact | 0.0000001 | PASS |
| FX-008 | FULL | LIS A->B | 0.496848 | `LIS_AB` | 0.496848 | exact | 0.0000002 | PASS |
| FX-008 | FULL | ipSAE_d0res B->A | 0.636590 | `ipsae_BA` | 0.636590 | exact | 0.0000004 | PASS |
| FX-008 | FULL | ipSAE_d0chn B->A | 0.763113 | `ipsae_d0chn_BA` | 0.763113 | exact | 0.0000003 | PASS |
| FX-008 | FULL | ipSAE_d0dom B->A | 0.705959 | (not published) | -- | -- | -- | n/a |
| FX-008 | FULL | ipTM_d0chn B->A | 0.520393 | `ipsae_iptm_d0chn_BA` | 0.520393 | exact | 0.0000001 | PASS |
| FX-008 | FULL | pDockQ2 B->A | 0.325409 | `pDockQ2_BA` | 0.325409 | exact | 0.0000001 | PASS |
| FX-008 | FULL | LIS B->A | 0.459825 | `LIS_BA` | 0.459826 | exact | 0.0000011 | PASS |
| FX-009 | FULL | ipSAE_d0res A->B | 0.711865 | `ipsae_AB` | 0.711865 | exact | 0.0000001 | PASS |
| FX-009 | FULL | ipSAE_d0chn A->B | 0.894745 | `ipsae_d0chn_AB` | 0.894745 | exact | 0.0000000 | PASS |
| FX-009 | FULL | ipSAE_d0dom A->B | 0.869095 | (not published) | -- | -- | -- | n/a |
| FX-009 | FULL | ipTM_d0chn A->B | 0.576863 | `ipsae_iptm_d0chn_AB` | 0.576863 | exact | 0.0000002 | PASS |
| FX-009 | FULL | pDockQ2 A->B | 0.262381 | `pDockQ2_AB` | 0.262381 | exact | 0.0000003 | PASS |
| FX-009 | FULL | LIS A->B | 0.440766 | `LIS_AB` | 0.440768 | exact | 0.0000024 | PASS |
| FX-009 | FULL | ipSAE_d0res B->A | 0.672218 | `ipsae_BA` | 0.672218 | exact | 0.0000003 | PASS |
| FX-009 | FULL | ipSAE_d0chn B->A | 0.879691 | `ipsae_d0chn_BA` | 0.879691 | exact | 0.0000001 | PASS |
| FX-009 | FULL | ipSAE_d0dom B->A | 0.854967 | (not published) | -- | -- | -- | n/a |
| FX-009 | FULL | ipTM_d0chn B->A | 0.418560 | `ipsae_iptm_d0chn_BA` | 0.418560 | exact | 0.0000003 | PASS |
| FX-009 | FULL | pDockQ2 B->A | 0.289935 | `pDockQ2_BA` | 0.289935 | exact | 0.0000003 | PASS |
| FX-009 | FULL | LIS B->A | 0.414833 | `LIS_BA` | 0.414834 | exact | 0.0000011 | PASS |
| FX-010 | FULL | ipSAE_d0res A->B | 0.688783 | `ipsae_AB` | 0.688783 | exact | 0.0000002 | PASS |
| FX-010 | FULL | ipSAE_d0chn A->B | 0.773648 | `ipsae_d0chn_AB` | 0.773649 | exact | 0.0000008 | PASS |
| FX-010 | FULL | ipSAE_d0dom A->B | 0.698266 | (not published) | -- | -- | -- | n/a |
| FX-010 | FULL | ipTM_d0chn A->B | 0.595829 | `ipsae_iptm_d0chn_AB` | 0.595829 | exact | 0.0000004 | PASS |
| FX-010 | FULL | pDockQ2 A->B | 0.238004 | `pDockQ2_AB` | 0.238003 | exact | 0.0000006 | PASS |
| FX-010 | FULL | LIS A->B | 0.421410 | `LIS_AB` | 0.421410 | exact | 0.0000004 | PASS |
| FX-010 | FULL | ipSAE_d0res B->A | 0.120251 | `ipsae_BA` | 0.120251 | exact | 0.0000002 | PASS |
| FX-010 | FULL | ipSAE_d0chn B->A | 0.859679 | `ipsae_d0chn_BA` | 0.859679 | exact | 0.0000003 | PASS |
| FX-010 | FULL | ipSAE_d0dom B->A | 0.799233 | (not published) | -- | -- | -- | n/a |
| FX-010 | FULL | ipTM_d0chn B->A | 0.418319 | `ipsae_iptm_d0chn_BA` | 0.418319 | exact | 0.0000003 | PASS |
| FX-010 | FULL | pDockQ2 B->A | 0.272080 | `pDockQ2_BA` | 0.272080 | exact | 0.0000001 | PASS |
| FX-010 | FULL | LIS B->A | 0.346870 | `LIS_BA` | 0.346871 | exact | 0.0000007 | PASS |

**70 of 70 comparisons pass** (84 table rows, of which 14 are the unpublished `ipSAE_d0dom`).
On the five full-precision entries the worst |delta| is
**2.4e-6** — about 400x inside tolerance. On the two 2 dp entries every truncated value matches
the stored value exactly.

`ipSAE_d0dom` is the one score AFDB does not publish. Its two intermediates are published and
both match exactly, per direction, which is the asymmetry R003 existed to fix:

| Fixture | n0dom A->B | AFDB | n0dom B->A | AFDB | d0dom A->B | AFDB | d0dom B->A | AFDB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| FX-001 | 342 | 342 | 342 | 342 | 6.74288 | 6.74 | 6.74288 | 6.74 |
| FX-002 | 220 | 220 | 220 | 220 | 5.51150 | 5.51 | 5.51150 | 5.51 |
| FX-006 | 250 | 250 | 252 | 252 | 5.85205 | 5.85205 | 5.87369 | 5.87369 |
| FX-007 | 475 | 475 | 532 | 532 | 7.77211 | 7.77211 | 8.15219 | 8.15219 |
| FX-008 | 249 | 249 | 268 | 268 | 5.84118 | 5.84118 | 6.04263 | 6.04263 |
| FX-009 | 821 | 821 | 843 | 843 | 9.73985 | 9.73985 | 9.84390 | 9.8439 |
| FX-010 | 513 | 513 | 476 | 476 | 8.02875 | 8.02875 | 7.77904 | 7.77904 |

`n0dom` differs between the two directions on all five heterodimers (250/252, 475/532, 249/268,
821/843, 513/476) and AFDB's per-direction values match ours integer for integer. `d0dom` matches
at every decimal place AFDB stores (5 or 6 significant figures).

### 5.3 Reference (b): rolled-up values vs AFDB, and the precision regimes

Each fixture falls in one of two storage regimes, detected from the data rather than assumed —
by asking whether AFDB's per-direction score fields carry more than two decimal places:

| Regime | Fixtures | Per-direction fields | Rolled-up fields |
|---|---|---|---|
| **TWO-DP** | FX-001, FX-002 (both homodimers) | stored at 2 dp, **truncated** | stored at 2 dp, **truncated** |
| **FULL** | FX-006, FX-007, FX-008, FX-009, FX-010 (all five heterodimers) | full precision (~1e-6) | stored at 2 dp, **rounded** |
| Fixture | Regime | Value (ours) | AFDB field | AFDB | Applied | \|delta\| | Other convention | \|delta\| | Verdict |
|---|---|---:|---|---:|---|---:|---|---:|---|
| FX-001 | TWO-DP | ipSAE_d0res (max) 0.914309 | `ipSAE` | 0.91 | truncate -> 0.91 | 0.0000 | round -> 0.91 | 0.0000 | PASS |
| FX-001 | TWO-DP | pDockQ 0.691274 | `pDockQ` | 0.69 | truncate -> 0.69 | 0.0000 | round -> 0.69 | 0.0000 | PASS |
| FX-001 | TWO-DP | pDockQ2 (max) 0.926932 | `pDockQ2` | 0.92 | truncate -> 0.92 | 0.0000 | round -> 0.93 | 0.0100 | PASS |
| FX-001 | TWO-DP | LIS (max, AFDB convention) 0.756452 | `LIS` | 0.75 | truncate -> 0.75 | 0.0000 | round -> 0.76 | 0.0100 | PASS |
| FX-002 | TWO-DP | ipSAE_d0res (max) 0.769888 | `ipSAE` | 0.76 | truncate -> 0.76 | 0.0000 | round -> 0.77 | 0.0100 | PASS |
| FX-002 | TWO-DP | pDockQ 0.687671 | `pDockQ` | 0.68 | truncate -> 0.68 | 0.0000 | round -> 0.69 | 0.0100 | PASS |
| FX-002 | TWO-DP | pDockQ2 (max) 0.728214 | `pDockQ2` | 0.72 | truncate -> 0.72 | 0.0000 | round -> 0.73 | 0.0100 | PASS |
| FX-002 | TWO-DP | LIS (max, AFDB convention) 0.659194 | `LIS` | 0.65 | truncate -> 0.65 | 0.0000 | round -> 0.66 | 0.0100 | PASS |
| FX-006 | FULL | ipSAE_d0res (max) 0.705718 | `ipSAE` | 0.71 | round -> 0.71 | 0.0000 | truncate -> 0.70 | 0.0100 | PASS |
| FX-006 | FULL | pDockQ 0.145233 | `pDockQ` | 0.15 | round -> 0.15 | 0.0000 | truncate -> 0.14 | 0.0100 | PASS |
| FX-006 | FULL | pDockQ2 (max) 0.705404 | `pDockQ2` | 0.71 | round -> 0.71 | 0.0000 | truncate -> 0.70 | 0.0100 | PASS |
| FX-006 | FULL | LIS (max, AFDB convention) 0.608776 | `LIS` | 0.61 | round -> 0.61 | 0.0000 | truncate -> 0.60 | 0.0100 | PASS |
| FX-007 | FULL | ipSAE_d0res (max) 0.711961 | `ipSAE` | 0.71 | round -> 0.71 | 0.0000 | truncate -> 0.71 | 0.0000 | PASS |
| FX-007 | FULL | pDockQ 0.273730 | `pDockQ` | 0.27 | round -> 0.27 | 0.0000 | truncate -> 0.27 | 0.0000 | PASS |
| FX-007 | FULL | pDockQ2 (max) 0.591426 | `pDockQ2` | 0.59 | round -> 0.59 | 0.0000 | truncate -> 0.59 | 0.0000 | PASS |
| FX-007 | FULL | LIS (max, AFDB convention) 0.417852 | `LIS` | 0.42 | round -> 0.42 | 0.0000 | truncate -> 0.41 | 0.0100 | PASS |
| FX-008 | FULL | ipSAE_d0res (max) 0.636590 | `ipSAE` | 0.64 | round -> 0.64 | 0.0000 | truncate -> 0.63 | 0.0100 | PASS |
| FX-008 | FULL | pDockQ 0.241786 | `pDockQ` | 0.24 | round -> 0.24 | 0.0000 | truncate -> 0.24 | 0.0000 | PASS |
| FX-008 | FULL | pDockQ2 (max) 0.325409 | `pDockQ2` | 0.33 | round -> 0.33 | 0.0000 | truncate -> 0.32 | 0.0100 | PASS |
| FX-008 | FULL | LIS (max, AFDB convention) 0.496848 | `LIS` | 0.50 | round -> 0.50 | 0.0000 | truncate -> 0.49 | 0.0100 | PASS |
| FX-009 | FULL | ipSAE_d0res (max) 0.711865 | `ipSAE` | 0.71 | round -> 0.71 | 0.0000 | truncate -> 0.71 | 0.0000 | PASS |
| FX-009 | FULL | pDockQ 0.719113 | `pDockQ` | 0.72 | round -> 0.72 | 0.0000 | truncate -> 0.71 | 0.0100 | PASS |
| FX-009 | FULL | pDockQ2 (max) 0.289935 | `pDockQ2` | 0.29 | round -> 0.29 | 0.0000 | truncate -> 0.28 | 0.0100 | PASS |
| FX-009 | FULL | LIS (max, AFDB convention) 0.440766 | `LIS` | 0.44 | round -> 0.44 | 0.0000 | truncate -> 0.44 | 0.0000 | PASS |
| FX-010 | FULL | ipSAE_d0res (max) 0.688783 | `ipSAE` | 0.69 | round -> 0.69 | 0.0000 | truncate -> 0.68 | 0.0100 | PASS |
| FX-010 | FULL | pDockQ 0.580514 | `pDockQ` | 0.58 | round -> 0.58 | 0.0000 | truncate -> 0.58 | 0.0000 | PASS |
| FX-010 | FULL | pDockQ2 (max) 0.272080 | `pDockQ2` | 0.27 | round -> 0.27 | 0.0000 | truncate -> 0.27 | 0.0000 | PASS |
| FX-010 | FULL | LIS (max, AFDB convention) 0.421410 | `LIS` | 0.42 | round -> 0.42 | 0.0000 | truncate -> 0.42 | 0.0000 | PASS |

**28 of 28 comparisons pass** under the regime-aware rule. The "other convention" columns are
printed deliberately: they are what a naive comparison would have produced.

### 5.4 A correction to the R060 methodology note

The R060 outcome note in `rework-plan.md` says, in bold:

> R092 must compare against the full-precision per-direction fields where they exist, and must
> **truncate rather than round** when only 2 dp is stored.

**The first half is right and load-bearing. The second half is right only for the TWO-DP
regime.** Applied uniformly it is wrong, and it fails loudest on exactly the entries R060 was
trying to protect.

The evidence, all of it in the tables above:

- **TWO-DP entries truncate.** FX-001 `pDockQ2` 0.926932 is stored as 0.92, not 0.93. FX-002
  `pDockQ` 0.687671 is stored as 0.68, not 0.69; `pDockQ2` 0.728214 as 0.72, not 0.73; `LIS`
  0.659194 as 0.65, not 0.66; `ipSAE` 0.769888 as 0.76, not 0.77. Five separate values where the
  two conventions differ, and all five say truncate. `pDockQ` is the decisive one, because it has
  no per-direction field, so "max over already-truncated per-direction values" cannot explain it.
- **FULL entries round their rolled-up fields.** FX-006 `ipSAE` 0.705718 is stored as 0.71, which
  truncation cannot produce. So are FX-006 `pDockQ` 0.145233 -> 0.15, FX-006 `pDockQ2` 0.705403
  -> 0.71, FX-006 `LIS` 0.608774 -> 0.61, FX-008 `ipSAE` 0.636590 -> 0.64, FX-008 `pDockQ2`
  0.325409 -> 0.33, FX-008 `LIS` 0.496848 -> 0.50, FX-009 `pDockQ` 0.719113 -> 0.72, FX-009
  `pDockQ2` 0.289935 -> 0.29, FX-010 `ipSAE` 0.688783 -> 0.69. Ten values where the conventions
  differ, and all ten say round.

Applying R060's rule uniformly across this fixture set produces **11 false failures**, each of
exactly 0.01: FX-006 on all four rolled-up values, FX-008 on three, FX-009 on two, FX-007 and
FX-010 on one each. Applying the regime-aware rule produces **zero**.

For completeness, the other two ways to manufacture a false failure against AFDB on this same
fixture set, both of which the first pass of this task hit and both of which are comparison
errors rather than numerical ones: comparing `ipTM_d0chn` against `complexPredictionAccuracy_ipTM`
fails on 6 of 7 fixtures (Section 5.7), and comparing our reported LIS mean against AFDB's
rolled-up LIS max fails on all 5 heterodimers (Section 5.5). Together with the 11 above, that is
22 apparent failures, none of them a disagreement about a number.

The rule that survives contact with all seven fixtures is:

> Detect the regime from whether the per-direction score fields carry more than 2 dp. On a
> TWO-DP entry, truncate everything. On a FULL entry, compare the per-direction fields exactly
> and **round** the rolled-up ones. Never compare a rolled-up field when a per-direction field
> for the same quantity exists.

This is a refinement of R060's finding, not a contradiction of it: R060 measured only FX-001 and
FX-006, and on FX-006 it only examined per-direction fields, where the question does not arise.
**R093 should carry the corrected rule.**

### 5.5 LIS is combined as the mean — verified, and the difference from AFDB is a convention

The table below carries both of the checks the task asked for explicitly: the pDockQ symmetry
flag (Section 5.6) and the two LIS directions beside the mean the notebook reports.

| Fixture | pDockQ symmetric (asserted, 1e-12) | LIS A->B | LIS B->A | LIS reported (MEAN) | max, for contrast |
|---|---|---:|---:|---:|---:|
| FX-001 | True | 0.756442 | 0.756452 | 0.756447 | 0.756452 |
| FX-002 | True | 0.659194 | 0.659179 | 0.659187 | 0.659194 |
| FX-006 | True | 0.608776 | 0.592070 | 0.600423 | 0.608776 |
| FX-007 | True | 0.417852 | 0.378096 | 0.397974 | 0.417852 |
| FX-008 | True | 0.496848 | 0.459825 | 0.478336 | 0.496848 |
| FX-009 | True | 0.440766 | 0.414833 | 0.427799 | 0.440766 |
| FX-010 | True | 0.421410 | 0.346870 | 0.384140 | 0.421410 |

The reported LIS is the arithmetic mean of the two directions on all seven fixtures, matching
`ipsae_v4.py:982`, which forms `LIS_Score = (LIS[c1][c2] + LIS[c2][c1]) / 2.0` for its combined
row. The notebook says so on the face of the output: *"LIS : 0.7564 [HIGH] (the MEAN of the two
directions, not the max)"*.

On a homodimer the two conventions are indistinguishable (FX-001 mean 0.756447 vs max 0.756452).
On the heterodimers the gap is real and large enough to change a printed 2 dp figure: FX-007 mean
0.397974 vs max 0.417852, FX-010 mean 0.384140 vs max 0.421410 — 0.037 apart, 37x the tolerance.

**AFDB's rolled-up `complexPredictionAccuracy_LIS` is the max.** In the table in Section 5.3 the
LIS row is therefore compared against our max, not our reported mean, and it agrees on all seven.
This is a documented convention difference between AFDB's roll-up and `ipsae.py`'s combined row,
**not a discrepancy**, and the notebook is on `ipsae.py`'s side of it, which is what D2 requires.
The per-direction `LIS_AB` / `LIS_BA` fields — the ones that carry no convention at all — agree
to <= 2.4e-6 (Section 5.2), which is the comparison that actually settles the question.

### 5.6 pDockQ is symmetric — verified, not assumed

`compute_pdockq` computes the value in both orientations and asserts they agree to 1e-12 before
returning (`complex_interface_utils.py`, `_pdockq_from(contacts.contact_mask.T, ...)` followed by
`assert abs(reverse[0] - score) < 1e-12 and reverse[1] == n_pairs`). That assertion executed and
held on all seven fixtures, including the 8.16:1 and 1:7.11 length ratios where an axis mix-up
would be most likely to show. `ipsae.py` agrees: its two asymmetric rows print an identical
`pDockQ` on every fixture, and its combined row prints it without a max
(`ipsae_v4.py:989`).

### 5.7 `complexPredictionAccuracy_ipTM` is not ipSAE's `ipTM_d0chn`

This is worth stating plainly because it is the single easiest way to manufacture a false failure
against AFDB, and the first pass of this task did exactly that.

| Fixture | AFDB `_ipTM` | AFDB `_iptm_af` | our ipTM_d0chn | AFDB `_ipsae_iptm_d0chn` max | gap `_ipTM` vs ours |
|---|---:|---:|---:|---:|---:|
| FX-001 | 0.95 | 0.95 | 0.952902 | 0.95 | 0.0029 |
| FX-002 | 0.82 | 0.82 | 0.813102 | 0.81 | 0.0069 |
| FX-006 | 0.78 | 0.7849 | 0.772954 | 0.772954 | 0.0070 |
| FX-007 | 0.7 | 0.7019 | 0.668721 | 0.668721 | 0.0313 |
| FX-008 | 0.85 | 0.8529 | 0.818349 | 0.818349 | 0.0317 |
| FX-009 | 0.61 | 0.615 | 0.576863 | 0.576863 | 0.0331 |
| FX-010 | 0.65 | 0.6497 | 0.595829 | 0.595829 | 0.0542 |

`complexPredictionAccuracy_ipTM` tracks `complexPredictionAccuracy_iptm_af` — AlphaFold's own
ipTM, the quantity `ipsae.py` prints in its `ipTM_af` column and leaves at 0.000 for AFDB inputs
because the AF3 summary file is not published. It is a **different measurement** from ipSAE's
`ipTM_d0chn`, which is derived from the PAE matrix with a chain-pair `d0`. On FX-010 they are
0.054 apart; on FX-009, 0.033; on FX-008, 0.032. Substituting `_ipTM` for `_ipsae_iptm_d0chn_*`
produces a false failure on 6 of the 7 fixtures; only FX-001 happens to agree at 2 dp.

The correct AFDB reference for `ipTM_d0chn` is `complexPredictionAccuracy_ipsae_iptm_d0chn_AB` /
`_BA`, which this notebook matches to <= 5e-7 on all five full-precision fixtures (Section 5.2).
**R093 should record this distinction** — nothing in the current spec pack warns about it, and
`fixture-manifest.md`'s "Independent Reference Agreement" section names `_ipsae_iptm_d0chn_*`
correctly but does not say why `_ipTM` must not be substituted for it.

### 5.8 Values outside +-0.001

**None.**

- vs `ipsae.py` v4: 49 reported-value comparisons + 84 per-direction comparisons = 133. Zero
  failures. Worst |delta| 4.8e-5, and every delta above 1e-6 is bounded by the reference's own
  print precision.
- vs AFDB, regime-aware: 70 per-direction + 28 rolled-up = 98. Zero failures. Worst |delta|
  2.4e-6 on the full-precision fields; exact on the 2 dp fields. (`ipSAE_d0dom` accounts for the
  14 per-direction slots AFDB does not publish.)

The two references disagree with each other nowhere either, which is the strongest single
statement available: two independent implementations and this notebook produce the same numbers.

---

## 6. Negative paths, verbatim

These are the messages a user sees. Each is the full text of the raised exception; in every case
execution halts at that cell under Run All, and no scoring cell runs.

### 6.1 FX-004, malformed accession (`AF-NOT_A_REAL_ACCESSION`) — halts at cell 6, `AccessionLookupError`

```
AFDB rejected the accession 'AF-NOT_A_REAL_ACCESSION' as malformed (HTTP 400).
  found     : the endpoint refused the identifier before looking anything up.
              the service said: Invalid identifier format. Please use a UniProt accession or a supported AlphaFold DB ID.
  expected  : an AFDB model accession -- 'AF-0000000065889468' (a two-chain complex),
              'AF-P0A6Q3-F1' (one UniProt entry), or a bare UniProt accession
              such as 'P0A6Q3'.
  what to do: check ACCESSION_ID = 'AF-NOT_A_REAL_ACCESSION' in the input cell for a typo. To find
              an accession this notebook supports:
                https://alphafold.ebi.ac.uk/api/search?q=oligomericState:dimer&type=main&rows=5
              To analyse a structure of your own instead, set USE_LOCAL_FILE = True.
```

Raised before any download. The service's own message is quoted rather than paraphrased.

### 6.2 FX-005, monomer (`AF-O15552-F1`) — halts at cell 6, `UnsupportedAssemblyError`

```
AF-O15552-F1 is a monomer; this notebook analyses two-chain dimers only.
  found     : the AFDB metadata declares isComplex=false,
              and describes 1 chain: A
  supported : exactly 2 chains (oligomericState 'dimer'), homodimer or
              heterodimer alike.
  why       : ipTM, ipSAE, pDockQ, pDockQ2 and LIS are all *inter-chain*
              measurements -- they read the PAE block between two
              different chains and the contacts across the interface.
              With one chain there is no such block and no interface,
              so there is no number to report, not even a bad one.
  what to do: choose an AFDB *complex* accession. Every complex in AFDB
              is a dimer, and you can list some with:
                https://alphafold.ebi.ac.uk/api/search?q=oligomericState:dimer&type=main&rows=5
              A known-good homodimer: AF-0000000065889468.
              For per-residue confidence of O15552 on its own, the AFDB
              entry page https://alphafold.ebi.ac.uk/entry/O15552 already shows
              pLDDT and PAE; this notebook adds nothing for a single chain.
```

Raised at the declared-assembly gate, before the three document downloads. This confirms the
correction R023 recorded: the `CLAUDE.md` edge-case table still claims a monomer "produces
zero-length chain B; scores will be 0.0 (no informative error)", which is false at this revision.
R093 already owns that fix.

### 6.3 Local-file mode, mmCIF only — halts at cell 8, `MissingLocalDocumentError`

```
Local-file mode needs all three documents; 2 were not uploaded.
  uploaded  : mmCIF structure
  missing   : PAE document, pLDDT document
  why       : six of the seven values are computed from the PAE matrix
              (ipSAE d0res/d0chn/d0dom, ipTM_d0chn, pDockQ2, LIS) and the
              seventh, pDockQ, needs per-residue pLDDT. Without them
              Sections 3, 4, 6 and 7 have nothing to show and the traffic
              light has nothing to colour, so there is no partial run
              worth offering -- only a contact count and a contact map.
  what to do: upload the missing files in the widget above, then re-run this
              cell. For an AFDB model, the download links are the
              `cifUrl`, `paeDocUrl` and `plddtDocUrl` fields of
                https://alphafold.ebi.ac.uk/api/prediction/AF-0000000065889468
              and the files are named
                AF-<id>-predicted_aligned_error_v<n>.json
                AF-<id>-confidence_v<n>.json
              For a model of your own, export the PAE matrix and the
              per-residue pLDDT in the AFDB JSON layout that `parse_pae`
              and `parse_plddt` document.
              To analyse an AFDB entry instead, set USE_LOCAL_FILE = False
              and put its accession in ACCESSION_ID.
```

### 6.4 Local-file mode, mmCIF + PAE, no pLDDT — halts at cell 8, `MissingLocalDocumentError`

```
Local-file mode needs all three documents; one was not uploaded.
  uploaded  : mmCIF structure, PAE document
  missing   : pLDDT document
  ...
  what to do: upload the missing file in the widget above, then re-run this
              cell. ...
              and the files are named
                AF-<id>-confidence_v<n>.json
```

The count, the "uploaded" list, the "missing" list, the singular/plural of "the missing file(s)"
and the named filenames all track what is actually absent. R025's promise — refuse at upload
time, naming every missing file at once, rather than advertising a skip that never happens — is
kept in both cases.

### 6.5 Local-file mode, monomer documents — halts at cell 10, `UnsupportedAssemblyError`

```
This model is a monomer; this notebook analyses two-chain dimers only.
  found     : 1 chain in the structure and the PAE/pLDDT documents:
                A (330 residues)
              Monomer -- 1 chain (A), chain identities unknown
  supported : exactly 2 chains (oligomericState 'dimer'), homodimer or
              heterodimer alike.
  ...
```

Not required by the task, run because it is cheap and it exercises the one gate local-file mode
cannot skip. It fires from `verify_chain_identity`, three cells later than the online monomer
refusal, because in local mode there is no metadata to gate on — exactly as cell 6 states.

### 6.6 One observation about the negative paths

Under Run All the halt is clean and nothing downstream executes. Under a "run every cell and
ignore errors" driver — which is not a mode Jupyter or Colab offer by default, but is what
`nbclient` does with `allow_errors=True` — the 27 cells after the halt each raise their own
`NameError` on the variable the failed cell would have defined. That is inherent to a linear
notebook and is not a defect in this one; it is recorded so that a future automated run does not
mistake a cascade of `NameError`s for 28 separate bugs.

---

## 7. Local-file mode reproduces the online path exactly

Run 10a used FX-001's three documents from disk, with the AFDB API never contacted. All seven
values are **bit-identical** to the online run, not merely equal to 4 dp:

| Value | Online (run 1) | Local file (run 10a) |
|---|---|---|
| ipSAE_d0res | 0.9143088732360954 | 0.9143088732360954 |
| ipSAE_d0chn | 0.9529017899515078 | 0.9529017899515078 |
| ipSAE_d0dom | 0.9526724165644239 | 0.9526724165644239 |
| ipTM_d0chn | 0.9529017899515078 | 0.9529017899515078 |
| pDockQ | 0.6912735609411721 | 0.6912735609411721 |
| pDockQ2 | 0.9269317982158931 | 0.9269317982158931 |
| LIS | 0.7564470270478529 | 0.7564470270478529 |

Contacts, interface residue counts and every figure matched too (116 pairs, 43 / 43). The only
difference is Section 6: no Mol* views, with the reason printed.

---

## 8. What was **not** verified

Stated plainly, because an unqualified "validated" would be wrong.

1. **Google Colab. Not tested at all.** Every run was local, on macOS/arm64, with
   `IN_COLAB = False`. That means three things went unexercised: the `git clone --depth 1 --branch
   homodimer-notebook-rework` path, the stale-clone refresh path (`git fetch` + `git reset --hard
   FETCH_HEAD`), and the conditional `pip install molviewspec`. The 60 s Colab install budget in
   `CLAUDE.md` is therefore **unmeasured**. Local wall time per fixture was 3.9-5.3 s including all
   downloads, which is consistent with the target but is not evidence for it. This is the largest
   remaining gap and it needs a human with a browser.
2. **The `TODO(merge)` in the bootstrap cell is still live.** `REPO_BRANCH` is pinned to
   `homodimer-notebook-rework`. A Colab user running the published notebook today would clone a
   branch, not a release tag. R093 owns this; it is called out here because a Colab test that
   passes before the merge proves nothing about a Colab run after it.
3. **Interactive widget behaviour.** The upload widgets were driven by a stub subclass, so the
   real click-and-upload path, the widget's own rendering and the browser file picker are
   untested. What was tested is everything downstream of the bytes arriving.
4. **The three-chain local case** from `fixture-manifest.md`'s local-file table. It needs a
   synthetic mmCIF and is not one of the ten registered fixtures this task was scoped to.
5. **Mol\* rendering.** Six `<iframe>` elements were emitted per online run and the builders ran
   without error, but nothing confirms the 3D scene inside them draws correctly. That needs eyes
   on a browser.
6. **`ipSAE_d0dom` against AFDB.** AFDB publishes the `d0dom` and `n0dom` intermediates but not
   the score, so this one value has only `ipsae.py` behind it, not two references. Its
   intermediates match AFDB exactly per direction, which is strong circumstantial support, but it
   is not a second independent check of the score itself.
7. **Other Python and numpy versions.** One environment only (3.11.9 / numpy 2.3.3).

---

## 9. Findings and outstanding risk

No bug was found in the notebook, the module, or `interface.py`. Per the task's constraints
nothing outside this file was modified. The items below are for later tasks.

| # | Finding | Severity | Owner |
|---|---|---|---|
| 1 | R060's "truncate, do not round" rule is correct only for the TWO-DP storage regime; the rolled-up 2 dp fields on full-precision AFDB entries are **rounded**. Uniform application produces 11 false failures on this fixture set. Corrected rule in Section 5.4. | methodology, not a code defect | R093 |
| 2 | `complexPredictionAccuracy_ipTM` is AlphaFold's own ipTM, not ipSAE's `ipTM_d0chn`; up to 0.054 apart on this set, and a false failure on 6 of 7 fixtures if substituted. Nothing in the spec pack warns against it. | documentation | R093 |
| 3 | FX-002 is scored at `ipSAE_d0res` 0.7699 and is **not** borderline. Its manifest status (`candidate-needs-scoring`) and its stated role are both stale. | documentation | manifest edit |
| 4 | The R080 Section 7 defect recorded as open in `fixture-manifest.md` does not reproduce; the prose now reads "mixed" and names the offending score. The manifest's "What is still open" entry is stale. | documentation | manifest edit |
| 5 | `CLAUDE.md`'s monomer edge case ("zero-length chain B; scores will be 0.0, no informative error") is contradicted by the measured behaviour in Section 6.2. | documentation | R093 (already listed) |
| 6 | Colab is entirely unverified, and the bootstrap still pins `REPO_BRANCH` to a feature branch. | **risk** | R093 + a human Colab run |
| 7 | `fixture-manifest.md` records `pDockQ` as "unreconciled at full precision" because AFDB publishes only 2 dp. **Closed by this report**: `ipsae.py` gives it at 4 dp and agrees to <= 3.3e-5 on all seven fixtures. | closes an open item | manifest edit |

### The one risk worth restating

Everything numerical here rests on a single environment and a single day's live AFDB responses.
The scores are reproducible from cached documents; the AFDB *reference* values are not — they are
whatever the search endpoint returned on 2026-09-08. If AFDB reprocesses an entry, reference (b)
moves and this report's deltas go stale while reference (a) does not. That is an argument for
treating `ipsae.py` as the primary reference, as D2 already does, and AFDB as the corroborating
one.

---

## 10. Verification steps, as run

1. **Execution matrix** — Section 3. Ten registered fixtures plus three local-file cases; 13 runs.
2. **Agreement tables against both references, with the precision regime per fixture** —
   Sections 5.1 (ipsae.py), 5.2 and 5.3 (AFDB, regime stated per fixture and per table row).
3. **Any value outside +-0.001, with diagnosis** — Section 5.8: none, out of 231 comparisons
   (133 against `ipsae.py`, 98 against AFDB). The 22 apparent failures a naive comparison
   produces are diagnosed in Sections 5.4, 5.5 and 5.7, and every one of them is an artefact of
   the comparison rule rather than of the notebook's numbers.
4. **Negative-path messages** — Section 6, verbatim, all five.
5. **`git status --porcelain`** — run after all work. The only entries are the
   `specs/nmr_restraints/*`, `notebooks/protein_model_chem.ipynb`,
   `notebooks/nmr_restraints_visualization.ipynb` and `.claude/` changes that were already
   present at the start of this session and belong to unrelated work. `git diff HEAD` over
   `notebooks/homodimer_diagnostic.ipynb`, `src/insightfold/`, `CLAUDE.md` and
   `specs/homodimer_diagnostic/` is **empty**: nothing this task touched, other than adding this
   file.

## 11. Provenance

- Notebook and module: `homodimer-notebook-rework @ e162059`.
- Reference (a): `specs/homodimer_diagnostic/references/ipsae_v4.py`, sha256
  `10cf9b08c68c91e06cb28526cf2026f47a3980c9048fd3226d13e3304eaf1c27`, invoked as
  `python ipsae_v4.py <af3_pae.json> <model.cif> 10 8`.
- Reference (b): `https://alphafold.ebi.ac.uk/api/search?q=modelEntityId:<accession>&type=main&rows=1`,
  retrieved 2026-09-08.
- Fixture documents: `cifUrl`, `paeDocUrl`, `plddtDocUrl` from
  `https://alphafold.ebi.ac.uk/api/prediction/<accession>`, retrieved 2026-09-08.
- Tolerance: +-0.001 throughout, per R092.
