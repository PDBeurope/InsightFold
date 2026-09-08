# Homodimer Diagnostic Notebook — Rework Plan

Source of truth for the April-2026 rework of `notebooks/homodimer_diagnostic.ipynb`.
Every task below traces to a note in the user's "Planned notebook fixes" list or to a
defect found while auditing the notebook against `ipsae.py` v4 and the live AFDB API.

Status legend: `[ ]` not started · `[~]` in progress · `[x]` done · `[!]` blocked

Related documents:
- `specs/homodimer_diagnostic/tasks.md` — the original build plan (Phases 1-7). Still valid
  history; this file supersedes it for all work from 2026-09-07 onward.
- `specs/homodimer_diagnostic/fixture-manifest.md` — fixture registry, extended by R024.
- `CLAUDE.md` — project skill file. Several of its formulas are wrong; corrected by R093.

---

## Locked decisions

| # | Decision | Source |
|---|----------|--------|
| D1 | All extracted utilities live in a single flat module named `complex_interface_utils.py`. | user |
| D2 | Numerical ground truth is `DunbrackLab/IPSAE` `ipsae.py` v4 (Jan 2026). Publications are used for the explanatory prose and threshold rationale only. | user |
| D3 | Homodimers and heterodimers are both supported automatically, with no user-facing parameter to switch between them. Heterodimer support applies everywhere in the notebook, not just the scoring cells. | user |
| D4 | Architecture stays dimer-scoped for now, but APIs take an explicit chain pair so an N-chain generalisation is a later extension, not a rewrite. | user |
| D5 | Each task is verified by running the affected cells locally against a live AFDB fetch. Outputs are stripped before commit. | user |
| D6 | Notebook filename stays `homodimer_diagnostic.ipynb`. | user |
| D7 | `complex_interface_utils.py` lives in `src/insightfold/`. | user |
| D8 | Modules are named by **domain**, not by notebook, so several notebooks can share one. A function is promoted out of a domain module into a shared one when its **second** consumer appears, not in advance. | user |
| D9 | Colab bootstrap is `git clone --depth 1` plus `sys.path.insert(0, root/'src')`. **Not** `pip install`, which would resolve `pyproject.toml` and drag in biopython, gemmi, torch and plotly, breaking both the prohibited-dependency rule and the 60 s install budget. | user's other notebook, adapted |
| D10 | Redundant modules and duplicated inline copies are identified now but **deleted only after all work is verified**, so an early deletion cannot masquerade as a later bug. | user |

**Verified facts behind D9.** `PDBeurope/InsightFold` is public: the GitHub API returns 200
unauthenticated and `raw.githubusercontent.com` serves from `main`. No token, no `getpass`,
no auth header. This matters beyond convenience, because `getpass` inside a Run-all notebook
blocks forever on input nobody is watching.

## Assumptions I made (flag if wrong)

- **A2 — Em dash scope.** The sweep covers markdown cells, plot titles/labels, and printed
  output. Code comments are left alone. En dashes inside numeric ranges (`70–90`) are not
  em dashes and stay.
- **A3 — Zero-contact defaults.** `ipsae.py` returns `0.0` for pDockQ and pDockQ2 when there
  are no contacts; the notebook and `CLAUDE.md` return the sigmoid minima (0.018 / 0.005).
  I will follow `ipsae.py` and document the difference. See R006.

## Reference publications

Supplied by the user, stored in `specs/homodimer_diagnostic/references/`, which is
**gitignored** — these are publisher PDFs and must not be redistributed from a public repo.

| File | Paper | Needed by |
|------|-------|-----------|
| `ipsae_dunbrack_2025.pdf` | ipSAE, bioRxiv 2025.02.10.637595v2 | R002, R061 (rationale for the three `d0` variants, thresholds) |
| `lis_kim_2024.pdf` | LIS, bioRxiv 2024.02.19.580970v1 | R002 (widest threshold disagreement in the repo) |
| `pdockq2_zhu_2023.pdf` | pDockQ2, Bioinformatics btad424 | R002, R061 |
| `pdockq_bryant_2022.pdf` | pDockQ, Nat Commun s41467-022-28865-w | R061 (the 0.23 cutoff is already agreed across sources) |
| `ipsae_confidence_bands_report.pdf` | **AFDB internal technical note, "ipSAE Confidence Bands for AlphaFold DB Homodimers", March 2026, 9 pp.** Supplied by the user 2026-09-08. Literature review behind the AFDB bands, with per-band rationale, a cross-reference table to ipTM and DockQ tiers, ready-made user-facing FAQ text, and a caveats section. Feeds R002b, R061, R072, R080. |
| `afdb_dimers_2026.pdf` | **Han, Tsenkov, Venanzi et al., "AlphaFold Database expands to proteome-scale quaternary structures", bioRxiv 10.64898/2026.03.27.714458v2, posted 3 Jul 2026** | R009, R080. The authoritative source for the AFDB confidence criteria. The user is second author. |

Also: EBI news item, "Millions of protein complexes added to AlphaFold Database shed light on
how proteins interact", 16 Mar 2026,
`https://www.ebi.ac.uk/about/news/technology-and-innovation/first-complexes-alphafold-database/`
— a lay-audience companion to the preprint. Cites 1.7M high-confidence homodimers and states
no thresholds. Usable as a friendly link for notebook readers, not as a technical citation.

## Open questions

**Q4 (2026-09-08) — PAE cutoff: RESOLVED by decision, caveat recorded.**

`PAE_CUTOFF` stays at **10.0**, as in `CLAUDE.md` and the notebook today (user, 2026-09-08:
use the original thresholds).

Recorded so R092 is not surprised by it. The AFDB internal report §7 notes Dunbrack recommends
15 Å for AlphaFold2 and 10 Å for AF3/Boltz-1, and AFDB dimers are AlphaFold-Multimer via
ColabFold. The preprint's Methods (p.19) do not state which cutoff the production C++ ipSAE
implementation used, and the API exposes no ipSAE value to infer it from. Measured effect of
the choice on the fixtures:

| Fixture | @ 10 | @ 15 | delta |
|---------|------|------|-------|
| `AF-0000000065889468` homodimer | 0.9143 | 0.9143 | 0.0000 |
| `AF-0000000211034637` heterodimer | 0.7057 | 0.6981 | 0.0076 |

Consequence to state once in the notebook, not to act on: a value computed here may differ
slightly from the one behind an AFDB entry's badge if the production pipeline used a different
cutoff, and near a band boundary that difference can change the band. R092 compares against
`ipsae.py` at cutoff 10, which is self-consistent and is what the tolerance applies to.

---

## W0 — Correctness ground truth

Blocking for W6 and W8. The audit is already done; these tasks are about applying it.

### `[x] R001 — Publish the verified formula reference` — DONE 2026-09-07

**Delivered.** `specs/homodimer_diagnostic/formula-reference.md`.

**That file is now the single source of truth for all seven formulas.** The draft audit table
that previously lived here has been deleted deliberately: keeping a second, less accurate copy
would recreate exactly the multiple-sources-of-truth problem this rework exists to fix.

**Verification.** The agent re-derived every claim from `references/ipsae_v4.py` rather than
transcribing the draft. The load-bearing citations were then independently confirmed against
source: `:977` (pDockQ2 max over directions), `:982` (LIS mean), `:746`/`:747` (ipTM vs
ipSAE_d0chn), `:116-136` (the two `d0` helpers), `:751-756` (the `n0dom` sets), `:669`/`:700`
(zero-contact returns), plus `:111-113`, `:385`, `:534`, `:851`, `:989`.

**Outcome.** All six flagged draft claims confirmed correct. One substantive omission found,
plus three citation defects and three completeness gaps, all recorded in the document's
"Corrections to the draft audit" section. The substantive omission is promoted to R007 below.

### `[ ] R002b — Use the AFDB bands report to fill prose gaps only`

**Scope, set by the user 2026-09-08: the thresholds are settled and do not change.** The
report is a source for the explanatory material that is currently missing, and for nothing
else. No number in `threshold-reference.md` or `THRESHOLDS` moves as a result of this task.

**The agreed thresholds, unchanged and final:**

| Value | Green | Amber |
|-------|-------|-------|
| ipSAE (all three variants) | 0.70 | 0.60, with AFDB's published bands >=0.8 / 0.7-0.8 / 0.6-0.7 / <0.6 |
| ipTM_d0chn | 0.70 | 0.30 |
| pDockQ | 0.23 | 0.12 |
| pDockQ2 | 0.23 | 0.10 |
| LIS | 0.21 | 0.10 |
| Joint AFDB rule | ipSAE_d0res >= 0.6 AND pDockQ2 >= 0.23 | |

**What to take from the report — prose only:**
- Per-band rationale for each ipSAE tier, with its citations (report §4).
- The user-facing FAQ text in AFDB's house style (report §6). R072 and R080 follow its
  wording rather than inventing parallel phrasing.
- The cross-reference table aligning ipSAE, ipTM and DockQ tiers (report §5), as explanatory
  context for a reader who knows one scale and not the others.
- **The homodimer caveat (report §3.9), which belongs beside the traffic light:** AlphaFold's
  false-positive rate is markedly worse for homodimers than heterodimers, with the
  true-positive rate at 1% FPR dropping from 63% to 18%, because proteins that do not
  homodimerise in vivo often have homologs that do.
- The report's own caveats (§7): the bands are "evidence-informed rather than statistically
  optimised on a dedicated homodimer validation set", ipSAE is not experimentally calibrated
  for homodimers, and near-boundary scores (0.59 vs 0.61) must not be over-interpreted.

**Recorded negative result, so it is not re-investigated.** Kim 2024 was re-searched in full
for a second LIS tier. There is none: one Youden-optimal cutoff per metric (best LIS 0.21,
ipTM 0.38, p.30) and no gradation below it. The nearby numbers 0.911 and 0.891 are **AUC
values, not thresholds** (p.5) and must never be quoted as cutoffs. The LIS amber of 0.10
stays a judgement call and is the only one in the table.

**Not adopted.** The DockQ CAPRI ladder (0.23 / 0.49 / 0.80, report §3.3) would also transfer
to pDockQ and pDockQ2 by the same sigmoid-fitting argument that licenses 0.23. It is recorded
here as context a reader may find useful, and explicitly **not** used to change any threshold.

**Done when.** The report's rationale, FAQ wording, cross-reference table and caveats are
recorded for R061/R072/R080 to use, and no threshold has changed.

### `[x] R007 — Fix pDockQ2 direction handling` — DONE 2026-09-08

**What.** Compute pDockQ2 for both directions and report `max(A→B, B→A)`.

**Why.** Found by the R001 independent check; missed by the original audit. `ipsae.py`
computes `mean_ptm` from the directional PAE block only (`ipsae_v4.py:686`) and reports
`max(pDockQ2[A][B], pDockQ2[B][A])` (`ipsae_v4.py:977, 990`). The notebook calls
`compute_pdockq2(dist_matrix, pae_AB, ...)` in cell 19 and reports the A→B value alone.

**Supporting notes.** Same class of defect as the `d0dom` asymmetry in R003, and it will be
invisible on a homodimer test while being wrong on a heterodimer, which is precisely the gap
W2 exists to close. `mean_plddt` is symmetric and needs no change; only `mean_ptm` is
directional. Land alongside R003 in the same scoring pass, and give R074 the per-residue
`mean_ptm` intermediate it needs while the function is open.

**Done when.** pDockQ2 matches `ipsae.py` within ±0.001 on both fixtures, including the
heterodimer where the two directions genuinely differ.

### `[x] R002 — Establish one canonical threshold table` — DONE 2026-09-07

**Delivered.** `specs/homodimer_diagnostic/threshold-reference.md`. **Single source of truth
for all thresholds.** The four-way conflict table that previously lived here is deleted; see
that file's "Changes this file makes to the repo's current numbers" section instead.

**Canonical values.** ipSAE (all three variants) 0.60/0.30 · ipTM_d0chn 0.70/0.30 ·
pDockQ 0.23/0.12 · pDockQ2 0.23/0.10 · LIS 0.21/0.10.

**Verification performed.** Both PUBLISHED citations confirmed verbatim from the PDFs:
Kim 2024 p.30 "best LIS values exceeding the optimal threshold (0.21)"; Zhu 2023 p.6
"9 had pTM > 0.5 and min pDockQ2 > 0.23". The negative claim that Dunbrack states no ipSAE
cutoff was confirmed by full-text search: **zero** matches for `thresh*`, `cutoff of 0.`,
`ipSAE [><] 0.`, `above 0.x`, `below 0.x` across all 24 pages. The empirical ipSAE table was
reproduced independently against the live AFDB PAE endpoints and matched to three decimals.

**Key findings.**
1. **Only two of the fourteen numbers are genuinely published** (pDockQ2 green, LIS green).
   Six are outright judgement calls, collected in the document's "Values with no literature
   basis" section. This is now visible rather than implied.
2. **`ipSAE >= 0.6` is not Dunbrack's number.** It is an AFDB *dataset release filter*.
   Because every AFDB homodimer has already passed it, the threshold is near-tautological on
   AFDB input and will almost never show red. The notebook must say so.
3. **`ipSAE_d0chn >= ipSAE_d0dom >= ipSAE_d0res` is a theorem, not an observation** (d0 is
   monotone in L, ptm is increasing in d0, and `n0res_i <= n0dom <= nA+nB`). Verified
   pointwise on both fixtures. Consequence: **cell 26's per-variant thresholds have the bias
   backwards** — they give `d0chn`, the provably most permissive variant, the *lowest* green.
   One shared threshold set is used instead, with the ordering caveat surfaced to the user.
4. LIS green moves 0.15/0.3 -> 0.21; pDockQ2 green moves 0.15/0.5 -> 0.23. Neither of the
   old numbers appears in either paper. `CLAUDE.md`'s low LIS green was likely compensating
   for its own max-instead-of-mean bug.
5. pDockQ2's sigmoid is bounded on [0.005, **1.315**] and can exceed 1.0. It must never be
   presented as a probability or a percentage.

### `[x] R009 — Adopt the AFDB joint criterion and published confidence bands` — DONE 2026-09-08

**What.** Fold the AFDB preprint into `threshold-reference.md` and the module's `THRESHOLDS`,
replacing invented numbers with AFDB's published scheme. Decided with the user 2026-09-08.

**Verified quotes from `references/afdb_dimers_2026.pdf`** (extracted and checked against the
PDF text; page numbers are PDF pages):

- p.5 — the criterion and its provenance:
  > "We adopted a combined high-confidence criterion requiring both community-established
  > cutoffs of ipSAEmax >= 0.6[32,34,35] and pDockQ2max >= 0.23, corresponding to the DockQ
  > 'acceptable' quality boundary[36]. This joint cutoff yielded a precision of 0.924 (False
  > Positive Rate (FPR)= 0.043) for homodimers and 0.958 (FPR=0.004) for heterodimers
  > (Fig. 2a,b), supporting its use as a quality filter that prioritises precision over recall
  > given the scale of the release."
- p.12 — the published confidence bands:
  > "surfaced entries are further categorised in 'very high-confidence' (ipSAEmax >= 0.8),
  > 'confident' (0.7 <= ipSAEmax < 0.8), and 'low-confidence' (0.6 <= ipSAEmax < 0.7)."
- p.11 — the headline scale:
  > "we compiled over 31 million candidate complexes and surfaced 1.81 million high-confidence
  > assemblies."
- p.12 — below-threshold models are still available:
  > "Dimers not passing the previously defined threshold, together with their interface scores,
  > are provided on the FTP page ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/."

**Changes to make.**
1. **Joint AFDB decision rule as a binary badge.** `ipSAE_d0res >= 0.6 AND pDockQ2 >= 0.23`
   -> PASS / FAIL, presented as "would this model qualify for AFDB high-confidence release?".
   This is the actual AFDB classifier; seven independent traffic lights are not.
2. **`ipSAE_d0res` adopts AFDB's four published bands** in place of the invented 0.60/0.30
   green/amber: `>= 0.8` very high-confidence, `0.7-0.8` confident, `0.6-0.7` low-confidence,
   `< 0.6` below the AFDB threshold. This removes the largest remaining judgement call.
3. **Provenance upgrade.** `ipSAE_d0res` green moves from `AFDB-EMPIRICAL` to `PUBLISHED`,
   cited to Han/Tsenkov 2026 rather than to Dunbrack. `pDockQ2` 0.23 is now double-sourced
   (Zhu 2023 directly; AFDB 2026 independently).
4. **Corroborations to record, not re-derive.** The paper's Supplementary Fig. 5 compares max-
   and min-chain-level aggregation of ipSAE and pDockQ2, confirming that the reported values
   are the **max over directions** — independent support for R007. The paper's own naming
   (`ipSAEmax`, `pDockQ2max`) says the same thing.

**Nuance that must reach the notebook, not be smoothed over.** The paper adopts 0.6 as
"community-established" by citation; it does not derive it. What it derives is the
*validation* (precision 0.924 / 0.958). Supplementary Figs. 1-2 report MCC-optimal cutoffs of
**0.104** for homodimers and **0.520** for heterodimers, far below 0.6, and the paper states
the joint cutoff "prioritises precision over recall given the scale of the release". So 0.6 is
a deliberately conservative operating point, not an optimum. That is a more accurate and more
useful thing to tell a notebook user than "the threshold is 0.6".

**Also correct the scale figures.** ~31 million candidate complexes compiled; 19,148,379
homodimers and 7,561,477 heterodimers analysed; 1,735,475 (9.1%) and 79,392 (1.0%)
high-confidence respectively; 1.81 million high-confidence assemblies in total.

**Citation policy — decided with the user 2026-09-08.** Separate the *metric* from the
*threshold*, and credit both.

- **ipSAE the metric is Dunbrack's, and must be cited as his** wherever the score is
  explained, defined or computed: Dunbrack, bioRxiv 2025.02.10.637595. Adopting a different
  threshold does not diminish authorship of the method. The earlier finding that the paper
  states no cutoff is a statement about *thresholds only*, and must never be presented as a
  criticism of the work or as a reason to under-credit it.
- **The 0.6 threshold is AFDB's**, cited to Han/Tsenkov et al. 2026. Do not attribute it to
  Dunbrack, which is the error `CLAUDE.md` currently makes by placing it under a
  DunbrackLab-attributed heading.
- **The EBI news item is cited too**, despite predating the v2 preprint, as the accessible
  entry point to the original work:
  `https://www.ebi.ac.uk/about/news/technology-and-innovation/first-complexes-alphafold-database/`
  Present it as background reading, with the preprint as the technical citation. Its 1.7M
  figure covers high-confidence homodimers only, so do not use it for the scale numbers.
- Same principle throughout: pDockQ cites Bryant, pDockQ2 cites Zhu, LIS cites Kim, and the
  `ipsae.py` implementation cites DunbrackLab/IPSAE v4, regardless of whose threshold is used.

**Confirmed for the notebook (user, 2026-09-08).** The conservative-not-optimal point is to be
stated explicitly: 0.6 is a deliberately conservative operating point chosen to prioritise
precision over recall at release scale, not an optimum. The MCC-optimal cutoffs on the paper's
own benchmarks are 0.104 (homodimer) and 0.520 (heterodimer).

**Done when.** `threshold-reference.md` and `THRESHOLDS` carry the AFDB bands and the joint
rule, every changed number cites the preprint by page, the metric-vs-threshold citation split
is written down, and the conservative-not-optimal nuance is recorded for R080 to surface.

### `[ ] R008 — Surface directional asymmetry as a diagnostic`

**What.** Wherever a score is computed per direction, show **both** directions and their
difference, rather than silently collapsing to one number. Requested by the user 2026-09-07.

**Conflict to resolve, and how.** The user asked for "A to B by default, with a parameter to
view B to A". Taken literally for the *reported* score this would break D2 and R092, because
`ipsae.py` does not report A→B — it reports `max(A→B, B→A)` for ipSAE and ipTM_d0chn and
pDockQ2, and the `mean` for LIS (`formula-reference.md` §§1-7). Reporting A→B would put every
headline value outside the ±0.001 agreement criterion on any asymmetric complex.

Resolution, which delivers the intent without breaking ground truth:
- **The headline score stays exactly as `ipsae.py` defines it** (max, or mean for LIS). This
  is the number in the summary table and the traffic light.
- **A `direction` parameter is added for inspection**, accepting `'AB'`, `'BA'`, or the
  default `None` meaning "the ipsae.py-reported combination". It changes what is *displayed
  and plotted*, never what is scored.
- **A directional-difference check runs automatically** and reports `|A→B - B→A|` per score,
  flagging any score whose two directions differ by more than a stated tolerance.

**Where it applies.** Directional: ipTM_d0chn, ipSAE_d0res, ipSAE_d0chn, ipSAE_d0dom,
pDockQ2, LIS. Also `n0dom`/`d0dom`, which differ per direction (R003). **Not applicable to
pDockQ**, which is provably symmetric — `npairs` and the interface residue set are both
invariant under swapping the chains (`formula-reference.md` §5), which is why `ipsae.py`
prints it unmaxed at `ipsae_v4.py:989`. The notebook should say pDockQ is symmetric rather
than showing an empty diff for it.

**Why this is worth doing.** Measured on the fixtures: the homodimer has `n0dom` 342 in both
directions, the heterodimer 250 vs 252. The asymmetry is real, it is invisible on homodimers,
and it is the mechanism behind both R003 and R007. A visible diff turns a class of silent bug
into something a user can see.

**Done when.** Every directional score displays both directions plus the delta; the headline
values still match `ipsae.py` within ±0.001; pDockQ is labelled symmetric rather than diffed.

### `[x] R003 — Fix d0dom directional asymmetry` — DONE 2026-09-08

**What.** Compute `n0dom` and `d0dom` separately for A→B and B→A, per `ipsae.py` L775-778.

**Supporting notes.** Current code computes one `n_dom` from `pae_ab` and passes it into
both `_per_res(pae_ab, nA)` and `_per_res(pae_ba, nB)`. Correct form per direction:
`n0dom(X→Y) = (block < cutoff).any(axis=1).sum() + (block < cutoff).any(axis=0).sum()`
where `block` is the X→Y PAE quadrant. Report the `d0dom` belonging to whichever direction
supplied the maximum, matching `ipsae.py`'s `n0dom_max` / `d0dom_max` logic (L866-880).

**Done when.** `ipSAE_d0dom` matches `ipsae.py` within ±0.001 on both fixtures.

### `[ ] R004 — Rename ipTM to ipTM_d0chn in all user-facing text`

**What.** The notebook's "ipTM" is `ipTM_d0chn`, a reimplementation. AlphaFold's own
`ipTM_af` is a different number that AFDB does not expose in these endpoints.

**Why.** A reader comparing this notebook's "ipTM" against an AlphaFold output will get a
mismatch and assume the notebook is broken.

**Done when.** Labels, the summary table, and the Section 4 prose say `ipTM_d0chn` and
explain the distinction.

### `[x] R005 — Match calc_d0_array clamping for d0res` — DONE 2026-09-08

**What.** Use `max(26, L)` clamping followed by the unconditional cubic form for the
per-residue `d0`, and keep the scalar form for `d0chn` / `d0dom`.

**Supporting notes.** `ipsae.py` genuinely uses two different `d0` functions. Only affects
residues with exactly 27 valid pairs, but the project tolerance is ±0.001 and the difference
there is 0.039.

### `[x] R006 — Settle the zero-contact return values` — DONE 2026-09-08

**What.** Follow `ipsae.py` (`0.0`) and document why it differs from the sigmoid minima that
Bryant's and Zhu's published code returns. Update the `CLAUDE.md` edge-case table.

---


**R012 verification, 2026-09-08.** All seven values agree with the real `ipsae_v4.py`, run as a
subprocess on both fixtures at cutoffs `10 8`. Worst delta 4.7e-5, and that is the reference's
own print precision: `ipsae.py` writes pDockQ/pDockQ2/LIS with `%8.4f` and the ipSAE/ipTM
columns with `%8.6f`, where agreement is ~4e-7. The ipSAE values were additionally cross-checked
against an independent orchestrator-side computation written before the module existed, and
match to six decimals on both fixtures.

**Evidence that the fixes were needed, not theoretical.**
- R003: heterodimer `n0dom` is 250 (A→B) vs 252 (B→A). Reusing the A→B `d0dom` for B→A, as the
  notebook does, corrupts the B→A value by **0.00126** — past the ±0.001 tolerance.
- R007: heterodimer pDockQ2 is 0.705404 (A→B) vs 0.685271 (B→A), a **0.020** spread.
- Both fixtures happen to have A→B win the max, so neither bug moves the headline number *here*.
  R007 is right by luck on these two inputs; swap the chain order and the notebook is wrong by
  0.020. That is precisely why the fix is not optional.
- R005: a residue with exactly 27 valid pairs scores 0.519068 with `d0_array` versus 0.500000
  with `d0_scalar` — 19x the tolerance.
- R006: verified `pDockQ` and `pDockQ2` both return exactly `0.0` at zero contacts.

**One agent claim rejected on verification.** The agent reported that `formula-reference.md` D1
and the R007 task text have the pDockQ2 directions swapped, citing figures 0.6853 and 0.7054.
Checked: neither number appears in any spec file, and D1 contains no numbers at all — it states
only that the notebook reports the A→B value alone, which the agent's own measurements confirm.
The claim appears to conflate the plan's Q4 PAE-cutoff table (0.7057 vs 0.6981) with pDockQ2
directions. No spec edit was made, and none is needed.

## W1 — Module extraction

> User note: *"Offload every utils and function into a module. It doesn't have to be perfect,
> just focus on decluttering the notebook with the purpose of decluttering it."*

Target: notebook code cells become orchestration and narrative only. Anything that is a
reusable function moves out.

### `[x] R010 — Create the module and the Colab bootstrap` — DONE 2026-09-07
### `[x] R010b — Fix the Colab bootstrap` — DONE 2026-09-08

**Why this exists.** R010 was verified locally but not on Colab, and Colab was the path that
mattered. The user ran it and it failed with a bare `ImportError`. Three defects, all real:

1. **Branch pin.** The clone hardcoded `main`, which does not carry
   `complex_interface_utils.py`. Now pinned to `homodimer-notebook-rework` behind a
   `TODO(merge)` marker; R093 reverts it.
2. **Stale clone.** An existing clone was reused without refreshing, so a checkout predating
   the module kept failing even after the branch was pushed. Now fetches and hard-resets,
   falls back to the existing checkout when offline, and moves a non-git directory aside
   rather than deleting it.
3. **Opaque failure.** The `ImportError` named no cause. It now reports repo root, revision in
   use, branch expected, and the remedy.

**Safety property, verified independently.** `git reset --hard` appears only inside
`ensure_colab_clone()`, which has exactly one call site, guarded by
`IN_COLAB and (REPO_ROOT is None or REPO_ROOT == COLAB_CLONE_DIR.resolve())` and passed only
`COLAB_CLONE_DIR`. A locally discovered repo root can never be reset. A Colab user with their
own checkout elsewhere is also left alone, because `REPO_ROOT` then matches neither arm.

**Lesson recorded.** Verifying only the convenient environment is not verification. Any task
whose deliverable runs in two environments must be exercised in both, or the untested one must
be declared untested. See the amended verification protocol below.

**What.** Create `src/insightfold/complex_interface_utils.py` (D7) with a clear section
layout, and rewrite notebook cell 1 as the D9 clone-and-path bootstrap.

**Supporting notes.**
- Current cell 1 shells out to `pip install -q molviewspec` on every run including local.
  The replacement installs only what is missing and stays inside the 60 s Colab budget from
  `CLAUDE.md`. Allowed dependencies remain numpy, matplotlib, seaborn, requests, molviewspec,
  ipywidgets.
- Bootstrap shape, adapted from the user's `latent-structure-diffusion` notebook. Keep:
  `IN_COLAB` detection via `try: import google.colab`, a `find_repo_root()` that walks up
  looking for `pyproject.toml` + `src/`, and `sys.path.insert(0, root / 'src')`. Drop:
  Drive mounting, the GitHub token / `getpass` / auth-header path (repo is public), torch and
  CUDA setup, `gdown` dataset download, pandas/pyarrow artifact directories, the
  `BACKBONE_DIFFUSION_*` env-flag surface, and `py3Dmol`.
- Add `--depth 1` to the clone; full history drags several MB of committed notebook outputs
  for no benefit.
- Once the rework lands, pin the clone to a tag rather than `main` so a Colab run a year
  later reproduces instead of picking up drift. Track as part of R093.
- The `find_repo_root()` pattern is what lets one cell work unchanged locally and on Colab,
  which is the property worth preserving from the original.

### `[x] R011 — Move AFDB fetch and parsing into the module` — DONE 2026-09-08

**Covers cells 5, 6, 7, 8, 9.** `fetch_metadata`, `download_structure`, `download_pae`,
`download_plddt`, `parse_mmcif_atoms`, `extract_chain_coords`, `parse_pae`, `parse_plddt`.

**Supporting notes.**
- The PAE and pLDDT documents are served **gzipped**; `requests` handles this transparently,
  raw `urllib` does not. Keep using `requests`.
- Cell 9 has dead code: `plddt_ch_sorted` is computed and never used, and both branches of
  the `if len(plddt_chains_info) >= 2` check are identical.
- The notebook's `parse_mmcif_atoms` strips the `_atom_site.` prefix from column names;
  `interface.py`'s keeps it. Pick one convention in the module (see R015).

### `[x] R012 — Move scoring into the module` — DONE 2026-09-08

**Covers cell 18.** `d0_scalar`, `d0_array`, `ptm_func`, `compute_iptm_d0chn`,
`compute_ipsae`, `compute_pdockq`, `compute_pdockq2`, `compute_lis`, plus `THRESHOLDS` and
`traffic_light`. Apply R003, R005, R006 while moving.

**Supporting notes.** Per D4, each function takes an explicit `(block_XY, block_YX, nX, nY)`
or a small `ChainPair` object rather than assuming "A and B".

### `[x] R013 — Move plotting into the module` — DONE 2026-09-08

**Covers cells 13, 15, 16, 20, 22, 27.** Each becomes a `plot_*` function returning `fig`,
so the notebook cell is one call plus a title. Palette handling from R051 lives here.

### `[ ] R014 — Move MolViewSpec builders into the module`

**Covers cell 24.** One builder function per view, plus the shared `show_mol_view` helper.
Interacts with everything in W7.

### `[ ] R015 — Absorb interface.py, defer its deletion`

**What.** `complex_interface_utils.py` absorbs the useful parts of `src/insightfold/interface.py`.
Per D10, `interface.py` itself stays on disk until R095.

**Supporting notes.**
- `interface.py` is imported by nothing in the repo; it is referenced only in `CLAUDE.md`,
  `specs/homodimer_notebook_extraction.md`, and the two agent-suite documents.
- There are currently **three** divergent copies of `extract_cb_coords`: `interface.py`,
  `homodimer_diagnostic.ipynb` cell 8, and `notebooks/analysis_template.ipynb`. They differ
  in whether the `_atom_site.` prefix is stripped from column names. Pick one convention.
- **Name collision to resolve, or it becomes a silent scoring bug.**
  `InterfaceResult.n_contacts` means *interface residue count* (`mask_A.sum() + mask_B.sum()`),
  while the notebook's `n_contacts` means *contact pair count*. pDockQ needs the **pair**
  count (`ipsae.py` L652, `npairs`). Rename to `n_interface_residues` and `n_contact_pairs`
  so neither can be passed where the other is meant.
- `interface.py`'s dataclasses (`ChainCoords`, `InterfaceResult`) are worth keeping; they are
  the natural carrier for the explicit chain-pair API in D4.

### `[ ] R016 — Switch the notebook onto the module and delete the inline code`

**This is the task that delivers the decluttering.** R011-R014 build module functions
*alongside* the notebook's own copies, so that the notebook never breaks mid-refactor. Until
R016 runs, every piece of logic exists twice and the notebook still executes its own copy.

**Plan defect, found by the user 2026-09-08 and corrected here.** R016 previously read "Slim
the notebook import cell. Covers cell 2." That covered the imports only. No task in the plan
owned removing the inline implementations, which is the user's original request ("Offload
every utils and function into a module... focus on decluttering"). Rescoped below.

**Measured at the time of writing (1-indexed cells).** 1058 code lines total: 268 already
duplicated by the module, 361 more once R013/R014 land, leaving ~429 lines of genuine
orchestration and narrative.

| Cell | Lines | Inline logic | Replaced by |
|------|-------|--------------|-------------|
| 3 | 26 | imports, plot style, colour constants | keep the style block and colours; drop the rest, and drop the blanket `warnings.filterwarnings('ignore')` which currently hides the deprecations R075 needs to see |
| 7 | 21 | AFDB metadata fetch | `fetch_afdb_metadata`, `download_*` |
| 9 | 86 | mmCIF parsing, CB/CA extraction | `parse_mmcif_atoms`, `parse_structure` |
| 10 | 37 | PAE and pLDDT parsing, quadrant slicing | `parse_pae`, `parse_plddt`, `ordered_pair` |
| 13 | 21 | CB-CB distance matrix, interface masks | `detect_interface` |
| 14, 16, 17, 21, 23, 28 | 210 | matplotlib figures | `plot_*` from R013 |
| 19 | 103 | all seven scoring functions | `compute_*` from R012 |
| 25 | 143 | MolViewSpec view builders | R014 builders |

**Blast radius — why this is one coherent pass, not a per-cell trickle.** Cell 19 defines the
`compute_*` functions and cell 20 calls them, building the `scores` dict and the `res_*`
result dicts that cells 21, 27, 28 and 29 then consume by key (`res_ipsae['per_res_AB_d0res']`
and similar). The module returns dataclasses, not dicts. So switching the producer forces
updating every consumer in the same task. Splitting that across tasks would leave the notebook
broken in between, violating the Invariant.

**Equivalence criterion.** Capture every printed number and every score from a full Run All
*before* the switch, on both fixtures. Repeat *after*. Any numerical difference is a bug, not
drift, because R011 and R012 both verified bit-identical behaviour against the inline code at
the time they landed. Figures are compared visually; numbers must match exactly.

**Sequencing.** R016 runs last in M2, after R013, R014 and R015, so the notebook is rewritten
once against a complete module rather than twice against a partial one.

---


**R013 verification, 2026-09-08.** Deliberately a behaviour-preserving move, not a redesign:
R050/R051/R052 own the visual changes, and mixing a move with a redesign would leave no way to
tell an intended change from a regression. Numerical fixes could be landed early in R012
because `ipsae.py` is an objective oracle; visual work has none.

Result: **byte-identical PNGs**. The full PAE heatmap and the contact map render to the same
SHA-256 as the notebook's original inline code on both fixtures, zero differing pixels, and the
score-mask panel and pLDDT figure match too. Independently re-verified by the orchestrator on
the heterodimer. All six figures render on both fixtures with correct per-chain axes (181 and
101 residues), so nothing assumed a square matrix.

**Seams built for later milestones.**
- `PAE_CMAP` + `resolve_pae_cmap()` read the default at call time, so R051 flips one constant
  rather than rewriting six functions. Default stays `RdBu_r` today.
- `apply_plot_style()` is explicit and never runs at import, so importing the module cannot
  mutate global matplotlib state. Verified: `figure.dpi` is unchanged by import.
- Figures are built from `matplotlib.figure.Figure`, not `plt.subplots`, so a returned figure
  is not double-displayed by the inline backend and nothing leaks into pyplot's registry.
  Verified: 0 open figures after a module call.

**Two figures deliberately NOT byte-identical, both correct.** The per-residue profiles now
carry the R003/R007 directional fixes, so their data genuinely differs. The agreement matrix
normalises by `THRESHOLDS[...].green` instead of the notebook's stale inline dict and labels
the score `ipTM_d0chn` per R004 — preserving the old numbers there would have meant preserving
thresholds R002/R009 corrected.

## W2 — Heterodimer support

> User note: *"Make sure to test notebook on heterodimers too, as it's only really been
> tested on homodimers at the time of writing it."* Plus D3: works automatically, everywhere.

### `[ ] R020 — Fix the metadata chain-selection bug`

**What.** Cell 6 does `meta = entries[0]`. The AFDB prediction endpoint returns **one entry
per chain**, and for the fixtures I checked **chain B comes first**.

**Supporting notes — CORRECTED 2026-09-08. The entry order is non-deterministic.**

The original note here recorded "chain B comes first", from two samples that both happened to
return B first. That was an over-generalisation. Re-tested with repeated calls:

```
AF-0000000065889468  try1 ['A','B']   try2 ['A','B']   try3 ['B','A']
AF-0000000211034637  try1 ['B','A']   try2 ['B','A']   try3 ['A','B']
```

**The same accession returns a different order on different calls.** So `meta = entries[0]` is
not merely wrong, it is *non-deterministically* wrong: for a heterodimer, the sequence, gene,
protein name, monomer length and UniProt ID shown in the cell-10 report can flip between chains
on consecutive runs of the same notebook with the same input. For a homodimer it stays
invisible, because both entries describe the same protein.

Never code against either order. Select by `chainId` explicitly:
`{entry['chainId']: entry for entry in entries}`, and report per chain. R011 already preserves
all entries in `AFDBPrediction` with `entry_for_chain()`, so this task is the selection logic,
not a rewrite.

`cifUrl` / `bcifUrl` / `paeDocUrl` / `plddtDocUrl` are identical across entries on both
fixtures, and `AFDBPrediction.document_url()` now enforces that agreement rather than assuming
it, so a future per-chain document scheme fails loudly instead of mis-slicing silently.

**This also strengthens the case for R008's directional diff:** a value that changes between
runs of an unchanged notebook is the hardest class of bug to notice by eye.

### `[ ] R021 — Chain identity and labelling`

**What.** Resolve structure chain IDs against the PAE `chains` field, and use real protein
names in labels instead of "Chain A" / "Chain B".

**Supporting notes.** Both the PAE and pLDDT JSON carry a `chains` array with a `name` field
that the notebook currently ignores:

```json
{"name": "Interferon-stimulated gene 20 kDa protein", "label_asym_id": "A",
 "sequenceStart": 1, "sequenceEnd": 181}
```

The `CLAUDE.md` edge-case table already warns that asymmetric chain labels can misalign PAE
quadrant slicing. Add an explicit assertion that structure-derived chain lengths equal the
PAE-derived `nA` / `nB`, and fail loudly rather than silently mis-slicing.

### `[ ] R022 — Remove nA == nB assumptions`

**What.** Sweep every cell for places that only work when the chains are the same length.

**Supporting notes — known sites.**
- Cell 16: `mask_iptm = np.ones((nA, nB))` uses PAE-derived lengths while `mask_pdockq2 =
  contact_mask` uses structure-derived lengths. Shape mismatch if they ever differ.
- Cell 13: `ax2.set_xlim(0, max(nA, nB))` puts both chains on one axis, so the shorter chain's
  bar looks truncated rather than shorter. Needs per-chain axes or an explicit shared scale.
- Cell 24 views 3 and 4: index per-residue score arrays (length = PAE `nA`) with structure
  residue lists. Fine only while the two agree.
- Cell 20: per-residue profile panels use `sharex=False` already, so unequal lengths are OK,
  but the axis labels need the chain name from R021.
- Cell 27 and 28: the diagnostic prose says "the two chains" and "identical copies" in places
  that are false for a heterodimer.

### `[ ] R023 — Detect and display the assembly type`

**What.** Read `assemblyType`, `oligomericState`, and `complexComposition` from the metadata
and state plainly at the top of the report whether this is a homodimer or a heterodimer.

**Supporting notes.** Available fields, verified live:

```
assemblyType: "Homo" | "Hetero"
oligomericState: "dimer"
complexComposition: [{"identifierType": "uniprotAccession",
                      "identifier": "P0A6Q3", "stoichiometry": 2}]
isComplex: true
```

None of these are documented in `CLAUDE.md`'s field table — add them in R093. Also add a
clear early failure when `oligomericState` is not `dimer`, replacing the current silent
zero-length chain B path noted in the `CLAUDE.md` edge cases.

### `[ ] R024 — Register heterodimer fixtures`

**What.** Add heterodimer fixtures to `fixture-manifest.md` as FX-006 onward.

**Supporting notes.** There is a working AFDB search endpoint the spec pack does not mention:

```
GET /api/search?q=assemblyType:Hetero AND oligomericState:dimer&type=main&rows=6
    -> numFound: 80240
```

Verified candidates:

| Candidate | Chains | Global pLDDT | Note |
|-----------|--------|--------------|------|
| `AF-0000000211034637` | Q96AZ6 (181 aa) + P63166 (101 aa) | 87.82 | primary heterodimer fixture; unequal lengths, cross-species (human + mouse), PAE 282×282, parses clean |
| `AF-0000000211619209` | Q96M98 + Q5JSS6 | 73.57 | mid-confidence human heterodimer |
| `AF-0000000211157965` | Q9P2Y5 + Q14457 | 68.91 | lower-confidence, useful for amber/red paths |

`AF-0000000211034637` is the one I will use for per-task verification under D5.

---

## W3 — Section 2, Interface Detection

> User note: *"we should maybe add a 3D structure of residues involved with interaction.
> Currently it nicely shows the 2D view, via a contact matrix and an annotation track."*

### `[ ] R030 — Add a 3D interface view to Section 2`

**What.** Insert a MolViewSpec view directly after the contact map showing the interface
residues in 3D, so the 2D contact matrix and the 3D reality sit side by side.

**Supporting notes.** This is essentially current cell 24 View 1, moved forward. Decide
whether Section 6's View 1 then becomes redundant or stays as the "here is the complex"
orientation shot. Depends on R014 and R071 (both-chains side-chain colouring).

---

## W4 — Prose pass

### `[ ] R040 — Em dash sweep`

> User note: *"Read over all the markdown text, making sure the em dashes disappear, unless
> you want them there explicitly."*

**Scope.** Per A2: markdown cells, plot titles and axis labels, printed output. Not code
comments. En dashes in numeric ranges stay.

**Supporting notes.** Em dashes appear in all eight markdown cells and in roughly a dozen
plot titles and print statements. Replace with commas, colons, parentheses or full stops
according to what each sentence actually needs, rather than a blanket substitution.

**Done last** in each section, so it catches text written by the other tasks. Better run
once at the end of W5-W8 than repeatedly.

---

## W5 — Section 3, PAE Matrix Decomposition

### `[ ] R050 — Fix the full PAE heatmap size`

> User note: *"reduce the size of the PAE matrix (currently you have an in-cell scroll, which
> makes it really hard to view it all at the same time)."*

**Supporting notes.** Cell 15 uses `figsize=(10, 9)` at `figure.dpi = 150` (set in cell 2),
producing a 1500×1350 px image that overflows the output area and triggers the scrollbar.
Fix by reducing the figure size and setting `aspect='equal'` so the matrix stays square.
`aspect='auto'` is currently used, which also distorts the matrix for a heterodimer.

### `[ ] R051 — Green default palette with a switchable alternative`

> User note: *"Change the colour to the standard green palette... Personally, green is hard
> to see against white. Also the issue of colour blindness. But the green palette should be
> the default regardless."*

**What.** A module-level `PAE_CMAP` setting. Default: the AFDB-style sequential green where
dark green = low PAE = confident and pale = high PAE = uncertain. Alternatives offered:
a colourblind-safe sequential option, and the current `RdBu_r` for continuity.

**Supporting notes.** Two things to get right. First, low PAE must map to the *dark,
saturated* end, so confident regions stay legible on a white background — that addresses the
"hard to see against white" point directly. Second, cell 14's markdown currently says
"Low PAE (blue) ... High PAE (red)", which becomes wrong the moment the default changes; it
has to be rewritten to describe the palette generically or to match the new default.

### `[ ] R052 — Rebuild the score-mask panel`

> User note, three parts: explain the quadrant provenance, apply the palette, fix the aspect
> ratio and colour bar.

**Supporting notes.** Cell 16 problems, in order:
1. **Provenance.** The panel shows the A→B quadrant but never says so beyond the axis labels.
   The user explicitly wants it stated that this is the top-right quadrant of the full matrix
   from cell 15, extracted as `pae_matrix[:nA, nA:nA+nB]`. Consider a small inset locator
   showing which quadrant of the full matrix is being displayed.
2. **Aspect ratio.** `plt.colorbar(im, ax=axes[-1], ...)` attaches a colour bar to the last
   axis only, which steals width from that panel alone, so the fourth matrix renders narrower
   than the other three. Fix with one shared figure-level colour bar
   (`fig.colorbar(im, ax=axes, ...)`) or a dedicated colour bar axis. The user's fallback of
   stacking the four vertically is acceptable if a 2×2 cannot be made clean.
3. **Palette.** Inherit `PAE_CMAP` from R051. The grey "not used" underlay needs to stay
   distinguishable from the palette's pale end, which it will not be against a green ramp.
4. Also `aspect='auto'` here distorts the block; use `aspect='equal'`.
5. The subtitle text should name the actual cut-offs in force rather than hard-coding them.

---

## W6 — Section 4, Score Computation

### `[ ] R060 — Apply the formula audit`

**What.** Land R003, R004, R005, R006 in the notebook and the module, then verify against
`ipsae.py` on both fixtures.

**Supporting notes.** The user's ask was *"a one-over sanity check of each formula, to make
sure it makes sense."* The audit is complete and is recorded in R001. Net result: the
notebook is correct on pDockQ, pDockQ2 and LIS (where `CLAUDE.md` is wrong), and has one real
bug of its own (d0dom asymmetry) plus one edge-case divergence (`L == 27`).

### `[ ] R061 — Restructure Section 4 into per-score subsections`

> User note: *"We should [have] a section for each score, explain the formula, explain the
> thinking behind it based on the publication, and then calculate it, before moving on to the
> next section."*

**What.** Replace the single "define everything, then run everything" pair of cells with
subsections 4.1 to 4.6, each being markdown (what it measures → the formula → why the
authors built it that way → how to read the number) followed by one compute-and-report cell.

**Proposed order.** ipTM_d0chn → ipSAE (all three d0 variants together, since the variants
only differ in `L`) → pDockQ → pDockQ2 → LIS → combined score table.

**Depends on.** Q1. The "thinking behind it based on the publication" cannot be written
honestly without the papers.

**Supporting notes.** The user considers the existing PAE-matrix walkthrough and mask panel
sufficient for the mechanical breakdown, so these subsections should carry the *reasoning*,
not re-derive the arithmetic.

### `[ ] R062 — Rework the per-residue score profile plot`

> User note: *"The plot in the per-residue score profiles is odd... Actually, no, it's only
> specific for ipTM and ipSAE d0res, d0chn, and d0dom. I like the lightly shaded interface
> residues."*

**What.** Keep the plot restricted to the four values that genuinely have a per-residue
decomposition (ipTM_d0chn, ipSAE_d0res, ipSAE_d0chn, ipSAE_d0dom), keep the interface
shading and the pLDDT underlay, and add markdown explaining *why* only these four appear:
pDockQ, pDockQ2 and LIS are pooled statistics over contact pairs, so they have no per-residue
value to plot. Also mark the argmax residue on each panel, since that single residue is what
the reported score actually equals.

**Supporting notes.** The current cell 20 already plots exactly the right four series — the
user's confusion was about whether that was intentional. The fix is mostly explanatory text
plus the argmax marker, not a rewrite. Chain labels come from R021.

---

## W7 — Section 6, 3D Structure Visualisation

### `[ ] R070 — Write supporting text for every view`

> User note: *"What each plot is showing, and what we are colouring by."*

**What.** Each view gets a markdown block above it: what is drawn, what the colours encode,
what to look for, and what a problem looks like.

### `[ ] R071 — Show side chains on both chains, coloured by chain`

> User note: *"because it's homodimer, we are only visualising the side-chains of only one
> chain. But I think we should visualise on both chains regardless... So currently it's tomato
> (chain A) and teal (chain B), but the side chains on the chain A would be cornflower blue
> and orange for chain B."*

**What.** Every interface view renders side chains for both chains, with a distinct
side-chain colour per chain that reads clearly against that chain's cartoon colour.

**Supporting notes.** Two corrections to the note, both cosmetic. The current cartoon colours
are teal `#009688` for A and coral `#FF7043` for B (the note has them the other way round).
And View 1 already draws chain A's interface side chains only — Views 3 and 4 also draw chain
A only in View 4's case, while View 3 does cover both chains. Proposed scheme, to check
against a colourblind simulation before landing:

| Element | Chain A | Chain B |
|---------|---------|---------|
| cartoon | teal `#009688` | coral `#FF7043` |
| interface side chains | cornflower `#6495ED` | orange `#F5A623` |

### `[ ] R072 — Explain View 3 (ipSAE d0res on the interface)`

> User note: *"why are we using the d0res score? Not that it's wrong, but nothing is explained
> about it. And also, what is low, mid, and high? Need to explain the value ranges."*

**What.** Two additions. First, state why `d0res` is the variant shown: it is the AFDB primary
classifier and the only variant whose `d0` is derived per residue, which makes it the one that
varies meaningfully *along the chain* and therefore the only one worth painting onto a
structure. Second, replace "red=low, yellow=mid, green=high" with the actual numbers, tied to
the R002 thresholds, and add a colour bar or an explicit legend so the mapping is readable.

### `[ ] R073 — Rebuild View 4 (disagreement)`

> User note: *"it currently mentions 'Disagreement (green=PAE+contact, blue=PAE confident/no
> contact, red=contact/low PAE)', but what does that mean? I see PAE + contact. Is it good PAE
> or bad PAE or what?"*

**Supporting notes.** The current categories, decoded from cell 24, use a hard-coded
`ipsae_thresh = 0.5` that is unrelated to the Section 7 thresholds:

| Colour | Condition | What it actually means |
|--------|-----------|------------------------|
| green `#4CAF50` | `per_res_ipsae > 0.5` **and** is an interface residue | both signals agree that this residue is a confident contact |
| blue `#2196F3` | `per_res_ipsae > 0.5` **and** *not* an interface residue | PAE is confident about the partner chain, but no CB atom is within 8 Å |
| red `#F44336` | `per_res_ipsae <= 0.5` **and** is an interface residue | there is a physical contact, but PAE is not confident about it |

Three fixes: name the categories in terms the reader can act on rather than as a colour key,
source the threshold from R002 instead of a magic 0.5, and colour both chains (only chain A is
coloured today, and the `colour_map_view4` tuples carry a fourth label element that is
discarded with `_`, so the legend was clearly intended and never built).

### `[ ] R074 — Add pDockQ2 views`

> User note: *"why did we only plot ipSAE? ... What about pDockQ2? You mentioned how strict it
> is. So let's give it space to shine as well. Maybe have two plots for it like the third and
> fourth Mol* plots for ipSAE."*

**What.** Two new views mirroring 3 and 4.
- **pDockQ2 contact quality.** Colour each interface residue by the mean `ptm(PAE, d0=10)`
  over that residue's contact pairs, which is the per-residue quantity pDockQ2 pools. This is
  a genuine per-residue decomposition of pDockQ2 and does not exist anywhere in the notebook
  yet.
- **pDockQ2 vs ipSAE disagreement.** Residues where the contact-PAE view and the ipSAE view
  reach different conclusions, which is exactly the "why do the scores disagree" question the
  notebook exists to answer.

**Supporting notes.** Requires the module to expose the per-residue `mean_ptm` intermediate,
which `compute_pdockq2` currently pools away. Fold that into R012.

### `[ ] R075 — Fix the MolViewSpec implementation issues`

**Supporting notes.**
- `cm.get_cmap('RdYlGn')` in View 3 is deprecated and slated for removal. Use
  `plt.get_cmap`. The blanket `warnings.filterwarnings('ignore')` in cell 2 is currently
  hiding this (see R016).
- Every view builds **one component per residue**. View 2 builds one for all 282 residues of
  the heterodimer fixture. Collapse contiguous residue runs of the same colour into a single
  `ComponentExpression` with `beg_label_seq_id` / `end_label_seq_id` — the API already
  supports ranges, the code just always passes `beg == end`.
- Views 3 and 4 hard-code `label_asym_id='A'` / `'B'`; they should use the resolved chain IDs
  from R021.
- `struct_fmt` is chosen with `'bcif' if 'bcifUrl' in meta else 'mmcif'`, but `struct_url` is
  chosen with `meta.get('bcifUrl', meta.get('cifUrl', cif_url))`. If `bcifUrl` is present but
  empty (the `USE_LOCAL_FILE` path builds exactly that), the URL is empty while the format
  says bcif.
- The `USE_LOCAL_FILE` path sets `bcifUrl: ''`, so Section 6 cannot work at all offline.
  Either disable the section with a clear message or resolve a usable URL.

---

## W8 — Section 7, Diagnostic Summary

> User note: *"the traffic light system is great. Just need to make sure we are using the
> right values or scores when describing each metric."*

### `[ ] R080 — Wire the summary to the canonical thresholds`

**What.** Cells 26 and 28 both read `THRESHOLDS` from the module (R002). Delete cell 28's
`green/2` amber recomputation, which is the direct cause of the table and the prose being
able to disagree.

**Also required by R002.** Surface threshold **provenance** in the summary table. A user
reading a coloured badge currently has no way to tell a published ROC-optimal cutoff from a
round number someone picked. Each row should declare PUBLISHED / DERIVED / AFDB-EMPIRICAL /
HEURISTIC. Source: `threshold-reference.md`.

Two specific statements the notebook must make (user, 2026-09-07):
1. **Where ipSAE 0.6 comes from:** a data-driven analysis of the ~31 million predicted dimers
   in AlphaFold DB. Unpublished but empirically grounded, and a far larger calibration set
   than any of the four papers. It is AFDB's number, not Dunbrack's — the paper states no
   ipSAE cutoff at all.
2. **Where the traffic light is and is not discriminating:** for an AFDB accession the model
   has already passed the 0.6 filter, so green is largely re-measuring the criterion that
   caused the entry to exist and red will essentially never appear. For a model uploaded via
   `USE_LOCAL_FILE` no such filtering has occurred, so the full range including red is
   reachable and meaningful. State both; the limitation is specific to the AFDB path.

**All five remaining judgement calls are amber thresholds.** Every green now has either a
publication or the AFDB calibration behind it. None of the four papers defines a middle band
at all, so the amber tier is this notebook's own construct and should be presented that way.

### `[ ] R081 — Audit the metric descriptions`

**Supporting notes.** Cell 26's `descriptions` dict is mostly sound but has issues to check
against R001: "ipTM: All PAE cells used, no cutoff" is correct but should say `ipTM_d0chn`
(R004); "pDockQ: Contact count x interface pLDDT" understates it (`log10` of the contact
count, and the pLDDT is a mean over the union of contacting residues in both chains);
"ipSAE_d0chn: Penalises small interfaces more" needs checking — a larger `L` gives a larger
`d0`, which is *more* lenient, so the direction of that claim may be backwards.

### `[ ] R082 — Review the agreement matrix`

**Supporting notes.** Cell 27 normalises each score by its own green threshold and caps at
1.5, then plots pairwise absolute differences. That makes the matrix depend entirely on the
R002 thresholds being right, and the cap means two scores that are both far above threshold
register as perfect agreement. Decide whether this is the comparison worth showing, and label
the axes with what the normalisation actually is.

---

## W9 — Validation and documentation

### `[ ] R090 — Homodimer end-to-end run`
Fixture FX-001 `AF-0000000065889468`. Full run, clean, no warnings.

### `[ ] R091 — Heterodimer end-to-end run`
Fixture `AF-0000000211034637` (R024). Full run, every section correct for unequal chains.

### `[ ] R092 — Numerical agreement with ipsae.py`
All seven values within ±0.001 of `ipsae.py` v4 on both fixtures. Record in
`specs/homodimer_diagnostic/validation-report.md`, which T038 in the original `tasks.md`
already calls for and which does not yet exist.

### `[ ] R093 — Sync CLAUDE.md and the spec pack`
- Correct the pDockQ2 sigmoid constants, the pDockQ contact count, the LIS combination rule,
  and the d0dom definition (R001).
- Replace `THRESHOLDS` with the R002 table.
- Document the AFDB fields the notebook now uses: `assemblyType`, `oligomericState`,
  `complexComposition`, `isComplex`, `chainId`, and the `name` field inside PAE `chains`.
- Document that the prediction endpoint returns one entry per chain, chain B first (R020).
- Document the `/api/search` endpoint used for fixture sourcing (R024).
- Update the module layout section for `complex_interface_utils.py`, and record D8 (name
  modules by domain; promote on the second consumer) so future sessions follow it.
- Document the D9 bootstrap and pin the clone to a release tag rather than `main`.
- Update the edge-case table for R006 and R023.
- **Reconcile the dependency policy.** `CLAUDE.md` states biopython, torch and
  torch-geometric are prohibited, but `pyproject.toml` lists `biopython>=1.83,<2` as a hard
  dependency (alongside gemmi, scipy, networkx, plotly, pandas, jupyter). Both cannot be
  right. This is exactly why D9 avoids `pip install` — but the contradiction should be
  settled rather than routed around, since it will mislead the next reader.

### `[ ] R095 — Delete the superseded code (last task, per D10)`

**What.** Only once R090-R093 pass. Remove `src/insightfold/interface.py` (absorbed by R015)
and the duplicated inline `extract_cb_coords` in `notebooks/analysis_template.ipynb`, and
fix the four documents that reference `interface.py` by path.

**Why deferred.** Deleting early makes any later failure ambiguous: a real bug and a missing
import look identical from the traceback. Nothing imports `interface.py` today, so leaving it
in place costs nothing.

### `[ ] R094 — Close out the original tasks.md`
Mark T004 (threshold bands), T033-T038 (validation), and T039-T041 (review) against the work
done here, or point them at this file.

---

## Milestones and execution order

45 tasks in 8 milestones. Execution model, agreed with the user 2026-09-07: each task goes to
a **fresh agent**, strictly **sequential**, and each result is verified against source and by
execution before the next is dispatched. At every milestone boundary, work **stops** for user
review.

| # | Milestone | Tasks | n | Status |
|---|-----------|-------|---|--------|
| M1 | Ground truth | R001, R002, R002b, R009 | 4 | R002b added 2026-09-08 |
| M2 | Module extracted; **R016 is where the notebook actually shrinks** | R010, R010b, R011-R016 | 8 | R010, R010b, R011, R012, R013 done |
| M3 | Heterodimer support | R020-R024 | 5 | |
| M4 | Scoring correctness + Section 4 | R003-R008, R060-R062 | 9 | R003, R005, R006, R007 landed early in R012 |
| M5 | PAE visuals | R050-R052 | 3 | |
| M6 | 3D views | R070-R075, R030 | 7 | |
| M7 | Summary + prose | R080-R082, R040 | 4 | |
| M8 | Validation, docs, cleanup | R090-R095 | 6 | |

Rationale for the ordering. M1 first because every number downstream depends on it, and it
touches no code so it is cheap to redo. M2 before any behavioural change, so later edits land
in one place rather than being made twice. M3 before M4 because heterodimer support is what
makes the directional scoring bugs (R003, R007) observable at all — both are invisible on a
homodimer. M4 before the visual work because the visuals plot the scores. M7's em dash sweep
runs late so it covers text written by earlier milestones. M8's deletions run last, per D10.

**Verification protocol.** Read the diff rather than the agent's report; execute the affected
cells against a live AFDB fetch; check the task's "done when"; for scoring tasks compare
numerically against `references/ipsae_v4.py` on both fixtures. A failed task is re-issued to a
fresh agent with the specific failure, not patched by the orchestrator.

**Environment coverage.** A deliverable that runs both locally and on Colab must be verified
in both, or the untested environment must be named as untested in the report. R010 passed
local verification and shipped a Colab-breaking bug; the gap was the orchestrator's, not the
agent's. Colab cannot be driven from here, so the standing substitute is: exercise the
underlying functions against a simulated remote, and state plainly which paths were executed
versus reasoned about.

**Invariant — the notebook stays runnable after every single task.** M2 moves code out of the
notebook into the module over seven tasks; at no point may the notebook be left in a state
that only works once a later task lands. Each task either leaves the inline code in place
alongside the new module function, or switches the call site over within the same task.
Verified by execution, not by inspection. This is what makes it safe for the user to open and
run the notebook at any time, not only at milestone boundaries.

**Note on diffing the notebook.** `notebooks/homodimer_diagnostic.ipynb` had uncommitted
changes before this rework began, so `git diff HEAD` shows deltas no agent made: a leading
blank line in cell 4, and stored outputs on cells 5, 6, 7 and 24. Compare against the previous
working state, not against HEAD, when attributing a change to an agent.
