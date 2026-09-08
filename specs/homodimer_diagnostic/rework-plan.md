# Homodimer Diagnostic Notebook — Rework Plan

Source of truth for the April-2026 rework of `notebooks/homodimer_diagnostic.ipynb`.
Every task below traces to a note in the user's "Planned notebook fixes" list or to a
defect found while auditing the notebook against `ipsae.py` v4 and the live AFDB API.

Status legend: `[ ]` not started · `[~]` in progress · `[x]` done · `[!]` blocked

Related documents:
- `specs/homodimer_diagnostic/tasks.md` — the original build plan (Phases 1-7). Still valid
  history; this file supersedes it for all work from 2026-09-07 onward.
- `specs/homodimer_diagnostic/fixture-manifest.md` — fixture registry, extended by R024.
- `CLAUDE.md` — project skill file. **Rewritten by R093 on 2026-09-08.** It no longer
  carries its own copies of the formulas or the thresholds, because those copies are what
  drifted; it points at `formula-reference.md` and `threshold-reference.md` instead.

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

**Q4 — PAE cutoff: ANSWERED BY EVIDENCE 2026-09-08, not by decision.**

The AFDB search endpoint publishes the parameter directly:
`complexPredictionAccuracy_ipsae_pae_cutoff = 10.0`. **AFDB's production pipeline used 10,
which is what this notebook uses.** The earlier concern that AF2-family models might warrant
Dunbrack's 15 does not apply to how AFDB actually computed the released scores. No caveat is
needed and the 0.6 threshold transfers directly to our numbers.

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

### `[x] R002b — Use the AFDB bands report to fill prose gaps only` — DONE 2026-09-08

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

### `[x] R008 — Surface directional asymmetry as a diagnostic` — DONE 2026-09-08

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

### `[x] R004 — Rename ipTM to ipTM_d0chn in all user-facing text` — DONE 2026-09-08

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

### `[x] R014 — Move MolViewSpec builders into the module` — DONE 2026-09-08

**Covers cell 24.** One builder function per view, plus the shared `show_mol_view` helper.
Interacts with everything in W7.

### `[x] R015 — Absorb interface.py, defer its deletion` — DONE 2026-09-08

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

### `[x] R016 — Switch the notebook onto the module and delete the inline code` — DONE 2026-09-08

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


**R014 verification, 2026-09-08.** Behaviour-preserving move, same discipline as R013. A
verbatim replica of notebook cell 24 was built alongside the module builders on identical
inputs: residue→colour maps and whole-structure components matched on all eight views across
both fixtures.

**Component-count reduction from collapsing contiguous same-colour runs.** The notebook emits
one `ComponentExpression` per residue even though the API supports ranges.

| Fixture | View | Before | After |
|---------|------|--------|-------|
| homodimer | pLDDT | 344 | **12** |
| homodimer | disagreement | 154 | 31 |
| homodimer | chain overview | 43 | 13 |
| heterodimer | pLDDT | **282** | **32** |
| heterodimer | disagreement | 26 | 13 |

Residue→colour mapping verified identical in every case, so this is payload only, not a render
change. View 3 barely collapses (86→83) because adjacent interface residues rarely share an
exact colour off a continuous ramp — which is itself evidence the collapse merges nothing it
should not.

**Implementation defects fixed (R075 subset).** Deprecated `cm.get_cmap` replaced; hard-coded
`label_asym_id='A'`/`'B'` in views 3 and 4 replaced with resolved chain ids; and the
url/format mismatch resolved — `resolve_structure_source` makes one decision for both, treats
a present-but-empty field as absent, and raises rather than handing Mol* an empty URL labelled
`bcif` (exactly what the `USE_LOCAL_FILE` path built).

**Seam for R074.** `build_interface_value_view` is deliberately score-agnostic: it takes any
per-residue value array plus a colormap. R074's two pDockQ2 views are therefore calls to
existing code, not new builders.

**Finding for M6 — two pLDDT ladders disagree.** `MVS_PLDDT_BANDS` (3D) uses half-open
intervals with inclusive-below edges, so a residue at exactly 100.00 matches no band and is
not drawn, and 90.0 lands in the top band; the 2D `plddt_band_colour` puts 90.0 in the second
band. Preserved rather than silently reconciled, because changing either shifts colours in a
figure the user has not reviewed. **R070/R072 must reconcile them when View 2 gains its
legend** — a legend that disagrees with the colours it explains would be worse than the
current inconsistency.


**R015 outcome, 2026-09-08. Framed as a gap analysis, and the honest answer was "almost
nothing to port".** Every behaviour in `interface.py` was already covered by
`complex_interface_utils.py`. One documentation gap was found and ported: the CB-CB versus
AFDB-production CA-CA warning, which had no counterpart in the new module and would have been
lost at R095. Total diff 22 lines, purely additive. `interface.py` itself untouched per D10.

**Behavioural equivalence proven on both fixtures.** Same records, same chains, same residue
counts, bitwise-equal coordinates, pLDDT, distance matrices, interface masks and contact
pairs. Homodimer 172+172, 116 contact pairs; heterodimer 181+101, 22.

**The rename mattered, with evidence.** On the homodimer `n_interface_residues` is 86 while
`n_contact_pairs` is 116. Two different numbers that the old bare name `n_contacts` conflated,
and pDockQ needs the 116. On the heterodimer they coincide at 22 by accident of geometry,
which is why the homodimer is the fixture that proves the point.

**Latent bug found in `interface.py`, reproduced independently.** Its parser treats a blank
line inside the `_atom_site` loop as a terminator and silently drops every atom after it. CIF
whitespace is not significant there; the loop ends at the next data name, `loop_`, or `#`.
Neither fixture triggers it, but any file with an internal blank line would have been parsed
with a truncated atom list and no error. The new module handles it correctly. This is a good
argument for R095 actually happening rather than the file lingering.


**R016 outcome, 2026-09-08. The decluttering is delivered.**

**1201 -> 556 code lines (-645, -53.7%).** Excluding the untouched bootstrap cell, -62.3%.
Cell count 29 -> 28: the scoring-functions cell was deleted outright rather than emptied,
since it contained nothing but definitions the module now owns. Biggest reductions: MolViewSpec
156->30, scoring functions 120->0, mmCIF parsing 97->7, pLDDT figure 64->10, contact map 40->1.

**Equivalence held.** Both fixtures executed end to end via `nbconvert --execute` before and
after. Every score and every printed statistic is unchanged: chain lengths, PAE shapes, quadrant
shapes, pLDDT means, contact pairs, interface counts and ranges, mask cell counts, pDockQ/
pDockQ2 intermediates, LIS valid counts, d0 values, argmax residues, and the full diagnostic
prose. Four of six figures are byte-identical by SHA-256; the two that differ are the
intended threshold and directional-fix changes. Independently re-verified by the orchestrator
on the heterodimer: 0 errors, 0 stderr, 15 figures, all seven scores matching.

**Removing the blanket `warnings.filterwarnings('ignore')` surfaced nothing.** Zero warnings on
either fixture. The agent verified the harness was not swallowing them (a probe notebook shows
warnings reaching stderr normally) and established that the expected molviewspec
`PydanticDeprecatedSince20` only appears under `simplefilter('always')`, since Python's default
filter hides third-party DeprecationWarnings anyway. No narrow suppression was needed. The
blanket filter was hiding nothing this notebook produces.

**`%matplotlib inline` was required in the imports cell.** The module's `plot_*` functions
return bare `Figure` objects built without pyplot, so `plt.show()` cannot render them and
`display(fig)` emits only `text/plain` unless the inline figure formatters are registered. The
magic registers them; the notebook no longer imports pyplot at all.

**Found and deliberately not fixed here.**
- **R020 confirmed live.** On the heterodimer, cell 11 reports chain B's UniProt ID, protein
  name and length as if they were the complex's. Left for M3 with a `TODO(R020)` at the call
  site.
- **Endpoint order returned B-first on both fixtures this run**, consistent with the recorded
  non-determinism. The notebook no longer indexes `entries[0]` anywhere; the only remaining
  positional index is `AFDBPrediction.primary_entry` in the module, already marked.
- **For R052:** `plot_pae_score_masks` hard-codes the panel title `'ipTM  (all inter-chain, no
  cutoff)'`, so that one figure still says `ipTM` while every other label now says
  `ipTM_d0chn` (R004). `complex_interface_utils.py:3925`.
- **For R061:** the summary-table description "ipSAE_d0chn: Conservative variant. Penalises
  small interfaces more" is backwards. `d0chn >= d0dom >= d0res` is a theorem, so d0chn is the
  most permissive variant. Carried over verbatim rather than silently reworded.


**R060 verification, 2026-09-08. No agent dispatched: the work was already done and the
acceptance criterion already met, so spawning one would have manufactured activity.**

R003, R004, R005, R006 and R007 all landed in R012, verified against `ipsae.py` run as a
subprocess (worst delta 4.7e-5, itself the reference's print precision). R024 then added an
independent confirmation from AFDB's own production values.

**A methodology finding for R092: AFDB's stored precision varies by entry, and it truncates.**
The heterodimer `AF-0000000211034637` stores per-direction scores at full precision (agreement
to ~1e-6). The homodimer `AF-0000000065889468` stores **every** field at 2dp, including the
per-direction ones. Comparing against the 2dp fields naively produced two apparent mismatches;
both dissolved once truncation was accounted for:

| Value | Ours | Truncated | Rounded | AFDB |
|-------|------|-----------|---------|------|
| ipSAE_d0res | 0.914309 | 0.91 | 0.91 | 0.91 |
| ipTM_d0chn | 0.952902 | 0.95 | 0.95 | 0.95 |
| pDockQ2 | 0.926932 | **0.92** | 0.93 | **0.92** |
| LIS_AB | 0.756447 | **0.75** | 0.76 | **0.75** |
| pDockQ | 0.691274 | 0.69 | 0.69 | 0.69 |

pDockQ2 and LIS distinguish the two conventions, and both say truncate.

**CORRECTED 2026-09-08 by R092, which found this rule over-generalised.** There are **two
regimes**, and the rounding behaviour differs between them:

| Regime | Entries | Per-direction fields | Rolled-up 2dp fields |
|--------|---------|----------------------|----------------------|
| TWO-DP | FX-001, FX-002 (both homodimers) | stored at 2dp, **truncated** | truncated |
| FULL | FX-006 to FX-010 | full precision | **rounded** |

Verified independently: FX-002 `pDockQ` 0.687671 -> AFDB 0.68 is decisive for truncation
(rounding gives 0.69, and pDockQ has no per-direction field to fall back on). FX-006 `ipSAE`
0.705718 -> AFDB 0.71 is decisive for rounding (truncation gives 0.70). Applying the original
blanket rule uniformly produces **11 false failures**.

**Correct method:** use the full-precision per-direction fields wherever they exist; fall back
to a 2dp field only when the entry has nothing better, and then truncate.

## W2 — Heterodimer support

> User note: *"Make sure to test notebook on heterodimers too, as it's only really been
> tested on homodimers at the time of writing it."* Plus D3: works automatically, everywhere.

### `[x] R020 — Fix the metadata chain-selection bug` — DONE 2026-09-08

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

### `[x] R021 — Chain identity and labelling` — DONE 2026-09-08

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

### `[x] R022 — Remove nA == nB assumptions` — DONE 2026-09-08

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

### `[x] R023 — Detect and display the assembly type` — DONE 2026-09-08

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

### `[x] R024 — Register heterodimer fixtures` — DONE 2026-09-08

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


**R020 outcome, 2026-09-08.** Fixed structurally rather than by convention: `AFDBPrediction`
now **sorts `entries` by `chainId` in `__post_init__`**, so arrival order is discarded before
any caller can observe it. `primary_entry` is gone entirely — grep across `src/` and
`notebooks/` returns no matches — along with its two notebook call sites.

**The bug, measured.** Eight consecutive endpoint calls per accession returned mixed
`['A','B']` / `['B','A']` orders for both fixtures. Before the fix, the heterodimer report
flipped five identity fields wholesale between runs of an unchanged notebook: UniProt Q96AZ6
vs P63166, organism *Homo sapiens* vs *Mus musculus*, protein ISG20 vs SUMO1, monomer length
181 vs 101.

**Determinism proven.** Ten live fetches per fixture across both arrival orders produced one
distinct report each, by SHA-256. Independently re-verified by the orchestrator over six more
fetches: one distinct hash, matching the agent's.

**Report design satisfies D3 with no user-facing switch.** Chains resolving to the same
identity collapse into one block, so a homodimer does not repeat itself; a heterodimer expands
to two. Same code path either way.

**New module API:** `format_metadata_report()` (a pure function, which is what made the
determinism assertion possible), `chain_field()`, `shared_field()` (raises when entries
disagree, rather than trusting entry 0), `describes_chain()`.

**Two incidental fixes.** `Gene` and `Model version` previously printed `N/A` on every run
because the live service sends `gene` and `latestVersion`, not the `geneNames` and
`modelVersion` that `CLAUDE.md` documents. Both now resolve. **R093 must correct the
`CLAUDE.md` field table**, which is wrong about both names.

**Missing-chain handling is now explicit**, exercised live: a structure chain the endpoint did
not describe prints a note under its own header and downgrades the assembly line; a metadata
chain absent from the structure is flagged; `entry_for_chain` on an unknown chain raises with
the list of chains that were described.


**R021 outcome, 2026-09-08.**

**Chain-id mismatch now fails rather than guessing, and says why.** `verify_chain_identity()`
reconciles all three sources (structure, PAE, pLDDT) in one call, and
`verify_document_agreement()` closes the PAE-vs-pLDDT edge directly instead of inferring it
through the structure. The refusal to map positionally is argued in the error text itself:
neither source carries a trustworthy order (the documents' `chains` arrays are re-sorted here
precisely because arrival order is not authoritative, and R020 measured the prediction
endpoint's non-determinism), and chain length is no tie-breaker because a homodimer's chains
are equal by definition. A positional guess would mis-slice every quadrant into
plausible-looking wrong scores rather than an error. Metadata-side gaps stay non-fatal notes,
since they cannot misalign a slice.

**Real protein names throughout.** Heterodimer axis labels went from `Chain A` / `Chain B` to
`ISG20 (A)` / `Sumo1 (B)`, with the full form
`Chain A — Interferon-stimulated gene 20 kDa protein (Q96AZ6)` where space allows. Three label
forms (`token`, `short`, `full`) let each context pick what fits. The chain id is retained in
every form, so a homodimer's two chains stay distinguishable. Fallback ladder verified across
eight cases down to a bare `Chain A` when nothing is known; the homodimer's service-truncated
`'3-hydroxydecanoyl-'` is marked with an ellipsis rather than discarded.

Scores unchanged on both fixtures; 191 doctests pass.

### `[x] R025 — Make USE_LOCAL_FILE mode honest about missing PAE/pLDDT` — DONE 2026-09-08

**Found during R021, pre-existing, not caused by it.** In `USE_LOCAL_FILE` mode the upload
cells print "PAE file not uploaded — PAE-dependent analyses will be skipped", but nothing
skips them: `pae_raw` / `plddt_raw` stay `None` and the parse cell raises `TypeError` before
any verification runs.

**Why it matters more than it looks.** `threshold-reference.md` §1 rests on the local-file path
being the case where the traffic light does real work: for an AFDB accession the model has
already passed the 0.6 filter, so red essentially never appears. The one route that exercises
the full range is the one that crashes.

**Options.** Either genuinely skip the PAE-dependent sections with a clear banner listing what
was skipped, or fail immediately at upload time with a message saying PAE and pLDDT are
required. Do not keep advertising a skip that does not happen.

**Done when.** Local-file mode with no PAE/pLDDT either completes with a clear statement of
what was skipped, or fails at the point of upload with an actionable message. Verified by
running it both ways.


**R022 outcome, 2026-09-08.** Audited 17 sites systematically rather than working the plan's
list, which had gone stale: **three of the four sites the plan named were already fixed** by
R011-R021, and four real issues the plan did not list were found.

**Already fixed incidentally:** the `score_masks` shape mismatch (R013 derives the mask from
`block.shape`, and R021's verification makes the two length sources a checked precondition);
views 3 and 4 indexing per-residue arrays with structure residue lists (R014 takes and
length-checks `values_x` and `values_y` separately); the profile axes (already `sharex=False`).

**Found and fixed, not in the plan:** the coverage bar clipped residue 0 (`xlim` started at 0,
not -0.5); cell 12 printed the interface residue range for one chain only; the Section 3
markdown said "for a homodimer, the matrix has four quadrants", implying four equal blocks;
and the notebook H1 said "Homodimer Confidence Metric Diagnostic Notebook", false for the
heterodimer path D3 requires. H1 is now "Dimer Confidence Metric Diagnostic Notebook" — **D6
respected, the filename is unchanged.**

**Coverage bars fixed for unequal chains.** Residue counts in the tick labels, an end cap on
each track, and an annotation on the shorter track naming where it ends, so a 69-residue chain
on a 563-wide axis no longer reads as truncated data.

**New stress fixture: `AF-0000000204661110`**, Sptlc3 (563 aa) + Gm6993 (69 aa), **8.16:1**
versus 1.79:1 for the standard heterodimer. Found via the search endpoint by ranking 100
hetero dimers by length ratio. Runs end to end with no crash, no wrong number and no shape
error. Verified independently: lengths 563/69, labels `Sptlc3 (A)` / `Gm6993 (B)`,
ipSAE_d0res 0.7120, 42 contact pairs. Register it in R024.

**Gap worth knowing:** both real heterodimer fixtures put the *long* chain first, so `nx < ny`
never occurs live. Only synthetic cases (30+400, 5+300) cover the reverse orientation. R024
should look for a short-first fixture.

**Residual, deliberately not fixed here.** At 8.16:1 the `aspect='auto'` distortion stops being
cosmetic: the 563x69 mask panels render as near-squares, hiding the block's true proportions.
Axis labels stay correct, so nothing is false. Two of the three sites are already owned
(R050 for the PAE heatmap, R052 item 4 for the mask panel). **The third, the contact map's
left panel, is owned by no task — fold it into R050.** Fixing it needs a figure-size change,
which R022 was forbidden from making.


**R023 outcome, 2026-09-08.**

**Measured fact worth keeping: every AFDB complex is a dimer.** Verified independently on the
search endpoint — `isComplex:true` and `oligomericState:dimer` both return **2,010,763**;
`trimer` and `tetramer` return nothing. Homo 1,930,523 + Hetero 80,240 sums exactly. So >2
chains is reachable only through a local mmCIF, which scopes D4's future N-chain work: it is a
local-file feature, not an AFDB one. R093 should record this in `CLAUDE.md`.

**Three failure gates, each where its evidence first exists.**
1. `fetch_afdb_metadata` raises `AccessionLookupError`, distinguishing malformed (HTTP 400,
   quoting the service's own message) from absent (404).
2. A cheap metadata gate in cell 6, before three downloads. Deliberately **not** authoritative:
   it refuses only when the declaration *and* the endpoint's own per-chain entry list agree.
   A declared `trimer` on a two-entry record is exactly the disagreement this task surfaces, so
   refusing on it would prefer metadata over the model on evidence not yet gathered.
3. The structural gate inside `verify_chain_identity`, which always fires, online or local.

The agent rejected a standalone `validate_input()` the notebook calls, on the grounds that for
a non-expert hitting Run All a check that can be skipped is a check that will be.
`verify_chain_identity` is the one call every downstream cell already depends on, so there is
no path to a score that bypasses it.

**Disagreement is shown, not resolved silently.** Both sources print, and the rule for labels
("the chains win, because that is what the scores are computed from") is stated inside the
conflict text. When chains are present but unidentifiable, the homo/hetero call is downgraded
to plain `Dimer` rather than repeating an unverifiable declaration.

**Namespace hazard fixed that R020 did not have:** chain identity keys now come from one
namespace for all chains, never mixed. Mixing would make two copies of one protein look like
two proteins, manufacturing the exact false "Hetero" this task exists to catch.

**Corrects a stale `CLAUDE.md` edge case.** Its table says a monomer accession "produces
zero-length chain B; scores will be 0.0 (no informative error)". Wrong before and after: the
real pre-fix behaviour was a bare `IndexError` at `pae.chain_ids[1]`. Now it is
`UnsupportedAssemblyError` at cell 6, before any download. Added to R093.

232 doctests pass, up from 191. All three fixtures execute end to end with scores unchanged.


**R024 + R025 outcome, 2026-09-08. Milestone M3 complete.**

**Major finding: AFDB publishes its own IPSAE-derived scores, and we match them.** The search
endpoint (`q=modelEntityId:<acc>`) returns per-direction production values:
`complexPredictionAccuracy_ipsae_AB/_BA`, `_ipsae_d0chn_*`, `_ipsae_iptm_d0chn_*`,
`_pDockQ2_*`, `_LIS_*`, `_ipsae_n0dom_*`, plus the parameters used.

Verified independently on FX-006, nine directional values:

| Value | Ours | AFDB | Delta |
|-------|------|------|-------|
| ipSAE d0res A->B | 0.555450 | 0.555450 | 5e-7 |
| ipSAE d0res B->A | 0.705718 | 0.705718 | 1e-7 |
| ipSAE d0chn A->B | 0.805532 | 0.805532 | 3e-7 |
| ipTM_d0chn A->B | 0.668714 | 0.668714 | 2e-7 |
| ipTM_d0chn B->A | 0.772954 | 0.772954 | 4e-7 |
| pDockQ2 A->B | 0.705404 | 0.705403 | 5e-7 |
| pDockQ2 B->A | 0.685271 | 0.685271 | 0 |
| LIS A->B | 0.608776 | 0.608774 | 2.2e-6 |
| LIS B->A | 0.592070 | 0.592071 | 6e-7 |

**Worst delta 2.2e-6, roughly 450x tighter than the ±0.001 tolerance.** And `n0dom` matches
per direction: **250 (A->B) and 252 (B->A)**, exactly the asymmetry R003 was written to fix.
AFDB's own pipeline confirms both R003 and R007 were right: it publishes pDockQ2 and n0dom per
direction because they genuinely are directional.

**This changes R092.** The plan assumed `ipsae.py` run locally was the only reference. AFDB's
production values are a second, independent one, available for every released entry with no
subprocess. R092 should use both. Only `pDockQ` remains unconfirmed beyond 2 dp, since that is
all AFDB publishes.

**Two conventions to document, not discrepancies:** AFDB rolls LIS up as the max while
`ipsae.py` (and this module) use the mean; AFDB's `ipsae_dist_cutoff` of 15.0 governs different
counts than our contact `DIST_CUTOFF` of 8.0.

**R022's short-first observation was an artefact, now corrected.** Sampling 100 consecutive
heterodimers: **41 short-first, 57 long-first, 2 equal**. Roughly 41% of AFDB heterodimers put
the shorter chain first, so the reverse orientation is common, not rare. `AF-0000000211026350`
(Rbx1 108 + CUL3 768, 1:7.11) is registered as FX-010 and runs clean.

**Fixture manifest grew 5 -> 10.** FX-006 primary heterodimer, FX-007 asymmetry stress
(8.16:1), FX-008 borderline (the only LOW-CONFIDENCE band value), FX-009 ipTM length dilution
(ipTM 0.5769 amber vs ipSAE 0.7119 green) and the only fixture with backbone clashes, FX-010
short-first and metric disagreement. **FX-003, the metric-disagreement fixture that had been
blocking since the original spec pack, is closed** — FX-010 supplies it from live data.

**R025 chose to fail at upload rather than skip.** Six of seven values need the PAE matrix and
the seventh needs pLDDT, so a "skip" would leave a contact count and a contact map: a different
notebook with a banner on it, not this one. The error names every missing file at once and
tells the user where to get them. Local-file mode with all three files reproduces the online
scores exactly.

**Defect recorded for R080:** FX-008 and FX-009 show Section 7 printing "consistently HIGH
confidence across all metrics" while a traffic light is amber.

## W3 — Section 2, Interface Detection

> User note: *"we should maybe add a 3D structure of residues involved with interaction.
> Currently it nicely shows the 2D view, via a contact matrix and an annotation track."*

### `[x] R030 — Add a 3D interface view to Section 2` — DONE 2026-09-08

**What.** Insert a MolViewSpec view directly after the contact map showing the interface
residues in 3D, so the 2D contact matrix and the 3D reality sit side by side.

**Supporting notes.** This is essentially current cell 24 View 1, moved forward. Decide
whether Section 6's View 1 then becomes redundant or stays as the "here is the complex"
orientation shot. Depends on R014 and R071 (both-chains side-chain colouring).

---

## W4 — Prose pass

### `[x] R040 — Em dash sweep` — DONE 2026-09-08

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


**R008 + R062 outcome, 2026-09-08.**

**The directional breakdown turned out to matter far more than the fixtures first suggested.**
The heterodimer's 0.15 ipSAE spread was not the extreme case. On FX-010 (`AF-0000000211026350`,
Rbx1 108 + CUL3 768) the spread is **0.5685**: ipSAE_d0res is 0.6888 read one way and 0.1203
read the other, and **5 of 6 directional scores are flagged**. Verified against AFDB's own
published fields: 0.688783 / 0.120251, matching to 4dp. A reader who does not know the score is
a max over directions would read 0.6888 as "this interface is fine" when one chain's view of
the other is 0.12.

The report states the consequence in band terms, generated from `traffic_light` rather than
canned: "The gap crosses a band edge: CONFIDENT one way, BELOW AFDB THRESHOLD the other."

**`DIRECTION` is a viewing control, not a scoring parameter.** It accepts `None`/`'xy'`/`'yx'`,
positional `'ab'`/`'ba'`, or real chain ids, and a typo raises rather than silently showing the
wrong direction. Verified: the `Score Results` block is byte-identical between
`DIRECTION=None` and `DIRECTION='yx'`, so the headline numbers keep matching `ipsae.py` and
AFDB.

**pDockQ is reported as symmetric with a dash, never a zero delta**, on the grounds that a zero
would suggest a difference had been measured. On the homodimer the symmetry of the other scores
is stated as measured rather than assumed, with both columns shown so it can be checked.

**R062: the profile plot now marks the argmax.** Machine-verified on all four fixtures that the
starred residue's value equals the reported score for every series. That was the missing link —
the reported score *is* one residue's value, and nothing on the plot said so.

**`mean_ptm_by_residue` deliberately left off**, with a good argument: every other series on the
figure obeys "reported score = max over the profile", which is exactly what the new star
asserts. pDockQ2 pools over contact *pairs*, so it is neither the max nor any function of the
per-residue array; starring it would reintroduce the "is this correct?" confusion R062 exists to
remove. It is also `NaN` for non-interface residues (86 of 344 residues are interface on the
homodimer, 22 of 282 on the heterodimer), so it would render as a mostly-absent broken line.
Reserved for R074's 3D view, where the contact set is visible.

250 doctests pass. Four fixtures execute end to end, plus a fifth run with `DIRECTION='yx'`.


**R061 outcome, 2026-09-08. Milestone M4 complete.**

Section 4 is now 4.1 ipTM_d0chn, 4.2 ipSAE (all three variants together, since they differ
only in what `L` is fed to `d0_func`), 4.3 pDockQ, 4.4 pDockQ2, 4.5 LIS, 4.6 combined. Each
subsection: what it measures in plain language, the formula, the reasoning from the
publication, how to read the number, then the computation.

**56 publication claims, every one page-cited.** Eight of the load-bearing ones were
independently spot-checked against the PDFs and all were found on the exact page claimed,
including Dunbrack's RAF1/RIPK1 walk-down (ipTM 0.290 -> 0.459 with the cutoff -> 0.044 with
per-residue d0), Bryant's "AUC of 0.95" for the product and "average error of 0.11", Zhu's
"large, highly confident incorrect interfaces", and the AFDB note's 63% -> 18% homodimer
false-positive figure.

**Seven claims could not be sourced and are labelled as such in the prose rather than
asserted**, including the LIS linear transform `(12 - PAE)/12`: Kim says only "inversely
mapping PAE values to a 0-to-1 scale", so the linear form comes from the reference
implementation, not the paper. The notebook says so.

The agent also declined to make six claims earlier drafts had suggested, among them that
pDockQ2 supersedes pDockQ, and any authorial rationale for LIS averaging rather than taking a
max. Refusing to source-dress an inference is the behaviour this task needed.

**Both known text defects fixed:** the backwards `ipSAE_d0chn` description ("Conservative
variant. Penalises small interfaces more" -> "Most permissive variant"), and the
`plot_pae_score_masks` panel title that still said `ipTM` after the R004 rename.

**Numerical equivalence at full float precision** on both fixtures: all seven values identical
to 16 significant figures. Four fixtures execute end to end, 26 code cells each, zero errors.

**Open judgement call for R080.** Each compute cell now prints its band label and threshold
provenance next to the number, which answers "how to read this" where the reader is looking.
It reads the same canonical `THRESHOLDS` as Section 7, so they cannot disagree, but the band
appears twice. If Section 7 should own it exclusively, remove the `traffic_light` lines in
cells `cell-016b/f/h/j` and the loop in `cell-016d`.

## W5 — Section 3, PAE Matrix Decomposition

### `[x] R050 — Fix the full PAE heatmap size, and the contact map's aspect` — DONE 2026-09-08

**Added by R022:** `plot_interface_contact_map`'s left panel also uses `aspect='auto'` and is
owned by no other task. At 8.16:1 (fixture `AF-0000000204661110`) it distorts badly. Fixing it
needs a figure-size change, which is why R022 left it. Handle it with the heatmap below, since
both are the same class of fix.


> User note: *"reduce the size of the PAE matrix (currently you have an in-cell scroll, which
> makes it really hard to view it all at the same time)."*

**Supporting notes.** Cell 15 uses `figsize=(10, 9)` at `figure.dpi = 150` (set in cell 2),
producing a 1500×1350 px image that overflows the output area and triggers the scrollbar.
Fix by reducing the figure size and setting `aspect='equal'` so the matrix stays square.
`aspect='auto'` is currently used, which also distorts the matrix for a heterodimer.

### `[x] R051 — Green default palette with a switchable alternative` — DONE 2026-09-08

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

### `[x] R052 — Rebuild the score-mask panel` — DONE 2026-09-08

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


**M5 outcome, 2026-09-08 (R050 + R051 + R052, batched as one visual pass).**

**Sizing.** Full PAE matrix 1500x1350 -> **900x810**, which needs no downscaling in a ~900px
notebook output column. Sized against dpi rather than lowering dpi, so text keeps its size.

**Aspect fixed, and it was badly wrong before.** Every mask panel previously rendered at ~0.757
*regardless of the data* — an 8.16:1 block and a 1:7.11 block came out the same near-square.
Now exact on all four fixtures: 1.000 / 1.792 / 8.159 / 0.141. Panel widths equal at
302/163/55/647 px respectively, against 929/929/929/**743** before, where the fourth panel was
20% narrower because the colour bar was attached to it alone. Fixed with a dedicated colour-bar
column in the gridspec.

**The grid now follows the data**: 2x2 at low asymmetry, a single row of tall slivers at
8.16:1, and the user's suggested vertical stack at 1:7.11. A 2x2 of 8.16:1 blocks is mostly
whitespace.

**Palette.** Default is `Greens_r`, dark green at PAE 0 (luminance 0.198) fading to near-white
at high PAE (0.982) — verified, and it is the "green is hard to see against white" requirement:
confident regions are the dark ones.

**A real bug found by checking rather than assuming:** `PAE_CMAP_CHOICES` had
`'colourblind_safe': 'viridis_r'`, which puts **bright yellow at PAE 0** and inverts the
dark-equals-confident rule the whole design rests on. Corrected to `'viridis'`.

**Colour-blindness verified, not asserted.** Machado-Oliveira-Fernandes 2009 at severity 1.0.
`Greens_r` and `viridis` keep monotone lightness under deuteranopia, protanopia and
tritanopia; `RdBu_r` does not, under any of them, which is why it is kept but documented as
not recommended. The simulation was also run over the rendered figures, not just the ramps.

**The grey underlay problem was solved properly.** `#90A4AE` chosen by search over OKLab dE,
scoring 14.03 normal / 9.13 deuteranopia / 10.97 protanopia against the green ramp. But a
single-hue ramp spans the whole lightness range, so **no** neutral can separate from all of it
under every palette (plain grey collapses to 4.8-6.0 against viridis and RdBu_r). Hence
`PAE_UNUSED_HATCH = '///'`: texture is the one channel independent of palette choice, and it
also covers monochrome printing.

**Two `dataviz` skill checks were overridden on domain grounds, stated rather than
suppressed:** the green ramp's light-end contrast (a PAE heatmap is a contiguous field in a
drawn frame, and fading to near-background at the uncertain end is semantically right —
truncating the ramp would compress real dynamic range), and viridis's single-hue rule (it is
the canonical perceptually-uniform scientific ramp).

**Provenance stated for the mask panel**, as asked: a thumbnail locator of the full matrix with
the displayed block shaded, beside a sentence naming the exact slice
(`pae_matrix[:172, 172:344]`). Note the user's original wording said "top left quadrant"; it is
the **upper-right** (`[:nA, nA:]`) — rows from the first chain, columns from the second. The
locator diagram shades the actual block so the wording is not the only thing carrying it.

**Orchestrator note on verification.** Three of my own measurement approaches disagreed with
the agent and with each other. A control against plain `imshow(aspect='equal')` — correct by
definition — showed my methods were the faulty ones; the agent's
(`draw_without_rendering()` -> `apply_aspect()` -> `get_window_extent()`) reproduces the known
answer exactly. Validate a measurement against a known-good case before using it to contradict
a result.

## W6 — Section 4, Score Computation

### `[x] R060 — Apply the formula audit` — DONE 2026-09-08 (satisfied by R012 + AFDB cross-check)

**What.** Land R003, R004, R005, R006 in the notebook and the module, then verify against
`ipsae.py` on both fixtures.

**Supporting notes.** The user's ask was *"a one-over sanity check of each formula, to make
sure it makes sense."* The audit is complete and is recorded in R001. Net result: the
notebook is correct on pDockQ, pDockQ2 and LIS (where `CLAUDE.md` is wrong), and has one real
bug of its own (d0dom asymmetry) plus one edge-case divergence (`L == 27`).

### `[x] R061 — Restructure Section 4 into per-score subsections` — DONE 2026-09-08

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

### `[x] R062 — Rework the per-residue score profile plot` — DONE 2026-09-08

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

### `[x] R070 — Write supporting text for every view` — DONE 2026-09-08

> User note: *"What each plot is showing, and what we are colouring by."*

**What.** Each view gets a markdown block above it: what is drawn, what the colours encode,
what to look for, and what a problem looks like.

### `[x] R071 — Show side chains on both chains, coloured by chain` — DONE 2026-09-08

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

### `[x] R072 — Explain View 3 (ipSAE d0res on the interface)` — DONE 2026-09-08

> User note: *"why are we using the d0res score? Not that it's wrong, but nothing is explained
> about it. And also, what is low, mid, and high? Need to explain the value ranges."*

**What.** Two additions. First, state why `d0res` is the variant shown: it is the AFDB primary
classifier and the only variant whose `d0` is derived per residue, which makes it the one that
varies meaningfully *along the chain* and therefore the only one worth painting onto a
structure. Second, replace "red=low, yellow=mid, green=high" with the actual numbers, tied to
the R002 thresholds, and add a colour bar or an explicit legend so the mapping is readable.

### `[x] R073 — Rebuild View 4 (disagreement)` — DONE 2026-09-08

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

### `[x] R074 — Add pDockQ2 views` — DONE 2026-09-08

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

### `[x] R075 — Fix the MolViewSpec implementation issues` — DONE 2026-09-08

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


**M6a outcome, 2026-09-08 (R071 + R073 + R074 + R075 remainder). Six views now, all on both
chains, all with a legend.**

**The user's suggested side-chain colours would have failed, and the agent measured it rather
than shipping them.** Cornflower on chain 1 + orange on chain 2 gives **dE 7.3** under
deuteranopia, because coral simulates to `#BBA83F` and orange to `#D4BE29` — chain 2's side
chains would be invisible against chain 2's own cartoon, the single most important contrast in
the view. Fix is structural: chain 2's cartoon is warm so its side chains must be cool; chain
1's teal desaturates to neutral grey so its side chains must be strongly chromatic. Cornflower
is **kept but moved to chain 2**, chain 1 takes gold `#FFD400`. Worst pair 7.3 -> **15.1**,
every within-chain contrast >= 29.9.

The category palette was worse still: the old green/blue/red triple scored **dE 6.0** under
deuteranopia between *agree* and *contact-only* — the two categories a reader most needs to
separate were the two that merged. Now Okabe-Ito plus a luminance-separated dark red, worst
pair 21.5, a 3.6x improvement.

**The disagreement threshold is now sourced, not invented.** `MVS_DISAGREEMENT_THRESHOLD` is
`THRESHOLDS['ipsae_d0res'].amber` = 0.60, provenance PUBLISHED, replacing the hard-coded 0.5.
Verified identical.

**pLDDT ladders reconciled by aliasing, so they cannot drift again:**
`MVS_PLDDT_BANDS is PLDDT_BANDS` is literally `True`. The 2D ladder won because it is
AlphaFold's published one (strict `>`). Verified: all nine edge values agree, including 100.00
and 90.00. Consequence — the old 3D ladder **did not draw a residue at exactly pLDDT 100.00**;
View 2 now paints 172/172 rather than 171/172.

**A genuine scientific finding from the new pDockQ2-vs-ipSAE view.** The `ipsae_only` category
is empty on all four fixtures, and not from a bug: across the 102 confident interface residues,
**not one has a mean contact PAE worse than 2.61 A**, against a 10 A cutoff (minimum contact-ptm
0.9362). *(Corrected 2026-09-08: this note first said "1.8-2.6 A", which was the range of the
per-fixture maxima, not the span over residues. The true span is 0.85-2.61 A. Caught by the M6b
agent, which declined to repeat the figure and wrote the notebook prose as a bound instead;
re-measured independently.)* **pDockQ2 is never stricter than ipSAE at the
residue level** — its reputation for strictness comes entirely from the pooling (mean pLDDT x
mean ptm through a steep sigmoid), not from contact quality. Worth saying in R070/R072's prose,
since it corrects a natural misreading.

**The directional spread is now localised in 3D.** On FX-010 (0.5685 spread) View 6 shows the
0.69 direction with 7 confirmed residues on Rbx1 and the 0.12 direction with **zero** on CUL3,
while 49 of CUL3's 60 interface residues are "pDockQ2 only".

**A trap avoided:** `build_interface_value_view` now drops non-finite values, because
`mean_ptm_by_residue` is `NaN` for non-interface residues and a colormap paints "bad" values
**opaque black** — the most emphatic colour in the scene would have marked the residues with no
measurement at all.

283 doctests (was 254). All seven scores bit-identical on all four fixtures. Notebook runs
clean with six Mol* views.

## W8 — Section 7, Diagnostic Summary

> User note: *"the traffic light system is great. Just need to make sure we are using the
> right values or scores when describing each metric."*

### `[x] R080 — Wire the summary to the canonical thresholds` — DONE 2026-09-08

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

### `[x] R081 — Audit the metric descriptions` — DONE 2026-09-08

**Supporting notes.** Cell 26's `descriptions` dict is mostly sound but has issues to check
against R001: "ipTM: All PAE cells used, no cutoff" is correct but should say `ipTM_d0chn`
(R004); "pDockQ: Contact count x interface pLDDT" understates it (`log10` of the contact
count, and the pLDDT is a mean over the union of contacting residues in both chains);
"ipSAE_d0chn: Penalises small interfaces more" needs checking — a larger `L` gives a larger
`d0`, which is *more* lenient, so the direction of that claim may be backwards.

### `[x] R082 — Review the agreement matrix` — DONE 2026-09-08

**Supporting notes.** Cell 27 normalises each score by its own green threshold and caps at
1.5, then plots pairwise absolute differences. That makes the matrix depend entirely on the
R002 thresholds being right, and the cap means two scores that are both far above threshold
register as perfect agreement. Decide whether this is the comparison worth showing, and label
the axes with what the normalisation actually is.

---


**M6b outcome, 2026-09-08 (R070 + R072).** Section 6 restructured from one header plus a loop
into header, setup, and six (markdown, code) pairs, so each view's explanation sits directly
above it. Cell count 42 -> 54. Module untouched.

**View 3's `d0res` choice is now argued, not asserted.** Two reasons, both checkable: it is the
variant AFDB publishes as `ipSAEmax`, and it is the only one whose `d0` varies per residue, so
`d0chn` painted on a structure would show the same shape at a more forgiving scale. The vague
"red=low, yellow=mid, green=high" is replaced by a seven-row table mapping value -> exact ramp
hex -> AFDB band, verified against `value_colours` rather than assumed.

**A caveat the agent added unprompted and was right to:** the ramp midpoint (0.50) is not a
band edge, so the colours are a continuous scale and the AFDB bands are not colour breaks. A
reader would otherwise reasonably assume the green/yellow transition meant something.

**Every factual claim carries a source**, and three are explicitly marked as reasoning rather
than citation: the residue-level `d0` ordering (derived from the definitions, since the notebook
states the theorem only at score level), "steep sigmoid" as qualitative, and the
what-to-look-for guidance.

**The M6b agent corrected an error in my own M6a note above**, declining to repeat a figure it
could not reproduce. Re-measured independently: 102 confident interface residues, span
0.85-2.61 A, minimum contact-ptm 0.9362. The note is fixed. Worth recording as a pattern: the
agent wrote the notebook prose as a bound ("not one worse than 2.7 A") precisely because it
could not verify the range it had been handed.

## W9 — Validation and documentation

### `[x] R090 — Homodimer end-to-end run` — DONE 2026-09-08
Fixture FX-001 `AF-0000000065889468`. Full run, clean, no warnings.

### `[x] R091 — Heterodimer end-to-end run` — DONE 2026-09-08
Fixture `AF-0000000211034637` (R024). Full run, every section correct for unequal chains.

### `[x] R092 — Numerical agreement with ipsae.py` — DONE 2026-09-08
All seven values within ±0.001 of `ipsae.py` v4 on both fixtures. Record in
`specs/homodimer_diagnostic/validation-report.md`, which T038 in the original `tasks.md`
already calls for and which does not yet exist.

### `[x] R093 — Sync CLAUDE.md and the spec pack` — DONE 2026-09-08
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

### `[x] R095 — Delete the superseded code (last task, per D10)` — DONE 2026-09-08

**What.** Only once R090-R093 pass. Remove `src/insightfold/interface.py` (absorbed by R015)
and the duplicated inline `extract_cb_coords` in `notebooks/analysis_template.ipynb`, and
fix the four documents that reference `interface.py` by path.

**Why deferred.** Deleting early makes any later failure ambiguous: a real bug and a missing
import look identical from the traceback. Nothing imports `interface.py` today, so leaving it
in place costs nothing.

### `[x] R094 — Close out the original tasks.md` — DONE 2026-09-08
Mark T004 (threshold bands), T033-T038 (validation), and T039-T041 (review) against the work
done here, or point them at this file.

---


**R030 outcome, 2026-09-08. Milestone M6 complete.**

**Decision: the chain overview MOVED to Section 2 rather than being duplicated there**, keeping
its number so cross-references and `MVS_VIEW_LABELS` stay valid. Section 6 now draws Views 2-6
and says in bold where View 1 went. The module was not modified.

**The argument that settled it was a real constraint, not tidiness.** Each view is a live Mol*
instance holding a WebGL context; browsers cap concurrent contexts (16 in Chrome and Firefox)
and silently blank the oldest past the cap. The notebook already sits at six. A seventh showing
a picture already on the page would spend headroom, an extra BinaryCIF download and an extra
Mol* boot for no new information. **Verified: 6 iframes before, 6 after.** The executed-notebook
size delta equals the source delta exactly (+3,776 B), which is arithmetic proof no iframe was
added.

The reader argument is as strong: the contact map and the 3D view are the *same two arrays* in
two coordinate systems, and the useful check — "several separated blocks in the map may still
be one physical patch" — is a comparison. Fifteen cells apart it is a memory exercise. Placing
the copy next to the map would have left the original doing nothing new.

**The agent corrected its own draft against the fixtures.** Its first version said a short
residue range on a long chain means a small patch. FX-007 falsifies that (Sptlc3: 21 interface
residues spanning residues 65-528 of 563), as does FX-010 (CUL3: 60 spanning 441-768). The text
now says the range is an outer bound, not a size, and points the reader at the count instead.
That is a claim that would have read as authoritative and been wrong on half the fixtures.

**Local-file mode gained a sentence Section 6 did not need:** that skipping the 3D views costs
nothing numerically. At Section 2, where the interface is being defined, a reader needs telling
that the contact map and every score are unaffected.

54 cells, 32 code cells, 0 errors on all four fixtures, 283 doctests.


**M7a outcome, 2026-09-08 (R080 + R081 + R082 + R002b/R004 remainders).**

**The self-contradiction is fixed structurally, not by patching a conditional.** FX-008 printed
"consistently HIGH confidence across all metrics" while its own table showed
`ipSAE_d0res 0.6366` in amber — twice in four lines. Cause: `if n_green >= 4` over five scores.
Now one `ConfidenceSummary` object is built by `summarise_scores()`, and the table, the margin
plot and the prose are three renderings of it. There is no second path that could reach a
different verdict. Verified on all six fixtures: table and prose agree, and the phrase is gone.

**The agent caught me repeating an error I had already corrected.** My brief told it to state
that 0.6 came from "a data-driven analysis of ~31 million candidate complexes".
`threshold-reference.md` §0.2 records that this is wrong — 0.6 was **adopted by citation** as
community-established, **validated** on 230/117 and 94/250 post-training benchmarks, and
**applied** at ~31M scale. The agent wrote the accurate version and said why. The plan was
right; my brief regressed against it.

**R082: the agreement matrix is replaced, with an argument.** `|ai - aj|` over scalars is a
distance matrix over points on a line: 25 cells re-encoding 5 numbers. Worse, `cap=1.5` made
two scores comfortably past threshold register as agreeing *perfectly* — so the figure was
closest to blank on the AFDB path, which is the notebook's main path. Replaced by
`plot_threshold_margins`: five bars on one axis, each normalised by its own green threshold,
where the gap between rows *is* the pairwise difference. It restores which side of the cutoff
each score falls on and how far past, with provenance in the tick label. The loss (one
glanceable cell per pair) is stated in the notebook rather than hidden.

**R081 found a second wrong description, not just the known one.** `LIS: "Density of
inter-chain PAE < 12"` is wrong — the mean runs only over sub-cutoff cells, so the count never
enters, and a single cell at PAE 0 scores 1.0. It measures how far below the cutoff the
confident cells sit, not how many there are. Also corrected: pDockQ ("+" should be a product of
mean pLDDT and log10 of contact **pairs**), pDockQ2 (a product, and the PAE term is a mean of
the TM transform, not a mean PAE), and `ipsae_d0res` ("primary AFDB classifier" -> one **half**
of the release criterion, which is the conjunction).

**An HTML injection bug found and fixed.** `PAE < 10` in a description and `0.55 < 0.60` in a
FAIL reason both contain a raw `<`, which a browser reads as an opening tag and silently eats
along with everything after it. Confirmed in rendered FX-008 output. All interpolated strings
in `format_summary_table_html` are now escaped.

299 doctests. Six fixtures execute clean. Scores numerically identical.


**R040 outcome, 2026-09-08. Milestone M7 complete.**

**80 em dashes in user-facing text, now 0, with zero deliberate keeps.** Replacements were
varied by what each sentence needed rather than substituted uniformly: colons where the clause
explains or enumerates, commas for parentheticals, parentheses for genuine asides, full stops
where two independent thoughts were being joined, and words where the glyph was doing a word's
job (`"— symmetric —"` -> `"symmetric"`, table `"—"` -> `"n/a"`).

**A `&mdash;` HTML entity that no Unicode grep would find.** The AFDB verdict banner in
`format_summary_table_html` carried `&nbsp;&mdash;&nbsp;`. It surfaced only by diffing rendered
HTML before against after. `&#8212;` and `&#x2014;` were also checked: none.

**Fixed-width alignment preserved.** The pDockQ row of the directional table had
`— symmetric —` inside `f"{sym:>{col+2}s}"` fields. Machine-checked: the `Reported` column ends
at character 74 before and after, and the row is now the same length as its neighbours rather
than two characters short.

**Numeric identity proven properly:** every stream, execute_result and HTML output captured for
all six fixtures before and after, MolViewSpec base64 blobs stripped (they embed a wall-clock
timestamp), then every numeric token compared in order. 483 to 507 tokens per fixture, all
identical.

**Follow-up the agent flagged and the orchestrator then did.** Five bullets in the intro used an
**en** dash as a prose dash (`- **LIS** – Density of…`), which A2's scope left alone because
they are not em dashes but which read identically to a user. Converted to colons; the 26
genuine numeric-range en dashes (`70–90`, `pp. 5–6`, `CB–CB`) are untouched. That edit initially
dropped the trailing newlines and ran the bullets together, caught by re-reading the rendered
cell and repaired.


**M8a outcome, 2026-09-08 (R090 + R091 + R092). Delivered
`specs/homodimer_diagnostic/validation-report.md`, 820 lines: the report the original spec
pack's T038 called for and never got.**

**13 runs across all ten registered fixtures. Zero errors, zero bytes to stderr on every
successful run.** Seven scoring fixtures at 3.9-5.3 s each, six figures and six Mol* views
apiece. Both negative fixtures halt cleanly at the intended gate. Local-file mode with all
three documents reproduces the online scores bit-identically; with documents missing it names
every missing file.

**231 numerical comparisons against two independent references, zero outside tolerance.**
vs `ipsae.py` v4 as a subprocess: worst |delta| **4.7e-5**, every delta above 1e-6 falling on
the three values the reference prints at 4dp; the 6dp ones agree to <= 6e-7. vs AFDB production
values: worst |delta| **2.4e-6** on full-precision fields, exact on 2dp, and `n0dom` matching
integer-for-integer per direction on all five heterodimers.

**Three findings for later tasks.**
1. The R060 methodology note was over-generalised; corrected above.
2. **`complexPredictionAccuracy_ipTM` is AlphaFold's own ipTM, not our `ipTM_d0chn`.** It
   tracks `_iptm_af` and is up to 0.054 away on FX-010. Substituting it produces a false
   failure on 6 of 7 fixtures. Nothing in the spec pack warns about this: R093 must record it.
3. **FX-002 is not borderline** at `ipSAE_d0res` 0.7699 (CONFIDENT); FX-008 (0.6366) holds that
   role. R094 should correct the manifest.

**Two manifest items closed:** the R080 self-contradiction no longer reproduces, and `pDockQ`
is reconciled beyond 2dp (<= 3.3e-5 on all seven), closing the last numerical gap.

**No bug found** in the notebook, module or `interface.py`.

**Principal outstanding risk: Colab is entirely unverified.** Every run had `IN_COLAB = False`,
so the clone, stale-clone-refresh and `pip install molviewspec` paths never executed and the
60 s budget is unmeasured. The bootstrap's `TODO(merge)` still pins the feature branch, so a
Colab test today would not predict post-merge behaviour anyway.


**M8b outcome, 2026-09-08 (R093 + R094 + R095). Milestone M8 complete; the rework is closed.**

**R093 — `CLAUDE.md` rewritten, not patched.** Fifteen verified errors were corrected. The
structural change matters more than any one of them: the file no longer carries its own copies
of the seven formulas or the threshold table. Duplication is what produced every error in the
list, so the fix is to make `formula-reference.md` and `threshold-reference.md` the named
authorities and leave `CLAUDE.md` as a map plus the traps. What it keeps is what a copy cannot
drift into being wrong about: which authority owns what, the *name* of each trap, and the facts
that live nowhere else (the API field names, the endpoint's non-determinism, the module layout,
D4/D8/D9/D10).

Corrections, all sourced:

| # | `CLAUDE.md` said | Now | Authority |
|---|---|---|---|
| 1 | pDockQ2 sigmoid `0.715 / -12.3 / 0.605` | `1.31 / -0.075 / 84.733`, in the module | `formula-reference.md` §6, `ipsae_v4.py:695` |
| 2 | pDockQ2 pLDDT weighted by contact multiplicity (`np.repeat`/`np.tile`) | unweighted mean over unique interface residues | `formula-reference.md` §6 (D6), `ipsae_v4.py:692` |
| 3 | pDockQ `n` = interface residue count | contact **pair** count, and the two names are separated (`n_interface_residues` vs `n_contact_pairs`) | `formula-reference.md` §5, `ipsae_v4.py:653, 664` |
| 4 | LIS combined with `max` | **mean** of the two directions | `formula-reference.md` §7, `ipsae_v4.py:982` |
| 5 | `d0dom` = rows of `pae_AB` + rows of `pae_BA`, one value for both directions | per direction: rows-with-any + cols-with-any of that direction's own block | `formula-reference.md` §4, `ipsae_v4.py:751-756` |
| 6 | `compute_ipsae` a generator, prose promising a flat dict, no per-residue arrays | removed; `ciu.compute_ipsae` returns `IPSAEResult` with all three variants, both directions and the profiles | `formula-reference.md` D7 |
| 7 | only the scalar `d0_func` | both helpers named, with the `L == 27` divergence (1.0 vs 1.038891) and which score uses which | `formula-reference.md` §0 |
| 8 | `THRESHOLDS` = `ipsae 0.6/0.4, iptm 0.7/0.5, pdockq 0.23/0.09, pdockq2 0.5/0.23, lis 0.15/0.09` | the canonical table with provenance per number, plus the four published ipSAE bands | `threshold-reference.md` |
| 9 | "Primary AFDB classifier: ipSAE_d0res >= 0.6" under a Dunbrack-attributed heading | the joint criterion `ipSAE_d0res >= 0.6 AND pDockQ2 >= 0.23`, attributed to AFDB 2026, with the metric-vs-threshold citation split | `threshold-reference.md` §0, "Citation policy" |
| 10 | metadata fields `geneNames`, `modelVersion` | `gene`, `latestVersion` (both printed `N/A` on every run for months) | R020, re-verified live |
| 11 | no `assemblyType` / `oligomericState` / `complexComposition` / `isComplex` / `chainId`, no `chains[].name` | all documented, `oligomericStateDescription` marked optional (absent on both fixtures) | R021, R023, re-verified live |
| 12 | "consume `result[0]`" | one entry per chain in **non-deterministic** order; key by `chainId` | R020 |
| 13 | no `/api/search` | documented, including AFDB's own `complexPredictionAccuracy_*` production scores as a second validation reference, the two rounding regimes, the LIS max-vs-mean convention, and `ipsae_pae_cutoff = 10.0` | R024, R092 |
| 14 | — (absent) | **`complexPredictionAccuracy_ipTM` is AlphaFold's own ipTM, not `ipTM_d0chn`**, up to 0.054 apart; substituting it fails 6 of 7 fixtures | `validation-report.md` §5.7 |
| 15 | monomer accession "produces zero-length chain B; scores will be 0.0 (no informative error)" | wrong before *and* after: it was a bare `IndexError`, and is now `UnsupportedAssemblyError` before any download | R023 |
| 16 | zero contacts -> sigmoid minima 0.018 / 0.005 | `0.0`, which is what `ipsae.py` returns and what the ±0.001 criterion requires | `formula-reference.md` §8 |

Also added, because they existed nowhere a future session would look: the module layout and D8,
the D9 bootstrap with the `TODO(merge)` branch pin flagged prominently, the measured fact that
every AFDB complex is a dimer (so D4's N-chain work is a local-file feature), the Colab-unverified
risk, and the conventions section (British spelling, no em dashes, bare `Figure` returns, palette).

**The biopython contradiction is resolved by scoping, and both sides are kept.** The prohibition
is per-consumer, not repo-wide: it binds `complex_interface_utils.py` and the two dimer
notebooks, and grep confirms none of the three imports any prohibited package.
`pyproject.toml`'s heavy set belongs to the other pipelines: `Bio` is imported only by
`notebooks/protein_model_chem.ipynb`, `pandas` and `scipy` only by `src/insightfold/variants/`
and `src/pdbe_interfaces/`, `networkx` only by that same notebook. `plotly` is imported by
nothing and is flagged as a removal candidate. Nothing was deleted from `pyproject.toml`,
because deleting a dependency another pipeline uses to satisfy a rule that does not bind it
would be the wrong repair. This is also, restated, exactly why D9 forbids `pip install`.

**R094 — `tasks.md` closed, and three of its tasks were closed by results it did not predict.**
All 41 tasks are now resolved. T002's borderline blocker was closed by discovering the premise
was wrong rather than by curation; T003 got expected outputs for seven fixtures against two
references rather than one fixture against one; T039's single review pass was replaced by the
per-task protocol, which is the stronger check. T041 records a **recommendation** of beta,
conditional on the `TODO(merge)` flip and a real Colab run, and says explicitly that the
decision is the user's rather than making it.

**The FX-002 correction, and why it is worth more than a relabel.** The manifest called FX-002 a
"provisional borderline candidate" and left it at `candidate-needs-scoring` from May to
September. It scores `ipSAE_d0res` **0.7699** (CONFIDENT) with all seven values green: it is not
borderline and never was. FX-008 (0.6366) holds that role. FX-002 is re-labelled as a second
high-confidence homodimer, and kept rather than retired for a reason the original entry could
not have known: with FX-001 it is one of only two entries in AFDB's all-2-dp storage regime, and
its `pDockQ` 0.687671 -> AFDB `0.68` is the single measurement that proves AFDB **truncates**
rather than rounds (rounding gives 0.69, and `pDockQ` has no per-direction field to fall back
on). Without a second entry in that regime the rule could not be checked at all. The manifest's
"What is still open" section is now empty, with each of its four items struck through and its
resolution named.

**R095 — the deletions, and one that turned into a repair.**

`src/insightfold/interface.py` is deleted. Nothing imported it: the only occurrence of the
module path anywhere in the repo was inside its own docstring, and the notebook's
`detect_interface` call resolves to `ciu.detect_interface`. Four documents referenced it by
path and all four are fixed (`CLAUDE.md`, `specs/homodimer_notebook_extraction.md`, and both
agent-suite documents), each pointed at `complex_interface_utils.py` with the deletion noted.
The two original spec documents that say "the attached `interface.py`" are left alone: they mean
AFDB's *production* file, which is a different thing and is still accurately described.

**`analysis_template.ipynb` could not have been run since it was written.** Removing the
duplicated inline `extract_cb_coords` meant checking whether the notebook still worked, and the
check found it did not work *before*: cell 8 carried a `SyntaxError`,
`f"...{len(meta.get('sequence', '')) or 'N/A')}..."`, with an unmatched `)`. Every code cell was
parsed to confirm it was the only one. So "do not break it" was not the constraint it looked
like, and pointing the template at the module was clearly better than leaving a third divergent
parser beside a corrected `CLAUDE.md`.

What changed, and why it is one coherent change rather than scope creep: making the module
importable is the prerequisite for removing `extract_cb_coords`, and once the module is imported
the template's other inline copies are strictly worse than a call. They carried **every** defect
R093 had just corrected in `CLAUDE.md`: `THRESHOLDS` at `0.6/0.4`, `0.7/0.5`, `0.5/0.23`,
`0.15/0.09`; pDockQ on the residue count; pDockQ2 on `0.715 / -12.3 / 0.605`; LIS on `max`;
`d0dom` single-valued; `meta_list[0]`; `geneNames` and `modelVersion`. Leaving them would have
left the corrected `CLAUDE.md` contradicted by a file in the repo that a future session is meant
to copy from.

Cells 2, 3, 8, 10, 11, 13, 14, 15, 16 and 20 now call the module; the D9 bootstrap replaces the
bare `pip install`, carrying the same `TODO(merge)` marker; the blanket
`warnings.filterwarnings('ignore')` is gone. Three markdown cells that described the old
inline design were rewritten rather than left describing code that is no longer there.
Incidental repairs: the `SyntaxError`, the non-deterministic `meta_list[0]`, the two wrong field
names, the hard-coded `'bcif'` format against a possibly-empty URL, and cell ids (the file
declared nbformat 4.5 and carried none, which is a warning that becomes a hard error).

**Verified after every deletion, not before.**

| Check | Result |
|---|---|
| `import insightfold.complex_interface_utils` with `src` on the path | OK |
| `import insightfold.interface` | `ModuleNotFoundError`, as intended |
| Module doctests | **299 passed, 0 failed** (unchanged) |
| `homodimer_diagnostic.ipynb` on FX-001, `nbclient` | 33 code cells, **0 errors, 0 bytes to stderr**, 4.6 s |
| FX-001 scores after the deletion | ipSAE_d0res 0.9143, d0chn 0.9529, d0dom 0.9527, ipTM_d0chn 0.9529, pDockQ 0.6913, pDockQ2 0.9269, LIS 0.7564, 116 contact pairs, AFDB PASS - **identical to `validation-report.md`** |
| `analysis_template.ipynb` on FX-001, `nbclient` | 15 code cells, **0 errors, 0 bytes to stderr**, 3.2 s - it now runs at all, and produces the same seven values as the reference notebook |

The template agreeing with the reference notebook to four decimals is the useful result: it is
an independent orchestration of the same module, so it would not have matched if the deletion
had disturbed anything.

**One thing this milestone did not close, and cannot.** Colab is still unverified and the
bootstrap still pins `homodimer-notebook-rework`. Both are recorded at the top of `CLAUDE.md`
under a `TODO(merge)` banner and in `tasks.md` T040/T041, because the flip has to happen at
merge time and testing Colab before it would not predict post-merge behaviour anyway.


**M8b outcome, 2026-09-08. PROJECT COMPLETE.**

**`CLAUDE.md` no longer carries its own copies of the formulas or thresholds.** That
duplication caused every one of the 16 errors corrected, so the file is now a map plus the
traps, naming `formula-reference.md` and `threshold-reference.md` as the authorities. Sixteen
verified errors fixed: five wrong formulas, the whole threshold block, the `geneNames` /
`modelVersion` field names that made two fields print `N/A` on every run, the "consume
result[0]" instruction that is unsafe under non-deterministic ordering, and the monomer
edge-case row that was wrong both before and after the rework.

**The biopython contradiction was resolved by scoping, not by deleting.** The prohibition is
per-consumer, not repo-wide: it binds the module and the two dimer notebooks, none of which
import anything prohibited. `pyproject.toml`'s heavy set belongs to the other pipelines
(`Bio` only in `protein_model_chem.ipynb`; `pandas`/`scipy` only in `variants/` and
`pdbe_interfaces/`). Deleting a dependency another pipeline uses, to satisfy a rule that does
not bind it, would have been the wrong repair. `plotly` is imported by nothing and is flagged.

**`src/insightfold/interface.py` deleted** after proving nothing imported it. Four documents
that referenced it by path corrected; two spec documents left alone because they mean AFDB's
production file, which still exists.

**`analysis_template.ipynb` did not work before this task.** Checking whether removing its
duplicated `extract_cb_coords` would break it found a pre-existing `SyntaxError` in cell 8, so
the file had been unrunnable. It also carried inline copies of **every** defect R093 had just
corrected, which would have left the corrected `CLAUDE.md` contradicted by the very file a
future session copies from. It now calls the module and runs clean.

**The load-bearing verification:** the template and the reference notebook are independent
orchestrations of the same module, and both produce ipSAE_d0res 0.9143 on FX-001. They would
not agree if the deletion had disturbed anything.

**FX-002 corrected and kept**: it scores 0.7699 (CONFIDENT), not borderline, so FX-008 holds
that role. It stays because with FX-001 it is one of only two entries in AFDB's all-2dp
regime, and its `pDockQ` 0.687671 -> 0.68 is the single measurement proving AFDB truncates.

**Two items deliberately left open, both banner-flagged:** the bootstrap still pins
`REPO_BRANCH = 'homodimer-notebook-rework'` behind `TODO(merge)`, and Colab is unverified.
T041 records a recommendation of beta conditional on both, and says the decision is the user's.

## Milestones and execution order

46 tasks in 8 milestones. Execution model, agreed with the user 2026-09-07: each task goes to
a **fresh agent**, strictly **sequential**, and each result is verified against source and by
execution before the next is dispatched. At every milestone boundary, work **stops** for user
review.

| # | Milestone | Tasks | n | Status |
|---|-----------|-------|---|--------|
| M1 | Ground truth | R001, R002, R002b, R009 | 4 | R002b added 2026-09-08 |
| M2 | Module extracted; notebook shrunk 54% | R010, R010b, R011-R016 | 8 | **COMPLETE** 2026-09-08 |
| M3 | Heterodimer support | R020-R025 | 6 | **COMPLETE** 2026-09-08 |
| M4 | Scoring correctness + Section 4 | R003-R008, R060-R062 | 9 | **COMPLETE** 2026-09-08 |
| M5 | PAE visuals | R050-R052 | 3 | **COMPLETE** 2026-09-08 |
| M6 | 3D views | R070-R075, R030 | 7 | **COMPLETE** 2026-09-08 |
| M7 | Summary + prose | R080-R082, R002b, R004, R040 | 6 | **COMPLETE** 2026-09-08 |
| M8 | Validation, docs, cleanup | R090-R095 | 6 | **COMPLETE** 2026-09-08 |

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
