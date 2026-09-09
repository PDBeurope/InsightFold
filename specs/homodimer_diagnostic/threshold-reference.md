# Canonical Threshold Reference — 7 Complex Confidence Values

**Purpose.** One authoritative threshold set per confidence value, with the provenance of every
number recorded, plus the one *joint* rule that AlphaFold DB actually uses to decide what gets
released. This file supersedes the four conflicting sources previously in the repo
(`CLAUDE.md` `THRESHOLDS`, notebook cell 25 markdown, notebook cell 26 `thresholds`, notebook
cell 28). It is transcribed into `src/insightfold/complex_interface_utils.py`; edit this file
first and re-transcribe, never the other way round.

**Revised 2026-09-08 (R009)** against the AFDB preprint, which was not available when the first
version of this document was written. Six things changed; see "Revision history" at the end.

**Definitions.** Every value below is the quantity defined in
`specs/homodimer_diagnostic/formula-reference.md`, i.e. the `max` row of `ipsae.py` v4. In
particular `ipTM` here means **`ipTM_d0chn`**, a PAE-derived reimplementation, *not* AlphaFold's
own `ipTM` read from the model summary. PAE cutoff for the three ipSAE variants is **10.0**;
LIS cutoff is hardcoded 12; contact cutoff is CB–CB ≤ 8 Å.

The AFDB preprint's `ipSAEmax` and `pDockQ2max` are **`ipsae_d0res`** and **`pdockq2`** as this
project computes them: the per-residue-`d0` ipSAE and pDockQ2, each maxed over the two chain
directions. The paper states the aggregation explicitly (p. 19: "we computed scores for both
chain orderings and defined minimum and maximum such as ipSAEmin = min(ipSAE_A→B, ipSAE_B→A)")
and justifies choosing `max` over `min` on benchmark performance (p. 5, Supplementary Fig. 5,
p. 26) — **independent confirmation of R007**, which fixes this project's pDockQ2 to report the
max over directions rather than A→B alone.

**Sources.** Page numbers refer to the PDFs in `specs/homodimer_diagnostic/references/`
(gitignored). Preprint page numbers are PDF pages.

| Short name | Full reference |
|---|---|
| **AFDB 2026** | **Han, Tsenkov, Venanzi et al., *AlphaFold Database expands to proteome-scale quaternary structures*, bioRxiv 10.64898/2026.03.27.714458v2 (posted 3 Jul 2026), 33 pp.** |
| Dunbrack 2025 | Dunbrack, *ipSAE*, bioRxiv 2025.02.10.637595v2 (posted 15 Dec 2025), 24 pp. |
| Kim 2024 | Kim et al., *LIS*, bioRxiv 2024.02.19.580970v1 (posted 21 Feb 2024), 35 pp. |
| Zhu 2023 | Zhu, Shenoy, Kundrotas, Elofsson, *Bioinformatics* 39(7):btad424, 7 pp. |
| Bryant 2022 | Bryant, Pozzati, Elofsson, *Nat Commun* 13:1265, 11 pp. |

Background reading, not a technical citation: EBI news item, *Millions of protein complexes
added to AlphaFold Database shed light on how proteins interact*, 16 Mar 2026,
<https://www.ebi.ac.uk/about/news/technology-and-innovation/first-complexes-alphafold-database/>.
It is the accessible companion to AFDB 2026 and is worth linking for notebook readers, but it
states no thresholds and its 1.7M figure covers high-confidence **homodimers only**, so it must
never be used as the source for a number. Recorded in the module as `AFDB_NEWS_URL`.

**Provenance labels.**

- `PUBLISHED` — the paper explicitly states or applies this cutoff as an operating criterion.
- `DERIVED` — obtained from a paper's own quantitative statements by a documented inferential
  step, but not stated by it as a recommendation. The step is written out.
- `HEURISTIC` — no literature basis and no large-scale calibration. A judgement call.

The `AFDB-EMPIRICAL` label used by the first version of this document is **retired**. It existed
solely because the ipSAE 0.6 cutoff was then an unpublished AlphaFold DB release filter. The
filter is now published (AFDB 2026), so the label had exactly zero remaining users and was
removed from the module rather than left dead.

---

## 0. The decision rule: AFDB's joint high-confidence criterion

**This, not seven independent traffic lights, is the classifier.**

> `ipSAE_d0res >= 0.6` **AND** `pDockQ2 >= 0.23` → **PASS** / **FAIL**

It answers one question: *would this model qualify for AlphaFold DB high-confidence release?*
Implemented as `afdb_high_confidence(ipsae_d0res, pdockq2)`, which returns the boolean, the two
per-side booleans, the names of the failing side(s), and a display-ready reason.

**Source.** AFDB 2026, p. 5:

> "We adopted a combined high-confidence criterion requiring both community-established cutoffs
> of ipSAEmax ≥ 0.6[32,34,35] and pDockQ2max ≥ 0.23, corresponding to the DockQ 'acceptable'
> quality boundary[36]. This joint cutoff yielded a precision of 0.924 (False Positive Rate
> (FPR) = 0.043) for homodimers and 0.958 (FPR = 0.004) for heterodimers (Fig. 2a,b), supporting
> its use as a quality filter that prioritises precision over recall given the scale of the
> release."

Restated as the release criterion on p. 12 ("selected using the criteria ipSAEmax ≥ 0.6,
pDockQ2max ≥ 0.23") and used throughout the paper's own analyses (pp. 4, 6, 10, 26, 27, 29).

**Benchmark it was validated on** (AFDB 2026, p. 5, Fig. 2a,b legend p. 4): structurally
non-redundant experimental homodimers (n = 230) and heterodimers (n = 94) released *after* the
AlphaFold-Multimer training cutoff with no detectable dimer-level structural similarity to the
training set, against 117 monomeric and 250 non-interacting chain-pair controls.

| | Precision | FPR |
|---|---|---|
| Homodimer benchmark | 0.924 | 0.043 |
| Heterodimer benchmark | 0.958 | 0.004 |

**Why the conjunction and not ipSAE alone.** ipSAEmax "showed the clearest distributional
separation and the most stable Matthews Correlation Coefficient (MCC) plateau (Supplementary
Figs. 1,2), as well as the highest average precision (Supplementary Fig. 3)"; adding pDockQ2max
"further reduced clash-prone predictions not captured by ipSAEmax alone (Supplementary Fig. 4)"
(p. 5). So pDockQ2 is in the rule as a clash-catcher, not as a second opinion on the same
signal.

### 0.1 The nuance that must not be smoothed over: 0.6 is conservative, not optimal

The paper **adopts** 0.6 as "community-established" by citation (refs 32 = Dunbrack 2025, 34 =
Overath et al. 2025, 35 = the Adaptyv Nipah design-competition blog post). It does **not** derive
it. What it derives is the *validation* — the precision/FPR figures above.

Its own MCC-optimal cutoffs for `ipSAEmax`, selected by `argmax(MCC)` on the post-training
benchmarks, are far lower (AFDB 2026, Supplementary Figs. 1 and 2, p. 24):

| Benchmark | MCC-optimal `ipSAEmax` | MCC-optimal `pDockQ2max` |
|---|---|---|
| Homodimer (230 pos / 117 neg) | **0.104** | 0.013 |
| Heterodimer (94 pos / 250 neg) | **0.520** | 0.130 |

The paper's own words for why it uses 0.6 anyway: the filter "prioritises precision over recall
given the scale of the release" (p. 5), and "This joint cutoff for high-confidence predictions
prioritises confidence for AFDB integration; below-threshold predictions remain available for
download and may still contain true interactions" (p. 12).

**Consequence for a notebook user, which R080 must state:** a FAIL means *"not selected for
AFDB high-confidence release"*. It does **not** mean *"not an interaction"*. At release scale a
9.1% pass rate on homodimers is the point of the filter, not a measure of how many predictions
are wrong. Below-threshold dimers are published with their interface scores at
`ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/` (p. 12).

### 0.2 Release scale

All figures AFDB 2026; use these and no others.

| Figure | Value | Page |
|---|---|---|
| Candidate complexes compiled | over **31 million** | p. 11 |
| Homodimers analysed | **19,148,379** | pp. 3, 29 (Supp. Table 1) |
| Heterodimers analysed | **7,561,477** | pp. 3, 29 |
| Homodimers available for bulk download | 21,456,712 | pp. 3, 29 |
| High-confidence homodimers | **1,735,475 (9.1%)** | pp. 3, 5, 12 |
| High-confidence heterodimers | **79,392 (1.0%)** | pp. 3, 5, 12 |
| High-confidence assemblies in total | **1,814,832**, quoted as **1.81 million** | pp. 10, 11 |

Band populations (AFDB 2026, p. 6): very high-confidence 982,188 homo / 23,228 hetero;
confident 440,072 / 31,728; low-confidence 313,215 / 24,436.

**Corrections this replaces.** The first version of this document said the 0.6 cutoff "was
derived data-drivenly from the ~31 million predicted dimers in AlphaFold DB". Two errors:
(a) 0.6 was **not** derived from the release set — the paper adopts it by citation and validates
it on a 230/117 and 94/250 benchmark; (b) ~31 million is the number of **candidate complexes
compiled** (p. 11), not the number of dimers analysed, which is 19,148,379 homo + 7,561,477
hetero. Any prose repeating the old claim is wrong and must be fixed.

---

## Decision table

| Value | Green | Amber | Green provenance | Amber provenance | Bands |
|---|---|---|---|---|---|
| `ipSAE_d0res` | **0.70** | **0.60** | PUBLISHED (AFDB 2026, p. 12) | PUBLISHED (AFDB 2026, p. 12) | **4, published** |
| `ipSAE_d0chn` | **0.70** | **0.60** | DERIVED (transferred; §5) | DERIVED (transferred; §5) | 3 |
| `ipSAE_d0dom` | **0.70** | **0.60** | DERIVED (transferred; §5) | DERIVED (transferred; §5) | 3 |
| `ipTM_d0chn` | **0.70** | **0.30** | DERIVED (Dunbrack 2025, p. 14) | DERIVED (Dunbrack 2025, p. 14) | 3 |
| `pDockQ` | **0.23** | **0.12** | DERIVED (Bryant 2022, pp. 1, 3, 10) | DERIVED (Bryant 2022, pp. 4, 7) | 3 |
| `pDockQ2` | **0.23** | **0.10** | PUBLISHED (Zhu 2023, p. 6; **and** AFDB 2026, p. 5) | HEURISTIC | 3 |
| `LIS` | **0.21** | **0.10** | PUBLISHED (Kim 2024, p. 30) | HEURISTIC | 3 |

Semantics: `value >= green` → green; `green > value >= amber` → amber; `value < amber` → red.
Comparisons are `>=` on both edges. The amber threshold is a stored number, **never** recomputed
as `green / 2` (this is what cell 28 currently does for four of the five values it prints, which
is why the summary table and the prose can disagree on the same run).

### Two band schemes, one colour vocabulary

AFDB publishes **four named bands** for ipSAE. The other six values have only the three-colour
scheme this project invented. Something had to give, and the resolution is:

1. **The colour vocabulary is invariant.** Every value, under either scheme, resolves to exactly
   one of green / amber / red, and `traffic_light()` always returns a 2-tuple. A caller that
   paints a badge never needs to know which scheme a score uses. Nothing about the return type
   or arity varies per score.
2. **The label vocabulary is per-scheme, and is the band's own name.** Six values label their
   bands HIGH / MODERATE / LOW. `ipSAE_d0res` uses AFDB's published wording, because inventing a
   synonym for a published label ("HIGH" for "very high-confidence") would be a worse outcome
   than carrying two label sets.
3. **The difference is explicit, never inferred.** `Threshold.scheme` is `'published_bands'` or
   `'traffic3'`, `Threshold.n_bands` is 4 or 3, and `Threshold.confidence_band(value)` returns
   the whole band record — key, label, colour, lower edge, provenance and page citation.
   `Threshold.__post_init__` refuses to construct a threshold whose named ladder disagrees with
   its own green/amber edges, so the two representations cannot drift.

**The colour mapping for `ipSAE_d0res`**, and why it is the honest one:

| AFDB band | Range | Colour | Reasoning |
|---|---|---|---|
| very high-confidence | `>= 0.80` | green | AFDB's top band |
| confident | `0.70 – 0.80` | green | AFDB calls it confident; nothing to caution about |
| low-confidence | `0.60 – 0.70` | **amber** | AFDB's own caution word. Passes the release filter, but AFDB itself flags it |
| below AFDB threshold | `< 0.60` | **red** | Does not qualify for high-confidence release |

The red edge is therefore *exactly* the ipSAE side of the joint criterion, which is the property
that makes the traffic light and the PASS/FAIL badge tell a consistent story instead of two
different ones.

**Rejected alternative:** collapsing the four bands to three by merging "low-confidence" into
red. That would colour a model AFDB actually released as if it had been rejected.

---

## 1. `ipSAE_d0res`

**Four published bands (AFDB 2026). Green 0.70 — PUBLISHED. Amber 0.60 — PUBLISHED.**

**The bands.** AFDB 2026, p. 12 (Availability):

> "1,735,475 high-confidence homodimer predictions and 79,392 high-confidence heterodimers,
> selected using the criteria ipSAEmax ≥ 0.6, pDockQ2max ≥ 0.23, are available as individual
> entry pages through AlphaFold Database at alphafold.ebi.ac.uk. To support AFDB users in
> interpreting the structure models, surfaced entries are further categorised in 'very
> high-confidence' (ipSAEmax ≥ 0.8), 'confident' (0.7 ≤ ipSAEmax < 0.8), and 'low-confidence'
> (0.6 ≤ ipSAEmax < 0.7)."

Restated with band populations on p. 6:

> "For user interpretation, AFDB entries are additionally labelled by ipSAEmax as very
> high-confidence (≥ 0.8, 982,188 homo- and 23,228 heterodimers), confident (0.7–0.8, 440,072
> homo- and 31,728 heterodimers), or low-confidence (0.6–0.7, 313,215 homo- and 24,436
> heterodimers)."

| Band | Range | Colour | Provenance |
|---|---|---|---|
| very high-confidence | `>= 0.80` | green | PUBLISHED — AFDB 2026, p. 12 |
| confident | `0.70 <= x < 0.80` | green | PUBLISHED — AFDB 2026, p. 12 |
| low-confidence | `0.60 <= x < 0.70` | amber | PUBLISHED — AFDB 2026, p. 12 |
| below AFDB threshold | `< 0.60` | red | PUBLISHED — AFDB 2026, pp. 5, 12 |

These are AFDB's own user-facing labels for AFDB's own entries. This project displays them
verbatim; it does not paraphrase them.

**The metric is Dunbrack's; the thresholds are AFDB's.** Both must be credited, and the
distinction must not be blurred in either direction — see the Citation policy section below.
The AFDB paper itself cites Dunbrack as ref. 32 when adopting the cutoff.

**Dunbrack 2025 states no ipSAE cutoff of any kind**, and that finding stands. It was checked
across the whole 24-page PDF: the Results, "Benchmark of recent PDB entries" (p. 14), the
actifpTM comparison (pp. 16–17), the Genz benchmark (p. 18) and the Discussion (pp. 19–20)
contain no recommended numeric threshold, no ROC-optimal cutoff, and no sensitivity/specificity
table. Full-text search returned zero matches for `thresh*`, `cutoff of 0.`, `ipSAE [><] 0.`,
`above 0.x`, `below 0.x`. **This is a statement about thresholds only.** It is not a criticism
of the paper, whose framing is deliberately comparative — ipSAE separates true from false dimers
*better than* ipTM — and it is not a reason to credit the metric to anyone else.

**Where 0.6 actually comes from, precisely.** AFDB 2026 adopts it as a *community-established*
cutoff by citation to refs 32 (Dunbrack 2025), 34 (Overath et al., *Predicting experimental
success in de novo binder design*, bioRxiv 2025.08.14.670059) and 35 (the Adaptyv Bio Nipah
protein design competition blog post), p. 5. It then validates the joint rule to precision 0.924
/ 0.958. So the number's status is: **published as an operating criterion by AFDB, adopted from
community practice, validated but not optimised.** It is *not* a number the AFDB paper derived,
and it is *not* Dunbrack's number. §0.1 records the MCC-optimal alternatives (0.104 homodimer,
0.520 heterodimer, p. 24), which R080 must surface so a user is not told 0.6 is an optimum.

**Scope caveat — narrower than it first appears.** The joint criterion is a *dataset inclusion*
rule: it decided which models were released, not which released models are good. So for an
**AFDB accession** the value has already passed 0.6, red will essentially never appear, and the
informative distinction is which of the three green/amber bands the entry falls in. For a
**user-uploaded model via `USE_LOCAL_FILE`** no such filtering has occurred and the full range
including red is reachable and meaningful. Both facts belong wherever the badge is presented;
the limitation applies to the AFDB path only, not to the notebook as a whole.

**Independent supporting evidence that the 0.6–0.8 region separates true from false**
(Dunbrack 2025, Table 1, p. 16). At a PAE cutoff of 10 Å, three complexes AlphaFold placed
correctly scored 0.684 / 0.551 / 0.733, while three non-interacting pairs scored 0.019 / 0.012 /
0.000:

| | 1ycr | 2a25 | 3zgc | 4h3b (wrong peptide) | RAF1kd/LysC (non-interacting) | RAF1/RIPK1 (non-interacting) |
|---|---|---|---|---|---|---|
| ipSAE (10 Å) | 0.684 | 0.551 | 0.733 | 0.019 | 0.012 | 0.000 |

n = 6, so this is an illustration, not a calibration; AFDB's 230/117 and 94/250 post-training
benchmarks (p. 5) are the real calibration. But the gap between true and false is wide and
empty, and the 0.6 boundary falls inside it. The paper's supporting prose (p. 11):

> "the overall value of ipSAE is now 0.044, indicating that the proteins are not likely to
> interact"

**The invented amber of 0.30 is gone.** It had no source and was described in the first version
of this document as "a guess". AFDB's published low-confidence band (0.6–0.7) replaces it, so
every edge on this value is now published.

**What the bands mean.** ≥ 0.8: AFDB's top user-facing category. 0.7–0.8: AFDB "confident".
0.6–0.7: passes the release filter but AFDB itself labels it low-confidence. < 0.6: would not be
surfaced as a high-confidence AFDB entry — it would be published on the FTP site with its scores
instead (p. 12), and it may still be a real interaction (§0.1).

---

## 2. `ipSAE_d0chn`

**Green 0.70 — DERIVED. Amber 0.60 — DERIVED.** The numbers are transferred from `ipSAE_d0res`;
see §5 for the justification and for the guaranteed inequality `ipSAE_d0chn >= ipSAE_d0res`.

No publication states a cutoff for this variant. AFDB 2026 calibrated and published on
`ipSAEmax`, which is the `d0res` score (formula-reference.md §1). Dunbrack 2025 p. 22 describes
`ipSAE_d0chn` only as a diagnostic emitted alongside the primary score; it carries no separate
calibration. Transferring AFDB's numbers is therefore one inferential step — hence `DERIVED`
rather than `PUBLISHED`.

**AFDB's band *names* are deliberately not reused here.** Only the numbers transfer. Calling a
`d0chn` value of 0.85 "very high-confidence" would put AFDB's published label on a quantity AFDB
did not measure, on the permissive side of a provable inequality. This variant keeps the
three-colour HIGH / MODERATE / LOW scheme.

---

## 3. `ipSAE_d0dom`

**Green 0.70 — DERIVED. Amber 0.60 — DERIVED.** Same transfer as §2; see §5.

No publication states a cutoff for this variant. Note the transfer is imperfect in the same way:
AFDB calibrated on `ipSAEmax` (`d0res`), not on this variant, and by §5 this variant is provably
`>= d0res`. Three-colour scheme, no published band labels.

---

## 4. `ipTM_d0chn`

**Green 0.70 — DERIVED. Amber 0.30 — DERIVED.**

**The evidence.** Dunbrack 2025, "Benchmark of recent PDB entries", p. 14, describing the two top
panels of Figure 7 (left: AlphaFold's own ipTM from the pTM matrix; right: ipTM computed from the
PAE matrix instead):

> "If we use the PAE values in the ipTM expression (but no PAE cutoff), instead of the AF2 pTM
> matrix, we get quite similar distributions (top right panel). **In both panels, there is overlap
> in the density between values of ipTM or ipSAE from 0.3 to 0.7 for the true dimers** (full-length
> Uniprot sequences, blue curves and data points) **and false dimers** (full-length Uniprot
> sequences, magenta curves and data points)."

The right-hand panel — "ipTM calculated from the PAE matrix instead of the pTM matrix" (Figure 7
legend, p. 16) — **is exactly this project's `ipTM_d0chn`**: same d0 from the summed chain
lengths, same absence of a PAE cutoff, PAE substituted for the pTM probability sum. So the
0.3–0.7 band is measured on the same quantity we compute, on a benchmark of 40 true and 70 false
heterodimers. Green = top of the overlap band, amber = bottom of it.

**Do AlphaFold's published ipTM thresholds transfer?** Partially, and this needs stating in the
notebook.

- *In favour:* Dunbrack measured the two side by side and calls the distributions "quite similar"
  (p. 14). `ipTM_d0chn` also inherits ipTM's defining failure mode — it averages over **all**
  partner residues with no PAE cutoff, so disorder and non-interacting accessory domains dilute
  it exactly as they dilute AlphaFold's ipTM. Dunbrack pp. 7–8 is an extended analysis of this,
  including the worked KRAS/RAF1-RBD case where adding 120 disordered residues to both chains
  drops ipTM from 0.85 to 0.56. A metric with the same artefact should be read on the same scale.
- *Against:* they are not the same number. AlphaFold's ipTM sums over the pTM probability
  distribution in 64 distance bins (Dunbrack p. 5, Eq. 5); `ipTM_d0chn` substitutes the single
  expected PAE value into the same expression. Agreement is empirical, not algebraic, and
  Dunbrack quantifies it only as "quite similar distributions". Treat `ipTM_d0chn` as
  *ipTM-like*, and never quote it to a user as "the AlphaFold ipTM".

**Independent cross-check.** Kim 2024, Figure 4K legend, p. 30, gives a ROC/Youden-optimal cutoff
of **0.38** for AlphaFold's own best ipTM ("best ipTM scores falling below the optimal threshold
(0.38)"). That sits inside Dunbrack's 0.3–0.7 grey zone and just above the amber floor, which is
consistent: a Youden-optimal single cutoff should land inside a two-sided ambiguity band, nearer
its lower edge when negatives outnumber positives.

**AFDB 2026 does not change this value.** Its benchmarks report MCC-optimal cutoffs for ipTM of
**0.391** (homodimer) and **0.471** (heterodimer) (Supplementary Figs. 1–2, p. 24), but ipTM is
not part of its release criterion and those numbers are benchmark diagnostics, not
recommendations. They also disagree with each other, and both refer to AlphaFold's own ipTM
rather than to `ipTM_d0chn`. Recorded, not adopted.

**Rejected alternatives.** `CLAUDE.md`'s 0.70/0.50 has the right green but an amber that is inside
the measured overlap zone, which mislabels genuinely ambiguous models as MODERATE-only. Notebook
cells 25/26/28 use 0.6/0.3; 0.6 is inside the overlap band and would call ambiguous models HIGH.

**What the threshold means.** Green: above the range where Dunbrack's true and false dimer
densities overlap. Red: below that range. Amber: inside it — the score genuinely cannot separate
a true dimer from a decoy here, which is the honest reading.

---

## 5. The three ipSAE variants — can one threshold set serve all three?

**Decision: yes, one set (0.70 / 0.60) for all three, with a documented caveat.** AFDB 2026's
published edges are calibrated on `d0res` only (§1); §2 and §3 transfer the *numbers* to the
other two variants but not AFDB's band *labels*.

**The variants are not on the same scale, and their ordering is a theorem, not an observation.**
All three average the same `ptm(PAE, d0)` over the same set of sub-cutoff pairs; they differ only
in the `L` fed to `d0_func`. `d0_func` is monotone non-decreasing in `L`, and
`ptm(x, d0) = 1/(1 + (x/d0)^2)` is increasing in `d0`. Therefore, per residue `i`:

- `L_chn = nA + nB` (largest possible),
- `L_dom = n0dom` = rows-with-any + cols-with-any ≥ `n0res_i` for every `i`,
- `L_res = n0res_i` (smallest),

giving `ipSAE_d0chn >= ipSAE_d0dom >= ipSAE_d0res` for **every** model, and hence for the
`max`-row value too. (The single exception is the `L == 27` scalar/array `d0` discrepancy
documented in `formula-reference.md` §0, which is far too small to reverse an ordering.)

**Measured offsets.** Computed directly from the AFDB PAE endpoints for the five accessions named
in the spec pack, using the `formula-reference.md` formulas at PAE cutoff 10:

| Accession | nA / nB | `ipTM_d0chn` | `d0chn` | `d0dom` | `d0res` | `d0chn − d0res` |
|---|---|---|---|---|---|---|
| AF-0000000065889468 | 172 / 172 | 0.953 | 0.953 | 0.953 | 0.914 | +0.039 |
| AF-0000000066503175 | 123 / 123 | 0.813 | 0.866 | 0.855 | 0.770 | +0.096 |
| AF-0000000211034637 | 181 / 101 | 0.773 | 0.806 | 0.789 | 0.706 | +0.100 |
| AF-0000000211157965 | 699 / 450 | 0.577 | 0.895 | 0.869 | 0.712 | +0.183 |
| AF-0000000211619209 | 296 / 88 | 0.818 | 0.831 | 0.774 | 0.637 | +0.194 |

(Columns 4–6 are `ipSAE_d0chn` / `ipSAE_d0dom` / `ipSAE_d0res`. These are AFDB-released models,
so all five have already passed the joint criterion of §0 — `ipSAE_d0res >= 0.6` **and**
`pDockQ2 >= 0.23`. This is a high-confidence sample and says nothing about the low end. Note that
under the new bands the last row, `d0res = 0.637`, is AFDB "low-confidence" / amber, not green.)

The offset is real and grows as the interface becomes a smaller fraction of the total length
(0.04 for a compact 172+172 homodimer, 0.19 for a 296+88 pair). Note also row 4, where
`ipTM_d0chn = 0.577` but `ipSAE_d0chn = 0.895` — length dilution of ipTM in a 1149-residue
complex, exactly the artefact Dunbrack designed ipSAE to remove.

**Why one threshold set anyway.**

1. **There is a literature cutoff for exactly one of the three, and it is `d0res`.** AFDB 2026
   publishes bands for `ipSAEmax`, which is the `d0res` score. Assigning *different* numbers to
   `d0chn` and `d0dom` would mean inventing two numbers to sit beside one published one, and
   would manufacture false precision — two decimals of apparent calibration where there is none.
   Transferring the published edges unchanged is the smaller claim.
2. **The variants exist to be compared with each other.** Their whole diagnostic value is that
   `d0chn` high + `d0res` low means "the interface is small relative to the chains". Different
   thresholds per variant would flatten exactly that signal into three independent-looking
   traffic lights and defeat the notebook's stated purpose.
3. **The offset is smaller than the uncertainty in the anchor.** The 0.60 release edge is an
   AFDB operating choice adopted from community practice, deliberately conservative rather than
   optimal (§0.1); its own MCC-optimal counterpart is 0.104 on the homodimer benchmark. A
   0.04–0.19 systematic offset does not justify fabricating per-variant corrections around a
   number whose own placement is that loose.

**Mandatory caveat to surface in the notebook.** Because `ipSAE_d0chn >= ipSAE_d0dom >=
ipSAE_d0res` always holds, a shared threshold is *strictly most permissive for `d0chn` and
strictly least permissive for `d0res`*. `d0chn` will turn green first and red last. Users must
not read three green lights as three independent confirmations — they are one measurement seen
through three progressively more forgiving normalisations. **`ipSAE_d0res` is the one to quote**,
and it is the only one AFDB's criterion and bands are defined on.

**Rejected alternative.** Notebook cell 26's per-variant set (`d0res` 0.6/0.3, `d0chn` 0.5/0.25,
`d0dom` 0.55/0.28) moves the thresholds in the **wrong direction**: it makes `d0chn` — the
provably most permissive variant — the *easiest* to turn green, compounding the bias instead of
compensating for it. If per-variant thresholds were ever adopted, `d0chn` would need a *higher*
green than `d0res`, not a lower one. No source in the repo justifies these three numbers.

---

## 6. `pDockQ`

**Green 0.23 — DERIVED. Amber 0.12 — DERIVED.**

**0.23 traces to Bryant 2022 and it is the DockQ acceptability cutoff, not a pDockQ
recommendation.** Bryant never writes "use pDockQ ≥ 0.23". What the paper states is:

- Abstract, p. 1: "generate models with **acceptable quality (DockQ ≥ 0.23)** for 63% of the
  dimers."
- p. 3: "We measure the separation between **correct (DockQ ≥ 0.23) and incorrect models**
  provided by several metrics using a receiver operating characteristic (ROC) curve."
- Methods, p. 10, Eq. 6: "`SR = Fraction of predicted models with DockQ ≥ 0.23`".

0.23 is the CAPRI "acceptable" boundary imported from Basu & Wallner's DockQ. The single
inferential step that turns it into a pDockQ threshold is that **pDockQ is a sigmoid fitted
directly to DockQ as the target**, so a pDockQ of 0.23 *is* a predicted DockQ of 0.23. Methods,
p. 10, Eq. 7–8, with `L = 0.724`, `x0 = 152.611`, `k = 0.052`, `b = 0.018` — the constants in
`ipsae.py` and in the notebook. Hence DERIVED rather than PUBLISHED: one step, and a short one,
but a step.

Zhu 2023 independently confirms the reading of 0.23 (p. 2, §2.4.1): "According to the study
(Basu and Wallner 2016), **if the DockQ score > 0.23, then the quality of the model is acceptable
by the CAPRI criteria**."

**Amber 0.12 — derived from the score's own stated error.** Bryant reports the calibration error
of pDockQ twice:

- p. 4: "This enables the prediction of the DockQ scores (pDockQ) in a continuous manner **with an
  overall average error of 0.11** on the test set."
- p. 7 (Discussion): "we can predict the DockQ score (pDockQ) **with an average error of 0.1**."

A model at `pDockQ = 0.23 − 0.11 = 0.12` is one average prediction error below the acceptability
line, i.e. it could still be an acceptable model that pDockQ has under-scored. Below 0.12 the
model is more than one stated error away from acceptable. This is a defensible amber floor from
the paper's own numbers rather than an invented one.

**Range note.** The sigmoid is bounded on `[b, L + b] = [0.018, 0.742]`, so 0.23 is not far above
the floor; a value of exactly 0.018 means *zero contacts*, not "slightly bad" (see
`formula-reference.md` §8 — `ipsae.py` returns 0.0 there instead, and that discrepancy must be
resolved separately from this document).

**Known weakness, worth surfacing.** Zhu 2023 p. 5 found pDockQ over-optimistic outside Bryant's
heterodimer training regime: "**More than 10% of the chains in all these sets have pDockQ > 0.5
and DockQ_i < 0.23.** pDockQ does not utilize the predicted average errors (PAEs) but only
considers the size of the interface and the predicted quality (pLDDT) of residues in the
interface. Therefore, it does not work if a method generates models with large, highly confident
incorrect interfaces." Since this project's models are AlphaFold-Multimer homodimers — precisely
the regime Zhu tested — a green pDockQ should never be read as confirmation on its own.

**All four current sources already agree on 0.23**; only the amber differs (`CLAUDE.md` 0.09,
notebook 0.115 = green/2). Neither has a recorded source. 0.12 replaces both.

---

## 7. `pDockQ2`

**Green 0.23 — PUBLISHED, and double-sourced. Amber 0.10 — HEURISTIC.**

**Two independent publications land on 0.23.** Zhu 2023 applies it as an operating cutoff, and
AFDB 2026 arrives at it separately as the second half of its release criterion (p. 5):

> "…and pDockQ2max ≥ 0.23, corresponding to the DockQ 'acceptable' quality boundary[36]."

AFDB's role for pDockQ2 in the joint rule is specifically to catch clash-prone predictions that
ipSAE alone lets through (p. 5, Supplementary Fig. 4), which is a different justification from
Zhu's and therefore genuine corroboration rather than a restatement. Note that AFDB's own
MCC-optimal cutoff for `pDockQ2max` is **0.013** on the homodimer benchmark and **0.130** on the
heterodimer benchmark (Supplementary Figs. 1–2, p. 24), so 0.23 is conservative here too. The
metric is Zhu's in both cases.

**Zhu et al. settle the number: the green is 0.23, not 0.15 and not 0.5.** pDockQ2 is fitted to the same
target as pDockQ — "As in pDockQ, we fit a sigmoid curve (Equation 1) (by the scipy package) **to
the actual DockQ_i values**, yielding the coefficients L = 1.31, x0 = 84.733, k = 0.075, and
b = 0.005" (§2.6, p. 3; constants match `ipsae.py`). So the same CAPRI acceptability boundary
applies. Unlike Bryant with pDockQ, Zhu then **actually applies it as an operating cutoff**
(p. 6):

> "Out of the 29 potentially novel complexes (Supplementary Table S8) that were predicted by
> AlphaFold-Multimer, **9 had pTM > 0.5 and min pDockQ2 > 0.23** (see Supplementary Fig. S10).
> **Two highly confident complexes (pTM > 0.75 and pDockQ2 > 0.23)** were obtained."

This is the paper selecting real biological conclusions using `pDockQ2 > 0.23`. That makes it
PUBLISHED, and — mildly ironically — better sourced than pDockQ's own 0.23.

**Where 0.15 and 0.5 came from: nowhere.** Notebook cells 25/26/28 use 0.15; `CLAUDE.md` uses 0.5.
Neither number appears anywhere in Zhu 2023 in connection with a pDockQ2 cutoff. The only 0.5 in
the paper attached to a pDockQ-family score is the *criticism* of pDockQ quoted in §6 above
("pDockQ > 0.5 and DockQ_i < 0.23") — a statement that high pDockQ can be wrong, which is the
opposite of a recommendation. Both are dropped.

Note also that `CLAUDE.md`'s pDockQ2 has the wrong sigmoid constants entirely
(`0.715 / (1 + exp(-12.3(x - 0.605))) + 0.005`, see `formula-reference.md` §6), so its 0.5/0.23
threshold pair was calibrated — if at all — against a different function. It cannot be carried
over.

**Amber 0.10 — HEURISTIC.** Zhu gives no second tier, no ROC table, and no calibration error for
pDockQ2, so there is nothing to derive from. 0.10 is chosen as a round number below green; the
only support is that 0.1 is the boundary Zhu themselves use when bucketing low-confidence
interfaces in the outlier analysis — "Using **pDockQ2 < 0.1** and DockQ_i > 0.6 resulted in 55
interfaces" (p. 6), and the worst-case hexamer example is reported as "the min pDockQ2 over all
the interfaces for this prediction is 0.011". That is a bucket boundary for error analysis, not
a recommended cutoff, and it must not be cited as one. **This is a judgement call.**

**Range note.** The sigmoid is bounded on `[0.005, 1.315]` — it can exceed 1.0, which no
DockQ-scaled quantity should. Do not present pDockQ2 as a probability or as a percentage.

---

## 8. `LIS`

**Green 0.21 — PUBLISHED. Amber 0.10 — HEURISTIC.**

**Kim et al. state a specific cutoff, and it is 0.21.** Figure 4K legend, p. 30:

> "The selected PPIs exhibit **best LIS values exceeding the optimal threshold (0.21)** and best
> ipTM scores falling below the optimal threshold (0.38)."

The thresholds are ROC-derived. Methods, "Receiver operating characteristic (ROC) analysis",
p. 15:

> "ROC analysis was conducted using the sklearn library in Python, with AUC values computed using
> positive sets (fly/human PRS) and negative control sets (fly/human RRS, PRS+GFP, PRS+Wg).
> **Optimal thresholds were established based on the highest Youden's Index.**"

**"Best" is the right variant for this project.** Kim reports two LIS flavours (Methods, p. 14):
"The 'best' LIS was derived from the **rank 1 model**, while 'average' LIS was computed from ranks
1–5." AFDB serves one model per accession, so the rank-1 ("best") threshold is the transferable
one. The paper quotes only the best-LIS cutoff numerically (0.21); the average-LIS optimal value
is stated to exist (Supplementary Figure 8) but is not printed in the main text or the figure
legends of this PDF. **Do not cite an average-LIS number — it is not in the source.**

**Three things that make 0.21 transfer cleanly to our implementation:**

1. **Same PAE cutoff.** "A cutoff PAE value of 12, determined to provide the highest AUC, was used
   for both average and best LIS as per ROC analysis" (Methods, p. 14; result on p. 5: "a PAE
   cutoff of 12 maximizes the area under the curve (AUC)"). `ipsae.py` hardcodes 12
   (`ipsae_v4.py:712`).
2. **Same directional combination — the mean, not the max.** Kim's own worked example (Figure 2A,
   p. 25) shows an A:B interface score of 0.457 and a B:A score of 0.309 combining as
   `(0.457 + 0.309) / 2 = 0.383 (LIS)`. `ipsae.py:982` takes the mean. The notebook takes the
   mean. **`CLAUDE.md` takes the `max`** (`formula-reference.md` §7), which is always ≥ the mean
   and therefore systematically inflates LIS — so `CLAUDE.md`'s low green of 0.15 was, at best,
   compensating for a bug. Once LIS is computed as the mean, 0.21 is the correct comparison.
3. **Comparable model provenance.** Kim ran ColabFold / AlphaFold-Multimer (LocalColabFold script,
   Methods p. 14). AFDB homodimer entries report `"toolUsed": "ColabFold v1.6.0 /
   AlphaFold-Multimer"`. Same generator family.

**Neither 0.3 nor 0.15 has a source.** Notebook cells 25/26/28 use green 0.3; `CLAUDE.md` uses
green 0.15. The number 0.3 does not appear in Kim 2024 as a LIS threshold anywhere. Both are
dropped in favour of the paper's 0.21.

**Calibration context for reading a green LIS.** Kim reports that at these thresholds
"over 70% of yeast positive PPIs and ELM-annotated PPIs were predicted to be direct interactions.
**Approximately 5% of yeast negative control sets had values that exceeded these thresholds**"
(p. 7). So green LIS ≈ a ~5% false-positive operating point on their reference sets. On the
genome-scale datasets the pass rates were far lower — 19% of FlyBi Y2H pairs, 13% of DPiM IP-MS
pairs, 24% of literature-curated pairs met the best-LIS threshold (pp. 10–11).

**AFDB 2026 does not change this value either.** It benchmarks `LISmax` alongside ipSAE and
pDockQ2 and reports MCC-optimal cutoffs of **0.085** (homodimer) and **0.347** (heterodimer)
(Supplementary Figs. 1–2, p. 24) — a factor of four apart, on a metric it did not adopt into its
release criterion. Kim's own ROC-optimal 0.21 remains the published operating cutoff and stays.

**Amber 0.10 — HEURISTIC.** A Youden-optimal cutoff is by construction a single boundary; Kim
publishes no second tier and no gradation below 0.21. 0.10 is chosen as a round number roughly
half-way to zero. LIS has no sigmoid floor — a pair with no inter-chain PAE below 12 scores
exactly 0.0 — so values in 0.10–0.21 do mean "there is a real sub-cutoff PAE patch, but a weak
or small one". That is the intended reading, but the specific number is **a judgement call.**

**What the threshold means.** Green: at or above Kim's ROC-optimal operating point for calling a
pair a *direct, physical* interactor — a statement about whether the proteins interact, **not**
about whether the docked pose is correct. This is a different question from pDockQ/pDockQ2, which
predict pose accuracy. A model can legitimately be green on LIS and red on pDockQ.

---

## Citation policy — metric vs threshold

Decided with the user 2026-09-08. **Separate the metric from the threshold, and credit both.**
This is the rule for every consumer of this document: the notebook, the module docstrings, the
summary table, and any prose R080 generates.

1. **A metric belongs to whoever published it, whatever threshold is later applied to it.**
   ipSAE is Dunbrack's (bioRxiv 2025.02.10.637595); pDockQ is Bryant's; pDockQ2 is Zhu's; LIS is
   Kim's; the reference implementation is DunbrackLab/IPSAE `ipsae.py` v4. Wherever a score is
   explained, defined or computed, that is the citation. Adopting a different threshold does not
   diminish authorship of the method.
2. **The finding that Dunbrack states no cutoff is about *thresholds only*.** It must never be
   presented as a criticism of the ipSAE paper, nor used as a reason to under-credit it. AFDB
   2026 itself cites Dunbrack (ref. 32) when adopting the 0.6 cutoff.
3. **The ipSAE 0.6 cutoff and the four confidence bands are AlphaFold DB's**, cited to Han,
   Tsenkov, Venanzi et al. 2026. Do not attribute them to Dunbrack — which is the error
   `CLAUDE.md` currently makes by placing "Primary AFDB classifier: `ipSAE_d0res ≥ 0.6`" under a
   heading whose formulas are all attributed to DunbrackLab/IPSAE. (Fixing `CLAUDE.md` is out of
   scope for R009 and belongs to the task that rewrites it.)
4. **pDockQ2's 0.23 is double-sourced and both sources are named** (§7): Zhu 2023 p. 6 applies
   it directly, and AFDB 2026 p. 5 arrives at it independently via the DockQ "acceptable"
   boundary. The metric is Zhu's in both cases.
5. **The EBI news item is cited as background**, not as a technical source:
   <https://www.ebi.ac.uk/about/news/technology-and-innovation/first-complexes-alphafold-database/>
   It predates the v2 preprint and is the accessible entry point to the same work. Present it as
   further reading with the preprint as the technical citation. Its 1.7M figure covers
   high-confidence homodimers only — never use it for the scale numbers (§0.2).

In the module these are carried in `Threshold.source` and `Threshold.note` per score, in
`AFDB_JOINT_CRITERION.source` / `.note` for the joint rule, and in `AFDB_NEWS_URL`.

---

## Values with no literature basis

Everything in this section is a judgement call. Nothing here is supported by a publication, and
none of it should be presented to a user as a literature-derived number.

| Value | Threshold | Number | Why it is a guess |
|---|---|---|---|
| `pDockQ2` | amber | 0.10 | No source. Zhu 2023 gives no second tier, no ROC table, no calibration error, and AFDB 2026 publishes only the 0.23 edge. 0.1 appears in Zhu only as an error-analysis bucket boundary (p. 6). |
| `LIS` | amber | 0.10 | No source. Kim 2024 publishes a single Youden-optimal cutoff with no gradation below it. |

**Two judgement calls, both amber, down from six.** R009 retired the three invented ipSAE ambers
(0.30) by adopting AFDB's published bands, and moved the three ipSAE greens off the retired
`AFDB-EMPIRICAL` label. Every remaining green has a publication behind it. The two surviving
ambers are where the guesswork lives, because neither Zhu nor Kim defines a middle band — both
answer a binary question and say nothing about an intermediate tier. That amber tier is this
notebook's own construct and should be presented as such.

### Transferred, not published

| Value | Threshold | Number | Basis |
|---|---|---|---|
| `ipSAE_d0chn` | green / amber | 0.70 / 0.60 | `DERIVED`. AFDB 2026's published edges, transferred from `d0res` by §5. AFDB did not calibrate this variant, and it is provably the more permissive one. |
| `ipSAE_d0dom` | green / amber | 0.70 / 0.60 | `DERIVED`. Same transfer, same caveat. |

**Papers that give no usable cutoff at all, stated plainly:**

- **Dunbrack 2025 gives no ipSAE threshold.** Not for `d0res`, not for `d0chn`, not for `d0dom`,
  not at any PAE cutoff. Its only numeric guidance is about the *PAE* cutoff, not the score:
  "Cutoffs of 10 or 15 Å may be most suitable" (p. 14). This is a fact about the paper's scope,
  not a shortcoming — see the Citation policy above.
- **Bryant 2022 never states a pDockQ cutoff.** 0.23 is DockQ's, imported one inferential step.
- **Kim 2024 does not print the average-LIS optimal threshold** in the main text or figure
  legends; only the best-LIS value 0.21 is given. No average-LIS number should be quoted from
  this source.
- **AFDB 2026 publishes cutoffs only for `ipSAEmax` and `pDockQ2max`.** Its Supplementary
  Figs. 1–2 (p. 24) report MCC-optimal cutoffs for ipTM (0.391 homo / 0.471 hetero) and LISmax
  (0.085 / 0.347) as benchmark diagnostics, **not** as recommendations, and they disagree with
  each other across the two benchmarks by more than a factor of four. They are recorded here so
  nobody re-derives them, and they are deliberately **not** adopted as thresholds.

---

## Changes this file makes to the repo's current numbers

| Value | `CLAUDE.md` | notebook cell 26 | pre-R009 | **Canonical** | Change |
|---|---|---|---|---|---|
| `ipSAE_d0res` | 0.6 / 0.4 | 0.6 / 0.3 | 0.60 / 0.30 | **0.70 / 0.60** + 4 published bands | adopts AFDB's published bands (R009) |
| `ipSAE_d0chn` | — | 0.5 / 0.25 | 0.60 / 0.30 | **0.70 / 0.60** | transferred with `d0res`; cell 26 had the bias backwards (§5) |
| `ipSAE_d0dom` | — | 0.55 / 0.28 | 0.60 / 0.30 | **0.70 / 0.60** | transferred; same reason |
| `ipTM_d0chn` | 0.7 / 0.5 | 0.6 / 0.3 | 0.70 / 0.30 | 0.70 / 0.30 | unchanged by R009 |
| `pDockQ` | 0.23 / 0.09 | 0.23 / 0.115 | 0.23 / 0.12 | 0.23 / 0.12 | unchanged by R009 |
| `pDockQ2` | 0.5 / 0.23 | 0.15 / 0.075 | 0.23 / 0.10 | 0.23 / 0.10 | green now double-sourced (Zhu **and** AFDB) |
| `LIS` | 0.15 / 0.09 | 0.3 / 0.15 | 0.21 / 0.10 | 0.21 / 0.10 | unchanged by R009 |

Plus the new first-class concept: **`ipSAE_d0res >= 0.6 AND pDockQ2 >= 0.23` → PASS/FAIL**
(§0), which is not a threshold on any single value and does not belong in the table above.

Additionally, when this table is implemented:

1. **Amber must be read from the table, never recomputed.** Cell 28 currently derives amber as
   `green / 2` for `ipSAE_d0res`, `ipTM`, `pDockQ` and `LIS` while cell 26 reads the dict, so the
   summary table and the prose statement can label the same score differently in a single run.
   Both must consume one shared structure.
2. **`ipTM` should be relabelled `ipTM_d0chn`** everywhere in the notebook. It is not AlphaFold's
   ipTM, AFDB does not expose AlphaFold's ipTM on these endpoints, and thresholds published for
   AlphaFold's ipTM transfer only approximately (§4).
3. **Provenance should be visible in the UI.** Each traffic light should be able to state whether
   its threshold is published, derived, or a judgement call, and `Threshold.confidence_band()`
   returns the page citation for the exact edge that produced the colour.
4. **The PASS/FAIL badge leads; the seven lights explain it.** The joint criterion is the
   headline verdict. The per-score lights are diagnostics that say *why*, and must not be
   presented as seven independent confirmations — see §5 for why three of them are one
   measurement seen three ways.
5. **`AFDB_RELEASE_SCALE` exists so prose never invents a scale figure.** Every number in §0.2 is
   in the module with its page.

---

## Revision history

**2026-09-08 — R009**, against `references/afdb_dimers_2026.pdf` (33 pp.), which was unavailable
when this document was first written. Six changes:

1. **Added §0, the AFDB joint criterion** `ipSAE_d0res >= 0.6 AND pDockQ2 >= 0.23`, as a
   first-class concept rather than a seventh traffic light, with its precision/FPR validation.
   Implemented as `afdb_high_confidence()`.
2. **`ipSAE_d0res` adopts AFDB's four published bands** (p. 12), replacing the invented 0.60/0.30
   green/amber. Colour edges become 0.70 / 0.60.
3. **`ipSAE_d0res` green moves to `PUBLISHED`**, cited to Han/Tsenkov/Venanzi et al. 2026, not to
   Dunbrack. `ipSAE_d0chn` / `d0dom` become `DERIVED` transfers. The **`AFDB-EMPIRICAL` label is
   retired** — it had no remaining users and was removed from the module rather than left dead.
   (It appears in this document only in this history and in the label list's retirement note.)
4. **`pDockQ2` 0.23 recorded as double-sourced** — Zhu 2023 p. 6 directly, AFDB 2026 p. 5
   independently.
5. **Citation policy written down** (metric vs threshold), including the EBI news URL as
   accessible background.
6. **Scale figures corrected** (§0.2). The previous claim that 0.6 was "derived data-drivenly
   from the ~31 million predicted dimers" was wrong on both counts and is corrected in §0.2, and
   the conservative-not-optimal nuance with the MCC-optimal cutoffs (0.104 / 0.520, p. 24) is
   recorded in §0.1 for R080 to surface.
