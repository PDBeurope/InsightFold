# Canonical Threshold Reference — 7 Complex Confidence Values

**Purpose.** One authoritative green/amber threshold per confidence value, with the provenance of
every number recorded. This file supersedes the four conflicting sources currently in the repo
(`CLAUDE.md` `THRESHOLDS`, notebook cell 25 markdown, notebook cell 26 `thresholds`, notebook
cell 28). Nothing in this file has been implemented yet — a later task does that.

**Definitions.** Every value below is the quantity defined in
`specs/homodimer_diagnostic/formula-reference.md`, i.e. the `max` row of `ipsae.py` v4. In
particular `ipTM` here means **`ipTM_d0chn`**, a PAE-derived reimplementation, *not* AlphaFold's
own `ipTM` read from the model summary. PAE cutoff for the three ipSAE variants is **10.0**;
LIS cutoff is hardcoded 12; contact cutoff is CB–CB ≤ 8 Å.

**Sources.** Page numbers refer to the PDFs in `specs/homodimer_diagnostic/references/`
(gitignored):

| Short name | Full reference |
|---|---|
| Dunbrack 2025 | Dunbrack, *ipSAE*, bioRxiv 2025.02.10.637595v2 (posted 15 Dec 2025), 24 pp. |
| Kim 2024 | Kim et al., *LIS*, bioRxiv 2024.02.19.580970v1 (posted 21 Feb 2024), 35 pp. |
| Zhu 2023 | Zhu, Shenoy, Kundrotas, Elofsson, *Bioinformatics* 39(7):btad424, 7 pp. |
| Bryant 2022 | Bryant, Pozzati, Elofsson, *Nat Commun* 13:1265, 11 pp. |

**Provenance labels.**

- `PUBLISHED` — the paper explicitly states or applies this cutoff as an operating criterion.
- `DERIVED` — obtained from the paper's own quantitative statements, but not stated as a
  recommendation. The derivation is written out.
- `AFDB-EMPIRICAL` — not in any publication, but established by AlphaFold DB from a
  large-scale internal analysis. Empirically grounded; simply unpublished.
- `HEURISTIC` — no literature basis and no large-scale calibration. A judgement call.

---

## Decision table

| Value | Green | Amber | Green provenance | Amber provenance |
|---|---|---|---|---|
| `ipSAE_d0res` | **0.60** | **0.30** | AFDB-EMPIRICAL (31M-dimer analysis; see §1) | HEURISTIC |
| `ipSAE_d0chn` | **0.60** | **0.30** | AFDB-EMPIRICAL (same set; see §5 for why one set) | HEURISTIC |
| `ipSAE_d0dom` | **0.60** | **0.30** | AFDB-EMPIRICAL (same set; see §5) | HEURISTIC |
| `ipTM_d0chn` | **0.70** | **0.30** | DERIVED (Dunbrack 2025, p. 14) | DERIVED (Dunbrack 2025, p. 14) |
| `pDockQ` | **0.23** | **0.12** | DERIVED (Bryant 2022, pp. 1, 3, 10) | DERIVED (Bryant 2022, pp. 4, 7) |
| `pDockQ2` | **0.23** | **0.10** | PUBLISHED (Zhu 2023, p. 6) | HEURISTIC |
| `LIS` | **0.21** | **0.10** | PUBLISHED (Kim 2024, p. 30) | HEURISTIC |

Semantics: `value >= green` → green/HIGH; `green > value >= amber` → amber/MODERATE;
`value < amber` → red/LOW. Comparisons are `>=` on both edges. The amber threshold is a stored
number, **never** recomputed as `green / 2` (this is what cell 28 currently does for four of the
five values it prints, which is why the summary table and the prose can disagree on the same run).

---

## 1. `ipSAE_d0res`

**Green 0.60 — AFDB-EMPIRICAL. Amber 0.30 — HEURISTIC.**

**The Dunbrack paper states no ipSAE cutoff of any kind.** This was checked across the whole
24-page PDF: the Results, the "Benchmark of recent PDB entries" section (p. 14), the
actifpTM comparison (pp. 16–17), the Genz benchmark (p. 18) and the Discussion (pp. 19–20)
contain no recommended numeric threshold, no ROC-optimal cutoff, and no sensitivity/specificity
table. The paper's own framing is comparative — ipSAE separates true from false dimers *better
than* ipTM — not classificatory.

**Where 0.60 actually comes from.** It is an **AlphaFold DB production release criterion**, not a
result in the paper. It enters this repo at
`specs/homodimer_notebook_extraction.md:325` and
`specs/homodimer_diagnostic/homodimer_diagnostic_notebook_spec.md:17`:

> "The AlphaFold Database released predicted homodimers filtered by a confidence band:
> `ipSAE >= 0.6`, `backbone clashes <= 10`, and `average pLDDT >= 70`."

`CLAUDE.md:199` restates this as "Primary AFDB classifier". **Flagged:** the phrasing in
`CLAUDE.md` sits directly under a heading whose formulas are all attributed to
DunbrackLab/IPSAE, which invites the reader to assume 0.6 is Dunbrack's number. It is not.
It is AFDB's.

**It is nonetheless empirically grounded, and that must be stated in the notebook**
(user, 2026-09-07): the 0.6 cutoff was **derived data-drivenly from the ~31 million predicted
dimers in AlphaFold DB**. It is unpublished, not unfounded. On sample size alone that is a
far broader calibration than anything in the four reference papers — Kim's ROC-optimal 0.21
comes from a few thousand curated pairs, Dunbrack's Table 1 illustration from six complexes.
The notebook should say where 0.6 comes from rather than presenting it as a bare number, and
should not imply it is Dunbrack's.

**Scope caveat — and it is narrower than it first appears.** The AFDB filter is a *dataset
inclusion* criterion: it decided which models were released, not which released models are
good. So for an **AFDB accession**, the green light is largely re-measuring the filter that
caused the entry to exist, and red will essentially never appear. For a **user-uploaded model
via `USE_LOCAL_FILE`**, no such filtering has occurred and the full range including red is
reachable and meaningful. Both facts should be stated where the traffic light is presented;
the limitation applies to the AFDB path only, not to the notebook as a whole.

**Supporting evidence that 0.60 is at least not absurd** (Dunbrack 2025, Table 1, p. 16). At a
PAE cutoff of 10 Å, three complexes AlphaFold placed correctly scored 0.684 / 0.551 / 0.733,
while three non-interacting pairs scored 0.019 / 0.012 / 0.000:

| | 1ycr | 2a25 | 3zgc | 4h3b (wrong peptide) | RAF1kd/LysC (non-interacting) | RAF1/RIPK1 (non-interacting) |
|---|---|---|---|---|---|---|
| ipSAE (10 Å) | 0.684 | 0.551 | 0.733 | 0.019 | 0.012 | 0.000 |

n = 6, so this is an illustration, not a calibration. But the gap between true and false is wide
and empty, and both 0.60 and 0.30 fall inside it. The paper's supporting prose (p. 11):

> "the overall value of ipSAE is now 0.044, indicating that the proteins are not likely to interact"

**Amber 0.30.** No basis in the paper. Chosen because Dunbrack's only quantitative statement about
an ambiguous band for a PAE-based ipTM-like score is 0.3–0.7 (p. 14, quoted in §4 below), and
because in Table 1 nothing genuine falls below 0.49 and nothing false rises above 0.02, so any
value in roughly 0.1–0.5 is equally defensible. 0.30 is a round number in the middle of that
range. **This is a guess.**

**What the threshold means.** Green: the model would have passed AFDB's release filter, and its
ipSAE is in the range Dunbrack's true complexes occupy. Red: the value is in the range Dunbrack's
non-interacting decoys occupy.

---

## 2. `ipSAE_d0chn`

**Green 0.60 — AFDB-EMPIRICAL (inherited). Amber 0.30 — HEURISTIC.** Same set as
`ipSAE_d0res`; see §5 for the justification and for the guaranteed inequality
`ipSAE_d0chn >= ipSAE_d0res`.

No publication states a cutoff for this variant. Dunbrack 2025 p. 22 describes `ipSAE_d0chn` only
as a diagnostic emitted alongside the primary score; it carries no separate calibration.

---

## 3. `ipSAE_d0dom`

**Green 0.60 — AFDB-EMPIRICAL (inherited). Amber 0.30 — HEURISTIC.** Same set; see §5.

No publication states a cutoff for this variant. Note the inheritance is imperfect: AFDB's
31M-dimer calibration was performed on `ipSAE` as `ipsae.py` reports it (the `d0res` max row),
not on this variant, and by §5 this variant is provably >= `d0res`.

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

**Rejected alternatives.** `CLAUDE.md`'s 0.70/0.50 has the right green but an amber that is inside
the measured overlap zone, which mislabels genuinely ambiguous models as MODERATE-only. Notebook
cells 25/26/28 use 0.6/0.3; 0.6 is inside the overlap band and would call ambiguous models HIGH.

**What the threshold means.** Green: above the range where Dunbrack's true and false dimer
densities overlap. Red: below that range. Amber: inside it — the score genuinely cannot separate
a true dimer from a decoy here, which is the honest reading.

---

## 5. The three ipSAE variants — can one threshold set serve all three?

**Decision: yes, one set (0.60 / 0.30) for all three, with a documented caveat.**

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
so all five have already passed the `ipSAE >= 0.6` filter — this is a high-confidence sample and
says nothing about the low end.)

The offset is real and grows as the interface becomes a smaller fraction of the total length
(0.04 for a compact 172+172 homodimer, 0.19 for a 296+88 pair). Note also row 4, where
`ipTM_d0chn = 0.577` but `ipSAE_d0chn = 0.895` — length dilution of ipTM in a 1149-residue
complex, exactly the artefact Dunbrack designed ipSAE to remove.

**Why one threshold set anyway.**

1. **There is no literature cutoff for *any* of the three.** Assigning three different numbers
   would be inventing a fourth, fifth and sixth number on top of an already-invented one, and
   would manufacture false precision — three decimals of apparent calibration where there is
   none.
2. **The variants exist to be compared with each other.** Their whole diagnostic value is that
   `d0chn` high + `d0res` low means "the interface is small relative to the chains". Different
   thresholds per variant would flatten exactly that signal into three independent-looking
   traffic lights and defeat the notebook's stated purpose.
3. **The offset is smaller than the uncertainty in the anchor.** 0.60 itself is an AFDB
   operational choice with no error bar; a 0.04–0.19 systematic offset does not justify
   fabricating per-variant corrections around it.

**Mandatory caveat to surface in the notebook.** Because `ipSAE_d0chn >= ipSAE_d0dom >=
ipSAE_d0res` always holds, a shared threshold is *strictly most permissive for `d0chn` and
strictly least permissive for `d0res`*. `d0chn` will turn green first and red last. Users must
not read three green lights as three independent confirmations — they are one measurement seen
through three progressively more forgiving normalisations. `ipSAE_d0res` is the one to quote.

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

**Green 0.23 — PUBLISHED. Amber 0.10 — HEURISTIC.**

**Zhu et al. settle this: the green is 0.23, not 0.15 and not 0.5.** pDockQ2 is fitted to the same
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

## Values with no literature basis

Everything in this section is a judgement call. Nothing here is supported by a publication, and
none of it should be presented to a user as a literature-derived number.

| Value | Threshold | Number | Why it is a guess |
|---|---|---|---|
| `ipSAE_d0res` | amber | 0.30 | No source. Middle of the empty gap between true (0.55–0.73) and false (0.00–0.02) in Dunbrack Table 1, n = 6. |
| `ipSAE_d0chn` | amber | 0.30 | No source. Inherited from `d0res` by the §5 argument. |
| `ipSAE_d0dom` | amber | 0.30 | No source. Inherited. |
| `pDockQ2` | amber | 0.10 | No source. Zhu 2023 gives no second tier, no ROC table, no calibration error. 0.1 appears only as an error-analysis bucket boundary (p. 6). |
| `LIS` | amber | 0.10 | No source. Kim 2024 publishes a single Youden-optimal cutoff with no gradation below it. |

**Five judgement calls, all of them amber.** Every green threshold now has either a
publication or a large-scale empirical calibration behind it. The ambers are where the
guesswork lives, because none of the four papers defines a second tier — they all answer a
binary question (acceptable / not, interacting / not) and say nothing about a middle band.
The amber tier is this notebook's own construct and should be presented as such.

### Not a guess, but not published either

| Value | Threshold | Number | Basis |
|---|---|---|---|
| `ipSAE_d0res` | green | 0.60 | `AFDB-EMPIRICAL`. Derived data-drivenly from the ~31M predicted dimers in AlphaFold DB. Unpublished, but a far larger calibration set than any of the four papers. See §1. |
| `ipSAE_d0chn` | green | 0.60 | `AFDB-EMPIRICAL`, inherited via §5. |
| `ipSAE_d0dom` | green | 0.60 | `AFDB-EMPIRICAL`, inherited via §5. |

**Not in either table** (these are sourced from publications): `ipTM_d0chn` 0.70/0.30
(DERIVED, Dunbrack p. 14 overlap band); `pDockQ` 0.23 (DERIVED, CAPRI/DockQ acceptability via
Bryant pp. 1/3/10) and 0.12 (DERIVED, Bryant's stated average error of 0.11, pp. 4/7);
`pDockQ2` 0.23 (PUBLISHED, Zhu p. 6); `LIS` 0.21 (PUBLISHED, Kim p. 30).

**Papers that give no usable cutoff at all, stated plainly:**

- **Dunbrack 2025 gives no ipSAE threshold.** Not for `d0res`, not for `d0chn`, not for `d0dom`,
  not at any PAE cutoff. Its only numeric guidance is about the *PAE* cutoff, not the score:
  "Cutoffs of 10 or 15 Å may be most suitable" (p. 14). Three of the seven values in this table
  therefore have no literature anchor whatsoever.
- **Bryant 2022 never states a pDockQ cutoff.** 0.23 is DockQ's, imported one inferential step.
- **Kim 2024 does not print the average-LIS optimal threshold** in the main text or figure
  legends; only the best-LIS value 0.21 is given. No average-LIS number should be quoted from
  this source.

---

## Changes this file makes to the repo's current numbers

| Value | `CLAUDE.md` | notebook cell 26 | **Canonical** | Change |
|---|---|---|---|---|
| `ipSAE_d0res` | 0.6 / 0.4 | 0.6 / 0.3 | 0.60 / 0.30 | amber unified to 0.30 |
| `ipSAE_d0chn` | — | 0.5 / 0.25 | 0.60 / 0.30 | raised; cell 26 had the bias backwards (§5) |
| `ipSAE_d0dom` | — | 0.55 / 0.28 | 0.60 / 0.30 | raised; same reason |
| `ipTM_d0chn` | 0.7 / 0.5 | 0.6 / 0.3 | 0.70 / 0.30 | green from `CLAUDE.md`, amber lowered to the measured overlap floor |
| `pDockQ` | 0.23 / 0.09 | 0.23 / 0.115 | 0.23 / 0.12 | amber now derived from the score's stated error |
| `pDockQ2` | 0.5 / 0.23 | 0.15 / 0.075 | 0.23 / 0.10 | green corrected to the published operating cutoff |
| `LIS` | 0.15 / 0.09 | 0.3 / 0.15 | 0.21 / 0.10 | green corrected to the published ROC-optimal cutoff |

Additionally, when this table is implemented:

1. **Amber must be read from the table, never recomputed.** Cell 28 currently derives amber as
   `green / 2` for `ipSAE_d0res`, `ipTM`, `pDockQ` and `LIS` while cell 26 reads the dict, so the
   summary table and the prose statement can label the same score differently in a single run.
   Both must consume one shared structure.
2. **`ipTM` should be relabelled `ipTM_d0chn`** everywhere in the notebook. It is not AlphaFold's
   ipTM, AFDB does not expose AlphaFold's ipTM on these endpoints, and thresholds published for
   AlphaFold's ipTM transfer only approximately (§4).
3. **Provenance should be visible in the UI.** Each traffic light should be able to state whether
   its threshold is published, derived, or a judgement call — six of the fourteen numbers in this
   table are judgement calls, and a user reading a coloured badge has no way to know that
   otherwise.
