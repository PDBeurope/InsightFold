# Cluster Quality Diagnostic Decision Capture

Date: 2026-09-21
Stage: Scope (lifecycle stage 3), produced by `$idea-scoping-interview` then `$scoping-decision-capture`
Next skill: `$concept-to-prd`
Working name: `cluster_quality_diagnostic`
Revision: v2. Three rounds of measurement during scoping overturned several v1 decisions. Superseded
decisions are marked rather than deleted, and the evidence is in the Measurement Log.

---

## 1. Core Problem Or Question

**A user has an AlphaFold model. They cannot tell whether its confidence is normal for its protein
family, what quality that family can achieve at all, or whether anything that looks better is even
the same protein.**

The AFDB entry page shows a single global pLDDT with no peer context. A user reading 71.4 has no way
to know whether that is roughly what this family always achieves, or whether they have landed on an
unusually poor member of an otherwise well-predicted family. Those two situations demand opposite
actions, and the website cannot currently distinguish them.

Measurement during scoping sharpened this considerably. Across a 40-cluster stratified random
sample, **79.8% of member-level pLDDT variance is between clusters, not within them** (ICC(1) =
0.798). Cluster means span 49.4 to 96.7 while the median within-cluster standard deviation is only
4.03. Protein families are internally homogeneous; the real variation is between families.

That makes the **family's achievable quality the primary finding**, and the user's position within
the family a secondary one that matters only in the tail. It also makes "you are typical for your
family, this is its ceiling, stop looking" a genuinely useful answer rather than a disappointing one.

Why it matters, and to whom:

- **Structural biologists and modellers** choose an input structure for docking, virtual screening,
  molecular dynamics or homology modelling. Knowing the family's ceiling tells them whether shopping
  around is worth any effort at all.
- **AFDB and PDBe internal analysts** need to see where prediction quality is systematically weak
  across family space, which is invisible from per-entry views.
- **Bench scientists** choosing which homologue to work on treat confidence as a rough proxy for
  tractability.
- **AFDB website users**, since this is a prototype for what the Similar proteins tab could show.

What this adds over the existing Similar proteins tab: that tab lists cluster members but does not
analyse them. It shows membership, not the quality distribution of that membership, does not locate
the user's own protein within it, and does not flag that some listed members are not the same
protein at all.

Broader InsightFold context: this is the second domain notebook after `homodimer_diagnostic`, and
the first to work at family scale. It reuses the repo's established shape (explicit data contracts,
pinned fixtures, bare-Figure plotting, MolViewSpec views, restart-and-run-all validation) against a
different AFDB API.

---

## 2. How The System Works

### Core operational logic

The user supplies one bare UniProt accession. The notebook fetches that protein's two AFDB
clusterings, reports what quality the family achieves, locates the user's protein within it,
characterises the low tail, then produces a screened and ranked shortlist of alternatives, and
finally shows the query, best and worst as structures, pairwise alignments and MSA coverage.

**One spine in three acts**, ordered so the family's level is established before any recommendation:

**Act A, the family.** Fetch both clusterings. Report the family's achievable quality as the
headline: central tendency, band composition, and the spread. Locate the query, with the honest
statement that position within the body of the distribution carries little information and that only
the tail is diagnostic. Characterise what the low tail is actually made of. Show taxonomic
composition as a descriptive and evolutionary reading of the family, carrying no explanatory claim.

**Act B, shortlist.** Apply the divergence screen, then four user-editable eligibility filters. Rank
survivors and present them with the length-deviation statistic alongside.

**Act C, the models themselves.** Render query, best and worst in Mol\*; align query against best and
against worst; show MSA coverage for those models.

### Key components and processes

| Component | Role |
|---|---|
| AFDB cluster members API | Both clusterings for the query accession. The only source for member lists. Requires ~2.5 s request spacing; see constraints |
| AFDB prediction API | Per-accession `sequence`, `cifUrl`, `globalMetricValue`, `msaUrl`. Called for the query plus highlighted models only, never for all members |
| Family quality summary | Central tendency and pLDDT band composition of the cluster. The headline output |
| Outlier check | Whether the query sits in the family's low tail, which is the only region where position is informative |
| Length deviation | Both raw length correlation and absolute deviation from cluster median. Members far from the median length in either direction score worse |
| Divergence screen | Flags members whose `uniprotDescriptions` diverge from the cluster's modal description, or whose length is atypical. Excluded from the shortlist by default |
| Taxonomic composition | Descriptive view of which taxa the family spans, split by pLDDT band. Genus level or higher |
| Eligibility filters | Comparable length, reviewed status, organism, reference proteome. All served directly by cluster API fields |
| MolViewSpec views | Query, best, worst, pLDDT-coloured. Collapses to two when the query is an extreme |
| Pairwise alignment | Biopython `PairwiseAligner`, rendered with pyMSAviz. Collapses to one alignment when the query is an extreme |
| MSA coverage | ColabFold-style occupancy plot in matplotlib, from the a3m at `msaUrl` or a user-supplied file |

### Inputs and outputs

**Input:** one bare UniProt accession, for example `Q9I1F6`. Optionally a user-supplied a3m per
highlighted model.

**Outputs:** the family's quality summary; whether the query is typical or an outlier; distribution
figures for both clusterings; length and length-versus-pLDDT figures; a taxonomic composition
figure; a characterisation of the low tail; a screened and ranked candidate table; two or three
Mol\* views; one or two pairwise alignments; MSA coverage where available.

### Metrics and signals

Average pLDDT per member, as reported by AFDB, is the single quality signal. Sequence length,
species and UniProt description are the covariates. **The notebook defines no new score.** Every
number shown traces back to a field AFDB already publishes, which is deliberate.

### Assumptions the system relies on

1. Average pLDDT is comparable across members of one cluster. Reasonable, since all members were
   predicted by the same pipeline.
2. Cluster membership implies homology strong enough that members are informative about one another.
   **Partly false, and the notebook addresses it directly**: clusters contain functionally divergent
   paralogues. See caveat C5 and decision D20.
3. The cluster members API returns the full cluster in one response. Verified to 93,793 members.

### Known limitations

- The notebook can say a model is better predicted. It cannot say a model is biologically
  appropriate for the user's purpose. See caveat C3.
- Average pLDDT is one number over a whole chain. A model with an excellent domain and a long
  disordered tail scores poorly and may still be exactly what the user needs.
- Within-family spread is small, so for most queries the answer is "you are typical". See caveat C9.

---

## 3. Constraints, Caveats, And Reliability

These are first-class notebook content, not footnotes.

### C1. Foldseek clusters are WIDER and LOWER-scoring than MMseqs2 clusters

**REVISED. The v1 version of this caveat asserted the opposite, and was wrong.**

**What:** Measured on a 40-cluster stratified random sample, restricted to clusters with at least 50
members on both flags (32 qualifying). Foldseek standard deviation exceeded MMseqs2 standard
deviation in **30 of 32** cases (sign test p = 2.5e-07), and the Foldseek mean was lower in **28 of
32** (p = 1.9e-05). Median standard-deviation ratio 0.66, median mean difference −1.48 pLDDT.

**Two mechanisms compete, and the widening one wins:**

*Tightening, representative selection.* Foldseek was run over the roughly 52.3M AFDB50
representatives only, and each sequence cluster's representative was chosen as its highest-pLDDT
member. This is real and verified from Barrio-Hernandez et al. 2023.

*Tightening, structural alignability.* Low-confidence models are extended and ribbon-like, align
poorly, and fail the E-value and 90% structural overlap criteria. Also real.

*Widening, superfamily aggregation.* A Foldseek cluster spans a structural superfamily across much
greater evolutionary distance, aggregating heterogeneous sequence families, whereas an MMseqs2
cluster is one homogeneous family at 50% identity. **Empirically this dominates the other two.**

**Why it matters:** this reverses the observation that motivated the project, and it reverses what
the v1 capture asserted. It is not a small-sample artefact: Foldseek clusters have the *smaller* n
(median 5.3 times smaller), and small-sample standard deviation is downward biased, so the bias
works against the finding.

**Communication:** the notebook states plainly that structure-based clusters are broader and more
heterogeneous, and does not present either clustering as the better one.

### C2. Length deviation predicts pLDDT, in both directions

**What:** On the random sample (31 clusters, n ≥ 200), median r(raw length, pLDDT) = −0.310 and
median r(absolute deviation from median length, pLDDT) = −0.250. Median length coefficient of
variation is 8.8%. Members more than 20% from the median length score far worse, with Cohen's d
reaching −8.41.

**Note on an intermediate finding:** a hand-picked seven-cluster round measured length CV at 2.3 to
6.3% and raw-length correlation near zero, and concluded the length confound was designed out by
AFDB50's 90% overlap criterion. The random sample contradicts this; those clusters were unusually
length-homogeneous. **The random-sample result stands.**

**Mitigation:** the comparable-length filter is on by default, the length-versus-pLDDT plots make the
trend visible, and the length-deviation statistic appears beside raw pLDDT in the candidate table.

### C3. Better predicted is not the same as biologically appropriate

**What:** A high-confidence bacterial homologue is useless to someone who needs the human protein.
The notebook has no knowledge of the user's downstream purpose.

**Why it matters:** the central overclaim risk of the notebook, and now compounded by C5, since a
"better" member may not even be the same protein.

**Communication:** a dedicated markdown cell immediately before the candidate table. The notebook
recommends in strong language, because that is what makes it useful, while stating plainly that
appropriateness is the user's judgement.

### C4. The cluster API does not reliably return a cluster containing the query

**What:** On AFDB/Foldseek the query was absent in 5 of 5 accessions tested, which follows from the
representative logic in C1. On AFDB50/MMseqs2 the behaviour is harder to explain: querying `Q9KM69`,
`A0A7X7SVK1` or `A0A536HFQ7` returns a byte-identical response to `Q9I1F6`'s 4,628-member cluster,
and none of those three appear in it. Same for `A0A8C8MD66` against `P69905`.

**Mitigation:** the notebook looks the query up explicitly and handles absence as a normal path. See
decision D9. The MMseqs2 behaviour is open question O2 and should be raised with the AFDB team.

### C5. Clusters contain functionally divergent paralogues

**What:** Cluster membership does not guarantee the same protein. Measured examples: a `P0A6F5`
GroEL cluster's low tail contains 1,535 T-complex protein 1 subunits, which are the eukaryotic CCT
chaperonin; a `Q13148` TDP-43 cluster's tail contains 400 RBFOX1 entries; a `P0AA25` thioredoxin
cluster's tail contains KaiB and glutaredoxin-like proteins; a `P69905` haemoglobin α cluster's tail
contains haemoglobin β.

**Why it matters, and this is a correctness risk rather than a cosmetic one:** Act B recommends
"use this better-predicted member instead". An unscreened shortlist can recommend a CCT subunit as a
substitute for GroEL. That is actively harmful advice.

**How strong is the signal:** weaker than first measured. On the random sample, modal-description
token divergence between the bottom 5% and the cluster body gives a median Cohen's d of only −0.33,
with d < −0.3 in 15 of 35 clusters and a *positive* d in 9 of 35. Median 57.7% of low-tail members
are flagged. Strong cases exist (`P9WN91` d = −1.08, `Q58847` −1.04, `P08204` −0.97). Divergent
paralogues are **one contributor among several, not the explanation for the low tail.**

**But the tail was the wrong place to look.** Round 4 measured the **top** of the distribution, which
is where the shortlist actually draws from, and found impurity there that is both severe and
harmful at rank 1: `Q13148`'s entire top 20 is yeast Prp24p at 18% identity to TDP-43. So the
correctness risk is real and measured, not inferred, and it is worse at the operating point that
matters than at the one first examined. Decision D20 is justified by that measurement.

**Mitigation:** decision D20, with percent identity to the query as the primary signal and
description divergence as the complement.

**Also:** AFDB50/MMseqs2 groups at a maximum of 50% identity with at least 90% bidirectional
sequence overlap against the representative's longest sequence; AFDB/Foldseek groups at E-value
below 0.01 with at least 90% bidirectional structural overlap against the representative's largest
structure. Neither is an InterPro family or homologous superfamily. The notebook uses InterPro's
vocabulary because that is how the audience thinks, and must state once that these clusters
approximate those concepts rather than instantiate them.

### C6. The MSA endpoint is currently unavailable

**What:** `msaUrl` returns HTTP 403 with a 134-byte HTML body, across multiple accessions and both
v4 and v6. Not 404, not 503, not JSON.

**Mitigation:** fallback detection keys on status and content type, never on parsing the body. The
user-supplied a3m path is the tested path and doubles as the workaround. The automatic path is
declared untested.

### C7. Taxonomic diversity does NOT explain within-family pLDDT variation

**REVISED. v1 proposed this as the notebook's explanatory spine. It does not survive measurement.**

**What:** Three independent tests, all null and sign-inconsistent across clusters. Environmental
versus named-isolate comparison gives Cohen's d of +0.70 on `P0A6F5`, 0.00 on `P0AA25` and −0.49 on
`Q9I1F6`; four of seven clusters have no environmental members at all; low-tail enrichment
contradicts itself in direction. Genus rarity correlates with pLDDT at r between −0.030 and +0.338,
mostly noise. Phylum-level band fractions are indistinguishable.

**Why it matters:** the notebook must not present taxonomic diversity as evidence for an MSA-depth
mechanism. The taxonomic view survives as description only. See decision D19.

### C8. pLDDT is not an experimental tractability predictor

High confidence reflects the predictor's certainty, not how a protein behaves in a tube. One
explicit sentence in Act B.

### C9. A percentile is nearly uninformative in the body of the distribution

**What:** On the random sample, the median p25 to p75 width is **4.2 pLDDT points** and the median
p10 to p90 width is 8.2. Moving from the 25th to the 75th percentile of your family is within
rounding of anything a user cares about. But the tail is far out: the median gap from the cluster
median to its 1st percentile is **11.8 points, or 3.4 standard deviations**, and a median of 1.8% of
members sit more than 10 points below the median (maximum 21.9%).

**Why it matters:** "am I better than average for my family" is close to meaningless. "Am I an
outlier" is answerable. The notebook must not dress up a 4-point spread as a meaningful ranking.

**Communication:** the query's position is reported as a category (typical, or in the low tail) with
the percentile shown but explicitly de-emphasised, and the band width stated so the reader can see
how little separates the body.

### Operational constraints

- Must run on Google Colab free tier, restart-and-run-all, within roughly 60 s of install time.
- **The two AFDB endpoints have different rate policies and must not share one.** The **cluster**
  endpoint rate-limits: 40 accessions at 0.8 s spacing produced 35 spurious HTTP 500 responses, while
  the same requests at **2.5 s spacing produced zero failures**. The **prediction** endpoint needs no
  spacing at all: 45 back-to-back requests at 0.0 s spacing gave zero failures at 0.15 s median
  latency. This is what makes decision D20's per-candidate identity check affordable, at 3 to 4
  seconds for 20 rows. The cluster-endpoint finding also casts doubt on the earlier conclusion that
  `P12345` has a reproducible per-accession 500; see open question O12.
- Cluster size is not a performance problem: the largest measured, `P0A6F5`, returned 93,793 members
  at 8.82 MB in 1.0 s. No cap, no pagination, no sampling.
- The cluster endpoint is public but unadvertised, and a rate-limited mirror will replace it. The URL
  will appear in plain text in a shared notebook, so "unadvertised" cannot mean "hidden". See D15.

---

## 4. Scope Decisions

### D1. REVISED. The core question is the family's ceiling, plus an outlier check

We decided the notebook's primary question is *"what quality can this family achieve, is my model
typical of it, and is anything that looks better actually the same protein"*, rather than v1's
*"is my model an outlier in its family"*. The reason is measured: ICC(1) = 0.798 means almost four
fifths of variance is between families, so the family's level is where the information is, and the
query's position within the family is informative only in the tail. **Tradeoff accepted:** the most
common honest output becomes "you are typical, this is your family's ceiling", which is a less
exciting headline than v1 promised but is what the data supports.

### D2. All four audiences are served, through one spine

We decided against nominating a single primary user. Act A serves the analyst and method-comparison
readings, Act B the modeller and bench scientist. **Tradeoff accepted:** a longer notebook, and a
modeller scrolls past the landscape to reach the table.

### D3. REVISED. Outlier status is tail-focused, and reported with both length statistics

We decided to report the query's percentile but de-emphasise it per caveat C9, to state the family's
band width so the reader sees how little separates the body, and to make the headline judgement
categorical: typical, or in the low tail. Both raw length correlation and absolute deviation from
median length are reported, since the random sample supports both (median r = −0.310 and −0.250).
**Tradeoff accepted:** two length statistics rather than one, and a headline that is a category
rather than a number.

### D4. Four eligibility filters, all served by the API

Comparable length, reviewed (SwissProt) status, organism, reference proteome membership, because
those are exactly the four discriminating fields the cluster API returns beyond pLDDT. **Tradeoff
accepted:** no functional or structural-coverage filtering in v1. Organism matching likely needs a
taxonomic level rather than exact species, since exact match will frequently return nothing.

### D5. Input is a single bare UniProt accession

One accession, both clusterings fetched automatically with no user-facing switch, matching the
`homodimer_diagnostic` precedent. Batch input, local files and non-AFDB proteins are out of scope.
**Tradeoff accepted:** AFDB internal triage across many proteins is unsupported in v1. Note the API
accepts bare accessions only; `AF-Q9I1F6-F1` returns 404.

### D6. Strong recommending language, with substitutability first-class

We decided to shortlist and rank rather than merely describe, because a purely descriptive notebook
does not answer the question users have. **Tradeoff accepted:** real overclaim risk, accepted
knowingly, mitigated by caveats C3 and C5 and by decision D20.

### D7. REVISED. MSA work covers the highlighted models only

We decided against per-member MSA retrieval, which is infeasible at cluster scale, and in favour of
genuine MSA coverage plots for the two or three highlighted models. Users may supply their own a3m.
**The v1 plan to use taxonomic diversity as a cluster-wide MSA-depth proxy is dropped**, per caveat
C7. **Tradeoff accepted:** the notebook offers no cluster-wide account of MSA depth at all.

### D8. REVISED. The two clusterings answer different questions, with the characterisation reversed

We decided AFDB50/MMseqs2 is presented as **one homogeneous sequence family**, tighter and
higher-scoring, and AFDB/Foldseek as **a structural superfamily spanning greater evolutionary
distance**, broader, lower-scoring and more heterogeneous. They are not plotted as a like-for-like
comparison. **Tradeoff accepted:** v1 asserted the reverse characterisation on the basis of a
plausible mechanism that measurement did not support, and the appealing "structure clustering
predicts quality better" reading is given up, as is the founding observation that sequence
clustering shows the wider spread.

### D9. Query absence is handled by marking the representative and explaining the hop

Where the query is absent, the notebook states that the protein belongs to sequence cluster X whose
representative is R, marks R on the Foldseek distribution, and draws the query as a separate
reference line. **Tradeoff accepted:** an extra concept for the reader.

### D10. Cluster-wide extremes drive the figures; filters drive the table only

Act A shows cluster-wide best and worst, and the Mol\* views and alignments use those same extremes,
while filters govern only the Act B table. Keeps figures stable when a filter changes. **Tradeoff
accepted:** a rendered structure may be one the user's filters excluded. **Note:** the divergence
screen of D20 is an exception and applies to the worst-model selection too, since rendering a CCT
subunit as "the worst GroEL" would be actively misleading.

### D11. Ranking is on raw average pLDDT, with the bias made visible

Rank on raw average pLDDT so "best" matches the AFDB website, and handle length bias by plotting
lengths explicitly and showing the length-deviation statistic alongside. **Tradeoff accepted:** the
top-ranked candidate is sometimes not the best-predicted-for-its-size.

### D12. Fewer than three eligible candidates warns and still shows the table

Act A runs in full; Act B names which filter removed the most candidates, suggests relaxing it, and
shows what survived. **Tradeoff accepted:** the user may see an unhelpful table. Rejected: silently
auto-relaxing filters.

### D13. Dependencies: keep pandas and seaborn, add biopython and pyMSAviz, drop BigQuery and ptitprince

numpy, pandas, matplotlib, seaborn, requests, molviewspec, biopython, pyMSAviz. The repo's pandas
prohibition is per-consumer and binds the dimer notebooks, not this one. BigQuery is dropped because
the cluster API carries every field needed. **Tradeoff accepted:** roughly 6 s install, within the
Colab budget; ptitprince retained as fallback only.

### D14. Visualisation uses pyMSAviz and static rendering, no JavaScript widget

pyMSAviz accepts an in-memory alignment, returns a bare `matplotlib.figure.Figure` matching repo
convention, and takes a custom residue-to-colour dict from the user's `cleancolors.json`. Navigation
is by paging, because Colab's custom widget manager is a documented failure point. Taxonomic
composition is a **horizontal stacked bar**, not a Sankey; see D19. **Tradeoff accepted:** no
interactive scrolling MSA.

### D15. The cluster endpoint lives behind a single constant and is not named in prose

Base URL defined once so the switch to the rate-limited mirror is a one-line change; no markdown
cell promotes the endpoint. Request spacing of at least 2.5 s is built in. **Tradeoff accepted:**
the URL is readable in shared source, so this achieves "not advertised", not "secret". The hardcoded
API key and test-server base from the source notebook are not carried over.

### D16. InterPro family-coherence annotation is deferred to v1.1

**Tradeoff accepted:** the clearest differentiator from the Similar proteins tab is not in v1.

**Correction, round 4.** An earlier revision of this decision claimed InterPro would be a far better
divergence detector than description token overlap. **That was wrong.** InterPro answers *family
membership*, not *same protein as the query*: GroEL and the CCT chaperonin plausibly share a
chaperonin superfamily entry and would both pass, which is precisely the case decision D20 must
catch. It is a better family detector and a worse substitutability detector, and it would be a third
API with its own limits, unaffordable across thousands of members. Percent identity to the query is
cheaper, query-anchored and strictly better suited. InterPro stays deferred, for cluster-coherence
description rather than for screening.

### D17. Lifecycle position: prototype, promotable to the Similar proteins tab

Internal planning note; must not appear as notebook prose.

### D18. REVISED. Act B ranks the MMseqs2 cluster; the Foldseek toggle is demoted

We decided the candidate pool is the AFDB50/MMseqs2 cluster. **The v1 rationale for offering the
Foldseek set as a "curated representatives" pool is withdrawn**: caveat C1 shows Foldseek clusters
are broader and more heterogeneous, which makes them *more* exposed to the divergent-paralogue risk
of caveat C5, not less. The toggle may remain as an exploratory option but must not be described as
curated or higher-quality. **Tradeoff accepted:** the appealing "pre-filtered pool of best
representatives" framing is lost.

### D19. REVISED. Taxonomic composition is descriptive only, and is a stacked bar

We decided Act A includes a taxonomic composition view split by pLDDT band, **as description and as
the light evolutionary reading the user asked for, carrying no explanatory claim** per caveat C7.
Rendering is a horizontal stacked bar in pure matplotlib, taxon rank on the y axis and pLDDT band as
the stack, returning a bare `Figure` and reusing the existing band palette. Rank is selectable,
phylum by default and class for eukaryotic clusters where phylum is uninformative.

A Sankey was the user's initial preference and is **rejected**: Sankeys show flow between levels,
whereas this is composition split by one ordinal variable; hundreds of genera collapse the labels
into noise; and Plotly's Sankey is a JavaScript widget that D14 already rules out.

Resolution is genus or higher. Long-tail rule: keep the top 6 taxa by member count, or the smallest
set covering 90% of members, whichever is fewer, and aggregate the rest into a single
`Other / unresolved` row, keeping unresolved genera visibly distinct from genuinely rare taxa.
**Tradeoff accepted:** a view that describes rather than explains, and a diagram form the user did
not initially ask for.

### D20. REVISED. The identity screen is query-anchored, two-signal, on by default and disableable

We decided the shortlist is screened for "is this candidate the same protein as the query", on by
default, with the removal count stated and a toggle to reveal what was removed. The user-facing
behaviour is as originally decided. **Everything about the mechanism changed in round 4**, and the
v1 rationale is withdrawn as motivated reasoning: it substituted the *severity* of the harm for
*evidence that the screen reduces it*. Round 4 supplies the missing evidence, and it is strong.

**Justified at the right operating point.** The harmful case occurs at rank 1, not in the tail: a
user holding TDP-43 (`Q13148`) would be handed yeast Prp24p at 18% identity as the single best
recommendation. Measured, and verified not to be the caveat C4 lookup anomaly.

**Anchored on the query, never on the cluster mode.** The target is "same protein as the query", not
"typical of this cluster". Mode anchoring inverts when the query is a minority member, which is the
common case: query and modal descriptions diverge in 4 of 6 clusters measured, and on `P69905` the
modal description is haemoglobin *beta*, so mode anchoring flags the genuine alpha chains. Cluster
mode is a fallback only when the query is absent from the returned cluster, per caveat C4.

**Two complementary signals, neither sufficient alone.** Percent identity to the query, computed with
Biopython `PairwiseAligner` over the shortlist candidates only, plus query-anchored description
divergence. Identity catches mixed clusters where impure and genuine members coexist; description
catches uniformly impure clusters where identity has no relative signal to work with. Percent
identity is shown as a first-class column beside pLDDT, so each row carries its own falsifiable
evidence rather than an assurance.

**No absolute identity cut.** A flat threshold cannot work: 25% correctly excludes `Q13148`'s entire
top 20 but kills 19 of 20 plausibly legitimate LacI-family regulators in `Q9I1F6`; 30% kills 8
genuine thioredoxins. Flagging is on **within-cluster relative** identity combined with description
divergence. Calibration is open question O11a.

**Applicability gate, because the screen can invert on unannotated clusters.** The modal description
is computed over *informative* descriptions only, excluding "Uncharacterized protein" and similar,
since otherwise at high unannotated fractions the mode becomes "Uncharacterized" and the screen flags
the well-annotated members. The screen self-disables below an informative-fraction threshold, prints
a named banner, falls back to identity alone, and marks the table's identity column unverified.
Uninformative members are treated as **unevidenced, not evidenced-same**.

**Two operating points, kept separate.** The shortlist screen runs at the top of the distribution.
The decision D10 exception that stops a divergent member being rendered as "the worst model" runs at
the bottom, where the only available signal is weak (median Cohen's d −0.33). At the bottom the
notebook says a member is *description-divergent*, which describes the output, and does not claim it
has identified a paralogue.

**Affordable.** The prediction endpoint requires no request spacing (45 back-to-back requests, zero
failures, 0.15 s median latency), so a 20-row identity check costs 3 to 4 seconds. The 2.5 s spacing
requirement applies to the cluster endpoint only.

**Tradeoff accepted:** two signals and a gate is more machinery than a single heuristic, and the
combined rule still fails on `Q9I1F6`, where neither signal separates GntR from LacI. Rejected: hard
non-disableable exclusion, which gates on a heuristic with no expert override; warn-only, which
leaves a harmful recommendation at the top of the table; and the v1 mode-anchored description-only
screen, which is actively wrong in two thirds of the clusters measured.

### Refusal gates decided

| Condition | Behaviour |
|---|---|
| Accession in `AF-{ACC}-F1` form | Normalise to bare accession, or reject naming the expected form |
| Accession absent or malformed | HTTP 404 with a service message. Report it, quoting the service |
| Multi-fragment protein, e.g. `Q8WZ42` | 404 on both flags. Fragmented proteins are absent from clustering. Refuse with a clear explanation |
| Isoform, e.g. `P04637-2` | Supported; resolves to the parent's cluster. Match the exact `AF-{ACC}-F1` form, since prefix matching wrongly catches isoforms |
| Singleton or fewer than three members | Refuse: the notebook needs at least three members. Name the cluster size and suggest another accession |
| HTTP 5xx | Retry with backoff before reporting, since 5xx is often rate limiting rather than a fault. See open question O12 |
| `msaUrl` returns 403 | Fall back to the user-supplied a3m; state the endpoint is under maintenance |

---

## 5. User Journey / Narrative Flow

**Stage 1, entry.** The user edits one parameter cell containing a UniProt accession and runs all.
The notebook validates and normalises it, and either proceeds or refuses with a named reason.
*Trust point:* refusals name the actual condition, never a bare traceback.

**Stage 2, orientation.** The notebook explains the two clusterings, quoting the AFDB definitions,
relating them to InterPro's family and homologous superfamily concepts, and stating that they
approximate rather than instantiate those concepts. *Trust point:* caveat C5's first half lands here.

**Stage 3, the family's ceiling.** The headline output: what quality this family achieves, as central
tendency and band composition, for the sequence cluster. *The user learns* whether they are working
in a family AlphaFold predicts well at all. *Interpretation happens here*, and this is the stage that
carries the most information, since between-family variation dominates.

**Stage 4, where the query sits.** The query is located in the distribution, reported categorically
as typical or in the low tail, with the percentile shown but de-emphasised and the band width stated.
*Trust point:* caveat C9. The notebook must not dress up a 4-point spread as a ranking.

**Stage 5, what the low tail is made of.** Length-atypical members and description-divergent members
are characterised. *The user learns* that some cluster members are not the same protein, which is
both interesting in itself and the justification for the screen they are about to meet. *Trust
point:* caveat C5, stated with its real effect size rather than oversold.

**Stage 6, the structural superfamily.** The Foldseek set is introduced as broader and more
heterogeneous, spanning greater evolutionary distance. The query's representative is marked and the
hop explained; the query is drawn as a reference line. *Trust point:* caveats C1 and C4.

**Stage 7, the family across the tree of life.** Taxonomic composition, split by pLDDT band, as
description. *The user learns* which taxa the family spans. *Trust point:* caveat C7. The notebook
must not claim this explains the quality variation.

**Stage 8, candidates.** The divergence screen runs, reporting what it removed. The user adjusts four
filters and reads a ranked table of accession, description, species, length, average pLDDT and length
deviation. *Decision-making happens here.* *Trust point:* caveats C3 and C8 immediately above the
table, and D20's removal count visible rather than silent. Fewer than three survivors triggers the
D12 warning.

**Stage 9, looking at the models.** Query, best and worst in Mol\*, pLDDT-coloured, collapsing to two
when the query is an extreme. *The user learns* where in the structure the differences sit.

**Stage 10, sequence-level comparison.** Pairwise alignments of query against best and worst,
collapsing to one when the query is an extreme, with percent identity reported. *The user learns* how
similar the recommended alternative actually is, which is the evidence caveat C3 hands them.

**Stage 11, MSA evidence.** Coverage plots for the highlighted models, from `msaUrl` or a supplied
a3m. A supplied file is accepted only when its gap-stripped target sequence matches the AFDB sequence
for that accession; the `AF-{ACC}-F1` header convention is a hint that speeds lookup, not the
validation itself.

**Stage 12, exit.** A short summary restating the family's ceiling, whether the query is typical, the
top candidate if any, and the caveats that bound the recommendation.

---

## 6. Open Questions And Deferred Decisions

**O1. RESOLVED, 2026-09-21. Does taxonomic diversity explain within-family pLDDT spread?**
No. Three independent tests, all null and sign-inconsistent. See caveat C7 and decision D19. The
taxonomic view survives as description only.

**O2. Is the MMseqs2 lookup behaviour a bug, and should it be raised with the AFDB team?**
Open. An API-side question, not a notebook decision; the notebook tolerates it either way per D9.
Evidence in caveat C4.

**O3. RESOLVED, 2026-09-21. Which pool does Act B rank?** The MMseqs2 cluster. See decision D18.

**O4. What is the default comparable-length tolerance?**
Now informed but not decided. Median length CV is 8.8%, and members more than 20% from the median
score far worse (d to −8.41), which suggests ±20% as a defensible default. Needs confirming against
the fixture set.

**O5. Which MSA file formats does the user-supplied path accept?**
Open. a3m is required since that is what `msaUrl` serves. Whether plain aligned FASTA is also
accepted is undecided. a3m lowercase characters denote insertions and must be stripped before the
alignment is rectangular, so the two formats need different handling.

**O6. Is the flat `clustal` scheme acceptable, or is true Clustal X colouring needed?**
Open. All twelve schemes in `cleancolors.json` are flat per-residue lookups, but the one named
`clustal` is the simplified Jalview variant rather than conservation-threshold Clustal X. Real
Clustal X means writing the column conservation rules by hand.

**O7. Colab has not been verified.**
Open, deferred to validation. The repo has been burned before: every prior `homodimer_diagnostic`
validation run had `IN_COLAB = False`, leaving the clone and install paths unexercised. Until an
actual Colab run happens, the Colab path is declared untested.

**O8. Does `clusterTotal` ever disagree with the returned array length?**
Open. Always equal so far. Matters because if the query is a member the response omits, true cluster
size and every percentile is off by one.

**O9. How much taxonomic resolution is worth fetching?**
Open, and lower priority now that D19 is descriptive only. Genus level is viable: batched UniProt
taxonomy queries at batch size 20 resolved 53 of 60 genera covering 75.4% of members in 1.82 s across
3 calls, with full phylum, class, order and family lineage in the same payload. Batch size 30
silently drops hits as results crowd out the `size` window, so 20 is the safe ceiling. Species level
is not viable (3,274 unique names for `Q9I1F6` at about 0.12 s each is roughly 390 s serial) and is
unsafe besides, since `"Streptomyces sp."` fuzzy-matches to an arbitrary named species. There is no
bulk taxId route from AFDB: the cluster endpoint ignores `fields`, `include`, `taxonomy` and
`format`, and the search endpoint carries `taxId` but rejects OR queries with a 404. An offline
taxdump is 60 MB or more and blows the Colab budget.

**O10. NEW. Is the Foldseek widening due to superfamily aggregation or to cluster-size disparity?**
Open. The effect scales with how much the sequence cluster outnumbers the structure cluster
(r(log(mm_n/fs_n), sd_ratio) = −0.688), which is consistent with the aggregation explanation but does
not separate it from a size confound. The notebook's *description* of the two clusterings depends on
which it is. Resolved by: a design that controls for cluster size, or by advisory judgement.

**O11. NEW. What divergence threshold, and how to handle uninformative descriptions?**
**RESOLVED in part, 2026-09-21, and replaced by O11a to O11c below.** The question as posed was
mis-framed: Cohen's d of −0.33 was measured between the bottom 5% and the cluster body, which is an
*explanatory* statistic at an operating point the shortlist never visits. Round 4 measured the top
of the distribution and found impurity that is severe and harmful at rank 1, so decision D20 is
justified on evidence rather than on severity. The surviving sub-questions are calibration, not
existence.

**O11a. NEW. How is the identity plus description flag calibrated?**
Open, and the last substantive thing blocking the notebook's recommending language. No absolute
identity cut works: 25% correctly excludes `Q13148`'s entire top 20 but kills 19 of 20 plausibly
legitimate LacI-family regulators in `Q9I1F6`; 30% kills 8 genuine thioredoxins; 40% kills
everything. The rule must combine within-cluster relative identity with query-anchored description
divergence. Resolved by: tuning against the fixture set once O11b's hand label exists.

**O11b. NEW. Is `Q9I1F6`'s GntR-versus-LacI top 20 genuine impurity or a naming collision?**
Open, and it decides whether O11a is tunable at all. It is the one measured case where neither
signal separates the classes, and the reading was from descriptions rather than expert labelling.
Resolved by: a hand label from someone who knows the GntR and LacI regulator families. **This is the
single highest-value piece of domain input outstanding.**

**O11c. NEW. What are the applicability-gate thresholds?**
Open. Decision D20 disables the description signal when too few members are informatively annotated,
but the informative-fraction and modal-coverage thresholds are not set. `C1C553` (60.0%
uninformative) is the fixture for this. Resolved by: tuning against that fixture.

**O11 original text, retained for context.**
Modal-description token divergence is weak on
average (median Cohen's d −0.33, d < −0.3 in only 15 of 35 clusters, positive in 9 of 35), and
captures a median 57.7% of low-tail members. Uninformative descriptions such as "Uncharacterized
protein" are a median 2.7% of members but reach 68.7% in the worst cluster measured, and the
reference-proteome stratification of that sample over-samples well-annotated organisms, so the true
rate is likely worse. Resolved by: choosing a threshold against the fixture set, and deciding what
the screen does when descriptions are mostly uninformative. Note that InterPro annotation (D16,
deferred) would be a far better detector.

**O12. NEW. Is `P12345`'s HTTP 500 a per-accession fault or rate limiting?**
Open. It was recorded as reproducible (2 of 2) on MMseqs2 while Foldseek succeeded, but the random
sample later showed that 0.8 s request spacing produces spurious 500s that vanish at 2.5 s. The
earlier observation may have been rate limiting. Resolved by: re-testing `P12345` in isolation with
generous spacing. Affects whether it remains a useful fixture.

### Explicitly deferred features

- InterPro family-coherence annotation, to v1.1 (D16), now with added value per caveat C5.
- Batch input across multiple accessions, which would serve AFDB internal triage.
- Local mmCIF or non-AFDB protein support.
- An interactive, navigable MSA viewer (D14).
- Any cluster-wide account of MSA depth (D7).

---

## 7. Success Criteria / Graduation Signal

### The notebook is successful at v1 if

1. A user supplying any pinned fixture accession gets a correct family summary, query position,
   distribution figures, screened candidate table, structures and alignments in one
   restart-and-run-all pass.
2. Every refusal gate fires with a named, human-readable reason rather than a traceback.
3. Caveats C1, C3, C5 and C9 appear as prose a reader cannot skip. In particular, a reader who looks
   only at the figures must not come away believing that a 4-point percentile spread is a meaningful
   ranking, nor that structure clustering predicts quality better.
4. It runs on Colab free tier within the install budget.
5. **Named negative controls pass:** no Prp24p appears in `Q13148`'s default shortlist, no
   *T-complex protein 1* or *Thermosome* in `P0A6F5`'s, no haemoglobin beta in `P69905`'s. These are
   string assertions with real biological claims behind them.
6. **The false-positive control passes:** `P0AA25` retains its genuine thioredoxins despite their
   27 to 39% identity to the query, and `P69905` retains its alpha chains. Without this, a screen
   that removed half of every cluster would satisfy every other criterion.
7. **The applicability gate fires on `C1C553`:** the screen self-disables, prints its banner, falls
   back to identity alone, and marks the identity column unverified.
8. **The honesty invariant holds:** the printed removal count equals the number of suppressed rows,
   and the toggle restores exactly those rows. This is what criterion 6 of v2 actually asserted.
9. **Cluster API responses are pinned as fixture files**, so none of the above depends on a live,
   rate-limited, soon-to-be-replaced endpoint.
10. A structural biologist reading the output can state, in their own words, what quality their
    family achieves and whether their model is typical of it. **This needs a protocol before it is a
    criterion**: how many readers, which fixtures, what counts as a pass. Until then it belongs in
    graduation signals rather than here.

**Not success criteria.** Do not assert Cohen's d values, low-tail flag rates or identity
distributions in tests. Those are properties of AFDB's data rather than of this code, and they will
drift with AFDB releases. Record them as monitored measurements instead.

### Graduation signals

- **Toward the Similar proteins tab:** the family-ceiling summary and the divergence screen are
  judged useful enough to prototype as a web component, and the numbers can be precomputed per entry
  rather than fetched live.
- **Toward a standing notebook:** repeat use by people other than the author, and requests for batch
  mode, which is the natural second version.
- **Toward further investment:** evidence that the divergence screen surfaces genuine cluster
  impurity that AFDB would want to act on, which would make this a data-quality tool as much as a
  user-facing one.

### Stop or reconsideration signals

- If the divergence screen proves too weak to be trusted (open question O11), Act B's safety argument
  collapses and decision D6's recommending language must weaken substantially.
- If the MMseqs2 lookup behaviour in caveat C4 is widespread rather than occasional, the family
  summary is unreliable for an unknown fraction of queries and the core claim is undermined.
- If users consistently read the percentile as a ranking despite caveat C9, the statistic should be
  removed rather than caveated.
- If users act on recommendations across organism boundaries where they should not, decision D6's
  language must weaken.

---

## Measurement Log

Three rounds during scoping, 2026-09-21. Rounds 1 and 2 used hand-picked clusters; round 3 used a
stratified random sample and is the authority where they disagree.

### Round 3 sampling procedure (the authoritative round)

Stratified random draw, seed 20260921. Eight reference proteomes spanning bacteria (*E. coli* K12,
*B. subtilis* 168, *M. tuberculosis* H37Rv), archaea (*M. jannaschii*) and eukaryotes
(*S. cerevisiae*, *A. thaliana*, *H. sapiens*, *D. melanogaster*). Five uniform-random index offsets
per organism via the AFDB search endpoint filtered to
`isComplex:false AND isUniProtReferenceProteome:true`, offsets capped at 40,000 to avoid deep-paging
504s. 40 of 40 accessions resolved; 80 of 80 cluster fetches returned HTTP 200. No overlap with the
hand-picked sets. **Known bias:** reference-proteome stratification over-samples well-annotated
organisms, so the uninformative-description rates in caveat C5 are likely optimistic.

### What each round established

| Finding | Round 1 (2 clusters) | Round 2 (7 clusters, hand-picked) | Round 3 (40 random) | Status |
|---|---|---|---|---|
| Taxonomic / MSA-depth proxy | null | null, sign-inconsistent | not re-tested | **Dead.** C7 |
| Raw length correlation | r ≈ 0 | r sign-inconsistent, CV 2.3–6.3% | **median r = −0.310, CV 8.8%** | **Real.** C2 |
| Absolute length deviation | not tested | r to −0.517, AUC to 0.82 | median r = −0.250 | **Real.** C2 |
| Within vs between variance | not tested | max sd 6.4, means 41.6–96.3 | **ICC(1) = 0.798** | **Confirmed.** C9, D1 |
| Foldseek wider and lower | not tested | 6 of 7 | **30 of 32, p = 2.5e-07** | **Confirmed.** C1 |
| Divergent paralogues in tail | not tested | 74–96% capture, by eye | **median d = −0.33, 15 of 35** | **Partial.** C5, O11 |
| API rate limiting | not observed | not observed | **0.8 s fails, 2.5 s clean** | **Confirmed.** Constraints |

Round 2 is superseded on length by round 3: its two length-homogeneous clusters produced a CV of
2.3 to 6.3% that does not generalise, and the conclusion drawn from it, that AFDB50's 90% overlap
criterion designs the length confound out, was wrong.

### Round 3 detail, Foldseek reversal

Restricted to n ≥ 50 on both flags, 32 of 40 qualify. Foldseek sd exceeds MMseqs2 sd in 30 of 32
(sign test p = 2.5e-07); Foldseek mean is lower in 28 of 32 (p = 1.9e-05). Median sd ratio (mm/fs)
0.66, IQR 0.52 to 0.81. Median mean delta −1.48 pLDDT, range −5.64 to +1.08. Mean sd: MMseqs2 4.12,
Foldseek 5.99. Foldseek has the smaller n (median 5.3 times smaller), and small-sample sd is
downward biased, so the bias works against the finding.

### Round 3 detail, variance structure

35 clusters, 826,697 members. ICC(1) = 0.798. Cluster means span 49.4 to 96.7 (sd 12.46); median
within-cluster sd 4.03. Median p25–p75 width 4.2 pLDDT points; p10–p90 8.2; median-to-p5 gap 7.1;
median-to-p1 gap 11.8, which is 3.4 sd. A median of 1.8% of members sit more than 10 points below
the median, maximum 21.9%.

### Round 4: impurity at the shortlist's operating point, and it is severe

An independent evaluation review identified that every impurity observation to that point came from
the **bottom 5%**, while the shortlist draws from the **top**, so the correctness risk justifying
decision D20 was asserted from a region the shortlist never touches. Round 4 measured the right
operating point: percent identity to the query across the top 20 by pLDDT, six fixture clusters,
Biopython `PairwiseAligner`, BLOSUM62, affine gaps −11/−1, identity over global alignment columns.
Ungapped and shorter-sequence conventions differ by only 2 to 6 points, so the conclusion is
convention-independent, and member lengths are well matched, so it is not a global-alignment artefact.

| Cluster | Query | Top-20 median identity | min | max | <25% | <30% | <40% |
|---|---|---|---|---|---|---|---|
| `Q13148` | TDP-43 | **18.0** | 16.8 | 19.0 | 20 | 20 | 20 |
| `P9WN91` | fumarate reductase | **19.1** | 18.1 | 53.3 | 15 | 15 | 15 |
| `Q9I1F6` | GntR regulator | **21.0** | 18.3 | 26.1 | 19 | 20 | 20 |
| `P0AA25` | thioredoxin 1 | 31.2 | 27.0 | 39.4 | 0 | 8 | 20 |
| `P0A6F5` | GroEL | 54.9 | 22.8 | 59.5 | 2 | 2 | 2 |
| `P69905` | haemoglobin α | 86.6 | 34.5 | 100.0 | 0 | 0 | 1 |

**The harmful scenario occurs at rank 1, not in the tail.** `Q13148`'s entire top 20 is yeast
Prp24p / U6 snRNP protein at 18% identity, a different protein with a different domain architecture
(four RRMs against TDP-43's two). A user holding TDP-43 would be handed Prp24p as the single best
recommendation. `P9WN91`'s top 20 mixes dihydrolipoyl dehydrogenase, dimethylglycine oxidase and
aminomethyltransferase at 18 to 20% with genuine fumarate reductase at 47 to 53%. `P0A6F5` carries
*Thermosome subunit alpha*, an archaeal group II chaperonin, at 22.8% at rank 10.

**Verified 2026-09-21 that this is genuine impurity and not the lookup anomaly of caveat C4.** The
query is present in the returned cluster in all six cases: `Q13148` at index 11,493 of 33,211 with
the head of the cluster being Prp24p at pLDDT 86.81 against the query's 65.19; `P9WN91` at 1,741 of
32,580; `Q9I1F6` at 1,650 of 4,628; `P0A6F5` at 11,255 of 93,793; `P69905` at 17 of 5,223; `P0AA25`
at 13,792 of 62,441.

**A finding for the AFDB team in its own right.** AFDB50/MMseqs2 clusters drift far beyond their
nominal 50% identity criterion at scale, because the criterion binds each member to the
representative and transitivity accumulates across large clusters. A 33,211-member cluster spanning
both TDP-43 and Prp24p is the clearest case measured. This is worth raising alongside the caveat C4
lookup anomaly.

**The inversion bug is confirmed, and is the common case, not an edge case.** The query's description
diverges from the cluster's modal description in **4 of 6** clusters.

| Cluster | Query-anchored median Jaccard / n<0.3 | Mode-anchored / n<0.3 | Mode versus query |
|---|---|---|---|
| `P69905` | **1.00** / 1 | 0.33 / 1 | α against **β** |
| `P9WN91` | 0.00 / 16 | 0.00 / **20** | fumarate reductase against **MnmG** |
| `Q13148` | 0.00 / 20 | 0.00 / 18 | TDP-43 against RRM-domain |
| `Q9I1F6` | 0.10 / 16 | 0.12 / 11 | GntR against LacI |
| `P0A6F5` | 1.00 / 4 | 1.00 / 4 | same |
| `P0AA25` | 1.00 / 1 | 1.00 / 1 | same |

`P69905` is the clean demonstration: the modal description is haemoglobin **beta** at 20.7% of the
cluster, so mode anchoring flags the genuine alpha chains. `P9WN91` mode anchoring flags 20 of 20,
including both legitimate fumarate reductases.

**No flat identity threshold works.** A 25% cut correctly excludes `Q13148` (20 of 20), `P9WN91`
(15 of 20) and the thermosome, and correctly excludes nothing from `P69905` or `P0AA25`, but kills
19 of 20 of `Q9I1F6`'s LacI-family regulators, which are plausibly legitimate distant homologues.
30% additionally kills 8 genuine thioredoxins; 40% kills all 20 plus all of `Q9I1F6`.

**The two detectors are complementary, not redundant.** Query-anchored description catches the
uniformly impure case (`Q13148`, 20 of 20 at Jaccard 0) that identity cannot separate, because when
everything is low-identity there is no relative signal. Identity catches the mixed cases (`P9WN91`,
`P0A6F5`) where impure and genuine members coexist. Both fail on `Q9I1F6`.

**Rate tolerance, measured per endpoint.** The **prediction endpoint needs no spacing at all**: 12
requests each at 0.2, 0.5 and 1.0 s spacing gave zero failures, then 45 back-to-back requests at
0.0 s spacing gave zero failures at a median latency of 0.15 s. A 20-row identity check costs about
3 to 4 seconds, which is trivially affordable inside the Colab budget. The **cluster endpoint still
requires 2.5 s spacing.** The two endpoints must not share a rate policy.

**Limits.** Six clusters, chosen as existing fixtures rather than at random, so this establishes that
top-of-distribution impurity **exists and is severe** but not its prevalence. The
legitimate-versus-different reading is from descriptions, not expert labelling; Prp24p against
TDP-43, thermosome against GroEL and dihydrolipoyl dehydrogenase against fumarate reductase are
unambiguous, but `Q9I1F6`'s GntR against LacI case is a naming collision that needs a hand label.

### Unverified

- Whether the Foldseek widening is superfamily aggregation or cluster-size disparity (open question
  O10). The two were not separated by this design.
- The **prevalence** of top-of-distribution impurity. Round 4 shows it exists and is severe on six
  hand-chosen clusters; it does not establish how often it happens across AFDB.
- Whether `Q9I1F6`'s GntR-versus-LacI top 20 is genuine impurity or a naming collision between
  legitimate homologues. Needs a hand label; it is the case that decides whether a relative-identity
  rule can be tuned at all.
- How the AFDB Cluster resource itself renders taxonomic composition. Both
  `afdb-cluster.steineggerlab.workers.dev` and `cluster.foldseek.com` are JavaScript single-page
  applications whose content did not render for inspection.
- Whether the divergent-paralogue pattern generalises beyond the clusters where descriptions were
  inspected by eye.

---

## Fixtures Identified During Scoping

| Accession | Purpose |
|---|---|
| `Q9I1F6` | Primary. 4,628 MMseqs2 members, 913 Foldseek; query present in MMseqs2, absent from Foldseek |
| `Q13148` | TDP-43. 33,211 members, mean 62.4, 87.1% below 70. The **low-ceiling family** case, where the honest answer is "this family is poorly predicted throughout" |
| `P0A6F5` | GroEL. 93,793 members, 8.82 MB, 1.0 s. Size stress case, **and** the clearest cluster-impurity case: 1,535 T-complex protein 1 subunits in the low tail |
| `P69905` | Haemoglobin α. 5,223 members, mean 96.3. High-ceiling family; low tail contains haemoglobin β |
| `P0AA25` | Thioredoxin. 62,441 members; low tail contains KaiB and glutaredoxin-like |
| `P9WN91` | **Mixed-impurity case.** Top 20 mixes dihydrolipoyl dehydrogenase, dimethylglycine oxidase and aminomethyltransferase at 18 to 20% identity with genuine fumarate reductase at 47 to 53%. Also a mode-anchoring inversion case: modal description is MnmG, not fumarate reductase |
| `Q13148` (second role) | **The harmful-recommendation case.** Entire top 20 is yeast Prp24p at 18% identity to TDP-43. The negative control: no Prp24p may appear in a default shortlist |
| `C1C553` | **Uninformative-annotation case.** n = 110, 60.0% uninformative descriptions, modal "Myofilin" at 26.4%, mean pLDDT 51.7. Triggers decision D20's applicability gate. Note `A0A654FPV7` did **not** reproduce its previously reported 68.7% (measured 6.3%) and must not be pinned without re-derivation |
| `P0AA25` (second role) | **False-positive control.** All 20 top members are genuine thioredoxins yet sit at 27 to 39% identity, so a badly calibrated identity cut removes them all. The sharpest test that the screen is not over-eager |
| `P69905` (second role) | **Clean control**, top-20 median identity 86.6%, none below 30%. Also the mode-anchoring inversion demonstration: modal description is haemoglobin beta |
| `P00533` | Foldseek `clusterTotal` of 1. Small-cluster refusal gate |
| `P04637-2` | Isoform handling; resolves to the parent's cluster |
| `Q8WZ42` | Titin. Multi-fragment, 404 on both flags. Refusal gate |
| `P12345` | HTTP 500 on MMseqs2. Retained pending open question O12, which may reclassify it as rate limiting |

Source material to recycle: `~/PycharmProjects/alphafold-db-data-analysis/statistics_page.ipynb`,
cells 51 to 63. Cell 54 is the raincloud figure, cells 57 and 59 the length histograms, cells 61 and
62 the length-versus-pLDDT joint and density plots. Cells 64 to 67 are BigQuery and FASTA export and
are not carried over. **Cell 52 contains a live hardcoded Google API key against a test server**; it
is not carried over, and the key may warrant rotation.

Colour schemes: `~/Downloads/msa-colorschemes-master.zip`. `cleancolors.json` is a flat
`scheme -> residue -> "#RRGGBBAA"` dictionary covering twelve schemes, all flat per-residue lookups,
usable after stripping the trailing alpha. The bundled `.woff2` colour fonts are a browser technique
and are not usable from matplotlib.

---

## Recommended Advisory Review Before `$concept-to-prd`

| Agent | What to review |
|---|---|
| `computational-structural-biologist` | Open question O10, since the notebook's description of what a Foldseek cluster *is* depends on it. Whether caveat C1's three-mechanism account is correctly stated. Whether caveat C9's handling of a 4-point spread is honest enough. Whether the cluster-impurity finding in C5 is a known phenomenon with existing literature we should cite |
| `bioinformatics-data-engineer` | Caveat C4 and open questions O2, O8 and O12; the rate-limiting constraint and retry policy; the refusal-gate table; decision D13's dependency set |
| `evaluation-benchmarking-specialist` | Open question O11 above all: whether a screen with median Cohen's d of −0.33 can carry the safety argument decision D20 rests on. Whether the fixture set covers the claim space. What a regression test asserts given there is no ground truth |
| `molecular-visualization-specialist` | Decision D14, the Mol\* collapse rule in D10 including its divergence-screen exception, decision D19's stacked bar, and open question O6 |

---

## Readiness Judgment For `$concept-to-prd`

**Ready.**

Twenty decisions are recorded with rationale and accepted tradeoffs, five of them revised against
measurement during scoping. Twelve open questions are isolated in section 6. The evidence behind
every revision is in the Measurement Log, including which round superseded which.

Two qualifications, neither blocking.

Open question **O11b** is now the highest-value piece of outstanding domain input: whether
`Q9I1F6`'s GntR-versus-LacI top 20 is genuine impurity or a naming collision between legitimate
homologues. It is the one measured case where neither the identity signal nor the description signal
separates the classes, so it decides whether the decision D20 flag can be calibrated at all
(O11a). A hand label from someone who knows those regulator families resolves it.

Open question **O10** affects how the notebook *describes* a Foldseek cluster, which is PRD-level
prose rather than implementation detail. The finding that Foldseek clusters are wider is confirmed
either way; only the explanation is uncertain.

**Superseded by round 4:** the earlier concern that decision D20's screen was too weak to carry its
safety argument. That concern rested on a statistic measured at the wrong operating point. Impurity
at the top of the distribution is measured, severe, and harmful at rank 1, so the screen is justified
on evidence. What remains is calibration, not existence.

A note for whoever writes the PRD: this capture is unusual in that three of its founding premises
did not survive scoping. The MSA-depth explanation is dead, the length model was twice mis-specified
before landing, and the motivating observation about sequence clusters being wider is empirically
backwards. The notebook is stronger for it, but the PRD should not reach back to the original idea
statement for framing.
