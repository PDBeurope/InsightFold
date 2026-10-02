# PRD: Cluster Quality Diagnostic Notebook

## 0. Document Overview

| Field | Value |
|---|---|
| Product / Notebook Name | `cluster_quality_diagnostic` |
| Version | 0.2 (scope reduced from 0.1; see section 4) |
| Date | 2026-09-21 |
| Owner | Maxim Tsenkov |
| Stakeholders | AFDB team (potential Similar proteins tab integration); InsightFold notebook users |
| Status | Draft |
| Related Artifacts | `prompts/templates/cluster-quality-diagnostic-decision-capture.md`; `notebooks/homodimer_diagnostic.ipynb`; `~/PycharmProjects/alphafold-db-data-analysis/statistics_page.ipynb` cells 51 to 63 |
| Regulatory Scope | RUO / exploratory |
| Intended Lifecycle Outcome | Prototype, promotable to AFDB Similar proteins tab |

---

## 1. Executive Summary

A user holding an AlphaFold model sees a single global pLDDT with no peer context. They cannot tell
whether that number is normal for the protein's family or unusual.

This notebook takes one UniProt accession, retrieves the protein's two AFDB clusterings, and shows
**how pLDDT is distributed across the family**, where the user's protein sits in that distribution,
and what the best and worst-predicted members of the cluster actually look like as structures and in
alignment against the user's protein.

Measurement during scoping established that between-cluster spread dwarfs within-cluster spread:
between-cluster standard deviation 12.46 against a median within-cluster standard deviation of 4.03,
with a median IQR of 4.2 pLDDT points (40-cluster stratified random sample). Families are internally
homogeneous, so the family's overall level is the primary finding and the user's rank within it
matters mainly at the extremes.

Success means a structural biologist can state what quality their family achieves, whether their
model is typical, and what the range looks like at each end.

---

## 2. Background And Context

**Scientific context.** AFDB clusters its predictions two ways. AFDB50/MMseqs2 groups sequences at a
maximum of 50% identity with at least 90% bidirectional sequence overlap against the
representative's longest sequence. AFDB/Foldseek groups structures at an E-value below 0.01 with at
least 90% bidirectional structural overlap. These approximate what protein evolution calls a family
and a homologous superfamily in InterPro's sense, but they are operational proxies rather than
curated InterPro entries.

**Current pain point.** The AFDB Similar proteins tab lists cluster members but does not analyse
them. It shows membership without the quality distribution of that membership, and does not locate
the user's own protein within it.

**Why a notebook first.** This is a prototype for a capability the Similar proteins tab could carry.
Building it as a notebook lets the scientific claims be tested and caveated before anything is
precomputed or exposed on the website.

**Prior art in this repo.** `homodimer_diagnostic` establishes the conventions this notebook
inherits: explicit data contracts, pinned fixtures, plot functions returning bare Figures,
MolViewSpec views, restart-and-run-all validation, and a per-consumer dependency policy.

**A note on this concept's history.** Several founding premises did not survive scoping. An
MSA-depth explanation for within-family quality variation was measured and rejected; the motivating
observation that sequence clustering shows wider quality spread than structure clustering is
empirically backwards on a random sample; and a planned candidate-recommendation feature was cut
(section 4). The decision capture records all of this with evidence.

---

## 3. Goals

| ID | Goal | Rationale |
|---|---|---|
| G-001 | Show how pLDDT is distributed across a protein's family, for both clusterings | The core of the notebook. Between-family variance dominates, so the family's level is where the information is |
| G-002 | Show where the user's protein sits in that distribution | Answers "is my model typical for its family" |
| G-003 | Show the best and worst-predicted members of the cluster as structures and in alignment against the user's protein | Turns the distribution's extremes into something inspectable rather than a number |
| G-004 | Show how sequence length relates to the distribution | Makes visible the main covariate that shapes it |
| G-005 | Run reproducibly on Colab free tier in one restart-and-run-all pass | InsightFold portability requirement; the notebook is intended to be shared |

---

## 4. Non-Goals

**Cut from v0.1 and deferred to a later version:**

- **No ranked candidate table and no recommendation** that the user switch to a different model.
  Measurement showed AFDB clusters can contain functionally divergent proteins near the top of the
  pLDDT distribution, and no reliable screen for this is available without family-level annotation,
  which is itself deferred. Rather than ship a recommendation the notebook cannot stand behind, v1
  shows the distribution and its extremes and leaves the judgement to the reader.
- **No eligibility filters** (length, reviewed status, organism, reference proteome). They existed to
  serve the candidate table.
- **No identity or divergence screening.**
- **No InterPro annotation.**

**Out of scope generally:**

- No new score. Every number traces to a field AFDB publishes.
- No batch input, local structure files, or non-AFDB proteins.
- No cluster-wide account of MSA depth; the proxy for it was measured and rejected.
- No interactive or navigable MSA viewer.
- No claim about experimental tractability.
- No assertion that either clustering method is better than the other.
- Not production. RUO and exploratory.

---

## 5. Target Users

| User Type | Needs | Expertise Assumptions |
|---|---|---|
| Structural biologist / modeller | To know whether their model's confidence is typical for its family, and what the family's range looks like | Comfortable with pLDDT and structural homology; will not read AFDB clustering methodology unprompted |
| AFDB / PDBe internal analyst | Where prediction quality sits across family space, and how the two clusterings differ | Deep familiarity with AFDB data; will scrutinise methodology |
| Bench scientist | Context for whether a protein is well predicted relative to its relatives | Understands protein families; may over-read confidence as a tractability signal and must be guarded against it |

---

## 6. Use Cases / User Stories

| ID | Use Case / Story | Priority |
|---|---|---|
| UC-001 | As a modeller, I want to see the pLDDT distribution of my protein's family, so I know what quality this family achieves | Must |
| UC-002 | As a modeller, I want to see where my protein sits in that distribution, so I know whether it is typical | Must |
| UC-003 | As a modeller, I want to see the best and worst-predicted members as 3D structures, so I can see where the confidence differences sit | Must |
| UC-004 | As a modeller, I want my protein aligned against those two, with percent identity shown, so I can judge how related they actually are | Must |
| UC-005 | As an analyst, I want to see both clusterings and understand why they differ | Should |
| UC-006 | As any user, I want to see how sequence length relates to pLDDT across the cluster | Should |
| UC-007 | As any user, I want the notebook to refuse clearly and name the reason when my protein cannot be analysed | Must |

---

## 7. Functional Requirements

| ID | Requirement | Priority | Notes |
|---|---|---|---|
| FR-001 | Accept one bare UniProt accession as the sole user input, validating and normalising it | Must | The `AF-{ACC}-F1` form must be normalised or rejected with guidance |
| FR-002 | Retrieve both the AFDB50/MMseqs2 and AFDB/Foldseek clusterings automatically, with no user-facing switch | Must | Matches the `homodimer_diagnostic` precedent |
| FR-003 | Show the pLDDT distribution of each clustering, with summary statistics | Must | G-001. The core output |
| FR-004 | Mark the user's protein on the distribution and report its position, categorically as typical or extreme, with the percentile shown but de-emphasised | Must | Position within the body carries little information; see section 12 |
| FR-005 | Present the two clusterings as answering different questions, never as a like-for-like comparison | Must | MMseqs2 is one sequence family; Foldseek aggregates families across a structural superfamily |
| FR-006 | Present the Foldseek distribution as one point per sequence family scored by that family's best member, and draw no marker for the user's protein on it | Must | Foldseek clusters AFDB50 representatives, each the highest-pLDDT member of its family, so the distribution is of per-family maxima. Marking an individual on it would be a category error |
| FR-007 | Show sequence-length distributions and the relationship between length and pLDDT | Should | G-004 |
| ~~FR-008~~ | **REMOVED 2026-09-23.** Taxonomic composition figure. The cluster API returns free-text species names with no lineage, so the figure depended entirely on a secondary UniProt lookup and was dropped rather than kept as a weakly-grounded extra | n/a | |
| FR-009 | Identify the best and worst-predicted members of the cluster | Must | The API already returns members sorted descending by pLDDT |
| FR-010 | Render the user's protein, the best and the worst member in Mol\*, pLDDT-coloured, collapsing to two views when the user's protein is itself an extreme | Must | G-003 |
| FR-010a | Plot the per-residue pLDDT profile for each rendered model | Must | Average pLDDT over a chain is substantially a function of disorder content, so the profile is what distinguishes a uniformly poor model from an ordered domain plus a disordered tail |
| FR-011 | Produce **local** pairwise alignments of the user's protein against the best and against the worst member, collapsing to one when it is an extreme, reporting identity, aligned length and query coverage. These numbers are mandatory; only the rendered figure is optional | Must | G-003, UC-004. Identity alone cannot separate a functionally divergent member from a legitimate distant homologue, so it is always reported with aligned length and coverage |
| FR-012 | Show MSA coverage for the highlighted models, degrading gracefully when the AFDB MSA endpoint is unavailable | Should | The endpoint currently returns 403 |
| FR-013 | Accept a user-supplied MSA for a highlighted model, validated by sequence match rather than filename alone | Could | |
| FR-014 | Refuse clearly, naming the condition, for proteins absent from clustering, singleton or near-singleton clusters, fragmented proteins, and malformed or absent accessions | Must | UC-008 |
| FR-014 (extended) | When the accession does not resolve, explain the likely causes and what to do | Must | Name both "not catalogued in AFDB" and "deleted from UniProtKB", and suggest trying an equivalent or closely related protein. This is the notebook's first possible blocker, so the message must be actionable |
| FR-015 | State the caveats that bound the analysis as prose the reader cannot skip past, including how MMseqs2 cluster membership works | Must | Specifically: the narrow within-family spread; the two clusterings not being comparable; and that members are within threshold of the cluster *representative* rather than of each other, which is why a distantly related member can appear and why percent identity is reported |
| FR-016 | Close with a summary restating the family's level, the user's position, and the caveats | Should | |

---

## 8. Non-Functional Requirements

| ID | Requirement | Acceptance Signal |
|---|---|---|
| NFR-001 | Runs on Google Colab free tier within roughly 60 s of install time | An evidenced Colab run, not an inference. The repo has previously shipped an unverified Colab path |
| NFR-002 | Restart-and-run-all completes with no hidden state | Clean-kernel run on every fixture accession |
| NFR-003 | Respects each AFDB endpoint's own rate policy | Cluster endpoint: at least 2.5 s between requests, since 0.8 s produced spurious 5xx. Prediction endpoint: no spacing needed, measured at 45 back-to-back requests with zero failures. The two must not share one policy |
| NFR-004 | Tolerates large clusters without special handling | 93,793 members return in about 1 s; no cap, pagination or sampling |
| NFR-005 | Plot functions return bare `matplotlib.figure.Figure` objects built without pyplot | Repo convention |
| NFR-006 | No JavaScript widgets | Colab's custom widget manager is a known failure point; static or IFrame rendering only |
| NFR-007 | Dependencies limited to numpy, pandas, matplotlib, seaborn, requests, molviewspec, biopython and pyMSAviz | The repo's pandas prohibition is per-consumer and does not bind this notebook |
| NFR-008 | API base URLs defined once so the forthcoming rate-limited mirror is a one-line change | No endpoint promoted in markdown prose |
| NFR-009 | British spelling; no em dashes in markdown, plot labels or printed output | Repo convention |
| NFR-010 | Notebook outputs stripped before commit | Repo convention |

---

## 9. Data Requirements

| Data Source | Required Fields | Provenance / Licensing | Risk |
|---|---|---|---|
| AFDB cluster members API | `clusterTotal`; `afdbAccessions`, `uniprotDescriptions`, `speciesNames`, `sequenceLength`, `averagePlddt` | AFDB, public but unadvertised; a rate-limited mirror is planned | The only source of member lists, so no fallback. Rate-limits under fast querying. Does not reliably return a cluster containing the query. Returns no sequences |
| AFDB prediction API | `sequence`, `cifUrl`, `globalMetricValue`, `msaUrl`, `latestVersion` | AFDB, public | Called for the user's protein plus the two highlighted members only. No rate spacing required |
| AFDB MSA files | a3m alignments via `msaUrl` | AFDB, public | Currently HTTP 403. The automatic path cannot be validated until the endpoint returns, and must be declared untested |
| User-supplied MSA | a3m | User's own file | Must be validated against the AFDB sequence, not trusted on filename |
| Colour schemes | `cleancolors.json` from the user's `msa-colorschemes` archive | Third-party, check before redistribution | Only the JSON is usable; the bundled colour fonts are browser-specific |

**Caching:** none required. The notebook makes a small number of calls per run.

**Provenance note:** the source notebook contains a live hardcoded Google API key against an AFDB
test server. It is not carried over, and the key should be rotated.

---

## 10. Notebook / UX Expectations

**Entry point.** A single parameter cell holding one UniProt accession. The user edits one line and
runs all.

**Primary outputs**, in order: the family's pLDDT distribution for both clusterings; the user's
position in it; length and taxonomy context; two or three Mol\* views; one or two pairwise
alignments with percent identity; MSA coverage where available.

**Interpretation burden.** Four explanations must land or the notebook misleads. That a percentile
spanning about four pLDDT points is not a meaningful ranking. That a Foldseek point is one family
scored by its best member, not one protein. That cluster members can be distant homologues with
different biological function, which is why the alignments report identity, aligned length and
coverage. And that average pLDDT is substantially a function of disorder content, so a low-mean
family is not necessarily a badly predicted one.

**Visualisation classes required:** distribution plots with the user's protein marked; length
histograms and length-versus-pLDDT joint plots; a taxonomic composition stacked bar split by pLDDT
band; pLDDT-coloured 3D structure views; and paged, colour-coded sequence alignment figures.

**User-editable parameters:** the accession, the alignment colour scheme, the taxonomic rank, and
optionally a path to a user-supplied MSA.

**Educational level.** Assumes familiarity with pLDDT and protein families. Does not assume
knowledge of how AFDB clustering was constructed, which must be explained because it determines how
the results should be read. Uses InterPro's vocabulary for family and homologous superfamily, with
one explicit statement that these clusters approximate rather than instantiate those concepts.

---

## 11. Validation And Success Criteria

| Criterion | Method | Target |
|---|---|---|
| Runs end to end on every fixture accession | Restart-and-run-all on a clean kernel | 100% pass, no manual intervention |
| Runs on Colab free tier | An actual Colab execution, recorded | One successful run within the install budget, explicitly evidenced |
| Refusal gates fire correctly | Run each refusal fixture | Every refusal names its condition in human-readable prose; no bare tracebacks |
| Distribution statistics are correct | Compare against directly computed statistics on pinned cluster snapshots | Exact agreement |
| Percent identity is correct | Compare against an independent alignment on pinned sequences | Agreement to within rounding, with the identity convention stated |
| Extremes are correctly identified | Check against the API's own descending sort on pinned snapshots | Best is the first member, worst is the last |
| Caveats are unavoidable | Review by someone who did not build it | A reader who looks only at the figures does not conclude that the percentile is a ranking, nor that structure clustering predicts quality better |
| Tests do not depend on a live endpoint | Cluster API responses pinned as fixture JSON | Every assertion runs offline against a rate-limited, soon-to-be-replaced service |
| Dependency budget respected | Inspect imports | No package outside NFR-007 |

---

## 12. Risks, Assumptions, And Mitigations

| Item | Type | Impact | Mitigation / Owner |
|---|---|---|---|
| Users read a narrow percentile spread as a meaningful ranking | Risk | Medium. Median p25 to p75 width is only 4.2 pLDDT points | FR-004 reports position categorically and de-emphasises the percentile; the band width is stated |
| The best-scoring cluster member is a distant homologue with different biological function | Risk | Medium in v1, since nothing is recommended. Measured: the top 20 of a TDP-43 cluster are all a yeast Prp24p protein at 18% identity | FR-011's percent identity makes it visible on the same figure; FR-015 caveat states it explicitly. This risk is why the candidate table was cut |
| Users read high pLDDT as experimental tractability | Risk | Medium | Explicit statement alongside the distribution |
| The cluster API returns a cluster not containing the query | Risk | Medium. Statistics could be computed against the wrong population | FR-006 handles absence explicitly and never assumes membership. Accessions sharing an amino acid sequence resolve to the same cluster, which accounts for the cases observed |
| Some cluster members no longer resolve in UniProt | Expected behaviour | Low. AFDB and UniProt drift between releases | FR-014 explains it at the entry point and suggests an equivalent protein; the alignment path degrades with a named message when a member's sequence is unavailable |
| A large cluster spans wider divergence than "50% identity" suggests | Expected behaviour | Low, once explained. Members are within threshold of the cluster representative, not of each other | FR-015 requires the notebook to explain the clustering mechanism before showing the extremes, so a distantly related member reads as expected rather than as an error |
| The MSA endpoint stays unavailable | Risk | Low. Degrades an optional section | FR-012 graceful degradation; FR-013 user-supplied path; the automatic path declared untested |
| Colab remains unverified at release | Risk | Medium. The repo has shipped this mistake before | NFR-001 requires evidenced Colab execution, not inference |
| Average pLDDT is comparable across members of one cluster | Assumption | Foundational | Reasonable: all members come from the same prediction pipeline |
| The cluster API returns the full cluster in one response | Assumption | Verified to 93,793 members | Q-005 tracks the residual `clusterTotal` question |

---

## 13. Lifecycle Readiness

**Expected next step:** `$prd-to-notebook-spec`, then spec review, fixture curation, implementation,
execution validation and notebook review.

**Human and domain review needed before release:**

- Computational structural biology review of how the notebook describes a Foldseek cluster (Q-003),
  and of whether cluster impurity is a known phenomenon with literature to cite rather than a novel
  observation.
- A domain reader who did not build the notebook confirming that the three counter-intuitive
  explanations in section 10 land.

**Graduation path.** Promotion of the family distribution view into the AFDB Similar proteins tab,
where it could be precomputed per entry. A standing-notebook outcome is also credible.

**Evidence needed before graduation:** repeat use by people other than the author, and an evidenced
Colab run.

**Deferred to a later version:** the candidate shortlist with its eligibility filters and identity
screening, which requires family-level annotation to be safe; and InterPro cluster-coherence
annotation.

---

## 14. Open Questions

| ID | Question | Blocking? | Owner / Next Step |
|---|---|---|---|
| Q-001 | RESOLVED. Four accessions return a byte-identical cluster they are not members of | No | Most likely they share the same amino acid sequence and resolve to the same cluster. Expected behaviour. FR-006 handles absence regardless |
| Q-002 | Is the Foldseek widening effect due to superfamily aggregation or to cluster-size disparity? | No. The finding holds either way; only the explanation is uncertain | Affects FR-005 prose. Computational structural biology review |
| Q-003 | Can `clusterTotal` disagree with the returned array length? | No | Always equal so far. If the query is a member the response omits, percentiles are off by one. Spec-level verification |
| Q-004 | Is `P12345`'s HTTP 500 a per-accession fault or rate limiting? | No | Re-test in isolation with generous spacing. Determines whether it remains a useful fixture |
| Q-005 | Which MSA formats does FR-013 accept beyond a3m? | No | Spec-level decision |
| Q-006 | Is the flat `clustal` colour scheme acceptable, or is true Clustal X colouring required? | No | User preference. Real Clustal X requires hand-written column conservation rules |
| Q-007 | How much taxonomic resolution is worth fetching for FR-008? | No | Genus level is viable and measured; species level is infeasible and unsafe |

---

## 15. Spec Handoff Notes

The following belong to `$prd-to-notebook-spec`, not this PRD:

- Notebook section order and cell blueprints.
- Data contracts: exact response schemas, field types, the query-absence resolution procedure,
  accession normalisation rules, and the refusal-gate behaviour table.
- Fixture manifest and expected output snapshots, via `$fixture-selection`. Pinning cluster API
  responses as files is a requirement, not an option, given the endpoint is rate-limited and due to
  be replaced.
- Validation plan, including restart-and-run-all checks and the Colab verification procedure.
- Requirement traceability matrix.
- Task decomposition and module breakdown, including whether logic belongs in `src/insightfold/`.
  Note the repo's D8 rule: a function is promoted when its second consumer appears, not in advance.
- User-facing copy drafts beyond the messaging intent in section 10.
- Plot specifications: figure types, palettes, the taxonomic long-tail aggregation rule, and the
  alignment rendering configuration.
