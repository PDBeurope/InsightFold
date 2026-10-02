# Requirement Traceability Matrix

Regenerated against `requirements.md` v2, after readiness and scientific spec review.

## PRD Goals To Requirements

| PRD Goal | Requirements |
|---|---|
| G-001 Show the pLDDT distribution across the family, both clusterings | REQ-002, REQ-003, REQ-005, REQ-006, REQ-V1, REQ-V2 |
| G-002 Show where the user's protein sits | REQ-004, REQ-E2 |
| G-003 Show best and worst as structures and alignments | REQ-009, REQ-010, REQ-010a, REQ-011, REQ-011a, REQ-017, REQ-017a, REQ-017b, REQ-017c, REQ-017d, REQ-V6, REQ-V7, REQ-V8, REQ-V11, REQ-V12, REQ-E4 |
| G-004 Show how length and taxonomy relate to the distribution | REQ-007, REQ-008, REQ-V4, REQ-V5 |
| G-005 Run reproducibly on Colab in one pass | NFR-001 to NFR-010 |

## Requirements To Sections, Cells, Fixtures And Validation

| Requirement | Goal | Section / Cell | Fixture | Validation Check |
|---|---|---|---|---|
| REQ-001 Accession normalised, not rejected | G-001 | S1, S2 / C005 to C008 | `AF-Q9I1F6-F1`, `P04637-2`, junk | Accession normalisation |
| REQ-001a Isoform substitution stated | G-001 | S2 / C007 | `P04637-2` | Isoform substitution is stated |
| REQ-002 Fetch both clusterings | G-001 | S4 / C010 | all pinned | `clusterTotal` agreement; Rate policy separation |
| REQ-003 Distribution and summary statistics | G-001 | S5 / C012, C013 | `Q13148`, `P69905` | Distribution statistics |
| REQ-004 Query position, categorical, p5/p95 rule | G-002 | S5 / C015 to C017 (inside the cluster-summary section) | `P69905`, `Q13148`, `C1C553` | Query position; Low-n caveat |
| REQ-005 Two clusterings, different questions | G-001 | S3 / C017 then C022, adjacent | `P69905`, `P0A6F5` | Foldseek figure carries no query marker |
| REQ-006 Structure cluster, no query marker | G-001 | S3 / C020 to C022 | `P69905`, `Q9I1F6` | Foldseek figure carries no query marker |
| REQ-007 Length against pLDDT | G-004 | S4 / C019 | `P69905`, `P0A6F5` | Restart-and-run-all |
| REQ-009 Identify best and worst | G-003 | S10 / C026 | all pinned | Extremes identified; Sort assumption guarded |
| REQ-010 Mol\* views with collapse rule | G-003 | S11 / C026, C030 | `SYN-COLLAPSE`, `P69905` | Collapse rule; Viewer budget; Mol\* renders or degrades |
| REQ-010a Per-residue pLDDT profile | G-003 | S11 / C030a | `Q13148` | Per-residue profile renders |
| REQ-011 Local alignment, identity, coverage and aligned range | G-003 | S12 / C027, C031 | FX-02, FX-03, FX-05, FX-12 | Local alignment convention; Distance is visible; Low identity is not overinterpreted; Missing sequence handled |
| REQ-011a Numbers mandatory, figure optional | G-003 | S12 / C032, C033 | `P69905` | Identity and coverage are mandatory |
| REQ-017 TM-align superposition, displayed in Mol\* | G-003 | S12a / C033a, C033b, C033d | `P69905` against `P0C0U8` | Transform orthogonal, determinant +1, reproduces the reported RMSD |
| REQ-017a TM-score primary, RMSD de-emphasised | G-003 | S12a / C033a, C033b | `O15552` | Numbers print outside every figure guard |
| REQ-017b Alignment rendered in scheme and in pLDDT bands, grouped per model | G-003 | S12a / C033d | `O15552` | Both panels produced from one alignment |
| REQ-017c SSE track from `_struct_conf` | G-003 | S12a / C033c, C033d | `O15552`, `Q13148` | Helix, strand and turn codes present; no external tool |
| REQ-017d Mol\* displays, TM-align computes | G-003 | S12a / C033a | none | Stated in prose |
| REQ-V11 Superposition view, query grey and fixed | G-003 | S12a / C033d | `O15552` | Transform applied to the member only |
| REQ-V12 Structural alignment panels | G-003 | S12a / C033d | `O15552` | Titled with TM-score, RMSD, aligned length, identity |
| REQ-012 MSA coverage, graceful degradation | G-003 | S13 / C034, C036 | `SYN-A3M` | MSA 403 handled; MSA path declared untested |
| REQ-013 User a3m validated by sequence | G-003 | S13 / C023b, C035 | `SYN-A3M` wrong-sequence variant | User MSA validated by sequence |
| REQ-014 Refusal gates | all | S4 / C007, C011 | `Q8WZ42`, `P12345`, junk | Refusal: fragmented / malformed / 5xx |
| REQ-014a Degrade per clustering | all | S4, S8 / C021 | `P00533` | Partial degradation |
| REQ-015 Caveats as unskippable prose | all | C012, C015, C020, C025, C037 | `P69905`, `Q13148` | Caveats precede their outputs; Caveats are unavoidable |
| REQ-016 Closing summary | all | S14 / C037 | `P69905` | Restart-and-run-all |
| REQ-017 Actionable failure at the entry point | all | S4 / C011 | junk, known-deleted accession | Refusal: malformed or absent |
| REQ-E1 Clustering definitions, InterPro framing | G-001 | C009 | n/a | Caveats are unavoidable |
| REQ-E2 Narrow spread is not a ranking; use sd and IQR, not ICC | G-002 | C015 | `P69905` | Caveats are unavoidable |
| REQ-E3 Remote homology detection by structure | G-001 | C020 | `P69905` | Foldseek figure carries no query marker |
| REQ-E4 Distant homologue, not "a different protein" | G-003 | C025 | `Q13148`, `Q9I1F6` | Distant homologue not called a different protein |
| REQ-E5 Not tractability; not experimentally validated | G-001 | C012 | n/a | Caveats are unavoidable |
| REQ-E6 British spelling, InterPro vocabulary | all | all markdown | n/a | Prose conventions |
| REQ-E7 Cascaded clustering, transitive linkage | G-003 | C025 | `Q13148`, `Q9I1F6` | Caveats are unavoidable |
| REQ-E8 Average pLDDT is disorder-weighted | G-001 | C012, C025 | `Q13148` | Per-residue profile renders |
| REQ-E9 Cite, do not claim discovery | G-003 | C025 | n/a | Caveats are unavoidable |
| REQ-V1 Distribution raincloud, exactly one | G-001 | C017 | `P69905` | Figures are bare Figures; exactly one raincloud per clustering |
| REQ-V1a Shared fixed x-axis across both rainclouds | G-001 | C017, C022 | `P69905`, `Q13148` | Both figures report identical xlim |
| REQ-V2 Foldseek distribution, labelled, no query line | G-001 | C022 | `P69905` | Foldseek figure carries no query marker |
| REQ-V4 Length versus pLDDT | G-004 | C019 | `P69905`, `P0A6F5` | Figures are bare Figures |
| REQ-V6 Mol\* via `$molviewspec-rendering` | G-003 | C030 | `P69905` | Mol\* renders or degrades; Viewer budget |
| REQ-V7 Alignment figures, built-in colour fallback | G-003 | C033 | `Q13148`, `P69905` | Missing colour file; Identity and coverage are mandatory |
| REQ-V8 Per-residue pLDDT profile | G-003 | C030a | `Q13148` | Per-residue profile renders |
| REQ-V9 MSA coverage plot | G-003 | C036 | `SYN-A3M` | MSA coverage plot is executed |
| NFR-001 Colab budget | G-005 | C002, C003 | `P69905` | Colab execution; `TODO(merge)` resolved |
| NFR-002 Restart-run-all, no hidden state | G-005 | all | all | Restart-and-run-all; Cell-order independence; No hidden state |
| NFR-003 Per-endpoint rate policy | G-005 | C008, C010, C027 | `P69905` | Rate policy separation |
| NFR-004 Large clusters | G-005 | C010 | `P0A6F5` | Large cluster |
| NFR-005 Bare Figures | G-005 | all visualization cells | n/a | Figures are bare Figures |
| NFR-006 No JS widgets | G-005 | C030, C033, C036 | n/a | No JS widgets |
| NFR-007 Dependency set, including transitive tree | G-005 | C003, C004 | n/a | Dependency budget; Prohibited imports absent; Lazy imports |
| NFR-008 Single API base constant | G-005 | module | n/a | source inspection |
| NFR-009 Prose conventions | G-005 | all | n/a | Prose conventions |
| NFR-010 Outputs stripped | G-005 | git | n/a | Outputs stripped |
| NFR-011 Full-width figures, unwrapped prose | G-005 | C004, all markdown | `P69905` | Figures are bare Figures; CSS emitted once |
| NFR-012 Resilient downloads | G-005 | module `_get`, `fetch_model_plddt` | `O15552` (311 KB mmCIF) | Retry on IncompleteRead; truncated mmCIF raises |

## What The Spec Review Closed

The readiness review correctly found that the previous version of the "no requirement is untested"
claim checked only that a validation *row existed*, not that the row named a runnable artifact.
Seven rows failed that stricter test. All seven are now closed:

| Previously failing | Now |
|---|---|
| Collapse-rule branch had no fixture | `SYN-COLLAPSE`, construction specified, T015 |
| Missing-sequence degradation had no fixture | FX-12 `SYN-NOSEQ`, T017, plus a prediction-endpoint failure table. Note T017 found **no real trigger**: AFDB serves a sequence even for UniProt-deleted members |
| MSA coverage plot had no check that renders it | `SYN-A3M`, T016, T071b |
| Taxonomic rank parameter was unreachable | C023b moved before C023c and C024 |
| Colour scheme had no fallback for a missing file | Built-in default, T014, T071d |
| REQ-003 statistics were both asserted and forbidden | Assert-or-monitor line drawn by input source |
| REQ-013 user a3m was synthetic and unspecified | `SYN-A3M` wrong-sequence variant |

The scientific review closed a further set:

| Finding | Resolution |
|---|---|
| Foldseek reference line was a category error | REQ-006, REQ-V2: no query marker; axis states one point per family |
| REQ-E7 described cascaded clustering incorrectly | Rewritten as transitive linkage through intermediate representatives |
| Default fixture showed the alarming case | Default changed to `P69905`; `Q9I1F6` retained as the distant-but-legitimate fixture |
| "Different protein" framing was an overclaim | REQ-E4: distant homologue with possibly different biological function |
| Global alignment was the wrong default | REQ-011: local alignment, reported with aligned length and query coverage |
| ICC was over-loaded and its assumptions violated | Dropped from user-facing prose; sd 12.46 against 4.03 and IQR 4.2 used instead |
| Average pLDDT conflates disorder with quality | REQ-E8 plus REQ-010a per-residue profiles |
| Cluster heterogeneity presented as a discovery | REQ-E9: cite Steinegger and Söding 2018, Barrio-Hernandez et al. 2023 |
| Foldseek circularity unstated | REQ-E3: alignability selection stated; also settles OQ-2 |
| Unobservable acceptance criteria | REQ-004 and REQ-005 rewritten as checks; former REQ-V9 deleted |
| Thresholds unreviewed and unsourced | New "Thresholds And Their Sources" table in `requirements.md` |

## Orphan Check

Every cell in the blueprint carries at least one requirement ID, and every requirement above appears
in at least one cell, at least one task and at least one validation check.

Requirements resting on human judgement (REQ-015 and the REQ-E series) route to named human review
in `validation.md`. The former REQ-V9, "figures cannot contradict prose", was **deleted** as
unobservable rather than left as an untestable requirement.

## Deferred Scope, Traced

| Cut item | Why | Where the reasoning lives |
|---|---|---|
| Ranked candidate table | No reliable screen for cluster heterogeneity without family-level annotation | PRD section 4; decision capture |
| Four eligibility filters | Existed only to serve the candidate table | PRD section 4 |
| Identity and divergence screening | Could not separate FX-02 (26.73%, functionally divergent) from FX-03 (25.09%, legitimate relatives). Identity ranks them the wrong way round, and coverage is parameter-sensitive | Decision capture, round 4 |
| InterPro annotation | Deferred by the owner | PRD section 13 |
