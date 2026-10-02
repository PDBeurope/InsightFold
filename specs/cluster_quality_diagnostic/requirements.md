# Requirements

Source PRD: `prd/cluster_quality_diagnostic_notebook_prd.md` (v0.2)
Revision: v2, after readiness and scientific spec review (2026-09-21). Corrections are marked.

## User Story

As a structural biologist, AFDB analyst or bench scientist holding one AlphaFold model, I want to see
how pLDDT is distributed across my protein's family, where my model sits in that distribution, and
what the best and worst-predicted members look like, so that I understand whether my model's
confidence is typical and what the family's range actually is.

### Workflow stories

| ID | Story |
|---|---|
| WS-1 | I paste one UniProt accession, run all, and get an answer without editing anything else |
| WS-2 | I see the pLDDT distribution of my protein's family |
| WS-3 | I see where my protein sits, without being misled by a spread too narrow to rank on |
| WS-4 | I see the best and worst-predicted members as 3D structures, with their per-residue confidence |
| WS-5 | I see my protein aligned against each of them, with identity and coverage, so I can tell how distant they are |
| WS-6 | I see how sequence length relates to the distribution |
| WS-7 | When my protein cannot be analysed, I am told why in plain language |

## Functional Requirements

| ID | Requirement | Acceptance Criteria | Priority |
|---|---|---|---|
| REQ-001 | Accept one bare UniProt accession as the sole user input; **normalise** the `AF-{ACC}-F1` form rather than rejecting it, printing what was normalised | `AF-Q9I1F6-F1` becomes `Q9I1F6` and the notebook prints "interpreted as Q9I1F6". Malformed input rejected before any network call. `P04637-2` accepted | Must |
| REQ-001a | When an isoform accession resolves to a parent cluster, print which sequence is used for alignment and which structure is rendered | For `P04637-2`, output names whether the isoform or canonical sequence and structure were used | Must |
| REQ-002 | Fetch both clusterings automatically, no user-facing switch | Two responses per run, at least 2.5 s apart. `clusterTotal` and array length both recorded | Must |
| REQ-003 | Show the pLDDT distribution of the **sequence** cluster with summary statistics: mean, median, standard deviation, band composition, p25 to p75 width. Printed labels are spelled out, matching the glossary in the markdown above them, and the prose carries **no constants derived from any sample other than the user's own cluster** | Against the pinned `Q13148` snapshot, the module returns mean 62.4 and 87.1% below 70 exactly. **Prose must not call this a "low-quality family"**; see REQ-E8 | Must |
| REQ-004 | Mark the user's protein on the distribution and classify its position against **that cluster's own distribution**: typical if between p5 and p95, otherwise high or low extreme. Presented **inside the cluster-summary section**, not as a section of its own | Output gives the protein's own mean pLDDT, then the descriptor, then the percentile, then the cluster's p25 to p75 width. The descriptor thresholds are stated in the notebook. When the width is narrow the output says so explicitly. When the protein is absent from the cluster it says it has no rank, without error | Must |
| REQ-005 | Present the two clusterings as answering different questions | No single `Axes` contains both distributions. The two rainclouds sit adjacent in the cluster-summary section so they can be compared, with the explanation of what a structure-cluster member is between them | Must |
| REQ-006 | Present the structure cluster as structurally similar proteins, grouped by shape rather than sequence, and state that the user's protein is usually not shown because the clustering is built from cluster representatives. Do **not** draw the user's pLDDT as a reference line. Keep the explanation short: it is orientation for a reader, not a methods note | Prose says these are structurally similar proteins and that the query is usually absent, in under four sentences. No query reference line appears | Must |
| REQ-007 | Show the relationship between sequence length and pLDDT | A joint density with marginal distributions on both axes. Standalone length histograms were removed; the marginals carry the same information in context | Should |
| REQ-009 | Identify the best and worst-predicted members of the sequence cluster | Best is index 0, worst index −1 of the descending-sorted array, **asserted not assumed** against pinned snapshots | Must |
| REQ-010 | Render the user's protein, the best and the worst member in Mol\*, coloured by **AFDB pLDDT confidence band** using the same palette as REQ-V1, collapsing to two views when the user's protein is itself an extreme | Three views normally, two when collapsed. At most three live viewers. Colour is applied as repeated scoped `.color(color=..., selector=[...])` calls on a **single** cartoon representation, grouped by band. **Never** as one extra representation per run over a base representation of the same atoms: that renders overlapping cartoon geometry which z-fights, and the model appears uniformly grey. A legend of the four bands is shown above the views | Must |
| REQ-010b | **NEW.** Let the user choose which clustering the rendered models are drawn from. **The shipped default is `structure`**; REQ-010b originally specified `sequence` and the notebook prose still claimed it, which was corrected 2026-10-02 by removing the claim rather than changing the value. Either is defensible and the choice is the user's | `MODEL_SOURCE_CLUSTER` accepts "sequence" or "structure". Choosing "structure" states that the models are family representatives rather than individual proteins, and falls back to the sequence family with a message when no usable structure cluster exists | Should |
| REQ-010a | **NEW.** Plot the per-residue pLDDT profile for each rendered model | One profile per view target, from the already-downloaded mmCIF. Makes visible whether a low average is a uniformly poor model or an ordered domain plus a disordered tail | Must |
| REQ-011 | **CORRECTED.** Align the user's protein against the best and worst member using **local** alignment, reporting identity, aligned length and query coverage. Collapse to one alignment when the user's protein is an extreme | Output reads "X% identity over N aligned columns, covering Y% of the query". Convention stated in the notebook. Matches an independent local alignment on pinned sequences | Must |
| REQ-011a | **NEW.** The identity, aligned length and coverage figures are **mandatory and must be printed as text**. Only the rendered alignment figure is optional | With pyMSAviz unavailable, the three numbers still print and the run continues. With the numbers unavailable, the run raises | Must |
| REQ-017 | **NEW.** Superpose the user's protein on each rendered extreme with **TM-align**, and display each superposition in Mol\*, the query held fixed and the member transformed onto it | `tmtools.tm_align` supplies the rotation, translation, RMSD, TM-score and its own structure-based alignment. The transform returned by the library maps chain 1 onto chain 2, so the **inverse** (`uᵀ`, `-uᵀt`) is applied to the member instead, leaving the query in its deposited frame. The matrix is handed to MolViewSpec `.transform()` **column-major**; row-major silently renders a valid-looking but wrongly oriented model. The emitted matrix is orthogonal with determinant +1 |
| REQ-017a | **NEW.** Report the **TM-score** as the primary number, with RMSD secondary and explicitly de-emphasised | TM-score is normalised to the query length, is scale-free, and carries published interpretation points: ≥ 0.5 implies a shared fold, < 0.3 is the expectation for an unrelated pair (Zhang and Skolnick 2004, 2005). RMSD has no such scale and TM-align maximises TM-score rather than minimising RMSD, so a large RMSD over a good alignment is ordinary. The notebook must not present a high RMSD alone as failure |
| REQ-017b | **NEW.** Render TM-align's **own** alignment twice, once in the residue colour scheme and once with every residue in the confidence band of its own model, **grouped per model**: each superposition is immediately followed by the two panels for that same pair, rather than all superpositions followed by all panels | Both panels show the same structure-derived columns, so a column asserts spatial equivalence rather than sequence similarity. The banded panel is the one that answers this notebook's question, because it shows whether a structural match rests on well-predicted residues on both sides. Implemented through pyMSAviz `set_custom_color_func`, which falls through to the scheme when it returns `None` |
| REQ-017c | **NEW.** Annotate each aligned sequence with an **SSE** (secondary structure element) track read from the model's own `_struct_conf` records. The row is labelled `SSE`, and the abbreviation is expanded in the markdown cell above the figures, since the label itself has no room for it | `H` helix (any subtype), `E` strand, `T` turn or bend, blank where unassigned or gapped. No external secondary structure tool: AFDB mmCIF files carry DSSP-style assignments already. Track colours sit outside the pLDDT ramp so a structure row is never read as a confidence row |
| REQ-017d | **NEW.** State plainly that the superposition is computed by TM-align and only **displayed** by Mol\* | Mol\* computes no alignment here. The notebook must not imply that it does |
| REQ-012 | Show MSA coverage for the highlighted models, degrading gracefully | On HTTP 403 the section prints a named message and continues. Detection by status and content type, never by parsing the body | Should |
| REQ-013 | Accept a user-supplied a3m, validated by gap-stripped sequence match | A correctly named file with the wrong sequence is rejected. The `AF-{ACC}-F1` header convention is a hint, never the validation | Could |
| REQ-018 | **NEW.** Every user input is a **Colab form field**, not an editable line of code | `#@param` on the assignment, `#@title ... { display-mode: "form" }` on the cell, so Colab shows widgets and hides the source. A reader changes a value without touching code, and cannot break the cell by editing around the literal. The annotations are plain comments, so local Jupyter and PyCharm runs are unaffected. Validated by parsing every `#@param`: the declared widget type must match the assigned literal, and a dropdown's option list must contain the assigned value |
| REQ-018a | **NEW.** A dropdown's option list must stay in step with the module constant it mirrors | Colab requires the options to be a literal inside the comment, so the alignment-scheme list is necessarily a second copy of `ALIGNMENT_SCHEME_CHOICES`. `resolve_colour_scheme()` therefore **names** an unrecognised scheme and says what the valid ones are, instead of the previous silent fallback to Clustal. Checked by comparing the notebook's literal against the module tuple |
| REQ-014 | Refuse clearly, naming the condition, for every unsupported input | See the edge-case table below. No bare tracebacks | Must |
| REQ-014a | **NEW.** Degrade **per clustering**, not per run. Refuse the whole run only when the **sequence** cluster is unusable | `P00533` has Foldseek `clusterTotal` 1 but a usable MMseqs2 cluster: the sequence analysis completes and the Foldseek section is skipped with a named message | Must |
| REQ-015 | Present the bounding caveats as prose the reader cannot skip past | Each caveat appears as a markdown cell before the output it bounds | Must |
| REQ-016 | Close with a summary restating the family's level, the user's position, the identity and coverage against each extreme, and the bounding caveats | Final cell reflects the run's actual values | Should |
| REQ-017 | When the accession does not resolve, explain the likely causes and what to do | Message names both "not catalogued in AFDB" and "deleted from UniProtKB", and suggests trying an equivalent protein | Must |

## Non-Functional Requirements

| ID | Requirement | Acceptance Criteria |
|---|---|---|
| NFR-001 | Runs on Colab free tier within roughly 60 s of install time | An evidenced Colab run, not an inference |
| NFR-002 | Restart-and-run-all completes with no hidden state | Clean-kernel run passes on every fixture |
| NFR-003 | Respects each AFDB endpoint's own rate policy | Cluster endpoint at least 2.5 s between requests; prediction endpoint needs no spacing. The two must not share one policy |
| NFR-004 | Tolerates large clusters without special handling | `P0A6F5` at 93,793 members completes. No cap, pagination or sampling |
| NFR-005 | Plot functions return bare `matplotlib.figure.Figure` objects built without pyplot | No pyplot registry leak; nothing double-displayed |
| NFR-006 | No JavaScript widgets | Alignments and MSAs render statically. Mol\* is embedded with an **`srcdoc` iframe, not a `data:` URI**: PyCharm's renderer treats a `data:text/html` iframe as a navigation to a separate document and hands it to the system browser, opening one Chrome window per viewer on top of rendering inline. `srcdoc` carries the document in an attribute so there is nothing to delegate. This is a deliberate divergence from `complex_interface_utils.py`, which uses the data-URI pattern |
| NFR-007 | **CANONICAL dependency list.** Allowed: numpy, pandas, matplotlib, seaborn, requests, molviewspec, biopython, pyMSAviz, **tmtools** (a thin wheel around the original TM-align C++; the alternative was a 91 MB headless Mol\* npm tree to reach the same algorithm). Prohibited: torch, gemmi, plotly, scipy, ptitprince, google-cloud-bigquery, networkx | No other third-party import, including in pyMSAviz's resolved transitive tree |
| NFR-013 | **NEW.** The notebook holds presentation and parameters, not implementation. Every cell below the bootstrap is a single call into `cluster_quality_utils`; no loops, no `try`/`except`, no formatting | Grep **with comments stripped** (a prose comment may legitimately contain the word "for"): no `for `, `while `, `try:`, `except` or `print(` in any cell except the bootstrap. Measured 2026-10-02: 29 code lines outside the bootstrap, down from about 200. The orchestration layer is `report_*` (prints, mandatory) and `show_*` (renders, optional), a split that keeps mandatory numbers out of figure guards per REQ-011a. Verified behaviour-preserving by byte-comparing printed output and the display sequence before and after the move |
| NFR-008 | Each API base URL defined once | **The forthcoming mirror arrived 2026-10-02** and the switch to `api/workbench/cluster-family` was indeed a one-line change, which is what this requirement existed to buy. That route is advertised, so the former bar on naming the cluster endpoint in prose is lifted; the prediction endpoint is unaffected |
| NFR-009 | British spelling; no em dashes in markdown, plot labels or printed output | Grep clean for em dash, `&mdash;` and `&#8212;` |
| NFR-010 | Notebook outputs stripped before commit | No `outputs` payload in the committed `.ipynb` |
| NFR-012 | **NEW.** Downloads survive a slow or unstable connection. Every request retries transient failures, and every streamed file is verified complete before use | `IncompleteRead`, socket timeouts, connection resets, HTTP 429 and 5xx retry with exponential backoff and jitter, 4 attempts. HTTP 4xx other than 429 fails immediately. The cluster endpoint's 2.5 s spacing is applied before **every** attempt, retries included. File downloads use a 180 s timeout, JSON APIs 30 s. mmCIF completeness is checked against the prediction's residue count and a short file triggers a fresh download, up to 3 times. Failures name the file and the cause, never a bare `IncompleteRead(163840 bytes read)` |
| NFR-011 | **NEW.** Figures fill the width of the notebook column and reflow with the window; markdown prose is not hard-wrapped | `apply_notebook_style()` injects CSS scoped with `:not(.cq-panel-img)`, so the Mol\* and PAE flex panels keep their own layout. Figures render at 144 dpi so they stay crisp when scaled up. Markdown cells hold one line per paragraph and let the renderer wrap |

## Non-Goals

**Cut from an earlier scope and deferred:** ranked candidate table; eligibility filters; identity or
divergence screening; InterPro annotation.

**Out of scope generally:** no new score; no batch input, local files or non-AFDB proteins; no
cluster-wide account of MSA depth; no interactive MSA viewer; no tractability claim; no assertion
that either clustering method is better.

## Thresholds And Their Sources

Every threshold the notebook uses, stated once so none is silent.

| Threshold | Value | Source |
|---|---|---|
| pLDDT bands | 50, 70, 90 | AFDB's own confidence bands; Tunyasuvunakool et al. 2021, Nature 596:590 |
| Typical versus extreme | outside p5 or p95 of **that cluster's** distribution | Chosen for this notebook. Relative to the cluster, never a global cutoff |
| Minimum usable cluster size | 3 members to run at all; below 30, print a low-n caveat because a percentile and an IQR are unstable | Chosen for this notebook |
| Comparable-length note | members beyond ±20% of cluster median flagged in the summary | Median length CV 8.8% measured across 31 clusters |
| Cluster request spacing | 2.5 s | Measured: 0.8 s produced spurious 5xx, 2.5 s produced zero |

## Edge Cases And Failure Behaviour

| Condition | Example | Required behaviour | Requirement |
|---|---|---|---|
| Accession in AF form | `AF-Q9I1F6-F1` | **Normalise**, printing what was interpreted | REQ-001 |
| Isoform | `P04637-2` | Supported; resolves to the parent cluster. Print which sequence and structure are used | REQ-001, REQ-001a |
| Absent or malformed accession | HTTP 404 | Report, quoting the service message, and name both likely causes with a suggestion to try an equivalent protein | REQ-014, REQ-017 |
| Multi-fragment protein | `Q8WZ42` | 404 on both flags. Explain that fragmented proteins are absent from clustering | REQ-014 |
| Sequence cluster too small | fewer than 3 members | Refuse the run, naming the size | REQ-014 |
| **Foldseek cluster too small, sequence cluster fine** | `P00533` | **Run the sequence analysis; skip the Foldseek section with a named message** | REQ-014a |
| Cluster between 3 and 30 members | `C1C553` at 110 is above; smaller clusters exist | Run, with a named low-n caveat on the percentile and IQR | REQ-004 |
| HTTP 5xx | `P12345` | Retry with backoff before reporting | REQ-014, NFR-003 |
| Query absent from returned cluster | Foldseek, always | Normal path. Do not assume membership | REQ-006 |
| `msaUrl` returns 403 | all accessions currently | Named message; section continues | REQ-012 |
| Member sequence unavailable | deleted UniProt entry | Alignment degrades with a named message; the 3D view may still render | REQ-011 |
| `cleancolors.json` absent | Colab, or licence unresolved | Fall back to the module's built-in default scheme. Never block | REQ-V7 |

| Connection cut mid-download | slow or remote link | Retry with backoff. If every attempt fails, name the file and say re-running will try again | NFR-012 |
| **Structure file arrives truncated** | measured: 160 KiB of a 311 KB mmCIF for `AF-O15552-F1` | **Must raise.** A short mmCIF parses happily into a partial structure: the measured case yielded 150 residues of 330 with no error. Completeness is checked against the prediction's residue count | NFR-012 |

## Educational And Documentation Requirements

| ID | Requirement |
|---|---|
| REQ-E1 | Explain both clustering definitions, quoting the AFDB FAQ criteria, and relate them to InterPro's family and homologous superfamily concepts, stating once that these clusters approximate rather than instantiate those concepts |
| REQ-E2 | Explain, in the cluster-summary section and before the percentile is shown, that a rank is only worth what the cluster's own p25 to p75 width makes it worth. **Derive this from the user's own cluster**, never from constants measured on some other sample. **Use the directly interpretable numbers**: between-cluster sd 12.46 against median within-cluster sd 4.03, and median IQR 4.2 points. Do not quote ICC in user-facing prose |
| REQ-E3 | Explain in plain language that these are proteins grouped by shape rather than sequence, and that the user's protein is usually not among them because the structure clustering is built from cluster representatives. Keep it to a few sentences; this is orientation, not methodology. **Measured context, for maintainers rather than the notebook**: the two member lists share no AFDB entry and no exact UniProt accession across 8 accessions, though the same protein can appear in both as different isoforms (18 base accessions shared in `Q13148`, including `P04637` itself for `P04637-2`). The API flags no representative, so which member that is cannot be confirmed from the response |
| REQ-E4 | **CORRECTED.** Explain that a cluster's members can be **distant homologues sharing a fold or domain family but differing in biological function**. TDP-43 and Prp24p are both RRM-domain proteins, so the honest claim is functional divergence, not "a different protein". Identity measures how distant, not whether the function is shared. Direct the reader to identity, coverage, `uniprotDescriptions` and `speciesNames` together |
| REQ-E5 | State that pLDDT is not an experimental tractability predictor, that pLDDT differences between cluster members have not been validated against experimental structures for these families, and that all values move with AFDB releases |
| REQ-E6 | Use British spelling and InterPro's vocabulary throughout |
| REQ-E7 | **CORRECTED.** Explain MMseqs2 cluster membership accurately: AFDB50 uses **cascaded** clustering, so members are linked to the final representative **transitively, through a chain of intermediate representatives**, each link within threshold. Identity to the final representative, and between arbitrary members, is **unbounded below**. Cite Steinegger and Söding 2018, Nat Commun 9:2542. **Do not write that every member is within threshold of the representative; that is false for cascaded clustering** |
| REQ-E8 | **NEW.** State that average pLDDT over a chain is substantially a function of **disorder content**, so differences between members partly measure intrinsic disorder rather than prediction quality. Do not describe a low-mean family as "low quality". Cite Tunyasuvunakool et al. 2021 and Akdel et al. 2022 |
| REQ-E9 | **NEW.** Cite cluster heterogeneity as a known phenomenon rather than presenting it as a discovery: Steinegger and Söding 2018 for cascaded and transitive clustering, and Barrio-Hernandez et al. 2023, Nature 621:637 for AFDB structural clustering merging sequence clusters across remote homology |

## Visualisation Requirements

| ID | Requirement |
|---|---|
| REQ-V1 | Sequence-cluster pLDDT distribution as **one** raincloud: half violin, jittered points coloured by pLDDT band, boxplot, with the user's protein marked. **Exactly one per clustering**; no unmarked duplicate of the same data. Reimplemented in matplotlib and seaborn, not ptitprince. Points subsampled above 2,000 with the subsampling stated in the title; the violin and box always use every member |
| REQ-V1b | **NEW.** Band legend order follows the spatial arrangement of the data it labels. **Horizontal** legends run ascending, `<50` first and `>90` last, reading left to right in increasing confidence: the rainclouds, the 3D band legend, the taxonomy stack. The **vertical** legend on the length-versus-pLDDT density runs descending, `>90` at the top, so that reading down the legend tracks reading down the pLDDT y axis beside it. Band names in HTML legends must be escaped: `<50` is a literal `<` that a browser reads as the start of a tag, silently swallowing the label. Applies to the raincloud, the joint density, the 3D band legend and the taxonomy figure |
| REQ-V1a | **NEW.** Both rainclouds share one fixed x-axis, computed from the pooled minimum and maximum across both clusterings, padded by 20 pLDDT points, rounded outward to the nearest 10 (half-up, so 45 gives 50), and clamped to 0 to 100. Neither figure may autoscale: two distributions drawn one above the other on different axes look alike when they are not. The padding always exceeds half a step, so the rounded bound cannot clip a data point |
| REQ-V2 | **CORRECTED.** Structure-cluster raincloud, same form as REQ-V1, axis labelled "one point = one AFDB50 family, scored by its highest-pLDDT member". **No query reference line.** Drawn on its own Axes, never beside REQ-V1 |
| REQ-V3 | **REMOVED.** Standalone length histograms are dropped; the marginal distributions of REQ-V4 carry the same information in context |
| REQ-V4 | Length-versus-pLDDT **joint density** with marginals on both axes, built on a bare Figure with gridspec rather than `seaborn.JointGrid`, which creates its own pyplot figure. Points are coloured by pLDDT band, matching REQ-V1 and REQ-010, and the pLDDT marginal is banded to match. The 2D histogram and its colourbar are **omitted when colouring by band**, since a density under coloured points is invisible and its colourbar is then dead space; KDE contours carry the density instead. `colour_by_band=False` restores the plain `mako` density. Points subsample above 6,000 |
| REQ-V6 | Mol\* pLDDT-coloured structure views via `$molviewspec-rendering`. No alternative 3D viewer library |
| REQ-V7 | Pairwise alignment figures via pyMSAviz, colour scheme from `cleancolors.json` **with a built-in default when that file is absent**. Identity, aligned length and coverage in the title |
| REQ-V8 | Per-residue pLDDT profile per rendered model, from the downloaded mmCIF |
| REQ-V9 | MSA coverage plot, ColabFold `plot_msa_v2` style, reimplemented in matplotlib |
| REQ-V10 | **NEW.** Each rendered model's **predicted aligned error** matrix, plotted beside its Mol\* view in one flex row. `Greens_r`, `aspect='equal'`, colourbar labelled in Angstrom. Reuses `parse_pae` from `complex_interface_utils`, passing `fallback_lengths` because monomer PAE documents omit the `chains` array |
| REQ-V11 | **NEW.** Structural superposition view: two cartoon structures in one Mol\* instance, the query in grey and held fixed, the member in green and transformed. Colour here encodes **identity of chain**, not confidence, so the palette is deliberately disjoint from the pLDDT ramp and a key states which is which |
| REQ-V12 | **NEW.** The structural alignment panels of REQ-017b and REQ-017c, wrapped at 80 columns, titled with the TM-score, RMSD, aligned length and identity |

*(The former REQ-V9, "figures cannot contradict prose", is removed as unobservable. It is covered by
the human reader review, task T093.)*

## Open Questions

| ID | Question | Status |
|---|---|---|
| OQ-1 | RESOLVED. Four accessions return a byte-identical cluster they are not in | Most likely shared amino acid sequence resolving to the same cluster. REQ-006 handles absence regardless |
| OQ-2 | **RESOLVED.** Why are Foldseek clusters lower-scoring? | Representative selection is a *positive* selection on pLDDT, so a **lower** Foldseek mean can only come from aggregating lower-scoring families. This is a stronger argument than the spread comparison and settles the question in favour of superfamily aggregation |
| OQ-3 | Can `clusterTotal` disagree with the returned array length? | **non-blocking**: REQ-002 records both, so it surfaces in output |
| OQ-4 | Is `P12345`'s 500 a per-accession fault or rate limiting? | **non-blocking**: affects fixture classification only |
| OQ-5 | Which MSA formats beyond a3m does REQ-013 accept? | **non-blocking**: a3m only in v1 |
| OQ-6 | Flat `clustal` scheme or true Clustal X colouring? | **non-blocking**: flat schemes ship |
| OQ-7 | RESOLVED, 2026-09-23. Taxonomic composition figure | **Removed.** The cluster API returns free-text species names with no lineage, so the figure depended on a bolt-on UniProt lookup to say anything at all. Dropped on the owner's call rather than kept as a weakly-grounded extra |
| OQ-8 | **NEW.** Is `cleancolors.json` redistributable? | **non-blocking**: the built-in default scheme makes this an enhancement, not a dependency |
