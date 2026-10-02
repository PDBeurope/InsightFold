# Jupyter Cell Blueprint

Kernel assumption: Python 3.11, fresh kernel, cells executed strictly top to bottom exactly once. No
cell may depend on having been run twice, and no cell may depend on a later cell.

> **Cell IDs are spec-side only.** They do **not** appear in the notebook and must not be
> reintroduced there: they mean nothing to a reader and nothing in the notebook explains them.
> Identify a cell by its purpose and position when tracing a requirement. The same applies to
> REQ, NFR, FX and OQ identifiers: they live in this spec pack, never in notebook comments or
> markdown.

| Cell ID | Type | Purpose | Inputs | Outputs | User Editable? | Hidden-State Risk | Validation Hook | Requirements |
|---|---|---|---|---|---|---|---|---|
**Every cell below the bootstrap is a single call.** The "Contents" column describes what the
cell *does*, which is unchanged; the logic itself lives in the orchestration layer of
`cluster_quality_utils.py` (see `notebook-design.md`). Two merges followed from this:

- **C002 and C003 folded into C001.** The optional-package install, the inline backend and the
  style are all reached through `insightfold.notebook_setup.setup()`, which the bootstrap cell
  ends by calling. The module owns its own optional dependency list rather than restating it in
  a cell, and the bootstrap is now a Colab form that collapses to a title bar.
- **C010 folded into C009.** Both refusal gates now run inside `load_clusters()`, beside the fetch
  they guard, so no cell can read a cluster that has not been gated.
- **C005, C006 and C007 folded into one Colab form cell.** The accession textbox, the
  normalisation and the identity card belong together: in form display-mode Colab shows the
  widget and the printed identity beneath it, which is the confirm-before-you-continue step
  REQ-001 asks for.

**The two parameter cells are Colab forms** (REQ-018), so the reader sets values through widgets
and never edits code. `#@param` and `#@title` are ordinary comments outside Colab, so local runs
are byte-identical; this was verified by running the notebook before and after and comparing the
printed output and the display sequence.

The three cells that stay as visible code are the interface, not implementation: the bootstrap
(it runs before the module is importable and carries the `TODO(merge)` pin), the `ACCESSION`
line, and the display-preferences cell.

| C001 | markdown | Title, purpose, what the user will get | none | none | no | none | prose review | REQ-E1 |
| C002 | code | Bootstrap: `find_repo_root()`, clone on Colab, `sys.path.insert`. Never `pip install` the package | none | `REPO_ROOT` | no | **`sys.path` mutation.** Guard against duplicate inserts | import succeeds | NFR-001, NFR-007 |
| C003 | code | `pip install molviewspec biopython pymsaviz`, guarded by an import probe | none | none | no | Colab-only path; unexercised locally | install under budget | NFR-001, NFR-007 |
| C004 | code | Imports, `%matplotlib inline`, `apply_notebook_style()` | C002 | module symbols | no | **Global matplotlib state.** The dpi change and the width CSS live in the explicit call, never at module import | no pyplot registry leak; CSS emitted exactly once | NFR-005, NFR-011 |
| C005 | markdown | Callout naming the input and giving example accessions. **No longer an "EDIT ME" banner**: the input is a Colab form widget (REQ-018), so the instruction lives on the form itself and the markdown would only duplicate it | none | none | no | none | prose review | REQ-001 |
| C006 | code | `ACCESSION = "P69905"` | none | `ACCESSION` | **yes** | none | value is a bare accession | REQ-001 |
| C007 | code | **Normalise** the AF form rather than rejecting it, printing what was interpreted; reject malformed input. For an isoform, record which sequence and structure will be used | C006 | `ACCESSION` normalised | no | none | `AF-Q9I1F6-F1` becomes `Q9I1F6` with the interpretation printed; junk rejected with a named error; `P04637-2` states its substitution | REQ-001, REQ-001a, REQ-014 |
| C008 | code | Fetch query metadata and sequence from the prediction endpoint (no spacing) | C007 | `query_meta`, `query_sequence` | no | none | **First visible confirmation.** Prints description, organism, length, pLDDT, `latestVersion`, run timestamp | REQ-001, NFR-003 |
| C009 | markdown | What AFDB50/MMseqs2 and AFDB/Foldseek are; the InterPro framing; approximate not instantiate | none | none | no | none | quotes the AFDB FAQ criteria | REQ-E1 |
| C010 | code | Fetch both clusterings, 2.5 s apart; print member counts and a patience note | C007 | `mm_cluster`, `fs_cluster` | no | **Rate limiting.** Spacing enforced inside the fetch helper, not by the caller | counts printed; `clusterTotal` compared to array length | REQ-002, NFR-003, NFR-004 |
| C011 | code | Refusal gates: fewer than three members, fragmented, 404, 5xx with backoff. A 404 must explain likely causes (not catalogued in AFDB, or deleted from UniProtKB) and suggest trying an equivalent protein | C010 | none, or raises | no | none | each refusal fixture produces its named message; the 404 message is actionable, not a bare restatement | REQ-014, REQ-017 |
| C012 | markdown | **Cluster summary**: a term-by-term glossary of what C013 prints and what to read from each number. Every abbreviation expanded except pLDDT. Also that pLDDT is not tractability, and that a chain average is substantially a measure of disorder content so a low mean is not necessarily poor prediction | none | none | no | none | prose review; **no sample-derived constants in the prose**, since the notebook must read sensibly for any clustered protein | REQ-003, REQ-E5, REQ-E8 |
| C013 | code | Compute `family_summary`: mean, median, sd, band composition, p25 to p75 width | C010 | `family_summary` | no | none | matches directly computed values on pinned JSON | REQ-003 |
| C015 | markdown | **Your protein within the cluster.** Defines the descriptor (low extreme / typical / high extreme) against this cluster's own p5 and p95, and explains that a rank is only worth what the p25 to p75 width printed above makes it worth. Sits inside section 3; there is no separate section for it | none | none | no | none | appears **before** C016 | REQ-E2, REQ-004 |
| C016 | code | Locate the user's protein: percentile, category, band width. Handle absence | C010, C008 | `position` | no | **Query may be absent** from the returned cluster. Must not assume membership | percentile and p5/p95 category correct against pinned JSON; low-n caveat fires below 30 members | REQ-004 |
| C017 | visualization | **The single** sequence-cluster raincloud, with the user's protein marked. There is deliberately no second unmarked copy: two near-identical plots of the same data taught nothing | C010, C016 | `Figure` | no | none | seaborn and matplotlib, not ptitprince; bare Figure; exactly one raincloud per clustering | REQ-V1, REQ-003, REQ-004 |
| C019 | visualization | Length-versus-pLDDT **joint density with marginals**. The only length figure: the standalone histograms were removed, since the marginals carry the same information in context | C010 | `Figure` | no | Optional section; must fail soft | bare Figure or named skip; built with gridspec, not `seaborn.JointGrid` | REQ-V4, REQ-007 |
| C020 | markdown | **Structurally similar proteins**, before the figure. Structure is more conserved than sequence, so Foldseek detects remote homologues whose sequence identity has fallen below the threshold for reliable inference. The user's protein is usually not shown, because the clustering is built from cluster representatives. Four sentences, not a methods note | none | none | no | none | appears **before** C021 | REQ-E3, REQ-006 |
| C021 | code | Resolve the AFDB50 representative; skip the section with a named message if the Foldseek cluster is unusable while the sequence cluster is fine | C010 | `representative` | no | **Query reliably absent from Foldseek**, and that is expected, not an error. Must not compute a query percentile here | representative named; `P00533` skips this section but completes the run | REQ-006, REQ-014a |
| C022 | visualization | Structure-cluster raincloud, placed **immediately after C017** so the two distributions can be compared. **No query reference line**: the query is usually not a member | C021 | `Figure` | no | The query is usually not a member of this cluster, so a marker for it would be misleading | bare Figure; no query marker present | REQ-V2, REQ-006 |
| C023a | markdown | Optional-settings callout for the display preferences form | none | none | no | none | prose review | REQ-010b, REQ-013 |
| C023b | code | `ALIGNMENT_COLOUR_SCHEME`, `USER_MSA_PATH`, `MODEL_SOURCE_CLUSTER`, `SHOW_PAE` | none | display settings | **yes** | **Must precede every cell that reads these**: C026, C030, C033 and C035 | defaults valid | REQ-010b, REQ-013 |
| C025 | markdown | **How cascaded clustering works, and what follows.** AFDB50 links members to the final representative *transitively, through a chain of intermediate representatives*, so identity to the representative and between members is unbounded below; hence a measured local identity of 26.73% to the best member over only residues 80-171 of 414. Consequence: a best or worst member can be a **distant homologue sharing a fold or domain family but differing in biological function** (TDP-43 and Prp24p are both RRM proteins). Identity and coverage measure how distant, not whether function is shared; read them with the description and species. Cite Steinegger and Söding 2018 and Barrio-Hernandez et al. 2023; this is a known phenomenon, not a discovery | none | none | no | none | appears **before** C026 | REQ-E7, REQ-E4, REQ-E9 |
| C026 | code | Identify best and worst from the sorted arrays; apply the collapse rule when the user's protein is an extreme | C010, C008 | `view_targets` | no | none | collapses to two targets when the query is an extreme | REQ-009, REQ-010 |
| C027 | code | Fetch `cifUrl`, `sequence` and `msaUrl` for the non-query view targets (no spacing) | C026 | target records enriched | no | **Must not** inherit the cluster endpoint's 2.5 s spacing. A deleted UniProt entry may lack a sequence | at most 2 extra calls; missing sequence degrades gracefully | REQ-011, NFR-003 |
| C030 | visualization | Mol\* views via `$molviewspec-rendering`, pLDDT-coloured, lazy import in `try/except ImportError` | C026, C027 | 2 or 3 viewers | no | **WebGL context budget.** At most 3 | degrades with a named message | REQ-V6, REQ-010 |
| C030a | visualization | Per-residue pLDDT profile for each rendered model, from the downloaded mmCIF | C027 | 2 or 3 `Figure` | no | Shows whether a low average is a uniformly poor model or an ordered domain plus a disordered tail, which is what REQ-E8 warns about | bare Figure; one per view target | REQ-010a, REQ-V8 |
| C031 | code | **Local** alignment of the user's protein against best and worst; compute identity, aligned length and query coverage | C027, C008 | `alignments`, identity, aligned length, coverage | no | Collapse rule must match C026 exactly. **Mandatory: must raise if the numbers cannot be produced** | convention stated; matches an independent local alignment on pinned sequences | REQ-011, REQ-011a |
| C032 | code | **Print** identity, aligned length and coverage as text, in the form "X% identity over N aligned columns, covering Y% of the query" | C031 | printed lines | no | **This cell is mandatory and must not be inside the figure's try/except** | numbers print with pyMSAviz uninstalled | REQ-011a |
| C033 | visualization | pyMSAviz alignment figures in the selected colour scheme, with identity, aligned length and coverage in the title | C031, C023b | 1 or 2 `Figure` | no | Lazy import; scheme from `cleancolors.json` **with a built-in module default when that file is absent** | degrades to a named message; C032's numbers already printed | REQ-V7, REQ-011 |
| C033a | markdown | Introduce the structural superposition: what it adds over the sequence alignment, how to read TM-score, and that **TM-align computes it while Mol\* only displays it** | none | none | no | Caveat before output. Must not imply Mol\* aligns anything | TM-score bands stated; RMSD explicitly de-emphasised | REQ-017, REQ-017a, REQ-017d |
| C033b | code | Install `tmtools` if absent; download each view target's mmCIF **once, verified complete**; TM-align the query against each extreme; print TM-score, RMSD, aligned length and identity | C026, C027 | `structures`, `superpositions` | no | **Self-contained install**, so the section can be removed whole without touching the shared install cell. Completeness matters more here than anywhere: a truncated mmCIF parses into a partial chain and TM-align then reports a confident transform for it | text summary prints outside any figure guard | REQ-017, REQ-017a, NFR-012 |
| C033c | markdown | Explain the per-model block, the two alignment panels, and expand **SSE** as secondary structure element | none | none | no | Caveat before output. SSE is expanded here because the track label inside the figure has no room for it | states that a column means spatial equivalence, not similarity | REQ-017b, REQ-017c |
| C033d | visualization | **One block per model**: superposition in Mol\*, then the same alignment twice, in the colour scheme and in pLDDT bands, each carrying the SSE track. Band and SSE legends shown once above the blocks | C033b, C023b | 1 or 2 viewers and 2 or 4 `Figure`, interleaved | no | **Grouped by model, not by output type**, so each superposition sits with the alignment that explains it. **WebGL context budget**: adds up to 2 on top of C030's 3, so 5 of the browser's 16. Matrix passed column-major. Both panels reuse the single verified download from C033b | viewer and panels degrade independently, each with a named message; C033b's numbers already printed | REQ-017b, REQ-017c, REQ-V11, REQ-V12 |
| C034 | code | Resolve `msaUrl` for the view targets; detect 403 by status and content type, never by parsing the body | C027 | `msa_sources` | no | **Endpoint currently 403.** This path is untested against a working service | named message on 403; continues | REQ-012 |
| C035 | code | Optional user-supplied a3m; validate by gap-stripped sequence match | C023b, C026 | `msa_sources` updated | no (path set in C023b) | Header convention is a hint, never the validation | correct name plus wrong sequence is rejected | REQ-013 |
| C036 | visualization | MSA coverage plot, ColabFold style, a3m insertions stripped | C034, C035 | 0 to 3 `Figure` | no | Lowercase a3m characters are insertions | bare Figure or named skip; **exercised offline by the synthetic a3m fixture**, so the plotting code is not shipped unexecuted | REQ-V9, REQ-012 |
| C037 | code | Summary: family level, the user's position, identity and coverage against each extreme, bounding caveats | C013, C016, C021, C031 | printed summary | no | none | values match the run | REQ-016 |
| C038 | markdown | What this notebook cannot tell you; what is deferred to a later version | none | none | no | none | states that candidate recommendation is out of scope and why | REQ-015 |

## Cell Order And Hidden-State Hazards

| Hazard | Where | Mitigation |
|---|---|---|
| `sys.path` inserted twice on re-run | C002 | Guard against duplicate entries before inserting |
| Matplotlib global state mutated on import | C004 | `apply_plot_style()` is explicit; importing the module must not mutate anything |
| Rate policy leaking between endpoints | C010 versus C008, C027 | Spacing lives inside each endpoint's own fetch helper, never in shared calling code |
| Query assumed present in cluster | C016, C021, C026 | All three must handle absence explicitly; absence is the normal Foldseek path |
| Display settings read before they are set | C023b | **Resolved by placement.** C023b sits before every cell that reads it: C023c, C024, C033 and C035. `TAXONOMIC_RANK` reaches both the lineage selection and the figure, so it is a live knob rather than a dead one |
| A mandatory output buried inside an optional section | C031 to C033 | Identity, aligned length and coverage are computed in C031 and printed in C032, both mandatory. Only the pyMSAviz figure in C033 is optional. C025's claim that identity is how the reader checks is therefore backed by a cell that cannot silently vanish |
| Foldseek section refusing the whole run | C021 | A small Foldseek cluster skips only that section. The run is refused only when the **sequence** cluster is unusable |
| Best or worst member has no sequence | C027, C031 | A deleted UniProt entry may lack one. Degrade the alignment with a named message; the 3D view may still render |
| Figures double-displayed | all visualization cells | Bare `Figure` built without pyplot; the notebook's `%matplotlib inline` registers formatters |
| Optional-section failure aborting the run | C019, C023c, C024, C030, C033, C033d, C034, C036 | Each wrapped to print a named message and continue. **C031 and C032 are excluded**: the alignment numbers are mandatory and must raise, per REQ-011a |
| Collapse rule diverging between sections | C026 versus C031 | Computed once in C026 and passed forward; never recomputed |

## Install And Kernel Assumptions

- C002 and C003 are the only cells that touch the environment. Neither may run implicitly from an
  import.
- C003 must be a no-op when the three packages already import cleanly, so local runs cost nothing.
- Colab preinstalls numpy, pandas, matplotlib, seaborn and requests, so C003 covers `molviewspec`,
  `biopython` and `pymsaviz` only.
- `tmtools` is imported lazily inside `superpose()`, and installed by C033b rather than by the
  shared install cell, so the whole superposition section can be deleted as one block.
- `molviewspec` and `pymsaviz` are imported lazily inside C030 and C033, so an install failure
  degrades the optional sections without breaking the mandatory ones.
- No cell caches to disk. Five API calls per run makes re-fetching cheaper than reasoning about
  cache invalidation.
