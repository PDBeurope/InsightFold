# Validation Plan

| Check | Method | Pass Criteria | Fixture / Artifact | Requirement |
|---|---|---|---|---|
| Restart-and-run-all | Restart kernel, run all cells in order, once | No errors, no manual intervention | `P69905` | NFR-002 |
| Restart-and-run-all, all fixtures | Same, per fixture accession | Passes on every non-refusal fixture | all | NFR-002 |
| Cell-order independence | Run all; then restart and run only C001 to C016 | No cell depends on a later cell. C023b's display settings are read only by C023c, C024, C033 and C035, all of which follow it | `P69905` | NFR-002 |
| No hidden state | Run all twice in one kernel | Second run produces identical values. `sys.path` is not duplicated | `Q9I1F6` | NFR-002, C002 |
| Colab execution | Actual Colab free-tier run, recorded with a link or screenshot | Completes within roughly 60 s of install time | `P69905` | NFR-001 |
| Dependency budget | Grep all imports in notebook and module | Nothing outside numpy, pandas, matplotlib, seaborn, requests, molviewspec, biopython, pyMSAviz | source | NFR-007 |
| Prohibited imports absent | Grep | No torch, gemmi, plotly, scipy, ptitprince, bigquery, networkx | source | NFR-007 |
| Lazy imports | Inspect | molviewspec and pymsaviz imported inside functions, not at module top | source | NFR-007 |
| Distribution statistics | Compare module output against directly computed numpy values on pinned JSON | Exact agreement on mean, median, sd, p25-p75 width, band fractions | all pinned | REQ-003 |
| Query position | Compare against a hand-computed percentile on pinned JSON | `Q9I1F6` at index 1,650 of 4,628; percentile matches | `Q9I1F6`, `P69905`, `Q13148` | REQ-004 |
| Query absence handled | Run the Foldseek path | Section completes, names the representative, does not raise, and computes no query percentile there | `P69905`, `Q9I1F6` | REQ-006 |
| Extremes identified | Compare against the API's descending sort on pinned JSON | Best is index 0, worst is index −1, **asserted not assumed** | all pinned | REQ-009 |
| Sort assumption guarded | Feed a deliberately unsorted pinned response | Module sorts defensively and logs that the assumption failed | synthetic | REQ-009 |
| **Distance is visible** | Run the alignment section | FX-02 reports **26.73% identity over 101 aligned columns, covering 21.3% of the query, residues 80-171 of 414**, printed and in the figure title. Assert the identity, columns and range; do **not** assert coverage as a threshold, since it moves between 21% and 55% across gap settings | FX-02 | REQ-011, REQ-E4 |
| Low identity is not overinterpreted | Run the alignment section | FX-05 (33.68%, 84.4% coverage) and FX-03 (25.09%, 74.6%, residues 18-277 of 343) are reported without any claim of a different protein. FX-03's identity is **below** FX-02's | FX-05, FX-03 | REQ-011, REQ-E4 |
| Missing sequence handled | Feed FX-12, a prediction response with the sequence fields removed | Alignment degrades with a named message; run continues. **No real accession is known to trigger this**; AFDB serves sequences even for UniProt-deleted members | FX-12 | REQ-011 |
| Collapse rule | Run FX-11 both ways: query at index 0, then at index −1 | Exactly 2 view targets and 1 alignment in both | FX-11 | REQ-010, REQ-011 |
| Identity and coverage are mandatory | Uninstall pyMSAviz and run | Identity, aligned length and coverage still **print**; only the figure is skipped. With the numbers unavailable, the run raises rather than continuing | `P69905` | REQ-011a |
| Local alignment convention | Compare against an independent **local** alignment on pinned sequences | Identity, aligned length and query coverage all agree to within rounding | `P69905`, `Q13148`, `Q9I1F6` | REQ-011 |
| **Coverage and range are evidence, not a rule** | Compare FX-02 and FX-03 output | FX-02: 26.73% over 101 columns, 21.3% coverage, residues 80-171. FX-03: 25.09% over 279 columns, 74.6% coverage, residues 18-277. Identity ranks them the wrong way round. Coverage separates them **at the default gap settings only**: measured across -6/-1 to -20/-1 the FX-02 coverage moves 21% to 55% and at -20/-1 the ordering inverts. The notebook must present these as evidence for the reader, never as a threshold | FX-02, FX-03 | REQ-011, REQ-E4 |
| Distant homologue not called a different protein | Read the FX-03 output | No prose claims a different protein. FX-02 is described as a distant homologue with different biological function sharing an RRM domain, not as an unrelated protein | FX-02, FX-03 | REQ-E4 |
| Foldseek figure carries no query marker | Inspect the figure and its axis label | Axis and caption state "one point = one AFDB50 family, scored by its highest-pLDDT member". No query reference line present | `P69905` | REQ-006, REQ-V2 |
| Partial degradation | Run FX-09 | The sequence analysis completes on its 380-member sequence cluster; only the structure section is skipped, with a named message | FX-09 (both flags pinned) | REQ-014a |
| Per-residue profile renders | Run | One pLDDT profile per view target. On `Q13148` it shows an ordered region plus a low-confidence tail rather than a flat low trace | `Q13148` | REQ-010a, REQ-V8 |
| MSA coverage plot is executed | Run against FX-13 | Coverage plot renders; lowercase insertions stripped before the alignment is treated as rectangular | FX-13 | REQ-V9, REQ-012 |
| Missing colour file | Remove `cleancolors.json` and run | Alignment figure renders with the module's built-in default scheme. No error | `P69905` | REQ-V7 |
| Low-n caveat | Run FX-14, n = 12 | Percentile and IQR printed with a named low-n caveat | FX-14 | REQ-004 |
| Low mean is not a messy cluster | Run FX-06 | Mean 51.72 with a best member at 96.15% identity and 99.7% coverage. Prose must not equate a low family mean with cluster heterogeneity | FX-06 | REQ-E8 |
| Isoform substitution is stated | Run `P04637-2` | Output names which sequence was aligned and which structure rendered | `P04637-2` | REQ-001a |
| Accession normalisation | Run `AF-Q9I1F6-F1` | Normalised to `Q9I1F6` and the interpretation printed. Not rejected | `AF-Q9I1F6-F1` | REQ-001 |
| Collapse rule consistency | Compare C026's targets against C031's alignment pairs | Identical; the rule is computed once and passed forward | `SYN-COLLAPSE` | REQ-010, REQ-011 |
| Refusal: fragmented | Run | Named message explaining fragmented proteins are absent from clustering. No traceback | `Q8WZ42` | REQ-014 |
| Refusal: too small | Truncate a pinned response below 3 members | Named message giving the cluster size | derived from FX-14 | REQ-014 |
| Refusal: malformed or absent | Run | 404 reported, quoting the service `detail`, **and naming both likely causes** with a suggestion to try an equivalent protein | FX-08 (404 on both endpoints), plus a junk string | REQ-014, REQ-017 |
| Refusal: 5xx | Run | Retries with backoff, then reports. Does not present rate limiting as a data error | `P12345` | REQ-014, NFR-003 |
| Isoform supported | Run | Resolves to the parent cluster | `P04637-2` | REQ-001 |
| Rate policy separation | Instrument request timing | Cluster calls at least 2.5 s apart; prediction calls not delayed | `P69905` | NFR-003 |
| Large cluster | Run | Completes without cap, pagination or sampling | `P0A6F5` (93,793 members) | NFR-004 |
| `clusterTotal` agreement | Run | Compared against array length and recorded in output | all | REQ-002, OQ-3 |
| MSA 403 handled | Run | Named message; detection by status and content type, never by parsing the body | any | REQ-012 |
| MSA path declared untested | Read the notebook | States the endpoint is under maintenance and this path is unverified | any | REQ-012 |
| User MSA validated by sequence | Supply FX-13's wrong-sequence variant, which carries the correct header | Rejected, naming the mismatch | FX-13 | REQ-013 |
| Mol\* renders or degrades | Run, then simulate ImportError | 2 or 3 viewers, or a named explanation | `Q9I1F6` | REQ-V6, REQ-010 |
| Viewer budget | Count live viewers | At most 3 | `Q9I1F6` | REQ-V6 |
| Figures are bare Figures | Inspect return types and pyplot registry | Every plot function returns `matplotlib.figure.Figure`; registry stays empty; nothing double-displayed | source | NFR-005 |
| No JS widgets | Grep and inspect | Mol\* uses the base64 data-URI IFrame pattern; alignments and MSAs are static | source | NFR-006 |
| Caveats precede their outputs | Read the notebook in order | Each caveat cell appears before the output it bounds | all | REQ-015 |
| Caveats are unavoidable | Review by someone who did not build it | A figures-only reader does not conclude the percentile is a ranking, nor that structure clustering predicts quality better | `Q9I1F6`, `Q13148` | REQ-015, REQ-E2, REQ-E3 |
| Prose conventions | Grep | British spelling; no em dash, `&mdash;` or `&#8212;` anywhere including generated HTML | source | NFR-009 |
| Outputs stripped | Inspect committed `.ipynb` | No `outputs` payload | git | NFR-010 |
| Tests run offline | Run the harness with the network disabled | Every assertion above that uses a fixture passes from pinned JSON | fixtures | fixture manifest |
| `TODO(merge)` resolved | Grep before release | Bootstrap pins `main` or a release tag, not a feature branch | source | NFR-001 |

## Reference Implementation And Tolerance

There is no external reference implementation for this notebook's outputs, unlike
`homodimer_diagnostic`'s comparison against `ipsae.py`. Statistics are compared against directly
computed numpy values on pinned inputs, which is exact agreement rather than a tolerance.

Percent identity is the one place a tolerance applies. Compare against an independent alignment
implementation and require agreement to within rounding, having first fixed the convention
(identical columns over the **local** aligned block, reported with aligned length and query
coverage). Note that ungapped and shorter-sequence
conventions differ by 2 to 6 percentage points, which is why the convention must be stated in the
notebook rather than only here.

## Assert Or Monitor: The Line Is Drawn By Input Source

**Assert exactly** any statistic computed from **pinned JSON**. It is a property of this notebook's
code: means, band fractions, percentiles, identities, aligned lengths, coverage, extreme selection.

**Monitor, never assert**, the same statistic recomputed from the **live API**. It is a property of
AFDB's data and will drift with AFDB releases. Record it in the fixture manifest with its capture
date and investigate changes, rather than failing a build over them.

This supersedes the earlier wording that forbade asserting band fractions outright, which conflated
the two cases and contradicted REQ-003's own acceptance criterion.

Figure images are never snapshotted.

## Human Review Needs

| Review | Who | What they judge |
|---|---|---|
| Scientific correctness | Computational structural biologist | Whether the Foldseek description (OQ-2) is right; whether cluster impurity has literature to cite rather than being presented as novel |
| Interpretability | A domain reader who did not build it | Whether the three counter-intuitive explanations land: the narrow spread, the Foldseek reversal, and cluster impurity |
| Visualisation quality | Whoever reviews figures | Whether each figure is legible, correctly coloured, and cannot support a conclusion the prose contradicts |
| Pedagogy and documentation | Technical reader | Whether a new user can run it and interpret the output without this spec |
