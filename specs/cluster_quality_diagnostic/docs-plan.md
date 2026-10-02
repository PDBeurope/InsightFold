# Documentation Plan

The notebook must stand alone. A reader who opens it on Colab with no other context should be able to
run it and interpret the result. This plan assigns each documentation type a home and an owner
section, so that documentation is written alongside the notebook rather than after it.

| Documentation Type | Artifact | Audience | Notes |
|---|---|---|---|
| Tutorial | Notebook cells C001, C005, C009 | New users | How to run it: edit one accession, run all. What they will get. What the two clusterings are |
| How-to | Notebook cells C023a and C023b plus a short README section | Returning users | All four are Colab form widgets: the clustering the 3D models come from, the alignment colour scheme, whether PAE is drawn, and a path to the reader's own MSA |
| Reference | `specs/cluster_quality_diagnostic/data-contracts.md` plus module docstrings | Developers and power users | Inputs, outputs, parameters, API contracts, rate policies, the percent-identity convention |
| Explanation | Notebook caveat cells C012, C015, C020, C025, C037 | All users | Scientific background, limitations, interpretation, and what the notebook deliberately does not do |

## In-Notebook Documentation Requirements

These are notebook content, not external docs, and they are the parts most likely to be skipped
under time pressure. Each maps to a requirement and a cell.

| Content | Cell | Requirement |
|---|---|---|
| What the notebook answers and what it does not | C001 | REQ-015 |
| How to run it: type an accession into the form, then Run all | C005, C006 | REQ-001, REQ-018 |
| What AFDB50/MMseqs2 and AFDB/Foldseek are, quoting the AFDB FAQ criteria, with the InterPro framing and the statement that these clusters approximate rather than instantiate those concepts | C009 | REQ-E1 |
| **Cluster summary.** A glossary of every term the summary prints (members, mean, median, standard deviation, p25 to p75, confidence bands), each with what it is and what to read from it. No abbreviations left unexpanded except pLDDT itself. **No sample-derived constants**: the reader applies the interpretation rules to their own cluster's numbers | C012 | REQ-003 |
| That pLDDT is not an experimental tractability predictor | C012 | REQ-E5 |
| **Your protein within the cluster.** The descriptor and its thresholds, and why a rank is only worth what this cluster's own p25 to p75 width makes it worth. No constants from any other sample | C015 | REQ-E2, REQ-004 |
| That structure-based clusters are broader and lower-scoring than sequence-based ones, which reverses a common expectation, with the mechanisms stated | C020 | REQ-E3 |
| How MMseqs2 cluster membership works: AFDB50 uses **cascaded** clustering, so members are linked to the final representative *transitively, through a chain of intermediate representatives*. Identity to the representative and between members is unbounded below, which is why the best cluster member can align to the query over only a single domain. Cite Steinegger and Soding 2018 | C025 | REQ-E7 |
| That a member can therefore be a distant homologue with different biological function. Worked example: TDP-43 against Prp24p, 26.73% identity over 101 columns covering residues 80-171 of 414, which is an RRM1-only match. Both are RRM proteins, so this is remote homology, not an unrelated protein. Framed as a consequence of the clustering mechanism and cited, not presented as a discovery | C025 | REQ-E4, REQ-E9 |
| What the notebook cannot tell you, and that candidate recommendation is deliberately out of scope because no single available number reliably separates a divergent member from a legitimate distant homologue | C037 | REQ-015 |
| That average pLDDT is substantially a function of disorder content, so a low-mean family is not necessarily a badly predicted one. The per-residue profiles show the difference | C012, C030a | REQ-E8 |
| That the superposition is computed by **TM-align** and only displayed by Mol\*, and that TM-score rather than RMSD is the number to read | C033a | REQ-017a, REQ-017d |
| That a column in the structural alignment means the two residues occupy the same position in space, not that they resemble each other, and what the SSE track codes mean, with the abbreviation expanded | C033c | REQ-017b, REQ-017c |
| That the MSA section's automatic path is untested while the endpoint returns 403 | C034 | REQ-012 |

## Tone And Convention Requirements

- British spelling throughout (NFR-009).
- **No internal identifiers in the notebook.** Cell IDs, requirement IDs and fixture IDs belong to the spec pack. A reader sees none of them, and nothing in the notebook would explain them.
- **Markdown prose is not hard-wrapped.** One line per paragraph; let the renderer wrap it to whatever width the reader's window is. Tables, list items, headings and fenced code keep their own line structure. **This includes prose inside raw HTML callouts and consecutive `#@markdown` lines in a form cell**, both of which are concatenated before rendering and so wrap just as artificially. Those two were the only places the rule had drifted, found 2026-10-02.
- No em dashes in markdown, plot titles, axis labels or printed output. En dashes inside numeric
  ranges are fine. Check generated HTML for `&mdash;` and `&#8212;`, which no Unicode grep will find.
- InterPro's vocabulary for family and homologous superfamily.
- Caveats are placed **before** the output they bound, never after. A caveat below a figure has
  already failed.
- Confident where the evidence supports it, explicitly uncertain where it does not. The notebook
  states what it measured and what it assumes.

## External Documentation

| Artifact | Change needed |
|---|---|
| `README.md` | Add `cluster_quality_diagnostic` to the notebook list with a one-line description and its default accession |
| `CLAUDE.md` | No change required. This notebook does not touch `complex_interface_utils.py`, and its dependency policy is per-consumer. If it later shares code with the dimer work, D8 applies: promote on the second consumer, not in advance |
| `specs/cluster_quality_diagnostic/` | This spec pack is the reference documentation. Keep `data-contracts.md` current if AFDB changes |
| `prompts/templates/cluster-quality-diagnostic-decision-capture.md` | The scoping record. Do not update it as the notebook evolves; it is a point-in-time decision log, and its value is that it records what was believed when |

## Documentation Debts To Record Rather Than Hide

| Debt | Where stated |
|---|---|
| The MSA automatic path is untested against a working endpoint | Notebook C033 and `validation.md` |
| Colab is unverified until T074 runs | `validation.md`. Must not be described as Colab-ready before then |
| Candidate recommendation is deferred, not forgotten | Notebook C037, PRD section 4, traceability matrix |
| The taxonomic composition figure was removed: the cluster API carries no lineage, so it depended on a secondary UniProt lookup to say anything | `data-contracts.md` section 4, `requirements.md` OQ-7 |
