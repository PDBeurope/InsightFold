# Cluster Quality Diagnostic Spec Pack Overview

| Field | Value |
|---|---|
| Source PRD | `prd/cluster_quality_diagnostic_notebook_prd.md` (v0.2) |
| Status | Draft v2, revised after readiness and scientific spec review (both returned `fail`; all blocking and major findings addressed). Ready for re-review |
| Target Notebook | `notebooks/cluster_quality_diagnostic.ipynb` |
| Target Module | `src/insightfold/cluster_quality_utils.py` |
| Decision record | `prompts/templates/cluster-quality-diagnostic-decision-capture.md` |
| Date | 2026-09-21 |

## Summary

The notebook takes one UniProt accession and shows **how pLDDT is distributed across that protein's
family**, where the user's protein sits in that distribution, and what the best and worst-predicted
members of the cluster look like as 3D structures and in pairwise alignment against the user's
protein.

It is descriptive. It does not rank candidates or recommend that the user switch models.

Three measured findings shape the design. All three were revised after scientific spec review.

1. **Between-cluster spread dwarfs within-cluster spread.** Between-cluster standard deviation is
   **12.46** against a median within-cluster standard deviation of **4.03**, with a median IQR of
   **4.2** pLDDT points (40-cluster stratified random sample). The family's level is therefore the
   headline and the user's rank within it matters mainly at the extremes. *An ICC(1) of 0.798 was
   computed from the same sample but is **not** used in user-facing prose: group sizes span 50 to
   93,793, pLDDT is ceiling-compressed near 100, and the stratified frame fixes between-cluster
   spread by construction. The three direct numbers carry the argument without those assumptions.*
2. **Foldseek clusters run lower than MMseqs2 clusters** (mean lower in 28 of 32). The decisive
   argument: Foldseek members are each their sequence family's highest-pLDDT representative, which
   is *positive* selection on pLDDT, so a **lower** mean can only come from aggregating
   lower-scoring families. A Foldseek distribution is therefore one point per family scored by its
   best member, further conditioned on structural alignability. The notebook says this before
   showing the figure, and draws no query marker on it.
3. **Cluster members can be distant homologues with different biological function.** The twenty
   best-scoring members of a TDP-43 cluster are yeast Prp24p at 26.73% local identity, aligning over
   residues 80-171 of a 414-residue query, which is an RRM1-only match. Both are RRM-domain proteins,
   so this is remote homology rather than an unrelated protein, and the notebook says so.
   **This is a known phenomenon, not a discovery**: cascaded clustering links members transitively
   through chains of intermediate representatives (Steinegger and Söding 2018, Nat Commun 9:2542),
   and AFDB structural clustering merges sequence clusters across remote homology
   (Barrio-Hernandez et al. 2023, Nature 621:637). The notebook reports identity, aligned length and
   query coverage and lets the reader judge, rather than filtering.

A fourth constraint, added on review: **average pLDDT over a chain is substantially a function of
disorder content**, so cross-member differences partly measure intrinsic disorder rather than
prediction quality. The notebook must not call a low-mean family "low quality", and it plots a
per-residue pLDDT profile for each rendered model to make the distinction visible.

## Scope Reduction From v0.1

An earlier revision specified a ranked candidate table with eligibility filters and an identity
screen. **All of it was cut.** Measurement showed that no available signal separates a confirmed
true positive (the TDP-43 case) from a confirmed false positive (a `Q9I1F6` cluster whose top 20 are
legitimate distant LacI-family homologues at similar identity). The signal that does separate them is
family-level annotation, which the owner deferred to a later version.

Rather than ship a recommendation the notebook cannot stand behind, v1 shows the distribution and
its extremes with the identity evidence visible, and leaves the judgement to the reader. The cut is
traced in `traceability-matrix.md` so a later version can pick it up without re-deriving the
reasoning.

## Blocking Questions

**None.** No open question prevents implementation starting.

## Assumptions Accepted For Build

| ID | Assumption | Basis | Revisit when |
|---|---|---|---|
| A-1 | Taxonomic rank defaults to phylum, with a top-6 plus `Other / unresolved` long-tail rule | Display default; phylum is uninformative for some eukaryotic clusters, hence the user override | A reviewer finds the default unreadable on a common case |
| A-2 | **REVISED.** Identity is computed by **local** alignment and reported as "X% identity over N aligned columns, covering Y% of the query" | Global identity penalises length and domain-architecture mismatch, which is exactly what this notebook surfaces, and nothing bounds query-to-member length ratio. BLOSUM62 with −11/−1 is also the BLASTP local default | Stays stated in the notebook. Revisit only if a reviewer prefers another convention |
| A-3 | `P12345` remains a 5xx fixture | Its 500 may be rate limiting rather than a per-accession fault | T003 re-tests it in isolation |
| A-5 | Typical versus extreme is defined as outside p5 or p95 **of that cluster's own distribution** | A relative rule, never a global cutoff. Stated in the notebook | A reviewer proposes a better rule |
| A-6 | Clusters below 30 members get a named low-n caveat rather than a higher refusal gate | A percentile and an IQR are unstable at n of 3 to 30, but refusing loses real information | Reader review finds the caveat insufficient |
| A-4 | The cluster API's descending sort holds | Confirmed on all 13 responses examined | The module asserts it and sorts defensively if it fails |

## Non-Blocking Open Questions

Carried from the PRD; none gate implementation.

| ID | Question | Next step |
|---|---|---|
| OQ-1 | RESOLVED. Four accessions return a byte-identical cluster they are not in | Most likely shared amino acid sequence resolving to the same cluster. Expected behaviour; see note 1 below |
| OQ-2 | Is the Foldseek widening superfamily aggregation or cluster-size disparity? | Affects REQ-005 prose only. Structural biology review (T092) |
| OQ-3 | Can `clusterTotal` disagree with the returned array length? | REQ-002 records both, so it surfaces in output rather than silently |
| OQ-4 | Is `P12345`'s 500 a fault or rate limiting? | T003 |
| OQ-5 | MSA formats beyond a3m | Spec-level; a3m only in v1 |
| OQ-6 | Flat `clustal` scheme or true Clustal X colouring? | User preference; flat schemes ship |
| OQ-7 | Taxonomic rank default and long-tail cutoff | Accepted as A-1 |

## Required Advisory Reviews

| Review | Role | What they judge | Task |
|---|---|---|---|
| Spec review | `$notebook-spec-review` | This pack, before implementation | T001 |
| Fixture curation | `$fixture-selection` | Confirm and freeze the manifest; pin the JSON | T002 |
| Scientific correctness | `computational-structural-biologist` | The Foldseek description (OQ-2); whether cluster impurity has existing literature to cite rather than being presented as a novel observation | T092 |
| Interpretability | Domain reader who did not build it | Whether the three counter-intuitive explanations land: the narrow spread, the Foldseek reversal, and cluster impurity | T093 |
| Execution validation | `$notebook-execution-validation` | Restart-run-all, hidden state, dependency budget, Colab | T090 |
| Final review | `$notebook-review` | Scientific quality, pedagogy, maintainability, lifecycle readiness | T091 |

## Notes On AFDB Data Behaviour

Observed during scoping. **None of these is treated as a defect or requires action.** Recorded so a
future reader does not re-investigate them.

1. **Identical cluster responses across different accessions.** `Q9KM69`, `A0A7X7SVK1`,
   `A0A536HFQ7` and `A0A8C8MD66` each return a byte-identical cluster response to another
   accession's. Most likely explanation: they share the same underlying amino acid sequence, so AFDB
   resolves them to the same cluster. Expected behaviour, not an anomaly.
2. **`referenceLabel` reflects the UniProt release AFDB was built against**, not the current one.
   AFDB and UniProt drift between releases. Not worth validating per accession; revisit only if it
   causes an observable problem.
3. **Some cluster members no longer resolve in UniProt.** Same cause. The notebook handles it where
   it matters: REQ-017 explains it at the entry point, and REQ-011's alignment path degrades with a
   named message when a member's sequence is unavailable.
4. **AFDB50 clusters span wider sequence divergence than "50% identity" suggests**, because MMseqs2
   binds each member to the cluster representative rather than to every other member. Expected for
   cascaded clustering. **The notebook explains this** (REQ-E7) rather than flagging it, since it is
   the reason distantly related members appear in a cluster at all.

## Pack Contents

| File | Purpose |
|---|---|
| `spec-pack-overview.md` | This file |
| `requirements.md` | REQ, NFR, educational and visualisation requirements; edge cases |
| `notebook-ux-contract.md` | User question, entry point, flow, trust messaging |
| `notebook-design.md` | Section outline, data flow, handoff, dependency policy, design decisions |
| `cell-blueprint.md` | Cell-by-cell plan with hidden-state hazards |
| `traceability-matrix.md` | Goals to requirements to cells to fixtures to validation |
| `data-contracts.md` | Every external API and file contract, measured |
| `fixture-manifest.md` | Fixtures, expected snapshots, provenance |
| `validation.md` | Validation plan and what is deliberately not asserted |
| `tasks.md` | Executable task list, phases A to G |
| `docs-plan.md` | Documentation types, in-notebook requirements, recorded debts |
