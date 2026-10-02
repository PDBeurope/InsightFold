# Notebook UX Contract

## Primary User Question

> **How is pLDDT distributed across my protein's family, where does my protein sit in it, and what do
> the best and worst-predicted members look like?**

The order is load-bearing: the cluster is summarised before the user's protein is placed in it,
because how much a rank is worth depends on how wide that cluster's middle is. Both live in one
section, and the notebook derives the judgement from the user's own cluster rather than quoting
constants measured on some other sample.

## Target User And Expertise

Primary: a structural biologist or modeller wanting context for their model's confidence. Secondary:
an AFDB analyst looking at family-level quality, and a bench scientist wanting to know whether a
protein is well predicted relative to its relatives.

**Assumed:** familiarity with pLDDT, protein families, and structural homology.

**Not assumed:** any knowledge of how AFDB clustering was constructed. This must be explained,
because it determines how every figure should be read.

**Actively guarded against:** reading confidence as experimental tractability; reading a narrow
percentile spread as a ranking; assuming a cluster member is the same protein.

## First Runnable User Input Cell

One parameter cell, near the top, visually separated:

```python
# ─── Colab form fields, not editable code (REQ-018) ──────────────────────
ACCESSION = "P69905"        # a bare UniProt accession, e.g. Q9I1F6, Q13148
# ─────────────────────────────────────────────────────────────────────────
```

Display preferences live in a second cell placed later, before the sections they affect, so the
entry point stays a single line:

```python
ALIGNMENT_COLOUR_SCHEME = "clustal2"
USER_MSA_PATH           = None         # optional a3m for a highlighted model
MODEL_SOURCE_CLUSTER    = "sequence"   # which clustering the 3D models come from
SHOW_PAE                = True         # plot each model's PAE beside its 3D view
```

## User-Editable Parameters Versus Internal Variables

| User-editable | Internal, must not be presented as a knob |
|---|---|
| `ACCESSION` | Cluster and prediction API base URLs (NFR-008) |
| `ALIGNMENT_COLOUR_SCHEME` | Request spacing per endpoint, which differs between the two endpoints (NFR-003) |
| `USER_MSA_PATH` | pLDDT band boundaries and palette |
| `MODEL_SOURCE_CLUSTER` | Mol\* viewer count budget |
| `SHOW_PAE` | The shared x-axis range across the two rainclouds |
| | The PAE colourmap and its scaling |


## Happy-Path Default Example

**`P69905`**, haemoglobin alpha. 5,223 MMseqs2 members: familiar, fast, and its top-20 median
identity to the query is 86.6%, so a first run shows a coherent family where the best member is
plainly the same protein.

Changed from `Q9I1F6` on scientific review. Its best member sits at 25.09% identity, which is
*lower* than the 26.73% of the case the notebook elsewhere describes as functionally divergent. Putting it on the default path would show a first-time reader the
ambiguous case before they have the vocabulary for it. `Q9I1F6` is retained as a fixture precisely
because it is the distant-but-legitimate case.

## User Flow

| Stage | User does | Notebook does | User sees | User learns |
|---|---|---|---|---|
| 1 | Edits `ACCESSION`, runs all | Validates, normalises, refuses if unsupported | Accession echoed with description and organism | That the notebook understood its input |
| 2 | Reads | Explains both clusterings | Two definitions plus the InterPro framing | What population is about to be shown |
| 3 | Waits | Fetches both clusters, 2.5 s apart | Member counts, patience note for large clusters | Scale of the family |
| 4 | Reads | Summarises the cluster, then locates the user's protein in it, then plots **one** raincloud with the protein marked | A term glossary, the statistics, the protein's own value and descriptor, and the figure | What quality this family achieves and whether their protein is typical of it |
| 5 | Reads | Length against pLDDT as a joint density | One figure with marginals | How length relates to confidence here |
| 4b | Reads | Structure-cluster raincloud, directly beneath the sequence one, no query marker | A second distribution to compare against the first | That structurally similar proteins are a different and broader set |
| 7 | Inspects | Renders 2 or 3 Mol\* views, each beside its PAE matrix | pLDDT-coloured structures and their predicted aligned error | Where in the structure the differences sit, and which parts are positioned confidently relative to each other |
| 8 | Inspects | Locally aligns against best and worst | Identity, aligned length, coverage and aligned range as text, plus 1 or 2 figures | How distant those extremes actually are |
| 9 | Reads | MSA coverage where available | Coverage plot or a named message | Whether depth tracks the differences |
| 10 | Reads | Summarises | Family level, position, caveats | What to take away |

**Interpretation happens at stages 4 and 8.** Trust messaging must land before each, not after.

## Expected First Visible Confirmation

Within a few seconds of running all, and **before** any call to the cluster endpoint, the user must
see their accession echoed with its UniProt description, organism, length and global pLDDT, taken
from the prediction endpoint, which needs no rate spacing. If that line names the wrong protein, they
should stop. This is the cheapest possible early failure signal.

The second confirmation is the member counts for both clusterings, which shows progress and gives the
scale of what follows.

## Fixtures Versus User Flow

Fixtures must not become the user flow:

- The default `ACCESSION` is a real, interesting protein, not a test artefact.
- No fixture list, test accession table or validation block appears in the notebook body.
- Pinned cluster JSON exists for the **test harness**, which lives outside the notebook.
- The notebook has no `TESTING` branch and no fixture-loading path. It always calls the live API. A
  notebook that behaves differently under test is not the notebook being tested.

## Trust, Caveats, And Limitations To Show In Notebook

Each is a markdown cell placed **before** the output it bounds, never after.

| Placement | Message | Requirement |
|---|---|---|
| Before stage 2 figures | These clusters approximate InterPro families; they are not curated entries | REQ-E1 |
| Before the percentile (stage 4) | A rank is only worth what this cluster's own p25 to p75 width makes it worth. The extremes are informative; the body may not be | REQ-E2, REQ-004 |
| Before the structure raincloud (stage 4b) | Proteins grouped by structural similarity. Structure is more conserved than sequence, so Foldseek detects remote homologues whose sequence identity has fallen below the threshold for reliable inference. The user's protein is usually not shown, because the clustering is built from cluster representatives | REQ-E3, REQ-006 |
| Before the 3D and alignment sections (stage 7) | Cascaded clustering links members transitively, so a member can be a distant homologue with different biological function. Worked example: a human TDP-43 query's best cluster member is yeast Prp24p, aligning at 26.73% identity over residues 80-171 of 414, an RRM1-only match. Both are RRM proteins. Read identity, coverage, aligned range, description and species together; no one of them decides it. A known phenomenon, cited, not a discovery | REQ-E4, REQ-E7, REQ-E9 |
| Alongside the distribution | pLDDT is not an experimental tractability predictor, has not been validated against experimental structures for these families, and moves with AFDB releases | REQ-E5 |
| Alongside the distribution | Average pLDDT over a chain is substantially a function of disorder content, so a low-mean family is not necessarily a badly predicted one. The per-residue profiles show the difference | REQ-E8 |
| Stage 11 | The MSA endpoint is under maintenance; this path is untested | REQ-012 |

**Tone.** The notebook is descriptive. It shows the distribution and its extremes and does not tell
the user to switch models. Where the best-predicted member is only distantly related, the notebook
does not hide it or filter it out; it reports identity, aligned length and coverage beside the
member's description and species, and leaves the judgement to the reader. It does not claim that a
low identity means a different protein, because identity alone cannot establish that: 26.73%
(divergent), 25.09% and 33.68% (both legitimate) across three measured clusters, so the divergent
case scores highest of the three.
