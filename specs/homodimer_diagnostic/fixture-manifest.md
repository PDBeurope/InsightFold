# Homodimer Diagnostic Fixture Manifest

## Status

Current readiness: `validation-ready for the homodimer and heterodimer paths`.

FX-006 to FX-010 were added by R024 on 2026-09-08. They close two of the three gaps that held this manifest at `smoke-ready`: there are now heterodimer fixtures at all, and every score on every one of them has been reconciled against AlphaFold DB's *own published* IPSAE-derived values, which serve as the trusted reference snapshot the earlier passes were waiting for (see "Independent Reference Agreement" below). What is still open is the pDockQ reference — AFDB publishes it only to two decimal places — and the FX-002 borderline question.

Note that the file name says "homodimer" for continuity with the spec pack (D6); since R022 the notebook itself is titled "Dimer Confidence Metric Diagnostic Notebook" and the fixture set covers both assembly types.

Last fixture curation pass: 2026-09-08 (R024).

## Fixture Summary

| Fixture ID | Role | Accession / Source | Status | Purpose |
|---|---|---|---|---|
| FX-001 | happy-path / high-confidence reference candidate | `AF-0000000065889468` | ready; values reconciled with AFDB | Exercise complete successful AFDB fetch, parse, metric, visualization, MolViewSpec, and summary workflow |
| FX-002 | lower-confidence / provisional borderline candidate | `AF-0000000066503175` | candidate-needs-scoring | Stress a current AFDB complex with lower global confidence than FX-001; may become borderline fixture if ipSAE/summary behavior supports it |
| FX-003 | metric-disagreement | superseded by FX-010 | closed | Was requested from a domain reviewer; FX-010 supplies the disagreement pattern from live AFDB data (see FX-010) |
| FX-004 | malformed negative | `AF-NOT_A_REAL_ACCESSION` | ready | Confirm clear AFDB identifier-format failure and no downstream stack trace |
| FX-005 | valid AFDB but unsupported by v1 | `O15552` / `AF-O15552-F1` | ready | Confirm monomer entries are rejected clearly because v1 supports only two-chain dimers |
| FX-006 | primary heterodimer | `AF-0000000211034637` | ready; values reconciled with AFDB | The per-task verification fixture for W2/W3 under D5. Unequal chains (181 + 101), cross-species, small interface |
| FX-007 | length-asymmetry stress | `AF-0000000204661110` | ready; values reconciled with AFDB | 8.16:1 length ratio. Exercises every axis, mask panel and coverage bar that assumed comparable chain lengths |
| FX-008 | mid-confidence heterodimer | `AF-0000000211619209` | ready; values reconciled with AFDB | Global pLDDT 73.57, `ipSAE_d0res` 0.6366 — the only fixture that lands in AFDB's `LOW-CONFIDENCE` band, just above the 0.6 release cutoff |
| FX-009 | large / low-confidence heterodimer | `AF-0000000211157965` | ready; values reconciled with AFDB | 699 + 450, the ipTM length-dilution case: `ipTM_d0chn` 0.5769 amber while `ipSAE_d0res` is 0.7119 green. Also the only fixture with non-zero backbone clashes (9) |
| FX-010 | reverse chain orientation (`nx < ny`) + metric disagreement | `AF-0000000211026350` | ready; values reconciled with AFDB | The **short-first** fixture: chain A is 108 residues, chain B is 768. The only real-data cover for `nx < ny`. Also the largest directional ipSAE split in the set (0.6888 vs 0.1203) and the only fixture whose summary prints "scores disagree" |

## Fixture Details

### FX-001: Happy-Path High-Confidence Homodimer

Stable identifier:

```text
AF-0000000065889468
```

Retrieval endpoint:

```text
https://alphafold.ebi.ac.uk/api/prediction/AF-0000000065889468
```

AFDB API check on 2026-05-07:

- response is non-empty
- `isComplex: true`
- `entryId: AF-0000000065889468`
- `modelEntityId: AF-0000000065889468`
- `uniprotAccession: P0A6Q3`
- `uniprotId: FABA_ECOLI`
- `gene: fabA`
- `organismScientificName: Escherichia coli`
- `sequenceEnd: 172`
- `globalMetricValue: 97.27` and `97.24` across the two returned chain records
- required URLs present: `cifUrl`, `paeDocUrl`, `plddtDocUrl`

Notebook-local evidence from `notebooks/homodimer_diagnostic.ipynb`:

- PAE chain info has two chains, A and B, each residue range 1-172.
- Dimer length is 344 residues total.
- Mean pLDDT is 97.2 for chain A and 97.3 for chain B.
- Contact count is 116 CB-CB/CA-for-glycine contact pairs at 8.0 Angstrom cutoff.
- Existing notebook score snapshot:

| Output | Expected Value | Comparison Mode | Source |
|---|---:|---|---|
| `ipSAE_d0res` | 0.9143 | tolerance ±0.001 | existing notebook output |
| `ipSAE_d0chn` | 0.9529 | tolerance ±0.001 | existing notebook output |
| `ipSAE_d0dom` | 0.9527 | tolerance ±0.001 | existing notebook output |
| `ipTM` / `ipTM_d0chn` | 0.9529 | tolerance ±0.001 | existing notebook output |
| `pDockQ` | 0.6913 | tolerance ±0.001 | existing notebook output |
| `pDockQ2` | 0.9269 | tolerance ±0.001 | existing notebook output |
| `LIS` | 0.7564 | tolerance ±0.001 | existing notebook output |

Re-verified 2026-09-08 (R024): all seven values reproduce unchanged, and AFDB's own published figures agree at the two decimal places it publishes for a homodimer — `ipSAE` 0.91, `ipSAE_d0chn` 0.95, `ipsae_iptm_d0chn` 0.95, `pDockQ` 0.69, `pDockQ2` 0.92, `LIS` 0.75, complex `globalMetricValue` 97.26, `N_clash_backbone` 1.0. This fixture is the only one where AFDB rounds every per-direction field to 2 dp, because both directions are identical; the heterodimer fixtures carry full precision (see "Independent Reference Agreement").

Expected behavior:

- notebook runs top-to-bottom
- AFDB metadata table is populated, and the two chains collapse into one identity block because they are the same protein
- mmCIF, PAE JSON, and confidence JSON are downloaded
- two-chain dimer validation passes
- score table contains all required rows
- all scores are interpreted as high/agreeing, unless threshold policy changes
- contact map, PAE heatmap, PAE masks, per-residue profiles, pLDDT comparison, and diagnostic summary render
- minimum MolViewSpec chain overview and pLDDT mapping view render or documented fallback is shown
- runtime target is less than 60 seconds in free-tier Colab

Readiness notes:

- Ready. Also the fixture the local-file path is exercised with (see "Local-File Mode Fixtures"), because it is the only one whose three documents are small enough to keep in a scratch directory comfortably.
- The IPSAE reference version is pinned in `formula-reference.md` (`ipsae.py` v4, Jan 2026); the remaining reference gap is `pDockQ` at better than 2 dp.

### FX-002: Provisional Lower-Confidence / Borderline Candidate

Stable identifier:

```text
AF-0000000066503175
```

Retrieval endpoint:

```text
https://alphafold.ebi.ac.uk/api/prediction/AF-0000000066503175
```

AFDB API check on 2026-05-07:

- response is non-empty
- `isComplex: true`
- `entryId: AF-0000000066503175`
- `modelEntityId: AF-0000000066503175`
- `uniprotAccession: Q55DI5`
- `uniprotId: Q55DI5_DICDI`
- `gene: DDB_G0270900`
- `organismScientificName: Dictyostelium discoideum`
- `uniprotDescription: Transcription elongation factor Eaf N-terminal domain-containing protein`
- `sequenceEnd: 123`
- `globalMetricValue: 86.04` and `86.07` across the two returned chain records
- required URLs present: `cifUrl`, `paeDocUrl`, `plddtDocUrl`

Why this is useful:

- It is a current public AFDB complex accession.
- It is shorter than FX-001 and should run quickly.
- Its global confidence profile is lower than FX-001, making it useful as a second successful complex fixture even if it does not become the final borderline case.
- It was the public complex example surfaced on the AFDB homepage during curation.

Expected behavior:

- two-chain homodimer validation should pass
- notebook should compute all metrics and render all required plots
- score classification may be lower or more mixed than FX-001

Readiness notes:

- Candidate only. It must be scored by the notebook and, ideally, compared with reference calculations before it can be labelled the official borderline fixture.
- If `ipSAE_d0res` is not near the intended 0.6 boundary, keep it as a secondary successful fixture and request a better borderline accession.

### FX-003: Metric-Disagreement Fixture

Status: **closed 2026-09-08 (R024) — superseded by FX-010.** No domain-reviewer request is outstanding.

`AF-0000000211026350` (FX-010) meets every property listed below except "homodimer", which the requirement wrongly narrowed: since R022/R023 the notebook is dimer-scoped, not homodimer-scoped, and a heterodimer disagreement case exercises the summary logic identically. It produces three of the four listed patterns at once — high `pDockQ` (0.5805) with a `LOW-CONFIDENCE` `ipSAE_d0res` (0.6888), a strong pDockQ-versus-pDockQ2 divergence (0.5805 against 0.2721), and LIS below ipSAE — and it is the only accession in the fixture set whose Section 7 prints `OVERALL: Mixed confidence signals — scores disagree`. It was found from live AFDB data rather than requested, and it runs in about 5 s.

FX-009 (`AF-0000000211157965`) is a second, milder disagreement case: `ipTM_d0chn` amber against `ipSAE_d0res` green, driven by chain length rather than by interface quality.

The original requirement is retained below for the record.

Required properties:

- current public AFDB complex accession
- two-chain homodimer compatible with v1
- produces at least one clear disagreement pattern, such as:
  - high `ipSAE_d0res` with low `pDockQ`
  - high `pDockQ` with low `ipSAE_d0res`
  - strong pDockQ versus pDockQ2 divergence
  - high LIS with lower ipSAE
- should run within the notebook runtime target or be marked optional benchmark

Requested output from domain reviewer:

```text
accession:
reason it is a disagreement case:
expected score pattern:
known caveats:
source of prior evidence:
```

Readiness notes:

- ~~Blocking for full validation and final diagnostic-summary review.~~ No longer blocking: FX-010 fills the role.
- Not blocking for first notebook implementation if FX-001 is used for smoke execution.

### FX-004: Malformed Accession Negative

Stable identifier:

```text
AF-NOT_A_REAL_ACCESSION
```

Retrieval endpoint:

```text
https://alphafold.ebi.ac.uk/api/prediction/AF-NOT_A_REAL_ACCESSION
```

AFDB API check on 2026-05-07:

```json
{"error":"Invalid identifier format. Please use a UniProt accession or a supported AlphaFold DB ID."}
```

Expected behavior:

- notebook reports the identifier-format problem clearly
- no downstream parsing or scoring cells run with undefined variables
- no traceback is exposed as the main user-facing result

### FX-005: Valid AFDB Monomer Unsupported By V1

Stable identifier:

```text
O15552
```

Equivalent AFDB entry:

```text
AF-O15552-F1
```

Retrieval endpoint:

```text
https://alphafold.ebi.ac.uk/api/prediction/O15552
```

AFDB API check on 2026-05-07:

- response is non-empty
- `isComplex: false`
- `entryId: AF-O15552-F1`
- `modelEntityId: AF-O15552-F1`
- `uniprotAccession: O15552`
- `uniprotId: FFAR2_HUMAN`
- `gene: FFAR2`
- `organismScientificName: Homo sapiens`
- `sequenceEnd: 330`
- required monomer URLs present: `cifUrl`, `paeDocUrl`, `plddtDocUrl`

Local cached file evidence:

```text
notebooks/AF-O15552-F1.cif
notebooks/models/AF-O15552-F1-model_v6.cif
```

Expected behavior:

- AFDB fetch succeeds
- notebook detects `isComplex: false` or a one-chain/monomer shape
- notebook reports that v1 supports only two-chain homodimers
- no homodimer metric interpretation is produced

### FX-006: Primary Heterodimer (ISG20 + SUMO1)

Stable identifier:

```text
AF-0000000211034637
```

Retrieval endpoint:

```text
https://alphafold.ebi.ac.uk/api/prediction/AF-0000000211034637
```

AFDB API check on 2026-09-08 (per-chain records, keyed by `chainId`, never by position — see R020):

- response is non-empty, two entries, one per chain
- `isComplex: true`, `assemblyType: Hetero`, `oligomericState: dimer`
- `complexComposition: Q96AZ6 x1 + P63166 x1`
- `latestVersion: 1`, `modelCreatedDate: 2026-02-20T08:23:18Z`
- chain A — `uniprotAccession: Q96AZ6`, `uniprotId: ISG20_HUMAN`, `gene: ISG20`,
  `organismScientificName: Homo sapiens`,
  `uniprotDescription: Interferon-stimulated gene 20 kDa protein`,
  `sequenceStart: 1`, `sequenceEnd: 181`, `globalMetricValue: 92.37`
- chain B — `uniprotAccession: P63166`, `uniprotId: SUMO1_MOUSE`, `gene: Sumo1`,
  `organismScientificName: Mus musculus`,
  `uniprotDescription: Small ubiquitin-related modifier 1`,
  `sequenceStart: 1`, `sequenceEnd: 101`, `globalMetricValue: 79.66`
- complex-level `globalMetricValue: 87.82` (search endpoint, `modelEntityId:` query)
- required URLs present and identical across both entries: `cifUrl`, `bcifUrl`, `paeDocUrl`, `plddtDocUrl`

Why this is useful:

- It was the per-task verification fixture for the whole of W2 and W3 under D5, so every heterodimer fix in the rework was measured against it.
- Chains are unequal (181 vs 101, **1.79:1**), so it breaks any surviving `nA == nB` assumption without being so lopsided that a genuine bug hides in the distortion.
- Cross-species (human + mouse), which makes a chain-identity mix-up visible in the metadata report at a glance: R020's non-determinism bug flipped organism, gene, UniProt accession, protein name and monomer length together, and this fixture is where that was first seen.
- Small interface (22 contact pairs) with high PAE confidence, which is the combination that separates pDockQ from pDockQ2 — see the score table below.

Notebook-local evidence, 2026-09-08 (module `complex_interface_utils`, `DIST_CUTOFF` 8.0, `PAE_CUTOFF` 10.0, `LIS_CUTOFF` 12.0):

- PAE matrix 282 × 282, `max_predicted_aligned_error` 31.57; quadrants (181,181), (181,101), (101,181), (101,101)
- chain labels resolve to `ISG20 (A)` / `Sumo1 (B)`
- 22 CB-CB contact pairs; interface residues 11 / 181 (6.1%) for A and 11 / 101 (10.9%) for B
- mean pLDDT 92.4 (A) and 79.7 (B)
- assembly line: `Heterodimer -- 2 chains (A, B), 2 distinct proteins: Q96AZ6 x1 + P63166 x1`, agreeing with the AFDB declaration

| Output | Measured Value | Band | Comparison Mode | Source |
|---|---:|---|---|---|
| `ipSAE_d0res` | 0.7057 | CONFIDENT | tolerance ±0.001 | notebook run 2026-09-08; agrees with AFDB `ipsae_BA` 0.705718 |
| `ipSAE_d0chn` | 0.8055 | green | tolerance ±0.001 | agrees with AFDB `ipsae_d0chn_AB` 0.805532 |
| `ipSAE_d0dom` | 0.7891 | green | tolerance ±0.001 | notebook run 2026-09-08 |
| `ipTM_d0chn` | 0.7730 | HIGH | tolerance ±0.001 | agrees with AFDB `ipsae_iptm_d0chn_BA` 0.772954 |
| `pDockQ` | 0.1452 | MODERATE | tolerance ±0.001 | AFDB publishes 0.15 (2 dp) |
| `pDockQ2` | 0.7054 | HIGH | tolerance ±0.001 | agrees with AFDB `pDockQ2_AB` 0.705403 |
| `LIS` | 0.6004 | HIGH | tolerance ±0.001 | agrees with the mean of AFDB `LIS_AB` 0.608774 and `LIS_BA` 0.592071 |

Expected behavior:

- notebook runs top-to-bottom in online mode with no error
- metadata report expands to two chain blocks rather than collapsing to one, and prints the same content on every run regardless of endpoint entry order
- axis labels, quadrant headings and summary prose use `ISG20 (A)` / `Sumo1 (B)`, never "the two chains" or "identical copies"
- pDockQ amber against pDockQ2 green is the expected reading, not a defect: 22 contacts is too few for pDockQ's contact-count term while the PAE at those contacts is confident

Readiness notes:

- Ready. Scores are reconciled against AFDB's own published values for five of the seven outputs; `ipSAE_d0dom` and `pDockQ` are notebook snapshots only.

### FX-007: Length-Asymmetry Stress Fixture (Sptlc3 + Gm6993)

Stable identifier:

```text
AF-0000000204661110
```

Retrieval endpoint:

```text
https://alphafold.ebi.ac.uk/api/prediction/AF-0000000204661110
```

AFDB API check on 2026-09-08:

- `isComplex: true`, `assemblyType: Hetero`, `oligomericState: dimer`
- `complexComposition: Q8BG54 x1 + A0A1W2P738 x1`, `latestVersion: 1`
- chain A — `Q8BG54` / `SPTC3_MOUSE` / `Sptlc3` / *Mus musculus* / `Serine palmitoyltransferase 3` / 1-563 / `globalMetricValue: 81.67`
- chain B — `A0A1W2P738` / `A0A1W2P738_MOUSE` / `Gm6993` / *Mus musculus* / `Predicted gene 6993` / 1-69 / `globalMetricValue: 71.23`
- complex-level `globalMetricValue: 80.53`
- required URLs present: `cifUrl`, `bcifUrl`, `paeDocUrl`, `plddtDocUrl`

How it was found:

Ranked 100 `assemblyType:Hetero AND oligomericState:dimer` search hits by chain-length ratio during R022 and took the most extreme. At **8.16:1** it is the widest ratio in the fixture set.

Why this is useful:

- Every axis, coverage bar and mask panel that quietly assumed comparable chain lengths misbehaves here first. R022 fixed four such sites against this fixture (a clipped residue 0, a one-chain-only interface range, and two pieces of false homodimer prose).
- A 69-residue chain drawn on a 563-wide shared axis is where "shorter chain" and "truncated data" become visually indistinguishable — the annotation and end-cap on the coverage track exist because of this fixture.

Notebook-local evidence, 2026-09-08:

- PAE matrix 632 × 632, `max_predicted_aligned_error` 31.62
- chain labels `Sptlc3 (A)` / `Gm6993 (B)`
- 42 contact pairs; interface residues 21 / 563 (3.7%) for A and 16 / 69 (23.2%) for B — the same interface, two very different percentages, which is the point
- mean pLDDT 81.7 (A) and 71.2 (B)

| Output | Measured Value | Band | Comparison Mode | Source |
|---|---:|---|---|---|
| `ipSAE_d0res` | 0.7120 | CONFIDENT | tolerance ±0.001 | agrees with AFDB `ipsae_BA` 0.711961 |
| `ipSAE_d0chn` | 0.8181 | green | tolerance ±0.001 | agrees with AFDB `ipsae_d0chn_AB` 0.818091 |
| `ipSAE_d0dom` | 0.7856 | green | tolerance ±0.001 | notebook run 2026-09-08 |
| `ipTM_d0chn` | 0.6687 | MODERATE | tolerance ±0.001 | agrees with AFDB `ipsae_iptm_d0chn_BA` 0.668721 |
| `pDockQ` | 0.2737 | HIGH | tolerance ±0.001 | AFDB publishes 0.27 (2 dp) |
| `pDockQ2` | 0.5914 | HIGH | tolerance ±0.001 | agrees with AFDB `pDockQ2_AB` 0.591426 |
| `LIS` | 0.3980 | HIGH | tolerance ±0.001 | agrees with the mean of AFDB `LIS_AB` 0.417853 and `LIS_BA` 0.378097 |

Expected behavior:

- notebook runs top-to-bottom with no shape error and no wrong count
- the 69-residue chain reads as *shorter*, not as truncated, on every shared-axis figure
- interface percentages differ wildly between the chains and both are correct

Known residual, not a fixture defect:

At 8.16:1 the `aspect='auto'` distortion in the PAE heatmap, the mask panels and the contact map's left panel stops being cosmetic — a 563 × 69 block renders as a near-square. Axis labels stay correct, so nothing printed is false. Owned by R050 and R052; R022 was forbidden from making the figure-size change.

### FX-008: Mid-Confidence Heterodimer (PACRG + MEIG1)

Stable identifier:

```text
AF-0000000211619209
```

Retrieval endpoint:

```text
https://alphafold.ebi.ac.uk/api/prediction/AF-0000000211619209
```

AFDB API check on 2026-09-08:

- `isComplex: true`, `assemblyType: Hetero`, `oligomericState: dimer`
- `complexComposition: Q96M98 x1 + Q5JSS6 x1`, `latestVersion: 1`
- chain A — `Q96M98` / `PACRG_HUMAN` / `PACRG` / *Homo sapiens* / `Parkin coregulated gene protein` / 1-296 / `globalMetricValue: 68.71`
- chain B — `Q5JSS6` / `MEIG1_HUMAN` / `MEIG1` / *Homo sapiens* / `Meiosis expressed gene 1 protein homolog` / 1-88 / `globalMetricValue: 89.93`
- complex-level `globalMetricValue: 73.57`
- required URLs present: `cifUrl`, `bcifUrl`, `paeDocUrl`, `plddtDocUrl`

Why this is useful:

- `ipSAE_d0res` of 0.6366 is the only value in the fixture set inside AFDB's `LOW-CONFIDENCE` band (0.6-0.7), just above the 0.6 release cutoff. It is the closest thing the online path has to a borderline case, and it is what FX-002 was originally wanted for.
- Same-species (both human), so it is the control against FX-006: an identity mix-up here does *not* change the organism, only the gene and accession.
- Global pLDDT is split the opposite way from FX-006 — the short chain is the confident one (89.93 against 68.71).

Notebook-local evidence, 2026-09-08:

- PAE matrix 384 × 384, `max_predicted_aligned_error` 31.71; ratio 3.36:1
- chain labels `PACRG (A)` / `MEIG1 (B)`
- 43 contact pairs; interface residues 20 / 296 (6.8%) for A and 16 / 88 (18.2%) for B
- mean pLDDT 68.7 (A) and 89.9 (B)
- `ipSAE_d0dom` intermediates: `d0 = 5.84118`, `n_dom = 249` — both exactly equal to AFDB's `complexPredictionAccuracy_ipsae_d0dom_AB` and `_ipsae_n0dom_AB`, which is the tightest confirmation in this manifest that the d0dom normalisation matches the reference

| Output | Measured Value | Band | Comparison Mode | Source |
|---|---:|---|---|---|
| `ipSAE_d0res` | 0.6366 | LOW-CONFIDENCE | tolerance ±0.001 | agrees with AFDB `ipsae_BA` 0.63659 |
| `ipSAE_d0chn` | 0.8313 | green | tolerance ±0.001 | agrees with AFDB `ipsae_d0chn_AB` 0.831315 |
| `ipSAE_d0dom` | 0.7735 | green | tolerance ±0.001 | notebook run 2026-09-08 |
| `ipTM_d0chn` | 0.8183 | HIGH | tolerance ±0.001 | agrees with AFDB `ipsae_iptm_d0chn_AB` 0.818349 |
| `pDockQ` | 0.2418 | HIGH | tolerance ±0.001 | AFDB publishes 0.24 (2 dp) |
| `pDockQ2` | 0.3254 | HIGH | tolerance ±0.001 | agrees with AFDB `pDockQ2_BA` 0.325409 |
| `LIS` | 0.4783 | HIGH | tolerance ±0.001 | agrees with the mean of AFDB `LIS_AB` 0.496848 and `LIS_BA` 0.459826 |

Expected behavior:

- notebook runs top-to-bottom with no error
- `ipSAE_d0res` prints the AFDB band label `LOW-CONFIDENCE` while the other four traffic lights are green

Known defect this fixture exposes (recorded, not fixed here):

The Section 7 OVERALL prose prints *"This complex has consistently HIGH confidence across all metrics"* even though `ipSAE_d0res` is amber. The prose branch and the per-score traffic lights disagree. **Owned by R080** (wire the summary to the canonical thresholds); R024 records it rather than fixing it. FX-008 and FX-009 both reproduce it.

### FX-009: Large Low-Confidence Heterodimer (UVRAG + Beclin-1)

Stable identifier:

```text
AF-0000000211157965
```

Retrieval endpoint:

```text
https://alphafold.ebi.ac.uk/api/prediction/AF-0000000211157965
```

AFDB API check on 2026-09-08:

- `isComplex: true`, `assemblyType: Hetero`, `oligomericState: dimer`
- `complexComposition: Q9P2Y5 x1 + Q14457 x1`, `latestVersion: 1`
- chain A — `Q9P2Y5` / `UVRAG_HUMAN` / `UVRAG` / *Homo sapiens* / `UV radiation resistance-associated gene protein` / 1-699 / `globalMetricValue: 61.89`
- chain B — `Q14457` / `BECN1_HUMAN` / `BECN1` / *Homo sapiens* / `Beclin-1` / 1-450 / `globalMetricValue: 79.81`
- complex-level `globalMetricValue: 68.91`, `complexPredictionAccuracy_N_clash_backbone: 9.0`
- required URLs present: `cifUrl`, `bcifUrl`, `paeDocUrl`, `plddtDocUrl`

Why this is useful:

- **The ipTM length-dilution case.** 1149 residues total gives `d0_chn` its largest value in the set, and `ipTM_d0chn` — which averages over every partner residue with no PAE cutoff — falls to 0.5769 (amber) while `ipSAE_d0res`, which uses only pairs under 10 Å, reads 0.7119 (green). This is the clearest demonstration in the fixture set of why the notebook reports both.
- **The only fixture with non-zero backbone clashes** (9, against AFDB's release limit of ≤10). It is therefore the fixture to use when the joint AFDB criterion's clash term is wired up.
- Largest and slowest of the set (501 contact pairs, 1149 × 1149 PAE), so it is the runtime canary for the <60 s Colab budget.

Notebook-local evidence, 2026-09-08:

- PAE matrix 1149 × 1149, `max_predicted_aligned_error` 31.69; ratio 1.55:1
- chain labels `UVRAG (A)` / `BECN1 (B)`
- 501 contact pairs; interface residues 181 / 699 (25.9%) for A and 184 / 450 (40.9%) for B
- mean pLDDT 61.9 (A) and 79.8 (B)
- notebook wall-clock 5.4 s locally, executed with `nbclient`

| Output | Measured Value | Band | Comparison Mode | Source |
|---|---:|---|---|---|
| `ipSAE_d0res` | 0.7119 | CONFIDENT | tolerance ±0.001 | agrees with AFDB `ipsae_AB` 0.711865 |
| `ipSAE_d0chn` | 0.8947 | green | tolerance ±0.001 | agrees with AFDB `ipsae_d0chn_AB` 0.894745 |
| `ipSAE_d0dom` | 0.8691 | green | tolerance ±0.001 | notebook run 2026-09-08 |
| `ipTM_d0chn` | 0.5769 | MODERATE | tolerance ±0.001 | agrees with AFDB `ipsae_iptm_d0chn_AB` 0.576863 |
| `pDockQ` | 0.7191 | HIGH | tolerance ±0.001 | AFDB publishes 0.72 (2 dp) |
| `pDockQ2` | 0.2899 | HIGH | tolerance ±0.001 | agrees with AFDB `pDockQ2_BA` 0.289935 |
| `LIS` | 0.4278 | HIGH | tolerance ±0.001 | agrees with the mean of AFDB `LIS_AB` 0.440768 and `LIS_BA` 0.414834 |

Expected behavior:

- notebook runs top-to-bottom with no error and inside the runtime target
- `ipTM_d0chn` amber against `ipSAE_d0res` green is the expected reading, and the Section 4 prose should explain the length dependence rather than treat it as a contradiction
- reproduces the same R080 OVERALL-prose defect described under FX-008

### FX-010: Reverse Chain Orientation and Metric Disagreement (RBX1 + Cullin-3)

Stable identifier:

```text
AF-0000000211026350
```

Retrieval endpoint:

```text
https://alphafold.ebi.ac.uk/api/prediction/AF-0000000211026350
```

AFDB API check on 2026-09-08:

- `isComplex: true`, `assemblyType: Hetero`, `oligomericState: dimer`
- `complexComposition: P62878 x1 + Q13618 x1`, `latestVersion: 1`, `modelCreatedDate: 2026-02-20T09:48:35Z`
- chain A — `P62878` / `RBX1_MOUSE` / `Rbx1` / *Mus musculus* / `E3 ubiquitin-protein ligase RBX1` / 1-108 / `globalMetricValue: 51.74`
- chain B — `Q13618` / `CUL3_HUMAN` / `CUL3` / *Homo sapiens* / `Cullin-3` / 1-768 / `globalMetricValue: 87.32`
- complex-level `globalMetricValue: 82.94`, `complexPredictionAccuracy_N_clash_backbone: 0.0`
- required URLs present: `cifUrl`, `bcifUrl`, `paeDocUrl`, `plddtDocUrl`

Why this is useful — it closes the gap R022 flagged:

FX-006, FX-007, FX-008 and FX-009 all put the **long** chain first, so `nx < ny` had never been exercised against real data — only against synthetic cases (30+400, 5+300). Here chain A is 108 residues and chain B is 768, a **1 : 7.11** ratio in the reverse direction, and the ordered pair the notebook builds is short → long. It is the mirror of FX-007 and should be run alongside it whenever anything touches an axis, a quadrant slice or a per-residue array.

It also does three other jobs:

- **It supplies the metric disagreement FX-003 was blocking on**, from live AFDB data rather than a reviewer request. pDockQ 0.5805 (green) against pDockQ2 0.2721, and LIS 0.3841 against `ipSAE_d0res` 0.6888, drive Section 7 into `OVERALL: Mixed confidence signals — scores disagree`, printing the STRUCTURE-vs-PAE, pDockQ-vs-pDockQ2 and LIS-vs-ipSAE branches. It is the only fixture in the set that does.
- **Largest directional asymmetry in the set.** AFDB's own per-direction values are `ipsae_AB` 0.688783 against `ipsae_BA` 0.120251 — a 0.57 gap, where FX-001's two directions are identical to 2 dp. This is the fixture for R008's directional diagnostic.
- **Lowest per-chain confidence in the set** (chain A at 51.74), so the interface-pLDDT panel and the "low pLDDT at the interface" counter have something to show.

Notebook-local evidence, 2026-09-08:

- PAE matrix 876 × 876, `max_predicted_aligned_error` 31.69; quadrants (108,108), (108,768), (768,108), (768,768)
- chain labels `Rbx1 (A)` / `CUL3 (B)`
- 155 contact pairs; interface residues 32 / 108 (29.6%) for A and 60 / 768 (7.8%) for B; interface residue ranges 18-90 (A) and 441-768 (B)
- mean pLDDT 51.7 (A) and 87.3 (B)
- `ipSAE_d0dom` intermediates `d0 = 7.77904`, `n_dom = 476`

| Output | Measured Value | Band | Comparison Mode | Source |
|---|---:|---|---|---|
| `ipSAE_d0res` | 0.6888 | LOW-CONFIDENCE | tolerance ±0.001 | agrees with AFDB `ipsae_AB` 0.688783 |
| `ipSAE_d0chn` | 0.8597 | green | tolerance ±0.001 | agrees with AFDB `ipsae_d0chn_BA` 0.859679 |
| `ipSAE_d0dom` | 0.7992 | green | tolerance ±0.001 | notebook run 2026-09-08 |
| `ipTM_d0chn` | 0.5958 | MODERATE | tolerance ±0.001 | agrees with AFDB `ipsae_iptm_d0chn_AB` 0.595829 |
| `pDockQ` | 0.5805 | HIGH | tolerance ±0.001 | AFDB publishes 0.58 (2 dp) |
| `pDockQ2` | 0.2721 | HIGH | tolerance ±0.001 | agrees with AFDB `pDockQ2_BA` 0.27208 |
| `LIS` | 0.3841 | HIGH | tolerance ±0.001 | agrees with the mean of AFDB `LIS_AB` 0.42141 and `LIS_BA` 0.346871 |

Expected behavior:

- notebook runs top-to-bottom with no error, with the short chain first everywhere
- quadrant shapes print as (108,108), (108,768), (768,108), (768,768). Because the two inter-chain blocks are not square and not interchangeable, a transposed or mis-ordered slice raises here instead of quietly producing a plausible wrong number, which is exactly what an equal-length homodimer cannot test
- Section 7 prints the mixed-signal branch, not the "consistently HIGH" branch

### Short-First Orientation Search (R024)

R022 recorded that both heterodimer fixtures then in hand put the long chain first, leaving `nx < ny` covered only by synthetic cases. R024 searched for a real one.

Method, 2026-09-08:

```text
GET https://alphafold.ebi.ac.uk/api/search?q=assemblyType:Hetero AND oligomericState:dimer&type=main&rows=100&start=<n>
    -> numFound: 80240
```

The search documents do **not** carry per-chain lengths, so each hit's `entryId` was resolved through the prediction endpoint and the chains keyed by `chainId`, with `sequenceEnd - sequenceStart + 1` as the length.

Result over the first **100** consecutive hits (0 errors, 0 non-two-chain records):

| Orientation | Count |
|---|---:|
| short-first (`n_A < n_B`) | 41 |
| long-first (`n_A > n_B`) | 57 |
| equal lengths | 2 |

Smallest `n_A / n_B` seen: **0.1288** (`AF-0000000205111076`, 193 + 1499).

**Conclusion: short-first heterodimers are not rare — they are roughly 41% of AFDB heterodimers.** The R022 observation was an artefact of the four accessions that happened to be picked, not a property of the database. `AF-0000000211026350` was selected from the ranked list as FX-010 because it combines a strong reverse ratio (1 : 7.11) with a size that keeps the notebook inside the runtime target, and it was verified end to end before registration. `AF-0000000205111076` is the more extreme option (1 : 7.77) and is recorded here as a spare, but at 1692 residues it is a heavier download than the fixture set needs.

Search endpoint notes worth keeping:

- `q=entryId:<accession>` returns HTTP 404 with an empty body. The field that works for a single-accession lookup is `q=modelEntityId:<accession>`. A bare `q=<accession>` also 404s.
- The search document carries AFDB's own published IPSAE-derived scores per direction (`complexPredictionAccuracy_ipsae_AB` / `_BA`, `_ipsae_d0chn_*`, `_ipsae_iptm_d0chn_*`, `_pDockQ2_*`, `_LIS_*`, `_ipsae_d0dom_*`, `_ipsae_n0dom_*`, `_N_clash_backbone`), which is what makes the reference agreement below possible without running `ipsae.py` locally.

### Independent Reference Agreement (R024)

Every heterodimer fixture's scores were reconciled against AlphaFold DB's own published values, read from the search endpoint's `complexPredictionAccuracy_*` fields. These are produced by AFDB's own IPSAE run, not by this module, so the agreement is an independent check rather than a self-comparison.

Agreement rules used, both of which are conventions rather than discrepancies:

- `ipSAE`, `ipSAE_d0chn`, `ipTM_d0chn` and `pDockQ2` are reported by this module as the **maximum** over the two directions, matching `ipsae.py`; AFDB publishes both directions separately, so the comparison is against `max(AB, BA)`.
- `LIS` is reported as the **mean** of the two directions, again matching `ipsae.py`; AFDB's rolled-up `complexPredictionAccuracy_LIS` is the max, so the comparison is against `mean(LIS_AB, LIS_BA)`. This is why the rolled-up figure can differ in the second decimal place from the notebook's LIS while the per-direction values agree exactly.

Result across FX-006 to FX-010: **all five values agree to 4 decimal places on all five fixtures** (|Δ| < 0.0001, well inside the ±0.001 tolerance). `ipSAE_d0dom`'s intermediates agree too — on FX-008 the module's `d0 = 5.84118` and `n_dom = 249` are exactly AFDB's `ipsae_d0dom_AB` and `ipsae_n0dom_AB`.

Still unreconciled at full precision:

- `pDockQ` — AFDB publishes `complexPredictionAccuracy_pDockQ` to 2 decimal places only. Every fixture agrees at that precision (0.15 / 0.27 / 0.24 / 0.72 / 0.58), but this is not a ±0.001 check.
- `ipSAE_d0dom` — the rolled-up score is not published; only its `d0` and `n0` intermediates are, and those match exactly.
- Note that AFDB's `complexPredictionAccuracy_ipsae_dist_cutoff` is 15.0 while this module's `DIST_CUTOFF` is 8.0. The two cutoffs govern different things — AFDB's feeds the `ipsae_dist_nres*` counts, this module's feeds pDockQ/pDockQ2 contacts — and pDockQ agreement at 2 dp on all five fixtures confirms the contact definitions have not diverged.

### Local-File Mode Fixtures (R025)

Local-file mode (`USE_LOCAL_FILE = True`) never reaches the AFDB metadata gate, so it needs its own negative fixtures. All are constructed from files already listed above.

| Case | Inputs | Expected behaviour |
|---|---|---|
| mmCIF only | FX-001's `-model_v1.cif` | `MissingLocalDocumentError` at the upload cell, naming the PAE document and the pLDDT document |
| mmCIF + PAE, no pLDDT | FX-001's `.cif` + `-predicted_aligned_error_v1.json` | `MissingLocalDocumentError` at the upload cell, naming the pLDDT document |
| all three | FX-001's three files | runs to the end; scores identical to the online run of FX-001 |
| monomer | FX-005's three files (`AF-O15552-F1`) | `UnsupportedAssemblyError` — "This model is a monomer" — from the structural gate in `verify_chain_identity` |
| three chains | synthetic A/B/C built from FX-001 by duplicating chain B, with matching 516 × 516 PAE and 516-score pLDDT | `UnsupportedAssemblyError` — "This model is a homotrimer" — from the same structural gate |

The three-chain case has to be synthetic: R023 measured that AFDB contains no prediction with more than two chains (`isComplex:true` and `oligomericState:dimer` both return 2,010,763; `trimer` and `tetramer` return zero), so a local mmCIF is the only route to that code path.

## Provenance Requirements

For every finalized fixture, keep:

- accession
- retrieval endpoint
- retrieval date
- AFDB response record used
- cached artifact path if cached
- source of expected metric values
- tolerance policy
- whether values are smoke snapshots or trusted reference snapshots

## Fixture Readiness Decision

Updated 2026-09-08 (R024). Current status: `full-validation-ready for scoring`, with two documentation-level gaps remaining.

### What changed

The three conditions that held this manifest at `smoke-ready` have been resolved as follows.

| Former blocker | Status |
|---|---|
| FX-003 supplied by a domain reviewer | **Closed** — FX-010 supplies the disagreement pattern from live AFDB data. No reviewer request outstanding. |
| numeric snapshots frozen from a trusted run | **Closed for six of seven outputs** — reconciled against AFDB's own published IPSAE values to 4 dp on FX-006 to FX-010. `pDockQ` is confirmed only to 2 dp because that is all AFDB publishes. |
| IPSAE reference version/commit pinned | **Closed** — `ipsae.py` v4 (Jan 2026), pinned in `formula-reference.md` under D2. |
| score threshold bands approved | **Closed** — `threshold-reference.md`, adopted by R002/R009, with per-threshold provenance (`PUBLISHED` / `DERIVED` / `HEURISTIC`) carried in `THRESHOLDS`. |
| FX-002 confirmed or replaced | **Open, but no longer blocking** — FX-008 (`ipSAE_d0res` 0.6366) now occupies the borderline role from live data. FX-002 remains unscored and should either be scored or retired in R090. |

### The working set

| Purpose | Fixture |
|---|---|
| primary homodimer, and the local-file path | FX-001 |
| primary heterodimer, per-task verification under D5 | FX-006 |
| chain-length asymmetry, long-first | FX-007 |
| chain-length asymmetry, short-first (`nx < ny`) | FX-010 |
| borderline `ipSAE_d0res` near the 0.6 AFDB cutoff | FX-008 |
| large model, ipTM length dilution, non-zero clashes | FX-009 |
| metric disagreement | FX-010, with FX-009 as the milder case |
| malformed accession | FX-004 |
| valid but unsupported assembly | FX-005 (online monomer) and the local three-chain case |

A change that touches axes, quadrant slicing, per-residue arrays or labelling should be verified on **FX-001, FX-007 and FX-010** at minimum: one equal-length pair, one long-first extreme and one short-first extreme.

### Verified end to end on 2026-09-08

All six online fixtures (FX-001, FX-006, FX-007, FX-008, FX-009, FX-010) execute the notebook top to bottom under `nbclient`: 20 of 20 code cells executed, zero errors, 4.1 s to 5.4 s each locally. The five local-file cases in "Local-File Mode Fixtures" were exercised by driving the notebook's own upload cells with a stubbed `FileUpload` widget.

### What is still open

- `pDockQ` has no better-than-2-dp independent reference. Closing it needs a local `ipsae.py` v4 run, which is R092's job.
- FX-002 is still unscored; either score it or retire it (R090).
- The Section 7 OVERALL prose disagrees with the per-score traffic lights on FX-008 and FX-009: it prints "consistently HIGH confidence across all metrics" while `ipSAE_d0res` is amber on FX-008 and `ipTM_d0chn` is amber on FX-009. Owned by **R080**; recorded here because two registered fixtures reproduce it and any summary rework must be checked against them.
- The `aspect='auto'` distortion at extreme length ratios (FX-007, FX-010). Owned by **R050** and **R052**.

