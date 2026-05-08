# Homodimer Diagnostic Fixture Manifest

## Status

Current readiness: `smoke-ready`.

The fixture set is now sufficient to start notebook implementation and smoke validation. It is not sufficient for full validation, beta, or scientific sign-off because a true metric-disagreement fixture and trusted reference score snapshots are still missing.

Last fixture curation pass: 2026-05-07.

## Fixture Summary

| Fixture ID | Role | Accession / Source | Status | Purpose |
|---|---|---|---|---|
| FX-001 | happy-path / high-confidence reference candidate | `AF-0000000065889468` | ready for smoke; reference values provisional | Exercise complete successful AFDB fetch, parse, metric, visualization, MolViewSpec, and summary workflow |
| FX-002 | lower-confidence / provisional borderline candidate | `AF-0000000066503175` | candidate-needs-scoring | Stress a current AFDB complex with lower global confidence than FX-001; may become borderline fixture if ipSAE/summary behavior supports it |
| FX-003 | metric-disagreement | requested from AFDB/PDBe domain reviewer | blocking | Validate summary logic when PAE-based and contact/contact-pLDDT scores disagree |
| FX-004 | malformed negative | `AF-NOT_A_REAL_ACCESSION` | ready | Confirm clear AFDB identifier-format failure and no downstream stack trace |
| FX-005 | valid AFDB but unsupported by v1 | `O15552` / `AF-O15552-F1` | ready | Confirm monomer entries are rejected clearly because v1 supports only two-chain homodimers |

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
| `ipSAE_d0res` | 0.9143 | tolerance, pending reference freeze | existing notebook output |
| `ipSAE_d0chn` | 0.9529 | tolerance, pending reference freeze | existing notebook output |
| `ipSAE_d0dom` | 0.9527 | tolerance, pending reference freeze | existing notebook output |
| `ipTM` / `ipTM_d0chn` | 0.9529 | tolerance, pending reference freeze | existing notebook output |
| `pDockQ` | 0.6913 | tolerance, pending reference freeze | existing notebook output |
| `pDockQ2` | 0.9269 | tolerance, pending reference freeze | existing notebook output |
| `LIS` | 0.7564 | tolerance, pending reference freeze | existing notebook output |

Expected behavior:

- notebook runs top-to-bottom
- AFDB metadata table is populated
- mmCIF, PAE JSON, and confidence JSON are downloaded
- two-chain homodimer validation passes
- score table contains all required rows
- all scores are interpreted as high/agreeing, unless threshold policy changes
- contact map, PAE heatmap, PAE masks, per-residue profiles, pLDDT comparison, and diagnostic summary render
- minimum MolViewSpec chain overview and pLDDT mapping view render or documented fallback is shown
- runtime target is less than 60 seconds in free-tier Colab

Readiness notes:

- Ready for smoke validation.
- Not yet ready as a reference fixture until the IPSAE reference version/commit is pinned and scores are regenerated from a trusted reference run.

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

Status: requested from AFDB/PDBe domain reviewer.

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

- Blocking for full validation and final diagnostic-summary review.
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

Current status: `smoke-ready`, not `full-validation-ready`.

Implementation may start using:

- FX-001 as the primary successful smoke fixture
- FX-004 as malformed negative fixture
- FX-005 as valid-but-unsupported fixture

Full validation remains blocked until:

- FX-002 is either confirmed as the borderline fixture or replaced
- FX-003 is supplied by a domain reviewer
- numeric snapshots are frozen from a trusted run
- IPSAE reference version/commit is pinned
- score threshold bands are approved

