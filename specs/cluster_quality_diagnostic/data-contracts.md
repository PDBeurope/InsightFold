# Data Contracts

All contracts below were measured against the live services on 2026-09-21. Where a field's behaviour
was verified, the evidence is named.

## 1. AFDB cluster family API

```
GET https://alphafold.ebi.ac.uk/api/workbench/cluster-family/{uniprot_accession}?cluster_flag={flag}
    flag ∈ {"AFDB50/MMseqs2", "AFDB/Foldseek"}   percent-encoded: AFDB50%2FMMseqs2
```

Public, no key, and **publishable**: this is the advertised route. It replaced an earlier
unadvertised internal route on 2026-10-02. That route is deliberately not named here or anywhere
else in the repository; it was never committed, and it should stay that way.

**Migration verified 2026-10-02, not assumed.** The two routes returned **byte-identical JSON**
across 10 accession-and-flag pairs (`P06965`, `O15552`, `Q13148`, `P04637`, `P00533`, both
clusterings), and identical status codes on every failure path: 404 with the same `detail` string
for absent and malformed accessions, 422 for a bad `cluster_flag`, and the same per-accession 500
on `P12345`. The slash in the flag works encoded or bare; the module encodes it.

**Rate policy: 2.5 s spacing retained, and its justification is now stale.** The 0.8 s failure
storm no longer reproduces on **either** route: 40 distinct accessions at 0.8 s against the new
endpoint returned 37 × 200, 2 × 404 (absent accessions) and the one known `P12345` 500, with no
rate-limit failures at all. The service appears to have relaxed the limit since 2026-09-21. The
spacing is kept because it is cheap insurance and nothing depends on removing it, but it is no
longer a measured requirement and should not be cited as one.

Response shape: `{clusterTotal, clusterMembers: {<seven parallel arrays>}}`.

| Field | Required | Type | Used for | Validation |
|---|---|---|---|---|
| `clusterTotal` | yes | int | Member count; compared against array length | Must equal `len(afdbAccessions)`. Log a warning if not (OQ-3) |
| `afdbAccessions` | yes | list[str] | Member identity, `AF-{ACC}-F1` form | All seven arrays must be the same length |
| `averagePlddt` | yes | list[float] | **Every distribution statistic** | Finite, 0 to 100. Verified sorted descending |
| `sequenceLength` | yes | list[int] | Length figures | Positive integers |
| `uniprotDescriptions` | yes | list[str] | Member labels in figures and the summary | May be empty or "Uncharacterized protein" |
| `speciesNames` | yes | list[str] | **Unused in v1** | Free text. **No taxId and no lineage**, which is why the taxonomic composition figure was dropped: it could say nothing without a secondary lookup |
| `reviewedStatus` | no | list[bool] | Display only | Unused for logic in v1 |
| `referenceLabel` | no | list[bool] | Unused in v1 | Reflects the UniProt release AFDB was built against rather than the current one. Not validated per accession |

**Verified behaviours:**

- Arrays are **sorted descending by `averagePlddt`**. Confirmed on all 13 responses examined, so best
  is index 0 and worst is index −1. REQ-009 depends on this and must assert it rather than assume it.
- Accepts **bare accessions only**. `AF-Q9I1F6-F1` returns 404.
- Isoforms work: `P04637-2` resolves to the parent's cluster.
- Returns the full cluster in one response, verified to 93,793 members at 8.82 MB in 1.0 s. No
  pagination.
- **The query is not reliably a member of the returned cluster.** Absent from the structure cluster
  in 7 of 8 accessions measured. *(An earlier revision attributed this to Foldseek clustering only
  AFDB50 representatives. That explanation is refuted; see the disjoint-member-sets entry below.)*
  On MMseqs2, four accessions
  (`Q9KM69`, `A0A7X7SVK1`, `A0A536HFQ7`, `A0A8C8MD66`) each return a byte-identical response to
  another accession's cluster, most likely because they share the same underlying amino acid
  sequence and resolve to the same cluster. Either way, REQ-006 must treat absence as a normal path
  and never assume membership.
- **Cluster breadth. CORRECTED.** AFDB50 uses **cascaded** clustering, so members are linked to the
  final representative **transitively, through a chain of intermediate representatives**, each link
  within threshold. Identity to the final representative, and between arbitrary members, is
  **unbounded below**; measured local identity to the best member reaches 25 to 27% with coverage
  from 21% to 75%. It is *not* true that every
  member is within threshold of the representative. Steinegger and Söding 2018, Nat Commun 9:2542.
  REQ-E7 requires the notebook to explain this correctly.
- **The two clusterings return disjoint member sets. MEASURED, mechanism unconfirmed.** For all 8
  accessions with both flags pinned, the sequence-cluster and structure-cluster member lists share
  **zero** accessions. The query is a member of exactly one of them, usually the sequence one
  (7 of 8; `P04637-2` is the exception, present in the structure cluster and absent from the
  sequence one). Querying any member of either set returns the same pair of clusters: `P01956`
  and `P0C0U8` both return `P69905`'s 5,223-member and 53-member clusters respectively.
  **An earlier revision claimed the sequence cluster's representative is what appears in the
  structure cluster. That is refuted**: the sequence cluster's top member is absent from the
  structure cluster in 8 of 8 cases. The mechanism is not determinable from the API alone and
  should be confirmed with the AFDB team. This is the same behaviour as OQ-1: the endpoint will
  return a cluster the queried accession does not belong to.

**Failure behaviour:**

| Status | Meaning | Notebook behaviour |
|---|---|---|
| 404 | Absent, malformed, or a fragmented protein | Refuse, quoting `detail` from the service |
| 422 | Invalid `cluster_flag` value | Unreachable from the notebook; the flag comes from `CLUSTER_FLAGS` |
| 5xx | A per-accession fault; rate limiting no longer observed | Retry with backoff before reporting |
| 200, `clusterTotal` < 3 | Too small to analyse | Refuse, naming the size |

**Provenance:** public and advertised as of 2026-10-02. Its predecessor was reachable but
unadvertised, which is why NFR-008 forbade naming the cluster endpoint in prose; that restriction no
longer applies to this route, and the predecessor's path stays unrecorded. The base URL is still
defined exactly once (NFR-008).

## 2. AFDB prediction API

```
GET https://alphafold.ebi.ac.uk/api/prediction/{accession}
```

Public, no key. **Needs no request spacing**: 45 back-to-back requests at 0.0 s produced zero
failures at 0.15 s median latency. Do not apply the cluster endpoint's 2.5 s policy here.

Returns a list with one entry per chain; for a monomer, one entry.

| Field | Required | Used for | Validation |
|---|---|---|---|
| `sequence` / `uniprotSequence` | yes | Pairwise alignment input | Non-empty; length must match `sequenceLength` from the cluster response to within the fragment caveat |
| `uniprotDescription` | yes | First visible confirmation | Printed in C008 |
| `organismScientificName` | yes | First visible confirmation | Printed in C008 |
| `globalMetricValue` | yes | The user's protein's pLDDT | May differ slightly from the cluster endpoint's `averagePlddt` rounding. Record which endpoint each number came from |
| `cifUrl` | yes | Mol\* structure source | Non-empty |
| `bcifUrl` | no | Mol\* alternative | **May be present but empty.** Treat empty as absent; an empty `bcifUrl` labelled `bcif` is a silent blank viewer |
| `msaUrl` | no | MSA coverage | e.g. `files/msa/AF-Q9I1F6-F1-msa_v6.a3m` |
| `latestVersion` | yes | Provenance line | Not `modelVersion`, which does not exist |

Called at most **three times per run**: once for the user's protein, once each for the best and worst
members. Never per cluster member.

**Failure behaviour.** Required, because a substantial share of cluster members are deleted from
UniProtKB (14 of 20 and 9 of 20 in two measured clusters) while their AFDB structures still exist.

| Status / shape | Meaning | Notebook behaviour |
|---|---|---|
| 404 | Accession not served by AFDB | For the user's protein: refuse per REQ-017. For a best or worst member: skip the alignment with a named message; the 3D view may still render if `cifUrl` was already obtained |
| 200 with empty or missing `sequence` | Record exists but carries no sequence | Same as 404 for alignment purposes. **Do not** treat an empty string as a valid sequence |
| 200 with empty `bcifUrl` | Known: the field can be present but empty | Treat as absent and use `cifUrl`. An empty `bcifUrl` labelled `bcif` is a silent blank viewer |
| 5xx | Transient | Retry with backoff. Unlike the cluster endpoint, this is not expected to be rate limiting |

**No real accession is known to trigger the missing-sequence path.** All 20 of `Q9I1F6`'s top members
were checked: 14 are deleted from UniProtKB, yet AFDB serves a sequence for every one. FX-12
`SYN-NOSEQ` is therefore a synthetic response with the sequence fields removed, and the manifest
records that no real case is known rather than implying one exists.

## 2a. AFDB structure and PAE file downloads

```
GET {cifUrl}        # e.g. .../AF-O15552-F1-model_v6.cif
GET {paeDocUrl}     # e.g. .../AF-O15552-F1-predicted_aligned_error_v6.json
```

**Measured 2026-09-23 on `AF-O15552-F1`:**

| Property | Value | Consequence |
|---|---|---|
| `Transfer-Encoding` | **chunked** | The body is streamed |
| `Content-Length` | **absent**, on both GET and HEAD | Truncation **cannot** be detected from the response |
| `Accept-Ranges` | **absent** | A partial download **cannot** be resumed; every retry restarts from zero |
| mmCIF size | 311,399 bytes for a 330-residue protein | Roughly 1 KB per residue |
| PAE size | 260,390 bytes for the same protein | Grows with the square of the length |

**Failure mode observed in the field.** A user on a slow link received `IncompleteRead(163840
bytes read)`: exactly 160 KiB, ten clean 16 KiB chunks of a 311 KB file. `IncompleteRead` is
neither an `HTTPError` nor a `URLError`, so a handler written for HTTP status codes will not
catch it.

**Why truncation must be checked rather than trusted.** Feeding those exact 163,840 bytes to
the mmCIF parser returns **150 residues of an expected 330, with no error raised**. Truncated
PAE and cluster JSON fail loudly on parse; a truncated mmCIF does not. Completeness is therefore
verified against the prediction's residue count, per NFR-012.

## 2b. mmCIF `_struct_conf` secondary structure

AFDB model files carry DSSP-style secondary structure already, so the notebook needs no external
assignment tool and no extra download. The loop appears once per model file:

```
_struct_conf.conf_type_id          # HELX_RH_AL_P, HELX_RH_3T_P, HELX_RH_PI_P, HELX_LH_PP_P, STRN, TURN_TY1_P, BEND
_struct_conf.beg_label_seq_id      # inclusive
_struct_conf.end_label_seq_id      # inclusive
```

Collapsed to three codes on a row labelled `SSE` (secondary structure element): any `HELX*`
to `H`, `STRN` to `E`, everything else to `T`.
Residues named by no record are unassigned and render blank.

**Measured 2026-09-23.** Present on every model checked. `AF-O15552-F1`: 111 helix and 17
turn residues of 330, no strand, consistent with a seven-transmembrane-helix receptor.
`AF-Q13148-F1` carries `STRN` over 21% of its residues. Coverage is partial by design: the
records describe regular secondary structure only, so a long disordered tail appears as a run
of blanks rather than as coil.

## 2c. TM-align via `tmtools`

`tmtools.tm_align(coords_a, coords_b, seq_a, seq_b)` takes two `(n, 3)` CA arrays with their
one-letter sequences and returns `u` (3x3 rotation), `t` (translation), `rmsd`,
`tm_norm_chain1`, `tm_norm_chain2`, and `seqxA` / `seqyA` / `seqM`, the structure-based
alignment as gapped strings.

**Two contracts that are silent when broken, and both were verified by measurement:**

| Contract | Measured | Why it matters |
|---|---|---|
| `u @ a + t` maps **chain 1 onto chain 2** | Applying it to `P69905` against `P0C0U8` reproduces the reported RMSD, 1.34 A against 1.33 A | The notebook holds the query fixed, so it needs the **inverse**: `uᵀ`, `-uᵀt`. Verified orthogonal, determinant +1 |
| MolViewSpec `.transform(rotation=...)` reads 9 values **column-major** | `u.flatten(order="F")` | Row-major produces a plausible-looking but wrongly oriented model, with no error anywhere |

The sequences passed in are truncated to the CA count, because a model's `sequence` field may be
longer than the residues actually present in its coordinates.

**Caution on comparing implementations.** Mol\*'s JavaScript TM-align and `tmtools` agree on easy
pairs and diverge on hard ones (`O15552`: 2.94 A against 5.35 A). They optimise the same objective
by different heuristics. The notebook reports whichever it ran, and says which.

## 3. AFDB MSA files

```
GET {msaUrl}      # e.g. https://alphafold.ebi.ac.uk/files/msa/AF-Q9I1F6-F1-msa_v6.a3m
```

Format: a3m. **Currently returns HTTP 403** with a 134-byte HTML body, across multiple accessions and
both v4 and v6.

| Check | Rule |
|---|---|
| Availability | Detect by HTTP status **and** content type. Never parse the body; it is HTML, not JSON |
| Insertions | Lowercase characters in a3m are insertions and must be stripped before the alignment is rectangular |
| Failure | Print a named message stating the endpoint is under maintenance, and continue |

**This path is untested against a working service** and must be declared untested in the notebook.

## 4. UniProt taxonomy REST (REMOVED)

**No longer used.** The taxonomic composition figure was removed on 2026-09-23: the cluster API
returns free-text `speciesNames` with no taxId and no lineage, so the figure depended entirely on
this secondary lookup. Retained here only so a future reader knows what was measured and why it was
dropped, rather than rediscovering it.

What was measured while it was in scope:

| Constraint | Value |
|---|---|
| Query form | `(scientific:"X" OR scientific:"Y" ...) AND rank:genus` |
| **Batch size** | **20. Never more.** At 30, results crowd out the `size` window and hits are silently dropped, giving 36% coverage instead of 75% |
| Measured performance | 53 of 60 genera covering 75.4% of members in 1.82 s across 3 calls |
| Resolution | **Genus or higher only.** Species level is infeasible (3,274 unique names at ~0.12 s each) and unsafe, since `"Streptomyces sp."` fuzzy-matches to an arbitrary named species |
| Failure | Skip the figure with a named message. Never block the run |

Unresolved genera must be kept **visually distinct** from genuinely rare taxa in the figure.
Conflating them previously produced a misleading low-confidence row that was a name-resolution
artefact rather than biology.

## 5. User-supplied MSA file

| Check | Rule |
|---|---|
| Format | a3m in v1 (OQ-5) |
| Identity validation | The gap-stripped target sequence must equal the AFDB `sequence` for that accession |
| Naming convention | `AF-{ACC}-F1` in the header is a **hint that speeds lookup, never the validation**. A correctly named file with the wrong sequence must be rejected |
| Failure | Reject, naming the mismatch |

## 6. Colour schemes

`cleancolors.json` from the user's `msa-colorschemes` archive. Flat two-level dict:
`scheme -> residue letter -> "#RRGGBBAA"`. Strip the trailing alpha before use.

Twelve schemes, all flat per-residue lookups: `buried`, `cinema`, `clustal`, `clustal2`, `helix`,
`hydrophobicity`, `lesk`, `mae`, `strand`, `taylor`, `turn`, `zappo`. Note the scheme named `clustal`
is the simplified flat Jalview variant, **not** true conservation-dependent Clustal X colouring
(OQ-6).

The bundled `.woff2` colour fonts are a browser technique and are **not usable from matplotlib**.
Only the JSON is used.

**`cleancolors.json` is an optional enhancement, not a dependency.** It is not in the repo, a Colab
run cannot reach the user's local copy, and its licence is unconfirmed. The module therefore ships a
**built-in default residue colour scheme** that needs no external file, and loads `cleancolors.json`
only when it is present and readable. A missing file is never an error.

## 7. Rate and batch measurements

These are measured values with no automatic guard, so they need an explicit re-verification trigger
in the same way fixture JSON does.

| Measurement | Value | Re-verify when |
|---|---|---|
| Cluster endpoint spacing | 2.5 s | An AFDB release, the switch to the mirror, or any recurrence of spurious 5xx |
| Prediction endpoint spacing | none required | Same triggers |
| UniProt taxonomy batch size | 20 | A UniProt API change, or coverage dropping below roughly 70% |
| Install time for the three pip packages | roughly 6 s | Any dependency version bump, or a Colab image change |

## 7. Cross-source validation

| Check | Rule | On failure |
|---|---|---|
| Array length agreement | All seven cluster arrays equal length | Raise; the response is malformed |
| `clusterTotal` agreement | Equals array length | Warn and record; percentiles may be off by one (OQ-3) |
| Sort order | `averagePlddt` descending | Sort defensively rather than trusting it, but log that the assumption failed |
| Sequence length agreement | Prediction `sequence` length matches the cluster's `sequenceLength` for the same accession | Warn; do not block the alignment |
| Query presence | Whether the accession appears in each returned cluster | Record explicitly. Absence is a normal path, not an error |
