# Fixture Manifest

**Status: pinned, measured and machine-verified, 2026-09-21.** Every numeric value below is
generated from `_expected.json` / `_alignments.json`, not transcribed by hand, and
`src/insightfold/cluster_quality_utils.py` reproduces all 61 of them exactly.

Files live in `specs/cluster_quality_diagnostic/fixtures/`, gzipped, 3.15 MB total.

> An earlier revision of this file carried hand-transcribed values that did not match the pinned
> data. They were caught by running the module against the fixtures. Do not hand-edit the numeric
> rows; regenerate them from `_expected.json` and `_alignments.json`.

| Field | Value |
|---|---|
| Captured | 2026-09-21 |
| Source | `https://alphafold.ebi.ac.uk/api/` |
| Cluster endpoint spacing used | 2.6 s |
| Prediction endpoint spacing used | none |
| Layout | `fixtures/cluster/{acc}_{mm\|fs}.json.gz`, `fixtures/prediction/{acc}.json.gz`, `fixtures/synthetic/`, plus `_expected.json`, `_prediction.json`, `_alignments.json`. All four directories and these three files are committed; the blanket `*.json` ignore rule that silently excluded them was corrected 2026-10-02. A capture receipt was also written at pin time and has since been deleted: it recorded only the request URL, status and byte sizes, none of which anything reads, and the keys already encode accession and clustering |
| Refresh trigger | An AFDB release, or the switch to the rate-limited mirror. Re-capture, diff `_expected.json`, investigate any change rather than accepting it |

Marker key: `exact` values must match; `presence` means the field or behaviour must exist;
`manual-review` means a human judges it.

---

## Happy path

Role: `happy-path`.

### FX-01 `P69905`, haemoglobin alpha. The notebook default.

| Property | Value | Mark |
|---|---|---|
| Source | `cluster/P69905_mm.json.gz`, `cluster/P69905_fs.json.gz` | |
| MMseqs2 | n 5,223, mean 96.3215, median 97.12, sd 3.467, IQR **0.76**, len CV 0.0226 | exact |
| Bands | >90 0.963431, 70-90 0.032548, 50-70 0.003446, <50 0.000574 | exact |
| Query | index **17**, percentile 99.5788 | exact |
| Best / worst | `AF-P01956-F1` 98.12 / `AF-A0A4X1SZ27-F1` 43.94 | exact |
| Foldseek | n **53**, mean 89.4721, sd **11.8155**, IQR 11.12 | exact |
| Alignments | best 90.07% over 141 cols, cov 99.3%; worst 32.85% over 137 cols, cov 85.21% | exact |
| Sorted descending | true | exact |

**Why it is the default.** Coherent family, best member plainly the same protein, fast. Its IQR of
**0.76 pLDDT points** is also the sharpest demonstration that a percentile within the body is
meaningless. Its Foldseek sd is **3.4 times** its MMseqs2 sd, which is the clearest single-fixture
evidence for REQ-006.

---

## Edge cases

Role: `edge-case`.

### FX-02 `Q13148`, TDP-43. Low-quality family, disorder confound, domain-level match.

| Property | Value | Mark |
|---|---|---|
| Source | `cluster/Q13148_mm.json.gz` | |
| MMseqs2 | n 33,211, mean 62.4395, median 61.88, sd 6.3655, IQR **9.22**, len CV 0.1255 | exact |
| Bands | >90 0.0, 70-90 0.129415, 50-70 0.85297, <50 0.017615 | exact |
| Query | index 11,493, percentile 65.3007 | exact |
| Best / worst | `AF-A0A0L8VJ73-F1` 86.81 (Prp24p) / `AF-A0A438IRX9-F1` 30.75 | exact |
| Alignment vs best | **26.73% over 101 cols, cov 21.26%** of a 414-residue query | exact |
| Per-residue profile | Ordered N-terminal region plus a low-confidence C-terminal tail, not a flat low trace | manual-review |

**Three jobs.** It is the low-mean family (REQ-003 must not call it "low quality"), the disorder case
(REQ-E8, REQ-010a), and the domain-level-match case. The best member aligns over **residues 80-171 of
414**, which sits on RRM1 (106-176), so this is a genuine single-domain match. The **aligned range is
the stable signal**; coverage moves between 21% and 55% across gap settings, so it is evidence for a
reader, not a threshold.

### FX-03 `Q9I1F6`, gluconate operon repressor. Distant but legitimate.

| Property | Value | Mark |
|---|---|---|
| Source | `cluster/Q9I1F6_mm.json.gz`, `cluster/Q9I1F6_fs.json.gz` | |
| MMseqs2 | n 4,628, mean 91.0368, median 91.25, sd 2.4982, IQR **2.0**, len CV 0.0631 | exact |
| Query | index 1,650, percentile 63.9585 | exact |
| Best / worst | `AF-A0A6J4IYQ6-F1` 94.69 / `AF-A0A415C607-F1` 49.81 | exact |
| Alignment vs best | **25.09% over 279 cols, cov 74.64%** | exact |
| Foldseek | n 913, mean 86.8859, sd 6.1435 | exact |

**The false-positive control.** Identity (25.09%) is *lower* than `Q13148`'s (26.73%) while the
members are genuine LacI-family relatives, which is direct proof that identity alone cannot separate
the cases. Its best member aligns over residues 18-277 of 343, a full-length match rather than a
single domain. Any rule that flags this cluster is wrong.

### FX-04 `P0A6F5`, GroEL. Size stress.

| Property | Value | Mark |
|---|---|---|
| Source | `cluster/P0A6F5_mm.json.gz` (8.82 MB raw, 1.49 MB gzipped) | |
| MMseqs2 | n 93,793, mean 88.4857, median 89.31, sd 2.5637, IQR **3.31**, len CV 0.0376 | exact |
| Query | index 11,255, percentile 87.5694 | exact |
| Best / worst | `AF-P51349-F1` 93.44 / `AF-A0AAE2BTP8-F1` 53.69 | exact |
| Runtime | Completes in about 1 s; no cap, pagination or sampling | range |

### FX-05 `P0AA25`, thioredoxin. Legitimate distant homologues, second control.

| Property | Value | Mark |
|---|---|---|
| Source | `cluster/P0AA25_mm.json.gz` | |
| MMseqs2 | n 62,441, mean 91.6856, median 92.81, sd 4.885, IQR **5.07**, len CV 0.0678 | exact |
| Query | index 13,792, percentile 77.2697 | exact |
| Alignment vs best | 33.68% over 95 cols, cov 84.4% | exact |

### FX-06 `C1C553`, sparse annotation, low pLDDT, coherent cluster.

| Property | Value | Mark |
|---|---|---|
| Source | `cluster/C1C553_mm.json.gz` | |
| MMseqs2 | n 110, mean 51.7245, median 51.515, sd 1.6678, IQR **2.5825**, len CV 0.0756 | exact |
| Query | index 8, percentile 91.8182 | exact |
| Alignment vs best | **96.15% over 312 cols, cov 99.68%** | exact |
| Low-n caveat | Does **not** fire; 110 is above the 30-member threshold. This fixture covers sparse annotation, not low n | presence |

**Teaching value.** Mean pLDDT 51.7 with 96% identity to its best member: a **coherent** cluster that
is simply hard to predict. The counterexample to "low pLDDT means a messy cluster".

### FX-07 `P04637-2`, isoform.

| Property | Value | Mark |
|---|---|---|
| Source | `cluster/P04637-2_mm.json.gz` | |
| MMseqs2 | n 1,024, mean 70.4885, median 71.0, sd 3.7344, IQR **3.56**, len CV 0.0795 | exact |
| Query index | **None.** `AF-P04637-2-F1` is **not** in the returned cluster | exact |
| Behaviour | Resolves to a cluster, but the isoform is absent from it. Must state which sequence is aligned and which structure rendered | presence |

**Unexpected result.** The isoform case also exercises query-absence, so REQ-001a and REQ-006 are
tested together on one fixture.

---

## Negative fixtures

Role: `negative`.

### FX-08 `Q8WZ42`, titin, fragmented. HTTP 404.

Source `cluster/Q8WZ42_mm.json.gz`, 81 bytes. Body: `{"detail":"No data found for the given protein
accession and cluster step flag."}`. Must refuse, naming fragmentation, and name both likely causes
per REQ-017. `exact` on the message being quoted, `presence` on the explanation.

### FX-09 `P00533`, EGFR Foldseek singleton.

Source `cluster/P00533_fs.json.gz`. Foldseek `clusterTotal` **1**, single member
`AF-A0A8J6H303-F1` at 76.62. **The partial-degradation case**: the MMseqs2 cluster is fine, so the
sequence analysis must complete and only the Foldseek section skip, per REQ-014a.

### FX-10 `P12345`, reproducible HTTP 500.

Source `cluster/P12345_mm.json.gz`, empty body. **OQ-4 is now RESOLVED**: this 500 reproduced at
**2.6 s spacing**, well inside the safe rate, so it is a genuine per-accession backend fault and not
rate limiting. Must retry with backoff, then report without presenting it as a data error.

---

## Synthetic fixtures

### FX-11 `SYN-COLLAPSE`, the collapse rule.

`synthetic/SYN-COLLAPSE-best.json.gz` and `-worst.json.gz`: the pinned `P69905` MMseqs2 response
truncated to its first 40 members. Treat the query as `AF-P01956-F1` (index 0) for the best case and
`AF-P08849-F1` (index −1) for the worst. Expected in both: **exactly 2 view targets and 1 alignment**,
not 3 and 2. Derivation is scripted and regenerates from the parent pin.

No real fixture has the query at an extreme, which is why this is synthetic.

### FX-12 `SYN-NOSEQ`, missing sequence.

`synthetic/SYN-NOSEQ.json.gz`: the pinned `A0A4X1SZ27` prediction response with `sequence` and
`uniprotSequence` removed.

**Important honesty note.** All 20 of `Q9I1F6`'s top members were checked: **14 are deleted from
UniProtKB, yet AFDB serves a sequence for every one of them.** No real accession is known to trigger
this path. The defensive code is kept because it costs nothing, but the manifest records that it has
**no known real-world trigger** rather than implying one exists.

### FX-13 `SYN-A3M`, MSA parsing and coverage plot.

`synthetic/SYN-A3M.a3m`: four sequences, the `P69905` query first under header `AF-P69905-F1`, with
a lowercase insertion run in one member and gap runs in others. Exercises REQ-V9's coverage plot
offline, since `msaUrl` returns 403 everywhere and the plotting code would otherwise ship unexecuted.

`synthetic/SYN-A3M-wrongseq.a3m`: identical header, reversed sequence. **Must be rejected** by
REQ-013's gap-stripped sequence match.

---

## What The Fixture Set Proves

| Claim | Fixture | Evidence |
|---|---|---|
| A percentile within the body is meaningless | FX-01 | IQR 0.76 points across 5,223 members |
| Foldseek clusters are wider | FX-01, FX-03 | sd ratios 0.29 and 0.41 |
| Identity alone cannot separate the cases | FX-02 vs FX-03 | 26.73% (divergent) is **above** 25.09% (legitimate) |
| The aligned range distinguishes them | FX-02 vs FX-03 | residues 80-171 of 414 (one domain) against 18-277 of 343 (full length). Coverage agrees at default gap settings but inverts at -20/-1, so it is not a rule |
| Low mean is not a messy cluster | FX-06 | mean 51.7, best member at 96% identity and 99.7% coverage |
| Arrays are sorted descending | all 10 | verified on every pinned response |
| Size is not a problem | FX-04 | 93,793 members in about 1 s |

### FX-14 `SYN-LOWN`, the low-n caveat.

`synthetic/SYN-LOWN.json.gz`: the pinned `C1C553` response truncated to 12 members, query still
present. Below the 30-member threshold, so the REQ-004 low-n caveat must fire on the percentile and
the IQR. Built because no real cluster under 30 members was located.

### Additional pins added to close the offline-coverage gap

Both clusterings are now pinned for every non-refusal fixture, and every view target has a
prediction response, so the whole validation plan runs offline.

| Pin | Result |
|---|---|
| `Q13148_fs` | `clusterTotal` 8,790 |
| `P0AA25_fs` | 5,902 |
| `P0A6F5_fs` | 1,567 |
| `P04637-2_fs` | 257 |
| `C1C553_fs` | 3 |
| `P00533_mm` | **380**, which is what makes FX-09 a partial-degradation case rather than a refusal: the sequence cluster is fine and only the structure cluster is unusable |
| prediction `Q8WZ42` | **HTTP 404.** Titin is absent from the prediction endpoint too, so this is a real 404 fixture for both endpoints |
| prediction `P12345` | **HTTP 200, 430 aa.** The protein exists in AFDB; only its MMseqs2 cluster call 500s, confirming a cluster-endpoint-specific fault |
| prediction `P51349`, `A0AAE2BTP8`, `P04637`, `P04637-2`, `P00533`, `A0A8J6H303`, `P0A6F5` | all 200 with sequences |

## Gaps And Blockers

| Item | Status |
|---|---|
| Real `msaUrl` a3m | **Blocked.** 403 everywhere. FX-13 covers the code path; the live path stays declared untested |
| Missing-sequence trigger | **No known real case.** FX-12 is synthetic by necessity, and this is recorded rather than hidden |
| Low-n caveat fixture (n below 30) | **Closed** by FX-14 `SYN-LOWN`, n = 12, verified to fire the caveat |
| Colab verification | **Open.** Requires an actual Colab run (T074) |
