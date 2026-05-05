# PRD: Homodimer Confidence Metric Diagnostic Notebook

**Author:** Maxim (AFDB / PDBe, EMBL-EBI)
**Status:** Draft
**Last updated:** 2026-03-25

---

## 1. Use Case

For a given predicted homodimer from the AlphaFold Database, produce an interactive diagnostic notebook that decomposes and visualises each confidence metric (ipSAE, ipTM, pDockQ, pDockQ2, LIS) down to its contributing residues, PAE regions, and structural features, so that a user with no structural biology background can understand why these scores agree or disagree and assess the quality of the predicted complex.

---

## 2. Problem Statement

The AlphaFold Database released predicted homodimers filtered by a confidence band: `ipSAE >= 0.6`, `backbone clashes <= 10`, and `average pLDDT >= 70`. However, the suite of confidence metrics (ipTM, ipSAE, pDockQ, pDockQ2, LIS) frequently disagree with one another for a given complex. These scores each interrogate a different aspect of the prediction using different inputs, different residue subsets, and different regions of the PAE matrix, but currently there is no tool that exposes *why* they diverge. Users (and the AFDB team internally) see five opaque numbers with no mechanistic interpretation.

---

## 3. Target Audience

A user who is **unfamiliar with structural bioinformatics**. The notebook must:

- Explain what each metric measures in plain language before showing results
- Annotate every visualisation with what it shows and how to interpret it
- Provide a summary panel that synthesises findings into a clear diagnostic statement
- Use consistent colour schemes with legends throughout
- Avoid jargon without definition

---

## 4. Inputs

### 4.1 User-Provided Input

A single AFDB homodimer accession ID (e.g., `AF-0000000065889468`).

### 4.2 Data Retrieved from the AFDB API

The notebook fetches all required data programmatically. The API endpoint:

```
GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}
```

Returns a JSON array of entries. Each entry contains:

| Field | Description | Used For |
|-------|-------------|----------|
| `cifUrl` | URL to mmCIF structure file | Atomic coordinates, chain assignments, B-factor (pLDDT) |
| `pdbUrl` | URL to PDB format structure file | Alternative coordinate source |
| `paeDocUrl` | URL to PAE JSON file | Full predicted aligned error matrix |
| `plddtDocUrl` | URL to confidence JSON file | Per-residue pLDDT scores |
| `paeImageUrl` | URL to PAE image (PNG) | Optional quick-look reference |
| `sequence` | Amino acid sequence | Residue mapping and annotation |
| `uniprotAccession` | UniProt ID | Cross-referencing |

### 4.3 PAE JSON Structure

```json
[
  {
    "predicted_aligned_error": [[...], ...],  // NxN matrix (N = total residues across all chains)
    "max_predicted_aligned_error": 27.47,
    "chains": [
      {"name": "...", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 172},
      {"name": "...", "label_asym_id": "B", "sequenceStart": 1, "sequenceEnd": 172}
    ]
  }
]
```

The PAE matrix is symmetric in dimensions but **not** in values (PAE[i][j] != PAE[j][i]). For a homodimer of length L, the matrix is (2L x 2L). Rows 0..(L-1) are chain A residues; rows L..(2L-1) are chain B residues.

### 4.4 Confidence JSON Structure

```json
{
  "residueNumber": [1, 2, ..., 344],
  "confidenceScore": [69.38, 76.88, ...],
  "confidenceCategory": ["Low", "Confident", ...],
  "chains": [
    {"name": "...", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 172},
    {"name": "...", "label_asym_id": "B", "sequenceStart": 1, "sequenceEnd": 172}
  ]
}
```

### 4.5 mmCIF Structure

Standard mmCIF with `_atom_site` records. Two chains (A and B) with identical sequences. CA (alpha carbon) and CB (beta carbon, CA for glycine) coordinates are required for distance calculations.

---

## 5. Score Definitions and Formulas

All scores follow the exact implementations in the [DunbrackLab/IPSAE repository](https://github.com/DunbrackLab/IPSAE) (`ipsae.py`, version 4, January 2026). The notebook must recompute each score from scratch using the same formulas, not call the script as a subprocess.

### 5.1 Shared Functions

#### TM-score function (ptm_func)

```
ptm(x, d0) = 1 / (1 + (x / d0)^2)
```

Where `x` is the PAE value and `d0` is a length-dependent normalisation factor.

#### d0 calculation

```
d0(L) = max(1.0, 1.24 * (L - 15)^(1/3) - 1.8)   for L > 27
d0(L) = 1.0                                         for L <= 27
```

Where `L` is the number of residues used for normalisation (varies by score variant).

### 5.2 Interface Detection

Contacts are identified using **CB-CB distances** (CA for glycine), with a cutoff of **8.0 Å**. This follows the `pDockQ_cutoff = 8.0` in the IPSAE code, which uses CB coordinates for the distance matrix.

Note: The AFDB team's production pipeline (see attached `interface.py`) uses CA-CA distances at 8.0 Å via a GPU-accelerated `radius_graph` implementation. The notebook should use CB-CB (CA for GLY) to match the IPSAE scoring code, but should document this distinction.

### 5.3 ipTM (ipTM_d0chn)

**What it measures:** Global inter-chain confidence based on the full PAE matrix (no PAE cutoff applied). Normalised by chain pair length.

**Formula:**
- For each residue `i` in chain1, compute the mean of `ptm(PAE[i][j], d0_chn)` over all residues `j` in chain2
- `d0_chn = d0(len(chain1) + len(chain2))`
- The per-chain-pair asymmetric value is the **maximum** per-residue mean across all residues in chain1
- The symmetric "max" value is `max(ipTM(A->B), ipTM(B->A))`

**No PAE cutoff** is applied. All inter-chain PAE values contribute.

**Key insight for users:** ipTM captures overall confidence in the relative positioning of the two chains. A high ipTM means AlphaFold is confident about where chain B sits relative to chain A, globally.

### 5.4 ipSAE (ipSAE_d0res)

**What it measures:** Inter-chain confidence with a PAE cutoff, normalised by the number of residues in chain2 that have good PAE values for a given residue in chain1.

**AFDB parameters:** PAE cutoff ≤ 10, distance cutoff ≤ 15.

**Formula:**
- For each residue `i` in chain1, identify residues `j` in chain2 where `PAE[i][j] < pae_cutoff` (the "valid pairs")
- `n0_res[i]` = count of valid pairs for residue `i`
- `d0_res[i] = d0(n0_res[i])`
- Per-residue score: mean of `ptm(PAE[i][j], d0_res[i])` over valid pairs only
- The chain-pair asymmetric value is the **maximum** per-residue score across all residues in chain1
- The symmetric value is `max(ipSAE(A->B), ipSAE(B->A))`

**Key insight for users:** ipSAE is stricter than ipTM because it only considers residue pairs where the PAE is already below a quality threshold. It asks: "Among the residues AlphaFold is somewhat confident about, how good is that confidence?" The per-residue d0 normalisation means small, tight interfaces can score well even if the overall chain is large.

### 5.5 Additional ipSAE Variants

The IPSAE code computes three d0 variants. The notebook should compute and display all three:

| Variant | d0 based on | Use |
|---------|-------------|-----|
| `ipSAE_d0res` | Num residues in chain2 with PAE < cutoff for each residue in chain1 | **Primary AFDB classifier** |
| `ipSAE_d0chn` | Total residues in chain pair (len(chain1) + len(chain2)) | More conservative; penalises small interfaces |
| `ipSAE_d0dom` | Total residues in both chains that have any inter-chain PAE < cutoff | Intermediate normalisation |

### 5.6 pDockQ

**What it measures:** Post-hoc structural plausibility based on contact count and pLDDT. Does **not** use the PAE matrix.

**Formula:**
- Count inter-chain contacts: CB-CB distance ≤ 8.0 Å (CA for GLY)
- Identify all interface residues (residues in either chain involved in at least one contact)
- `x = mean_pLDDT_of_interface_residues * log10(number_of_contacts)`
- `pDockQ = 0.724 / (1 + exp(-0.052 * (x - 152.611))) + 0.018`

**pLDDT source:** CB pLDDT values (from the B-factor column of the structure file or from the confidence JSON).

**Key insight for users:** pDockQ asks "Does this interface look physically reasonable?" A high pDockQ means there are many contacts between well-resolved residues. It can diverge from PAE-based scores when the structure looks plausible but AlphaFold wasn't confident about the chain arrangement, or vice versa.

**Reference:** Bryant, Pozzati, and Elofsson (2022). Nature Communications. https://www.nature.com/articles/s41467-022-28865-w

### 5.7 pDockQ2

**What it measures:** Combines contact-based assessment with PAE confidence at the interface.

**Formula:**
- For each residue `i` in chain1, find residues `j` in chain2 within 8.0 Å (CB-CB)
- For each such contact pair, compute `ptm(PAE[i][j], 10.0)` (fixed d0 of 10.0)
- `mean_ptm = sum of ptm values / total number of contact pairs`
- `x = mean_pLDDT_of_interface_residues * mean_ptm`
- `pDockQ2 = 1.31 / (1 + exp(-0.075 * (x - 84.733))) + 0.005`

**Key insight for users:** pDockQ2 bridges the gap between structure-based and PAE-based scores. It considers both whether contacts exist (structural) and whether AlphaFold is confident about those contacts (PAE). When pDockQ is high but pDockQ2 is low, it means the contacts exist but AlphaFold has high uncertainty about them.

**Reference:** Zhu, Shenoy, Kundrotas, Elofsson (2023). Bioinformatics. https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714

### 5.8 LIS (Local Interaction Score)

**What it measures:** How "intertwined" the interface is, based on a linear transform of PAE values below a threshold of 12.

**Formula:**
- Select all PAE values for the inter-chain block (chain1 rows, chain2 columns)
- Filter to values where `PAE < 12`
- For each valid value: `score = (12 - PAE) / 12`
- `LIS = mean(scores)` over all valid values
- The symmetric LIS is the average of `LIS(A->B)` and `LIS(B->A)`

**Key insight for users:** LIS captures how densely the two chains interact at the PAE level. A high LIS with low ipTM can indicate that the chains are entangled in the prediction (many low-PAE inter-chain pairs) but AlphaFold wasn't confident about the overall arrangement. LIS is more of a "quantity of interaction" measure than a "quality of interaction" measure.

**Reference:** Kim, Hu, Comjean, Rodiger, Mohr, Perrimon (2024). bioRxiv. https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1

---

## 6. Notebook Structure

The notebook is divided into sequential sections. Each section must include a markdown explanation cell before any code, written for the target audience (Section 3).

### Section 1: Setup and Data Loading

**Purpose:** Install dependencies, fetch data from the AFDB API.

**Steps:**
1. Accept user input: AFDB accession ID
2. Call `GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}`
3. Download the mmCIF structure file (`cifUrl`)
4. Download the PAE JSON (`paeDocUrl`)
5. Download the confidence JSON (`plddtDocUrl`)
6. Parse the mmCIF to extract:
   - CA and CB coordinates per residue per chain
   - Residue names and numbers
   - Chain assignments (A, B)
7. Parse the PAE JSON into an NxN numpy matrix
8. Parse the confidence JSON into per-residue pLDDT arrays
9. Display basic metadata: protein name, UniProt ID, organism, sequence length, number of chains

**Dependencies:** `numpy`, `matplotlib`, `seaborn`, `requests`, `molviewspec` (`pip install molviewspec`). No GPU-dependent libraries. No BioPython requirement (parse mmCIF manually to avoid dependency weight).

### Section 2: Interface Detection

**Purpose:** Identify which residues form the interface between the two chains.

**Steps:**
1. Compute the CB-CB distance matrix (CA-CA for GLY) between all residue pairs
2. Identify inter-chain contacts using 8.0 Å cutoff
3. Mark interface residues in each chain
4. Report: number of interface residues per chain, number of contact pairs, list of contacting residue ranges

**Visualisations:**
- **Contact map** (matplotlib): 2D scatter or heatmap showing which residue pairs in chain A vs chain B are in contact. Color by distance (darker = closer).
- **Interface coverage bar**: A horizontal bar for each chain showing which residues are at the interface (colored) vs. not (grey). This gives an immediate visual of how much of each chain participates.

### Section 3: PAE Matrix Decomposition

**Purpose:** Show the full PAE matrix and highlight which regions each score uses.

**Steps:**
1. Display the full PAE matrix as a heatmap
2. Annotate the four quadrants: intra-chain A (top-left), inter-chain A->B (top-right), inter-chain B->A (bottom-left), intra-chain B (bottom-right)
3. Overlay masks showing:
   - **ipTM region:** The entire inter-chain quadrants (no PAE cutoff)
   - **ipSAE region:** Only cells in the inter-chain quadrants where `PAE < 10`
   - **LIS region:** Only cells in the inter-chain quadrants where `PAE < 12`
   - **pDockQ2 region:** Only cells in the inter-chain quadrants where the corresponding residue pair is also within 8.0 Å CB-CB distance

**Visualisations:**
- **Full PAE heatmap** (seaborn/matplotlib): Standard heatmap with chain boundaries marked. Color scale: 0 (dark blue, confident) to 31.75 (red, uncertain).
- **Side-by-side PAE quadrant panels** (2x2 grid): Each panel shows the same inter-chain PAE block (A->B quadrant) but with a different mask applied (ipTM, ipSAE, LIS, pDockQ2), highlighting which cells contribute to each score. Non-contributing cells should be greyed out or made transparent.

### Section 4: Score Computation and Decomposition

**Purpose:** Recompute all five metrics from scratch, showing intermediate values at each step. The user sees exactly how each number is derived.

**Steps:**
1. Compute pDockQ: show contact count, interface residues, mean pLDDT, the `x` variable, and the sigmoid output
2. Compute pDockQ2: show mean_ptm at the interface, the `x` variable, and the sigmoid output
3. Compute LIS: show count of valid PAE values, mean score
4. Compute ipTM: show d0_chn, per-residue scores, the maximum residue
5. Compute ipSAE (all three d0 variants): show d0 values, per-residue scores, the maximum residue for each variant

**Visualisations:**
- **Per-residue score profiles** (matplotlib line plots): For each chain, plot the per-residue contribution to ipTM, ipSAE_d0res, ipSAE_d0chn, and ipSAE_d0dom as a function of residue number. Overlay on the same axes with a secondary y-axis for pLDDT. This shows where the high-scoring residues are and whether they correlate with high pLDDT.
- **Score intermediate values table** (displayed in notebook): A table showing, for each score, the key intermediate values (contact count, d0, n0, mean_pLDDT, etc.) that feed into the final number.

### Section 5: Per-Residue pLDDT at the Interface

**Purpose:** Bridge the PAE-based scores and the contact-based scores by examining backbone quality at the interface.

**Steps:**
1. Extract pLDDT for all interface residues (both chains)
2. Compare pLDDT distribution of interface residues vs. non-interface residues
3. Identify low-pLDDT interface residues (< 70) that may be dragging down pDockQ

**Visualisations:**
- **pLDDT distribution comparison** (matplotlib): Overlaid histograms or violin plots showing pLDDT distribution for interface vs. non-interface residues.
- **Per-residue pLDDT profile with interface annotation** (matplotlib): Line plot of pLDDT along the sequence, with interface residues highlighted (e.g., colored background band or markers).

### Section 6: 3D Structure Visualisation

**Purpose:** Map diagnostic information onto the 3D structure using MolViewSpec for interactive viewing.

**Technology:** [MolViewSpec](https://github.com/molstar/mol-view-spec) Python library (`pip install molviewspec`). MolViewSpec generates JSON view descriptions that render in Mol*, the same viewer used by PDBe and RCSB. In a Colab environment, the notebook should generate MolViewSpec JSON (`.mvsj`) and provide either:
- An inline Mol* widget (if ipywidget integration is available)
- A URL to open the view in the Mol* web viewer (fallback)

**Views to generate:**

1. **Overview:** Cartoon representation, chains colored distinctly (chain A: teal, chain B: coral), interface residues highlighted in a contrasting colour
2. **pLDDT mapping:** Cartoon representation colored by per-residue pLDDT using the standard AlphaFold colour scheme (dark blue > 90, light blue 70-90, yellow 50-70, orange < 50)
3. **Interface diagnostic:** Interface residues only, shown as ball-and-stick or surface, colored by the per-residue ipSAE_d0res contribution (gradient from red [low] to green [high])
4. **Disagreement regions:** Residues where PAE-based scores (ipSAE) suggest high confidence but contact-based scores (pDockQ) suggest low confidence, or vice versa, highlighted in distinct colours

**MolViewSpec implementation notes:**
- The Python `molviewspec` library uses a builder API to construct a tree of nodes (download -> parse -> structure -> component -> representation -> color)
- Per-residue coloring is achieved through MolViewSpec annotation files or inline color nodes with residue selectors (`label_asym_id`, `beg_label_seq_id`, `end_label_seq_id`)
- The structure can be loaded directly from the AFDB `cifUrl` or `bcifUrl` (BinaryCIF preferred for performance)
- Reference: Midlik et al. (2025), NAR. https://doi.org/10.1093/nar/gkaf370
- Example notebook: https://mybinder.org/v2/gh/molstar/mol-view-spec/master?labpath=test-data/notebooks/01_kras_structure_visualization.ipynb

### Section 7: Diagnostic Summary

**Purpose:** Synthesise all findings into a human-readable interpretation.

**Steps:**
1. Display a summary table of all five scores with traffic-light colouring (green = high confidence, amber = moderate, red = low confidence)
2. Compute a pairwise agreement matrix: for each pair of scores, report whether they agree (both high or both low) or disagree
3. Generate a plain-text diagnostic statement that answers: "Why do these scores agree or disagree for this complex?"

**Diagnostic logic (heuristics):**
- If all scores are high: "This complex has consistently high confidence across all metrics. The interface is well-resolved (high pLDDT), structurally plausible (high pDockQ), and the PAE matrix shows strong inter-chain confidence."
- If ipSAE is high but pDockQ is low: "AlphaFold is confident about the relative chain positioning (PAE), but there are few physical contacts at the interface. This can occur when chains interact via a small, tight interface or when the predicted distance between chains is slightly too large for contacts to form."
- If pDockQ is high but ipSAE is low: "The interface has many contacts between well-resolved residues, but AlphaFold's PAE indicates uncertainty about the relative chain arrangement. This can occur when the local structure of each chain is well-predicted but the docking orientation is uncertain."
- If pDockQ2 diverges from pDockQ: "pDockQ and pDockQ2 disagree, which indicates that although contacts exist, the PAE confidence at those specific contact points is low (or high). pDockQ2 incorporates PAE at the interface, so it's the more informative of the two."
- If LIS is high but ipSAE is low: "Many inter-chain PAE values are below 12 (dense interaction), but when the PAE cutoff is tightened to 10 and the TM-score formula is applied, the confidence drops. This suggests a broad but diffuse interaction rather than a tight, well-defined interface."

These are starting heuristics. The notebook should present them as interpretive guidance, not definitive classification.

---

## 7. Visualisation Specifications

### 7.1 Colour Schemes

| Element | Colour scheme | Library |
|---------|--------------|---------|
| PAE heatmap | Sequential: dark blue (0) -> white (15) -> red (31.75) | matplotlib `RdBu_r` or custom |
| pLDDT | AlphaFold standard: dark blue (>90), light blue (70-90), yellow (50-70), orange (<50) | Custom discrete cmap |
| Chain A | `#009688` (teal) | Consistent across all plots |
| Chain B | `#FF7043` (coral) | Consistent across all plots |
| Interface residues | `#FFC107` (amber) when shown as a highlight | |
| Contact map | Sequential: light -> dark by distance | `viridis_r` |
| Per-residue score profiles | Distinct line colours per score | matplotlib default or tab10 |

### 7.2 Figure Sizing and Style

- All matplotlib figures: `figsize` proportional to content, minimum `(10, 8)` for heatmaps
- Font size: minimum 12pt for labels, 10pt for tick labels
- All axes labelled with units
- All heatmaps include a colour bar
- `seaborn` style: `whitegrid` or `white`
- DPI: 150 for inline display

### 7.3 3D Visualisation (MolViewSpec)

- Builder pattern using the `molviewspec` Python API
- Structure loaded from AFDB BinaryCIF URL for performance
- Multiple views generated as separate MolViewSpec states
- Each view saved as `.mvsj` JSON and displayed via Mol* viewer link or inline widget

---

## 8. Acceptance Criteria

### 8.1 Functional Requirements

1. **Single input:** User provides only an AFDB accession ID. Everything else is fetched automatically.
2. **Score accuracy:** All recomputed scores must match the output of `ipsae.py` (DunbrackLab) when run with the same PAE cutoff (10) and distance cutoff (15) on the same structure and PAE data. Tolerance: ±0.001 for floating point.
3. **All five scores computed:** ipSAE (all three d0 variants), ipTM_d0chn, pDockQ, pDockQ2, LIS.
4. **Per-residue decomposition:** For ipSAE and ipTM, the notebook must show per-residue contribution profiles, not just aggregate scores.
5. **PAE matrix visualisation:** Full heatmap with chain boundaries and per-score region overlays.
6. **3D visualisation:** At least two MolViewSpec views (chain overview + pLDDT mapping). Interface and diagnostic views are stretch goals.
7. **Diagnostic summary:** Plain-text interpretation generated from score values.

### 8.2 Non-Functional Requirements

1. **Runtime:** < 60 seconds on a free-tier Google Colab instance (no GPU).
2. **Dependencies:** Only `numpy`, `matplotlib`, `seaborn`, `requests`, `molviewspec`. No BioPython, no PyTorch, no torch-geometric.
3. **Self-contained:** All code in a single `.ipynb` file. No external scripts or data files required (data is fetched from AFDB API).
4. **Educational:** Every code cell preceded by a markdown explanation cell. Target audience is a newcomer (Section 3).

### 8.3 Verification

The engineering session should validate the notebook against at least three test cases:

| Accession | Expected behaviour |
|-----------|--------------------|
| `AF-0000000065889468` | High-confidence homodimer (E. coli FabA). All scores should be high and in agreement. |
| (TBD: a borderline case) | Moderate ipSAE near the 0.6 cutoff. Scores may disagree. |
| (TBD: a disagreement case) | Known case where ipSAE is high but pDockQ is low, or vice versa. |

---

## 9. Technical Reference

### 9.1 IPSAE Repository

- **URL:** https://github.com/DunbrackLab/IPSAE
- **Script:** `ipsae.py` (version 4, January 2026)
- **License:** MIT
- **Paper:** https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2

### 9.2 MolViewSpec

- **Python package:** `pip install molviewspec` (PyPI: https://pypi.org/project/molviewspec/)
- **GitHub:** https://github.com/molstar/mol-view-spec
- **Documentation:** https://molstar.org/mol-view-spec-docs/
- **NAR paper:** Midlik et al. (2025). https://doi.org/10.1093/nar/gkaf370
- **Example notebook:** https://mybinder.org/v2/gh/molstar/mol-view-spec/master?labpath=test-data/notebooks/01_kras_structure_visualization.ipynb

### 9.3 Score References

| Score | Reference |
|-------|-----------|
| ipSAE | Dunbrack Lab (2025). bioRxiv. https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2 |
| pDockQ | Bryant, Pozzati, Elofsson (2022). Nat Commun. https://www.nature.com/articles/s41467-022-28865-w |
| pDockQ2 | Zhu, Shenoy, Kundrotas, Elofsson (2023). Bioinformatics. https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714 |
| LIS | Kim, Hu, Comjean, Rodiger, Mohr, Perrimon (2024). bioRxiv. https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1 |
| ipTM / pTM | Yang and Skolnick (2004). Proteins. TM-score d0 formula origin |

### 9.4 AFDB API

- **Prediction endpoint:** `GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}`
- **Example homodimer accession:** `AF-0000000065889468`
- **Homodimer files served:** mmCIF, PDB, PAE JSON, confidence JSON, PAE image

### 9.5 Interface Detection Reference

The attached `interface.py` (Maciej Majewski, Apache 2.0) provides the AFDB production implementation using PyTorch/torch-geometric. The notebook does **not** use this implementation but should match its interface definition logic (8.0 Å CA-CA cutoff). The notebook uses CB-CB distances (CA for GLY) to match the IPSAE scoring code's contact definition.

---

## 10. Out of Scope

- Batch processing of multiple homodimers (this is a single-complex diagnostic tool)
- Integration with the InsightFold MCP layer or Nextflow pipelines
- Comparison with experimental structures from the PDB
- Computation of backbone clashes (one of the AFDB filter criteria, but outside the scope of confidence metric decomposition)
- Support for heteromeric complexes (homodimers only for this use case)
- Any modifications to the IPSAE scoring formulas

---

## 11. Future Extensions (Not In Scope for v1)

- **Batch pre-computation:** Compute a "disagreement score" (rank variance across metrics) in BigQuery for all homodimers, surface outliers for drill-down with this notebook
- **InsightFold flywheel integration:** Users submit diagnostic results back to enrich the pre-computed layer
- **Heteromer support:** Extend to non-symmetric complexes where asymmetric ipSAE values (A->B vs B->A) differ substantially
- **Comparative view:** Side-by-side comparison of two homodimers from the same protein family
