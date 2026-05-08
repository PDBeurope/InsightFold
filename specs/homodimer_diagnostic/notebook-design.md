# Homodimer Diagnostic Notebook Design

## Notebook Target

Path: `notebooks/homodimer_diagnostic.ipynb`

Runtime target: Google Colab free tier, CPU-only, less than 60 seconds for the happy-path fixture.

Dependency policy:

- Allowed: `numpy`, `matplotlib`, `seaborn`, `requests`, `json`, standard library, `molviewspec`.
- Disallowed for v1: BioPython, PyTorch, torch-geometric, GPU-only dependencies, subprocess calls to external scoring scripts, local-only data files.

## Notebook Sections

| Section | Purpose | Key Outputs | Requirements |
|---|---|---|---|
| 1. Purpose and RUO framing | Explain the notebook, audience, and research-use limits | Intro markdown | REQ-001, REQ-021 |
| 2. Setup and accession input | Imports, optional `molviewspec` install/import, accession parameter | `ACCESSION_ID`, imports | REQ-001, REQ-020 |
| 3. Fetch AFDB data | Retrieve metadata and download mmCIF/PAE/confidence JSON | `metadata`, `cif_text`, `pae_json`, `confidence_json` | REQ-002, REQ-003 |
| 4. Validate data contracts | Check required fields and basic schema | validation summary | REQ-002, REQ-003, REQ-005 |
| 5. Parse structure and confidence | Parse mmCIF atoms, chain residues, pLDDT | `chains`, `residues`, coordinates, pLDDT arrays | REQ-004, REQ-006, REQ-007 |
| 6. Interface detection | Compute CB/CA distances and contacts | contact pairs, interface residues, contact count | REQ-008 |
| 7. PAE decomposition | Display full PAE matrix and score masks | heatmaps, masks | REQ-016 |
| 8. Metric computation | Compute ipTM, ipSAE variants, pDockQ, pDockQ2, LIS | final score table, intermediate tables | REQ-009 to REQ-014 |
| 9. Per-residue and interface profiles | Show score drivers per residue and interface/non-interface pLDDT | profiles, distribution plots | REQ-015, REQ-017 |
| 10. MolViewSpec views | Create chain overview and pLDDT mapping | `.mvsj` states or fallback message | REQ-018 |
| 11. Diagnostic summary | Explain why scores agree or disagree | plain-language summary panel | REQ-019, REQ-021 |
| 12. Validation snapshot | Show fixture checks, runtime notes, versions | validation summary | REQ-020 |

Every code section must be preceded by explanatory markdown written for non-specialists.

## Data Flow

```text
ACCESSION_ID
  -> AFDB prediction endpoint
  -> metadata record
  -> cifUrl, paeDocUrl, plddtDocUrl
  -> mmCIF parser -> residue table, CA/CB coordinates, B-factor pLDDT fallback
  -> PAE parser -> PAE matrix, chain ranges
  -> confidence parser -> pLDDT arrays
  -> homodimer validation
  -> interface detection
  -> metric computations
  -> tables, figures, MolViewSpec states, diagnostic summary
```

## Function Candidates

Reusable computations should be implemented as notebook-local functions first and may later move to `src/insightfold/` if reused.

| Function | Responsibility |
|---|---|
| `fetch_afdb_prediction(accession_id)` | Return validated metadata record and source URLs |
| `download_text(url)` / `download_json(url)` | Retrieve files with clear network errors |
| `parse_mmcif_atoms(cif_text)` | Extract atom records needed for coordinates and B-factors |
| `build_residue_table(atom_records)` | Produce one residue row per chain/residue with CA/CB coordinates |
| `parse_pae_payload(pae_json)` | Return PAE matrix and chain metadata |
| `parse_confidence_payload(confidence_json)` | Return per-residue pLDDT and confidence categories |
| `validate_homodimer(chains, pae_matrix)` | Confirm two-chain homodimer compatibility |
| `compute_interchain_contacts(residue_table, cutoff=8.0)` | Return contact pairs and interface residue masks |
| `d0(length)` | IPSAE/TM-score normalization |
| `ptm(pae, d0_value)` | PAE-to-confidence transform |
| `compute_iptm(...)` | ipTM asymmetric and symmetric scores |
| `compute_ipsae_variants(...)` | ipSAE d0res/d0chn/d0dom scores and intermediates |
| `compute_pdockq(...)` | pDockQ score and intermediates |
| `compute_pdockq2(...)` | pDockQ2 score and intermediates |
| `compute_lis(...)` | LIS asymmetric and symmetric scores |
| `build_diagnostic_summary(score_table, intermediates)` | Plain-language score interpretation |

## Visualization Plan

- Contact map: chain A residue index versus chain B residue index, color-coded by distance.
- Interface coverage: one horizontal track per chain, interface residues highlighted.
- Full PAE heatmap: quadrant labels and chain boundary lines.
- PAE masks: score-specific overlays for ipTM, ipSAE, pDockQ2 contacts, and LIS.
- Per-residue profiles: ipTM/ipSAE contribution curves for A->B and B->A.
- pLDDT comparison: interface versus non-interface distribution.
- MolViewSpec: chain overview and pLDDT mapping view.

## Error Handling

The notebook should stop early with a readable message for:

- AFDB response is empty.
- Required URL field is missing.
- PAE matrix shape does not match chain lengths.
- Structure is not a two-chain homodimer.
- Required coordinates are absent for too many residues.
- No inter-chain contacts are found.

Optional MolViewSpec failures should not block core metric computation.

## Assumptions

- Confidence JSON is preferred for pLDDT; mmCIF B-factors are fallback.
- CB-CB contact definition with CA fallback for glycine is used to match IPSAE-style scoring.
- AFDB production CA-CA interface definition is documented but not used for metric recomputation.
- At least one happy-path accession can be used for smoke validation; full validation needs three curated accessions.

