# Homodimer Diagnostic Notebook Requirements

Source PRD: `prd/homodimer_diagnostic_notebook_prd.md`

## Goal

Build a self-contained InsightFold notebook that accepts one AFDB homodimer accession, retrieves public AFDB data, recomputes homodimer confidence metrics, visualizes the residue/PAE/interface evidence behind those metrics, and explains score agreement or disagreement for users without structural bioinformatics expertise.

## Users

- Non-specialist AFDB/PDBe users who need plain-language interpretation.
- AFDB/PDBe bioinformaticians who need reproducible score calculations and intermediate values.
- Platform engineers who need lightweight, deterministic notebook behavior.
- Research collaborators who need interface and confidence visualizations.

## Functional Requirements

| ID | Requirement | Acceptance Criteria |
|---|---|---|
| REQ-001 | Accept one AFDB homodimer accession as input | Notebook exposes one clearly labelled accession variable or input cell; invalid or unsupported inputs produce a clear message |
| REQ-002 | Retrieve AFDB metadata | Calls `GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}`; validates non-empty JSON response; displays accession, UniProt accession, sequence length, and available file URLs |
| REQ-003 | Download required data files | Downloads `cifUrl`, `paeDocUrl`, and `plddtDocUrl` where available; reports missing required fields before downstream computation |
| REQ-004 | Parse mmCIF without BioPython | Extracts chain IDs, residue numbers, residue names, CA coordinates, CB coordinates, and B-factor pLDDT fallback using lightweight parsing |
| REQ-005 | Parse PAE JSON | Loads `predicted_aligned_error`, `max_predicted_aligned_error`, and chain metadata; validates matrix shape against chain lengths |
| REQ-006 | Parse confidence JSON or fallback pLDDT | Uses `confidenceScore` from confidence JSON when available; falls back to mmCIF B-factors with a visible provenance note |
| REQ-007 | Validate homodimer assumption | Confirms exactly two chains with identical sequence length or records a clear unsupported-input error |
| REQ-008 | Detect inter-chain interface | Computes CB-CB distances with CA fallback for glycine; uses 8.0 Angstrom cutoff; reports contact count and interface residues |
| REQ-009 | Compute ipTM_d0chn | Computes asymmetric A->B and B->A values and symmetric max value using all inter-chain PAE values |
| REQ-010 | Compute ipSAE variants | Computes `ipSAE_d0res`, `ipSAE_d0chn`, and `ipSAE_d0dom` using PAE cutoff 10; reports intermediate d0 and residue-subset behavior |
| REQ-011 | Compute pDockQ | Computes pDockQ from contact count and interface mean pLDDT using the published logistic formula |
| REQ-012 | Compute pDockQ2 | Computes pDockQ2 from interface pLDDT and PAE-transformed contact confidence |
| REQ-013 | Compute LIS | Computes asymmetric A->B/B->A LIS over inter-chain PAE values below 12 and reports symmetric average |
| REQ-014 | Show intermediate values | Displays calculation tables for contact count, interface residues, d0 values, PAE subset sizes, pLDDT summaries, and final scores |
| REQ-015 | Visualize interface evidence | Produces contact map and interface coverage visualization with labelled axes, legends, and interpretation text |
| REQ-016 | Visualize PAE evidence | Produces full PAE heatmap with quadrant labels and score-specific masks for ipTM, ipSAE, pDockQ2, and LIS |
| REQ-017 | Visualize per-residue evidence | Produces per-residue contribution profiles for ipTM/ipSAE and pLDDT interface versus non-interface comparisons |
| REQ-018 | Generate MolViewSpec views | Produces at least chain overview and pLDDT mapping MolViewSpec states or a documented fallback when unavailable |
| REQ-019 | Explain score agreement/disagreement | Generates a plain-language diagnostic summary linking score patterns to PAE, contacts, pLDDT, and interface evidence |
| REQ-020 | Stay lightweight | Runs in less than 60 seconds on free-tier Colab for the happy-path fixture; uses no BioPython, PyTorch, torch-geometric, GPU-only dependencies, external scripts, or local data assumptions |
| REQ-021 | Preserve RUO framing | States that outputs are interpretive research-use guidance, not clinical or biological proof |

## Non-Goals

- Batch processing.
- Heteromer support.
- Experimental PDB structure comparison.
- Backbone clash calculation.
- Production AFDB/PDBe pipeline integration.
- Formula modification.
- Clinical decision support.

## Open Questions

| ID | Question | Status |
|---|---|---|
| OQ-001 | Which AFDB accessions should represent borderline and disagreement validation cases? | blocking for full validation |
| OQ-002 | Should pLDDT prefer confidence JSON or mmCIF B-factors when both are present? | assumption: prefer confidence JSON, fallback to mmCIF |
| OQ-003 | What score bands define high, moderate, and low for each metric? | blocking for final user-facing interpretation thresholds |
| OQ-004 | Should MolViewSpec interface/disagreement views be required in v1? | assumption: chain overview and pLDDT mapping are required; richer views are stretch |
| OQ-005 | Should generated `.mvsj` files be saved or displayed only inline? | assumption: save local session files and display links where possible |
| OQ-006 | How should non-homodimer inputs be handled? | assumption: fail early with explanation rather than attempting partial analysis |

