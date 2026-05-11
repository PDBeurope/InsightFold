# Product Requirements Document (PRD)
## Life Sciences / Bioinformatics / Research Platforms

---

# 1. Document Overview

| Field | Details |
|---|---|
| Product Name | Homodimer Confidence Metric Diagnostic Notebook |
| Version | v1 Draft |
| Date | 2026-03-25 |
| Author(s) | Maxim, AFDB / PDBe, EMBL-EBI |
| Stakeholders | AFDB team, PDBe users, structural bioinformatics researchers, bioinformaticians, platform engineers |
| Status | Draft |
| Related Documents | specs/homodimer_diagnostic_notebook_spec.md; agent-skills/specs/homodimer_diagnostic_notebook_spec.md |
| Regulatory Scope | RUO |

---

# 2. Executive Summary

## Problem Statement
The AlphaFold Database provides predicted homodimers filtered by global confidence criteria, but the major confidence metrics used to assess a complex can disagree. Users currently see opaque values for ipSAE, ipTM, pDockQ, pDockQ2, and LIS without a clear explanation of which residues, PAE regions, contacts, or structural features drive those scores.

---

## Product Vision
Create a self-contained interactive diagnostic notebook that recomputes and decomposes homodimer confidence metrics from first principles, then explains score agreement or disagreement in plain language for users without structural biology expertise.

---

## Success Criteria

- Recompute ipSAE, ipTM, pDockQ, pDockQ2, and LIS from AFDB data.
- Match DunbrackLab IPSAE reference calculations within +/- 0.001 where equivalent inputs are used.
- Show per-residue and PAE-region contributions rather than only aggregate scores.
- Provide annotated 2D and 3D visualisations that explain how to interpret each result.
- Generate a clear diagnostic summary describing why metrics agree or diverge.
- Run in less than 60 seconds on free-tier Google Colab without GPU dependencies.

---

# 3. Background & Context

## Scientific Context

- Domain: AlphaFold Database predicted homodimer confidence assessment.
- Biological object: two-chain homodimeric protein complexes with identical sequences.
- Core data: mmCIF coordinates, PAE matrix, per-residue pLDDT confidence, chain metadata.
- Score families: PAE-based metrics, contact-based metrics, and hybrid interface metrics.
- Current challenge: confidence scores use different residue subsets and formulas, causing disagreement that is difficult for non-specialists to interpret.

---

## Business / Institutional Context

- Supports AFDB and PDBe users evaluating predicted homodimer quality.
- Helps internal AFDB review of complexes passing filters such as `ipSAE >= 0.6`, backbone clashes <= 10, and average pLDDT >= 70.
- Improves transparency and educational value around confidence metrics.
- Provides a reusable notebook pattern for future InsightFold diagnostic workflows.

---

## Existing Solutions

| Solution | Strengths | Weaknesses |
|---|---|---|
| AFDB prediction page | Provides structure files, PAE data, confidence data, and quick visual reference | Does not explain why multiple interface confidence scores agree or disagree |
| DunbrackLab IPSAE script | Reference implementation for ipSAE-related scoring | Script-oriented; not designed as a teaching or diagnostic notebook for non-specialists |
| AFDB production interface pipeline | Production-scale interface detection implementation | Uses GPU-oriented dependencies and a CA-CA contact definition; not suitable for lightweight educational notebooks |
| MolViewSpec / Mol* | Rich interactive molecular visualisation | Requires notebook integration and generated view specifications |

---

# 4. Goals & Objectives

## Product Goals

| Goal | Metric |
|---|---|
| Provide single-accession diagnostic workflow | User supplies only one AFDB homodimer accession ID |
| Recompute all required confidence metrics | ipSAE variants, ipTM_d0chn, pDockQ, pDockQ2, and LIS generated in notebook |
| Explain score drivers | Per-residue profiles, PAE masks, contact maps, and intermediate value tables displayed |
| Support non-specialist interpretation | Every code section preceded by explanatory markdown |
| Keep execution lightweight | Runtime under 60 seconds in free-tier Colab; no GPU-only dependencies |

---

## Non-Goals

Clearly define what this product will **not** address.

- Batch processing of multiple homodimers.
- Heteromeric complex support.
- Comparison against experimental PDB structures.
- Backbone clash computation.
- Nextflow, MCP, or InsightFold production pipeline integration.
- Modification of the IPSAE scoring formulas.

---

# 5. Target Users

## Primary Users

| User Type | Goals | Pain Points |
|---|---|---|
| Bioinformaticians | Diagnose confidence metric agreement and disagreement for an AFDB homodimer | Need reproducible metric calculations and residue-level explanations |
| Wet-lab Scientists | Understand whether a predicted homodimer is credible enough to inspect further | May lack structural bioinformatics background |
| Clinical Researchers | Interpret confidence signals without treating the notebook as clinical evidence | Need plain-language limitations and RUO framing |
| Data Stewards | Review data provenance, formats, and confidence metadata | Need automatic retrieval and clear mapping of AFDB API fields |
| Platform Engineers | Maintain a lightweight notebook that avoids heavy dependencies | Need deterministic formulas and simple data flow |

---

## Secondary Users

- AFDB internal reviewers.
- PDBe support and documentation teams.
- Consortium partners.
- Research collaborators.
- Public users inspecting AFDB homodimers.

---

# 6. User Personas

## Persona 1 - Non-specialist Research User

### Background
- Works with protein function or experimental follow-up.
- Has limited structural bioinformatics training.

### Needs
- Plain-language explanation of each metric.
- Visual cues showing which residues and PAE regions drive confidence.
- A final diagnostic statement that explains agreement or disagreement.

### Frustrations
- Multiple scores can be high, low, or mixed without an obvious reason.
- PAE matrices and interface metrics are difficult to interpret without guidance.

---

## Persona 2 - AFDB / PDBe Bioinformatician

### Background
- Understands AFDB data products and confidence metrics.
- Needs to inspect homodimer predictions passing production filters.

### Needs
- Recomputed scores that match reference formulas.
- Intermediate values for debugging and validation.
- Visual comparisons between PAE-based and contact-based evidence.

### Frustrations
- Existing scripts produce scores but not an explanatory diagnostic view.
- Different interface definitions can create confusing differences between pipelines.

---

# 7. Product Scope

## In Scope

- Single AFDB homodimer accession input.
- AFDB API retrieval of prediction metadata, mmCIF, PAE JSON, and confidence JSON.
- Manual mmCIF parsing without BioPython.
- CB-CB interface detection with CA fallback for glycine.
- Computation of ipSAE, ipTM, pDockQ, pDockQ2, and LIS.
- Per-residue metric decomposition.
- PAE heatmaps and metric-region masks.
- pLDDT and interface quality visualisation.
- MolViewSpec views for chain overview and pLDDT mapping.
- Plain-text diagnostic summary.

---

## Out of Scope

- Clinical diagnosis or clinical decision support.
- EHR storage or clinical workflow integration.
- Wet-lab automation.
- Billing systems.
- Large-scale homodimer batch analysis.
- Heteromeric complexes.
- Experimental structure comparison.

---

# 8. Use Cases

| ID | Use Case | Priority |
|---|---|---|
| UC-01 | Enter a single AFDB homodimer accession ID | High |
| UC-02 | Fetch all required structure, PAE, pLDDT, and metadata files automatically | High |
| UC-03 | Detect interface residues and contact pairs | High |
| UC-04 | Recompute all confidence metrics from scratch | High |
| UC-05 | Visualise PAE regions used by each score | High |
| UC-06 | Display per-residue score contribution profiles | High |
| UC-07 | Map confidence and interface data onto a 3D structure | Medium |
| UC-08 | Generate a plain-language diagnostic summary | High |
| UC-09 | Validate notebook behavior against reference accessions | Medium |

---

# 9. User Stories

## User Story Template

> As a [user type], I want to [action], so that [benefit].

---

## Example

> As a non-specialist AFDB user, I want to enter one homodimer accession and see why confidence metrics disagree, so that I can judge whether the predicted complex deserves further inspection.

### Acceptance Criteria

- The notebook accepts a single accession ID.
- All data is fetched automatically from the AFDB API.
- Each score is explained before it is computed.
- The output includes score values, intermediate calculations, visualisations, and a diagnostic statement.

---

## Additional User Stories

> As an AFDB bioinformatician, I want score calculations to match the IPSAE reference implementation, so that notebook results can be trusted during internal review.

> As a platform engineer, I want the notebook to avoid GPU-specific and heavy bioinformatics dependencies, so that it runs reliably in Google Colab.

> As a research collaborator, I want interface residues highlighted in 2D and 3D, so that I can connect confidence metrics to structural features.

---

# 10. Functional Requirements

| ID | Requirement | Priority |
|---|---|---|
| FR-01 | System shall accept one AFDB homodimer accession ID as user input | High |
| FR-02 | System shall call `GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}` | High |
| FR-03 | System shall download mmCIF, PAE JSON, and confidence JSON resources from AFDB response fields | High |
| FR-04 | System shall parse chain assignments, residue numbers, residue names, CA coordinates, and CB coordinates from mmCIF | High |
| FR-05 | System shall parse the PAE matrix and chain metadata from PAE JSON | High |
| FR-06 | System shall parse per-residue pLDDT values from confidence JSON or B-factors | High |
| FR-07 | System shall compute CB-CB inter-chain contacts using an 8.0 A cutoff, using CA for glycine | High |
| FR-08 | System shall compute ipTM_d0chn without a PAE cutoff | High |
| FR-09 | System shall compute ipSAE_d0res, ipSAE_d0chn, and ipSAE_d0dom using PAE cutoff 10 | High |
| FR-10 | System shall compute pDockQ from contact count and interface mean pLDDT | High |
| FR-11 | System shall compute pDockQ2 from interface pLDDT and PAE-transformed contact confidence | High |
| FR-12 | System shall compute LIS using inter-chain PAE values below 12 | High |
| FR-13 | System shall display intermediate calculation values for each score | High |
| FR-14 | System shall show per-residue ipTM and ipSAE contribution profiles | High |
| FR-15 | System shall generate full PAE heatmap and per-score PAE region overlays | High |
| FR-16 | System shall compare interface and non-interface pLDDT distributions | Medium |
| FR-17 | System shall generate at least chain overview and pLDDT MolViewSpec views | Medium |
| FR-18 | System shall generate a plain-language score agreement and disagreement summary | High |

---

# 11. Non-Functional Requirements

## Performance

- Upload throughput: not applicable; all data is fetched from AFDB URLs.
- Max concurrent workflows: one accession per notebook run.
- API latency: notebook should tolerate normal AFDB request latency with clear failure messages.
- Runtime target: less than 60 seconds on free-tier Google Colab without GPU.

---

## Scalability

- Expected users: individual researchers, AFDB reviewers, and collaborators running ad hoc notebook sessions.
- Expected datasets: one homodimer prediction per notebook execution.
- Storage growth/year: none required by the notebook; generated `.mvsj` files are local session artifacts.

---

## Reliability

- Uptime target: dependent on AFDB API and Colab runtime availability.
- Backup strategy: notebook is self-contained; no persistent output store required.
- Disaster recovery objectives: not applicable for v1 exploratory notebook.
- Error handling: missing API fields, malformed JSON, unavailable downloads, and non-homodimer structures should produce clear notebook messages.

---

## Security

- No user authentication required.
- No secrets or credentials stored in the notebook.
- Data fetched only from AFDB-provided public URLs.
- Notebook should avoid executing external scripts from downloaded data.
- Generated MolViewSpec links should not include sensitive information.

---

## Compliance

- RUO framing.
- No clinical diagnosis claims.
- FAIR Principles supported through transparent data provenance and reusable notebook outputs.
- GDPR, HIPAA, GxP, MDR, and IVDR are not in scope for v1 unless future deployment context changes.

---

# 12. Data Requirements

## Data Types

- Predicted protein structure coordinates.
- Predicted aligned error matrix.
- Per-residue pLDDT confidence scores.
- AFDB prediction metadata.
- Protein sequence and chain ranges.
- Derived interface contact and per-residue score data.

---

## Supported Formats

| Data Type | Formats |
|---|---|
| Structure | mmCIF, optional PDB or BinaryCIF for visualisation |
| PAE | AFDB PAE JSON |
| Confidence | AFDB confidence JSON, mmCIF B-factor pLDDT fallback |
| Metadata | AFDB prediction API JSON |
| Visualisation State | MolViewSpec `.mvsj` JSON |

---

## Metadata Standards

- AFDB prediction API response fields.
- UniProt accession cross-reference.
- mmCIF `_atom_site` records.
- Chain labels using `label_asym_id`.
- Residue ranges using `sequenceStart` and `sequenceEnd`.
- FAIR Principles for provenance and interpretability.

---

# 13. Architecture Overview

## System Components

- Notebook input cell for AFDB accession ID.
- AFDB API client using `requests`.
- Data download and validation layer.
- mmCIF parser for coordinate, residue, and chain extraction.
- PAE and pLDDT parsers.
- Interface detection module.
- Metric computation functions.
- Plotting and table generation functions.
- MolViewSpec state generation functions.
- Diagnostic summary heuristic layer.

---

## Integrations

| System | Purpose |
|---|---|
| AlphaFold Database API | Retrieve prediction metadata and file URLs |
| AFDB mmCIF / BinaryCIF files | Provide structural coordinates for contacts and 3D views |
| AFDB PAE JSON | Provide PAE matrix for ipTM, ipSAE, pDockQ2, and LIS |
| AFDB confidence JSON | Provide per-residue pLDDT values |
| MolViewSpec / Mol* | Render interactive 3D molecular views |
| Google Colab | Target notebook runtime environment |

---

# 14. Workflow & UX

## Key Workflows

1. User enters a homodimer accession ID.
2. Notebook fetches AFDB prediction metadata.
3. Notebook downloads mmCIF, PAE JSON, and confidence JSON.
4. Notebook parses residues, chains, coordinates, PAE, and pLDDT.
5. Notebook identifies interface contacts and interface residues.
6. Notebook displays contact map and interface coverage.
7. Notebook shows full PAE matrix and score-specific PAE masks.
8. Notebook computes and decomposes all scores.
9. Notebook visualises pLDDT at interface and non-interface regions.
10. Notebook generates MolViewSpec views.
11. Notebook displays a diagnostic summary and score agreement matrix.

---

## UX Requirements

- Every code cell must be preceded by a markdown explanation cell.
- Explanations must define structural bioinformatics terms before using them.
- Visualisations must include titles, legends, labelled axes, units, and interpretation notes.
- Colour schemes must remain consistent across plots.
- The summary panel must state what is high, moderate, low, and why.
- The notebook must avoid jargon without definition.
- The notebook must be usable by a newcomer without reading source code.

---

# 15. AI / ML Components (Optional)

## AI Use Cases

- No generative AI component is required for v1.
- The notebook interprets existing AlphaFold confidence outputs using deterministic formulas and rule-based heuristics.
- Future extensions may use ranking or summarisation methods to identify metric disagreement outliers across many homodimers.

---

## Model Requirements

| Requirement | Details |
|---|---|
| Explainability | Deterministic formulas and visible intermediate values are required |
| Validation datasets | At least three AFDB homodimer accessions, including high-confidence, borderline, and disagreement cases |
| Monitoring | Not applicable for notebook-only v1 |
| Bias mitigation | Not applicable for deterministic metric decomposition; interpretive text must avoid definitive biological or clinical claims |

---

# 16. Risks & Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Score mismatch with IPSAE reference implementation | High | Implement formulas directly from IPSAE v4 and validate within +/- 0.001 |
| Confusion between CA-CA and CB-CB interface definitions | Medium | Document that notebook uses CB-CB with CA for glycine to match IPSAE scoring, while AFDB production pipeline uses CA-CA |
| AFDB API response changes | Medium | Validate required fields and produce clear error messages |
| MolViewSpec integration fails in a notebook runtime | Medium | Generate `.mvsj` files and provide Mol* viewer links as fallback |
| Non-specialist users overinterpret confidence as biological truth | High | Use RUO framing and describe diagnostic statements as interpretive guidance |
| Runtime exceeds Colab target | Medium | Avoid BioPython, PyTorch, torch-geometric, and GPU-only dependencies |
| PAE matrix asymmetry is mishandled | High | Treat A->B and B->A quadrants separately where formulas require asymmetry |

---

# 17. Dependencies

- Python notebook environment, preferably Google Colab.
- `numpy`.
- `matplotlib`.
- `seaborn`.
- `requests`.
- `molviewspec`.
- AFDB API and public AFDB file URLs.
- DunbrackLab IPSAE formulas for reference validation.
- No BioPython, PyTorch, torch-geometric, GPU, external scripts, or local data files.

---

# 18. Milestones & Timeline

| Milestone | Owner | Target Date |
|---|---|---|
| Discovery complete | AFDB / InsightFold team | 2026-03-25 |
| MVP ready | Notebook implementation owner | TBD |
| Reference validation complete | Notebook implementation owner | TBD |
| Beta release | AFDB / PDBe review group | TBD |
| Production launch | TBD | TBD |

---

# 19. KPIs & Success Metrics

| Metric | Target |
|---|---|
| Score recomputation accuracy | Matches IPSAE reference within +/- 0.001 where comparable |
| Notebook runtime | Less than 60 seconds on free-tier Colab |
| Required scores computed | 100% of ipSAE variants, ipTM_d0chn, pDockQ, pDockQ2, and LIS |
| Required visualisations generated | Contact map, interface coverage, PAE heatmap, PAE masks, per-residue profiles, pLDDT comparison |
| MolViewSpec minimum views | At least 2 views: chain overview and pLDDT mapping |
| Educational coverage | Every code cell has preceding explanatory markdown |
| Test case coverage | At least 3 AFDB accessions validated |
| Dependency compliance | No BioPython, PyTorch, torch-geometric, or GPU-only libraries |

---

# 20. Open Questions

- Which AFDB accessions should be used as the borderline and disagreement validation cases?
- Should pLDDT be sourced preferentially from confidence JSON or mmCIF B-factors when both are available?
- What threshold bands should define high, moderate, and low for each displayed score?
- Should MolViewSpec interface and disagreement views be required for v1 or remain stretch goals?
- Should generated `.mvsj` files be saved alongside the notebook or only displayed inline?
- How should non-homodimer or malformed accession inputs be handled in the user interface?

---

# Appendix

## Glossary

| Term | Definition |
|---|---|
| AFDB | AlphaFold Database |
| PAE | Predicted aligned error; an AlphaFold estimate of positional uncertainty between residue pairs |
| pLDDT | Per-residue confidence score for local structure quality |
| ipTM | Inter-chain predicted TM-score style confidence metric |
| ipSAE | Interface predicted self-aligned error metric used for inter-chain confidence assessment |
| pDockQ | Contact and pLDDT-based score for predicted complex interface plausibility |
| pDockQ2 | Interface score combining contact evidence with PAE confidence |
| LIS | Local Interaction Score based on low-PAE inter-chain regions |
| mmCIF | Macromolecular Crystallographic Information File format |
| MolViewSpec | Specification and Python tooling for generating Mol* molecular views |
| RUO | Research Use Only |
| FAIR | Findable, Accessible, Interoperable, Reusable |

---

## References

- Source summary: `specs/homodimer_diagnostic_notebook_spec.md`
- AlphaFold Database prediction endpoint: `GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}`
- DunbrackLab IPSAE repository: https://github.com/DunbrackLab/IPSAE
- IPSAE paper: https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2
- Bryant, Pozzati, and Elofsson (2022), pDockQ: https://www.nature.com/articles/s41467-022-28865-w
- Zhu, Shenoy, Kundrotas, and Elofsson (2023), pDockQ2: https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714
- Kim, Hu, Comjean, Rodiger, Mohr, and Perrimon (2024), LIS: https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1
- MolViewSpec GitHub: https://github.com/molstar/mol-view-spec
- MolViewSpec documentation: https://molstar.org/mol-view-spec-docs/
- Midlik et al. (2025), MolViewSpec NAR paper: https://doi.org/10.1093/nar/gkaf370
- Yang and Skolnick (2004), TM-score d0 formula origin.
