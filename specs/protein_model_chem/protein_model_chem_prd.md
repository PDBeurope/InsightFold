# Product Requirements Document (PRD)
## Protein Model Chemistry Structural Reasoning Notebook

---

# 1. Document Overview

| Field | Details |
|---|---|
| Product Name | Protein Model Chemistry Structural Reasoning Notebook |
| Version | 0.1 |
| Date | 2026-05-08 |
| Author(s) | Codex, based on `protein_model_chem_summary.txt` |
| Stakeholders | Structural biologists, protein scientists, translational researchers, mutation-focused experimentalists, InsightFold maintainers |
| Status | Draft |
| Related Documents | `specs/protein_model_chem/protein_model_chem_summary.txt` |
| Regulatory Scope | RUO; non-clinical; not intended for diagnosis, pathogenicity classification, treatment decisions, or regulated clinical reporting |

---

# 2. Executive Summary

## Problem Statement

Structural biologists and protein scientists routinely inspect AlphaFold models, PDBe experimental structures, residue contacts, and mutation sites, but current workflows are often manual, fragmented, visualization-heavy, and weakly interpretable. Many tools show contacts or prediction scores without explaining the structural chemistry behind those signals, while others risk overclaiming mechanistic certainty.

This notebook should help users answer:

> How do mutations perturb the structural interaction and energetic constraint network underlying a protein's shape, and what mechanistic structural consequences may result?

The desired output is not merely structure visualization. The notebook must provide interpretable mechanistic reasoning grounded in structural chemistry, while clearly separating directly computed observations from cautious interpretations and higher-level hypotheses.

---

## Product Vision

The notebook will provide an InsightFold workflow for analyzing protein structures as structural constraint networks. It will retrieve or load AlphaFold and PDBe structures, detect chemically meaningful residue interactions, model mutation-induced perturbations to local compatibility and interaction topology, and explain likely structural consequences without claiming exact thermodynamic or clinical truth.

The long-term goal is to occupy a defensible niche between structure viewers, mutation predictors, and black-box pathogenicity tools: an educational and research-oriented system for interpretable structural reasoning.

---

## Success Criteria

- Users can load an AlphaFold or PDBe structure and inspect provenance, local coordinate confidence, and structural caveats.
- Users can identify residue-level interaction networks, including hydrogen bonds, salt bridges, hydrophobic contacts, aromatic interactions, disulfides, steric clashes, and local packing contacts.
- Users can introduce mutations and compare wild-type versus mutant interaction environments.
- The notebook produces layered outputs: observation, interpretation, hypothesis, structural disruption likelihood, and confidence.
- Expert users judge generated explanations as chemically reasonable and useful for experimental planning.
- The notebook avoids clinical, pathogenicity, exact ddG, and deterministic folding claims.
- The notebook can later be converted into an implementation-ready notebook spec without inventing core data sources, outputs, risks, or acceptance criteria.

---

# 3. Background & Context

## Scientific Context

- Biological domain: protein structural biology, structural chemistry, mutation interpretation, and residue interaction networks.
- Structural inputs: AlphaFold predicted structures and PDBe experimental structures.
- Core abstraction: a protein fold is treated as a network of residue-level structural constraints stabilized by cooperative energetic compatibility.
- Mutation model: a mutation is interpreted as a perturbation to local compatibility, packing, sterics, interaction topology, and constraint density.
- Mechanistic categories: secondary structure disruption, interface disruption, flexibility-related constraint changes, and structural destabilization signals.
- Scientific limits: static structures do not capture full thermodynamics, entropy, kinetics, allostery, ensemble dynamics, folding pathways, or biological pathogenicity.

---

## Business / Institutional Context

- Aligns with InsightFold's notebook-driven life-sciences analysis lifecycle.
- Supports research use cases where users need transparent structural reasoning rather than opaque mutation predictions.
- Provides a foundation for future structural chemistry notebooks, fixture selection, validation against known mutation effects, and educational workflows.
- Should remain positioned as a research and hypothesis-generation workflow until expert validation supports broader operationalization.

---

## Existing Solutions

| Solution | Strengths | Weaknesses |
|---|---|---|
| Molecular viewers such as PyMOL, ChimeraX, Mol* | Strong visual inspection and manual structural analysis | Require expertise, provide limited automated mechanistic interpretation, and can be time-consuming for mutation comparison |
| AlphaFold and AFDB confidence views | Provide broad structure availability and useful local confidence signals | Confidence is not equivalent to biological certainty; weak for disorder, interfaces, heterogeneity, ligands, and assemblies |
| PDBe experimental structure records | Provide experimentally derived structural evidence and metadata | May include missing residues, engineered constructs, crystal-packing artifacts, low-resolution sidechains, or incomplete biological assemblies |
| Mutation/pathogenicity predictors | Useful for prioritization in some biological contexts | Often mix evidence modalities and can be opaque; not suitable as a structural-mechanism explanation layer |
| Force-field or molecular dynamics workflows | Can model physical behavior in more detail | Expensive, complex, sensitive to setup, and outside the intended defensible MVP scope |

---

# 4. Goals & Objectives

## Product Goals

| Goal | Metric |
|---|---|
| Enable interpretable structural reasoning from AF and PDBe structures | Users can trace each conclusion to interactions, local confidence, and graph-derived signals |
| Reduce manual residue and mutation inspection burden | Expert users report faster hypothesis generation than manual viewer-only workflows |
| Support mutation perturbation analysis without overclaiming | Outputs consistently avoid exact ddG, pathogenicity, clinical severity, and deterministic folding language |
| Make confidence and provenance visible throughout the workflow | Every interaction, interpretation, and severity ranking carries source and confidence context |
| Prepare for implementation-ready notebook specification | PRD includes target users, data needs, functional requirements, risks, non-goals, and open questions |

---

## Non-Goals

Clearly define what this product will **not** address.

- It will not predict clinical pathogenicity, disease severity, diagnosis, treatment impact, or regulatory-grade variant interpretation.
- It will not compute rigorous thermodynamic free energies, validated ddG values, or total protein stability scores.
- It will not perform molecular dynamics, full ensemble modeling, folding pathway simulation, or allosteric propagation prediction.
- It will not produce authoritative mutant structures or claim exact mutant conformations.
- It will not include V1 ligand energetics, cofactors, PTMs, nucleic acids, membrane energetics, or explicit solvent models.
- It will not rely on conservation, genomics, clinical databases, or evolutionary scores for severity ranking in the structural-only MVP.

---

# 5. Target Users

## Primary Users

| User Type | Goals | Pain Points |
|---|---|---|
| Structural biologists | Understand residue interactions, local packing, interfaces, and mutation perturbations | Manual inspection is slow; confidence and mechanistic interpretation are often spread across tools |
| Protein scientists | Identify structurally important residues and generate testable mutation hypotheses | Need interpretable chemistry rather than opaque prediction scores |
| Translational researchers | Explore structural mechanisms behind mutations for research prioritization | Need non-clinical mechanistic context without overclaiming pathogenicity |
| Mutation-focused experimentalists | Prioritize constructs, mutagenesis experiments, or validation assays | Need actionable structural hypotheses without coding or raw force-field interpretation |
| InsightFold maintainers | Provide reproducible notebook workflows and validation fixtures | Need clear scope boundaries, data contracts, and scientific risk tracking |

---

## Secondary Users

- Educators teaching protein structure, mutation effects, and interaction networks.
- Research collaborators reviewing figures or mechanistic explanations.
- Bioinformatics developers extending structural chemistry workflows.
- Domain reviewers validating scientific interpretation and confidence language.

---

# 6. User Personas

## Persona 1 - Structural Biology Researcher

### Background
- Understands protein folds, sidechain chemistry, interfaces, and structure quality.
- Uses molecular viewers and databases, but may not want to write custom analysis code for every mutation.

### Needs
- Rapid comparison of wild-type and mutant interaction environments.
- Clear provenance, confidence, and structural caveats.
- Exportable figures and explanations for research discussions.

### Frustrations
- Manual contact inspection is repetitive.
- Existing tools often show contacts without ranking their mechanistic importance.
- Mutation outputs can sound more certain than the structure supports.

---

## Persona 2 - Protein Engineering Experimentalist

### Background
- Designs or tests protein variants and understands basic structural biology.
- Uses structural hypotheses to plan experiments.

### Needs
- Identification of residues that stabilize helices, beta sheets, buried cores, or interfaces.
- Mutation disruption likelihood that remains structural and non-clinical.
- Plain-language interpretation tied to visible evidence.

### Frustrations
- Raw force-field outputs are difficult to interpret.
- Black-box predictors provide scores without explaining local chemistry.
- Structural data may be low confidence, incomplete, or assembly-dependent.

---

## Persona 3 - InsightFold Notebook Developer

### Background
- Builds reproducible scientific notebooks from PRDs and specs.
- Needs explicit data contracts, fixtures, validation criteria, and implementation boundaries.

### Needs
- Clear required inputs, outputs, visualizations, and success metrics.
- Explicit unresolved choices for scoring, graph weighting, and calibration.
- Risks and non-goals that prevent unsafe scope expansion.

### Frustrations
- Scientific ambiguity can lead to hidden assumptions in implementation.
- Validation is difficult without fixture expectations and confidence rules.

---

# 7. Product Scope

## In Scope

- AlphaFold and PDBe structure retrieval or loading.
- Structure normalization across chains, residues, numbering, missing regions, alternate conformations, and provenance.
- Structure visualization with confidence and provenance indicators.
- Detection of residue-level interaction types: hydrogen bonds, salt bridges, hydrophobic contacts, aromatic interactions, disulfides, steric clashes, and local packing contacts.
- Residue-level interaction graph construction and network metrics.
- Mutation introduction and local environment recomputation.
- Wild-type versus mutant comparison of interactions, packing, clashes, compatibility, and network topology.
- Layered mechanistic explanation: observation, interpretation, and hypothesis.
- Structural disruption likelihood ranking: low, moderate, or high.
- Confidence propagation from structure source, local structural confidence, interaction robustness, and inference depth.
- Warnings for unsupported contexts such as disorder, low-confidence AF regions, membrane proteins, metal-dependent proteins, missing ligands, flexible loops, and uncertain assemblies.
- Exportable results for figures, tables, and explanatory summaries.

---

## Out of Scope

- Clinical diagnosis, clinical reporting, pathogenicity prediction, disease severity, or patient-specific interpretation.
- Molecular dynamics, ensemble thermodynamics, exact conformational prediction, folding free-energy prediction, or rigorous ddG calculation.
- Ligand, cofactor, PTM, nucleic acid, solvent, metal coordination, and membrane energetics in V1.
- Automated integration of conservation, genomic evidence, clinical databases, literature evidence, or evolutionary scores into severity ranking.
- Production-scale multi-user platform features such as authentication, RBAC, audit logging, and billing.
- Wet-lab automation or experimental protocol generation.

---

# 8. Use Cases

| ID | Use Case | Priority |
|---|---|---|
| UC-01 | Load an AlphaFold or PDBe protein structure and inspect provenance and confidence | High |
| UC-02 | Compute and visualize residue-level interaction networks | High |
| UC-03 | Inspect a residue's local interaction environment and network role | High |
| UC-04 | Introduce a mutation and compare wild-type versus mutant local compatibility | High |
| UC-05 | Identify lost interactions, steric clashes, packing disruption, and centrality shifts | High |
| UC-06 | Generate layered mechanistic explanations and structural hypotheses | High |
| UC-07 | Rank structural disruption likelihood as low, moderate, or high | Medium |
| UC-08 | Export figures, interaction tables, confidence notes, and written explanations | Medium |
| UC-09 | Warn users when structure quality or biological context limits interpretation | High |
| UC-10 | Use notebook outputs for experimental planning or educational explanation | Medium |

---

# 9. User Stories

## User Story Template

> As a [user type], I want to [action], so that [benefit].

---

## Story 1 - Structure Loading And Provenance

> As a structural biologist, I want to load an AF or PDBe structure and see provenance and confidence annotations, so that I can judge whether downstream interpretations are reliable.

### Acceptance Criteria

- User can provide a structure accession or supported structure file.
- The notebook records whether the structure came from AlphaFold, PDBe, or user upload.
- PDBe is preferred by default when a high-quality experimental structure exists, while manual override remains available.
- pLDDT or experimental quality metadata is surfaced where available.
- Low-confidence or incomplete regions are visibly marked and carried into downstream confidence.

---

## Story 2 - Interaction Network Inspection

> As a protein scientist, I want to inspect residue-level interaction networks, so that I can understand which interactions stabilize a protein's local and global shape.

### Acceptance Criteria

- The notebook detects required interaction classes and exposes the definitions or geometric thresholds used.
- Interactions are summarized as residue-level tables and graph edges.
- Users can inspect local residue neighborhoods, interaction density, and centrality-like signals.
- Interaction confidence is separated from structure confidence and mechanistic confidence.

---

## Story 3 - Mutation Perturbation Analysis

> As a mutation-focused experimentalist, I want to introduce a residue substitution and compare it to wild type, so that I can generate mechanistic hypotheses for experimental testing.

### Acceptance Criteria

- User can specify at least one mutation with chain, residue, original residue, and substituted residue.
- The notebook recomputes or estimates the local mutant interaction environment.
- Outputs identify lost interactions, gained clashes, packing changes, compatibility shifts, and local network perturbations.
- The notebook does not claim to predict the exact mutant structure.

---

## Story 4 - Mechanistic Explanation

> As a translational researcher, I want explanations that distinguish observations from hypotheses, so that I can use structural reasoning without mistaking it for clinical interpretation.

### Acceptance Criteria

- Every result is labeled as observation, interpretation, or hypothesis.
- Language uses cautious terms such as "suggests", "consistent with", "may require", and "possible destabilization".
- Language avoids deterministic terms such as "causes", "abolishes", and "guarantees".
- Structural disruption likelihood is explicitly non-clinical and non-pathogenicity.

---

## Story 5 - Export And Review

> As an InsightFold user, I want to export figures and concise explanations, so that I can include them in lab discussions, presentations, or follow-up notebook review.

### Acceptance Criteria

- Export includes structure provenance, interaction definitions, mutation input, key observations, confidence notes, and caveats.
- Exported summaries preserve the distinction between computed facts and speculative hypotheses.
- Figures or tables include enough context for domain review.

---

# 10. Functional Requirements

| ID | Requirement | Priority |
|---|---|---|
| FR-01 | Notebook shall support AlphaFold and PDBe structure inputs. | High |
| FR-02 | Notebook shall prefer suitable PDBe experimental structures when available, with manual override. | Medium |
| FR-03 | Notebook shall normalize chains, residue numbering, missing residues, alternate conformations, and structure provenance before analysis. | High |
| FR-04 | Notebook shall display structure visualization with confidence and provenance indicators. | High |
| FR-05 | Notebook shall detect hydrogen bonds, salt bridges, hydrophobic contacts, aromatic interactions, disulfides, steric clashes, and local packing contacts. | High |
| FR-06 | Notebook shall represent interactions as residue-level graph edges with interaction type, confidence, and optional weight fields. | High |
| FR-07 | Notebook shall compute network-derived signals such as interaction density, local connectivity loss, centrality, fragmentation, and constraint redundancy where scientifically defensible. | Medium |
| FR-08 | Notebook shall allow users to inspect a residue's local interaction environment, burial or packing context, compatibility signals, and network role. | High |
| FR-09 | Notebook shall accept mutation input with chain, residue index, original residue, and substituted residue. | High |
| FR-10 | Notebook shall estimate mutation-induced local compatibility perturbation without claiming exact mutant conformations. | High |
| FR-11 | Notebook shall compare wild-type and mutant interaction topology, including lost interactions, gained clashes, altered packing, centrality shifts, and compatibility maps. | High |
| FR-12 | Notebook shall categorize results as observation, interpretation, or hypothesis. | High |
| FR-13 | Notebook shall rank structural disruption likelihood as low, moderate, or high, based only on structural interactions, energetic approximations, and network signals. | Medium |
| FR-14 | Notebook shall propagate separate structure confidence, interaction confidence, and mechanistic confidence through outputs. | High |
| FR-15 | Notebook shall warn users when low-confidence regions, missing biological context, disorder, uncertain assemblies, or unsupported protein classes limit interpretation. | High |
| FR-16 | Notebook shall expose interaction definitions, thresholds, scoring assumptions, and caveats. | High |
| FR-17 | Notebook shall support export of figures, interaction tables, mutation comparison summaries, and confidence/caveat notes. | Medium |
| FR-18 | Notebook shall avoid exact energy, ddG, pathogenicity, disease severity, and clinical language in all outputs. | High |

---

# 11. Non-Functional Requirements

## Performance

- Structure retrieval and parsing target: TBD; should be suitable for interactive notebook use on representative single-chain and small multimer examples.
- Interaction detection target: TBD; should complete without long-running infrastructure for selected fixtures.
- Mutation comparison target: TBD; should be fast enough to support iterative residue inspection in notebook workflows.
- Export generation target: TBD.

---

## Scalability

- Expected users: research notebook users running local or hosted InsightFold workflows.
- Expected datasets: individual proteins or protein assemblies, not proteome-scale batch processing in V1.
- Storage growth/year: TBD; notebook should favor cacheable external structure retrieval and reproducible fixture data.

---

## Reliability

- Notebook must be restart-and-run-all compatible once implemented.
- Results should be reproducible for pinned structure fixtures and explicit thresholds.
- External API failures should be handled with clear messages and cached fixture fallbacks where available.
- Downstream interpretation must degrade when confidence is low rather than silently producing authoritative output.

---

## Security

- No sensitive human subject data is expected in V1.
- Notebook should not require clinical, protected health, or regulated variant datasets.
- If user-uploaded structures are supported later, storage and sharing assumptions must be revisited.

---

## Compliance

- Regulatory scope is RUO.
- Notebook outputs must not be framed as clinical decision support.
- FAIR-style provenance should be preserved for structures, accessions, versions, fixtures, and generated outputs.

---

# 12. Data Requirements

## Data Types

- Protein structures from AlphaFold and PDBe.
- Optional user-provided structure files.
- Optional mutation inputs with chain, residue, wild-type residue, and mutant residue.
- Structure provenance and confidence metadata.
- Residue-level interaction records.
- Residue graph and network metric outputs.
- Mechanistic interpretation and confidence summaries.

---

## Supported Formats

| Data Type | Formats |
|---|---|
| Protein structures | mmCIF and/or PDB; exact supported formats TBD during notebook specification |
| AlphaFold confidence | pLDDT and available model confidence metadata |
| PDBe metadata | Experimental method, resolution where applicable, chain/residue mapping, assembly metadata where available |
| Mutation inputs | Structured notebook fields; optional CSV/TSV batch format TBD |
| Interaction tables | DataFrame-compatible tabular format; CSV export desirable |
| Visual outputs | Notebook-rendered structure views, heatmaps, network views, and exportable static figures |
| Explanation outputs | Markdown and/or structured JSON-like summary TBD |

---

## Metadata Standards

- Structure accession, source database, retrieval date, model version where available, and file checksum where practical.
- Chain identifiers, residue numbering, insertion codes, missing residues, alternate conformations, and biological assembly assumptions.
- AlphaFold local confidence such as pLDDT where available.
- PDBe experimental metadata such as method and resolution where applicable.
- Interaction thresholds, scoring parameters, graph weighting assumptions, and confidence propagation rules.
- Mutation notation and validation status.

---

# 13. Architecture Overview

## System Components

- Notebook user interface and narrative sections.
- Structure retrieval layer for AlphaFold and PDBe.
- Structure parser and normalization layer.
- Provenance and confidence annotation layer.
- Interaction detection layer.
- Residue-level graph construction and network analysis layer.
- Mutation application or perturbation estimation layer.
- Wild-type versus mutant comparison layer.
- Mechanistic explanation and confidence classification layer.
- Visualization layer for structures, interactions, networks, and compatibility shifts.
- Export layer for figures, tables, and written summaries.

---

## Integrations

| System | Purpose |
|---|---|
| AlphaFold / AFDB | Retrieve predicted structures and local confidence metadata |
| PDBe | Retrieve experimental structures, metadata, and assembly context |
| Notebook runtime | Execute reproducible analysis and display educational narrative |
| Structure visualization library | Render structures, residue selections, interactions, and confidence overlays |
| Scientific Python stack | Parse structures, compute geometry, build graphs, tabulate outputs, and plot figures |

---

# 14. Workflow & UX

## Key Workflows

1. Load structure from AlphaFold, PDBe, or supported file input.
2. Normalize structure, annotate provenance, and show confidence indicators.
3. Compute residue-level interaction types and interaction graph.
4. Inspect global interaction network and local residue environments.
5. Select residue or enter mutation.
6. Estimate mutant local compatibility and interaction-network perturbation.
7. Compare wild-type versus mutant observations.
8. Generate layered mechanistic interpretation and structural disruption likelihood.
9. Review confidence, caveats, and unsupported-context warnings.
10. Export figures, tables, and concise explanations.

---

## UX Requirements

- Notebook should lead with the usable analysis workflow, not a marketing-style overview.
- Inputs should be explicit and validated before analysis.
- Results should be organized by evidence depth: observation, interpretation, hypothesis.
- Confidence language and caveats should appear near the relevant result, not only in a final disclaimer.
- Visualization should emphasize interaction perturbation maps, network changes, compatibility shifts, and local structural context.
- The notebook should avoid visual or textual framing that implies clinical or deterministic conclusions.
- Users should be able to inspect how an interaction or severity label was derived.
- Warnings should be specific: low pLDDT, missing residues, uncertain assembly, missing ligands, unsupported protein context, or weak interaction confidence.

---

# 15. AI / ML Components (Optional)

## AI Use Cases

- V1 should rely on transparent structural chemistry heuristics, graph analysis, and confidence propagation rather than opaque ML prediction.
- Natural-language mechanistic explanation may be generated from structured observations and interpretations, but it must be constrained by computed evidence.
- Future ML use could include calibrated structural disruption scoring if validated datasets and expert review support it.

---

## Model Requirements

| Requirement | Details |
|---|---|
| Explainability | Every explanation must trace back to structure source, interaction evidence, thresholds, and confidence classification |
| Validation datasets | TBD; should include known structural mutation effects and expert-reviewed examples before calibration claims |
| Monitoring | Not applicable for static notebook MVP; later hosted workflows would need output quality and misuse monitoring |
| Bias mitigation | Avoid training or claims based on clinical labels in V1; avoid mixing structural-only reasoning with pathogenicity evidence |

---

# 16. Risks & Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Users overinterpret static structures as biological certainty | High | Use layered confidence, cautious language, provenance display, and explicit structural-hypothesis framing |
| AlphaFold confidence is mistaken for biological truth | High | Propagate pLDDT, downweight low-confidence regions, visually mark uncertainty, and warn for interfaces/disorder/heterogeneity |
| PDBe structures are treated as absolute truth | Medium | Surface resolution, missing residues, constructs, assemblies, and experimental caveats |
| Heuristic energetic estimates imply false precision | High | Avoid exact energy and ddG language; use qualitative compatibility and destabilization terminology |
| Network metrics look convincing but lack biological meaning | High | Expose graph definitions and weights; validate against expert examples; avoid unsupported centrality claims |
| Mutation interpretation overclaims exact mutant conformation | High | Frame outputs as local compatibility perturbations and structural adaptation pressure |
| Missing ligands, cofactors, PTMs, metals, membranes, or assemblies distort interpretation | High | Detect or warn when context may be missing; exclude unsupported claims in V1 |
| Structure normalization errors produce incorrect residue mapping | High | Build robust normalization, validate mutation notation, and include fixture tests for numbering and chains |
| Explanations feel decorative rather than actionable | Medium | Tie every explanation to computed observations and expert-reviewable evidence |
| Interaction visualizations become cognitively overwhelming | Medium | Prioritize residue-centric summaries, filters, and focused comparison views |

---

# 17. Dependencies

- AlphaFold / AFDB structure availability and confidence metadata.
- PDBe structure and metadata access.
- Reliable structure parsing and normalization libraries.
- Geometry and chemistry rules for interaction detection.
- Network analysis tooling.
- Structure visualization support in notebook environments.
- Fixture set covering AF structures, PDBe structures, missing residues, low-confidence regions, interfaces, and representative mutations.
- Domain expert review for interaction definitions, confidence wording, and mechanistic interpretation.
- Future validation datasets for known mutation effects and disruption calibration.

---

# 18. Milestones & Timeline

| Milestone | Owner | Target Date |
|---|---|---|
| PRD draft complete | Codex | 2026-05-08 |
| Notebook spec pack drafted | TBD | TBD |
| Fixtures selected and documented | TBD | TBD |
| Interaction detection prototype complete | TBD | TBD |
| Mutation perturbation workflow prototype complete | TBD | TBD |
| Visualization and explanation layers implemented | TBD | TBD |
| Restart-and-run-all validation complete | TBD | TBD |
| Domain review complete | TBD | TBD |

---

# 19. KPIs & Success Metrics

| Metric | Target |
|---|---|
| Expert-rated chemical reasonableness of explanations | TBD; should be assessed during domain review |
| Successful restart-and-run-all execution for pinned fixtures | 100% for release candidate |
| Required interaction classes detected in fixture outputs | 100% where fixture biology supports the class |
| Outputs preserving observation/interpretation/hypothesis labels | 100% of generated mechanistic summaries |
| Outputs avoiding clinical/pathogenicity/exact-energy claims | 100% of notebook text and generated summaries |
| User ability to trace severity to structural evidence | TBD; evaluate in usability review |
| Manual inspection effort reduction | TBD; compare expert workflow with and without notebook |
| Correlation with known structural mutation effects | TBD; requires validation dataset |

---

# 20. Open Questions

- Which energetic framework should power interaction scoring: geometric heuristics, statistical potentials, empirical force fields, or a hybrid?
- How should interaction types contribute to graph edge weighting?
- Should mutant analysis use sidechain replacement only, local minimization, repacking, or a simpler compatibility heuristic?
- How should low, moderate, and high structural disruption thresholds be calibrated?
- Which benchmark datasets should be used to compare predictions against known structural mutation effects?
- Should the notebook compare multiple structures or models automatically when available?
- How should biological assemblies be reconstructed or selected by default?
- Should interaction persistence across multiple structures or models be computed?
- How detailed should natural-language mechanistic explanations be, and what guardrails prevent unsupported claims?
- What precomputation versus live computation split is acceptable for interactive notebook use?
- What exact structure file formats and parser libraries should be supported in V1?
- How should unsupported contexts such as metal coordination, ligands, membrane proteins, and PTMs be detected or flagged?

---

# Appendix

## Glossary

| Term | Definition |
|---|---|
| AF | AlphaFold, a predicted protein structure resource and model family |
| AFDB | AlphaFold Database |
| PDBe | Protein Data Bank in Europe, a source of experimental macromolecular structures and metadata |
| RUO | Research Use Only |
| pLDDT | AlphaFold per-residue confidence score |
| Structural constraint network | A residue-level graph representation of stabilizing interactions and compatibility constraints within a protein fold |
| Compatibility perturbation | A mutation-induced change in local structural fit, packing, sterics, or interaction favorability |
| Structural disruption likelihood | A non-clinical ranking of how strongly structural evidence suggests local or broader disruption |
| Observation | A directly computed result such as a lost hydrogen bond or detected clash |
| Interpretation | Mechanistic reasoning derived from observations, such as weakened hydrophobic packing |
| Hypothesis | A higher-level, more speculative consequence, such as possible local flexibility increase |
| ddG | Change in folding or binding free energy; exact ddG prediction is out of scope for this PRD |

---

## References

- `specs/protein_model_chem/protein_model_chem_summary.txt`
- AlphaFold / AFDB structure and confidence metadata
- PDBe experimental structure and metadata resources
- Relevant structural biology literature and expert review to be selected during fixture and notebook specification

---

## Assumptions

- Static structures are informative for structural hypothesis generation when confidence and provenance are handled explicitly.
- Local residue interactions can approximate stabilizing structural constraints for educational and research use.
- Interaction-network perturbation is a useful proxy for structural sensitivity, pending validation.
- Target users value interpretable mechanistic explanations more than opaque prediction scores.
- The first implementation should prioritize scientific defensibility over breadth of biological context.

---

## Human / Domain Review Needs

- Review interaction definitions, thresholds, and graph weighting assumptions.
- Review mutation consequence categories and severity wording.
- Review examples where missing biological context could invalidate interpretation.
- Review final notebook outputs for overclaiming, false precision, or clinical/pathogenicity framing.
- Approve validation fixtures before broader notebook lifecycle progression.
