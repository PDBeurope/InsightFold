# Protein Model Chemistry Notebook Requirements

Source PRD: `specs/protein_model_chem/protein_model_chem_prd.md`

## Goal

Build a self-contained InsightFold notebook that loads an AlphaFold or PDBe protein structure, detects chemically meaningful residue interactions, models a user-specified mutation as a local compatibility perturbation, compares wild-type versus mutant interaction networks, and produces layered mechanistic explanations with explicit confidence and RUO limits.

The notebook must remain an interpretable structural chemistry workflow. It must not become a pathogenicity predictor, exact stability engine, molecular dynamics workflow, or clinical decision-support tool.

## Users

- Structural biologists who need faster residue and mutation inspection.
- Protein scientists who need chemically interpretable mutation hypotheses.
- Mutation-focused experimentalists who need structural hypotheses for experimental planning.
- Translational researchers who need non-clinical mechanistic context.
- InsightFold maintainers and reviewers who need reproducible notebook behavior, fixtures, and scientific caveats.

## User Stories

| ID | User Story |
|---|---|
| US-001 | As a structural biologist, I want to load an AlphaFold or PDBe structure and see provenance and confidence annotations, so that I can judge whether downstream interpretations are reliable. |
| US-002 | As a protein scientist, I want to inspect residue-level interaction networks, so that I can understand which interactions stabilize a protein's local and global shape. |
| US-003 | As a mutation-focused experimentalist, I want to introduce a residue substitution and compare it to wild type, so that I can generate mechanistic hypotheses for experimental testing. |
| US-004 | As a translational researcher, I want explanations that distinguish observations from hypotheses, so that I can use structural reasoning without mistaking it for clinical interpretation. |
| US-005 | As an InsightFold user, I want to export figures and concise explanations, so that I can include them in lab discussions, presentations, or follow-up notebook review. |

## Functional Requirements

| ID | Requirement | Acceptance Criteria |
|---|---|---|
| REQ-001 | Accept a structure source input | Notebook exposes a clear input cell for AFDB accession, PDBe/PDB identifier, or supported local structure path; unsupported combinations fail before analysis with a readable message. |
| REQ-002 | Retrieve AlphaFold structures and confidence metadata | AFDB inputs fetch structure coordinates and pLDDT where available; source URL, accession, retrieval date, and confidence provenance are displayed. |
| REQ-003 | Retrieve PDBe structures and experimental metadata | PDBe inputs fetch coordinates and show experimental method, resolution where available, chain IDs, residue ranges, missing-residue caveats, and assembly assumptions. |
| REQ-004 | Prefer PDBe when a suitable experimental structure is available | Default mode records why PDBe or AFDB was selected; users can override source preference; unresolved PDBe suitability is reported rather than hidden. |
| REQ-005 | Parse and normalize protein structure data | Produces a residue table with chain, residue number, insertion code, residue name, atom coordinates, missing atoms, pLDDT or experimental-quality fields, and provenance. |
| REQ-006 | Validate mutation input | Mutation input includes chain, residue identifier, expected wild-type residue, and mutant residue; mismatches, missing residues, unsupported residues, or ambiguous mappings stop mutation analysis clearly. |
| REQ-007 | Detect residue-level interactions | Notebook detects hydrogen bonds, salt bridges, hydrophobic contacts, aromatic interactions, disulfides, steric clashes, and local packing contacts using documented geometric rules. |
| REQ-008 | Build an interaction graph | Interactions are represented as residue-level graph edges with source residue, target residue, interaction type, geometric evidence, confidence, and optional weight. |
| REQ-009 | Compute network-derived signals | Notebook computes local degree, interaction density, connectivity loss, centrality-like values, fragmentation, and constraint redundancy only where definitions are documented. |
| REQ-010 | Inspect a selected residue | User can inspect local interaction environment, burial/packing proxy, interaction partners, confidence, and network role for the mutation site or another residue. |
| REQ-011 | Estimate mutant local compatibility | Notebook applies a conservative mutation model and reports lost compatible interactions, newly incompatible contacts, likely clashes, packing/cavity signals, and local constraint redistribution. |
| REQ-012 | Compare wild type and mutant | Outputs include WT and mutant interaction tables, delta summaries, network-change summaries, compatibility maps, and side-by-side visualizations. |
| REQ-013 | Generate layered mechanistic explanation | Each conclusion is labeled as observation, interpretation, or hypothesis; generated text is traceable to computed evidence. |
| REQ-014 | Rank structural disruption likelihood | Notebook produces low/moderate/high structural disruption likelihood only from structural interactions, qualitative energetic approximations, and network signals; threshold policy is visible. |
| REQ-015 | Propagate confidence | Outputs distinguish structure confidence, interaction confidence, and mechanistic confidence; low-confidence regions downweight or qualify conclusions. |
| REQ-016 | Warn for unsupported biological context | Notebook warns for low-confidence AF regions, missing residues, uncertain assemblies, disorder, missing ligands/cofactors/PTMs/metals, membrane context, flexible loops, and unsupported protein classes. |
| REQ-017 | Provide visualizations | Notebook renders structure overview, confidence overlay, interaction network or neighborhood view, WT/mutant interaction perturbation view, compatibility map, and severity/confidence summary. |
| REQ-018 | Export analysis artifacts | Notebook can export interaction tables, mutation comparison summaries, confidence notes, caveats, and static figures in notebook-friendly formats. |
| REQ-019 | Preserve RUO and non-clinical framing | All notebook text and generated summaries avoid clinical, pathogenicity, disease severity, exact ddG, and deterministic folding claims. |

## Non-Functional Requirements

| ID | Requirement | Acceptance Criteria |
|---|---|---|
| NFR-001 | Runtime target | Happy-path fixture runs top-to-bottom in local Jupyter or Colab CPU in under 2 minutes after dependencies are available; larger structures may be marked optional. |
| NFR-002 | Dependency budget | Uses standard scientific Python and lightweight visualization libraries; no GPU-only, MD, PyTorch, or heavyweight force-field engines in v1. |
| NFR-003 | Reproducibility | Notebook is restart-and-run-all compatible and uses pinned fixtures or cached files for validation. |
| NFR-004 | Transparency | Interaction thresholds, graph weights, mutation model assumptions, and confidence rules are visible in the notebook. |
| NFR-005 | Privacy and compliance | Notebook does not require human subject, clinical, protected health, or regulated variant datasets. |
| NFR-006 | Graceful degradation | External API failures, missing optional metadata, or unavailable visualization backends produce clear warnings or fallbacks without corrupting downstream variables. |

## Non-Goals

- Clinical diagnosis, clinical reporting, pathogenicity prediction, disease severity, or patient-specific interpretation.
- Rigorous thermodynamic free energies, validated ddG values, total protein stability scores, or exact energetic claims.
- Molecular dynamics, folding pathways, ensemble thermodynamics, exact mutant conformational prediction, or allosteric propagation prediction.
- V1 support for ligand energetics, cofactors, PTMs, nucleic acids, explicit solvent, metal coordination, or membrane energetics.
- Automated integration of conservation, genomics, clinical databases, or evolutionary scores into severity ranking.
- Proteome-scale batch processing or production platform features.

## Edge Cases And Failure Behavior

| Case | Required Behavior |
|---|---|
| Invalid AFDB or PDBe identifier | Stop before retrieval-dependent cells; show identifier and expected format. |
| Network or API failure | Show source URL, exception summary, and fixture/cache fallback instructions. |
| Structure has no requested chain | Stop mutation analysis; display available chains. |
| Mutation wild-type residue does not match structure | Stop mutation analysis; show observed residue and requested residue. |
| Residue is missing coordinates or key atoms | Skip unsupported local calculations for that residue and explain why. |
| Low pLDDT or weak experimental quality near mutation | Continue only with downgraded confidence and visible warning. |
| No interactions detected near mutation site | Report zero detected interactions as an observation and avoid inferring stability from absence alone. |
| Structure likely requires missing ligand/cofactor/PTM/metal/membrane context | Warn that structural interpretation may be incomplete and suppress strong severity claims. |
| Visualization backend unavailable | Still produce tables and text summaries; show visualization fallback message. |

## Open Questions

| ID | Question | Type | Owner | Resolution |
|---|---|---|---|---|
| Q-001 | Which exact happy-path, edge-case, and negative fixtures should be pinned for this notebook? | blocking | Fixture selection / domain reviewer | Needed before implementation validation. |
| Q-002 | Is BioPython acceptable as the v1 parser for mmCIF/PDB, or should a more mmCIF-specialized parser replace it later? | human-review | Implementation lead | Provisional v1 decision: use BioPython parser wrappers; revisit only if fixture parsing exposes numbering, alternate-location, or assembly failures. |
| Q-003 | Are the provisional geometric thresholds and residue sets scientifically acceptable? | human-review | Domain reviewer | Adopt the interaction configuration in `notebook-design.md` for v1; require domain review before scientific sign-off. |
| Q-004 | Is sidechain-substitution-only mutation modeling sufficient for v1? | human-review | Domain reviewer / implementation lead | Adopt sidechain substitution only, no global minimization, no conformational sampling, and 8 Angstrom local recomputation radius. |
| Q-005 | Are provisional relative interaction weights acceptable as heuristic constraint weights? | human-review | Domain reviewer | Adopt baseline weights and relative-weight formula in `notebook-design.md`; label as non-thermodynamic. |
| Q-006 | How should low/moderate/high structural disruption likelihood thresholds be calibrated? | human-review | Domain reviewer | Use provisional descriptive rubric only; do not treat severity as calibrated until fixture/domain review. |
| Q-007 | What exact PDBe resolution and residue-coverage cutoffs qualify as a suitable experimental structure? | human-review | Domain reviewer / implementation lead | Provisional v1 policy: prefer PDBe when resolution and residue coverage are acceptable and document caveats; exact cutoffs remain domain-review items. |
| Q-008 | Should local user-uploaded structures be included in v1? | non-blocking | Product owner | Keep optional; AFDB/PDBe are primary. |
| Q-009 | How detailed should generated natural-language explanations be? | assumption | Product owner / reviewer | Use template-based summaries from computed evidence in v1. |
| Q-010 | Which visualization backend should be used for 3D structure rendering? | assumption | Implementation lead | Use MolViewSpec/Mol* for v1 3D structure rendering with Mol* viewer URL fallback; use Plotly/NetworkX for network views and Matplotlib for static fallback figures. |
