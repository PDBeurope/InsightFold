# Mapping Bond-Type/Energy for AF and PDBe Models
# Spec-Driven Development Specification

---

# 1. Project Overview

## 1.1 Purpose

This project implements a notebook-based structural interaction analysis system for:

- AlphaFold (AF) protein models
- PDBe experimental structures

The system analyzes proteins as:

> structural constraint networks.

The notebook detects:
- chemically meaningful interactions,
- compatibility relationships,
- interaction-network topology,
- mutation-induced perturbations,
- and local structural destabilization signals.

The primary objective is not prediction of biological truth.

The system is explicitly designed to:
- explain structural chemistry,
- interpret interaction perturbations,
- visualize compatibility changes,
- and generate mechanistic structural hypotheses.

The notebook must remain:
- deterministic,
- explainable,
- computationally lightweight,
- and scientifically defensible.

The system must not:
- train models,
- perform molecular dynamics,
- estimate rigorous ΔΔG,
- or claim pathogenicity prediction.

---

# 2. Scientific Model

## 2.1 Structural Philosophy

The system models proteins as:

> networks of cooperative structural constraints.

Protein stability is treated as an emergent property of:
- residue compatibility,
- interaction density,
- packing quality,
- electrostatic balance,
- and local geometric coherence.

Mutations are interpreted as:

> perturbations to the structural compatibility network.

The notebook should analyze:
- interaction loss,
- clash introduction,
- packing disruption,
- cavity formation,
- centrality changes,
- and altered local constraint density.

---

## 2.2 Scientific Positioning

The notebook is:

- a structural interpretation engine,
- a mechanistic reasoning system,
- and an interaction perturbation analyzer.

The notebook is NOT:

- a pathogenicity classifier,
- a protein folding simulator,
- a thermodynamic engine,
- or a machine-learning prediction system.

---

## 2.3 Structural Confidence Philosophy

The notebook must distinguish:

| Confidence Layer | Meaning |
|---|---|
| Structural confidence | trust in coordinates |
| Interaction confidence | trust in interaction geometry |
| Mechanistic confidence | trust in interpretation |

This distinction must propagate through:
- outputs,
- visualizations,
- scoring,
- and user explanations.

---

# 3. Notebook Scope

## 3.1 Included Features

The notebook must support:

### Structure retrieval
- AlphaFold structures
- PDBe experimental structures

### Interaction analysis
- hydrogen bonds
- salt bridges
- hydrophobic contacts
- aromatic interactions
- disulfides
- steric clashes
- local packing density

### Interaction network analysis
- weighted residue graphs
- interaction centrality
- local density
- network fragmentation
- local redundancy

### Mutation analysis
- sidechain substitution
- local interaction recomputation
- compatibility perturbation analysis
- local structural interpretation

### Visualization
- 3D structure rendering
- interaction overlays
- perturbation maps
- network graphs
- confidence overlays

### Interpretation generation
- observation layer
- mechanistic interpretation layer
- structural hypothesis layer

---

## 3.2 Explicitly Out of Scope

The notebook must NOT implement:

- molecular dynamics
- free-energy simulation
- ΔΔG prediction
- pathogenicity prediction
- ensemble modeling
- membrane energetics
- ligand energetics
- allostery prediction
- exact mutant structure prediction
- machine-learning models
- training pipelines

---

# 4. Technology Stack

## 4.1 Runtime

Required:

- Python 3.11+
- Jupyter Notebook or JupyterLab

---

## 4.2 Core Dependencies

### Required

```python
biopython
numpy
pandas
scipy
networkx
plotly
molviewspec
requests
```

---

### Recommended

```python
mdtraj
freesasa
matplotlib
seaborn
```

---

## 4.3 Optional

```python
biopandas
ipywidgets
cytoscape
```

---

# 5. Notebook Architecture

The notebook must be divided into deterministic execution sections.

Recommended notebook section ordering:

```text
1. Configuration
2. Structure Retrieval
3. Structure Normalization
4. Interaction Detection
5. Interaction Weighting
6. Network Construction
7. Mutation Analysis
8. Perturbation Analysis
9. Mechanistic Interpretation
10. Visualization
11. Export
```

---

# 6. Configuration Layer

## 6.1 Global Configuration Object

The notebook must expose a single configuration dictionary.

Example:

```python
CONFIG = {
    "structure_source_priority": ["pdbe", "alphafold"],
    "mutation_radius": 8.0,
    "hbond_distance_cutoff": 3.5,
    "salt_bridge_cutoff": 4.0,
    "hydrophobic_cutoff": 5.0,
    "aromatic_cutoff": 6.0,
    "disulfide_cutoff": 2.2,
    "af_low_confidence_threshold": 50,
    "af_high_confidence_threshold": 80,
    "interaction_weights": {
        "disulfide": 1.0,
        "salt_bridge": 0.75,
        "aromatic": 0.65,
        "hydrogen_bond": 0.60,
        "hydrophobic": 0.50,
        "vdw_packing": 0.40,
        "steric_clash": -1.0
    }
}
```

---

# 7. Input Specification

# 7.1 Supported Inputs

## AlphaFold

Accepted:

- UniProt accession
- AFDB accession
- local PDB/mmCIF

Examples:

```text
P04637
AF-P04637-F1-model_v4
```

---

## PDBe

Accepted:

- PDB accession
- local PDB/mmCIF

Examples:

```text
1TUP
7XYZ
```

---

## Mutation Input

Supported formats:

```text
R273H
A:W128A
CHAIN_A:W128A
```

---

# 7.2 Validation Rules

Reject if:

- invalid structure ID
- invalid mutation syntax
- missing residue
- wild-type mismatch
- unsupported file format
- empty structure

Warn if:

- low AF confidence
- missing residues
- alternate conformations
- low experimental resolution
- incomplete assembly

---

# 8. Structure Retrieval Module

## 8.1 Responsibilities

The retrieval module must:

- fetch structures
- retrieve metadata
- identify provenance
- handle retries
- normalize source selection

---

## 8.2 AlphaFold Retrieval

Base endpoint:

```text
https://alphafold.ebi.ac.uk/files/
```

Retrieve:

```text
AF-{UNIPROT}-F1-model_v4.pdb
```

---

## 8.3 PDBe Retrieval

Base endpoint:

```text
https://files.rcsb.org/download/
```

Preferred format:

```text
mmCIF
```

---

## 8.4 Structure Selection Logic

Default policy:

1. Use PDBe if:
   - acceptable resolution
   - sufficient residue coverage

2. Otherwise use AF

Manual override required.

---

# 9. Structure Normalization

## 9.1 Required Normalization Steps

The notebook must:

1. parse structure
2. normalize chain IDs
3. normalize residue numbering
4. remove unsupported heteroatoms
5. resolve alternate locations
6. identify unresolved residues
7. annotate missing atoms
8. build biological assembly if available

---

## 9.2 Confidence Annotation

The notebook must compute:

### AF confidence
From:
- pLDDT

### Experimental confidence
From:
- resolution
- atom completeness
- residue occupancy

---

# 10. Interaction Detection Engine

# 10.1 General Principles

Interactions represent:

> relative structural compatibility constraints.

They do NOT represent:
- rigorous physical energies,
- or experimentally validated energetic truth.

---

# 10.2 Hydrogen Bonds

## Detection Criteria

Distance:

```text
≤ 3.5 Å
```

Angle:

```text
≥ 120°
```

## Output Fields

```python
{
    "type": "hydrogen_bond",
    "residue_1": "A:45",
    "residue_2": "A:80",
    "distance": 2.9,
    "angle": 145,
    "geometry_confidence": 0.91
}
```

---

# 10.3 Salt Bridges

## Detection Criteria

Opposite charge groups.

Distance:

```text
≤ 4.0 Å
```

---

# 10.4 Hydrophobic Contacts

## Detection Criteria

Hydrophobic sidechain centroid distance:

```text
≤ 5.0 Å
```

Hydrophobic residues:

```text
A, V, I, L, M, F, W, Y
```

---

# 10.5 Aromatic Interactions

## Detection Criteria

Ring centroid distance:

```text
≤ 6.0 Å
```

Aromatic residues:

```text
F, W, Y, H
```

---

# 10.6 Disulfides

## Detection Criteria

Sulfur-sulfur distance:

```text
≤ 2.2 Å
```

---

# 10.7 Steric Clashes

## Detection Criteria

van der Waals overlap threshold exceeded.

Clashes should generate:
- negative compatibility weighting.

---

# 10.8 Packing Density

Packing density should estimate:

- local neighbor density
- local buriedness
- cavity introduction

Optional:

- SASA integration using FreeSASA

---

# 11. Interaction Weighting Model

# 11.1 Weighting Philosophy

Interaction weights represent:

> relative contributions to structural constraint stability.

They do NOT represent:

- kcal/mol,
- free energies,
- thermodynamic measurements.

---

# 11.2 Baseline Interaction Weights

| Interaction | Weight |
|---|---|
| disulfide | 1.00 |
| salt_bridge | 0.75 |
| aromatic | 0.65 |
| hydrogen_bond | 0.60 |
| hydrophobic | 0.50 |
| vdw_packing | 0.40 |
| steric_clash | -1.00 |

---

# 11.3 Relative Weight Formula

```text
relative_weight =
base_interaction_weight
× geometry_quality
× burial_factor
× confidence_factor
× cooperativity_factor
```

---

# 11.4 Factor Definitions

## Geometry Quality

Derived from:

- distance quality
- angle quality
- steric validity

Range:

```text
0.0 → 1.0
```

---

## Burial Factor

Approximates:

- solvent exclusion
- packing depth

Buried interactions should contribute more strongly.

---

## Confidence Factor

Derived from:

- AF pLDDT
- experimental resolution
- atom completeness

---

## Cooperativity Factor

Approximates:

- local interaction density
- network reinforcement
- cooperative stabilization

---

# 12. Network Construction

# 12.1 Graph Model

## Nodes

Residues.

---

## Edges

Interactions.

Multi-edge support required.

Reason:

Two residues may share multiple interaction classes.

---

# 12.2 Edge Schema

```python
{
    "source": "A:W128",
    "target": "A:F91",
    "interaction_type": "aromatic",
    "geometry_confidence": 0.91,
    "burial_factor": 0.82,
    "cooperativity_factor": 0.74,
    "relative_weight": 0.63
}
```

---

# 12.3 Network Metrics

Required metrics:

- weighted degree centrality
- betweenness centrality
- interaction density
- clustering coefficient
- local redundancy
- local fragmentation

Avoid:

- opaque graph metrics without interpretation value.

---

# 13. Mutation Engine

# 13.1 Mutation Philosophy

The notebook must:

- analyze compatibility perturbation,
- not predict exact mutant structures.

---

# 13.2 Mutation Strategy

V1 implementation:

- sidechain substitution only
- no global minimization
- no conformational sampling

---

# 13.3 Mutation Radius

Default:

```text
8 Å
```

Used for:

- local interaction recalculation
- compatibility analysis
- perturbation scoring

---

# 13.4 Mutation Analysis Outputs

Required outputs:

- lost interactions
- gained interactions
- clash introduction
- cavity formation
- packing reduction
- centrality changes
- compatibility shifts

---

# 14. Perturbation Interpretation Engine

# 14.1 Interpretation Layers

The notebook must generate:

## Observation Layer

Directly computed facts.

Examples:

- hydrogen bond lost
- clash introduced
- aromatic contact removed

---

## Interpretation Layer

Mechanistic structural reasoning.

Examples:

- hydrophobic packing weakened
- electrostatic compensation reduced
- β-sheet compatibility disrupted

---

## Hypothesis Layer

Higher-level speculative consequence.

Examples:

- local flexibility may increase
- local destabilization possible
- local backbone adaptation may be required

---

# 14.2 Language Constraints

The notebook must use:

- “consistent with”
- “may weaken”
- “suggests”
- “possible destabilization”

The notebook must avoid:

- “causes unfolding”
- “guarantees instability”
- “abolishes function”

---

# 15. Severity Ranking

# 15.1 Severity Philosophy

The notebook ranks:

> structural disruption likelihood.

The notebook does NOT rank:

- disease severity
- pathogenicity
- biological outcome

---

# 15.2 Severity Contributors

Severity should derive from:

| Signal | Interpretation |
|---|---|
| clash introduction | strong destabilization |
| buried polarity mismatch | strong concern |
| interaction loss | destabilization |
| cavity formation | packing reduction |
| centrality loss | structural importance |
| network fragmentation | cooperative destabilization |
| secondary structure incompatibility | moderate-high concern |

---

# 15.3 Severity Classes

```text
LOW
MODERATE
HIGH
```

These classes are:

- heuristic
- structure-centric
- non-clinical

---

# 16. Visualization Requirements

# 16.1 3D Structure Viewer

Recommended:

```python
molviewspec
```

Required features:

- residue selection
- chain filtering
- mutation highlighting
- confidence coloring
- interaction overlays
- `.mvsj` state generation when export is requested
- inline Mol* rendering or Mol* viewer URL fallback

MolViewSpec/Mol* rules:

- Keep MolViewSpec visualization optional; no data-loading, scoring, validation, or interpretation cell may depend on MolViewSpec-only variables.
- Verify the installed MolViewSpec API before generating views.
- Prefer AFDB `bcifUrl` when available; otherwise use mmCIF.
- Use parsed chain labels and residue numbering.
- Skip individual views with clear messages when required inputs are missing.

---

# 16.2 Interaction Overlay

Required interaction colors:

| Interaction | Color |
|---|---|
| hydrogen bond | blue |
| salt bridge | red |
| hydrophobic | gold |
| aromatic | purple |
| disulfide | green |
| steric clash | orange |

---

# 16.3 Perturbation Visualization

The notebook must visualize:

- interaction loss
- interaction gain
- clash zones
- local packing changes
- compatibility perturbation

---

# 16.4 Network Graph

Framework:

```python
plotly
networkx
```

Required features:

- node centrality scaling
- interaction filtering
- local neighborhood highlighting
- weighted edge rendering

---

# 16.5 Confidence Overlay

Must display:

- AF pLDDT confidence
- low-confidence dimming
- interpretation confidence labels

---

# 17. Notebook Workflow

# Step 1 — Structure Input

User enters:

- AF accession
- PDB accession
- local file

Notebook validates input.

---

# Step 2 — Structure Retrieval

Notebook:

- retrieves structure
- retrieves metadata
- identifies provenance

---

# Step 3 — Structure Normalization

Notebook:

- normalizes structure
- annotates missing regions
- resolves chain numbering

---

# Step 4 — Interaction Analysis

Notebook computes:

- interactions
- weighted interaction graph
- local packing density
- centrality metrics

---

# Step 5 — Structural Visualization

Notebook renders:

- structure
- interactions
- confidence overlays

---

# Step 6 — Mutation Submission

User enters mutation.

Notebook validates:

- residue existence
- WT compatibility
- chain validity

---

# Step 7 — Perturbation Analysis

Notebook recomputes:

- local interaction environment
- compatibility shifts
- centrality changes
- structural perturbations

---

# Step 8 — Mechanistic Interpretation

Notebook generates:

- observations
- interpretations
- structural hypotheses

---

# Step 9 — Severity Assessment

Notebook outputs:

- LOW / MODERATE / HIGH
- confidence explanation
- uncertainty notes

---

# Step 10 — Export

Notebook exports:

- JSON
- CSV
- PNG/SVG figures
- graph objects

---

# 18. Output Schema

# 18.1 Interaction Record

```python
{
    "residue_1": "A:W128",
    "residue_2": "A:F91",
    "interaction_type": "aromatic",
    "distance": 5.2,
    "geometry_confidence": 0.92,
    "relative_weight": 0.63,
    "burial_factor": 0.81,
    "confidence": "high"
}
```

---

# 18.2 Mutation Perturbation Record

```python
{
    "mutation": "W128A",
    "lost_interactions": 4,
    "gained_clashes": 1,
    "packing_change": "reduced",
    "centrality_delta": -0.42,
    "severity": "HIGH"
}
```

---

# 18.3 Interpretation Record

```python
{
    "observation": "Two aromatic interactions lost.",
    "interpretation": "Hydrophobic packing network weakened.",
    "hypothesis": "Local flexibility may increase.",
    "confidence": "moderate"
}
```

---

# 19. Error Handling

# 19.1 Retrieval Errors

Examples:

```text
Unable to retrieve structure from PDBe.
Attempted fallback to AlphaFold.
```

---

# 19.2 Validation Errors

Examples:

```text
Mutation residue does not match wild-type structure residue.
```

---

# 19.3 Confidence Warnings

Examples:

```text
Mutation occurs in low-confidence AlphaFold region.
Interpretation reliability reduced.
```

---

# 20. User-Facing Explanations

# 20.1 Onboarding

> This notebook analyzes proteins as interaction-based structural constraint networks. Mutations are interpreted as perturbations to structural compatibility rather than deterministic biological outcomes.

---

# 20.2 AlphaFold Warning

> AlphaFold models are structural predictions. Interactions involving low-confidence regions should be interpreted cautiously.

---

# 20.3 Structural Interpretation Warning

> This notebook estimates structural compatibility and interaction perturbation. It does not predict pathogenicity, folding thermodynamics, or exact mutant conformations.

---

# 20.4 Flexibility Warning

> Flexibility-related interpretations are inferred from weakened constraints and reduced packing, not directly measured dynamics.

---

# 21. Logging and Observability

# 21.1 Required Logging

Log:

- structure retrieval failures
- parsing failures
- mutation validation failures
- interaction engine warnings
- confidence downgrades
- export failures

---

# 21.2 Intermediate Artifact Storage

Retain:

- normalized structures
- interaction graphs
- perturbation graphs
- severity breakdowns
- interpretation traces

---

# 21.3 Telemetry

Track:

- structure source usage
- interaction-type frequency
- mutation analysis frequency
- confidence-warning frequency

---

# 22. Performance Requirements

# 22.1 Runtime Targets

## Structure analysis

Target:

```text
< 30 seconds
```

for average monomer.

---

## Mutation analysis

Target:

```text
< 10 seconds
```

for local perturbation recomputation.

---

# 22.2 Computational Constraints

The notebook must:

- run on CPU
- avoid GPU dependency
- avoid large-scale simulation
- remain interactive in Jupyter

---

# 23. Deliverables

The notebook implementation must provide:

## Required

- executable Jupyter notebook
- modular Python functions
- interaction graph generation
- mutation perturbation analysis
- visualization components
- JSON export support
- CSV export support

---

## Recommended

- reusable Python package structure
- notebook examples
- sample datasets
- unit tests for interaction detection

---

# 24. Future Expansion (Post-V1)

Potential future additions:

- local minimization
- interaction robustness across structures
- ensemble-aware interpretation
- ligand-aware interactions
- multimer-aware interface scoring
- mutation batch processing
- web application deployment

These are explicitly excluded from V1.
