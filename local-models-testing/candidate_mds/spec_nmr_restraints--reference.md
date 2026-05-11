# Project Specification

# NMR Restraints Visualization Notebook

---

## Overview

### Purpose

The NMR restraints visualization notebook is an interactive exploratory analysis system for understanding how experimentally derived NMR restraints relate to deposited macromolecular structures.

The system addresses the question:

> “Which regions of an NMR structure are experimentally supported, weakly constrained, or locally inconsistent with the deposited NMR restraints?”

The notebook transforms raw restraint and coordinate data into an interpretable spatial evidence model that allows users to explore experimental support directly within a 3D structural context.

---

### Intended Users

Primary users include:

- structural biologists
- protein biochemists
- computational biology researchers
- Protein Data Bank users
- researchers familiar with protein structure interpretation but not necessarily expert in NMR refinement

The notebook is intentionally designed for users who:
- understand protein structures,
- understand residue-level interpretation,
- but may not understand NMR restraint semantics or refinement workflows.

---

### Product Description

The notebook accepts a Protein Data Bank (PDB) identifier and generates:

- a reproducible exploratory analysis report,
- interactive 3D structural visualizations,
- restraint density maps,
- restraint violation analyses,
- residue-level navigation tools,
- and local evidence inspection interfaces.

The system integrates:
- structure coordinates,
- NMR restraint data,
- geometric analysis,
- and interactive visualization.

Visualization is state-driven and reproducible using MolViewSpec and Mol*.

---

### Primary Outcome

The notebook enables users to:

- identify experimentally well-supported structural regions,
- identify weakly restrained regions,
- inspect violated restraints spatially,
- understand whether violations are isolated or clustered,
- explore local experimental evidence networks,
- and interpret deposited NMR evidence without specialist tooling.

The notebook prioritizes:
- evidence transparency,
- contextual interpretation,
- and exploratory reasoning.

The system does not produce:
- automated quality scores,
- validation verdicts,
- or global confidence metrics.

---

### Ecosystem Context

The notebook fits within:
- structural biology analysis workflows,
- Protein Data Bank exploratory tooling,
- NMR validation and interpretation workflows,
- and reproducible notebook-based scientific reporting.

The project complements:
- geometric validation systems,
- refinement software,
- ensemble analysis packages,
- and specialist NMR interpretation tools.

Its strategic differentiation is that it visualizes experimental support directly instead of reducing evidence to aggregate validation metrics.

---

# Domain Context

## Scientific Context

NMR-derived structures are inferred structural models constrained by experimentally measured restraints.

Unlike X-ray crystallography or cryo-EM:
- the final structure is not directly observed,
- but inferred through optimization against experimental constraints.

The deposited structure therefore represents:
- one possible model,
- or one representative conformer,
- satisfying a collection of restraints.

Understanding the relationship between:
- restraints,
- violations,
- and local geometry
is essential for meaningful interpretation.

---

## Key Domain Concepts

### Distance Restraints

Distance restraints specify expected spatial proximity between atoms.

Representation:

| Field | Meaning |
|---|---|
| atom_1 | first restrained atom |
| atom_2 | second restrained atom |
| lower_bound | minimum allowed distance |
| upper_bound | maximum allowed distance |

Distance restraints commonly arise from NOE-derived experiments.

---

### Dihedral Restraints

Dihedral restraints constrain torsion angles defined by four atoms.

Representation:

| Field | Meaning |
|---|---|
| atom_1–4 | torsion-defining atoms |
| lower_bound | minimum allowed angle |
| upper_bound | maximum allowed angle |

---

### Violation

A restraint violation occurs when the measured geometric value lies outside the allowed range.

Violation magnitude is defined as:

```text
distance_to_nearest_allowed_boundary
```

Violations indicate:
- local inconsistency,
- competing restraints,
- model limitations,
- or experimental ambiguity.

Violations are not automatically evidence of incorrect structures.

---

### Restraint Density

Restraint density is defined as:

```text
number_of_restraints_associated_with_residue
```

Density is interpreted as:
- a proxy for experimental support,
- not a confidence or correctness score.

All restraints associated with a residue contribute equally:
- backbone restraints,
- sidechain restraints,
- distance restraints,
- dihedral restraints,
- and supported ambiguous restraints.

---

## Interpretation Risks

The notebook must avoid several misleading interpretations.

### Misconception: “No violations means the structure is correct”

False because:
- restraints may be sparse,
- refinement may overfit,
- or incompatible restraints may have been omitted.

---

### Misconception: “High restraint density equals certainty”

False because:
- restraints vary in informational value,
- restraint redundancy may occur,
- and refinement artifacts may exist.

---

### Misconception: “Violations are errors”

False because:
- violations may reflect:
  - conformational heterogeneity,
  - local dynamics,
  - incompatible experimental constraints,
  - or refinement tradeoffs.

The notebook must present violations as:
- evidence of inconsistency,
- not automatic failure.

---

# Input Specification

## Primary User Input

### Required Input

| Field | Type | Required | Description |
|---|---|---|---|
| pdb_id | string | yes | Four-character PDB accession |

Validation regex:

```regex
^[A-Za-z0-9]{4}$
```

Input must be normalized to uppercase before retrieval.

---

## Optional Runtime Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| model_index | integer | 0 | Conformer index to analyze |
| violation_threshold_distance | float | 0.5 | Minimum distance violation to highlight |
| violation_threshold_dihedral | float | 5.0 | Minimum angular violation to highlight |
| local_context_radius | float | 6.0 Å | Spatial neighborhood radius |
| max_visible_restraints | integer | 250 | Rendering clutter limit |
| cache_enabled | boolean | true | Enable local cache |

---

## External Data Sources

### Structure Coordinates

Primary source:
- RCSB PDB

Preferred format:
- Binary CIF (`bcif`)

Fallback:
- mmCIF

Endpoints:

```text
https://models.rcsb.org/{PDB_ID}.bcif
https://files.rcsb.org/download/{PDB_ID}.cif
```

---

### NMR Restraints

Source:
- BMRB/PDB NMR-STAR deposition files

Supported categories:

| Restraint Type | Category |
|---|---|
| Distance restraints | `_Gen_dist_constraint_list` |
| Dihedral restraints | `_Torsion_angle_constraint` |

The parser must support:
- multiple restraint loops,
- multiple saveframes,
- and multiple tables within a single file.

All compatible tables are merged into a unified restraint collection.

---

## Ambiguous Distance Restraints

### Detection Rule

A distance restraint is ambiguous when:

| Column | Value |
|---|---|
| `Member_logic_code` | `OR` |

---

### Evaluation Semantics

For ambiguous restraints:

```text
evaluate all candidate atom-pair distances
→ select smallest distance
→ compare against allowed bounds
```

The smallest distance:
- determines violation status,
- is the only member visualized,
- and is the only member reported in tables.

Alternative members are not visualized in v1.

---

## Atom Mapping Strategy

The implementation assumes deposited restraints are already annotation-compatible with deposited structure atoms.

No external atom mapping library is required.

Canonical atom identity:

```text
(
    auth_asym_id,
    auth_seq_id,
    auth_comp_id,
    auth_atom_id
)
```

---

### Distance Restraint Mapping

| Restraint Field | Model Field |
|---|---|
| `_Gen_dist_constraint.Auth_asym_ID_1` | `_atom_site.auth_asym_id` |
| `_Gen_dist_constraint.Auth_seq_ID_1` | `_atom_site.auth_seq_id` |
| `_Gen_dist_constraint.Auth_comp_ID_1` | `_atom_site.auth_comp_id` |
| `_Gen_dist_constraint.Auth_atom_ID_1` | `_atom_site.auth_atom_id` |

Equivalent mapping applies for `_2` fields.

---

### Dihedral Restraint Mapping

For atoms 1–4:

| Restraint Field | Model Field |
|---|---|
| `_Torsion_angle_constraint.Auth_asym_ID_n` | `_atom_site.auth_asym_id` |
| `_Torsion_angle_constraint.Auth_seq_ID_n` | `_atom_site.auth_seq_id` |
| `_Torsion_angle_constraint.Auth_comp_ID_n` | `_atom_site.auth_comp_id` |
| `_Torsion_angle_constraint.Auth_atom_ID_n` | `_atom_site.auth_atom_id` |

Matching must be exact.

---

## Mapping Coverage Thresholds

| Mapping Coverage | Behaviour |
|---|---|
| >95% | normal execution |
| 70–95% | continue with warning |
| <70% | abort analysis |

Coverage formula:

```text
mapping_coverage =
mapped_restraints / total_parsed_restraints
```

---

## Alternate Conformations

Alternate conformations (`altloc`) are ignored in v1.

Only:
- primary/default conformers,
- blank altlocs,
- or `.` altloc records
are used.

---

## Retrieval Behaviour

### Retry Logic

For network failures:
- retry up to 3 times,
- exponential backoff:
  - 1 s
  - 2 s
  - 4 s

---

## Failure Conditions

Hard failures:
- structure unavailable,
- restraint file unavailable,
- malformed coordinate file,
- no parseable restraints,
- mapping coverage below threshold.

---

## Warning Conditions

Warnings:
- partially unmapped restraints,
- unsupported ambiguity semantics,
- missing restraint categories,
- ignored alternate conformations,
- unsupported records.

---

# System Execution / Processing Logic

## Execution Flow

```text
Input PDB ID
→ Retrieve structure
→ Retrieve restraints
→ Parse coordinates
→ Parse restraint tables
→ Normalize atom identifiers
→ Map restraint atoms
→ Compute geometry
→ Compute violations
→ Compute residue density
→ Build MolViewSpec states
→ Generate interactive notebook
```

---

## Runtime Environment

### Recommended Environment

| Component | Requirement |
|---|---|
| Python | 3.11+ |
| Notebook | JupyterLab |
| OS | Linux/macOS/Windows |

---

## Core Dependencies

### Scientific Stack

| Library | Purpose |
|---|---|
| numpy | geometry computation |
| pandas | tabular processing |
| biopython | structure parsing |
| pynmrstar | NMR-STAR parsing |
| scipy | optional geometry utilities |

---

### Visualization Stack

| Library | Purpose |
|---|---|
| molviewspec | declarative visualization |
| Mol* | rendering |
| ipywidgets | controls |
| plotly | optional plots |

Canonical MolViewSpec version:

```text
1.8.1
```

---

## MolViewSpec Integration Strategy

### Visualization Philosophy

Mol* is treated as:
- a rendering engine,
- not the primary interaction controller.

MolViewSpec is the canonical:
- scene-definition layer,
- reproducibility layer,
- and visualization abstraction layer.

---

### Canonical Builder Pattern

```python
builder = mvs.create_builder()

structure = (
    builder
    .download(url=structure_url)
    .parse(format=structure_format)
    .model_structure()
)
```

---

### Supported Parse Formats

| Format | `parse(format=...)` value |
|---|---|
| Binary CIF | `'bcif'` |
| mmCIF | `'mmcif'` |

`bcif` is preferred when available.

---

### State Ownership

Notebook Python state is authoritative.

Mol* frontend state is ephemeral.

Selections originating from:
- sequence view,
- tables,
- filters,
- or notebook controls
must regenerate deterministic visualization states.

---

### Residue Selection Expressions

Canonical fields:

| Field | Purpose |
|---|---|
| `label_asym_id` | chain selection |
| `beg_label_seq_id` | residue range start |
| `end_label_seq_id` | residue range end |

---

### Per-Residue Coloring Strategy

Preferred:
- annotation-file-based coloring if supported in MolViewSpec 1.8.1.

Fallback:
- per-residue `ComponentExpression` loops.

---

### Large-State Guard

| Residues | Strategy |
|---|---|
| ≤ 500 | per-residue coloring |
| > 500 | grouped contiguous ranges |

---

### Viewer Wrapper

Canonical rendering helper:

```python
show_mol_view(state, label)
```

Preferred rendering path:

```python
state.molstar_html()
```

Fallback:
- serialized `.mvsj`
- external Mol* viewer URL.

---

### Failure Isolation

Each visualization panel must be isolated:

```python
try:
    ...
except Exception as e:
    ...
```

Failure in one panel must not prevent other views from rendering.

---

## Geometry Algorithms

### Distance Computation

Euclidean distance in Ångström units.

---

### Dihedral Computation

Torsion angles:
- computed from four atoms,
- normalized to `[-180°, +180°]`.

Allowed intervals must follow restraint-file semantics directly, including wrapped intervals.

Example:

```text
170° → -170°
```

represents a wrapped interval around ±180°.

---

### Distance Violation Logic

```python
if measured < lower:
    violation = lower - measured
elif measured > upper:
    violation = measured - upper
else:
    violation = 0
```

---

### Dihedral Violation Logic

Violation magnitude:
- computed as minimum angular distance to allowed interval,
- using circular angular geometry.

---

## Residue Density Computation

### Absolute Density

```text
absolute_density(residue) =
count(unique_restraints_associated_with_residue)
```

---

### Normalized Density

Structure-wide mean density:

```text
mean_density =
N_restraints / N_residues
```

Normalized density:

```text
normalized_density(residue) =
absolute_density(residue) / mean_density
```

---

## Caching

### Cached Artifacts

Cache:
- downloaded files,
- parsed structures,
- parsed restraints,
- computed violation tables.

---

### Cache Key

```text
{pdb_id}_{model_index}_{config_hash}
```

---

## Performance Targets

Target structure size:
- up to ~500 residues.

Target restraint count:
- up to ~10,000 restraints.

Expected runtime:
- <30 seconds on commodity hardware.

---

# Output Parsing and Interpretation

## Core Outputs

### Residue Density Table

| Field | Description |
|---|---|
| chain_id | chain identifier |
| residue_number | residue position |
| residue_name | residue type |
| absolute_density | raw restraint count |
| normalized_density | relative density |
| structure_mean_density | baseline mean |

---

### Violation Table

| Field | Description |
|---|---|
| restraint_id | unique identifier |
| restraint_type | distance/dihedral |
| ambiguity_type | NONE/OR |
| selected_member | evaluated atom pair |
| involved_residues | participating residues |
| measured_value | measured geometry |
| lower_bound | allowed minimum |
| upper_bound | allowed maximum |
| violation_magnitude | deviation magnitude |

---

## Violation Prioritization

Default sort order:
1. descending violation magnitude
2. restraint type
3. residue number

---

## Empty Result Behaviour

### No Violations

Display:

> No violations exceeding current thresholds were detected.

Additional explanatory text must clarify that this:
- does not imply correctness,
- and is not a validation verdict.

---

### Missing Restraints

Display:

> No compatible restraint data were available for this structure.

Notebook exits gracefully.

---

# Visualisations / UI Components

# 1. Global Structure View

## Purpose

Provide overview of experimental support distribution.

---

## Rendering

Representation:
- cartoon backbone.

Residues colored by:
- restraint density.

Default mode:
- normalized density.

Users may toggle:
- normalized density,
- absolute density.

---

## Global Restraint Rendering Policy

Only violated restraints are rendered globally.

Satisfied restraints are hidden.

If violated restraint count exceeds limit:
- render highest violations first.

---

## Interaction

Supported:
- rotate,
- zoom,
- hover,
- click.

Hover tooltip:
- residue identifier,
- absolute density,
- normalized density,
- local violation count.

---

# 2. Local Evidence View

## Purpose

Inspect local experimental evidence around a selected residue.

---

## Neighborhood Definition

Neighborhood is spatially defined.

A residue is included if:

```text
distance(
    selected_residue.CA,
    candidate_residue.CA
) ≤ local_radius
```

Default radius:

```text
6.0 Å
```

Radius must be user-configurable.

---

## Local Rendering Policy

Render:
- all restraints involving selected residue.

This includes:
- satisfied restraints,
- violated restraints,
- ambiguous restraints.

---

## Transparency Rules

| Element | Rendering |
|---|---|
| violated restraints | opaque/high saliency |
| satisfied restraints | semi-transparent |
| non-violated neighboring residues | reduced opacity |

---

# 3. Distance Restraint Overlay

## Purpose

Visualize pairwise spatial constraints.

---

## Encoding

| Violation Magnitude | Color |
|---|---|
| none | translucent gray |
| low | yellow |
| medium | orange |
| high | red |

---

# 4. Dihedral Restraint View

## Purpose

Visualize torsional constraints separately from distance restraints.

---

## Rendering

Highlight:
- restrained bond,
- local backbone region.

Distance and dihedral views are separated by default to reduce clutter.

---

# 5. Sequence View

## Purpose

Provide linear residue navigation.

---

## Behaviour

Selecting residues:
- updates local structure view,
- filters violation table.

Hover state synchronized with 3D view.

---

## Accessibility

Sequence visualization must not rely exclusively on color.

Additional encodings:
- saturation,
- borders,
- tooltips.

---

# 6. Violation Table

## Purpose

Provide sortable/filterable violation inspection.

---

## Behaviour

Capabilities:
- sort,
- search,
- filter,
- select.

Selecting a row:
- centers local visualization,
- highlights involved atoms,
- updates local evidence view.

---

# Workflow / User Journey

## Stage 1 — Notebook Launch

User provides:
- PDB ID,
- optional runtime parameters.

System validates inputs and begins retrieval.

---

## Stage 2 — Data Retrieval

System retrieves:
- structure coordinates,
- restraint files.

Progress indicators displayed.

Failures terminate gracefully with diagnostics.

---

## Stage 3 — Parsing and Mapping

System:
- parses coordinates,
- parses restraints,
- maps atoms.

Notebook reports:
- parsed restraints,
- mapping coverage,
- ignored records,
- excluded restraints.

---

## Stage 4 — Geometry Analysis

System computes:
- distances,
- torsions,
- violations,
- residue densities.

---

## Stage 5 — Global Exploration

Notebook displays:
- global structure,
- density mapping,
- violation overlays,
- synchronized tables.

---

## Stage 6 — Local Exploration

User selects:
- residues,
- sequence regions,
- or violation rows.

Notebook updates:
- local neighborhood,
- restraint overlays,
- local evidence context.

---

## Stage 7 — Interpretation

Users identify:
- dense support regions,
- sparse regions,
- isolated inconsistencies,
- clustered violations.

---

# Component / Module Breakdown

# 1. Retrieval Module

## Responsibility

Download:
- structure files,
- restraint files.

---

# 2. Structure Parser

## Responsibility

Parse coordinate structures into normalized atom tables.

---

# 3. Restraint Parser

## Responsibility

Extract:
- distance restraints,
- dihedral restraints,
- ambiguity metadata.

---

# 4. Atom Mapping Module

## Responsibility

Resolve restraint atoms against coordinate atoms using exact `auth_*` matching.

---

# 5. Geometry Engine

## Responsibility

Compute:
- distances,
- torsions,
- violations.

---

# 6. Density Engine

## Responsibility

Compute:
- absolute density,
- normalized density.

---

# 7. Visualization State Builder

## Responsibility

Generate deterministic MolViewSpec states.

---

# 8. Notebook UI Layer

## Responsibility

Coordinate:
- widgets,
- state synchronization,
- rendering,
- filtering,
- user interactions.

---

# User-Facing Explanations and Copy

## Notebook Introduction

> This notebook visualizes how deposited NMR restraints relate to the structural model of a protein.
>
> The goal is to help users explore where experimental restraint support is dense, sparse, or locally inconsistent.
>
> Restraint density is shown as a proxy for experimental support, not as a confidence or accuracy score.

---

## Density Tooltip

> Restraint density reflects the number of deposited restraints associated with a residue, including both backbone and sidechain restraints.
>
> Density is intended as a proxy for experimental coverage and should not be interpreted as a confidence or accuracy score.

---

## Normalized Density Tooltip

> Normalized density reflects restraint coverage relative to the structure-wide average restraint density.

---

## Violation Tooltip

> A restraint violation occurs when the measured geometry falls outside the deposited allowed range.
>
> Violations may reflect experimental ambiguity, conformational heterogeneity, refinement tradeoffs, or local inconsistency.

---

## No Violations Message

> No violations exceeding the current visualization thresholds were detected.
>
> This does not imply that the structure is fully validated or error-free.

---

## Mapping Warning

> Some restraint atoms could not be mapped onto the deposited structure coordinates.
>
> Results may represent only a subset of deposited restraints.

---

## Unsupported Ambiguity Warning

> Some ambiguous restraint semantics are not supported in v1 and were excluded from analysis.

---

## Empty Restraint State

> No compatible restraint data were available for this structure.

---

## Trust Guidance

> This notebook is intended as an exploratory evidence viewer.
>
> It does not produce automated quality scores or validation verdicts.

---

# Data, Logging, and Observability

## Stored Data

Optional persistence:
- downloaded files,
- parsed structures,
- computed tables,
- visualization states.

---

## Logging Requirements

Required logs:
- retrieval success/failure,
- parser diagnostics,
- mapping statistics,
- excluded restraint counts,
- runtime duration.

---

## Metrics

Collect:
- runtime duration,
- restraint counts,
- violation counts,
- rendering time.

---

## Debug Diagnostics

Diagnostic outputs:
- parse statistics,
- mapping coverage,
- ignored records,
- geometry failures.

---

## Security Considerations

The notebook must:
- avoid arbitrary code execution from downloaded files,
- use safe parsing practices,
- avoid executing remote content.

No authentication required for public PDB retrieval.

---

# Integration Points

## Upstream Dependencies

- RCSB structure retrieval
- BMRB/PDB restraint retrieval
- MolViewSpec 1.8.1
- Mol* rendering framework

---

## Export Formats

Supported exports:
- CSV violation tables,
- notebook HTML export.

Formal visualization-state export schemas are deferred from v1.

---

## Reusability

Outputs should remain reusable for:
- future web applications,
- batch processing pipelines,
- automated reporting workflows.

---

# Non-Functional Requirements

## Performance

Target runtime:
- <30 seconds for typical structures.

---

## Reliability

Notebook must fail gracefully with:
- partial results where possible,
- clear diagnostics,
- isolated visualization failures.

---

## Maintainability

Implementation must separate:
- parsing,
- geometry,
- visualization,
- and UI orchestration.

---

## Portability

Notebook should run:
- locally,
- in Binder/JupyterHub,
- in cloud notebook environments.

---

## Accessibility

Visual encodings must not rely exclusively on color.

Tables and controls should remain keyboard accessible where possible.

---

## Usability

Default visualization state must prioritize:
- low clutter,
- readability,
- orientability.

---

# Out of Scope for v1

| Feature | Rationale |
|---|---|
| Ensemble-wide analysis | Simplifies interpretation |
| Weighted ambiguity semantics | Excessive complexity |
| Pseudoatom handling | Complex mapping semantics |
| Dynamic averaging semantics | Beyond exploratory scope |
| Automated quality scoring | Avoid overclaiming certainty |
| Machine-learning interpretation | Not required for MVP |
| Refinement recalculation | Outside visualization scope |
| Cross-structure comparison | Deferred feature |
| Collaborative state sharing | Not required for notebook MVP |
| Formal visualization-state interchange schema | Deferred |
| Custom Mol* camera choreography | Default Mol* behaviour sufficient |

---

# Final v1 Scientific Constraints

The notebook intentionally:
- emphasizes evidence transparency,
- avoids global confidence scores,
- avoids binary validation labels,
- and exposes underlying restraint evidence directly.

The system is an exploratory experimental evidence viewer,
not a formal validation engine.
