# NMR Restraint Evidence Viewer – V1 Implementation Specification

## 1. Introduction and Purpose

This document defines a complete, implementation-ready specification for the NMR Restraint Evidence Viewer. The system produces an interactive exploratory report that maps experimental NMR restraints (distance and dihedral) onto a 3D protein structure, enabling users to visually distinguish strongly supported, weakly constrained, and inconsistently restrained regions.

The target audience is an engineer building the system. The specification is:
- **Implementation-ready** – all data structures, algorithms, and UI behaviors are defined.
- **Technically precise** – fields, types, formulas, and tolerances are explicit.
- **Operationally clear** – user interaction flows and system states are fully described.
- **Internally consistent** – every component’s output feeds deterministically into the next.
- **Detailed enough to build from without ambiguity.**

---

## 2. System Overview

The system runs as a reproducible computational notebook (e.g., Jupyter). A user supplies a PDB ID. The pipeline:

```
PDB ID → Data retrieval → Restraint parsing → Atom mapping →
Geometric evaluation → Violation computation → Interactive visualization
```

The output is an interactive dashboard containing:
- A 3D molecular viewer (Mol* driven by MolViewSpec)
- A sequence-based residue navigator
- A sortable, filterable violation table
- Control toggles (global/local, distance/dihedral modes)

All visualisation is **state-driven**; Mol* acts as a rendering surface, not the interaction controller.

---

## 3. Data Models

### 3.1 Canonical Atom Identifier

Used in both restraints and model atoms after mapping.

| Field            | Type   | Description                          | Example    |
|------------------|--------|--------------------------------------|------------|
| `chain_id`       | string | Chain identifier (auth_asym_id)      | "A"        |
| `residue_type`   | string | 3-letter amino acid code             | "ALA"      |
| `residue_number` | int    | Author residue sequence number       | 42         |
| `atom_id`        | string | Standard atom name (no quotes/spaces)| "CA"       |

### 3.2 Model Atom

Extends canonical identifier with Cartesian coordinates.

| Field   | Type  | Description       |
|---------|-------|-------------------|
| `x`     | float | Coordinate (Å)    |
| `y`     | float |                   |
| `z`     | float |                   |

### 3.3 Distance Restraint

```json
{
  "id": "dist_1",
  "atom_1": { "chain_id": "A", "residue_type": "ALA", "residue_number": 10, "atom_id": "CA" },
  "atom_2": { "chain_id": "A", "residue_type": "GLY", "residue_number": 25, "atom_id": "CA" },
  "min_dist": 2.0,
  "max_dist": 5.0
}
```

- `min_dist`, `max_dist` in Å.
- Restraint is satisfied if measured distance \( d \) satisfies `min_dist ≤ d ≤ max_dist`.

### 3.4 Dihedral Restraint

```json
{
  "id": "dih_1",
  "atom_1": { ... },
  "atom_2": { ... },
  "atom_3": { ... },
  "atom_4": { ... },
  "min_angle": -60.0,
  "max_angle":  30.0
}
```

- Four atoms define the torsion angle.
- `min_angle`, `max_angle` in degrees.
- Allowed set: angles θ such that `(θ - min_angle) mod 360 ≤ (max_angle - min_angle) mod 360`, with modulo yielding `[0,360)`.
- Special case: if `(max_angle - min_angle) mod 360 == 0`, the range is the whole circle (always satisfied).

### 3.5 Violation Object

```json
{
  "violation_id": "v_dist_1",
  "restraint_type": "distance",            // "distance" or "dihedral"
  "restraint_id": "dist_1",
  "atoms": [ <model atom references> ],
  "measured_value": 5.8,
  "expected_min": 2.0,
  "expected_max": 5.0,
  "violation_magnitude": 0.8
}
```

- `violation_magnitude` ≥ 0.
- Zero magnitude means the restraint is satisfied (no violation).
- Non-zero magnitude = shortest distance by which measured value lies outside allowed range, using the appropriate metric (linear distance for distance restraints, angular distance for dihedrals).

### 3.6 Residue-Level Density

Computed once per residue after mapping:

```json
{
  "chain_id": "A",
  "residue_number": 42,
  "residue_type": "ALA",
  "restraint_count": 17       // integer = number of restraints involving this residue
}
```

Interpretation: a proxy for experimental support, not a confidence score.

---

## 4. Data Pipeline

### 4.1 Data Retrieval

- Obtain model coordinates: download the first model (conformer 1) from the PDB in mmCIF format (e.g., `https://files.rcsb.org/download/<PDB_ID>.cif`).
- Obtain NMR restraints: download the NMR-STAR file from the PDB (e.g., `https://files.rcsb.org/download/<PDB_ID>.mr`).  
  V1 only parses:
  - Distance restraints (tag `_Dist_constraint` family or equivalent).
  - Dihedral angle restraints (tag `_Torsion_angle_constraint` family).
- Fallback: if NMR-STAR file unavailable or contains no parsable restraints, the system operates with zero restraints and displays an informational message.

### 4.2 Atom Mapping

The mapping layer reconciles restraint atom identifiers with model atom identifiers.  
**Assumptions for V1**:
- Model atom names are taken from `label_atom_id` in mmCIF.
- Restraint atom names follow standard PDB nomenclature; minimal cleanup: trim whitespace, remove quotes.
- Matching algorithm:
  1. Direct match: `chain_id`, `residue_number` (as integer), and `atom_id` (case‑insensitive) all identical.
  2. If direct match fails, attempt common alias substitutions (e.g., `"HN"` → `"H"`, `"1HA"` → `"HA1"`). A static alias table is provided in Appendix A.
  3. Unresolvable restraint atoms are logged and that restraint is excluded from further analysis.
- Mapping is performed once, and the resulting list of `MappedRestraint` objects (restraint + model atom references) is used throughout.

### 4.3 Geometric Evaluation

**Distance**  
Euclidean distance between two mapped atom coordinates:
```
d = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
violation_magnitude = max(0, min_dist - d, d - max_dist)
```

**Dihedral (torsion) angle**  
Standard atan2 computation using coordinates of four atoms:
```
θ = atan2( |b2| * a1 · (a2 × a3),  (a1 × a2) · (a2 × a3) )
```
where a1 = bond vector (atom1→atom2), a2 = (atom2→atom3), a3 = (atom3→atom4), and all vectors normalised.

After obtaining θ in degrees in range `[-180, 180)`, applied periodicity handling:

```
min_a = min_angle mod 360
max_a = max_angle mod 360
if (max_a - min_a) < 0: max_a += 360   // ensure interval width ≤ 360
θ_mod = θ mod 360
if θ_mod < min_a: θ_mod += 360
if min_a ≤ θ_mod ≤ max_a:
    violation_magnitude = 0
else:
    dist_to_min = min( |θ_mod - min_a|, |θ_mod - (min_a+360)| )
    dist_to_max = min( |θ_mod - max_a|, |θ_mod - (max_a+360)| )
    violation_magnitude = min(dist_to_min, dist_to_max)
```
This yields an angular deviation in degrees.

### 4.4 Violation Output

All restraints with `violation_magnitude > 0` become `Violation` objects. Violations are assigned stable IDs derived from restraint IDs.

### 4.5 Residue Density Computation

For each residue (chain, number, type), count distinct restraints where the residue appears in any mapped atom. Store as `restraint_count`.

---

## 5. State and Interaction Model

### 5.1 Application State

The entire UI behaviour is driven by a single, serialisable state object:

```json
{
  "pdb_id": "1ABC",
  "global_mode": "density",        // "density" | "distance" | "dihedral"
  "selected_residues": [           // list of residue keys
    {"chain_id": "A", "residue_number": 42}
  ],
  "selected_violation_id": null,   // or violation ID string
  "neighbourhood_radius": 5.0,     // Å
  "camera": null                   // optional saved camera orientation (quat/pos)
}
```

Default state: `global_mode = "density"`, no selections, PDB ID as input.

### 5.2 State Transitions & Events

| User Action                     | State Change                                                                          |
|---------------------------------|---------------------------------------------------------------------------------------|
| Input new PDB ID                | Reset state, reload data, set defaults.                                               |
| Click residue in Sequence View  | Toggle selection: if already selected, deselect; else add to `selected_residues`.    |
| Click row in Violation Table    | Set `selected_violation_id`, clear residue selections.                                |
| Toggle mode (Global/Local)      | When any residue selected → local view; empty selection → global view. Mode automatically derived. |
| Switch mode toggle              | Cycle `global_mode`: "density" → "distance" → "dihedral" → "density". In local view, mode fil‐ ters which restraints are shown. |
| Camera change (user rotates)    | `camera` updated (optional, for reproducible snapshot).                                |

**Combined logic**:  
- “Local view” is active when `selected_residues` is non‑empty **or** `selected_violation_id` is non‑null.  
- In local view, the 3D scene zooms to the relevant atoms; only nearby residues and restraints are shown (see 5.3).  
- When a violation row is selected, the scene focuses tightly on the atoms of that violation.

### 5.3 Local View Scope Definition

Given a set of **focus residues** (from selected residues or violations):

1. **Focus atoms**: all atoms of focus residues, plus (if a violation is selected) all atoms of that violation.
2. **Neighbourhood**: all residues that have at least one atom within `neighbourhood_radius` of any focus atom.
3. The displayed subset = focus residues ∪ neighbourhood residues.
4. Restraints shown in local view = those where **all** involved atoms belong to the displayed subset. (For dihedral restraints, this means all four atoms are in shown residues; for distance, both atoms.)
5. Opacity of non-focus residues may be reduced for visual emphasis (see Section 7).

---

## 6. UI Component Specifications

### 6.1 3D Molecular Viewer (Mol* via MolViewSpec)

- The viewer is a Mol* widget that receives a MolViewSpec (MVS) JSON scene definition.
- The scene is regenerated **completely** on every state change (deterministic, declarative).
- No built‑in Mol* selections or clicks; all interaction originates from external Sequence View and Violation Table.

Detailed scene assembly rules are in Section 7.

### 6.2 Sequence View

- Visual: a horizontal bar where each residue is represented as a rectangle (coloured by density).
- X‑axis: residue number (single chain assumed from the first model chain; multi‑chain handling out of scope for V1).
- Interaction: click/drag to select one or multiple residues. Shift‑click extends selection.
- Colour mapping: linear gradient based on `restraint_count`.  
  – Min → colour at zero count (light grey).  
  – Max → colour at highest density in protein (saturated blue to yellow‑orange from a viridis‑like scale; see Appendix B).  
  – Exact scale: `colour = interpolate( density_colors, count / max_count )` (clamped).
- Optional overlay: small red tick marks above residues that participate in at least one violation.

### 6.3 Violation Table

- Tabular display with columns:

| Column           | Description                                                    |
|------------------|----------------------------------------------------------------|
| Restraint ID     | Original restraint identifier (from file)                      |
| Type             | “Distance” or “Dihedral”                                       |
| Atoms            | List of involved atom identifiers (chain:resnum:atom)          |
| Measured         | Floating‑point measured value (Å or °)                         |
| Expected Range   | `[min, max]`                                                   |
| Violation Mag.   | Non‑negative number; sortable (default descending).            |

- Interaction: clicking a row sets `selected_violation_id` and triggers local view.
- Filtering: ability to show all violations, only distance, only dihedral, or those involving selected residues.

### 6.4 Control Bar

- PDB ID input (text box) with “Load” button.
- Global mode toggle buttons: “Density”, “Distance”, “Dihedral”.  
  – “Density”: shows residue‑coloured structure only (global view defaults).  
  – “Distance”: overlays distance restraint lines.  
  – “Dihedral”: shows dihedral violations via bond colours.  
  Active mode is highlighted.
- Local/Global indicator (auto‑determined; not a button).
- Reset selection button to return to global view.

---

## 7. MolViewSpec Scene Assembly Rules

MolViewSpec (MVS) is a JSON‑based declarative language for describing molecular scenes. The system generates MVS trees for each state.

### 7.1 Common Settings

- **Model**: only the first model (conformer) from the mmCIF file.
- **Representation**: cartoon for protein backbone; atoms displayed as small spheres only when needed for local detail.
- **Colouring**: residue‑level property mapping via `colorTheme` with per‑residue values.
- **Camera**: when local view active, include `focus` specifying the set of atoms (by selection expression) and optional `focusRadius`.

### 7.2 Global View MVS

When `selected_residues` empty and `selected_violation_id` null:

- **Density mode**:  
  - Cartoon coloured by `restraint_count` (property mapped).  
  - No distance lines, no bond emphasis.
- **Distance mode**:  
  - Cartoon coloured by density.  
  - All mapped distance restraints rendered as **Lines** (or thin cylinders) between the two atom positions.  
    - Colour: if violation magnitude == 0 → light grey with opacity 0.3; if violated → red with opacity linearly scaled from 0.3 at magnitude 0.01 to 1.0 at magnitude ≥ 2.0.  
    - Opacity clamping: max opacity 0.9.
- **Dihedral mode**:  
  - Cartoon coloured by density.  
  - Violated dihedral restraints only: colour the **central bond** (atom2–atom3) of the four atoms.  
    - Bond colour: red with opacity scaled similarly to distance mode based on violation magnitude.  
    - Non‑violated dihedrals not shown.

### 7.3 Local View MVS

Active when local view triggered by selection.

- **Visible residues**: only those in the displayed subset (focus + neighbours).
- **Cartoon**: colour focus residues with a highlight colour (e.g., yellow), neighbour residues with density colour but lowered opacity (0.5).
- **Ligand/water**: hidden.
- **Restraints**: only those satisfying the displayed subset condition.
  - Distance mode: lines as above, but all lines in the subset are drawn (both violated and non‑violated). Use the same colour/opacity scheme.
  - Dihedral mode: show violated dihedral central bonds as coloured lines. If a dihedral restraint is fully within subset but not violated, optionally show central bond in very faint grey (opacity 0.2) for context.
- **Camera**: `focus` set to the bounding sphere of all focus atoms with an additional 3 Å margin. `focusRadius` = radius of that sphere.

### 7.4 Geometry Generation Details

- **Distance lines**: For each restraint, generate a `line` primitive.  
  Coordinates: `[atom1.x, atom1.y, atom1.z]`, `[atom2.x, atom2.y, atom2.z]`.  
  Radius: 0.15 Å.
- **Dihedral central bonds**: For each violated dihedral, generate a line between atom2 and atom3 with radius 0.25 Å and colour red with appropriate opacity.
- All custom geometry is placed in a separate MVS component under a group labelled “Restraints”.

---

## 8. Operational Assumptions & Limitations

- **First conformer only**: all analysis is based on model 1 of the deposited ensemble.
- **Single chain**: V1 assumes the protein has one predominant chain. Multi‑chain targets will appear but may require manual chain selection (future enhancement).
- **Standard atom naming**: Atom mapping relies on reasonable atom name conventions; non‑standard names may cause incomplete mapping.
- **Restraints availability**: If no NMR‑STAR file or no parsable distance/torsion restraints, the system gracefully displays zero density, no violations.
- **No ensemble averaging**: Violations ignore ensemble spread; they represent a single snapshot.
- **No ambiguous/pseudoatom handling**: Such restraints are discarded if not resolvable by direct mapping.
- **No automated scoring or ML** – purely visual exploratory tool.

---

## 9. Implementation Plan (V1 Steps)

1. **Data retrieval module**  
   Download mmCIF and NMR‑STAR; parse out distance and dihedral restraint tables using a robust NMR‑STAR parser (e.g., `pynmrstar`).

2. **Atom mapping component**  
   Build lookup table; attempt direct match + alias expansion; output `MappedRestraint` list.

3. **Geometric evaluation**  
   Compute distances, dihedral angles, violations; produce list of `Violation` objects.

4. **Density computation**  
   Count restraints per residue.

5. **State management**  
   Implement state object and reducer/event handler.

6. **UI scaffolding**  
   Jupyter notebook with `ipywidgets` for sequence view (using `bqplot` bars or custom widget), violation table (`ipydatagrid` or `pandas` display + event), controls.

7. **MolViewSpec generation**  
   Functions to convert state into MVS JSON, handling all combinations (global density/distance/dihedral, local variations).

8. **Mol* integration**  
   Use Mol* Jupyter widget that accepts MVS updates on state change.

9. **Iteration and testing**  
   Validate with deposited NMR structures (e.g., 1UBQ, 2K39) where restraint data exists.

---

## Appendix A – Atom Name Alias Table

| Restraint Name | Model Alias  |
|----------------|--------------|
| HN             | H            |
| 1HA            | HA1          |
| 2HA            | HA2          |
| 1HB            | HB1          |
| 2HB            | HB2          |
| ...            | ...          |

(Expand as needed based on encountered NMR‑STAR files; V1 ships with a minimal but extensible dictionary.)

## Appendix B – Density Colour Scale

Uses a perceptually uniform sequential colormap (viridis).
- Low density (0): `#440154`
- Medium density: `#21918c`
- High density (max): `#fde725`

Linear interpolation in RGB space between defined knots.

---

*Document version 1.0 – ready for development.*
