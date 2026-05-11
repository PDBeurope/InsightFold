# Technical Specification: NMR Experimental Support & Violation Viewer

## 1. Project Overview
This system is an exploratory report pipeline designed to visualize the relationship between NMR structural models and their corresponding experimental restraints. It transforms raw PDB/mmCIF coordinates and NMR-STAR restraint data into an interactive 3D environment, allowing researchers to identify regions of high experimental support, sparse constraints, and local structural inconsistencies.

### 1.1 Objectives
- **Quantify Experimental Support:** Map restraint density to residues.
- **Identify Inconsistencies:** Calculate and visualize geometric violations of distance and dihedral restraints.
- **Spatially Contextualize Data:** Provide a seamless link between textual violation reports, sequence positions, and 3D coordinates.
- **Ensure Reproducibility:** Use a state-driven visualization approach via `MolViewSpec`.

### 1.2 Scope & Constraints
- **Model Selection:** Analysis is performed exclusively on the **first conformer (Model 1)** of the deposited ensemble.
- **Exclusions:** This version (V1) excludes ensemble statistics, pseudoatoms, ambiguous restraints, and automated quality scoring.
- **Input:** Requires a valid PDB ID with associated NMR-STAR files.

---

## 2. System Architecture

### 2.1 High-Level Pipeline
`PDB ID` $\rightarrow$ `Data Retrieval` $\rightarrow$ `Restraint Parsing` $\rightarrow$ `Atom Mapping` $\rightarrow$ `Geometric Evaluation` $\rightarrow$ `Violation Computation` $\rightarrow$ `Interactive Visualization`.

### 2.2 Component Decomposition

#### 2.2.1 Data Retrieval Layer
- **Coordinates:** Fetch `.cif` (mmCIF) files from the PDB.
- **Restraints:** Fetch NMR-STAR files containing distance and dihedral restraint data.
- **Output:** Raw coordinate arrays and restraint lists.

#### 2.2.2 Atom Annotation & Mapping Layer
The system must resolve the mapping between the restraint file identifiers and the coordinate file identifiers.
- **Canonical Atom Key:** A unique tuple consisting of `(chain_id, residue_type, residue_number, atom_id)`.
- **Mapping Process:** 
    1. Parse NMR-STAR restraint atom indices.
    2. Resolve indices to Canonical Atom Keys.
    3. Cross-reference keys with mmCIF coordinate entries to retrieve Cartesian coordinates $(x, y, z)$.

#### 2.2.3 Restraint Data Models
**Distance Restraint:**
- `id`: Unique identifier.
- `atoms`: Pair of Canonical Atom Keys $[A, B]$.
- `range`: Tuple of $(\text{distance}_{min}, \text{distance}_{max})$ in Angstroms ($\text{\AA}$).

**Dihedral Restraint:**
- `id`: Unique identifier.
- `atoms`: Sequence of four Canonical Atom Keys $[A, B, C, D]$.
- `range`: Tuple of $(\text{angle}_{min}, \text{angle}_{max})$ in degrees.

---

## 3. Computational Logic

### 3.1 Geometric Evaluation
For every restraint, the system must calculate the observed value from the structural coordinates.

- **Distance Calculation:** 
  $$d_{obs} = \sqrt{(x_2-x_1)^2 + (y_2-y_1)^2 + (z_2-z_1)^2}$$
- **Dihedral Calculation:** 
  Compute the torsion angle $\phi$ defined by atoms $A \rightarrow B \rightarrow C \rightarrow D$ using the standard vector cross-product method.

### 3.2 Violation Computation
A violation occurs if the observed value falls outside the specified range.

- **Distance Violation Magnitude ($V_d$):**
  $$V_d = \max(0, \text{distance}_{min} - d_{obs}, d_{obs} - \text{distance}_{max})$$
- **Dihedral Violation Magnitude ($V_\phi$):**
  $$V_\phi = \max(0, \text{angle}_{min} - \phi_{obs}, \phi_{obs} - \text{angle}_{max})$$
- **Violation Object:** Contains `restraint_id`, `atoms`, `measured_value`, `expected_range`, and `magnitude`.

### 3.3 Restraint Density Calculation
Restraint density is a residue-level metric used to proxy experimental support.
- **Definition:** For each residue $R_i$, count the number of distance and dihedral restraints where at least one atom belongs to $R_i$.
- **Formula:** $\text{Density}(R_i) = \sum (\text{all restraints associated with } R_i)$.

---

## 4. Visualization Specification

### 4.1 Framework & State Management
- **Engine:** Mol* (rendering surface).
- **State Controller:** `MolViewSpec` (declarative state).
- **Logic:** The notebook controls the state; Mol* simply renders the specified state.

### 4.2 View Modes

#### 4.2.1 Global View (Support Mapping)
- **Purpose:** Overview of experimental coverage.
- **Visual Encoding:** 
    - **Residue Coloring:** Gradient scale based on **Restraint Density**.
    - **Default Palette:** (e.g., Low Density $\rightarrow$ Light Blue; High Density $\rightarrow$ Dark Blue).
    - **Filtering:** Option to overlay high-magnitude violations.

#### 4.2.2 Local View (Evidence Inspection)
- **Purpose:** Detailed inspection of specific restraints.
- **Activation:** Triggered via Sequence View selection or Violation Table click.
- **Visual Encoding:**
    - **Distance Restraints:** Rendered as lines/cylinders between atom pairs.
    - **Distance Color:** Color-coded by `Violation Magnitude` (e.g., Green = Valid, Yellow = Weak Violation, Red = Strong Violation).
    - **Distance Saliency:** Adjust transparency based on magnitude to prevent clutter.
    - **Dihedral Restraints:** Highlight the specific bond/residue geometry associated with the torsion.

---

## 5. User Interface Components

### 5.1 Sequence View
- **Function:** Primary navigation tool.
- **Interactions:**
    - Clicking a residue $\rightarrow$ Centers 3D view on residue $\rightarrow$ Triggers Local View for that residue.
    - Visual encoding (Optional): Color bar indicating restraint density per residue.

### 5.2 Violation Table
- **Columns:** `Restraint ID`, `Type`, `Atoms`, `Measured Value`, `Range`, `Violation Magnitude`.
- **Interactions:**
    - Sorting by `Violation Magnitude`.
    - Clicking a row $\rightarrow$ Visualizes that specific restraint in 3D $\rightarrow$ Zooms to involved atoms.

### 5.3 3D Structure View
- **Function:** Spatial interpretation.
- **Controls:** Toggles for "Distance Restraints" vs "Dihedral Restraints" (separated modes to avoid visual noise).

---

## 6. Operational Summary Table

| Feature | Implementation Detail | Requirement |
| :--- | :--- | :--- |
| **Data Input** | PDB ID | Must have NMR-STAR data |
| **Model** | Model 1 Only | Ignore subsequent conformers |
| **Atom Mapping** | Canonical Key | `(chain, res_type, res_num, atom_name)` |
| **Support Metric** | Count per residue | Inclusive of all restraint types |
| **Violation** | $\Delta$ outside range | Calculated for distance and torsion |
| **Viz Engine** | Mol* $\rightarrow$ MolViewSpec | State-driven, reproducible |
| **Local View** | Line/Color encoding | Color $\propto$ Violation Magnitude |