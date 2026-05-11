This specification is structured to guide an engineering team through the development of a feature-rich, modular, and reproducible exploratory scientific tool. It adheres strictly to the principle of "specification first," ensuring that implementation decisions are driven by these explicit requirements.

***

# 🧬 NMR Structural Evidence Viewer Specification (V1.0)

**Document Status:** Specification Draft
**Target Audience:** Software Engineers, Scientific Computational Developers
**Goal:** To enable interactive, spatially-aware analysis of the experimental restraints supporting a deposited NMR structural model.

***

## 1. Goals and Scope

### 1.1. Core Problem Statement
The system addresses the gap in structural bioinformatics by transforming complex, raw NMR restraint data (e.g., distance and dihedral constraints) into an easily interpretable, localized visual context relative to the 3D atomic structure.

**Core Question:** For a given PDB structure, "Which regions are strongly supported, weakly supported, or locally inconsistent with the deposited NMR restraints?"

### 1.2. Success Criteria (Use Cases)
The implemented system must allow the user to achieve the following capabilities:

1.  **Global Support Mapping:** Visualize the overall density and quality of experimental support across the entire protein structure.
2.  **Local Inspection:** Select a residue or region to focus the view, revealing all associated restraints and violations in context.
3.  **Violation Deep Dive:** Select a specific violation (e.g., a single distance restraint) to view the involved atoms, the measured value, and the expected acceptable range.
4.  **Reproducible Reporting:** Generate a consistent analysis report derived from the pipeline state.

### 1.3. Operational Boundaries (Scope & Limitations)
The V1 system **must** operate under the following constraints:

*   **Conformer Limitation:** Only the coordinates and data associated with the **first deposited conformer** (or the primary model) will be analyzed. Ensemble statistics are explicitly out of scope.
*   **Data Source Reliance:** The system is entirely dependent on the quality and availability of external inputs (PDB/mmCIF coordinates and associated NMR-STAR restraint files).
*   **Model Type:** The analysis is designed to interpret the structure as an *inferred model*, never as a direct observation.
*   **Exclusions (V1):** Ensemble averaging, dynamic averaging, ambiguous restraints, pseudo-atoms, and automated quality scoring beyond explicit violation computation are **not** implemented.

***

## 2. System Architecture and Data Flow

The system must execute as a deterministic, linear pipeline triggered by a single primary input (`PDB ID`).

### 2.1. Workflow Diagram (State Machine Flow)
$$\text{PDB ID} \xrightarrow{1} \text{Data Retrieval} \xrightarrow{2} \text{Annotation/Mapping} \xrightarrow{3} \text{Evaluation Engine} \xrightarrow{4} \text{Violation Data} \xrightarrow{5} \text{Visualization View}$$

### 2.2. Core Components and Responsibilities

| Component | Responsibility | Inputs | Outputs |
| :--- | :--- | :--- | :--- |
| **Data Retrieval Layer** | Fetches geometry and raw restraints. | `PDB ID` | Coordinates ($\mathbf{C}$), Raw Restraint Files ($\mathbf{R}_{raw}$) |
| **Atom Annotation Layer** | Maps raw restraint IDs to canonical structure atoms. | $\mathbf{C}$, $\mathbf{R}_{raw}$ | Mapped Atomic Coordinates ($\mathbf{C}'$), Mapped Restraint Set ($\mathbf{R}_{mapped}$) |
| **Evaluation Engine** | Computes geometric deviations and densities. | $\mathbf{C}'$, $\mathbf{R}_{mapped}$ | Violation Set ($\mathbf{V}$), Density Map ($\mathbf{D}$) |
| **Visualization Layer** | Interprets $\mathbf{V}$ and $\mathbf{D}$ to generate interactive graphical states. | $\mathbf{V}$, $\mathbf{D}$, User Selections | 3D Scene State, Sequence View State, Table State |

***

## 3. Data Model Specification

To ensure technical precision, all data structures must be strictly defined.

### 3.1. Canonical Atomic Coordinate Object (State)
Every atom involved in the model or a restraint must conform to this structure.

| Field | Data Type | Description | Constraints |
| :--- | :--- | :--- | :--- |
| `chain_id` | String | Single character chain identifier (e.g., 'A', 'B'). | Mandatory |
| `residue_type` | String | Three-letter amino acid type (e.g., 'ALA', 'GLY'). | Mandatory |
| `residue_number` | Integer | Sequence position (1-indexed). | Mandatory |
| `atom_id` | String | Specific atom name (e.g., 'CA', 'N', 'C'). | Mandatory |
| `coordinates` | Vector(3) | Cartesian coordinates $(\text{X}, \text{Y}, \text{Z})$ in Å. | Mandatory (Real) |

### 3.2. Restraint Data Model ($\mathbf{R}_{mapped}$)
This structure encapsulates the necessary information for a single restraint.

#### 3.2.1. Distance Restraint ($\mathbf{R}_{dist}$)
*   **Definition:** A pair of atoms ($A, B$) expected to be within a defined range.
*   **Fields:**
    *   `type`: "Distance"
    *   `atom_A`: Canonical Atomic Object (Reference ID)
    *   `atom_B`: Canonical Atomic Object (Reference ID)
    *   `min_dist`: Float (Minimum allowed distance, Å)
    *   `max_dist`: Float (Maximum allowed distance, Å)

#### 3.2.2. Dihedral Restraint ($\mathbf{R}_{dihed}$)
*   **Definition:** A sequence of four atoms ($A, B, C, D$) expected to form a torsion angle ($\theta$) within a defined range.
*   **Fields:**
    *   `type`: "Dihedral"
    *   `atom_1`, `atom_2`, `atom_3`, `atom_4`: Canonical Atomic Object (Reference IDs)
    *   `min_angle`: Float (Minimum allowed angle, degrees)
    *   `max_angle`: Float (Maximum allowed angle, degrees)

### 3.3. Violation Data Model ($\mathbf{V}$)
A violation is generated when the measured geometry deviates from the restraint's allowed range.

| Field | Data Type | Description | Source |
| :--- | :--- | :--- | :--- |
| `associated_restraint` | Object | The original $\mathbf{R}_{dist}$ or $\mathbf{R}_{dihed}$ object. | $\mathbf{R}_{mapped}$ |
| `involved_atoms` | List[Atomic Object] | Atoms contributing to the measured geometry. | $\mathbf{C}'$ |
| `measured_value` | Float | The computed geometric value (distance or angle). | Evaluation Engine |
| `expected_range` | Tuple(Float, Float) | The allowed $[min, max]$ range. | $\mathbf{R}_{mapped}$ |
| `violation_magnitude` | Float | The absolute deviation magnitude beyond the range limits. | Calculated |

### 3.4. Support Density Model ($\mathbf{D}$)
Density is a **residue-level metric**.

$$\text{Density}(R) = \text{Count}(\text{All } \mathbf{R} \in \mathbf{R}_{mapped} \text{ associated with residue } R)$$

***

## 4. Functional Requirements (The Engine Logic)

The Evaluation Engine must perform the following calculations to produce $\mathbf{V}$ and $\mathbf{D}$.

### 4.1. Geometry Measurement
The system must calculate the following measurements based on the canonical coordinates $\mathbf{C}'$:

1.  **Measured Distance:** Standard Euclidean distance calculation between two points.
2.  **Measured Dihedral Angle:** Calculation of the angle formed by four consecutive atoms.

### 4.2. Violation Computation Logic
For every restraint $\mathbf{R}$, the violation calculation must be executed:

#### A. Distance Violation ($\mathbf{R}_{dist}$)
Let $d_{meas}$ be the measured distance and $[d_{min}, d_{max}]$ be the expected range.
*   **Condition for Violation:** $d_{meas} < d_{min}$ OR $d_{meas} > d_{max}$.
*   **Violation Magnitude:**
    $$\text{Magnitude} = \begin{cases} |d_{meas} - d_{min}| & \text{if } d_{meas} < d_{min} \\ |d_{meas} - d_{max}| & \text{if } d_{meas} > d_{max} \\ 0 & \text{otherwise} \end{cases}$$

#### B. Dihedral Violation ($\mathbf{R}_{dihed}$)
Let $\theta_{meas}$ be the measured angle and $[\theta_{min}, \theta_{max}]$ be the expected range.
*   **Note:** Angles must be handled modulo $360^\circ$.
*   **Violation Magnitude:**
    $$\text{Magnitude} = \begin{cases} |\theta_{meas} - \theta_{min}| & \text{if } \theta_{meas} < \theta_{min} \\ |\theta_{meas} - \theta_{max}| & \text{if } \theta_{meas} > \theta_{max} \\ 0 & \text{otherwise} \end{cases}$$

***

## 5. User Interface and Visualization Specification

The Visualization Layer must adhere to a **state-driven pattern**. Changes in one UI component (e.g., clicking a residue in the Sequence View) must update the other components (3D View, Violation Table) without re-running the entire calculation pipeline.

### 5.1. Global View (Default State)
*   **Goal:** Provide a rapid assessment of global structural integrity.
*   **Visualization:**
    *   **Residue Coloring:** Color the backbone (or C$\alpha$) of each residue based on its calculated **Restraint Density ($\mathbf{D}$)**.
        *   *Color Scheme:* Implement a monotonic gradient (e.g., low density = blue, high density = red).
    *   **Violation Markers:** Optionally render small, non-interactive markers on residues with a high density of violations ($\mathbf{V}$).
    *   **Clarity Priority:** The initial view must be low-clutter, prioritizing the overall structure's orientation.

### 5.2. Local View (Activated State)
This view is activated by user interaction (e.g., selecting a residue $R_s$ in the Sequence View).

*   **Goal:** Detail the local experimental evidence for the selected area.
*   **Visualization:**
    *   **Target Residues:** $R_s$ and neighboring residues (e.g., $\pm 3$ residues) must be visibly highlighted (e.g., slightly larger spheres, specific border color).
    *   **Restraint Rendering:** Display all $\mathbf{R}$ associated with $R_s$.
        *   **Distance:** Render as a line connecting the two restrained atoms.
        *   **Dihedral:** Render by coloring the two bonds defining the geometry (e.g., $\text{A} \rightarrow \text{B}$ and $\text{C} \rightarrow \text{D}$).
    *   **Encoding (Crucial):** The line/bond rendering must encode the **Violation Magnitude**:
        *   **Color:** Scale color intensity/hue based on $\text{Violation Magnitude}$. (e.g., $\text{Magnitude} \approx 0 \rightarrow \text{Green/Gray}$; $\text{Magnitude} \rightarrow \text{Red}$).
        *   **Opacity/Transparency:** Use transparency to modulate saliency. Highly violating restraints should be opaque; perfectly supported restraints should be highly transparent or light gray.

### 5.3. Violation Table Component
*   **Role:** A filterable, sortable list of all $\mathbf{V}$.
*   **Row Data:** Must display the associated restraint type, the involved atoms, the measured value, the expected range, and the violation magnitude.
*   **Interaction:** Selecting a row (a specific $\mathbf{V}$) must trigger the **Local View** to automatically jump to and highlight the involved atoms and restraint geometry in the 3D structure.

### 5.4. Sequence View Component
*   **Role:** Primary navigation tool.
*   **Visualization:** Each residue must be represented as a navigable element.
*   **Encoding (Future):** The visual state of the residue element should optionally encode:
    *   **Density:** Background color gradient (reflecting $\mathbf{D}$).
    *   **Violation Burden:** A small symbol (e.g., a count or icon) indicating the total number of associated violations.

***

## 6. Implementation Notes and Assumptions

1.  **Coordinate System:** All calculations assume a consistent right-handed Cartesian coordinate system (Ångströms).
2.  **Angle Normalization:** Dihedral angles must be handled cyclically (e.g., $-180^\circ \equiv +180^\circ$).
3.  **Separation of Modes:** The rendering of Distance restraints and Dihedral restraints must be mutually exclusive (i.e., they cannot be displayed simultaneously in the same manner). A mode toggle is required for switching between these distinct visualization types.
4.  **Performance:** Since the evaluation engine is computationally intensive, the system must utilize efficient data structures (e.g., cached intermediate results) to ensure that UI updates are near-instantaneous once the initial $\mathbf{V}$ and $\mathbf{D}$ sets are calculated.