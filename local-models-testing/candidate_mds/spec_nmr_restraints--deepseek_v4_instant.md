# NMR Restraint Visualization Specification

## 1. Introduction

### 1.1 Problem Statement
Users need to understand how well an NMR structural model is supported by, and consistent with, the experimental restraints used during structure determination. Specifically, the system answers: *“Which regions of an NMR structure are experimentally supported, weakly constrained, or locally inconsistent with the deposited NMR restraints?”*

### 1.2 Target Audience
- Protein biochemists, structural biology users, and researchers familiar with protein structure but not necessarily NMR refinement experts.
- Users exploring NMR structures in the Protein Data Bank (PDB).

### 1.3 Core Capabilities
The system enables users to:
- Visually identify experimentally supported vs. weakly restrained regions.
- Inspect violated restraints in structural context.
- Determine if violations are isolated or locally clustered.
- Explore local experimental support networks around residues.
- Relate restraint data directly to 3D structure.

**Scope:** The notebook is an *exploratory experimental evidence viewer*, not a full automated validation system.

### 1.4 Guiding Principles
- **Evidence transparency** over geometric validation scoring.
- **Contextual interpretation** over aggregate scores.
- **Separation of support from conflict** – no global confidence scores.
- **Reproducibility** through declarative, state-driven visualization.

---

## 2. System Architecture & Workflow

### 2.1 High-Level Pipeline
```
PDB ID → Data Retrieval → Restraint Parsing → Atom Mapping → Geometric Evaluation → Violation Computation → Interactive Visualization
```

### 2.2 Output Artifacts
The notebook produces:
- A reproducible analysis report (HTML/notebook format).
- Interactive 3D visualization (Mol* viewer).
- Sequence-based navigation widget.
- Violation inspection table.

---

## 3. Data Retrieval Specification

### 3.1 Input
- **PDB ID** (4-character alphanumeric string, e.g., `2K39`).

### 3.2 Data Sources
| Data Type | Source Format | Access Method |
|-----------|---------------|----------------|
| Structure coordinates | PDB/mmCIF | Direct download from PDB (e.g., `https://files.rcsb.org/download/{PDB_ID}.cif`) |
| Distance restraints | NMR-STAR | Retrieved via PDB or BMRB (e.g., `https://bmrb.io/ftp/pub/bmrb/entries/{PDB_ID}/nmr-star/{PDB_ID}_nmr-star_v3.txt`) |
| Dihedral restraints | NMR-STAR | Same file as above |

### 3.3 Required Data Elements
- **Coordinates:** Atom identifier (chain, residue number, residue type, atom name), Cartesian coordinates (x,y,z).
- **Restraints:** For distance: atom pair identifiers, min distance, max distance. For dihedral: four atom identifiers, min angle, max angle.

### 3.4 Retrieval Failure Handling
- If any required data source is unavailable → show descriptive error and abort.
- If NMR-STAR file has no distance or dihedral restraints → show warning and continue with empty restraint sets.

---

## 4. Atom Mapping Layer

### 4.1 External Dependency
The atom mapping layer is considered **external to the visualization logic**. It must provide normalized atom mappings with the following canonical representation.

### 4.2 Canonical Atom Representation
| Field | Type | Description |
|-------|------|-------------|
| `chain_id` | string (1-2 chars) | Chain identifier (e.g., "A", "B") |
| `residue_type` | string (3-letter code) | Standard amino acid code (e.g., "ALA") |
| `residue_number` | integer | Sequence position |
| `atom_id` | string | Standard atom name (e.g., "CA", "N", "C") |
| `coordinates` | (float, float, float) | Cartesian coordinates (Å) |

### 4.3 Mapping Requirements
The mapping layer must:
- Reconcile restraint identifiers (from NMR-STAR) with structure atom identifiers (from mmCIF).
- Resolve potential naming differences (e.g., "H" vs "HN") using standard residue atom name tables.
- Report unresolvable atoms as warnings; skip restraints containing unresolvable atoms.
- Ensure both atoms in a distance restraint map uniquely to model atoms; otherwise discard that restraint.

### 4.4 Interface Expectation
The visualization system assumes a function or module with the signature:
```python
def map_restraints_to_model(restraints: RestraintSet, model: StructureModel) -> MappedRestraintSet
```
Output includes mapped atom references (indices or full identifiers) and the original restraint bounds.

---

## 5. Restraint Data Models

### 5.1 Distance Restraint
| Field | Type | Description |
|-------|------|-------------|
| `id` | string (optional) | Unique restraint identifier from source, if any |
| `atom1` | CanonicalAtom | First restrained atom |
| `atom2` | CanonicalAtom | Second restrained atom |
| `d_min` | float (Å) | Lower bound of allowed distance |
| `d_max` | float (Å) | Upper bound of allowed distance |

**Interpretation:** The distance is expected to lie within `[d_min, d_max]`. Both bounds are inclusive.

### 5.2 Dihedral Restraint
| Field | Type | Description |
|-------|------|-------------|
| `id` | string (optional) | Unique identifier |
| `atom1`, `atom2`, `atom3`, `atom4` | CanonicalAtom | Atoms defining the torsion angle (in order) |
| `angle_min` | float (degrees) | Lower bound of allowed torsion |
| `angle_max` | float (degrees) | Upper bound of allowed torsion |

**Interpretation:** The torsion angle (0–360°) is expected to lie within `[angle_min, angle_max]` (circular). For simplicity V1 uses linear comparison after normalizing measured angle to [0,360).

---

## 6. Geometric Evaluation Layer

### 6.1 Distance Evaluation
For each distance restraint:
- Compute Euclidean distance: `d_measured = ||pos(atom1) - pos(atom2)||`
- Violation exists if `d_measured < d_min` or `d_measured > d_max`.

**Violation magnitude:**
```
if d_measured < d_min:
    magnitude = d_min - d_measured   (positive, units Å)
elif d_measured > d_max:
    magnitude = d_measured - d_max   (positive, units Å)
else:
    magnitude = 0.0
```

### 6.2 Dihedral Evaluation
For each dihedral restraint:
- Compute torsion angle from four atom coordinates using standard formula (atan2 of dot products of normals).
- Normalize measured angle to [0, 360) degrees.
- Determine if within `[angle_min, angle_max]`. Because bounds may cross 0°, the condition is:
  - If `angle_min <= angle_max` → within if `angle_min <= measured <= angle_max`
  - Else (wrap-around) → within if `measured >= angle_min` or `measured <= angle_max`
- Violation exists if not within the interval.

**Violation magnitude (angular deviation):**
```
def angular_distance(a, b):
    diff = abs(a - b) % 360
    return min(diff, 360 - diff)

if measured is outside interval:
    # distance to nearest bound
    dist_to_min = angular_distance(measured, angle_min)
    dist_to_max = angular_distance(measured, angle_max)
    magnitude = min(dist_to_min, dist_to_max)   # degrees
else:
    magnitude = 0.0
```

### 6.3 Output
For each evaluated restraint, store:
- `measured_value` (distance or angle)
- `violation_magnitude` (float, ≥0)
- `is_violated` (boolean)

---

## 7. Violation Data Model

A violation is a derived entity containing:

| Field | Type | Description |
|-------|------|-------------|
| `restraint_id` | string | Reference to original restraint (if exists) |
| `restraint_type` | enum | `"distance"` or `"dihedral"` |
| `atoms_involved` | list[CanonicalAtom] | List of atoms (2 for distance, 4 for dihedral) |
| `measured_value` | float | Distance (Å) or angle (degrees) |
| `expected_min` | float | Lower bound |
| `expected_max` | float | Upper bound |
| `violation_magnitude` | float | As defined in Section 6 |
| `residue_impact` | set[int] | Set of residue numbers directly involved (all atoms' residue numbers) |

---

## 8. Restraint Density Computation

### 8.1 Definition
For each residue position `r`:
```
restraint_density[r] = number of restraints (distance or dihedral) that involve at least one atom in residue r
```

**Counting rules:**
- A distance restraint between atoms in residues `i` and `j` contributes +1 to density of *both* residues.
- A dihedral restraint involving atoms in residues `i, j, k, l` contributes +1 to density of each distinct residue among them.
- Duplicate counting: each restraint counted at most once per residue.

### 8.2 Interpretation
- **Proxy for experimental support** – not a direct confidence score.
- Higher density → more experimental information anchoring that residue.
- Zero density → no direct restraints, weak experimental support.

### 8.3 Normalization for Visualization (Optional)
For color mapping, density values may be normalized across the protein:
```python
density_max = max(density values over all residues)
if density_max > 0:
    normalized = density / density_max
else:
    normalized = 0 for all
```

---

## 9. Visualization System Specification

### 9.1 Technology Stack
- **Mol***: Provides 3D rendering surface and scene management.
- **MolViewSpec**: Declarative specification language to define reproducible visual states. Used to generate scenes without direct Mol* API calls.

### 9.2 Architectural Constraint
Mol* acts **only as a rendering surface**, not the central interaction controller. All UI state (selection, filtering, coloring) is managed externally by the notebook (Python/IPywidgets) and passed to MolViewSpec.

### 9.3 State-Driven Model
The visualization is a pure function of a `VisualizationState` object that includes:
- Current PDB ID
- Restraint density per residue
- List of violations
- Selected residues (set of residue numbers)
- Selected violations (by ID)
- Active filter (distance vs dihedral, violation magnitude threshold)
- Visualization mode (`"global"` or `"local"`)

Changes to state regenerate the MolViewSpec JSON and re-render.

---

## 10. Global View Specification

### 10.1 Purpose
Display experimental support distribution across the entire structure.

### 10.2 Default Visual Encoding
- **Protein backbone** (cartoon representation) colored by restraint density:
  - Color map: sequential (e.g., blue = low density → red = high density).
  - Normalization as in Section 8.3.
- **Side chains** optionally shown in wireframe/low opacity (default off).
- **Violations** (optional, default only top N major violations to reduce clutter):
  - Show distance violations as dashed lines (if lines are enabled).
  - Show no dihedral violations in global view (defer to local view).
- Transparency: None for cartoon. For distance restraint lines, transparency determined by violation magnitude (lower opacity for small violations).

### 10.3 Clutter Control
- Do not show all restraints in global view.
- Show only violations exceeding a user‑settable threshold (default: magnitude > 0.5Å for distance, > 10° for dihedral).
- Maximum number of simultaneously displayed violation lines: 200 (configurable).

---

## 11. Local View Specification

### 11.1 Activation Triggers
Local view activates when:
- User selects a residue in the sequence view.
- User selects a row in the violation table.
- User applies a filtering control that isolates a region (e.g., residue range slider).

### 11.2 Visualization Content
When a residue is selected:
- Show the selected residue and all residues within a **spherical radius of 8 Å** (CA–CA) as opaque cartoon.
- All other residues shown as transparent cartoon or hidden (user preference).
- **Show all restraints** (both satisfied and violated) involving any atom in the selected residue(s):
  - Distance restraints: lines between atoms.
  - Dihedral restraints: color the central bond (between atoms 2 and 3) with a hue representing violation severity.
- Labels: optionally show residue numbers and restraint distances.
- Violations are highlighted with increased line thickness (e.g., 3× normal) and a pulsating effect (if supported).

### 11.3 Local Violation Context
When a violation is selected from the table:
- Center view on the geometric centroid of involved atoms.
- Show the mutant as in local view, plus all restraints involving those atoms.
- Display a tooltip or side panel with:
  - Restraint bounds and measured value.
  - Violation magnitude.
  - Atom identifiers and residue numbers.

---

## 12. Distance Restraint Visualization Encoding

| Property | Encoding |
|----------|----------|
| **Geometry** | Cylinder/line between atom centers. |
| **Color (satisfied)** | Gray (opacity 0.3). |
| **Color (violated)** | Gradient from yellow (small violation) to red (large violation). |
| **Line width** | 0.05 Å for satisfied; 0.15 Å for violated. |
| **Transparency** | Alpha = max(0.2, 1.0 - (violation_magnitude / max_expected_violation)). Or fixed 0.5 for satisfied. |
| **Dashed** | Violated restraints displayed as dashed lines (pattern: 5px dash, 2px gap). |

**Color mapping for violation magnitude:**
- Distance violation (Å): 0 → yellow `#FFFF00`, 5Å or more → red `#FF0000` (clamped).
- Dihedral violation (degrees): 0 → yellow, 180° → red.

---

## 13. Dihedral Restraint Visualization

Because dihedral restraints involve four atoms, they are visualized by **coloring the bond between the middle two atoms** (atoms 2 and 3 in the torsion definition).

| Property | Encoding |
|----------|----------|
| **Geometry** | Cylinder between atom2 and atom3. |
| **Color (satisfied)** | Gray (opacity 0.3). |
| **Color (violated)** | Same gradient as distance (yellow to red) based on violation magnitude. |
| **Line width** | 0.1 Å. |
| **Mode separation** | Dihedral restraints are not shown simultaneously with distance restraints. A mode toggle selects which type is visible in the 3D view. Default: distance restraints. |

---

## 14. Interaction Model (External Controllers)

### 14.1 Principles
- Mol* does **not** manage selection or filtering state.
- All interactions go through notebook widgets (ipywidgets, Plotly, etc.) that update the `VisualizationState` and regenerate the MolViewSpec.

### 14.2 Sequence View
- Component: Interactive residue sequence (e.g., using `ipywidgets.HBox` of buttons or a custom slider/glyph).
- **Encoding:**
  - Background color of each residue = restraint density (normalized).
  - Optional star/icon for residues with ≥1 violation.
- **Actions:**
  - Click on residue → add to selected residues set (local view). Multiple selection possible (Shift+click).
  - Clear selection button.

### 14.3 Violation Table
- Component: HTML table or `ipywidgets.DataGrid`.
- Columns:
  | Type | Magnitude | Involved Residues | Expected (range) | Measured |
  |------|-----------|------------------|------------------|-----------|
  | Distance | 1.2Å | 45-47 | [3.0, 5.0] Å | 6.2 Å |
- **Actions:**
  - Click row → select that violation, trigger local view.
  - Sort by any column (magnitude, residue).
  - Filter by type, magnitude threshold, residue range.

### 14.4 Filtering Controls
- **Violation magnitude slider** (distance: 0–10Å, dihedral: 0–180°).
- **Restraint type toggle** (distance / dihedral / both).
- **Residue range slider** (min, max residue numbers).
- **Density threshold slider** (show only residues with density ≥ threshold in global view – optional).

### 14.5 Mode Toggles
- **Global vs Local view** (radio button).
- **Show satisfied restraints** (checkbox, default off in global view, on in local view).
- **Show side chains** (checkbox, default off).

---

## 15. Main UI Components (Detailed)

### 15.1 Layout
The notebook should present a three‑panel layout:

| Panel | Content |
|-------|---------|
| Left (30% width) | Sequence view + Violation table + Filter controls |
| Right (70% width) | 3D Mol* view |
| Bottom (collapsible) | Log/output messages, restraint counts, metadata |

### 15.2 Sequence View Detailed Spec
- Represented as a horizontal or vertical strip of residue boxes.
- Dimensions: each residue box 30px wide, 20px high (scalable).
- Tooltip on hover: residue number, type, density, violation count.
- On selection: visual highlight (black border, glow).

### 15.3 Violation Table Detailed Spec
- Pagination if >100 rows (default 20 rows per page).
- Search box to filter by residue number or atom name.
- Export to CSV button.

### 15.4 3D View Initialization
- Mol* viewer embedded via `ipymolstar` or similar.
- Default camera: center of mass of the protein, distance such that entire protein fits in view.
- Background: white or light gray.

---

## 16. Operational Assumptions

The system assumes:
1. Deposited restraints are available in NMR-STAR format for the given PDB ID.
2. Atom mappings provided by the external mapping layer are correct.
3. Restraint files are compatible with model files (same chemical components).
4. The first conformer (model 1) is representative of the ensemble; analysis is done only on model 1.
5. The structure contains at least one chain with amino acids (no nucleic acids in V1).
6. Input PDB ID is valid and publicly accessible.

---

## 17. Known Limitations (V1)

The following features are **explicitly excluded** from V1:
- Ensemble statistics or multi‑model averaging.
- Ambiguous restraints (e.g., distance restraints with multiple possible atom pairings).
- Pseudoatom handling (e.g., methyl group averaging).
- Ensemble averaging semantics for restraint satisfaction.
- Dynamic averaging or mobility analysis.
- Automated quality scoring or machine learning interpretation.
- Nucleic acid or ligand restraints.
- Manual adjustment of atom mappings.
- Export of restraint lists or violation reports in custom formats (beyond CSV for violation table).

Only the **first model/conformer** is analyzed.

---

## 18. Output Specification

The notebook produces:
1. **A rendered HTML notebook** containing:
   - Interactive widgets (sequence view, table, sliders).
   - Embedded Mol* viewer with MolViewSpec generated state.
2. **In‑notebook text summary**:
   - Total number of distance restraints parsed.
   - Total number of dihedral restraints parsed.
   - Count of violated distance restraints.
   - Count of violated dihedral restraints.
   - Number of residues with zero restraints (warning if >10% of residues).
3. **Console/Log output** (optional for debugging): warnings for unresolvable atoms, skipped restraints.

---

## 19. Error Handling & Edge Cases

| Scenario | Behavior |
|----------|----------|
| PDB ID not found | Display error message, stop execution. |
| NMR-STAR file missing | Display warning, continue with empty restraints (all density = 0, no violations). |
| No distance restraints | Show "No distance restraints available" in UI, disable distance‑related features. |
| No dihedral restraints | Similar symmetric handling. |
| Atom mapping fails for a restraint | Log warning, skip that restraint (do not count it in density). |
| Residue numbering mismatch (e.g., insertion codes) | Normalize to integer numbering; ignore insertion codes (warn). |
| Structure has >1 chain | Analyze all chains independently; sequence view concatenates chains with chain letter headers. |
| Large protein (>500 residues) | Use downsampling of satisfied restraints in global view; keep all violations. Provide user warning. |

---

## 20. Performance Considerations

- **Restraint geometry generation**: For >5000 distance restraints, limit visualization in global view to violations and a random subset of satisfied restraints (user configurable).
- **Violation table**: Lazy loading; only rows visible in viewport are fully computed (use pagination).
- **Density computation**: O(N_restraints) preprocessing – acceptable for <50k restraints.
- **Coordinate evaluation**: Caching of computed distances/angles; recompute only when structure changes (never in V1).
- **Mol* scene updates**: Regenerate only when state changes; avoid animation loops.

---

## 21. Example Usage Flow

1. User enters PDB ID (`2K39`) in a text field → click "Analyze".
2. System retrieves mmCIF and NMR-STAR, runs atom mapping, computes violations and density.
3. Default global view appears: cartoon colored by density, no restraints shown.
4. User sees a residue with low density (blue) in a loop region.
5. User clicks that residue in sequence view → switches to local view: shows all restraints involving that residue, residues within 8Å.
6. User notices a violated distance restraint (red dashed line) between two side chains.
7. User clicks the violation in the table → view centers, side panel shows details: expected 3.0–5.0Å, measured 6.2Å.
8. User toggles view back to global, adjusts violation magnitude filter to >1Å, sees clustered violations in the same loop.
9. User exports violation table as CSV for further analysis.

---

## 22. Appendix: MolViewSpec State Example (Illustrative)

The following pseudo‑JSON illustrates the declarative state for global view:

```json
{
  "version": "0.1",
  "structures": [{ "url": "https://files.rcsb.org/download/2K39.cif" }],
  "components": [
    {
      "representation": "cartoon",
      "color": {
        "type": "residue_property",
        "property": "restraint_density",
        "palette": ["blue", "white", "red"]
      },
      "visible": true
    },
    {
      "representation": "distance_restraint",
      "data": { "source": "computed", "restraints": [...] },
      "color": { "type": "violation_magnitude", "palette": ["yellow","red"] },
      "style": { "lineWidth": 0.05, "dashed": true }
    }
  ]
}
```

*Note: Actual MolViewSpec syntax may differ; implementation must adapt to the library's official schema.*

---

**End of Specification**