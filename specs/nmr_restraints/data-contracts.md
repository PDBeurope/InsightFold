# NMR Restraints Visualization Data Contracts

Source PRD: `specs/nmr_restraints/nmr_restraints_prd.md`

## External Sources

| Source | Access Pattern | Required Fields | Optional Fields | Validation | Failure Behavior |
|---|---|---|---|---|---|
| PDBe model file | `https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif` | mmCIF atom records for selected model; chain, residue, atom names; Cartesian coordinates | altloc, insertion code, model metadata, entity metadata | HTTP success; non-empty mmCIF; Gemmi parse succeeds; selected model exists | Stop structure-dependent workflow; show URL/status/exception |
| PDBe NMR restraint file | `https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_nmr-data.str` | NMR-STAR loops containing compatible distance and/or dihedral restraints | multiple saveframes, multiple loops, ambiguity metadata | HTTP success; non-empty text; parser succeeds; at least one compatible restraint category | Show missing/empty restraint state; skip mapping/geometry |
| RCSB fallback model file | `https://files.rcsb.org/download/{pdb_id}.cif` | Same as model file | same | Used only if configured fallback is enabled and PDBe model retrieval fails | If fallback also fails, stop structure-dependent workflow |
| Local model file | user path ending `.cif` or `.mmcif` | Same as model file | same | Path exists; suffix allowed; Gemmi parse succeeds | Stop before parsing-dependent workflow |
| Local restraint file | user path ending `.str` | Same as restraint file | same | Path exists; suffix allowed; NMR-STAR parse succeeds | Show missing/empty restraint state or parser diagnostic |

## Input Contracts

### Remote Input

| Field | Type | Required | Validation | Notes |
|---|---|---|---|---|
| `pdb_id` | string | yes in remote mode | `^[A-Za-z0-9]{4}$` | Normalize to uppercase for display and lowercase/uppercase as required by endpoint construction. |
| `model_index` | integer | yes | `>= 0`; selected model exists | Default `0`; no ensemble-wide analysis in v1. |
| `violation_threshold_distance` | float | yes | `>= 0` | Default `0.5` Angstrom. |
| `violation_threshold_dihedral` | float | yes | `>= 0` and `<= 180` | Default `5.0` degrees. |
| `local_context_radius` | float | yes | `> 0` | Default `6.0` Angstrom. |
| `max_visible_restraints` | integer | yes | `>= 0` | Default `250`; `0` means hide restraint overlays. |
| `cache_enabled` | boolean | yes | boolean | Default `true`. |

### Local Input

| Field | Type | Required | Validation | Notes |
|---|---|---|---|---|
| `model_path` | path string | yes in local mode | existing `.cif` or `.mmcif` file | Contents parsed as data only. |
| `restraint_path` | path string | yes in local mode | existing `.str` file | Contents parsed as data only. |
| `pdb_id_label` | string | no | optional four-character ID or free text label | Used for provenance and cache key only. |

## Structure Parsing Contract

### `atom_table`

One row per atom in selected model after altloc filtering.

| Column | Type | Required | Description |
|---|---|---|---|
| `model_index` | int | yes | Selected model index. |
| `chain_id` | string | yes | Display chain identifier. |
| `auth_asym_id` | string | yes | Deposited author asym ID. |
| `auth_seq_id` | string/int | yes | Deposited author residue ID; preserve insertion semantics where needed. |
| `auth_comp_id` | string | yes | Deposited residue/component name. |
| `auth_atom_id` | string | yes | Deposited atom name. |
| `label_asym_id` | string | recommended | Label asym ID for MolViewSpec selection. |
| `label_seq_id` | int/null | recommended | Label sequence ID for MolViewSpec residue ranges. |
| `x` | float | yes | X coordinate in Angstroms. |
| `y` | float | yes | Y coordinate in Angstroms. |
| `z` | float | yes | Z coordinate in Angstroms. |
| `element` | string/null | recommended | Chemical element if available. |
| `altloc` | string/null | optional | Alternate-location identifier before filtering. |
| `included_altloc` | bool | yes | Whether atom passed v1 altloc filter. |

### `residue_table`

One row per residue/component present in the selected model after atom filtering.

| Column | Type | Required | Description |
|---|---|---|---|
| `residue_key` | tuple/string | yes | Stable key based on deposited identifiers. |
| `chain_id` | string | yes | Display chain. |
| `auth_asym_id` | string | yes | Deposited chain. |
| `auth_seq_id` | string/int | yes | Deposited residue number. |
| `auth_comp_id` | string | yes | Residue/component name. |
| `label_asym_id` | string | recommended | MolViewSpec selection chain. |
| `label_seq_id` | int/null | recommended | MolViewSpec sequence ID. |
| `ca_x`, `ca_y`, `ca_z` | float/null | recommended | CA coordinates for proteins; fallback centroid for non-protein if implemented. |
| `centroid_x`, `centroid_y`, `centroid_z` | float | yes | Residue/component centroid for fallback spatial queries. |
| `atom_count` | int | yes | Included atom count. |
| `residue_type` | string | optional | protein, nucleic_acid, ligand, water, unknown if classified. |

## Restraint Parsing Contract

### Observed Source Tags For V1 Parser

The parser must support these NMR-STAR tags, observed in the pinned `9L1V` fixture. Any compatible-looking loop missing required tags should be reported in `unsupported_records` with source saveframe/loop context.

| Normalized Field | Required Source Tag |
|---|---|
| distance logical ID | `_Gen_dist_constraint.ID` |
| distance member ID | `_Gen_dist_constraint.Member_ID` |
| distance member logic | `_Gen_dist_constraint.Member_logic_code` |
| distance atom 1 auth chain | `_Gen_dist_constraint.Auth_asym_ID_1` |
| distance atom 1 auth residue | `_Gen_dist_constraint.Auth_seq_ID_1` |
| distance atom 1 auth component | `_Gen_dist_constraint.Auth_comp_ID_1` |
| distance atom 1 auth atom | `_Gen_dist_constraint.Auth_atom_ID_1` |
| distance atom 2 auth chain | `_Gen_dist_constraint.Auth_asym_ID_2` |
| distance atom 2 auth residue | `_Gen_dist_constraint.Auth_seq_ID_2` |
| distance atom 2 auth component | `_Gen_dist_constraint.Auth_comp_ID_2` |
| distance atom 2 auth atom | `_Gen_dist_constraint.Auth_atom_ID_2` |
| distance lower bound | `_Gen_dist_constraint.Distance_lower_bound_val` |
| distance upper bound | `_Gen_dist_constraint.Distance_upper_bound_val` |
| distance source list | `_Gen_dist_constraint.Gen_dist_constraint_list_ID` |
| torsion logical ID | `_Torsion_angle_constraint.ID` |
| torsion name | `_Torsion_angle_constraint.Torsion_angle_name` |
| torsion atom n auth chain | `_Torsion_angle_constraint.Auth_asym_ID_n` |
| torsion atom n auth residue | `_Torsion_angle_constraint.Auth_seq_ID_n` |
| torsion atom n auth component | `_Torsion_angle_constraint.Auth_comp_ID_n` |
| torsion atom n auth atom | `_Torsion_angle_constraint.Auth_atom_ID_n` |
| torsion lower bound | `_Torsion_angle_constraint.Angle_lower_bound_val` |
| torsion upper bound | `_Torsion_angle_constraint.Angle_upper_bound_val` |
| torsion target | `_Torsion_angle_constraint.Angle_target_val` |
| torsion source list | `_Torsion_angle_constraint.Torsion_angle_constraint_list_ID` |

Source-tag requirements:

- V1 distance restraints require `_Gen_dist_constraint.Distance_lower_bound_val` and `_Gen_dist_constraint.Distance_upper_bound_val`.
- V1 dihedral restraints require `_Torsion_angle_constraint.Angle_lower_bound_val` and `_Torsion_angle_constraint.Angle_upper_bound_val`.
- Missing required bound tags or missing row values make the row unsupported for v1.
- `_Gen_dist_constraint.Member_logic_code = OR` identifies ambiguous distance candidate members when rows share the same `_Gen_dist_constraint.ID`.

### `distance_restraints_raw`

One row per parsed distance restraint member before ambiguity collapse.

| Column | Type | Required | Description |
|---|---|---|---|
| `restraint_id` | string/int | yes | Logical `_Gen_dist_constraint.ID` where available. |
| `member_id` | string/int | recommended | Row/member identifier if separate from logical ID. |
| `member_logic_code` | string/null | recommended | Used to detect `OR` ambiguity. |
| `auth_asym_id_1`, `auth_seq_id_1`, `auth_comp_id_1`, `auth_atom_id_1` | mixed | yes | First atom identity. |
| `auth_asym_id_2`, `auth_seq_id_2`, `auth_comp_id_2`, `auth_atom_id_2` | mixed | yes | Second atom identity. |
| `lower_bound` | float | yes | Minimum allowed distance. V1 expects both lower and upper bounds; rows missing either bound are unsupported. |
| `upper_bound` | float | yes | Maximum allowed distance. |
| `source_saveframe` | string/null | recommended | Provenance inside NMR-STAR file. |
| `source_loop` | string/null | recommended | Provenance inside NMR-STAR file. |

### `dihedral_restraints_raw`

One row per parsed dihedral restraint.

| Column | Type | Required | Description |
|---|---|---|---|
| `restraint_id` | string/int | yes | Logical `_Torsion_angle_constraint` ID where available. |
| `auth_asym_id_1..4` | mixed | yes | Four atom chain identifiers. |
| `auth_seq_id_1..4` | mixed | yes | Four atom residue identifiers. |
| `auth_comp_id_1..4` | mixed | yes | Four atom residue/component names. |
| `auth_atom_id_1..4` | mixed | yes | Four atom names. |
| `lower_bound` | float | yes | Minimum allowed angle in degrees. |
| `upper_bound` | float | yes | Maximum allowed angle in degrees. |
| `source_saveframe` | string/null | recommended | Provenance inside NMR-STAR file. |
| `source_loop` | string/null | recommended | Provenance inside NMR-STAR file. |

## Mapping Contract

Canonical atom key:

```text
(auth_asym_id, auth_seq_id, auth_comp_id, auth_atom_id)
```

| Output | Required Fields | Validation | Failure Behavior |
|---|---|---|---|
| `mapped_distance_restraints` | raw fields, logical restraint ID, member ID, atom keys, coordinates for atoms 1 and 2, `mapping_status`, `unmapped_reason` | Each mapped atom key exists in `atom_lookup`; mapping coverage is computed at logical-restraint level. A non-ambiguous logical restraint maps as `1/1` when both atoms map. An `OR` ambiguous logical restraint maps as `1/1` when at least one candidate member maps and can be evaluated. | Unmapped members retained with reason; member-level counts reported separately; low logical-restraint coverage can abort geometry |
| `mapped_dihedral_restraints` | raw fields, atom keys, coordinates for atoms 1 to 4, `mapping_status`, `unmapped_reason` | All four atom keys exist for mapped rows | Unmapped rows retained with reason; low coverage can abort geometry |
| `mapping_report` | parsed logical restraints, mapped logical restraints, unmapped logical restraints, excluded logical restraints, raw member counts, member-level mapped/unmapped counts, coverage, threshold status, examples | `coverage = mapped_logical_restraints / parsed_logical_restraints`; ambiguous `OR` groups count as one logical restraint | Abort below 70%; warn 70% to 95%; normal above 95% |

## Computed Output Contracts

### `evaluated_restraints`

| Column | Type | Required | Description |
|---|---|---|---|
| `restraint_id` | string/int | yes | Logical restraint ID. |
| `restraint_type` | string | yes | `distance` or `dihedral`. |
| `ambiguity_type` | string | yes | `NONE`, `OR`, or unsupported excluded before evaluation. |
| `selected_member` | string/null | yes | Selected atom pair for ambiguous distance restraints. |
| `involved_residues` | list/string | yes | Participating residue keys. |
| `measured_value` | float | yes | Distance in Angstroms or angle in degrees. |
| `lower_bound` | float | yes | Allowed minimum. |
| `upper_bound` | float | yes | Allowed maximum. |
| `violation_magnitude` | float | yes | Distance or angular deviation to nearest allowed boundary. |
| `is_violation` | bool | yes | Whether magnitude exceeds zero. |
| `passes_display_threshold` | bool | yes | Whether magnitude exceeds configured display threshold. |

### `residue_density_table`

| Column | Type | Required | Description |
|---|---|---|---|
| `chain_id` | string | yes | Display chain identifier. |
| `residue_number` | string/int | yes | Residue position. |
| `residue_name` | string | yes | Residue/component name. |
| `residue_key` | tuple/string | yes | Stable residue key. |
| `absolute_density` | int | yes | Count of unique restraints associated with residue. |
| `normalized_density` | float | yes | Absolute density divided by structure-wide mean density. |
| `structure_mean_density` | float | yes | `N_restraints / N_residues`. |
| `local_violation_count` | int | yes | Number of violations involving residue. |

### `violation_table`

Must be sorted by default:

1. descending `violation_magnitude`
2. `restraint_type`
3. first involved residue number

Additional required columns:

- `display_color_class`: `none`, `low`, `medium`, `high`
- `interpretation_note`: short caveat string for local inconsistency, not error
- `source_saveframe`
- `source_loop`

## Visualization State Contract

| State | Required Inputs | Required Behavior | Failure Behavior |
|---|---|---|---|
| `global_mvs_state` | model URL or notebook-embedded/served model content, residue density table, thresholded violation table | Shows whole structure as cartoon colored by density; violated restraints only; cap by `max_visible_restraints`; renders inline in the notebook cell using `state.molstar_html()` encoded into a data-URI `IFrame` | Show density table and warning if state build/render fails |
| `local_mvs_state` | selected residue, local context table, evaluated restraints, model URL or notebook-embedded/served model content | Shows selected residue, neighbors, local restraints, local violations, reduced-opacity context; renders inline in the notebook cell using the same embedded `IFrame` helper | Show local tables and warning if render fails |

## Provenance And Cache Contract

Cache key:

```text
{pdb_id_or_label}_{model_index}_{config_hash}
```

Cached artifacts may include:

- downloaded mmCIF and NMR-STAR files
- parsed `atom_table` and `residue_table`
- parsed restraint tables
- computed `evaluated_restraints`, `violation_table`, and `residue_density_table`
- generated MVS JSON state

Each cached artifact must record:

- source URL or local path
- retrieval/load timestamp
- cache creation timestamp
- config hash
- parser package names and versions where available
