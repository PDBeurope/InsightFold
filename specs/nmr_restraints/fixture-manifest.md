# NMR Restraints Visualization Fixture Manifest

| Field | Value |
|---|---|
| Notebook / Feature | NMR Restraints Visualization Notebook |
| Source Spec | `specs/nmr_restraints/requirements.md` |
| Owner | InsightFold maintainers / fixture curator TBD |
| Last Updated | 2026-05-08 |

## Fixture Policy

This notebook has no mutation input and no AlphaFold input path. Mutation and AF fixtures are therefore out of scope for this feature. The representative biological examples are PDB/PDBe NMR entries with deposited NMR-STAR restraints, plus negative and synthetic fixtures for expected failure and edge-case behavior.

Network fixtures may be used during development, but regression validation should use repo-local cached copies under `specs/nmr_restraints/fixtures/cache/`. That cache directory is intentionally ignored by git to avoid committing large coordinate/restraint files. Keep checksums and expected-output snapshots in tracked spec files.

## Fixture Summary

| Fixture ID | Role | Source | Stable Identifier | Why Chosen | Runtime Notes |
|---|---|---|---|---|---|
| `happy-path-9l1v` | happy-path / reference | PDBe entry files | PDB `9L1V` | Real solution NMR entry with model coordinates, distance restraints, dihedral restraints, and many `OR` ambiguous distance groups | Expected under 30 seconds after dependencies on local Jupyter CPU |
| `local-copy-9l1v` | local-file parity | Cached PDBe files | PDB `9L1V` local copies | Verifies that local-file mode matches remote mode without network dependency | Uses ignored cache files; same expected outputs as remote fixture |
| `edge-or-9l1v` | edge-case | Same PDBe restraint file as `9L1V` | `_Gen_dist_constraint.ID` groups in `9L1V` | Exercises logical-restraint mapping coverage and smallest-distance selection for ambiguous `OR` restraints | Reuses happy-path data; no extra runtime |
| `edge-dihedral-9l1v` | edge-case | Same PDBe restraint file as `9L1V` | `_Torsion_angle_constraint` loop in `9L1V` | Exercises dihedral parsing and torsion-angle violation workflow on a real deposited file | Reuses happy-path data; no extra runtime |
| `negative-missing-restraints-1crn` | negative | PDBe entry files | PDB `1CRN` | X-ray structure with valid model endpoint and missing PDBe `_nmr-data.str` endpoint; validates graceful non-NMR behavior | Retrieval-only negative fixture; no geometry expected |
| `edge-wrapped-dihedral-synthetic` | edge-case / unit fixture | Synthetic in-memory or tiny local fixture | `wrapped-dihedral-minimal` | Tests circular angle interval logic not covered by `9L1V`, whose torsion bounds are not wrapped | Should run as a fast unit-style validation, not a full notebook story |
| `edge-low-mapping-synthetic` | edge-case / negative | Synthetic edited restraint rows against `9L1V` model or tiny local fixture | `partial-mapping-warning` and `partial-mapping-abort` | Tests 70%-95% warning and below-70% abort thresholds with controlled unmapped atoms | Should run as parser/mapping validation, not user-facing story |

## Fixture Details

### `happy-path-9l1v`

| Field | Value |
|---|---|
| Role | happy-path / reference |
| Source Database / File | PDBe entry-file service |
| Stable Identifier | PDB `9L1V` |
| Retrieval Endpoint / Command | `curl -L -o specs/nmr_restraints/fixtures/cache/9l1v.cif https://www.ebi.ac.uk/pdbe/entry-files/download/9l1v.cif`; `curl -L -o specs/nmr_restraints/fixtures/cache/9l1v_nmr-data.str https://www.ebi.ac.uk/pdbe/entry-files/download/9l1v_nmr-data.str` |
| Retrieval Date | 2026-05-08 |
| Cached Path | `specs/nmr_restraints/fixtures/cache/9l1v.cif`; `specs/nmr_restraints/fixtures/cache/9l1v_nmr-data.str` |
| License / Use Constraints | Public PDB/PDBe data; preserve accession and source provenance |
| Data Format / Schema Version | mmCIF model file; NMR-STAR restraint file from PDBe entry files |
| Expected Runtime | Under 30 seconds after dependencies are available |

#### Expected Outputs

| Output | Expected Value | Comparison Mode | Tolerance / Notes |
|---|---|---|---|
| Model endpoint status | `200` | exact | Verified 2026-05-08 |
| Restraint endpoint status | `200` | exact | Verified 2026-05-08 |
| Model file byte length | `2362171` | exact | Repo-local cached copy |
| Restraint file byte length | `413575` | exact | Repo-local cached copy; PDBe HEAD returned same content length |
| Model file SHA-256 | `e2ec15597e973bed4fa2673a53bf6756b9ee9bac805988020bdf881f25b47e74` | exact | Repo-local cached copy |
| Restraint file SHA-256 | `b1c1535358df907f732e9946b8148ab9d4f0a130cd445f1ff3f1bda1cac6865b` | exact | Repo-local cached copy |
| Experimental method | `SOLUTION NMR` | exact | `_exptl.method` in mmCIF |
| Submitted NMR models | `20` | exact | `_pdbx_nmr_ensemble.conformers_submitted_total_number`; also 20 atom-site model IDs |
| Default selected model | model index `0`, corresponding to PDB model number `1` | exact | Notebook default |
| Model 1 atom count | `1330` | exact | Count of `ATOM` rows with `_atom_site.pdbx_PDB_model_num = 1` |
| Model 1 residue count | `82` | exact | Unique `(auth_asym_id, auth_seq_id, auth_comp_id)` in model 1 |
| Model 1 chains | `A` | exact | Single chain |
| Distance restraint loops | `2` loops | exact | `_Gen_dist_constraint` loops with 35 headers |
| Distance restraint member rows | `2219` | exact | Rows across both `_Gen_dist_constraint` loops |
| Logical distance restraints | `1150` | exact | Unique `_Gen_dist_constraint.ID` values |
| `OR` ambiguous logical distance restraints | `451` | exact | Unique IDs with `_Gen_dist_constraint.Member_logic_code = OR` |
| Torsion restraint loops | `1` loop | exact | `_Torsion_angle_constraint` loop with 53 headers |
| Torsion restraint rows | `149` | exact | Rows in `_Torsion_angle_constraint` loop |
| Logical torsion restraints | `149` | exact | Unique `_Torsion_angle_constraint.ID` values |
| Distance rows missing lower/upper bound | `0` | exact | Confirms v1 assumption that both bounds are present for this fixture |
| Torsion rows missing lower/upper bound | `0` | exact | Confirms bounds are present |
| Wrapped torsion bounds | `0` | exact | Wrapped interval behavior requires synthetic fixture |
| Mapping coverage | `TBD after implementation parser run` | manual-review | Must be recorded before validation sign-off |
| Top 5 violations | `TBD after implementation parser run` | manual-review | Must be recorded before validation sign-off |
| Density summary | `TBD after implementation parser run` | manual-review | Mean density, max-density residue, selected local-view residue |
| Global MolViewSpec render | Inline embedded IFrame exists | presence | Uses `state.molstar_html()` base64 IFrame pattern |

#### Validation Checks

- [ ] Remote retrieval succeeds for both model and restraint endpoints.
- [ ] Local cached files match recorded byte lengths and SHA-256 checksums.
- [ ] Parser reports 20 submitted models and selects model number 1 by default.
- [ ] Parser produces 1,330 atoms and 82 residues for model 1.
- [ ] NMR-STAR parser finds 2 distance loops, 2,219 distance member rows, 1,150 logical distance restraints, and 451 `OR` logical groups.
- [ ] NMR-STAR parser finds 1 torsion loop and 149 torsion restraints.
- [ ] Mapping report records logical-restraint coverage and separate member-level diagnostics.
- [ ] Global density view and local evidence view render inline or produce documented fallback.

#### Known Risks

- Mapping coverage and violation outputs are not frozen until the implementation parser and geometry code run.
- `9L1V` is suitable for `OR` ambiguity and dihedral coverage, but not wrapped-dihedral interval validation.

### `local-copy-9l1v`

| Field | Value |
|---|---|
| Role | local-file parity |
| Source Database / File | Repo-local ignored fixture cache |
| Stable Identifier | PDB `9L1V` |
| Retrieval Endpoint / Command | Same as `happy-path-9l1v`; files loaded from `specs/nmr_restraints/fixtures/cache/` |
| Retrieval Date | 2026-05-08 |
| Cached Path | `specs/nmr_restraints/fixtures/cache/9l1v.cif`; `specs/nmr_restraints/fixtures/cache/9l1v_nmr-data.str` |
| License / Use Constraints | Public PDB/PDBe data; cache ignored by git |
| Data Format / Schema Version | mmCIF and NMR-STAR |
| Expected Runtime | Under 30 seconds after dependencies are available |

#### Expected Outputs

| Output | Expected Value | Comparison Mode | Tolerance / Notes |
|---|---|---|---|
| Parsed model/restraint counts | Same as `happy-path-9l1v` | exact | Must match remote fixture for same cached files |
| Provenance mode | `local-file` | exact | Must display local paths and checksums |
| Network requirement | none | exact | Validation should pass with network disabled after files are cached |
| Inline MolViewSpec render | embedded notebook cell view or documented fallback | presence | Local mode must not rely on arbitrary browser access to filesystem paths |

#### Validation Checks

- [ ] Local file inputs load without network access.
- [ ] Parsed counts match `happy-path-9l1v`.
- [ ] Provenance identifies local paths and SHA-256 checksums.
- [ ] Tables are identical to remote mode for the same selected model and config.

#### Known Risks

- Live Mol* rendering from local model content requires embedding or notebook-accessible serving. The notebook must use the embedded IFrame pattern or provide a clear fallback.

### `edge-or-9l1v`

| Field | Value |
|---|---|
| Role | edge-case |
| Source Database / File | Same cached NMR-STAR file as `happy-path-9l1v` |
| Stable Identifier | PDB `9L1V`, `_Gen_dist_constraint` logical IDs with `Member_logic_code = OR` |
| Retrieval Endpoint / Command | Same as `happy-path-9l1v` |
| Retrieval Date | 2026-05-08 |
| Cached Path | `specs/nmr_restraints/fixtures/cache/9l1v_nmr-data.str` |
| License / Use Constraints | Public PDB/PDBe data |
| Data Format / Schema Version | NMR-STAR `_Gen_dist_constraint` |
| Expected Runtime | No additional runtime beyond happy path |

#### Expected Outputs

| Output | Expected Value | Comparison Mode | Tolerance / Notes |
|---|---|---|---|
| `OR` logical group count | `451` | exact | Unique `_Gen_dist_constraint.ID` with `Member_logic_code = OR` |
| First ambiguous logical ID | `1` | exact | Three members: MET 1 HA to THR 80 HG21/HG22/HG23 |
| Logical mapping rule | `1/1` mapped if at least one candidate member maps and can be evaluated | exact | Member-level mapping diagnostics still reported |
| Selected member | smallest measured candidate distance | exact after geometry | Must be recorded after implementation run |

#### Validation Checks

- [ ] Parser groups member rows by logical `_Gen_dist_constraint.ID`.
- [ ] Coverage counts the ambiguous group as one logical restraint.
- [ ] Candidate distances are all evaluated when atoms map.
- [ ] Selected member is the smallest measured candidate and is the only member shown in evaluated restraint output for v1.

#### Known Risks

- Some ambiguous members may fail exact atom mapping if atom naming differs between restraint and coordinate file. This should affect member diagnostics, not logical coverage if at least one member remains evaluable.

### `edge-dihedral-9l1v`

| Field | Value |
|---|---|
| Role | edge-case |
| Source Database / File | Same cached NMR-STAR file as `happy-path-9l1v` |
| Stable Identifier | PDB `9L1V`, `_Torsion_angle_constraint` loop |
| Retrieval Endpoint / Command | Same as `happy-path-9l1v` |
| Retrieval Date | 2026-05-08 |
| Cached Path | `specs/nmr_restraints/fixtures/cache/9l1v_nmr-data.str` |
| License / Use Constraints | Public PDB/PDBe data |
| Data Format / Schema Version | NMR-STAR `_Torsion_angle_constraint` |
| Expected Runtime | No additional runtime beyond happy path |

#### Expected Outputs

| Output | Expected Value | Comparison Mode | Tolerance / Notes |
|---|---|---|---|
| Torsion loop count | `1` | exact | 53 headers |
| Torsion restraint rows | `149` | exact | Also 149 logical IDs |
| First torsion restraint | PHI for MET 1 C / LEU 2 N / LEU 2 CA / LEU 2 C | exact | Bounds `-139.4` to `-63.18`, target `-101.29` |
| Wrapped intervals | `0` in this fixture | exact | Use synthetic wrapped fixture for circular-boundary validation |

#### Validation Checks

- [ ] Parser extracts all four atom identities and angle bounds.
- [ ] Geometry engine computes torsion values from mapped atoms.
- [ ] Violation engine evaluates unwrapped intervals.
- [ ] Wrapped interval unit test remains separate.

#### Known Risks

- This fixture validates real dihedral parsing but not wrapped-boundary behavior.

### `negative-missing-restraints-1crn`

| Field | Value |
|---|---|
| Role | negative |
| Source Database / File | PDBe entry-file service |
| Stable Identifier | PDB `1CRN` |
| Retrieval Endpoint / Command | Model: `https://www.ebi.ac.uk/pdbe/entry-files/download/1crn.cif`; restraint: `https://www.ebi.ac.uk/pdbe/entry-files/download/1crn_nmr-data.str` |
| Retrieval Date | 2026-05-08 endpoint check |
| Cached Path | none required |
| License / Use Constraints | Public PDB/PDBe data |
| Data Format / Schema Version | mmCIF model available; NMR-STAR restraint endpoint absent |
| Expected Runtime | Fast retrieval failure path |

#### Expected Outputs

| Output | Expected Value | Comparison Mode | Tolerance / Notes |
|---|---|---|---|
| Model endpoint status | `200` | exact | Verified with HEAD request 2026-05-08 |
| Model endpoint content length | `69506` | exact/range | HEAD value observed 2026-05-08; may change if PDBe updates file representation |
| Restraint endpoint status | `404` | exact | Verified with HEAD request 2026-05-08 |
| Notebook behavior | missing-restraint state | exact | No mapping or geometry attempted |
| User-facing message | says no compatible NMR restraint data are available | presence | Must not imply structure error |

#### Validation Checks

- [ ] Notebook retrieves or recognizes model availability.
- [ ] Restraint retrieval failure is caught and displayed.
- [ ] Geometry and visualization sections depending on restraints are skipped.
- [ ] No uncaught traceback appears.

#### Known Risks

- This is intentionally not an NMR science fixture. It validates non-NMR user input behavior.

### `edge-wrapped-dihedral-synthetic`

| Field | Value |
|---|---|
| Role | edge-case / unit fixture |
| Source Database / File | Synthetic fixture to be created during implementation |
| Stable Identifier | `wrapped-dihedral-minimal` |
| Retrieval Endpoint / Command | n/a |
| Retrieval Date | n/a |
| Cached Path | `specs/nmr_restraints/fixtures/synthetic/wrapped_dihedral_minimal.*` once created |
| License / Use Constraints | Repo-authored synthetic test data |
| Data Format / Schema Version | Minimal atom table or mmCIF plus NMR-STAR torsion row |
| Expected Runtime | Sub-second unit-style check |

#### Expected Outputs

| Output | Expected Value | Comparison Mode | Tolerance / Notes |
|---|---|---|---|
| Wrapped interval example | lower `170`, upper `-170` | exact | Interval crosses +/-180 degrees |
| Measured angle at `180` or `-180` | violation `0` | tolerance | `1e-5` degrees |
| Measured angle at `0` | violation `170` | tolerance | `1e-5` degrees |
| Notebook behavior | circular angular geometry used | exact | No false 340-degree violation |

#### Validation Checks

- [ ] Circular interval helper handles wrapped interval.
- [ ] Boundary equivalence at `-180/+180` is accepted.
- [ ] Fixture runs independently of external endpoints.

#### Known Risks

- This fixture is synthetic and validates algorithm correctness, not deposition complexity.

### `edge-low-mapping-synthetic`

| Field | Value |
|---|---|
| Role | edge-case / negative |
| Source Database / File | Synthetic edited restraint rows against `9L1V` or minimal local fixture |
| Stable Identifier | `partial-mapping-warning`; `partial-mapping-abort` |
| Retrieval Endpoint / Command | n/a |
| Retrieval Date | n/a |
| Cached Path | `specs/nmr_restraints/fixtures/synthetic/partial_mapping_*.str` once created |
| License / Use Constraints | Repo-authored synthetic test data |
| Data Format / Schema Version | NMR-STAR-like restraint rows with controlled invalid atom IDs |
| Expected Runtime | Sub-second parser/mapping check |

#### Expected Outputs

| Output | Expected Value | Comparison Mode | Tolerance / Notes |
|---|---|---|---|
| Warning fixture coverage | `>= 0.70` and `< 0.95` logical coverage | exact/range | Analysis proceeds with partial warning |
| Abort fixture coverage | `< 0.70` logical coverage | exact/range | Geometry skipped |
| Unmapped examples | at least one invalid atom identifier shown | presence | Must include logical restraint ID and atom key |

#### Validation Checks

- [ ] Logical-restraint coverage threshold is applied.
- [ ] Member-level counts are separate from logical coverage.
- [ ] Warning path and abort path are both exercised.

#### Known Risks

- Synthetic mapping failures should not be mistaken for real deposition behavior.

## Source Tag Mapping From `9L1V`

These observed tags should be added to parser/data-contract implementation references.

| Normalized Field | Observed NMR-STAR Tag |
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

## Deferred Fixture Candidates

| Candidate | Reason Deferred | Condition to Revisit |
|---|---|---|
| `2K39` | PDBe `_nmr-data.str` endpoint returned `404` during endpoint check, so it is not suitable as a required PDBe-entry-file fixture | Revisit only if another source supplies a stable NMR-STAR restraint file and the notebook adds that source |
| Additional real wrapped-dihedral entry | Not needed for v1 because synthetic wrapped fixture gives deterministic circular-geometry validation | Revisit if domain reviewer wants all scientific edge cases represented by deposited entries |
| Large NMR restraint entry near performance boundary | Current v1 needs a fast routine validation set first | Add after core workflow is implemented and runtime is understood |

## Fixture Selection Rationale

- `9L1V` is the required primary fixture because it exercises the notebook's main claim: deposited NMR coordinates plus real NMR-STAR distance and dihedral restraints can be parsed, mapped, summarized, and visualized.
- `9L1V` also covers ambiguous `OR` distance restraints at meaningful scale, with 451 logical `OR` groups.
- `1CRN` is intentionally an X-ray negative fixture. It tests a common user error: entering a valid PDB ID that has no NMR restraint file.
- Synthetic edge fixtures are appropriate for wrapped dihedral and mapping-threshold behavior because they require deterministic edge conditions that are hard to guarantee in deposited entries.
