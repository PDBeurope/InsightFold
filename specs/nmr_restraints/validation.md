# NMR Restraints Visualization Validation Plan

Source requirements: `specs/nmr_restraints/requirements.md`

## Validation Scope

Validation must prove that the notebook can be restarted and run top-to-bottom, handles pinned fixtures deterministically, computes documented geometry and density values, degrades gracefully on expected failures, and presents scientifically cautious interpretation.

Implementation is not scientifically complete until fixture expected snapshots are recorded and reviewed by an NMR-aware domain reviewer.

## Validation Matrix

| Check | Method | Pass Criteria |
|---|---|---|
| Top-to-bottom execution | Restart kernel and run all cells on happy-path fixture | No unhandled exceptions; final validation snapshot produced. |
| Local-file execution | Run notebook using cached/local `9L1V` files | Same key outputs as remote fixture, except provenance mode. |
| Input validation | Test invalid PDB IDs, missing local paths, unsupported extensions | Notebook stops before retrieval/parsing with readable messages. |
| Retrieval contract | Inspect provenance and cache behavior for remote fixture | URLs, timestamps, cache key, status, and byte lengths are recorded. |
| Structure parsing | Compare atom/residue counts and selected model metadata to expected snapshot | Exact match to snapshot for pinned fixture and parser version. |
| Restraint parsing | Compare parsed distance, dihedral, unsupported, saveframe, and loop counts | Exact match to snapshot for pinned fixture and parser version. |
| Atom mapping | Compare mapped/unmapped/excluded logical-restraint counts, raw member diagnostics, and coverage | Exact logical numerator/denominator match; `OR` groups with at least one evaluable member count as `1/1` logical restraint mapped; threshold decision correct. |
| Mapping thresholds | Run warning and abort fixtures | 70% to 95% fixture proceeds with partial warning; below 70% fixture skips geometry. |
| Ambiguous `OR` semantics | Run ambiguous fixture and inspect selected member | All candidates measured; smallest distance selected; non-selected members excluded from evaluated output. |
| Distance violation formula | Use fixture rows with known distances and bounds | Both lower and upper bounds are present for evaluated v1 distance restraints; magnitudes match expected values within `1e-6` Angstrom. Rows missing either bound are excluded as unsupported with diagnostics. |
| Dihedral formula | Use wrapped and unwrapped dihedral fixtures | Angles and violation magnitudes match expected values within `1e-5` degrees. |
| Residue density | Compare absolute and normalized density snapshot | Counts exact; normalized values match expected float tolerance. |
| Violation prioritization | Inspect default violation table order | Sorted by descending magnitude, type, then residue number. |
| Empty restraint behavior | Run missing/empty restraint fixture | Clear missing-restraint message; no geometry attempt; no misleading traceback. |
| No-violation behavior | Use fixture or threshold configuration yielding no displayed violations | No-violation message includes caveat that this is not a validation verdict. |
| Visualization state build | Inspect global and local MVS state objects and rendered notebook output | Objects exist and reference expected selections/colors/restraint caps; rendered Mol* views are embedded inline in notebook cells using the `state.molstar_html()` plus base64 `IFrame` pattern. |
| Visualization fallback | Disable or simulate MolViewSpec/Mol* render failure | Tables and summaries still render; warning is visible. |
| Hidden-state hazard | Run clean kernel twice and compare key output hashes/counts | Key tables and states are deterministic for same inputs/config. |
| Dependency/runtime budget | Record dependency versions and elapsed time | Happy-path fixture under 30 seconds after dependencies; exceptions documented with environment. |
| Export artifacts | Run export section | CSV violation export includes provenance; optional outputs listed in export manifest. |
| Documentation completeness | Review notebook against `docs-plan.md` | Tutorial, how-to, reference, and explanation content present at required depth. |
| Scientific language review | Search notebook/export text for overclaims | No density-as-confidence, violation-as-error, pass/fail, or clinical/validation verdict language. |

## Algorithm Validation Details

### Distance Violation

For each evaluated distance restraint:

```text
if measured < lower:
    violation = lower - measured
elif measured > upper:
    violation = measured - upper
else:
    violation = 0
```

Pass criteria:

- Fixture hand checks match within `1e-6` Angstrom.
- Bounds are preserved from parsed restraint source.
- Rows missing lower or upper bounds are excluded as unsupported in v1 and counted in diagnostics.

### Dihedral Violation

Validation must cover:

- four-atom torsion calculation
- normalization to `[-180, +180]`
- unwrapped intervals such as `-60` to `60`
- wrapped intervals such as `170` to `-170`
- nearest-boundary angular distance

Pass criteria:

- Fixture hand checks match within `1e-5` degrees.
- Circular equivalence at `-180`/`+180` is handled without false large violations.

### Ambiguous Distance Restraints

Validation must prove:

1. Rows sharing `_Gen_dist_constraint.ID` and `Member_logic_code = OR` are grouped.
2. All mapped candidate atom-pair distances are computed.
3. The candidate with the smallest measured distance is selected.
4. The selected candidate alone drives violation magnitude, reporting, visualization, and local inspection.
5. Alternative candidates remain traceable in diagnostics if possible, but are not rendered as evaluated members in v1.

### Residue Density

Validation must prove:

```text
absolute_density(residue) = count(unique_logical_restraints_associated_with_residue)
structure_mean_density = N_logical_restraints / N_residues
normalized_density(residue) = absolute_density / structure_mean_density
```

Pass criteria:

- Unique logical restraints, not ambiguous member rows, drive density counts.
- Residues with no associated restraints receive absolute density `0`.
- Normalized density handles zero-restraint cases without divide-by-zero errors.

## Fixture Validation Requirements

| Fixture | Required Checks |
|---|---|
| `happy-path-9l1v` | Retrieval, parser source-tag support, expected file checksums, structure counts, restraint counts, mapping, geometry, density, global view, local view, exports, runtime. |
| `local-copy-9l1v` | Local-file mode parity with remote fixture, provenance mode, cache-disabled behavior, inline embedded MolViewSpec rendering or documented fallback. |
| `edge-or-9l1v` | Ambiguous grouping, logical-restraint mapping coverage, member-level diagnostics, and smallest-distance selected member. |
| `edge-dihedral-9l1v` | Real deposited dihedral parser support and unwrapped torsion-interval evaluation. |
| `negative-missing-restraints-1crn` | Valid model endpoint with missing NMR restraint endpoint; missing-restraint state and graceful termination. |
| `edge-wrapped-dihedral-synthetic` | Circular dihedral interval behavior. |
| `edge-low-mapping-synthetic` | Warning and abort threshold behavior. |

## Snapshot Artifacts

After first validated run, store or document:

- dependency versions
- raw source checksums if files are cached
- structure counts
- restraint parsing counts
- mapping counts and coverage
- top violation rows
- density summary rows
- selected local residue snapshot
- generated MVS state presence and size
- runtime timings

Recommended machine-readable file:

```text
specs/nmr_restraints/fixtures/expected_snapshots.json
```

## Manual Review Checklist

- Density language says "experimental coverage" or equivalent, not "confidence" or "certainty".
- Violation language says "local inconsistency" or equivalent, not "error" or "wrong structure".
- No-violation message says absence of thresholded violations is not a validation verdict.
- Mapping warnings make partial coverage visible before interpretation.
- Ambiguity limitations are stated near ambiguous restraint outputs.
- RUO/non-clinical/non-validation scope appears in notebook introduction and exported summaries.
- Local evidence views are interpretable without needing raw NMR-STAR expertise.

## Validation Blockers

- No pinned happy-path fixture snapshot.
- No edge validation for ambiguous `OR` restraints.
- No wrapped-dihedral validation if dihedral restraints are implemented.
- Mapping coverage thresholds not exercised.
- Notebook depends on manual cell ordering or frontend-only state.
- Visualization failure prevents table outputs.
- Scientific caveats are missing or contradicted by generated summaries.
