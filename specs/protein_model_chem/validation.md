# Protein Model Chemistry Validation Plan

Source PRD: `specs/protein_model_chem/protein_model_chem_prd.md`

## Validation Status

Current status: `fixture-selected`.

Pinned fixtures are available in `fixture-manifest.md`. Exact local interaction counts, perturbation counts, graph deltas, and severity labels must be frozen after the first trusted implementation run and domain review.

## Validation Matrix

| Check | Method | Pass Criteria |
|---|---|---|
| Top-to-bottom execution | Restart kernel and run all cells on FX-001 | Notebook completes without errors in target runtime and produces all required tables and figures. |
| Hidden-state safety | Run cells only in notebook order after restart | No section depends on variables created by manual out-of-order execution. |
| Dependency audit | Inspect imports and install cells | No GPU-only, MD, heavy force-field, PyTorch, or unapproved parser dependency appears. |
| Configuration audit | Inspect displayed `CONFIG` and implementation constants | Source priority, 8 Angstrom mutation radius, interaction thresholds, AF confidence thresholds, and baseline weights match `notebook-design.md`. |
| Structure retrieval contract | Run AFDB and PDBe retrieval fixtures | Coordinates, metadata, provenance, and confidence/quality fields are present or clearly warned. |
| Structure parsing contract | Compare parsed residue and atom tables against fixture expectations | Chain IDs, residue counts, mutation site, missing residues, and coordinates match snapshots. |
| Mutation validation | Run valid and mismatch mutation fixtures | Valid mutation passes; mismatch stops with observed/requested residue message. |
| Interaction detection | Compare interaction table against fixture snapshots | Required interaction classes are present where expected; thresholds and evidence are recorded. |
| Interaction weighting | Compare interaction records against config and factor bounds | `geometry_confidence`, `burial_factor`, `confidence_factor`, `cooperativity_factor`, and `relative_weight` are present, bounded, and labeled heuristic/non-thermodynamic. |
| Graph construction | Check graph edge/node consistency against interaction table | Every graph edge maps to an interaction row and every mutation-site node is present. |
| Confidence propagation | Run low-confidence fixture | Structure, interaction, and mechanistic confidence are downgraded near low-quality regions. |
| Mutation perturbation output | Run FX-001 after fixture pinning | WT/mutant comparison reports reviewed perturbation signals without exact mutant structure claims. |
| Unsupported-context warnings | Run missing-context fixture when available | Ligand/cofactor/PTM/metal/membrane caveats appear and strong claims are suppressed. |
| Visualization existence | Inspect rendered notebook or generated artifacts | Structure overview, confidence overlay, local neighborhood, interaction summary, and perturbation view render or documented fallback appears. |
| MolViewSpec safety | Run with `molviewspec` available and unavailable when practical | MolViewSpec API is verified, views render or skip with clear messages, no required downstream cell depends on MolViewSpec-only variables, and Mol* viewer fallback or `.mvsj` export is emitted when inline rendering is unavailable. |
| Explanation guardrails | Scan generated summaries | Every statement has observation/interpretation/hypothesis label; no clinical/pathogenicity/exact-ddG/deterministic terms are used. |
| Export completeness | Inspect export bundle | Tables, figures, provenance, caveats, and summaries are present and traceable to input fixture. |
| Export formats | Inspect generated artifacts | JSON summary, CSV tables, PNG or SVG figures, and graph edge/node data or graph object export are present unless a documented fallback is triggered. |
| Documentation completeness | Review notebook and docs plan | Tutorial, how-to, reference, and explanation needs are covered. |

## Fixture-Specific Validation

### FX-001: Happy Path

Pinned input:

- Source: PDB/PDBe `2LZM`
- Mutation: `A:L99A`
- Expected source metadata: X-ray diffraction, 1.70 Angstrom resolution, chain `A`, sequence length 164.

Expected validation:

- retrieval succeeds
- mutation validates
- residue table includes mutation site
- WT interaction table includes local interactions near mutation site
- graph metrics are produced
- mutant compatibility delta includes at least one reviewed signal
- summary includes observation, interpretation, hypothesis, confidence, and structural disruption likelihood
- runtime is below target

Snapshot fields to freeze:

- structure source and accession
- mutation notation
- residue count
- mutation-site confidence
- local WT interaction count
- interaction types near mutation site
- perturbation signal counts
- severity label and rationale

### FX-002: AFDB Confidence Caveat

Pinned input:

- Source: AFDB `AF-P04637-F1` / UniProt `P04637`
- Mutation: `A:M1A`
- Expected source metadata: model version 6, sequence length 393, global metric 75.06, `fractionPlddtVeryLow` 0.298, `fractionPlddtLow` 0.104.

Expected validation:

- pLDDT or local confidence is parsed
- low-confidence region is detected near mutation or selected residue
- interaction and mechanistic confidence are downgraded
- strong mechanism and severity claims are suppressed

### FX-003: Interface Perturbation

Pinned input:

- Source: PDB/PDBe `6M0J`
- Mutation: `B:K417N` using label chain `B` / auth chain `E`
- Expected source metadata: X-ray diffraction, 2.45 Angstrom resolution, ACE2 chain `A`, Spike RBD label chain `B` / auth chain `E`.

Expected validation:

- interface or inter-chain contacts are detected
- assembly assumptions are shown
- interface disruption language appears only when supported by contact evidence

### FX-004: Missing Biological Context

Pinned input:

- Source: PDB/PDBe `6M0J`
- Mutation: `A:H374A`
- Expected context metadata: NAG, ZN, CL, membrane-associated Spike entity, interface context.

Expected validation:

- ligand/cofactor/PTM/metal/membrane caveat appears when metadata or structure indicates missing context
- summary uses cautious language and avoids complete-mechanism claims

### FX-005: Invalid Identifier

Pinned input:

- Source: PDB/PDBe `XXXX`
- Mutation: `A:L99A`

Expected validation:

- user-facing error identifies invalid/unavailable input
- downstream sections are skipped or gated
- no undefined-variable traceback is the primary output

### FX-006: Mutation Mismatch

Pinned input:

- Source: PDB/PDBe `2LZM`
- Mutation: `A:W99A`
- Expected observed residue: chain `A` residue 99 is leucine.

Expected validation:

- user-facing error reports observed residue and requested wild-type residue
- mutation perturbation and severity sections do not run

### FX-007: Missing Residue

Pinned input:

- Source: PDB/PDBe `6M0J`
- Mutation: `A:M1A`
- Status: optional until parser chain/residue mapping is confirmed.

Expected validation:

- missing or unavailable residue coordinate message is displayed
- available residue range is shown where parser can report it
- mutation perturbation and severity sections do not run

## Scientific Validation

Required before scientific sign-off:

- Domain reviewer approves interaction thresholds and definitions.
- Domain reviewer approves mutation model assumptions and examples.
- Domain reviewer reviews graph metrics for interpretability and false-authority risk.
- Domain reviewer approves severity bands or confirms they must remain provisional.
- At least one fixture has expert-reviewed expected observations and interpretation.

No exact ddG, total stability score, pathogenicity score, or clinical outcome validation is expected or allowed for v1.

## Data Contract Validation

Validate that each derived artifact follows `data-contracts.md`:

- `structure_request`
- AFDB or PDBe metadata
- `atom_table`
- `residue_table`
- mutation input and validation result
- `wt_interactions`
- interaction graph
- `mutant_compatibility_delta`
- `mechanistic_summary`
- export bundle

## Runtime And Dependency Budget

Target:

- FX-001 under 2 minutes after dependencies are installed.
- Structure analysis under 30 seconds for an average monomer fixture.
- Local mutation perturbation recomputation under 10 seconds after structure normalization.
- No external services beyond AFDB/PDBe retrieval unless fixture cache is used.
- No GPU requirement.
- No MD engine, force-field engine, or heavy ML dependency.

Record during validation:

- runtime
- package versions
- MolViewSpec version and accepted parse format, when installed
- platform
- structure-analysis runtime
- local mutation recomputation runtime
- whether network or cached fixture mode was used

## Documentation Validation

Notebook must include:

- clear run instructions
- input format examples
- interaction definitions and thresholds
- mutation model assumptions
- confidence model explanation
- RUO and non-clinical limitations
- unsupported-context caveats
- interpretation vocabulary explaining observation, interpretation, and hypothesis

## Failure Conditions

Validation fails if:

- implementation invents hidden fixture values not recorded in the fixture manifest
- generated summaries make clinical, pathogenicity, exact-energy, or deterministic folding claims
- low-confidence structure regions still produce high-confidence mechanism claims
- graph metrics are presented without definitions
- mutation analysis runs after a residue mismatch
- output cannot be traced back to structure, interaction, mutation, and confidence evidence
