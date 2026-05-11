# Protein Model Chemistry Data Contracts

Source PRD: `specs/protein_model_chem/protein_model_chem_prd.md`

## Structure Request

Required fields:

- `source`: one of `afdb`, `pdbe`, or `local`
- `identifier`: AFDB accession, UniProt accession, PDBe/PDB ID, or local path

Optional fields:

- `source_preference`: `pdbe_first`, `afdb_first`, or `manual`
- `chain_filter`
- `assembly_policy`: `asymmetric_unit`, `biological_assembly`, or `manual`
- `format_hint`: `mmcif` or `pdb`

Validation checks:

- source is supported
- identifier is non-empty and compatible with source
- local path exists if local files are enabled
- source preference does not silently override an explicit source

Failure behavior:

- invalid source or identifier stops retrieval
- ambiguous source preference produces a clear prompt-like message in notebook output, not hidden fallback behavior

## AlphaFold / AFDB Contract

Expected retrieval:

```text
GET https://alphafold.ebi.ac.uk/api/prediction/{accession_or_uniprot}
```

Required fields from selected metadata record:

- structure URL, preferably mmCIF
- accession or entry ID
- sequence or sequence length when available
- confidence source, such as pLDDT JSON or pLDDT in structure B-factors

Optional fields:

- model version
- organism and protein names
- PAE URLs if available
- global confidence fields
- complex or chain metadata

Validation checks:

- response is valid JSON and non-empty
- selected record contains a usable coordinate URL
- coordinate URL uses HTTP(S)
- pLDDT can be mapped to residues or a visible missing-confidence warning is emitted

Failure behavior:

- missing coordinates stops analysis
- missing confidence metadata permits analysis with downgraded confidence only
- low pLDDT near the mutation site downgrades mechanistic confidence and warns the user

## PDBe Contract

Expected retrieval:

```text
PDBe/PDB identifier -> mmCIF coordinate file + PDBe metadata endpoints
```

Required fields:

- coordinate file or coordinate text
- PDB ID
- chain IDs
- residue identifiers
- atom coordinates

Optional fields:

- experimental method
- resolution
- R-free / R-work where available
- biological assembly metadata
- missing residues
- ligands, metals, cofactors, PTMs, engineered construct notes

Validation checks:

- coordinate file is parseable
- requested chain exists
- requested mutation residue maps to one residue
- experimental quality metadata is displayed when available
- missing residues and alternate conformations are recorded

Failure behavior:

- missing coordinates or requested chain stops mutation analysis
- missing experimental metadata does not stop analysis but lowers provenance completeness
- ligand/cofactor/PTM/metal context near mutation triggers unsupported-context warning

## Local Structure Contract

Status: optional v1.

Supported formats:

- mmCIF, pending parser decision
- PDB, pending parser decision

Required fields:

- atom coordinates
- residue names
- chain IDs
- residue sequence identifiers

Optional fields:

- B-factors as confidence proxy
- REMARK records
- biological assembly records

Validation checks:

- file exists and is below notebook-friendly size
- parser supports the format
- residue table can be constructed

Failure behavior:

- unsupported local format stops analysis
- missing confidence metadata continues with explicit low-provenance note

## Mutation Input Contract

Required fields:

- `chain_id`
- `residue_id` or `residue_number`
- `wild_type_residue`
- `mutant_residue`

Optional fields:

- insertion code
- alternate conformer selection
- residue label for display
- mutation notation string such as `R273H`, `A:W128A`, or `CHAIN_A:W128A`

Parsed mutation notation fields:

- `raw_notation`
- `parsed_chain_id`
- `parsed_wild_type_residue`
- `parsed_residue_number`
- `parsed_mutant_residue`
- `parse_status`
- `parse_error`

Validation checks:

- chain exists
- residue maps uniquely
- observed residue matches requested wild-type residue
- mutant residue is a standard amino acid supported by v1
- mutation site has enough coordinates for local analysis

Failure behavior:

- mismatched wild-type residue stops mutation analysis
- unsupported mutant residue stops mutation analysis
- missing local coordinates stops mutation perturbation but preserves source/WT inspection outputs

## Derived Residue Table

Required columns:

- `residue_uid`
- `source`
- `chain_id`
- `residue_number`
- `insertion_code`
- `residue_name`
- `one_letter_code`
- `has_backbone`
- `has_sidechain`
- `ca_x`, `ca_y`, `ca_z`
- representative sidechain or centroid coordinates, if available
- `confidence_value`
- `confidence_source`
- `quality_flags`

Validation checks:

- residue identifiers are unique after normalization
- coordinates are numeric
- confidence values are in expected range where present
- missing coordinates are explicit flags, not null surprises downstream

## Interaction Table

Required columns:

- `interaction_id`
- `residue_uid_a`
- `residue_uid_b`
- `chain_id_a`
- `chain_id_b`
- `residue_label_a`
- `residue_label_b`
- `interaction_type`
- `geometry_metric`
- `threshold`
- `passes_threshold`
- `geometry_confidence`
- `burial_factor`
- `confidence_factor`
- `cooperativity_factor`
- `relative_weight`
- `structure_confidence`
- `interaction_confidence`
- `evidence_note`

Optional columns:

- `directionality`
- `atom_ids`
- `distance_angstrom`
- `angle_degrees`
- `burial_proxy`

Validation checks:

- every residue UID exists in residue table
- interaction type is one of the configured supported types
- threshold value is recorded
- confidence fields are present
- relative weights are present and are labeled as heuristic/non-thermodynamic

## Interaction Graph Contract

Required representation:

- NetworkX graph or edge/node tables with equivalent information.

Required node attributes:

- `residue_uid`
- chain and residue labels
- residue name
- confidence summary
- mutation-site flag

Required edge attributes:

- interaction type
- interaction confidence
- `geometry_confidence`
- `burial_factor`
- `confidence_factor`
- `cooperativity_factor`
- `relative_weight`
- evidence note

Validation checks:

- graph edges match interaction table rows
- graph metric outputs record formulas or library functions used
- multi-edge behavior preserves multiple interaction classes between the same residue pair

## Mutant Compatibility Delta Contract

Required columns:

- `delta_id`
- `residue_uid`
- `mutation`
- `signal_type`: `lost_interaction`, `weakened_interaction`, `new_clash`, `packing_change`, `cavity_signal`, `network_change`, or `confidence_warning`
- `wt_evidence`
- `mutant_estimate`
- `lost_interactions`
- `gained_interactions`
- `gained_clashes`
- `packing_change`
- `centrality_delta`
- `severity`
- `observation_text`
- `interpretation_text`
- `confidence`

Validation checks:

- every signal is traceable to WT interaction, geometry, mutation rule, or confidence rule
- exact mutant conformations and exact energies are not emitted
- confidence is downgraded for missing atoms or low-quality regions
- mutation perturbation is limited to sidechain substitution and local recomputation within the configured 8 Angstrom radius

## Mechanistic Summary Contract

Required fields:

- `layer`: `observation`, `interpretation`, or `hypothesis`
- `statement`
- `supporting_records`
- `confidence`
- `caveat`

Optional grouped interpretation record fields:

- `observation`
- `interpretation`
- `hypothesis`

Validation checks:

- observation statements describe computed facts only
- interpretation statements use cautious language
- hypothesis statements are explicitly marked speculative
- clinical/pathogenicity terms are absent

## Export Bundle Contract

Expected artifacts:

- `residue_table`
- `wt_interactions`
- `mutant_compatibility_delta`
- `network_metrics`
- `mechanistic_summary`
- `severity_summary`
- JSON summary
- CSV tables
- PNG or SVG figures
- NetworkX/graph object export or serialized edge/node tables
- static figures or notebook-rendered figures
- provenance and caveats

Validation checks:

- exported summaries preserve RUO wording
- filenames or in-notebook artifacts include fixture/accession/mutation identifiers where safe
- no hidden dependency on notebook execution order
