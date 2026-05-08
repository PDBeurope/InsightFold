# Homodimer Diagnostic Data Contracts

## AFDB Prediction Endpoint

Endpoint:

```text
GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}
```

Expected response: JSON array with at least one record.

Required fields:

- `cifUrl`
- `paeDocUrl`
- `plddtDocUrl`
- `sequence`
- `uniprotAccession`

Optional fields:

- `pdbUrl`
- `bcifUrl`
- `paeImageUrl`
- protein names, organism, model metadata, release metadata

Validation checks:

- response is valid JSON
- response is a non-empty list
- selected record contains required fields
- required URL fields use HTTP(S)

Failure behavior:

- missing required fields should stop the notebook before parsing
- optional missing fields should produce warnings only

## PAE JSON

Supported shape:

```json
[
  {
    "predicted_aligned_error": [[...]],
    "max_predicted_aligned_error": 27.47,
    "chains": [
      {"label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 172},
      {"label_asym_id": "B", "sequenceStart": 1, "sequenceEnd": 172}
    ]
  }
]
```

Required fields:

- `predicted_aligned_error`
- `chains`
- for each chain: `label_asym_id`, `sequenceStart`, `sequenceEnd`

Validation checks:

- PAE matrix is square
- matrix dimension equals total chain residue count
- exactly two chains are present for v1
- both chain ranges have equal length
- values are numeric

Failure behavior:

- invalid matrix or chain metadata stops metric computation

## Confidence JSON

Expected shape:

```json
{
  "residueNumber": [1, 2],
  "confidenceScore": [69.38, 76.88],
  "confidenceCategory": ["Low", "Confident"],
  "chains": [...]
}
```

Required fields:

- `residueNumber`
- `confidenceScore`

Optional fields:

- `confidenceCategory`
- `chains`

Validation checks:

- residue and score arrays have equal length
- score count matches total residues or can be mapped to structure residues
- scores are numeric and within expected pLDDT range 0-100

Fallback:

- if confidence JSON is missing or unusable, use mmCIF B-factors as pLDDT with visible provenance note

## mmCIF Atom Site Records

Required atom-level fields:

- atom name
- residue name
- chain label, preferably `label_asym_id`
- residue sequence index
- x/y/z coordinates
- B-factor or equivalent pLDDT field

Required atoms:

- CA for every residue used in distance fallback
- CB for non-glycine residues where available

Validation checks:

- exactly two chains after filtering
- chain residue counts match PAE chain ranges
- chain lengths are equal
- residue coordinates are numeric
- missing CB can fall back to CA only for glycine; otherwise record warning or skip affected residue pair

## Derived Data Contracts

### Residue Table

Required columns:

- `global_index`
- `chain_id`
- `chain_local_index`
- `residue_number`
- `residue_name`
- `ca_x`, `ca_y`, `ca_z`
- `cb_or_ca_x`, `cb_or_ca_y`, `cb_or_ca_z`
- `plddt`

### Contact Pairs

Required columns:

- `chain_a_residue_index`
- `chain_b_residue_index`
- `distance_angstrom`
- `chain_a_residue_number`
- `chain_b_residue_number`

### Score Table

Required rows:

- `ipTM_d0chn`
- `ipSAE_d0res`
- `ipSAE_d0chn`
- `ipSAE_d0dom`
- `pDockQ`
- `pDockQ2`
- `LIS`

Required columns:

- `score`
- `value`
- `directionality`
- `primary_inputs`
- `interpretation_note`

