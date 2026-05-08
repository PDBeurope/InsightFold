---
name: homodimer-confidence-scoring
description: Implement and review AlphaFold Database homodimer interface confidence scoring. Use when Codex needs AFDB homodimer PAE quadrant extraction, ipTM, ipSAE variants, pDockQ, pDockQ2, LIS, traffic-light thresholds, API contracts, mmCIF parsing rules, interface detection, or validation against reference scoring behavior.
---

# Homodimer Confidence Scoring

## Purpose

Compute interface quality scores from AlphaFold Database homodimer predictions with formulas and data contracts preserved. Use this skill with `$afdb-api-fetching`, `$ipsae-implementation`, `$pdockq-implementation`, `$lis-implementation`, `$pae-visualization`, and `$validation-against-reference` as needed.

## Inputs

Use:

- AFDB accession or already fetched AFDB metadata
- PAE JSON
- pLDDT JSON
- mmCIF or parsed chain coordinates
- chain lengths and chain IDs
- fixture expected outputs or a reference implementation when available

## Outputs

Produce:

- parsed PAE quadrants: `pae_AB`, `pae_BA`, `pae_AA`, `pae_BB`
- interface masks and contact counts
- score dictionaries for ipTM, ipSAE variants, pDockQ, pDockQ2, and LIS
- intermediate values needed for visualization and validation
- validation notes with tolerance, edge cases, and unresolved assumptions

## Core Rules

- Preserve PAE asymmetry. Never symmetrize the PAE matrix unless the user explicitly asks for a separate exploratory view.
- Sort PAE `chains` by `label_asym_id` before deriving chain lengths when chain order is uncertain.
- Use CB-CB contact distance <= 8.0 Angstrom; fall back to CA for glycine.
- Keep formulas inline and readable in notebooks.
- Return intermediate values in score dictionaries so later notebook sections do not recompute hidden state.
- Validate final numeric scores against reference behavior within `+/- 0.001` when a reference is available.

## Primitive Functions

Use these primitives consistently:

```python
def d0_func(L: int) -> float:
    """TM-score length-dependent normalisation."""
    if L <= 27:
        return 1.0
    return max(1.0, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)


def ptm_func(pae_vals, d0: float):
    """Vectorised TM-score transform."""
    return 1.0 / (1.0 + (pae_vals / d0) ** 2)
```

## PAE Quadrants

Extract quadrants from AFDB-style PAE JSON:

```python
chains = sorted(pae_json[0]["chains"], key=lambda c: c["label_asym_id"])
nA = chains[0]["sequenceEnd"] - chains[0]["sequenceStart"] + 1
nB = chains[1]["sequenceEnd"] - chains[1]["sequenceStart"] + 1
pae_matrix = np.array(pae_json[0]["predicted_aligned_error"], dtype=np.float32)
pae_AB = pae_matrix[:nA, nA:nA + nB]
pae_BA = pae_matrix[nA:nA + nB, :nA]
pae_AA = pae_matrix[:nA, :nA]
pae_BB = pae_matrix[nA:nA + nB, nA:nA + nB]
```

If the PAE `chains` field is absent, fall back to structure-derived chain lengths and record the assumption.

## Score Formulas

### ipTM

- Use all inter-chain PAE values.
- Use `d0_func(nA + nB)`.
- Compute both directions and take the maximum directional score.

### ipSAE

- Use PAE cutoff `< 10.0`.
- Preserve directional `pae_AB` and `pae_BA` scoring.
- Implement these variants:
  - `d0res`: `d0_func(n_valid)` per residue
  - `d0chn`: `d0_func(nA + nB)`
  - `d0dom`: `d0_func(n_dom)`, where `n_dom` is residues with any inter-chain PAE below cutoff
- Residues with no valid PAE pairs score `0.0`.

### pDockQ

- Interface contacts use CB-CB <= 8.0 Angstrom, CA fallback for glycine.
- If there are no contacts, return score `0.018`.
- Use:

```text
0.724 / (1 + exp(-0.052 * (x - 152.611))) + 0.018
```

where `x = mean_interface_plddt * log10(max(n_interface_residues, 1))`.

### pDockQ2

- Use contact-pair PAE values with fixed `d0 = 10.0`.
- If there are no contacts, return score `0.005`.
- Use the project reference constants from the relevant spec or reference implementation. If constants conflict across artifacts, stop and flag a human-review question.

### LIS

- Use PAE cutoff `< 12.0`.
- Score each direction as the mean of `(12.0 - PAE) / 12.0` over valid low-PAE inter-chain pairs.
- Return the maximum directional score.
- Explain LIS as density of low-PAE inter-chain interactions, not global docking confidence.

## Traffic-Light Thresholds

Use these defaults unless the spec overrides them:

```python
THRESHOLDS = {
    "ipsae_d0res": (0.6, 0.4),
    "iptm": (0.7, 0.5),
    "pdockq": (0.23, 0.09),
    "pdockq2": (0.5, 0.23),
    "lis": (0.15, 0.09),
}
```

## AFDB Data Contracts

Metadata endpoint:

```text
GET https://alphafold.ebi.ac.uk/api/prediction/{accession_id}
```

Consume the first JSON array item. Required or expected fields:

- `cifUrl`
- `bcifUrl`
- `paeDocUrl`
- `plddtDocUrl`
- `uniprotAccession`
- `sequence`
- `organismScientificName` or `organism`
- `proteinFullName` or `uniprotDescription`
- `geneNames`
- `modelVersion`

PAE JSON should contain `predicted_aligned_error`, `max_predicted_aligned_error`, and preferably `chains`.

pLDDT JSON should contain `residueNumber`, `confidenceScore`, `confidenceCategory`, and preferably `chains`.

## mmCIF Rules

Required atom-site fields:

- `group_PDB`
- `label_atom_id`
- `label_asym_id`
- `label_seq_id`
- `label_comp_id`
- `Cartn_x`, `Cartn_y`, `Cartn_z`
- `B_iso_or_equiv`
- `pdbx_PDB_model_num` when present

Filtering:

- keep only `group_PDB == "ATOM"`
- keep only `label_atom_id` in `("CA", "CB")`
- skip `label_seq_id` values `.` and `?`
- if `pdbx_PDB_model_num` exists, keep only model `1`

## Known Edge Cases

- zero inter-chain contacts: pDockQ `0.018`, pDockQ2 `0.005`
- residue with no valid PAE pairs for ipSAE: score `0.0`
- glycine has no CB atom: use CA
- missing PAE `chains`: use structure chain lengths and record limitation
- single-monomer accession: chain B is absent, scores should fail clearly or return documented non-informative values
- asymmetric chain labels such as A and C: sort and verify against PAE `chains`

## Quality Gate

Reject a scoring implementation if:

- PAE asymmetry is lost
- chain lengths and quadrant slices are undocumented
- contact definition differs from the spec without approval
- formula constants are invented from memory
- numeric validation is absent when a reference is available
- educational notebook prose fails to explain what each score measures and does not measure
