"""
CB/CA-distance interface detection for AlphaFold homodimer structures.

Contact definition: CB–CB Euclidean distance ≤ 8.0 Å.
GLY residues use CA because they have no CB atom.

NOTE — intentional divergence from AFDB production code:
  AFDB production interface.py (Majewski, Apache 2.0) uses CA–CA distances
  via PyTorch radius_graph on GPU. This module uses CB–CB to match the IPSAE
  scoring code's contact definition. Do not change the atom-type selection
  without updating dependent scoring functions (pDockQ, pDockQ2).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

@dataclass
class ChainCoords:
    """Coordinate arrays for one chain, indexed by residue position (1-based)."""
    chain_id: str
    coords: np.ndarray          # shape (N, 3), float32
    res_ids: np.ndarray         # shape (N,), int  — 1-based residue numbers
    res_names: np.ndarray       # shape (N,), str  — three-letter codes
    plddt: np.ndarray           # shape (N,), float32 — B-factor column


@dataclass
class InterfaceResult:
    """Output of detect_interface."""
    mask_A: np.ndarray          # shape (nA,) bool — True for interface residues in chain A
    mask_B: np.ndarray          # shape (nB,) bool — True for interface residues in chain B
    dist_matrix: np.ndarray     # shape (nA, nB) float32 — all pairwise CB/CA distances (Å)
    n_contacts: int             # total interface residue count (mask_A.sum() + mask_B.sum())
    contact_pairs: List[Tuple[int, int]] = field(default_factory=list)
    # ^ (i, j) index pairs where dist_matrix[i, j] <= cutoff


# ---------------------------------------------------------------------------
# mmCIF parsing
# ---------------------------------------------------------------------------

def parse_mmcif_atoms(cif_text: str) -> Tuple[List[List[str]], Dict[str, int]]:
    """
    Tokenise the _atom_site loop from an mmCIF text.

    Returns:
        records  — list of token lists, one per ATOM/HETATM line
        col_idx  — mapping from lowercase column name to position in each record

    The parser is column-order agnostic: field indices are built dynamically from
    the loop header, so it is robust to AFDB mmCIF version differences.
    """
    lines = cif_text.splitlines()
    loop_start = None

    for i, line in enumerate(lines):
        if line.strip().lower() == "loop_":
            j = i + 1
            headers: List[str] = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                headers.append(lines[j].strip())
                j += 1
            if any(h.lower().startswith("_atom_site.") for h in headers):
                loop_start = i
                break

    if loop_start is None:
        return [], {}

    j = loop_start + 1
    headers = []
    while j < len(lines) and lines[j].strip().startswith("_"):
        headers.append(lines[j].strip().lower())
        j += 1

    col_idx: Dict[str, int] = {h: idx for idx, h in enumerate(headers)}

    records: List[List[str]] = []
    while j < len(lines):
        stripped = lines[j].strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("_"):
            break
        if stripped.startswith("loop_"):
            break
        records.append(stripped.split())
        j += 1

    return records, col_idx


def extract_cb_coords(
    records: List[List[str]],
    col_idx: Dict[str, int],
) -> Dict[str, ChainCoords]:
    """
    Extract per-chain CB coordinates (CA fallback for GLY) from parsed atom records.

    Filtering applied:
        - group_PDB == 'ATOM' only
        - label_atom_id in ('CA', 'CB')
        - label_seq_id not in ('.', '?') and must be numeric
        - pdbx_PDB_model_num == '1' when that column is present

    Returns a dict keyed by chain ID (e.g. {'A': ChainCoords, 'B': ChainCoords}).
    Within each chain residues are sorted by residue number.
    """
    def _col(name: str) -> int | None:
        return col_idx.get(f"_atom_site.{name}")

    ix_group  = _col("group_pdb")
    ix_atom   = _col("label_atom_id")
    ix_chain  = _col("label_asym_id")
    ix_seq    = _col("label_seq_id")
    ix_resn   = _col("label_comp_id")
    ix_x      = _col("cartn_x")
    ix_y      = _col("cartn_y")
    ix_z      = _col("cartn_z")
    ix_bfact  = _col("b_iso_or_equiv")
    ix_model  = _col("pdbx_pdb_model_num")

    required = (ix_group, ix_atom, ix_chain, ix_seq, ix_resn, ix_x, ix_y, ix_z, ix_bfact)
    if any(c is None for c in required):
        raise ValueError("mmCIF is missing one or more required _atom_site columns.")

    # chain_id → {res_id → (atom_id, resname, x, y, z, plddt)}
    per_chain: Dict[str, Dict[int, Tuple]] = {}

    for toks in records:
        if len(toks) <= max(c for c in required if c is not None):
            continue
        if toks[ix_group] != "ATOM":
            continue
        if ix_model is not None and toks[ix_model] != "1":
            continue

        atom_id = toks[ix_atom]
        if atom_id not in ("CA", "CB"):
            continue

        res_str = toks[ix_seq]
        if res_str in (".", "?") or not res_str.lstrip("-").isdigit():
            continue

        chain_id = toks[ix_chain]
        res_id   = int(res_str)
        resname  = toks[ix_resn]

        try:
            x, y, z = float(toks[ix_x]), float(toks[ix_y]), float(toks[ix_z])
            plddt   = float(toks[ix_bfact])
        except ValueError:
            continue

        per_chain.setdefault(chain_id, {})
        existing = per_chain[chain_id].get(res_id)
        # Prefer CB over CA (CB overwrites an earlier CA record for the same residue)
        if existing is None or (atom_id == "CB" and existing[0] == "CA"):
            per_chain[chain_id][res_id] = (atom_id, resname, x, y, z, plddt)

    result: Dict[str, ChainCoords] = {}
    for chain_id, residues in per_chain.items():
        sorted_ids = sorted(residues)
        coords   = np.array([[residues[r][2], residues[r][3], residues[r][4]] for r in sorted_ids], dtype=np.float32)
        resnames = np.array([residues[r][1] for r in sorted_ids])
        plddt    = np.array([residues[r][5] for r in sorted_ids], dtype=np.float32)
        result[chain_id] = ChainCoords(
            chain_id=chain_id,
            coords=coords,
            res_ids=np.array(sorted_ids, dtype=np.int32),
            res_names=resnames,
            plddt=plddt,
        )
    return result


# ---------------------------------------------------------------------------
# Interface detection
# ---------------------------------------------------------------------------

def detect_interface(
    chain_a: ChainCoords,
    chain_b: ChainCoords,
    dist_cutoff: float = 8.0,
) -> InterfaceResult:
    """
    Detect interface residues between two chains using CB–CB Euclidean distances.

    GLY residues have no CB atom; extract_cb_coords already falls back to CA for those,
    so callers do not need to handle GLY specially.

    Args:
        chain_a:      Coordinates for chain A (from extract_cb_coords).
        chain_b:      Coordinates for chain B (from extract_cb_coords).
        dist_cutoff:  Contact distance threshold in Å (default 8.0).

    Returns:
        InterfaceResult with boolean interface masks, the full distance matrix,
        and a list of (i, j) index pairs for residues within the cutoff.

    Example
    -------
    >>> import numpy as np
    >>> from insightfold.interface import ChainCoords, detect_interface
    >>> a = ChainCoords('A',
    ...     coords=np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], dtype=np.float32),
    ...     res_ids=np.array([1, 2], dtype=np.int32),
    ...     res_names=np.array(['ALA', 'ALA']),
    ...     plddt=np.array([90.0, 85.0], dtype=np.float32))
    >>> b = ChainCoords('B',
    ...     coords=np.array([[3.0, 0.0, 0.0], [20.0, 0.0, 0.0]], dtype=np.float32),
    ...     res_ids=np.array([1, 2], dtype=np.int32),
    ...     res_names=np.array(['GLY', 'ALA']),
    ...     plddt=np.array([88.0, 80.0], dtype=np.float32))
    >>> result = detect_interface(a, b, dist_cutoff=8.0)
    >>> result.n_contacts      # A[0]+A[1] in interface (3 Å, 7 Å); only B[0] in interface (B[1] is 10 Å away)
    3
    >>> list(result.contact_pairs)
    [(0, 0), (1, 0)]
    """
    coords_a = chain_a.coords   # (nA, 3)
    coords_b = chain_b.coords   # (nB, 3)

    # Vectorised pairwise distance: broadcast over (nA, 1, 3) − (1, nB, 3)
    diff        = coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
    dist_matrix = np.sqrt((diff ** 2).sum(axis=-1)).astype(np.float32)   # (nA, nB)

    contact_mat   = dist_matrix <= dist_cutoff
    mask_a        = contact_mat.any(axis=1)   # residues in A with any B contact
    mask_b        = contact_mat.any(axis=0)   # residues in B with any A contact
    contact_pairs = list(zip(*np.where(contact_mat)))

    return InterfaceResult(
        mask_A=mask_a,
        mask_B=mask_b,
        dist_matrix=dist_matrix,
        n_contacts=int(mask_a.sum()) + int(mask_b.sum()),
        contact_pairs=[(int(i), int(j)) for i, j in contact_pairs],
    )
