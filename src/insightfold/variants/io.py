from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import requests

# Public endpoints (kept here so they are easy to update in one place)
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
UNIPROT_FIELDS_URL = (
    "https://rest.uniprot.org/uniprotkb/{acc}.json?fields=ft_variant%2Cft_mutagen%2Csequence"
)
AFDB_PREDICTION_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc_or_iso}"


# ---------- minimal mmCIF pLDDT reader (no external deps) ----------

def _plddt_from_mmcif_text(
    mmcif_text: str,
    chain_preference: Tuple[str, ...] = ("A",),
) -> Tuple[Dict[int, float], int]:
    """
    Minimal AlphaFold mmCIF parser:
      - locate the atom_site loop
      - select CA atoms
      - read _atom_site.B_iso_or_equiv as pLDDT
      - use auth_seq_id (or label_seq_id) as 1-based residue index
    Returns (plddt_by_pos, max_position).
    """
    lines = mmcif_text.splitlines()

    # Find an 'atom_site' loop
    loop_start = None
    for i, line in enumerate(lines):
        if line.strip().lower() == "loop_":
            j = i + 1
            headers: list[str] = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                headers.append(lines[j].strip())
                j += 1
            if any(h.lower().startswith("_atom_site.") for h in headers):
                loop_start = i
                break
    if loop_start is None:
        return {}, 0

    # Build header → index map
    j = loop_start + 1
    headers: list[str] = []
    while j < len(lines) and lines[j].strip().startswith("_"):
        headers.append(lines[j].strip())
        j += 1
    name_to_ix = {h: ix for ix, h in enumerate(headers)}

    def col(exact_lower_name: str) -> Optional[int]:
        for k, ix in name_to_ix.items():
            if k.lower() == exact_lower_name.lower():
                return ix
        return None

    ix_atom  = col("_atom_site.label_atom_id")
    ix_chain = col("_atom_site.label_asym_id")
    ix_seq   = col("_atom_site.auth_seq_id") or col("_atom_site.label_seq_id")
    ix_bfact = col("_atom_site.B_iso_or_equiv")
    if None in (ix_atom, ix_chain, ix_seq, ix_bfact):
        return {}, 0

    # Data lines follow until the loop ends (no longer ATOM/HETATM)
    data_end = j
    while data_end < len(lines) and (lines[data_end].startswith("ATOM") or lines[data_end].startswith("HETATM")):
        data_end += 1

    def collect(prefer_chains: Optional[Tuple[str, ...]]) -> Tuple[Dict[int, float], int]:
        out: Dict[int, float] = {}
        for r in range(j, data_end):
            toks = lines[r].split()
            # Be defensive about malformed rows
            try:
                if toks[ix_atom] != "CA":
                    continue
                if prefer_chains and toks[ix_chain] not in prefer_chains:
                    continue
                pos = int(float(toks[ix_seq]))
                b = float(toks[ix_bfact])
                out[pos] = b
            except Exception:
                continue
        return out, (max(out) if out else 0)

    plddt, maxpos = collect(chain_preference)
    if not plddt:
        plddt, maxpos = collect(None)  # fallback to any chain
    return plddt, maxpos


# ---------- HTTP helpers ----------

def _http_json(url: str, *, timeout: int = 30) -> dict:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.json()


def _http_text(url: str, *, timeout: int = 60) -> str:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.text


# ---------- AFDB ID building ----------

def build_afdb_id(uniprot_acc: str, isoform: Optional[str]) -> str:
    """Return AFDB accession for canonical or a specific isoform (F1 by design)."""
    iso = (isoform or "").strip()
    return f"AF-{uniprot_acc}-{iso}-F1" if iso else f"AF-{uniprot_acc}-F1"


# ---------- AFDB model metadata & mmCIF ----------

@dataclass
class ModelInfo:
    afdb_id: str
    mmcif_url: str
    model_version: Optional[str]
    af_crc64: Optional[str]   # from sequenceChecksum only
    mmcif_text: str
    seq_from_cif: str         # length synthesised from max residue index
    plddt_by_pos: Dict[int, float]  # 1-based position → pLDDT


def afdb_prediction_meta(uniprot_acc: str, isoform: Optional[str]) -> dict:
    """
    Fetch AFDB prediction metadata. The API sometimes returns a list; we normalise to a dict.
    """
    key = f"{uniprot_acc}-{isoform}" if (isoform or "").strip() else uniprot_acc
    js = _http_json(AFDB_PREDICTION_URL.format(acc_or_iso=key))
    if isinstance(js, list):
        return js[0] if js else {}
    return js


# ---------- CRC64 (ISO) for checksum matching ----------

def _crc64_table_iso():
    POLY = 0xD800000000000000
    tbl = []
    for i in range(256):
        crc = i << 56
        for _ in range(8):
            crc = (crc << 1) ^ POLY if (crc & (1 << 63)) else (crc << 1)
            crc &= 0xFFFFFFFFFFFFFFFF
        tbl.append(crc)
    return tbl


_CRC64_TBL = _crc64_table_iso()


def crc64_iso(s: str) -> str:
    """Return uppercase 16-hex CRC64 (ISO) string compatible with UniProt checksums."""
    crc = 0
    for ch in s.encode("utf-8"):
        idx = ((crc >> 56) ^ ch) & 0xFF
        crc = _CRC64_TBL[idx] ^ ((crc << 8) & 0xFFFFFFFFFFFFFFFF)
    return f"{crc:016X}"


# ---------- High level fetch ----------

def fetch_uniprot_entry(uniprot_acc: str) -> dict:
    """
    Fetch a trimmed UniProtKB JSON that includes features and sequence.
    Falls back to the full entry if fields filtering is unavailable.
    """
    try:
        return _http_json(UNIPROT_FIELDS_URL.format(acc=uniprot_acc))
    except Exception:
        return _http_json(UNIPROT_ENTRY_URL.format(acc=uniprot_acc))


def fetch_model_info(uniprot_acc: str, isoform: Optional[str]) -> ModelInfo:
    """
    Fetch AFDB metadata, enforce presence of 'sequenceChecksum', download mmCIF,
    and extract per-residue pLDDT (from CA B-factors) with a tiny built-in parser.
    """
    meta = afdb_prediction_meta(uniprot_acc, isoform)
    lower = {k.lower(): v for k, v in meta.items()}

    mmcif_url = (
        meta.get("cifUrl") or lower.get("cifurl") or
        meta.get("mmCifUrl") or lower.get("mmcifurl") or
        meta.get("modelUrl") or lower.get("modelurl")
    )
    if not (mmcif_url and str(mmcif_url).endswith(".cif")):
        raise RuntimeError("Could not locate a valid mmCIF URL in AFDB metadata.")

    version = meta.get("modelVersion") or lower.get("modelversion")

    # STRICT: use sequenceChecksum only
    af_crc64 = meta.get("sequenceChecksum") or lower.get("sequencechecksum")
    if not af_crc64:
        raise RuntimeError(
            "AlphaFold DB metadata does not contain 'sequenceChecksum'. "
            "Cannot perform the required strict CRC64 comparison."
        )
    af_crc64 = str(af_crc64).upper()

    afdb_id = build_afdb_id(uniprot_acc, isoform)

    # Fetch mmCIF and read pLDDT
    mmcif_text = _http_text(str(mmcif_url))
    plddt_by_pos, maxpos = _plddt_from_mmcif_text(mmcif_text, chain_preference=("A",))
    # Synthetise sequence length from max position observed (for mapping guards)
    seq_from_cif = "X" * max(0, maxpos)

    return ModelInfo(
        afdb_id=afdb_id,
        mmcif_url=str(mmcif_url),
        model_version=str(version) if version is not None else None,
        af_crc64=af_crc64,
        mmcif_text=mmcif_text,
        seq_from_cif=seq_from_cif,
        plddt_by_pos=plddt_by_pos,
    )
