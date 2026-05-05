from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple, Any

@dataclass
class UniProtFeature:
    feature_type: str          # "VARIANT" (Proteins API does not expose MUTAGEN)
    pos_uniprot: Optional[int]
    ref: Optional[str]
    alt: Optional[str]
    label: Optional[str]
    xrefs: list[str]
    frequency: Optional[float]
    sift: Optional[str]
    polyphen: Optional[str]

def _pos_from_feature(f: dict) -> Optional[int]:
    # Proteins API: begin/end are strings; single-site is begin==end
    b = f.get("begin"); e = f.get("end")
    try:
        if b is None or e is None:
            return None
        if str(b) == str(e):
            return int(str(b))
    except Exception:
        return None
    return None  # ranges ignored in v1

def _collect_xrefs(f: dict) -> list[str]:
    xs: list[str] = []
    for xr in f.get("xrefs", []) or []:
        name = xr.get("name")
        xid  = xr.get("id")
        if name:
            xs.append(f"{name}:{xid}" if xid else name)
    # also capture evidence sources (often pubmed / ClinVar tags)
    for ev in f.get("evidences", []) or []:
        src = ev.get("source") or {}
        name = src.get("name")
        sid  = src.get("id")
        if name:
            xs.append(f"{name}:{sid}" if sid else name)
    # dedupe preserve order
    seen, out = set(), []
    for t in xs:
        if t not in seen:
            seen.add(t); out.append(t)
    return out

def _maybe_frequency(f: dict) -> Optional[float]:
    # Proteins API sometimes has population frequency blocks; keep it simple
    # by scanning for numeric keys containing 'freq' in known places.
    cand_keys = ("frequency", "alleleFrequency", "af")
    for k in cand_keys:
        v = f.get(k)
        if isinstance(v, (int, float)):
            return float(v)
    # nested lists like populationFrequencies: [{"frequency": 0.001, ...}]
    for k in ("populationFrequencies", "frequencies"):
        arr = f.get(k)
        if isinstance(arr, list):
            for item in arr:
                if isinstance(item, dict):
                    v = item.get("frequency")
                    if isinstance(v, (int, float)):
                        return float(v)
    return None

def _maybe_predictions(f: dict) -> tuple[Optional[str], Optional[str]]:
    # Some payloads include predictions/scores; pass through if present.
    sift = f.get("siftPrediction") or f.get("siftScore")
    poly = f.get("polyphenPrediction") or f.get("polyphenScore")
    return (str(sift) if sift is not None else None,
            str(poly) if poly is not None else None)

def _label_from_feature(f: dict) -> Optional[str]:
    # Compose a short, human-friendly label:
    # take first description value if present, else consequence type,
    # else join clinicalSignificances types.
    if f.get("descriptions"):
        val = f["descriptions"][0].get("value")
        if val:
            return val
    if f.get("consequenceType"):
        return str(f["consequenceType"])
    if f.get("clinicalSignificances"):
        try:
            return ", ".join(d.get("type") for d in f["clinicalSignificances"] if d.get("type"))
        except Exception:
            pass
    return None

def _aa1(x: Optional[str]) -> Optional[str]:
    if not x:
        return x
    s = str(x).strip().upper()
    return s[:1] if s else s

def parse_uniprot_features(entry: dict, allow_mutagen: bool) -> tuple[str, str, List[UniProtFeature]]:
    """
    Parse the Proteins Variation API payload into a normalised feature list.
    Returns (sequence, crc64_or_empty, features).
    - We return an empty CRC64 so the caller can compute CRC64( sequence ) with crc64_iso.
    - The Proteins API does not expose MUTAGEN; 'allow_mutagen' is ignored for now.
    """
    seq = str(entry.get("sequence") or "")
    crc64 = ""  # compute later with crc64_iso(seq) in the notebook

    feats: List[UniProtFeature] = []
    for f in entry.get("features", []) or []:
        ftype = f.get("type") or "VARIANT"
        if ftype != "VARIANT":
            # Proteins Variation API is variant-focussed; ignore anything else for v1.
            continue

        pos = _pos_from_feature(f)
        ref = _aa1(f.get("wildType"))
        alt = _aa1(f.get("mutatedType") or f.get("alternativeSequence"))
        xrefs = _collect_xrefs(f)
        freq = _maybe_frequency(f)
        sift, poly = _maybe_predictions(f)
        label = _label_from_feature(f)

        feats.append(UniProtFeature(
            feature_type="VARIANT",
            pos_uniprot=pos,
            ref=ref, alt=alt,
            label=label,
            xrefs=xrefs,
            frequency=freq,
            sift=sift,
            polyphen=poly,
        ))

    return seq, crc64, feats
