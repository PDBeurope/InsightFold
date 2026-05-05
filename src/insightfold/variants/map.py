from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .normalise import UniProtFeature


def tidy_from_features(
    uniprot_acc: str,
    afdb_id: str,
    feats: List[UniProtFeature],
    plddt_by_pos: Dict[int, float],
    model_len: int,
    feature_selection: str,  # "Variants" | "Mutagenesis" | "Both"
    isoform: Optional[str] = None,
) -> pd.DataFrame:
    """
    Construct the tidy table expected by the notebook, with out-of-range rows flagged.
    """
    allow_variant = feature_selection in ("Variants", "Both")
    allow_mut = feature_selection in ("Mutagenesis", "Both")

    rows = []
    for f in feats:
        if f.feature_type == "VARIANT" and not allow_variant:
            continue
        if f.feature_type == "MUTAGEN" and not allow_mut:
            continue
        pos = f.pos_uniprot
        if pos is None:
            continue
        in_model = 1 <= pos <= model_len
        rows.append(
            {
                "accession": uniprot_acc,
                "afdb_id": afdb_id,
                "isoform": isoform or None,
                "feature_type": f.feature_type,
                "pos_uniprot": int(pos),
                "pos_model": int(pos) if in_model else None,
                "in_model": bool(in_model),
                "ref": f.ref,
                "alt": f.alt,
                "label": f.label,
                "source_xrefs": ";".join(f.xrefs) if f.xrefs else "",
                "frequency": float(f.frequency) if f.frequency is not None else None,
                "sift": f.sift,
                "polyphen": f.polyphen,
                "plddt": float(plddt_by_pos.get(int(pos), np.nan)) if in_model else np.nan,
                "notes": None if in_model else "out_of_range",
            }
        )

    cols = [
        "accession",
        "afdb_id",
        "isoform",
        "feature_type",
        "pos_uniprot",
        "pos_model",
        "in_model",
        "ref",
        "alt",
        "label",
        "source_xrefs",
        "frequency",
        "sift",
        "polyphen",
        "plddt",
        "notes",
    ]
    df = pd.DataFrame(rows, columns=cols)
    return df
