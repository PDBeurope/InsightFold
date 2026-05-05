from __future__ import annotations

from typing import Optional

import pandas as pd

try:
    import py3Dmol
except Exception as _e:  # pragma: no cover
    py3Dmol = None  # allows import of this module without the viewer installed


def show_structure_with_pins(
    mmcif_text: str,
    df: pd.DataFrame,
    *,
    colour_by_plddt: bool = True,
    label_pins: bool = True,
    width: int = 700,
    height: int = 500,
):
    """
    Render an mmCIF model with optional pins at the selected positions.
    Returns the py3Dmol view object. Caller should call .show() in notebooks.
    """
    if py3Dmol is None:
        raise RuntimeError("py3Dmol is not installed. Please `pip install py3Dmol`.")

    view = py3Dmol.view(width=width, height=height)
    view.addModel(mmcif_text, "mmcif")

    # Cartoon with a standard scheme. pLDDT colouring can be added later if desired.
    view.setStyle({"cartoon": {"colorscheme": "bPyMOL"}})

    if df is not None and not df.empty:
        for _, r in df[df["in_model"] & df["pos_model"].notna()].iterrows():
            pos = int(r["pos_model"])
            selection = {"chain": "A", "resi": pos}
            view.addStyle(selection, {"stick": {}})
            if label_pins:
                label = f"{r['feature_type']} {r['ref'] or ''}{pos}{r['alt'] or ''}"
                view.addLabel(label, {"fontSize": 10, "backgroundOpacity": 0.5}, selection)

    view.zoomTo()
    return view


def centre_on_position(view, pos: int, chain: str = "A"):
    """Recentre the camera on a given residue position."""
    if pos is None:
        return
    view.zoomTo({"chain": chain, "resi": int(pos)})
    view.update()
