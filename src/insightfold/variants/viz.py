from __future__ import annotations

import pandas as pd


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
    Legacy structure-pin viewer entrypoint.

    InsightFold now standardizes 3D protein structure visualization on
    MolViewSpec/Mol*. Reimplement this view with the molviewspec-rendering
    skill before using it in notebooks.
    """
    raise NotImplementedError(
        "Structure pin rendering must be implemented with MolViewSpec/Mol*."
    )


def centre_on_position(view, pos: int, chain: str = "A"):
    """Recentre the camera on a given residue position."""
    raise NotImplementedError(
        "Structure camera controls must be implemented with MolViewSpec/Mol*."
    )
