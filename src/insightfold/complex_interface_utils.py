"""
Shared utilities for AlphaFold DB protein-complex interface analysis.

Scope
-----
This module backs `notebooks/homodimer_diagnostic.ipynb`. It covers the whole
pipeline for a **predicted protein complex** served by the AlphaFold Database:
fetch metadata and the mmCIF / PAE / pLDDT documents, parse them, detect the
interface, compute the seven confidence values, and render the plots and 3D views.

Two architectural decisions shape every signature here (see
`specs/homodimer_diagnostic/rework-plan.md`, "Locked decisions"):

- **D3 — homodimers and heterodimers are both supported automatically.** Nothing
  assumes the two chains are identical, or even the same length.
- **D4 — the architecture is dimer-scoped, but the API is chain-pair-generic.**
  Functions take an explicit ordered chain pair (a directional PAE block plus the
  two chain lengths, or a chain-pair object) rather than assuming "A and B". An
  N-chain generalisation is then a later extension over all ordered pairs, not a
  rewrite.

Because PAE is asymmetric, "the chain pair" always means an *ordered* pair.
Where `ipsae.py` reports a combination over both directions (`max`, or the mean
for LIS), that combination is performed explicitly and is documented per score.

Numerical ground truth
----------------------
`DunbrackLab/IPSAE` `ipsae.py` v4 (Jan 2026), as transcribed and verified in
`specs/homodimer_diagnostic/formula-reference.md` (D2). Tolerance for agreement
with the reference is +/-0.001. Publications are used for explanatory prose and
threshold rationale only; the thresholds themselves come from
`specs/homodimer_diagnostic/threshold-reference.md`.

Dependency policy
-----------------
Allowed at module top level: `numpy`, `requests`, `matplotlib`, `seaborn`
(imports are added by the task that first needs them, so this file stays cheap
to import). `molviewspec` is imported **lazily**, inside the function that needs
it, so the notebook degrades gracefully when it is not installed and so a plain
`import complex_interface_utils` never pays for it.

Prohibited anywhere in this module: `biopython` / `Bio`, `torch`,
`torch-geometric`, `gemmi`, `pandas`, `scipy`, `plotly`. The notebook must run on
Colab free tier within a 60 s install budget, which is also why the Colab
bootstrap clones the repo and extends `sys.path` instead of `pip install`-ing the
package (D9).

Status
------
This file is being filled in over several tasks. Sections below carry a note
naming the task that populates them. Only the genuinely shared primitives
(the two `d0` helpers, `ptm_func`, the threshold table and the traffic light)
are implemented here, because everything else depends on them.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Literal, Tuple

import numpy as np

__all__ = [
    # constants
    "PAE_CUTOFF",
    "DIST_CUTOFF",
    "LIS_CUTOFF",
    # scoring primitives
    "d0_scalar",
    "d0_array",
    "ptm_func",
    # thresholds
    "Provenance",
    "Band",
    "Threshold",
    "THRESHOLDS",
    "traffic_light",
]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
# Cutoffs are module-level so that plot subtitles and prose can name the value
# actually in force rather than hard-coding a duplicate copy of it.

PAE_CUTOFF: float = 10.0
"""Inter-chain PAE cutoff for the three ipSAE variants, in Angstrom.

`ipsae.py` exposes this as a CLI argument defaulting to 10.0 (`ipsae_v4.py:58`)
and tests it strictly: a pair contributes when `pae < PAE_CUTOFF`
(`ipsae_v4.py:737`). Does not apply to ipTM_d0chn, which uses every partner
residue with no cutoff at all.
"""

DIST_CUTOFF: float = 8.0
"""CB-CB contact distance cutoff, in Angstrom, for pDockQ and pDockQ2.

`ipsae.py` `pDockQ_cutoff` (`ipsae_v4.py:644`), tested inclusively:
a pair is a contact when `dist <= DIST_CUTOFF` (`ipsae_v4.py:652`).
GLY has no CB atom and uses CA instead.
"""

LIS_CUTOFF: float = 12.0
"""PAE cutoff for LIS, in Angstrom.

Hardcoded to 12 in `ipsae.py` (`ipsae_v4.py:712`) and *independent* of
`PAE_CUTOFF`: Kim 2024 chose 12 as the value maximising AUC in their ROC
analysis, and the published LIS threshold is only valid at that cutoff.
Tested strictly (`pae < LIS_CUTOFF`).
"""


# ---------------------------------------------------------------------------
# AFDB access  --  filled by R011
# ---------------------------------------------------------------------------
# `fetch_metadata`, `download_structure`, `download_pae`, `download_plddt`.
# Uses `requests` (not `urllib`): the PAE and pLDDT documents are served gzipped
# and `requests` decompresses them transparently.


# ---------------------------------------------------------------------------
# Structure parsing  --  filled by R011
# ---------------------------------------------------------------------------
# `parse_mmcif_atoms`, `extract_chain_coords`. One column-naming convention for
# the `_atom_site.` prefix, chosen here and used everywhere (see R015, which
# absorbs the three divergent copies of `extract_cb_coords`).


# ---------------------------------------------------------------------------
# PAE / pLDDT parsing  --  filled by R011
# ---------------------------------------------------------------------------
# `parse_pae`, `parse_plddt`, and the ordered-chain-pair quadrant extraction that
# D4 makes explicit. Chain lengths derived here must be asserted equal to the
# structure-derived lengths rather than silently mis-slicing (R021).


# ---------------------------------------------------------------------------
# Interface detection  --  filled by R011
# ---------------------------------------------------------------------------
# CB-CB (CA for GLY) contact detection at `DIST_CUTOFF`, plus the `ChainCoords` /
# interface-result carriers that D4's explicit chain pair rides on. R015 absorbs
# the useful parts of `src/insightfold/interface.py` here and resolves its
# `n_contacts` name collision into `n_interface_residues` (residues) versus
# `n_contact_pairs` (pairs; this is the one pDockQ needs).


# ---------------------------------------------------------------------------
# Scoring primitives
# ---------------------------------------------------------------------------
# Implemented now: every score below depends on these, so they cannot wait for
# R012. Ground truth `ipsae_v4.py:111-136`, see formula-reference.md section 0.

def d0_scalar(L: float | int) -> float:
    """
    Scalar TM-score length normalisation, matching `ipsae.py` `calc_d0`.

    Used for `d0chn` (`ipsae_v4.py:732`) and `d0dom` (`ipsae_v4.py:778`), i.e.
    wherever `L` is a single whole-chain or whole-domain residue count.

    `1.0` for `L <= 27`; otherwise `max(1.0, 1.24 * (L - 15) ** (1/3) - 1.8)`.

    Args:
        L: Number of residues in the scoring domain.

    Returns:
        The `d0` normalisation constant.

    Note:
        **This is not interchangeable with `d0_array`.** They differ at exactly
        one point, `L == 27`, where this returns `1.0` and `d0_array` returns
        `1.038891`. The difference comes from the `else: d0 = 1.0` branch at
        `ipsae_v4.py:122-123`, which `calc_d0_array` has no counterpart for.
        `ipsae.py` genuinely uses both, so both are kept (R005).

    Example
    -------
    >>> d0_scalar(10)
    1.0
    >>> d0_scalar(27)
    1.0
    >>> round(d0_scalar(28), 6)
    1.115655
    >>> round(d0_scalar(344), 6)
    6.760261
    """
    L = float(L)
    if L > 27:
        return max(1.0, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)
    return 1.0


def d0_array(L: np.ndarray | float | int) -> np.ndarray:
    """
    Vectorised TM-score length normalisation, matching `ipsae.py` `calc_d0_array`.

    Used for `d0res` only (`ipsae_v4.py:787`), where `L` is a per-residue count
    of valid inter-chain pairs and therefore an array.

    Clamps `L` to at least 26, then applies the cubic form unconditionally with a
    `max(1.0, ...)` floor. The floor, not the clamp, is what makes `L < 26`
    return `1.0` (the clamped cubic at `L = 26` is 0.9578).

    Args:
        L: Residue count(s). Scalar or array; always returns an ndarray
           (0-d for a scalar input) so downstream broadcasting is uniform.

    Returns:
        `d0` normalisation constant(s), same shape as `L`.

    Note:
        Differs from `d0_scalar` at `L == 27` only, and deliberately so; see
        `d0_scalar`. Do not "simplify" the two into one function.

    Example
    -------
    >>> float(d0_array(10))
    1.0
    >>> round(float(d0_array(27)), 6)
    1.038891
    >>> round(float(d0_array(28)), 6)
    1.115655
    >>> d0_array(np.array([0, 27, 28])).shape
    (3,)
    """
    L_arr = np.maximum(26.0, np.asarray(L, dtype=float))
    return np.maximum(1.0, 1.24 * (L_arr - 15.0) ** (1.0 / 3.0) - 1.8)


def ptm_func(x: np.ndarray | float, d0: np.ndarray | float) -> np.ndarray | float:
    """
    TM-score transform of a PAE value, `1 / (1 + (x / d0) ** 2)`.

    `ipsae_v4.py:111-112`. Broadcast-safe: `x` and `d0` may each be a scalar or
    an ndarray, in any combination, so a whole PAE block can be transformed with
    one shared `d0` or with a per-row `d0` column vector.

    Args:
        x:  PAE value(s), in Angstrom.
        d0: Normalisation constant(s) from `d0_scalar` or `d0_array`.

    Returns:
        Value(s) in (0, 1]; 1.0 at `x == 0`, decreasing monotonically in `x`.

    Example
    -------
    >>> ptm_func(0.0, 10.0)
    1.0
    >>> ptm_func(10.0, 10.0)
    0.5
    >>> ptm_func(np.array([0.0, 10.0]), 10.0).tolist()
    [1.0, 0.5]
    """
    return 1.0 / (1.0 + (x / d0) ** 2)


# ---------------------------------------------------------------------------
# Score functions  --  filled by R012
# ---------------------------------------------------------------------------
# `compute_iptm_d0chn`, `compute_ipsae` (d0res / d0chn / d0dom), `compute_pdockq`,
# `compute_pdockq2`, `compute_lis`. Each takes an ordered chain pair per D4 and
# returns the per-residue intermediates alongside the headline value, so the
# plots and 3D views never recompute. R003 (per-direction d0dom), R005 (d0_array
# for d0res), R006 (zero-contact returns) and R007 (max over both pDockQ2
# directions) land with this section.


# ---------------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------------
# Implemented now: it is the single source of truth the summary table, the prose
# and the 3D colour maps all read from, so nothing downstream can invent its own
# copy. Transcribed verbatim from the Decision table in
# `specs/homodimer_diagnostic/threshold-reference.md` (R002). Do not edit the
# numbers here; edit that document and re-transcribe.

Provenance = Literal["PUBLISHED", "DERIVED", "AFDB-EMPIRICAL", "HEURISTIC"]
"""How much a threshold is actually worth.

- ``PUBLISHED``       the paper explicitly states or applies this cutoff as an
                      operating criterion.
- ``DERIVED``         obtained from the paper's own quantitative statements, but
                      not stated as a recommendation.
- ``AFDB-EMPIRICAL``  not in any publication, but established by AlphaFold DB
                      from a large-scale internal analysis (~31M predicted
                      dimers). Empirically grounded; simply unpublished.
- ``HEURISTIC``       no literature basis and no large-scale calibration. A
                      judgement call, and must be presented as one.
"""

Band = Literal["green", "amber", "red"]
"""Traffic-light band. Labels are HIGH / MODERATE / LOW respectively."""

_BAND_LABELS: Dict[str, str] = {
    "green": "HIGH",
    "amber": "MODERATE",
    "red": "LOW",
}


@dataclass(frozen=True)
class Threshold:
    """
    One score's green/amber cutoffs together with the provenance of each.

    Provenance is carried alongside the numbers deliberately: six of the fourteen
    values in the canonical table are judgement calls, and a user reading a
    coloured badge otherwise has no way to tell a published ROC-optimal cutoff
    from a round number someone picked (R080).

    Attributes:
        green:            `value >= green` is green/HIGH.
        amber:            `green > value >= amber` is amber/MODERATE;
                          `value < amber` is red/LOW.
        green_provenance: Provenance label for `green`.
        amber_provenance: Provenance label for `amber`.
        source:           Short citation for the green threshold.
        note:             One-line caveat to surface next to the badge.
    """

    green: float
    amber: float
    green_provenance: Provenance
    amber_provenance: Provenance
    source: str = ""
    note: str = ""

    def band(self, value: float) -> Band:
        """Return the band `value` falls in. Both edges compare with `>=`."""
        if value >= self.green:
            return "green"
        if value >= self.amber:
            return "amber"
        return "red"

    def provenance(self, band: Band) -> Provenance | None:
        """Provenance of the cutoff that defines `band`; `None` for red."""
        if band == "green":
            return self.green_provenance
        if band == "amber":
            return self.amber_provenance
        return None


THRESHOLDS: Dict[str, Threshold] = {
    "ipsae_d0res": Threshold(
        green=0.60,
        amber=0.30,
        green_provenance="AFDB-EMPIRICAL",
        amber_provenance="HEURISTIC",
        source="AlphaFold DB release filter, derived data-drivenly from ~31M "
               "predicted dimers. Not Dunbrack's number: the ipSAE paper states "
               "no cutoff of any kind.",
        note="For an AFDB accession the model has already passed this filter, so "
             "green is largely re-measuring the criterion that caused the entry "
             "to exist and red will essentially never appear. For a locally "
             "uploaded model the full range is reachable. Quote this variant, "
             "not d0chn or d0dom.",
    ),
    "ipsae_d0chn": Threshold(
        green=0.60,
        amber=0.30,
        green_provenance="AFDB-EMPIRICAL",
        amber_provenance="HEURISTIC",
        source="Inherited from ipsae_d0res; no publication states a cutoff for "
               "this variant.",
        note="ipSAE_d0chn >= ipSAE_d0dom >= ipSAE_d0res is a theorem, not an "
             "observation, so a shared threshold is strictly most permissive "
             "here. This variant turns green first and red last.",
    ),
    "ipsae_d0dom": Threshold(
        green=0.60,
        amber=0.30,
        green_provenance="AFDB-EMPIRICAL",
        amber_provenance="HEURISTIC",
        source="Inherited from ipsae_d0res; no publication states a cutoff for "
               "this variant.",
        note="Lies between d0res and d0chn by the same theorem. Three green "
             "lights are one measurement seen through three progressively more "
             "forgiving normalisations, not three confirmations.",
    ),
    "iptm_d0chn": Threshold(
        green=0.70,
        amber=0.30,
        green_provenance="DERIVED",
        amber_provenance="DERIVED",
        source="Dunbrack 2025, p. 14: true and false dimer densities overlap "
               "from 0.3 to 0.7 for ipTM computed from the PAE matrix.",
        note="This is ipTM_d0chn, a PAE-derived reimplementation, NOT AlphaFold's "
             "own ipTM, which AFDB does not expose on these endpoints. Amber "
             "means the score genuinely cannot separate a true dimer from a "
             "decoy in that range.",
    ),
    "pdockq": Threshold(
        green=0.23,
        amber=0.12,
        green_provenance="DERIVED",
        amber_provenance="DERIVED",
        source="Bryant 2022, pp. 1, 3, 10: the CAPRI/DockQ acceptability "
               "boundary, transferred because pDockQ is a sigmoid fitted "
               "directly to DockQ. Amber from the stated average error of 0.11 "
               "(pp. 4, 7).",
        note="Symmetric under swapping the chains, so it has no directional "
             "split. Zhu 2023 p. 5 found it over-optimistic on large, confident, "
             "incorrect interfaces, so a green pDockQ is never confirmation on "
             "its own.",
    ),
    "pdockq2": Threshold(
        green=0.23,
        amber=0.10,
        green_provenance="PUBLISHED",
        amber_provenance="HEURISTIC",
        source="Zhu 2023, p. 6: 'min pDockQ2 > 0.23' applied as an operating "
               "cutoff to select real biological conclusions.",
        note="The sigmoid is bounded on [0.005, 1.315] and can exceed 1.0. Never "
             "present pDockQ2 as a probability or a percentage.",
    ),
    "lis": Threshold(
        green=0.21,
        amber=0.10,
        green_provenance="PUBLISHED",
        amber_provenance="HEURISTIC",
        source="Kim 2024, p. 30: 'best LIS values exceeding the optimal "
               "threshold (0.21)', ROC/Youden-optimal.",
        note="Valid only at a PAE cutoff of 12 and for the mean of the two "
             "directions, which is what LIS_CUTOFF and compute_lis use. Answers "
             "whether the pair interacts at all, not whether the pose is right, "
             "so green LIS with red pDockQ is a legitimate combination.",
    ),
}
"""Canonical green/amber thresholds for the seven confidence values.

Single source of truth. Transcribed from the Decision table in
`specs/homodimer_diagnostic/threshold-reference.md`; that document carries the
full derivation for every number and supersedes the four conflicting copies that
previously lived in `CLAUDE.md` and in notebook cells 25, 26 and 28.

**Amber is a stored number and must never be recomputed as `green / 2`** — that
recomputation is exactly why the summary table and the prose could disagree on a
single run.

All five remaining judgement calls are ambers. Every green now has either a
publication or the AFDB calibration behind it; none of the four reference papers
defines a middle band at all, so the amber tier is this notebook's own construct.
"""


def traffic_light(value: float, score_name: str) -> Tuple[Band, str]:
    """
    Classify a score against its canonical threshold.

    Args:
        value:      The computed score.
        score_name: Key into `THRESHOLDS`, e.g. `'ipsae_d0res'`, `'iptm_d0chn'`,
                    `'pdockq'`, `'pdockq2'`, `'lis'`.

    Returns:
        `(band, label)` where band is `'green'` / `'amber'` / `'red'` and label
        is `'HIGH'` / `'MODERATE'` / `'LOW'`.

    Raises:
        KeyError: If `score_name` is not one of the seven canonical values. This
            is deliberate: a typo must fail loudly rather than silently defaulting
            to a band.

    Note:
        The threshold's provenance and caveat are not returned here; read them
        from `THRESHOLDS[score_name]`, which any consumer showing a badge should
        do (R080).

    Example
    -------
    >>> traffic_light(0.75, 'ipsae_d0res')
    ('green', 'HIGH')
    >>> traffic_light(0.45, 'ipsae_d0res')
    ('amber', 'MODERATE')
    >>> traffic_light(0.05, 'lis')
    ('red', 'LOW')
    """
    if score_name not in THRESHOLDS:
        raise KeyError(
            f"Unknown score name {score_name!r}. "
            f"Expected one of: {', '.join(sorted(THRESHOLDS))}."
        )
    band = THRESHOLDS[score_name].band(value)
    return band, _BAND_LABELS[band]


# ---------------------------------------------------------------------------
# Plotting  --  filled by R013
# ---------------------------------------------------------------------------
# One `plot_*` function per figure, each returning a `fig`, so a notebook cell is
# one call plus a title. Adds the `matplotlib` / `seaborn` top-level imports and
# the `PAE_CMAP` palette setting (R051): default is a sequential green with dark
# = low PAE = confident, plus a colourblind-safe alternative and `RdBu_r` for
# continuity.


# ---------------------------------------------------------------------------
# MolViewSpec views  --  filled by R014
# ---------------------------------------------------------------------------
# One builder function per 3D view, plus the shared `show_mol_view` helper.
#
# `molviewspec` MUST be imported lazily inside each builder, never at module top
# level, so that importing this module stays free and so the notebook degrades
# gracefully with a clear message when the package is absent.
