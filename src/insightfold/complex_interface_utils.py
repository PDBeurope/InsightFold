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

Citation policy: metric vs threshold
------------------------------------
The two are credited separately, and both are always credited (R009).

- **A metric belongs to whoever published it**, whatever threshold is applied to
  it later. ipSAE is Dunbrack's (bioRxiv 2025.02.10.637595); pDockQ is Bryant's;
  pDockQ2 is Zhu's; LIS is Kim's; the reference implementation is
  DunbrackLab/IPSAE `ipsae.py` v4. Adopting someone else's cutoff does not
  diminish authorship of the method, and the finding that the ipSAE paper states
  no cutoff is a statement about *thresholds only* - never a criticism of the
  work, and never a reason to under-credit it.
- **The ipSAE 0.6 cutoff and the four ipSAE confidence bands are AlphaFold DB's**
  (Han, Tsenkov, Venanzi et al., bioRxiv 10.64898/2026.03.27.714458v2, 3 Jul
  2026). They must not be attributed to Dunbrack.
- The EBI news item at `AFDB_NEWS_URL` is the accessible entry point to the same
  work and is worth linking as background, but it quotes high-confidence
  homodimers only and states no thresholds, so every number cites the preprint.

Two band schemes, one colour vocabulary
---------------------------------------
`ipsae_d0res` carries AlphaFold DB's four *published* confidence bands; the other
six values carry this project's own three-colour scheme. `traffic_light()` still
returns a 2-tuple whose first element is always `'green'` / `'amber'` / `'red'`
for every score - only the band *label* differs, and `Threshold.scheme`,
`.n_bands` and `.confidence_band()` expose the difference explicitly. See the
"Two band schemes coexist" comment above `Provenance`.

The classifier is `afdb_high_confidence()`
------------------------------------------
`ipSAE_d0res >= 0.6 AND pDockQ2 >= 0.23` is the joint rule AlphaFold DB actually
used to select which complexes to release. The seven traffic lights are
diagnostics that explain a verdict; they are not seven independent verdicts.

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
(the two `d0` helpers, `ptm_func`, the threshold table, the traffic light and the
AFDB joint criterion) are implemented here, because everything else depends on
them.
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
    "ConfidenceBand",
    "Threshold",
    "THRESHOLDS",
    "traffic_light",
    # the AFDB joint high-confidence criterion
    "AFDB_IPSAE_D0RES_MIN",
    "AFDB_PDOCKQ2_MIN",
    "AFDB_NEWS_URL",
    "AFDB_BELOW_THRESHOLD_FTP",
    "JointCriterion",
    "AFDB_JOINT_CRITERION",
    "AFDB_RELEASE_SCALE",
    "AFDBHighConfidence",
    "afdb_high_confidence",
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
# `specs/homodimer_diagnostic/threshold-reference.md` (R002, revised by R009).
# Do not edit the numbers here; edit that document and re-transcribe.
#
# Two band schemes coexist, deliberately (R009)
# ---------------------------------------------------------------------------
# AlphaFold DB publishes **four** named confidence bands for ipSAE, and this
# project's other six values have only a green/amber/red scheme invented here.
# Rather than flatten AFDB's published labels into three, or hand callers a
# tuple whose length changes per score, the two schemes are reconciled like
# this:
#
#   * The **colour vocabulary is invariant**: every score, under either scheme,
#     resolves to exactly one of `'green'` / `'amber'` / `'red'`, and
#     `traffic_light()` always returns a 2-tuple. A caller that only paints a
#     badge never has to know which scheme a score uses.
#   * The **label vocabulary is per-scheme**, and is the band's own name.
#     Six scores label their bands HIGH / MODERATE / LOW. `ipsae_d0res` labels
#     them with AFDB's published names, because inventing a synonym for a
#     published label would be worse than carrying it.
#   * `Threshold.scheme` and `Threshold.n_bands` say which scheme is in force,
#     and `Threshold.confidence_band()` returns the whole band record. Nothing
#     about the difference is implicit.
#
# The colour mapping for `ipsae_d0res` follows AFDB's own reading of its bands:
# "very high-confidence" and "confident" are both unambiguously good (green);
# "low-confidence" is AFDB's own caution (amber); below 0.6 the model does not
# qualify for high-confidence release at all (red). The red edge is therefore
# exactly the joint criterion's ipSAE cutoff, see `afdb_high_confidence`.

Provenance = Literal["PUBLISHED", "DERIVED", "HEURISTIC"]
"""How much a threshold is actually worth.

- ``PUBLISHED``  the paper explicitly states or applies this cutoff as an
                 operating criterion.
- ``DERIVED``    obtained from a paper's own quantitative statements by a
                 documented inferential step, but not stated by it as a
                 recommendation. The step is written out in `Threshold.source`.
- ``HEURISTIC``  no literature basis and no large-scale calibration. A
                 judgement call, and must be presented as one.

There is no ``AFDB-EMPIRICAL`` label any more. It existed for the ipSAE 0.6
cutoff while that number was an unpublished AlphaFold DB release filter; the
filter is now published (Han, Tsenkov, Venanzi et al. 2026), so the three ipSAE
entries moved to ``PUBLISHED`` / ``DERIVED`` and the label was retired rather
than left dead.
"""

Band = Literal["green", "amber", "red"]
"""Traffic-light colour. Invariant across both band schemes; the *label*
attached to a colour is not (see the section comment above)."""

_BAND_LABELS: Dict[str, str] = {
    "green": "HIGH",
    "amber": "MODERATE",
    "red": "LOW",
}
"""Default labels for the three-colour scheme used by six of the seven scores."""


@dataclass(frozen=True)
class ConfidenceBand:
    """
    One named band on a score's ladder, with the provenance of its lower edge.

    Attributes:
        key:        Stable machine identifier, e.g. `'very_high'`, `'high'`.
        label:      Display name. AFDB's published wording where one exists.
        colour:     `'green'` / `'amber'` / `'red'`.
        minimum:    Inclusive lower edge; `-inf` for the bottom band.
        provenance: Provenance of `minimum`.
        source:     Short citation for `minimum`, with a page number.
    """

    key: str
    label: str
    colour: Band
    minimum: float
    provenance: Provenance
    source: str = ""


@dataclass(frozen=True)
class Threshold:
    """
    One score's band ladder together with the provenance of every cutoff.

    Provenance is carried alongside the numbers deliberately: several values in
    the canonical table are judgement calls, and a user reading a coloured badge
    otherwise has no way to tell a published operating cutoff from a round
    number someone picked (R080).

    `green` and `amber` are the two **colour** edges and are always present, so
    every score can be rendered by a caller that knows nothing about named
    bands. `bands`, when non-empty, is the published named ladder; it is
    validated in `__post_init__` to agree with `green` and `amber`, so the two
    representations cannot drift apart.

    Attributes:
        green:            `value >= green` is green.
        amber:            `green > value >= amber` is amber; below is red.
        green_provenance: Provenance label for `green`.
        amber_provenance: Provenance label for `amber`.
        source:           Short citation for the green threshold, with a page.
        note:             One-line caveat to surface next to the badge.
        bands:            Published named ladder, ordered high to low. Empty
                          when the score has only the three-colour scheme.
    """

    green: float
    amber: float
    green_provenance: Provenance
    amber_provenance: Provenance
    source: str = ""
    note: str = ""
    bands: Tuple[ConfidenceBand, ...] = ()

    def __post_init__(self) -> None:
        if self.amber > self.green:
            raise ValueError(
                f"amber ({self.amber}) must not exceed green ({self.green})."
            )
        if not self.bands:
            return
        minima = [b.minimum for b in self.bands]
        if minima != sorted(minima, reverse=True):
            raise ValueError(
                "bands must be ordered from the highest band to the lowest; "
                f"got minima {minima}."
            )
        if self.bands[-1].minimum != float("-inf"):
            raise ValueError(
                "the lowest band must have a minimum of -inf so that every "
                "value falls in exactly one band."
            )
        greens = [b.minimum for b in self.bands if b.colour == "green"]
        ambers = [b.minimum for b in self.bands if b.colour == "amber"]
        if not greens:
            raise ValueError("a named ladder must contain at least one green band.")
        green_edge = min(greens)
        amber_edge = min(ambers) if ambers else green_edge
        if abs(green_edge - self.green) > 1e-12 or abs(amber_edge - self.amber) > 1e-12:
            raise ValueError(
                "named bands disagree with the colour edges: bands imply "
                f"green={green_edge}, amber={amber_edge}, but the Threshold "
                f"declares green={self.green}, amber={self.amber}."
            )

    # -- ladder -------------------------------------------------------------

    @property
    def scheme(self) -> str:
        """`'published_bands'` when a named ladder is carried, else `'traffic3'`."""
        return "published_bands" if self.bands else "traffic3"

    @property
    def ladder(self) -> Tuple[ConfidenceBand, ...]:
        """The band ladder, ordered high to low. Synthesised for `'traffic3'`."""
        if self.bands:
            return self.bands
        return (
            ConfidenceBand(
                key="high",
                label=_BAND_LABELS["green"],
                colour="green",
                minimum=self.green,
                provenance=self.green_provenance,
                source=self.source,
            ),
            ConfidenceBand(
                key="moderate",
                label=_BAND_LABELS["amber"],
                colour="amber",
                minimum=self.amber,
                provenance=self.amber_provenance,
                source=self.source,
            ),
            ConfidenceBand(
                key="low",
                label=_BAND_LABELS["red"],
                colour="red",
                minimum=float("-inf"),
                # The red band is defined by the amber edge, so it inherits that
                # edge's provenance.
                provenance=self.amber_provenance,
                source=self.source,
            ),
        )

    @property
    def n_bands(self) -> int:
        """Number of named bands. 3 for `'traffic3'`, 4 for `ipsae_d0res`."""
        return len(self.ladder)

    def confidence_band(self, value: float) -> ConfidenceBand:
        """Return the whole band record `value` falls in. Edges compare `>=`."""
        for band in self.ladder:
            if value >= band.minimum:
                return band
        raise AssertionError("unreachable: the lowest band has minimum -inf.")

    def band(self, value: float) -> Band:
        """Return the colour `value` falls in. Both edges compare with `>=`."""
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


# -- the AFDB high-confidence release criterion -----------------------------
# Han, Tsenkov, Venanzi et al., "AlphaFold Database expands to proteome-scale
# quaternary structures", bioRxiv 10.64898/2026.03.27.714458v2 (3 Jul 2026).
# Page numbers below are PDF pages of
# `specs/homodimer_diagnostic/references/afdb_dimers_2026.pdf`.

AFDB_IPSAE_D0RES_MIN: float = 0.60
"""ipSAE side of the AFDB joint high-confidence criterion (p. 5, p. 12)."""

AFDB_PDOCKQ2_MIN: float = 0.23
"""pDockQ2 side of the AFDB joint high-confidence criterion (p. 5, p. 12)."""

AFDB_NEWS_URL: str = (
    "https://www.ebi.ac.uk/about/news/technology-and-innovation/"
    "first-complexes-alphafold-database/"
)
"""EBI news item, 16 Mar 2026. Accessible background reading, not a technical
citation: it quotes 1.7M high-confidence *homodimers* only and states no
thresholds. Cite the preprint for every number."""

AFDB_BELOW_THRESHOLD_FTP: str = (
    "ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/"
)
"""Where models failing the joint criterion are still published, with their
interface scores (preprint p. 12). A FAIL verdict means "not surfaced as
high-confidence", not "not available"."""


@dataclass(frozen=True)
class JointCriterion:
    """The published AFDB high-confidence rule and its validation figures."""

    ipsae_d0res_min: float
    pdockq2_min: float
    source: str
    precision_homodimer: float
    fpr_homodimer: float
    precision_heterodimer: float
    fpr_heterodimer: float
    mcc_optimal_homodimer: float
    mcc_optimal_heterodimer: float
    note: str


AFDB_JOINT_CRITERION = JointCriterion(
    ipsae_d0res_min=AFDB_IPSAE_D0RES_MIN,
    pdockq2_min=AFDB_PDOCKQ2_MIN,
    source=(
        "Han, Tsenkov, Venanzi et al. 2026, p. 5: 'We adopted a combined "
        "high-confidence criterion requiring both community-established cutoffs "
        "of ipSAEmax >= 0.6 and pDockQ2max >= 0.23, corresponding to the DockQ "
        "‘acceptable’ quality boundary.' Restated as the release "
        "criterion on p. 12."
    ),
    precision_homodimer=0.924,
    fpr_homodimer=0.043,
    precision_heterodimer=0.958,
    fpr_heterodimer=0.004,
    mcc_optimal_homodimer=0.104,
    mcc_optimal_heterodimer=0.520,
    note=(
        "0.6 is a deliberately conservative operating point, not an optimum. "
        "The paper adopts it as 'community-established' by citation (p. 5) and "
        "derives only the validation: precision 0.924 (FPR 0.043) on the "
        "homodimer benchmark and 0.958 (FPR 0.004) on the heterodimer "
        "benchmark, 'supporting its use as a quality filter that prioritises "
        "precision over recall given the scale of the release' (p. 5). Its own "
        "MCC-optimal cutoffs for ipSAEmax are far lower: 0.104 on the homodimer "
        "and 0.520 on the heterodimer post-training benchmark (Supplementary "
        "Figs. 1-2, p. 24). A model below 0.6 is therefore 'not selected for "
        "high-confidence release', not 'not an interaction'."
    ),
)
"""The AFDB joint criterion, its validation, and the conservative-not-optimal
nuance, stored so R080 can surface it verbatim rather than paraphrase it."""


AFDB_RELEASE_SCALE: Dict[str, float] = {
    # p. 11: "we compiled over 31 million candidate complexes and surfaced 1.81
    # million high-confidence assemblies."
    "candidate_complexes": 31_000_000,
    # p. 3 and Supplementary Table 1, p. 29.
    "homodimers_analysed": 19_148_379,
    "heterodimers_analysed": 7_561_477,
    # p. 3, p. 5, p. 12.
    "homodimers_high_confidence": 1_735_475,
    "homodimer_high_confidence_fraction": 0.091,
    "heterodimers_high_confidence": 79_392,
    "heterodimer_high_confidence_fraction": 0.010,
    # p. 10 / p. 11: 1,814,832 high-confidence assemblies, quoted as 1.81M.
    "high_confidence_assemblies": 1_814_832,
    # p. 6: population of each published ipSAE band.
    "homodimers_very_high_confidence": 982_188,
    "heterodimers_very_high_confidence": 23_228,
    "homodimers_confident": 440_072,
    "heterodimers_confident": 31_728,
    "homodimers_low_confidence": 313_215,
    "heterodimers_low_confidence": 24_436,
}
"""Release-scale figures from the preprint, for prose that would otherwise
invent them. Every entry carries its page in the comment above it."""


_IPSAE_D0RES_BANDS: Tuple[ConfidenceBand, ...] = (
    ConfidenceBand(
        key="very_high",
        label="VERY HIGH-CONFIDENCE",
        colour="green",
        minimum=0.80,
        provenance="PUBLISHED",
        source="Han, Tsenkov, Venanzi et al. 2026, p. 12: 'very high-confidence' "
               "(ipSAEmax >= 0.8). Band populations on p. 6.",
    ),
    ConfidenceBand(
        key="confident",
        label="CONFIDENT",
        colour="green",
        minimum=0.70,
        provenance="PUBLISHED",
        source="Han, Tsenkov, Venanzi et al. 2026, p. 12: 'confident' "
               "(0.7 <= ipSAEmax < 0.8). Band populations on p. 6.",
    ),
    ConfidenceBand(
        key="low_confidence",
        label="LOW-CONFIDENCE",
        colour="amber",
        minimum=0.60,
        provenance="PUBLISHED",
        source="Han, Tsenkov, Venanzi et al. 2026, p. 12: 'low-confidence' "
               "(0.6 <= ipSAEmax < 0.7). This edge is also the ipSAE side of "
               "the joint release criterion (p. 5).",
    ),
    ConfidenceBand(
        key="below_threshold",
        label="BELOW AFDB THRESHOLD",
        colour="red",
        minimum=float("-inf"),
        provenance="PUBLISHED",
        source="Han, Tsenkov, Venanzi et al. 2026, p. 5 and p. 12: below "
               "ipSAEmax 0.6 a model is not surfaced as high-confidence. Such "
               "models remain downloadable with their interface scores (p. 12).",
    ),
)


THRESHOLDS: Dict[str, Threshold] = {
    "ipsae_d0res": Threshold(
        green=0.70,
        amber=0.60,
        green_provenance="PUBLISHED",
        amber_provenance="PUBLISHED",
        source="Han, Tsenkov, Venanzi et al. 2026, p. 12: AFDB's four published "
               "bands for ipSAEmax (>=0.8 very high-confidence, 0.7-0.8 "
               "confident, 0.6-0.7 low-confidence, <0.6 not surfaced). The "
               "*metric* is Dunbrack's (bioRxiv 2025.02.10.637595, cited as "
               "ref. 32 by the AFDB paper); the *thresholds* are AFDB's. Both "
               "must be credited.",
        note="Quote this variant, not d0chn or d0dom: AFDB's ipSAEmax is the "
             "per-residue-d0 score, maxed over the two chain directions. 0.6 is "
             "a conservative release cutoff, not an optimum (see "
             "AFDB_JOINT_CRITERION.note). For an AFDB accession the entry has "
             "already passed 0.6, so red will essentially never appear; for a "
             "locally uploaded model the full range is reachable.",
        bands=_IPSAE_D0RES_BANDS,
    ),
    "ipsae_d0chn": Threshold(
        green=0.70,
        amber=0.60,
        green_provenance="DERIVED",
        amber_provenance="DERIVED",
        source="Transferred from ipsae_d0res (Han, Tsenkov, Venanzi et al. 2026, "
               "pp. 5, 12). No publication states a cutoff for this variant, and "
               "AFDB calibrated on the d0res score, so the transfer is the one "
               "inferential step that makes this DERIVED rather than PUBLISHED.",
        note="ipSAE_d0chn >= ipSAE_d0dom >= ipSAE_d0res is a theorem, not an "
             "observation, so a shared threshold set is strictly most permissive "
             "here: this variant turns green first and red last. AFDB's band "
             "*names* are deliberately not reused for it - only the numbers "
             "transfer, not the published labels.",
    ),
    "ipsae_d0dom": Threshold(
        green=0.70,
        amber=0.60,
        green_provenance="DERIVED",
        amber_provenance="DERIVED",
        source="Transferred from ipsae_d0res (Han, Tsenkov, Venanzi et al. 2026, "
               "pp. 5, 12), same single inferential step as ipsae_d0chn.",
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
             "its own. Not part of the AFDB joint criterion; pDockQ2 is.",
    ),
    "pdockq2": Threshold(
        green=0.23,
        amber=0.10,
        green_provenance="PUBLISHED",
        amber_provenance="HEURISTIC",
        source="Double-sourced. Zhu 2023, p. 6: 'min pDockQ2 > 0.23' applied as "
               "an operating cutoff to select real biological conclusions. "
               "Independently, Han, Tsenkov, Venanzi et al. 2026, p. 5: "
               "'pDockQ2max >= 0.23, corresponding to the DockQ ‘acceptable’ "
               "quality boundary', used as the second half of the AFDB release "
               "criterion. The *metric* is Zhu's in both cases.",
        note="The sigmoid is bounded on [0.005, 1.315] and can exceed 1.0. Never "
             "present pDockQ2 as a probability or a percentage. AFDB's own "
             "MCC-optimal cutoff for pDockQ2max is 0.013 on the homodimer "
             "benchmark (Supplementary Fig. 1, p. 24), so 0.23 is conservative "
             "here too.",
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
"""Canonical thresholds for the seven confidence values.

Single source of truth. Transcribed from the Decision table in
`specs/homodimer_diagnostic/threshold-reference.md`; that document carries the
full derivation for every number and supersedes the four conflicting copies that
previously lived in `CLAUDE.md` and in notebook cells 25, 26 and 28.

**Amber is a stored number and must never be recomputed as `green / 2`** - that
recomputation is exactly why the summary table and the prose could disagree on a
single run.

Only `ipsae_d0res` carries a published named ladder (`scheme == 'published_bands'`,
four bands); the other six use the three-colour scheme. `traffic_light()` returns
a colour from the same three-value vocabulary either way.

Two judgement calls remain, both ambers (`pdockq2`, `lis`); R009 retired the three
invented ipSAE ambers by adopting AFDB's published bands.

The seven traffic lights are *diagnostics*, not the classifier. The classifier is
`afdb_high_confidence()`.
"""


def traffic_light(value: float, score_name: str) -> Tuple[Band, str]:
    """
    Classify a score against its canonical threshold.

    Args:
        value:      The computed score.
        score_name: Key into `THRESHOLDS`, e.g. `'ipsae_d0res'`, `'iptm_d0chn'`,
                    `'pdockq'`, `'pdockq2'`, `'lis'`.

    Returns:
        `(colour, label)`. `colour` is always one of `'green'` / `'amber'` /
        `'red'`, for every score, under either band scheme. `label` is the name
        of the band the value falls in: `'HIGH'` / `'MODERATE'` / `'LOW'` for the
        six scores on the three-colour scheme, and AFDB's own published wording
        for `ipsae_d0res`, which has four bands rather than three.

    Raises:
        KeyError: If `score_name` is not one of the seven canonical values. This
            is deliberate: a typo must fail loudly rather than silently defaulting
            to a band.

    Note:
        The return arity and the colour vocabulary never vary, so a caller that
        only paints a badge can ignore the scheme difference entirely. A caller
        that renders the label should read `THRESHOLDS[score_name].scheme` /
        `.n_bands`, or call `THRESHOLDS[score_name].confidence_band(value)` for
        the full record including the provenance and page citation of the edge
        (R080).

    Example
    -------
    >>> traffic_light(0.85, 'ipsae_d0res')
    ('green', 'VERY HIGH-CONFIDENCE')
    >>> traffic_light(0.75, 'ipsae_d0res')
    ('green', 'CONFIDENT')
    >>> traffic_light(0.65, 'ipsae_d0res')
    ('amber', 'LOW-CONFIDENCE')
    >>> traffic_light(0.55, 'ipsae_d0res')
    ('red', 'BELOW AFDB THRESHOLD')
    >>> traffic_light(0.75, 'iptm_d0chn')
    ('green', 'HIGH')
    >>> traffic_light(0.05, 'lis')
    ('red', 'LOW')
    """
    if score_name not in THRESHOLDS:
        raise KeyError(
            f"Unknown score name {score_name!r}. "
            f"Expected one of: {', '.join(sorted(THRESHOLDS))}."
        )
    band = THRESHOLDS[score_name].confidence_band(value)
    return band.colour, band.label


@dataclass(frozen=True)
class AFDBHighConfidence:
    """
    Verdict of the AFDB joint high-confidence criterion for one complex.

    Attributes:
        passed:           `True` when both sides pass.
        verdict:          `'PASS'` or `'FAIL'`, for display.
        ipsae_d0res:      The ipSAE_d0res value tested.
        pdockq2:          The pDockQ2 value tested.
        ipsae_d0res_pass: Whether the ipSAE side passed.
        pdockq2_pass:     Whether the pDockQ2 side passed.
        failed:           Names of the failing side(s), in table order. Empty on
                          a PASS; one entry when a single side failed; both when
                          neither passed.
        reason:           One-sentence explanation naming the values and cutoffs.
    """

    passed: bool
    verdict: str
    ipsae_d0res: float
    pdockq2: float
    ipsae_d0res_pass: bool
    pdockq2_pass: bool
    failed: Tuple[str, ...]
    reason: str


def afdb_high_confidence(ipsae_d0res: float, pdockq2: float) -> AFDBHighConfidence:
    """
    Apply AlphaFold DB's joint high-confidence release criterion.

    `ipSAE_d0res >= 0.6 AND pDockQ2 >= 0.23`. This single conjunction, not the
    seven per-score traffic lights, is the classifier AFDB actually used to
    decide which of ~31 million candidate complexes were surfaced. It answers
    one question: *would this model qualify for AFDB high-confidence release?*

    Both values must be the `max` over the two chain directions, matching the
    paper's `ipSAEmax` / `pDockQ2max` naming (p. 5) and `ipsae.py`'s own
    reporting; the paper's Supplementary Fig. 5 (p. 26) compares max- against
    min-aggregation and retains max.

    Args:
        ipsae_d0res: ipSAE with per-residue d0, maxed over both directions.
        pdockq2:     pDockQ2, maxed over both directions.

    Returns:
        An `AFDBHighConfidence` carrying the boolean, the per-side booleans, the
        names of any failing side(s), and a display-ready reason.

    Note:
        A FAIL means "not surfaced as an AFDB high-confidence entry", not "not
        an interaction". Below-threshold dimers remain published with their
        interface scores at `AFDB_BELOW_THRESHOLD_FTP` (p. 12), and 0.6 is a
        deliberately conservative operating point: see `AFDB_JOINT_CRITERION`
        for the precision/FPR figures and for the paper's own far lower
        MCC-optimal cutoffs (0.104 homodimer, 0.520 heterodimer, p. 24).

    Example
    -------
    >>> r = afdb_high_confidence(0.81, 0.64)
    >>> r.verdict, r.failed
    ('PASS', ())
    >>> afdb_high_confidence(0.55, 0.64).failed
    ('ipsae_d0res',)
    >>> afdb_high_confidence(0.81, 0.11).failed
    ('pdockq2',)
    >>> afdb_high_confidence(0.55, 0.11).failed
    ('ipsae_d0res', 'pdockq2')
    """
    ipsae_ok = ipsae_d0res >= AFDB_IPSAE_D0RES_MIN
    pdockq2_ok = pdockq2 >= AFDB_PDOCKQ2_MIN
    passed = ipsae_ok and pdockq2_ok

    failed = tuple(
        name
        for name, ok in (("ipsae_d0res", ipsae_ok), ("pdockq2", pdockq2_ok))
        if not ok
    )

    ipsae_clause = (
        f"ipSAE_d0res {ipsae_d0res:.3f} "
        f"{'>=' if ipsae_ok else '<'} {AFDB_IPSAE_D0RES_MIN:.2f}"
    )
    pdockq2_clause = (
        f"pDockQ2 {pdockq2:.3f} "
        f"{'>=' if pdockq2_ok else '<'} {AFDB_PDOCKQ2_MIN:.2f}"
    )
    if passed:
        reason = (
            f"PASS - meets AFDB's high-confidence release criterion "
            f"({ipsae_clause} and {pdockq2_clause})."
        )
    else:
        failing = " and ".join(
            clause
            for clause, ok in ((ipsae_clause, ipsae_ok), (pdockq2_clause, pdockq2_ok))
            if not ok
        )
        passing = " and ".join(
            clause
            for clause, ok in ((ipsae_clause, ipsae_ok), (pdockq2_clause, pdockq2_ok))
            if ok
        )
        tail = f"; {passing} passes." if passing else "; both sides fail."
        reason = (
            f"FAIL - would not be surfaced as an AFDB high-confidence entry: "
            f"{failing}{tail}"
        )

    return AFDBHighConfidence(
        passed=passed,
        verdict="PASS" if passed else "FAIL",
        ipsae_d0res=float(ipsae_d0res),
        pdockq2=float(pdockq2),
        ipsae_d0res_pass=ipsae_ok,
        pdockq2_pass=pdockq2_ok,
        failed=failed,
        reason=reason,
    )


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
