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
to import). Importing this module never mutates global matplotlib state: the
notebook's style block is applied only by an explicit `apply_plot_style()` call.
`molviewspec` is imported **lazily**, inside the function that needs
it, so the notebook degrades gracefully when it is not installed and so a plain
`import complex_interface_utils` never pays for it.

Prohibited anywhere in this module: `biopython` / `Bio`, `torch`,
`torch-geometric`, `gemmi`, `pandas`, `scipy`, `plotly`. The notebook must run on
Colab free tier within a 60 s install budget, which is also why the Colab
bootstrap clones the repo and extends `sys.path` instead of `pip install`-ing the
package (D9).

Status
------
This file is being filled in over several tasks; the sections still empty carry
a note naming the task that populates them. Implemented so far: AFDB access,
structure parsing, PAE / pLDDT parsing and interface detection (R011), the shared
scoring primitives (the two `d0` helpers and `ptm_func`), the seven score
functions (R012, landing R003 / R005 / R006 / R007), the threshold table with
its traffic light and the AFDB joint criterion, and the six matplotlib figures
(R013), and the four MolViewSpec 3D views with their shared display helper
(R014). Nothing is now left unimplemented.

The plotting section began as a behaviour-preserving move of the notebook's
inline figures (R013) and has since taken milestone M5: figure sizing and true
aspect ratios (R050), a sequential green default with switchable alternatives
through the `PAE_CMAP` seam (R051), and the rebuilt score-mask panel with its
shared colour bar and stated quadrant provenance (R052).

The MolViewSpec section is likewise a behaviour-preserving move: R070-R075 and
R030 redesign the views, and the four defects R075 lists that change no pixel
(the deprecated `cm.get_cmap`, one component per residue, hard-coded chain
letters, and the structure URL and format being resolved separately) are fixed
here while everything visual is left exactly as the notebook draws it.

The notebook still carries its own inline copies of the scoring code and still
runs off them; R016 switches the call sites over.
"""

from __future__ import annotations

import base64
import html as _html
import importlib
import json
import math
import subprocess
import sys
import textwrap
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import (
    Any, Dict, List, Literal, Mapping, NamedTuple, Optional, Sequence, Tuple,
)

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpatheffects
import matplotlib.pyplot as plt
import numpy as np
import requests
import seaborn as sns
from matplotlib.colors import Colormap
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator

__all__ = [
    # constants
    "PAE_CUTOFF",
    "DIST_CUTOFF",
    "LIS_CUTOFF",
    "AFDB_PREDICTION_URL",
    "DEFAULT_TIMEOUT",
    "BACKBONE_ATOMS",
    # AFDB access
    "AFDBPrediction",
    "fetch_afdb_metadata",
    "download_text",
    "download_json",
    "download_structure",
    "download_pae",
    "download_plddt",
    "AccessionLookupError",
    # local-file mode (R025)
    "REQUIRED_LOCAL_DOCUMENTS",
    "MissingLocalDocumentError",
    "require_local_documents",
    # assembly detection (R023)
    "SUPPORTED_N_CHAINS",
    "SUPPORTED_OLIGOMERIC_STATE",
    "AFDB_SEARCH_URL",
    "AFDB_DIMER_SEARCH",
    "EXAMPLE_DIMER_ACCESSION",
    "UnsupportedAssemblyError",
    "AssemblyDescription",
    "describe_assembly",
    "format_composition",
    "format_metadata_report",
    # structure parsing
    "ChainCoords",
    "parse_mmcif_atoms",
    "extract_chain_coords",
    "parse_structure",
    # PAE / pLDDT parsing
    "ChainSpan",
    "ChainPairPAE",
    "PAEMatrix",
    "parse_pae",
    "PLDDTScores",
    "parse_plddt",
    "verify_chain_lengths",
    "verify_document_agreement",
    # chain labelling (R021)
    "ChainLabel",
    "ChainIdentity",
    "verify_chain_identity",
    # interface detection
    "InterfaceContacts",
    "detect_interface",
    # scoring primitives
    "d0_scalar",
    "d0_array",
    "ptm_func",
    # score functions
    "ResidueProfile",
    "DirectionalPair",
    "compute_iptm_d0chn",
    "IPSAEResult",
    "compute_ipsae",
    "PDockQResult",
    "compute_pdockq",
    "PDockQ2Direction",
    "PDockQ2Result",
    "compute_pdockq2",
    "LISDirection",
    "LISResult",
    "compute_lis",
    # directional transparency (R008)
    "DIRECTIONAL_DELTA_TOLERANCE",
    "DirectionalDelta",
    "directional_deltas",
    "DIRECTION_FORWARD",
    "DIRECTION_REVERSE",
    "resolve_direction",
    "describe_direction",
    "format_directional_report",
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
    # Section 7: one summary object for the table and the prose (R080)
    "BAND_COLOUR_HEX",
    "SCORE_DESCRIPTIONS",
    "SummaryRow",
    "ConfidenceSummary",
    "summarise_scores",
    "format_afdb_verdict",
    "format_summary_table_html",
    "format_diagnostic_report",
    # plotting: style and palette
    "NOTEBOOK_RC_PARAMS",
    "apply_plot_style",
    "PAE_CMAP",
    "PAE_CMAP_CHOICES",
    "DIST_CMAP",
    "AGREEMENT_CMAP",
    "resolve_pae_cmap",
    # plotting: colours
    "COLOUR_A",
    "COLOUR_B",
    "COLOUR_IF",
    "COLOUR_NON_IF",
    "COLOUR_NON_IF_HIST",
    "CHAIN_COLOURS",
    "COLOUR_SIDE_CHAIN_X",
    "COLOUR_SIDE_CHAIN_Y",
    "SIDE_CHAIN_COLOURS",
    "SCORE_PROFILE_COLOURS",
    "PROFILE_SERIES",
    "PEAK_MARKER_COLOUR",
    "PLDDT_BANDS",
    "PLDDT_BAND_COLOURS",
    "plddt_band",
    "plddt_band_colour",
    "SCORE_DISPLAY_NAMES",
    "AGREEMENT_SCORES",
    # plotting: figures
    "score_masks",
    "plot_interface_contact_map",
    "plot_pae_matrix",
    "plot_pae_score_masks",
    "plot_residue_score_profiles",
    "plot_plddt_distribution",
    "plot_score_agreement",
    "THRESHOLD_MARGIN_CAP",
    "plot_threshold_margins",
    # MolViewSpec: capability check
    "MOLVIEWSPEC_MISSING_MESSAGE",
    "molviewspec_available",
    # MolViewSpec: view constants
    "MVS_VIEW_WIDTH",
    "MVS_VIEW_HEIGHT",
    "MVS_CONTEXT_COLOUR",
    "MVS_FAINT_COLOUR",
    "MVS_VALUE_CMAP",
    "MVS_PLDDT_BANDS",
    "MVS_DISAGREEMENT_THRESHOLD",
    "MVS_DISAGREEMENT_CATEGORIES",
    "MVS_CONTACT_PTM_THRESHOLD",
    "MVS_PDOCKQ2_AGREEMENT_CATEGORIES",
    "MVS_UNLIT_LABELS",
    "MVS_VIEW_LABELS",
    "format_view_label",
    # MolViewSpec: structure source
    "StructureSource",
    "resolve_structure_source",
    # MolViewSpec: legends
    "LegendEntry",
    "contact_ptm_to_pae",
    "format_category_labels",
    "category_legend",
    "plddt_legend",
    "value_ramp_legend",
    "chain_overview_legend",
    "legend_text",
    "legend_html",
    # MolViewSpec: per-residue colouring
    "ColourRun",
    "colour_runs",
    "value_colours",
    "add_residue_colours",
    # MolViewSpec: display and builders
    "mol_view_html",
    "show_mol_view",
    "build_chain_overview_view",
    "build_plddt_view",
    "build_interface_value_view",
    "build_category_view",
    "disagreement_categories",
    "build_disagreement_view",
    "disagreement_legend",
    "pdockq2_ipsae_categories",
    "build_pdockq2_agreement_view",
    "pdockq2_agreement_legend",
    # notebook orchestration (R100)
    "describe_checkout",
    "prepare_environment",
    "UploadPanel",
    "UPLOAD_INSTRUCTIONS",
    "local_upload_panel",
    "LOCAL_MODE_GATE_NOTE",
    "fetch_prediction",
    "SourceDocuments",
    "load_documents",
    "format_chain_lengths",
    "format_chain_report",
    "format_interface_report",
    "format_score_mask_usage",
    "format_iptm_report",
    "format_ipsae_report",
    "format_pdockq_report",
    "format_pdockq2_report",
    "format_lis_report",
    "format_score_values",
    "format_domain_sizes",
    "InterfacePLDDT",
    "interface_plddt_stats",
    "format_plddt_stats",
    "MOLVIEWSPEC_NO_URL_MESSAGE",
    "ViewRenderer",
    "view_renderer",
    "format_view_conventions",
    "format_interface_statistics",
    "format_references",
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
# AFDB access
# ---------------------------------------------------------------------------
# `requests`, not `urllib`: the PAE and pLDDT documents are served gzipped and
# `requests` decompresses them transparently, while `urllib.request.urlopen`
# hands back the raw deflate stream and `json.loads` then fails on byte 0x1f.

AFDB_PREDICTION_URL: str = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
"""AFDB prediction (metadata) endpoint. Returns a JSON array, one entry per chain."""

DEFAULT_TIMEOUT: float = 60.0
"""Seconds before a download is abandoned. Generous because Colab's egress is slow."""

_SHARED_DOCUMENT_FIELDS: Tuple[str, ...] = ("cifUrl", "bcifUrl", "paeDocUrl", "plddtDocUrl")
"""Metadata fields that address the *complex* and are therefore identical in
every per-chain entry. Read through `AFDBPrediction.document_url`, which checks
the agreement rather than assuming it."""

_UNSET: Any = object()
"""Sentinel distinguishing "no default given" from a default of `None`, so that
`shared_field('modelVersion')` can raise while `shared_field(..., default=None)`
returns `None`."""


@dataclass(frozen=True, eq=False)
class AFDBPrediction:
    """
    The whole prediction-endpoint response, with every per-chain entry preserved.

    The endpoint returns **one entry per chain**, and the order it returns them in
    is *non-deterministic*: the same accession answers `['A', 'B']` on one call and
    `['B', 'A']` on the next. Positional access is therefore not merely wrong for a
    heterodimer, it is wrong *intermittently* -- the UniProt accession, organism,
    protein name, gene and monomer length shown for an unchanged notebook could
    flip between chains from one run to the next. Two defences (R020):

    - `entries` is sorted by `chainId` at construction, so arrival order is thrown
      away before any caller can accidentally depend on it;
    - chain-specific fields are reachable only through `entry_for_chain` /
      `chain_field`, which raise `KeyError` for a chain the endpoint did not
      describe rather than quietly handing back another chain's identity.

    Whole-complex fields (the document URLs, the assembly type, the model version)
    go through `shared_field`, which checks that the entries agree instead of
    trusting whichever one happened to arrive first.

    Attributes:
        accession: The `AF-...` accession that was requested.
        entries:   Every entry the endpoint returned, sorted by `chainId`.

    Example
    -------
    >>> shuffled = AFDBPrediction('AF-1', ({'chainId': 'B', 'gene': 'Sumo1'},
    ...                                    {'chainId': 'A', 'gene': 'ISG20'}))
    >>> shuffled.chain_ids
    ('A', 'B')
    >>> shuffled.chain_field('A', 'geneNames', 'gene')
    'ISG20'
    >>> shuffled.entry_for_chain('C')
    Traceback (most recent call last):
    KeyError: "No metadata entry for chain 'C' in AF-1; the endpoint described chains A, B."
    """

    accession: str
    entries: Tuple[Dict[str, Any], ...]

    def __post_init__(self) -> None:
        # Sorting here, rather than at each call site, is what makes "never code
        # against the endpoint's order" structural instead of a convention.
        object.__setattr__(
            self,
            "entries",
            tuple(sorted(self.entries, key=lambda entry: str(entry.get("chainId", "")))),
        )

    @property
    def chain_ids(self) -> Tuple[str, ...]:
        """`chainId` of each entry, sorted -- never the endpoint's arrival order."""
        return tuple(str(entry.get("chainId", "")) for entry in self.entries)

    def describes_chain(self, chain_id: str) -> bool:
        """
        Whether the endpoint returned an entry for this chain.

        The structure may carry chains the metadata does not describe (and the
        reverse), so a caller that wants to report the gap rather than fail on it
        asks this first.
        """
        return chain_id in self.chain_ids

    def entry_for_chain(self, chain_id: str) -> Dict[str, Any]:
        """
        The metadata entry describing one chain.

        Args:
            chain_id: Structure chain label, e.g. `'A'`.

        Returns:
            That chain's entry.

        Raises:
            KeyError: If no entry carries that `chainId`. Failing is the point:
                the alternative is reporting some other chain's identity, or empty
                fields, as if they belonged to this one.
        """
        for entry in self.entries:
            if str(entry.get("chainId", "")) == chain_id:
                return entry
        described = ", ".join(self.chain_ids) or "(none)"
        raise KeyError(
            f"No metadata entry for chain {chain_id!r} in {self.accession}; "
            f"the endpoint described chains {described}."
        )

    def chain_field(self, chain_id: str, *names: str, default: Any = None) -> Any:
        """
        One chain's value for the first of `names` that entry actually carries.

        The fallback chain exists because the endpoint's field names are not
        stable across records: `proteinFullName` is documented but the live
        service sends `uniprotDescription`, and `geneNames` arrives as `gene`.

        Args:
            chain_id: Structure chain label.
            names:    Field names, tried in order.
            default:  Returned when the entry carries none of them.

        Returns:
            The value, or `default`.

        Raises:
            KeyError: If the endpoint described no such chain.
        """
        entry = self.entry_for_chain(chain_id)
        for name in names:
            value = entry.get(name)
            if value is not None and value != "":
                return value
        return default

    def shared_field(self, *names: str, default: Any = _UNSET) -> Any:
        """
        A whole-complex value, checked for agreement across the per-chain entries.

        Args:
            names:   Field names, tried in order within each entry.
            default: Returned when no entry carries any of them. Omit it to raise
                     `KeyError` instead.

        Returns:
            The single agreed value.

        Raises:
            ValueError: If the entries disagree, which means the field is
                per-chain after all and reporting it once would be a lie.
            KeyError: If absent from every entry and no `default` was given.

        Example
        -------
        >>> pred = AFDBPrediction('AF-1', ({'chainId': 'A', 'assemblyType': 'Homo'},
        ...                                {'chainId': 'B', 'assemblyType': 'Homo'}))
        >>> pred.shared_field('assemblyType')
        'Homo'
        >>> pred.shared_field('oligomericState', default='N/A')
        'N/A'
        """
        found: Dict[str, Any] = {}
        for entry in self.entries:
            for name in names:
                value = entry.get(name)
                if value is not None and value != "":
                    found.setdefault(str(value), value)
                    break
        if not found:
            if default is _UNSET:
                raise KeyError(
                    f"{' / '.join(names)} is absent from every entry of {self.accession}."
                )
            return default
        if len(found) > 1:
            raise ValueError(
                f"{' / '.join(names)} differs between chains of {self.accession}: "
                f"{sorted(found)}. It was assumed to describe the whole complex."
            )
        return next(iter(found.values()))

    def document_url(self, field: str) -> str:
        """
        A whole-complex document URL, checked for agreement across entries.

        Args:
            field: One of `'cifUrl'`, `'bcifUrl'`, `'paeDocUrl'`, `'plddtDocUrl'`.

        Returns:
            The URL.

        Raises:
            ValueError: If the entries disagree, which would mean the endpoint
                serves per-chain documents and the whole download path needs
                rethinking. Better to fail loudly than to silently download one
                chain's document and slice it as if it covered the complex.
            KeyError: If the field is absent from every entry.
        """
        return str(self.shared_field(field))

    @property
    def cif_url(self) -> str:
        """mmCIF download URL (coordinates)."""
        return self.document_url("cifUrl")

    @property
    def bcif_url(self) -> str:
        """BinaryCIF download URL (MolViewSpec only; never parsed here)."""
        return self.document_url("bcifUrl")

    @property
    def pae_url(self) -> str:
        """PAE JSON download URL."""
        return self.document_url("paeDocUrl")

    @property
    def plddt_url(self) -> str:
        """pLDDT JSON download URL."""
        return self.document_url("plddtDocUrl")


def download_text(url: str, timeout: float = DEFAULT_TIMEOUT) -> str:
    """
    GET a URL and return its decoded body.

    Args:
        url:     Absolute URL.
        timeout: Seconds.

    Returns:
        The response body as text.

    Raises:
        requests.HTTPError: On a non-2xx status.
    """
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    return response.text


def download_json(url: str, timeout: float = DEFAULT_TIMEOUT) -> Any:
    """
    GET a URL and parse its body as JSON.

    Used for the PAE and pLDDT documents, which AFDB serves gzip-encoded;
    `requests` decompresses them transparently on the basis of the
    `Content-Encoding` header, which is why this module does not use `urllib`.

    Args:
        url:     Absolute URL.
        timeout: Seconds.

    Returns:
        The parsed JSON: an array for PAE, an object for pLDDT.

    Raises:
        requests.HTTPError: On a non-2xx status.
    """
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    return response.json()


def fetch_afdb_metadata(accession: str, timeout: float = DEFAULT_TIMEOUT) -> AFDBPrediction:
    """
    Fetch the AFDB prediction metadata for one accession.

    Args:
        accession: An AFDB accession such as `'AF-0000000065889468'`.
        timeout:   Seconds.

    Returns:
        An `AFDBPrediction` holding *all* per-chain entries.

    Raises:
        AccessionLookupError: If the accession is malformed, unknown, or the
            endpoint is unreachable. A `ValueError` subclass. This is the *only*
            failure mode: `requests.HTTPError` is deliberately not allowed to
            escape (R023), because `400 Client Error: Bad Request for url: ...`
            tells a user pasting an accession into a Run-all notebook neither
            what a valid accession looks like nor what to do next.

    Example
    -------
    The message names the input, quotes the service, and says what to try:

    >>> try:                                       # doctest: +SKIP
    ...     fetch_afdb_metadata('AF-NOT_A_REAL_ACCESSION')
    ... except AccessionLookupError as exc:
    ...     print(str(exc).splitlines()[0])
    AFDB rejected the accession 'AF-NOT_A_REAL_ACCESSION' as malformed (HTTP 400).
    """
    url = AFDB_PREDICTION_URL.format(accession=accession)
    try:
        response = requests.get(url, timeout=timeout)
    except requests.RequestException as exc:
        raise AccessionLookupError(
            f"Could not reach the AFDB prediction endpoint for {accession!r}.\n"
            f"  found     : the request to {url}\n"
            f"              failed with {type(exc).__name__}: {exc}\n"
            f"  what to do: this is a network problem, not a problem with the "
            f"accession.\n"
            f"              Check the connection and re-run this cell. On Colab, "
            f"re-running\n"
            f"              the cell is usually enough; the notebook keeps no "
            f"partial state."
        ) from exc

    if not response.ok:
        said = _service_error_text(response)
        if response.status_code == 404:
            headline = (f"AFDB has no prediction for {accession!r} "
                        f"(HTTP 404).")
            found = ("  found     : the accession is well formed, but the database "
                     "holds no model\n              under it.")
        elif response.status_code == 400:
            headline = (f"AFDB rejected the accession {accession!r} as malformed "
                        f"(HTTP {response.status_code}).")
            found = ("  found     : the endpoint refused the identifier before "
                     "looking anything up.")
        else:
            headline = (f"AFDB could not return a prediction for {accession!r} "
                        f"(HTTP {response.status_code}).")
            found = ("  found     : the endpoint answered with an error status.")
        if said:
            found += f"\n              the service said: {said}"
        raise AccessionLookupError(
            f"{headline}\n{found}\n{_accession_advice(accession)}")

    try:
        payload = response.json()
    except ValueError as exc:
        raise AccessionLookupError(
            f"AFDB answered for {accession!r} with something that is not JSON.\n"
            f"  found     : {response.text[:120]!r}\n"
            f"{_accession_advice(accession)}"
        ) from exc

    entries = payload if isinstance(payload, list) else [payload]
    entries = [entry for entry in entries if isinstance(entry, dict)]
    if not entries:
        raise AccessionLookupError(
            f"AFDB returned an empty result for {accession!r}.\n"
            f"  found     : the endpoint answered 200 but described no model, "
            f"so there is\n              nothing to download or score.\n"
            f"{_accession_advice(accession)}")
    return AFDBPrediction(accession=accession, entries=tuple(entries))


def download_structure(prediction: AFDBPrediction, timeout: float = DEFAULT_TIMEOUT) -> str:
    """Download the mmCIF text for a prediction. See `download_text`."""
    return download_text(prediction.cif_url, timeout=timeout)


def download_pae(prediction: AFDBPrediction, timeout: float = DEFAULT_TIMEOUT) -> Any:
    """Download the raw PAE JSON for a prediction. Feed it to `parse_pae`."""
    return download_json(prediction.pae_url, timeout=timeout)


def download_plddt(prediction: AFDBPrediction, timeout: float = DEFAULT_TIMEOUT) -> Any:
    """Download the raw pLDDT JSON for a prediction. Feed it to `parse_plddt`."""
    return download_json(prediction.plddt_url, timeout=timeout)


# --- local-file mode (R025) --------------------------------------------------
# `USE_LOCAL_FILE = True` bypasses this whole section, so the notebook has to be
# told separately what a *complete* upload looks like. It used to announce that a
# missing PAE or pLDDT file would make it "skip PAE-dependent analyses", and then
# skip nothing: `pae_raw` stayed `None`, and two cells later `parse_pae(None)`
# raised `TypeError: 'NoneType' object is not subscriptable` -- a message with no
# hint that an upload was the cause, from a cell the user was not looking at.
#
# The advertised skip is not worth building, because there is nothing left to
# skip *to*. Six of the seven values come off the PAE matrix (ipSAE d0res, d0chn
# and d0dom, ipTM_d0chn, pDockQ2, LIS); the seventh, pDockQ, needs per-residue
# pLDDT. Drop those two documents and Sections 3, 4, 6 and 7 are empty and the
# traffic light -- the notebook's entire output -- has nothing to colour. What
# survives is a contact count and a contact map: a different, much smaller
# notebook, not this one with a banner on it.
#
# That matters more here than the crash does. `threshold-reference.md` s1 rests
# on the local-file path being the one place the traffic light does real work:
# every AFDB accession has already passed the 0.6 ipSAE filter, so red is
# essentially unreachable online. The single route that exercises the full range
# was the route that crashed.
#
# So: refuse at the point of upload, while the user is still looking at the
# widget, and name every missing file at once with somewhere to get it.

REQUIRED_LOCAL_DOCUMENTS: Tuple[Tuple[str, str, str], ...] = (
    ("structure", "mmCIF structure", "AF-<id>-model_v<n>.cif"),
    ("PAE", "PAE document", "AF-<id>-predicted_aligned_error_v<n>.json"),
    ("pLDDT", "pLDDT document", "AF-<id>-confidence_v<n>.json"),
)
"""The three uploads local-file mode requires: `(short name, description, AFDB filename)`.

All three are required, not two-plus-optionals. See the comment above for why a
partial upload is refused rather than partially analysed."""


class MissingLocalDocumentError(ValueError):
    """
    Local-file mode was run without all three of mmCIF, PAE and pLDDT.

    A `ValueError` subclass, like `AccessionLookupError` and
    `UnsupportedAssemblyError`, so a caller already catching `ValueError` keeps
    working. Raised at upload time rather than allowing a `TypeError` from
    `parse_pae(None)` further down, and rather than the skip the notebook used
    to promise and never perform (R025).
    """


def _uploaded(value: Any) -> bool:
    """True when an upload slot holds something. Empty text or bytes do not count."""
    if value is None:
        return False
    if isinstance(value, (str, bytes, bytearray, list, tuple, dict)):
        return len(value) > 0
    return True


def require_local_documents(
    structure: Any,
    pae: Any,
    plddt: Any,
    accession: str = "",
) -> None:
    """
    Refuse a local-file run that is missing any of the three documents.

    Every missing document is named in one message, so a user who uploaded only
    the mmCIF is not sent round the loop three times.

    Args:
        structure: The uploaded mmCIF text or bytes; `None` or empty if absent.
        pae:       The uploaded PAE document, raw or parsed; `None` if absent.
        plddt:     The uploaded pLDDT document, raw or parsed; `None` if absent.
        accession: `ACCESSION_ID`, quoted in the advice so that the AFDB URLs
                   which would supply the missing files are copy-pasteable.

    Returns:
        `None`, when all three are present.

    Raises:
        MissingLocalDocumentError: If any is missing, naming which, why it is
            needed and where to get it.

    Example
    -------
    A complete upload passes silently:

    >>> require_local_documents('data_AF...', [{'predicted_aligned_error': [[0.0]]}],
    ...                         {'confidenceScore': [90.0]})

    A partial one names everything that is missing, at once:

    >>> try:
    ...     require_local_documents('data_AF...', None, None, accession='AF-123')
    ... except MissingLocalDocumentError as exc:
    ...     print(chr(10).join(str(exc).splitlines()[:3]))
    Local-file mode needs all three documents; 2 were not uploaded.
      uploaded  : mmCIF structure
      missing   : PAE document, pLDDT document

    An empty file is not an upload:

    >>> try:
    ...     require_local_documents('', b'[]', {})
    ... except MissingLocalDocumentError as exc:
    ...     print(str(exc).splitlines()[2])
      missing   : mmCIF structure, pLDDT document
    """
    present = {"structure": _uploaded(structure),
               "PAE": _uploaded(pae),
               "pLDDT": _uploaded(plddt)}
    missing = [desc for key, desc, _ in REQUIRED_LOCAL_DOCUMENTS if not present[key]]
    if not missing:
        return

    uploaded = [desc for key, desc, _ in REQUIRED_LOCAL_DOCUMENTS if present[key]]
    filenames = [f"{name}" for key, _, name in REQUIRED_LOCAL_DOCUMENTS
                 if not present[key]]
    count = ("one was" if len(missing) == 1 else f"{len(missing)} were")
    metadata_url = AFDB_PREDICTION_URL.format(accession=accession or "<accession>")

    raise MissingLocalDocumentError(
        f"Local-file mode needs all three documents; {count} not uploaded.\n"
        f"  uploaded  : {_wrap_inline(uploaded, 14) if uploaded else '(nothing)'}\n"
        f"  missing   : {_wrap_inline(missing, 14)}\n"
        "  why       : six of the seven values are computed from the PAE matrix\n"
        "              (ipSAE d0res/d0chn/d0dom, ipTM_d0chn, pDockQ2, LIS) and the\n"
        "              seventh, pDockQ, needs per-residue pLDDT. Without them\n"
        "              Sections 3, 4, 6 and 7 have nothing to show and the traffic\n"
        "              light has nothing to colour, so there is no partial run\n"
        "              worth offering -- only a contact count and a contact map.\n"
        f"  what to do: upload the missing file{'' if len(missing) == 1 else 's'} "
        "in the widget above, then re-run this\n"
        "              cell. For an AFDB model, the download links are the\n"
        "              `cifUrl`, `paeDocUrl` and `plddtDocUrl` fields of\n"
        f"                {metadata_url}\n"
        "              and the files are named\n"
        + "".join(f"                {name}\n" for name in filenames) +
        "              For a model of your own, export the PAE matrix and the\n"
        "              per-residue pLDDT in the AFDB JSON layout that `parse_pae`\n"
        "              and `parse_plddt` document.\n"
        "              To analyse an AFDB entry instead, set USE_LOCAL_FILE = "
        "False\n"
        "              and put its accession in ACCESSION_ID."
    )



# --- assembly detection (R023) ----------------------------------------------
# Two independent kinds of evidence describe what a prediction *is*:
#
#   declared -- `assemblyType`, `oligomericState`, `complexComposition` and
#               `isComplex`, asserted by the AFDB metadata endpoint;
#   observed -- the chains actually present in the mmCIF / PAE / pLDDT documents
#               and the UniProt accession each of them resolves to.
#
# They are produced by different parts of the pipeline and can disagree, so this
# module reconciles them and *shows both* rather than picking a winner behind the
# reader's back. Nothing downstream branches on the answer -- every score is
# computed from the chains themselves, D3 -- so a disagreement is a reporting
# problem, not a scoring one, and is reported as such.
#
# The dimer restriction (D4) is enforced here too, because a valid prediction the
# notebook cannot analyse should say so in one sentence rather than crash eight
# cells later on an index that does not exist.

SUPPORTED_N_CHAINS: int = 2
"""Number of chains this notebook analyses. D4: it scores one ordered chain pair."""

SUPPORTED_OLIGOMERIC_STATE: str = "dimer"
"""`oligomericState` value corresponding to `SUPPORTED_N_CHAINS`."""

AFDB_SEARCH_URL: str = "https://alphafold.ebi.ac.uk/api/search"
"""AFDB search endpoint, used to point a user at an accession that *will* work.

Undocumented in `CLAUDE.md`; R093 should add it. Measured 2026-09-08:
`isComplex:true` returns 2 010 763 hits and `oligomericState:dimer` returns the
same 2 010 763, while `oligomericState:trimer` and `oligomericState:tetramer`
return zero. **Every complex in AFDB is a dimer**, which is why the >2-chain
path below can only be reached by a local mmCIF."""

AFDB_DIMER_SEARCH: str = (
    f"{AFDB_SEARCH_URL}?q=oligomericState:dimer&type=main&rows=5"
)
"""A ready-made search for accessions this notebook does support."""

EXAMPLE_DIMER_ACCESSION: str = "AF-0000000065889468"
"""A homodimer, quoted in error messages so the fix is copy-pasteable."""

_OLIGOMERIC_NAMES: Dict[int, str] = {
    1: "monomer", 2: "dimer", 3: "trimer", 4: "tetramer",
    5: "pentamer", 6: "hexamer", 7: "heptamer", 8: "octamer",
}
"""Chain count to the word AFDB's `oligomericState` would use for it."""


class AccessionLookupError(ValueError):
    """
    The accession could not be turned into a prediction.

    A `ValueError` subclass so that a caller already catching `ValueError` keeps
    working. Raised instead of letting `requests.HTTPError` reach the user,
    whose text (`400 Client Error: Bad Request for url: ...`) says nothing about
    what a valid accession looks like or what to do next.
    """


class UnsupportedAssemblyError(ValueError):
    """
    The prediction is real and readable, but not something this notebook scores.

    A monomer, or anything with more than two chains. Every metric here is an
    *inter-chain* measurement, so the honest response is a refusal that names the
    limit, not a zero or an `IndexError`.
    """


def _accession_advice(accession: str) -> str:
    """The shared 'what a valid accession looks like, and what to do' footer."""
    return (
        "  expected  : an AFDB model accession -- "
        f"'{EXAMPLE_DIMER_ACCESSION}' (a two-chain complex),\n"
        "              'AF-P0A6Q3-F1' (one UniProt entry), or a bare UniProt "
        "accession\n"
        "              such as 'P0A6Q3'.\n"
        f"  what to do: check ACCESSION_ID = {accession!r} in the input cell for a "
        "typo. To find\n"
        "              an accession this notebook supports:\n"
        f"                {AFDB_DIMER_SEARCH}\n"
        "              To analyse a structure of your own instead, set "
        "USE_LOCAL_FILE = True."
    )


def _service_error_text(response: "requests.Response") -> str:
    """The service's own `error` string, when it sent one, else `''`."""
    try:
        body = response.json()
    except ValueError:
        return ""
    if isinstance(body, dict):
        return str(body.get("error") or body.get("message") or "")
    return ""


def _wrap_inline(parts: Sequence[str], indent: int, width: int = 74) -> str:
    """A comma-joined list, wrapped onto continuation lines at `indent` spaces."""
    lines: List[str] = []
    current = ""
    for i, part in enumerate(parts):
        piece = part + ("," if i < len(parts) - 1 else "")
        if current and indent + len(current) + 1 + len(piece) > width:
            lines.append(current)
            current = piece
        else:
            current = f"{current} {piece}".strip()
    if current:
        lines.append(current)
    return ("\n" + " " * indent).join(lines)


def _normalise_composition(raw: Any) -> Tuple[Tuple[str, int], ...]:
    """
    `complexComposition` as a sorted `((identifier, stoichiometry), ...)` tuple.

    Sorted so that two compositions listing the same proteins in different orders
    compare equal -- the endpoint's ordering is no more trustworthy here than its
    entry ordering was (R020).

    Example
    -------
    >>> _normalise_composition([{'identifierType': 'uniprotAccession',
    ...                          'identifier': 'P63166', 'stoichiometry': 1},
    ...                         {'identifierType': 'uniprotAccession',
    ...                          'identifier': 'Q96AZ6', 'stoichiometry': 1}])
    (('P63166', 1), ('Q96AZ6', 1))
    >>> _normalise_composition(None)
    ()
    """
    if not isinstance(raw, (list, tuple)):
        return ()
    items: Dict[str, int] = {}
    for member in raw:
        if not isinstance(member, dict):
            continue
        identifier = str(member.get("identifier") or "").strip()
        if not identifier:
            continue
        try:
            count = int(member.get("stoichiometry") or 1)
        except (TypeError, ValueError):
            count = 1
        items[identifier] = items.get(identifier, 0) + count
    return tuple(sorted(items.items()))


def format_composition(composition: Sequence[Tuple[str, int]]) -> str:
    """
    A composition as a reader's phrase: `'P0A6Q3 x2'`, `'Q96AZ6 x1 + P63166 x1'`.

    Args:
        composition: `((identifier, stoichiometry), ...)`.

    Returns:
        The phrase, or `'unknown'` when the composition is empty.

    Example
    -------
    >>> format_composition((('P0A6Q3', 2),))
    'P0A6Q3 x2'
    >>> format_composition((('P63166', 1), ('Q96AZ6', 1)))
    'P63166 x1 + Q96AZ6 x1'
    >>> format_composition(())
    'unknown'
    """
    if not composition:
        return "unknown"
    return " + ".join(f"{identifier} x{count}" for identifier, count in composition)


def _declared_field(
    prediction: Optional[AFDBPrediction], *names: str
) -> Tuple[Any, Optional[str]]:
    """
    A whole-complex metadata value, plus a conflict string if the entries differ.

    `AFDBPrediction.shared_field` *raises* when the per-chain entries disagree,
    which is right for a document URL (a disagreement there breaks the download)
    but wrong for an assembly label (a disagreement there is exactly the thing
    this task exists to surface). So this reports instead of raising.
    """
    if prediction is None:
        return None, None
    seen: Dict[str, Tuple[Any, List[str]]] = {}
    for chain_id in prediction.chain_ids:
        value = prediction.chain_field(chain_id, *names, default=None)
        if value is None or value == "":
            continue
        key = repr(value)
        seen.setdefault(key, (value, []))[1].append(chain_id)
    if not seen:
        return None, None
    if len(seen) == 1:
        return next(iter(seen.values()))[0], None
    detail = "; ".join(
        f"{', '.join(chains)} say {value!r}" for value, chains in sorted(
            seen.values(), key=lambda item: item[1])
    )
    return None, (
        f"the per-chain metadata entries disagree about "
        f"{' / '.join(names)}: {detail}. Neither is taken as the complex's."
    )


_ASSEMBLY_ROW_WIDTH: int = 76
"""Wrap width for the assembly block. Wider than `_REPORT_WIDTH`'s rule because
`complexComposition` can carry two 10-character accessions and a URL must not be
broken across lines."""


def _wrap_row(label: str, value: str, width: int = _ASSEMBLY_ROW_WIDTH,
              label_width: int = 14) -> List[str]:
    """
    `'Label         : value'`, wrapped onto hanging-indented continuation lines.

    URLs and accessions are never broken mid-token, so a wrapped line can still
    be copied and pasted.

    Example
    -------
    >>> for line in _wrap_row('Assembly', 'alpha beta gamma delta', width=36):
    ...     print(line)
    Assembly      : alpha beta gamma
                    delta
    """
    head = f"{label:<{label_width}}: "
    pieces = textwrap.wrap(value, width=max(width - len(head), 20),
                           break_long_words=False, break_on_hyphens=False) or [""]
    return [head + pieces[0]] + [" " * len(head) + piece for piece in pieces[1:]]


@dataclass(frozen=True)
class AssemblyDescription:
    """
    What this prediction is, from both kinds of evidence, with the gaps named.

    Built by `describe_assembly`. Holds the AFDB metadata's *declaration* and the
    chains' own *observation* side by side, never collapsed into one number:
    they are produced independently and, when they disagree, the reader is the
    only one who can judge which to believe. `conflicts` is the list of those
    disagreements, in plain English.

    Attributes:
        accession:            The accession this describes.
        chain_ids:            The chains, sorted. Observed from the structure and
                              documents when `chains_observed`, otherwise the
                              chains the metadata endpoint happened to describe.
        chain_lengths:        `{chain_id: n_residues}`, empty when not yet parsed.
        chains_observed:      Whether `chain_ids` comes from the structure/PAE
                              (authoritative) or from the metadata alone (not).
        observed_composition: `((uniprot, count), ...)` resolved from the chains
                              themselves; `()` when a chain's identity is unknown.
        declared_type:        `assemblyType`: `'Homo'`, `'Hetero'` or `None`.
        declared_state:       `oligomericState`, e.g. `'dimer'`, or `None`.
        declared_description: `oligomericStateDescription`, e.g. `'Heterodimer'`.
        declared_composition: `complexComposition`, normalised and sorted.
        declared_is_complex:  `isComplex`, or `None` when absent.
        conflicts:            Disagreements between the two, as sentences.

    Example
    -------
    >>> homo = AFDBPrediction('AF-1', (
    ...     {'chainId': 'A', 'uniprotAccession': 'P0A6Q3', 'assemblyType': 'Homo',
    ...      'oligomericState': 'dimer', 'isComplex': True,
    ...      'complexComposition': [{'identifier': 'P0A6Q3', 'stoichiometry': 2}]},
    ...     {'chainId': 'B', 'uniprotAccession': 'P0A6Q3', 'assemblyType': 'Homo',
    ...      'oligomericState': 'dimer', 'isComplex': True,
    ...      'complexComposition': [{'identifier': 'P0A6Q3', 'stoichiometry': 2}]}))
    >>> assembly = describe_assembly(homo, {'A': 172, 'B': 172})
    >>> assembly.noun
    'homodimer'
    >>> assembly.headline
    'Homodimer -- 2 chains (A, B), one protein in 2 copies: P0A6Q3 x2'
    >>> assembly.conflicts
    ()
    """

    accession: str
    chain_ids: Tuple[str, ...] = ()
    chain_lengths: Dict[str, int] = field(default_factory=dict)
    chains_observed: bool = False
    observed_composition: Tuple[Tuple[str, int], ...] = ()
    declared_type: Optional[str] = None
    declared_state: Optional[str] = None
    declared_description: Optional[str] = None
    declared_composition: Tuple[Tuple[str, int], ...] = ()
    declared_is_complex: Optional[bool] = None
    conflicts: Tuple[str, ...] = ()
    metadata_present: bool = False

    @property
    def n_chains(self) -> int:
        """How many chains there are. Trustworthy only when `chains_observed`."""
        return len(self.chain_ids)

    @property
    def n_distinct_proteins(self) -> int:
        """Distinct proteins among the chains; `0` when their identity is unknown."""
        return len(self.observed_composition)

    @property
    def observed_type(self) -> Optional[str]:
        """
        `'Homo'`, `'Hetero'` or `None` -- from the chains, not from the metadata.

        `None` means at least one chain could not be resolved to a protein, which
        is the local-file case; it is *not* the same as "the chains differ".
        """
        if not self.observed_composition:
            return None
        return "Homo" if self.n_distinct_proteins == 1 else "Hetero"

    @property
    def observed_state(self) -> Optional[str]:
        """`'dimer'`, `'trimer'`, ... from the chain count, or `None` if unobserved."""
        if not self.chains_observed:
            return None
        return _OLIGOMERIC_NAMES.get(self.n_chains, f"{self.n_chains}-mer")

    @property
    def assembly_type(self) -> Optional[str]:
        """
        The homo/hetero call, observation first.

        The chains are what every score is actually computed from, so when the
        two sources disagree the chains are what the notebook's *labels* follow.
        The disagreement is never hidden: it is in `conflicts` and printed.

        When the chains *are* observed but their proteins are not identifiable,
        this is `None` rather than the declared value: an unverifiable claim is
        downgraded to "dimer", not repeated as if the chains had confirmed it.
        The declaration still appears verbatim in `declared_phrase`.
        """
        if self.observed_type:
            return self.observed_type
        if not self.chains_observed:
            return self.declared_type
        return None

    @property
    def noun(self) -> str:
        """
        What to call this thing: `'homodimer'`, `'heterodimer'`, `'monomer'`,
        `'3-chain complex'`.

        Example
        -------
        >>> base = dict(accession='AF-1', chains_observed=True)
        >>> AssemblyDescription(chain_ids=('A',), **base).noun
        'monomer'
        >>> AssemblyDescription(chain_ids=('A', 'B'),
        ...                     observed_composition=(('P1', 1), ('P2', 1)),
        ...                     **base).noun
        'heterodimer'
        >>> AssemblyDescription(chain_ids=('A', 'B', 'C'), **base).noun
        'trimer'
        >>> AssemblyDescription(chain_ids=tuple('ABCDEFGHI'), **base).noun
        '9-chain complex'
        """
        state = self.observed_state or self.declared_state
        if not state and self.declared_is_complex is False:
            # AFDB omits `oligomericState` entirely on single-chain entries, so
            # `isComplex: false` is the only thing that says "monomer" there.
            state = "monomer"
        if not state:
            return "complex of unknown size"
        if state == "monomer":
            return "monomer"
        if state not in _OLIGOMERIC_NAMES.values():
            return f"{self.n_chains}-chain complex"
        kind = self.assembly_type
        if kind in ("Homo", "Hetero"):
            return f"{kind.lower()}{state}"
        return state

    @property
    def composition_phrase(self) -> str:
        """`'one protein in 2 copies: P0A6Q3 x2'`, or `'2 proteins: A x1 + B x1'`."""
        if not self.observed_composition:
            return "chain identities unknown"
        listing = format_composition(self.observed_composition)
        if self.n_chains == 1:
            return listing
        n = self.n_distinct_proteins
        head = (f"one protein in {self.n_chains} copies" if n == 1
                else f"{n} distinct proteins")
        return f"{head}: {listing}"

    @property
    def headline(self) -> str:
        """
        The one-line answer to "what am I looking at?".

        Example
        -------
        >>> hetero = AssemblyDescription(
        ...     'AF-2', ('A', 'B'), {'A': 181, 'B': 101}, True,
        ...     (('P63166', 1), ('Q96AZ6', 1)))
        >>> hetero.headline
        'Heterodimer -- 2 chains (A, B), 2 distinct proteins: P63166 x1 + Q96AZ6 x1'
        """
        line = self._headline_plain
        if self.conflicts:
            line += "  [!] the AFDB metadata disagrees -- see below"
        return line

    @property
    def _headline_plain(self) -> str:
        """`headline` without the "see below" flag, for contexts with no below."""
        chains = ", ".join(self.chain_ids) or "none"
        noun = self.noun[:1].upper() + self.noun[1:]
        return (f"{noun} -- {self.n_chains} chain"
                f"{'' if self.n_chains == 1 else 's'} ({chains}), "
                f"{self.composition_phrase}")

    @property
    def declared_phrase(self) -> str:
        """What the metadata asserts, verbatim, or why there is nothing to quote."""
        if not self.metadata_present:
            return "not fetched (local file mode)"
        parts: List[str] = []
        if self.declared_type:
            parts.append(f"assemblyType={self.declared_type}")
        if self.declared_state:
            parts.append(f"oligomericState={self.declared_state}")
        if self.declared_description:
            parts.append(f"oligomericStateDescription={self.declared_description}")
        if self.declared_is_complex is not None:
            parts.append(f"isComplex={str(self.declared_is_complex).lower()}")
        if self.declared_composition:
            # Same members, same order on the page: when the declaration matches
            # the chains, show it in the chains' order so the reader is not left
            # comparing two differently-sorted lists of the same thing.
            shown = self.declared_composition
            if sorted(shown) == sorted(self.observed_composition):
                shown = self.observed_composition
            parts.append(f"complexComposition {format_composition(shown)}")
        if not parts:
            return "the endpoint asserted none of assemblyType, oligomericState, " \
                   "complexComposition or isComplex"
        return ", ".join(parts)

    @property
    def agreement_phrase(self) -> str:
        """`'agrees with the chains'`, or the count of disagreements."""
        if not self.metadata_present:
            return ""
        if self.conflicts:
            n = len(self.conflicts)
            return f"DISAGREES with the chains ({n} point{'' if n == 1 else 's'})"
        if self.declared_type or self.declared_state or self.declared_composition:
            return "agrees with the chains"
        return ""

    @property
    def is_supported(self) -> bool:
        """Whether `require_dimer` would pass."""
        try:
            self.require_dimer()
        except UnsupportedAssemblyError:
            return False
        return True

    def describe(self) -> str:
        """
        The multi-line assembly block for the metadata report.

        Example
        -------
        >>> lying = AFDBPrediction('AF-3', (
        ...     {'chainId': 'A', 'uniprotAccession': 'P0A6Q3',
        ...      'assemblyType': 'Hetero', 'oligomericState': 'dimer'},
        ...     {'chainId': 'B', 'uniprotAccession': 'P0A6Q3',
        ...      'assemblyType': 'Hetero', 'oligomericState': 'dimer'}))
        >>> block = describe_assembly(lying, {'A': 172, 'B': 172}).describe()
        >>> print(block.splitlines()[0])
        Assembly      : Homodimer -- 2 chains (A, B)
        >>> flat = ' '.join(block.split())      # the block wraps; the text is one string
        >>> 'DISAGREES with the chains (1 point)' in flat
        True
        >>> "AFDB says assemblyType 'Hetero'" in flat
        True
        >>> 'both chains resolve to the same protein, P0A6Q3 x2' in flat
        True
        """
        chains = ", ".join(self.chain_ids) or "none"
        noun = self.noun[:1].upper() + self.noun[1:]
        # No "see below" flag here: unlike `headline`, this block is followed by
        # the conflicts themselves, so pointing at them would be noise.
        lines = _wrap_row(
            "Assembly",
            f"{noun} -- {self.n_chains} chain"
            f"{'' if self.n_chains == 1 else 's'} ({chains})")
        lines += _wrap_row("Composition", self.composition_phrase)
        agreement = self.agreement_phrase
        lines += _wrap_row(
            "AFDB declares",
            self.declared_phrase + (f" -- {agreement}" if agreement else ""))
        for conflict in self.conflicts:
            wrapped = textwrap.wrap(conflict, width=_ASSEMBLY_ROW_WIDTH - 6,
                                    break_long_words=False,
                                    break_on_hyphens=False)
            lines.append(f"  [!] {wrapped[0]}")
            lines += [f"      {piece}" for piece in wrapped[1:]]
        return "\n".join(lines)

    def require_dimer(self) -> "AssemblyDescription":
        """
        Return self, or refuse this input with a message a non-expert can act on.

        Called twice on the online path, deliberately, at the two moments the two
        kinds of evidence first exist:

        1. straight after the metadata fetch, on the declaration alone. This
           costs nothing and saves three downloads, but the declaration is not
           authoritative about chain count -- the endpoint can describe fewer
           chains than the structure has (R020) -- so it refuses only when the
           declaration *and* the endpoint's own per-chain entry list agree that
           this is not a two-chain complex. A declaration the entries
           contradict is reported as a conflict, not acted on.
        2. inside `verify_chain_identity`, on the observed chains. That one is
           authoritative and always fires, online or in local-file mode.

        Returns:
            `self`, so it can be chained.

        Raises:
            UnsupportedAssemblyError: If this is not a two-chain complex.

        Example
        -------
        >>> monomer = AFDBPrediction('AF-O15552-F1', (
        ...     {'chainId': 'A', 'uniprotAccession': 'O15552', 'isComplex': False},))
        >>> try:
        ...     describe_assembly(monomer).require_dimer()
        ... except UnsupportedAssemblyError as exc:
        ...     print(str(exc).splitlines()[0])
        AF-O15552-F1 is a monomer; this notebook analyses two-chain dimers only.

        A declaration the endpoint's own entries contradict does not refuse
        anything here -- two chain entries are two chains, whatever the label
        says, and the structural gate will settle it:

        >>> mislabelled = AFDBPrediction('AF-6', (
        ...     {'chainId': 'A', 'uniprotAccession': 'P1',
        ...      'oligomericState': 'trimer'},
        ...     {'chainId': 'B', 'uniprotAccession': 'P2',
        ...      'oligomericState': 'trimer'}))
        >>> describe_assembly(mislabelled).require_dimer().declared_state
        'trimer'

        A dimer passes and hands itself back:

        >>> ok = AFDBPrediction('AF-1', ({'chainId': 'A', 'uniprotAccession': 'P1'},
        ...                              {'chainId': 'B', 'uniprotAccession': 'P1'}))
        >>> describe_assembly(ok, {'A': 5, 'B': 5}).require_dimer().noun
        'homodimer'
        """
        if self.chains_observed:
            if self.n_chains == SUPPORTED_N_CHAINS:
                return self
            raise UnsupportedAssemblyError(self._refusal(observed=True))
        # Metadata-only. Refuse solely when the declaration and the endpoint's
        # own chain list *both* say this is not a dimer. A declaration that the
        # entries themselves contradict -- `oligomericState: 'trimer'` on a
        # record with two chain entries -- is exactly the disagreement this
        # class exists to report, and refusing on it would be preferring the
        # metadata over the model on evidence not yet gathered. Such a case
        # falls through to the structural gate, which reads the real chains.
        declared_non_dimer = (
            self.declared_is_complex is False
            or (self.declared_state is not None
                and self.declared_state != SUPPORTED_OLIGOMERIC_STATE)
        )
        if declared_non_dimer and self.n_chains != SUPPORTED_N_CHAINS:
            raise UnsupportedAssemblyError(self._refusal(observed=False))
        return self

    def _refusal(self, observed: bool) -> str:
        """The refusal text: what was found, what is supported, why, what to do."""
        if observed:
            parts = [f"{cid} ({self.chain_lengths[cid]} residues)"
                     if cid in self.chain_lengths else cid
                     for cid in self.chain_ids]
            chains = _wrap_inline(parts, indent=16) or "(none)"
            found = (f"  found     : {self.n_chains} chain"
                     f"{'' if self.n_chains == 1 else 's'} in the structure and the "
                     f"PAE/pLDDT documents:\n"
                     f"                {chains}\n"
                     f"              {self._headline_plain}")
        else:
            described = ", ".join(self.chain_ids) or "(none)"
            found = (f"  found     : the AFDB metadata declares "
                     f"{self.declared_phrase},\n"
                     f"              and describes {self.n_chains} chain"
                     f"{'' if self.n_chains == 1 else 's'}: {described}")
        if self.metadata_present and observed:
            declared = textwrap.wrap(
                f"AFDB metadata: {self.declared_phrase}", width=60,
                break_long_words=False, break_on_hyphens=False)
            found += "".join(f"\n              {piece}" for piece in declared)

        if self.n_chains > SUPPORTED_N_CHAINS and observed:
            why = (
                "  why       : every score here is a property of one *ordered pair* of\n"
                f"              chains. {self.n_chains} chains give "
                f"{self.n_chains * (self.n_chains - 1)} ordered pairs and nothing in\n"
                "              the input says which one you meant, so a single ipSAE or\n"
                "              pDockQ for the whole assembly would be a claim about an\n"
                "              object this notebook never examined.\n"
                "  note      : AFDB itself contains no predictions with more than two\n"
                "              chains -- every isComplex entry is oligomericState\n"
                f"              '{SUPPORTED_OLIGOMERIC_STATE}' -- so this is almost "
                "certainly a local mmCIF."
            )
            todo = (
                "  what to do: extract the two chains you want to score into their own\n"
                "              mmCIF and run that, or run the notebook once per pair.\n"
                f"              A supported AFDB accession: {EXAMPLE_DIMER_ACCESSION}."
            )
        else:
            why = (
                "  why       : ipTM, ipSAE, pDockQ, pDockQ2 and LIS are all *inter-chain*\n"
                "              measurements -- they read the PAE block between two\n"
                "              different chains and the contacts across the interface.\n"
                "              With one chain there is no such block and no interface,\n"
                "              so there is no number to report, not even a bad one."
            )
            single = self.observed_composition or self.declared_composition
            entry = single[0][0] if single else ""
            todo = (
                "  what to do: choose an AFDB *complex* accession. Every complex in AFDB\n"
                "              is a dimer, and you can list some with:\n"
                f"                {AFDB_DIMER_SEARCH}\n"
                f"              A known-good homodimer: {EXAMPLE_DIMER_ACCESSION}."
            )
            if entry:
                todo += (
                    f"\n              For per-residue confidence of {entry} on its own, "
                    "the AFDB\n"
                    f"              entry page https://alphafold.ebi.ac.uk/entry/{entry} "
                    "already shows\n"
                    "              pLDDT and PAE; this notebook adds nothing for a "
                    "single chain."
                )
        supported = (
            f"  supported : exactly {SUPPORTED_N_CHAINS} chains "
            f"(oligomericState '{SUPPORTED_OLIGOMERIC_STATE}'), homodimer or\n"
            "              heterodimer alike."
        )
        headline = (
            f"{self.accession} is a {self.noun}; this notebook analyses "
            f"two-chain dimers only."
        )
        return "\n".join([headline, found, supported, why, todo])


def describe_assembly(
    prediction: Optional[AFDBPrediction] = None,
    chain_lengths: Optional[Mapping[str, int]] = None,
    labels: Optional[Mapping[str, "ChainLabel"]] = None,
    accession: Optional[str] = None,
) -> AssemblyDescription:
    """
    Reconcile what AFDB *says* the assembly is against what the chains *are*.

    Reads the four assembly fields the notebook previously ignored --
    `assemblyType`, `oligomericState`, `complexComposition` and `isComplex` --
    and compares each against evidence derived independently from the chains:
    their count, and the UniProt accession each resolves to. Every disagreement
    becomes a sentence in `conflicts`; none of them is resolved silently, because
    the two sources are independent and either can be the wrong one.

    Args:
        prediction:    Fetched metadata, or `None` in local-file mode.
        chain_lengths: `{chain_id: n_residues}` observed in the structure / PAE.
                       Omit it to describe the metadata alone, before download.
        labels:        `{chain_id: ChainLabel}`, used to identify a chain's
                       protein when the metadata does not describe that chain.
        accession:     Overrides the accession in messages.

    Returns:
        An `AssemblyDescription`.

    Example
    -------
    A heterodimer, agreeing:

    >>> hetero = AFDBPrediction('AF-4', (
    ...     {'chainId': 'A', 'uniprotAccession': 'Q96AZ6', 'assemblyType': 'Hetero',
    ...      'oligomericState': 'dimer', 'isComplex': True,
    ...      'complexComposition': [{'identifier': 'Q96AZ6', 'stoichiometry': 1},
    ...                             {'identifier': 'P63166', 'stoichiometry': 1}]},
    ...     {'chainId': 'B', 'uniprotAccession': 'P63166', 'assemblyType': 'Hetero',
    ...      'oligomericState': 'dimer', 'isComplex': True,
    ...      'complexComposition': [{'identifier': 'Q96AZ6', 'stoichiometry': 1},
    ...                             {'identifier': 'P63166', 'stoichiometry': 1}]}))
    >>> describe_assembly(hetero, {'A': 181, 'B': 101}).headline
    'Heterodimer -- 2 chains (A, B), 2 distinct proteins: Q96AZ6 x1 + P63166 x1'

    A stoichiometry that does not match the chains present:

    >>> wrong = AFDBPrediction('AF-5', (
    ...     {'chainId': 'A', 'uniprotAccession': 'P1', 'oligomericState': 'dimer',
    ...      'complexComposition': [{'identifier': 'P1', 'stoichiometry': 4}]},
    ...     {'chainId': 'B', 'uniprotAccession': 'P1', 'oligomericState': 'dimer',
    ...      'complexComposition': [{'identifier': 'P1', 'stoichiometry': 4}]}))
    >>> for conflict in describe_assembly(wrong, {'A': 5, 'B': 5}).conflicts:
    ...     print(conflict)
    AFDB says complexComposition 'P1 x4' but the chains present are 'P1 x2'. Only the chains in this model are analysed.

    Without a structure it describes the metadata alone, and says so:

    >>> describe_assembly(hetero).chains_observed
    False
    """
    accession = str(
        accession
        or (prediction.accession if prediction is not None else "")
        or "This model"
    )

    conflicts: List[str] = []

    declared_type, note = _declared_field(prediction, "assemblyType")
    if note:
        conflicts.append(note)
    declared_state, note = _declared_field(prediction, "oligomericState")
    if note:
        conflicts.append(note)
    declared_description, note = _declared_field(
        prediction, "oligomericStateDescription")
    if note:
        conflicts.append(note)
    raw_composition, note = _declared_field(prediction, "complexComposition")
    if note:
        conflicts.append(note)
    declared_is_complex, note = _declared_field(prediction, "isComplex")
    if note:
        conflicts.append(note)

    declared_composition = _normalise_composition(raw_composition)

    if chain_lengths is not None:
        lengths = {str(cid): int(n) for cid, n in chain_lengths.items()}
        chain_ids = tuple(sorted(lengths))
        chains_observed = True
    else:
        lengths = {}
        chain_ids = prediction.chain_ids if prediction is not None else ()
        chains_observed = False

    # Observed composition: what protein is each chain, from whichever source
    # knows. Two rules keep this honest.
    #
    # 1. Every chain is keyed from the *same* namespace. Comparing one chain's
    #    UniProt accession against another's protein name would make two copies
    #    of one protein look like two different ones -- exactly the false
    #    "Hetero" this function exists to catch.
    # 2. Incomplete knowledge yields `()` rather than a guess, so "unknown"
    #    never masquerades as "different".
    def _key(chain_id: str, namespace: str) -> str:
        if namespace == "uniprot":
            if prediction is not None and prediction.describes_chain(chain_id):
                value = prediction.chain_field(
                    chain_id, "uniprotAccession", default="")
                if value:
                    return str(value)
            if labels is not None and chain_id in labels and labels[chain_id].uniprot:
                return str(labels[chain_id].uniprot)
            return ""
        if labels is not None and chain_id in labels:
            label = labels[chain_id]
            return str(label.gene or label.protein_name or label.entry_name or "")
        return ""

    # Ordered by first appearance in chain order, not alphabetically: chain
    # order is already deterministic (chain ids are sorted), and a heterodimer
    # reads far better when the protein named first is chain A's. Comparisons
    # against the declared composition sort both sides, so the display order
    # cannot create a spurious disagreement.
    observed_composition: Tuple[Tuple[str, int], ...] = ()
    for namespace in ("uniprot", "name"):
        keys = [_key(chain_id, namespace) for chain_id in chain_ids]
        if keys and all(keys):
            counts: Dict[str, int] = {}
            for key in keys:
                counts[key] = counts.get(key, 0) + 1
            observed_composition = tuple(counts.items())
            break

    description = AssemblyDescription(
        accession=accession,
        chain_ids=tuple(chain_ids),
        chain_lengths=lengths,
        chains_observed=chains_observed,
        observed_composition=observed_composition,
        declared_type=str(declared_type) if declared_type is not None else None,
        declared_state=str(declared_state) if declared_state is not None else None,
        declared_description=(str(declared_description)
                              if declared_description is not None else None),
        declared_composition=declared_composition,
        declared_is_complex=(bool(declared_is_complex)
                             if declared_is_complex is not None else None),
        conflicts=(),
        metadata_present=prediction is not None,
    )

    # --- reconciliation ----------------------------------------------------
    observed_type = description.observed_type
    if declared_type and observed_type and str(declared_type) != observed_type:
        if observed_type == "Homo":
            conflicts.append(
                f"AFDB says assemblyType {str(declared_type)!r} (two different "
                f"proteins) but both chains resolve to the same protein, "
                f"{format_composition(observed_composition)}. The chains are what "
                f"every score is computed from, so the labels above follow them; "
                f"the metadata may be describing a different model, or may be wrong."
            )
        else:
            conflicts.append(
                f"AFDB says assemblyType {str(declared_type)!r} (one protein in "
                f"several copies) but the chains resolve to different proteins, "
                f"{format_composition(observed_composition)}. The chains are what "
                f"every score is computed from, so the labels above follow them; "
                f"the metadata may be describing a different model, or may be wrong."
            )

    observed_state = description.observed_state
    if declared_state and observed_state and str(declared_state) != observed_state:
        conflicts.append(
            f"AFDB says oligomericState {str(declared_state)!r} but the structure "
            f"and documents contain {description.n_chains} chain"
            f"{'' if description.n_chains == 1 else 's'} "
            f"({observed_state}). The chain count is taken from the files actually "
            f"parsed."
        )

    if declared_composition and observed_composition and \
            sorted(declared_composition) != sorted(observed_composition):
        conflicts.append(
            f"AFDB says complexComposition "
            f"{format_composition(declared_composition)!r} but the chains present "
            f"are {format_composition(observed_composition)!r}. Only the chains in "
            f"this model are analysed."
        )

    if declared_is_complex is not None and chains_observed:
        if not declared_is_complex and description.n_chains > 1:
            conflicts.append(
                f"AFDB says isComplex=false but the structure has "
                f"{description.n_chains} chains."
            )
        if declared_is_complex and description.n_chains < 2:
            conflicts.append(
                f"AFDB says isComplex=true but the structure has only "
                f"{description.n_chains} chain."
            )

    return replace(description, conflicts=tuple(conflicts))


# --- metadata reporting (R020) ---------------------------------------------
# The report is a pure function returning a string rather than a pile of
# `print` calls, so that "the same input renders the same text" is something a
# caller can assert instead of something a reader has to eyeball.

_REPORT_WIDTH: int = 59
"""Rule width for `format_metadata_report`; sized to the longest protein name
seen on the fixtures."""

_REPORT_LABEL_WIDTH: int = 14
"""Column width of the label before the `:` in the report."""


def _chain_identity(prediction: AFDBPrediction, chain_id: str) -> Tuple[Tuple[str, str], ...]:
    """
    One chain's display identity, as ordered `(label, value)` rows.

    Chains whose rows are equal are the same protein, which is what lets the
    report collapse a homodimer into one block without a homo/hetero switch (D3).
    """
    accession = prediction.chain_field(chain_id, "uniprotAccession")
    entry_name = prediction.chain_field(chain_id, "uniprotId")
    if accession and entry_name:
        uniprot = f"{accession} ({entry_name})"
    else:
        uniprot = str(accession or entry_name or "N/A")
    return (
        ("UniProt", uniprot),
        ("Protein", str(prediction.chain_field(
            chain_id, "proteinFullName", "uniprotDescription", default="N/A"))),
        ("Gene", str(prediction.chain_field(
            chain_id, "geneNames", "gene", default="N/A"))),
        ("Organism", str(prediction.chain_field(
            chain_id, "organismScientificName", "organism", default="N/A"))),
    )


def format_metadata_report(
    prediction: Optional[AFDBPrediction],
    chain_lengths: Mapping[str, int],
    accession: Optional[str] = None,
    width: int = _REPORT_WIDTH,
    assembly: Optional[AssemblyDescription] = None,
    labels: Optional[Mapping[str, "ChainLabel"]] = None,
) -> str:
    """
    Render the complex's metadata, with every identity field attributed to a chain.

    The prediction endpoint describes one chain per entry, so UniProt accession,
    protein name, gene, organism and monomer length are *chain* facts, not complex
    facts, and this report never presents one chain's as the complex's (R020).
    Chains that share an identity are collapsed into a single block, so a
    homodimer reads as one protein in two chains while a heterodimer shows both --
    the same code path either way, with no user-facing switch (D3). Nothing is
    keyed on chain position or on the endpoint's entry order (D4).

    Args:
        prediction:    The fetched metadata, or `None` in local-file mode, where
                       the complex-level geometry is still reported and the
                       missing identity is stated rather than shown as `N/A`.
        chain_lengths: `{chain_id: n_residues}` from the structure or PAE, in the
                       order the chains should be reported.
        accession:     Overrides the accession line; defaults to the prediction's.
        width:         Rule width.
        assembly:      A pre-computed `AssemblyDescription`, so that the report
                       and the check that let the notebook get this far cannot
                       drift apart. Computed here when omitted.
        labels:        `{chain_id: ChainLabel}`, used only to identify a chain's
                       protein when the metadata does not describe it.

    Returns:
        The report as a newline-joined string, with no trailing newline.

    Example
    -------
    >>> hetero = AFDBPrediction('AF-1', (
    ...     {'chainId': 'B', 'uniprotAccession': 'P63166', 'gene': 'Sumo1'},
    ...     {'chainId': 'A', 'uniprotAccession': 'Q96AZ6', 'gene': 'ISG20'}))
    >>> report = format_metadata_report(hetero, {'A': 181, 'B': 101})
    >>> [line for line in report.splitlines() if line.startswith('Chain')]
    ['Chain A       : 181 residues', 'Chain B       : 101 residues']

    Entry order cannot reach the output, because `AFDBPrediction` sorts entries
    at construction:

    >>> jumbled = AFDBPrediction('AF-1', tuple(reversed(hetero.entries)))
    >>> format_metadata_report(jumbled, {'A': 181, 'B': 101}) == report
    True

    Identical chains collapse into one block:

    >>> homo = AFDBPrediction('AF-2', ({'chainId': 'A', 'uniprotAccession': 'P0A6Q3'},
    ...                                {'chainId': 'B', 'uniprotAccession': 'P0A6Q3'}))
    >>> [line for line in format_metadata_report(homo, {'A': 172, 'B': 172}).splitlines()
    ...  if line.startswith('Chain')]
    ['Chains A, B   : 172 residues each']

    A structure chain the endpoint never described is named as such, not blanked:

    >>> partial = AFDBPrediction('AF-3', ({'chainId': 'A', 'uniprotAccession': 'P1'},))
    >>> [line for line in format_metadata_report(partial, {'A': 10, 'B': 10}).splitlines()
    ...  if 'no AFDB' in line]
    ['  no AFDB metadata entry for this chain (endpoint described: A)']

    The assembly is stated in the reader's words, from both kinds of evidence
    (R023). What AFDB declares is printed beside what the chains show, so the two
    can be compared instead of one being quietly preferred:

    >>> for line in format_metadata_report(homo, {'A': 172, 'B': 172}).splitlines()[4:7]:
    ...     print(line)
    Assembly      : Homodimer -- 2 chains (A, B)
    Composition   : one protein in 2 copies: P0A6Q3 x2
    AFDB declares : the endpoint asserted none of assemblyType, oligomericState,

    With a chain the endpoint never described, the homo/hetero call is downgraded
    rather than asserted from metadata the chains could not confirm:

    >>> [line for line in format_metadata_report(partial, {'A': 10, 'B': 10}).splitlines()
    ...  if line.startswith(('Assembly', 'Composition'))]
    ['Assembly      : Dimer -- 2 chains (A, B)', 'Composition   : chain identities unknown']
    """
    lengths = {str(chain_id): int(n) for chain_id, n in chain_lengths.items()}
    chain_ids = list(lengths)

    def row(label: str, value: str, indent: int = 0) -> str:
        pad = " " * indent
        return f"{pad}{label:<{max(_REPORT_LABEL_WIDTH - indent, 1)}}: {value}"

    # Group chains by identity. Insertion order follows `chain_lengths`, so the
    # report order is the caller's chain order and never the endpoint's.
    groups: List[Tuple[Tuple[Tuple[str, str], ...], List[str]]] = []
    undescribed: List[str] = []
    if prediction is not None:
        for chain_id in chain_ids:
            if not prediction.describes_chain(chain_id):
                undescribed.append(chain_id)
                continue
            identity = _chain_identity(prediction, chain_id)
            for key, members in groups:
                if key == identity:
                    members.append(chain_id)
                    break
            else:
                groups.append((identity, [chain_id]))

    # The assembly is stated from *both* kinds of evidence -- what AFDB declares
    # and what the chains are -- and any disagreement between them is printed
    # rather than resolved out of sight (R023).
    if assembly is None:
        assembly = describe_assembly(prediction, lengths, labels=labels,
                                     accession=accession)
    if prediction is None:
        version = "not fetched (local file mode)"
    else:
        version = str(prediction.shared_field("modelVersion", "latestVersion",
                                              default="N/A"))

    lines: List[str] = ["=" * width, "COMPLEX METADATA REPORT", "=" * width]
    lines.append(row("Accession", str(
        accession or (prediction.accession if prediction is not None else "N/A"))))
    lines.extend(assembly.describe().splitlines())
    lines.append(row("Model version", version))
    lines.append(row("Total length", f"{sum(lengths.values())} residues ("
                     + ", ".join(f"{c}: {lengths[c]}" for c in chain_ids) + ")"))
    lines.append("-" * width)

    if prediction is None:
        lines.append(row("Chains", ", ".join(chain_ids) or "(none)"))
        lines.append("  no AFDB metadata was fetched, so the protein identity of "
                     "each chain")
        lines.append("  is unknown here rather than defaulted")
    else:
        for identity, members in groups:
            head = (f"Chain {members[0]}" if len(members) == 1
                    else "Chains " + ", ".join(members))
            sizes = {lengths[c] for c in members}
            if len(members) == 1:
                size = f"{lengths[members[0]]} residues"
            elif len(sizes) == 1:
                size = f"{sizes.pop()} residues each"
            else:
                size = ", ".join(f"{c}: {lengths[c]}" for c in members) + " residues"
            lines.append(row(head, size))
            for label, value in identity:
                lines.append(row(label, value, indent=2))
        for chain_id in undescribed:
            lines.append(row(f"Chain {chain_id}", f"{lengths[chain_id]} residues"))
            lines.append("  no AFDB metadata entry for this chain (endpoint "
                         f"described: {', '.join(prediction.chain_ids) or '(none)'})")
        extra = [c for c in prediction.chain_ids if c not in lengths]
        if extra:
            lines.append(row("Note", "AFDB metadata also describes chains absent "
                                     f"from the structure: {', '.join(extra)}"))

    lines.append("=" * width)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Structure parsing
# ---------------------------------------------------------------------------
# Column-naming convention, chosen once here and used everywhere (R015 folds the
# three divergent copies of this parser into this one): keys in `col_idx` have
# the `_atom_site.` prefix **stripped** and are **lower-cased**, so the mmCIF's
# `_atom_site.B_iso_or_equiv` is looked up as `'b_iso_or_equiv'`. Lower-casing
# is what makes the choice safe rather than merely a coin flip -- mmCIF tags are
# case-insensitive by specification, and the AFDB files mix cases within one
# loop (`group_PDB`, `Cartn_x`, `label_seq_id`).

_ATOM_SITE_PREFIX = "_atom_site."

_REQUIRED_ATOM_SITE_COLUMNS: Tuple[str, ...] = (
    "group_pdb",
    "label_atom_id",
    "label_asym_id",
    "label_seq_id",
    "label_comp_id",
    "cartn_x",
    "cartn_y",
    "cartn_z",
    "b_iso_or_equiv",
)
"""The columns every downstream computation needs. `pdbx_pdb_model_num` is
optional and filtered on only when present."""

BACKBONE_ATOMS: Tuple[str, str] = ("CA", "CB")
"""The only atoms retained. CB defines a contact; CA is the substitute for
glycine, which has no CB."""


@dataclass(eq=False)
class ChainCoords:
    """
    One chain's per-residue coordinates and confidence, sorted by residue number.

    Every array is parallel and has length `n_residues`, so a boolean interface
    mask indexes all four consistently.

    Attributes:
        chain_id:  `label_asym_id`, e.g. `'A'`.
        coords:    `(N, 3)` float32. CB, or CA where the residue has no CB.
        res_ids:   `(N,)` int32, `label_seq_id` values (1-based).
        res_names: `(N,)` str, three-letter residue codes.
        plddt:     `(N,)` float32, from the `B_iso_or_equiv` column, which AFDB
                   repurposes to carry pLDDT.
    """

    chain_id: str
    coords: np.ndarray
    res_ids: np.ndarray
    res_names: np.ndarray
    plddt: np.ndarray

    @property
    def n_residues(self) -> int:
        """Number of residues carried, i.e. the chain length as parsed."""
        return int(self.coords.shape[0])

    def __len__(self) -> int:
        return self.n_residues


def parse_mmcif_atoms(cif_text: str) -> Tuple[List[List[str]], Dict[str, int]]:
    """
    Tokenise the `_atom_site` loop of an mmCIF file. No BioPython, no gemmi.

    Column indices are built from the loop header rather than hard-coded, so the
    parser is agnostic to column order and to AFDB mmCIF version differences.

    Args:
        cif_text: The full text of an mmCIF file.

    Returns:
        `(records, col_idx)`. `records` is one token list per atom line;
        `col_idx` maps a lower-cased, prefix-stripped column name to its position
        within a record.

    Raises:
        ValueError: If the file contains no `_atom_site` loop.

    Example
    -------
    >>> cif = '''loop_
    ... _atom_site.group_PDB
    ... _atom_site.label_atom_id
    ... _atom_site.label_asym_id
    ... _atom_site.label_seq_id
    ... _atom_site.label_comp_id
    ... _atom_site.Cartn_x
    ... _atom_site.Cartn_y
    ... _atom_site.Cartn_z
    ... _atom_site.B_iso_or_equiv
    ... ATOM CA A 1 GLY 0.0 0.0 0.0 90.0
    ... ATOM CA A 2 ALA 4.0 0.0 0.0 80.0
    ... ATOM CB A 2 ALA 5.0 0.0 0.0 80.0
    ... #
    ... '''
    >>> records, col_idx = parse_mmcif_atoms(cif)
    >>> len(records)
    3
    >>> sorted(col_idx)[:3]
    ['b_iso_or_equiv', 'cartn_x', 'cartn_y']
    >>> records[0][col_idx['label_comp_id']]
    'GLY'
    """
    lines = cif_text.splitlines()
    n_lines = len(lines)

    headers: List[str] = []
    row_start: int | None = None

    index = 0
    while index < n_lines:
        if lines[index].strip().lower() != "loop_":
            index += 1
            continue

        cursor = index + 1
        candidate: List[str] = []
        while cursor < n_lines and lines[cursor].strip().startswith("_"):
            candidate.append(lines[cursor].strip())
            cursor += 1
        if candidate and candidate[0].lower().startswith(_ATOM_SITE_PREFIX):
            headers = candidate
            row_start = cursor
            break
        index = cursor if cursor > index else index + 1

    if row_start is None:
        raise ValueError("No _atom_site loop found in the mmCIF text.")

    col_idx = {
        header.lower()[len(_ATOM_SITE_PREFIX):].strip(): position
        for position, header in enumerate(headers)
    }
    n_cols = len(headers)

    records: List[List[str]] = []
    cursor = row_start
    while cursor < n_lines:
        stripped = lines[cursor].strip()
        if not stripped:
            cursor += 1
            continue
        if stripped.startswith("#"):
            break
        if stripped.startswith("_") or stripped.lower() == "loop_":
            break
        tokens = stripped.split()
        if len(tokens) >= n_cols:
            records.append(tokens[:n_cols])
        cursor += 1

    return records, col_idx


def extract_chain_coords(
    records: List[List[str]],
    col_idx: Dict[str, int],
) -> Dict[str, ChainCoords]:
    """
    Build per-chain `ChainCoords` from tokenised `_atom_site` records.

    CB is the contact atom; glycine has none, so its CA is substituted. A CB
    record always wins over a CA record for the same residue regardless of the
    order they appear in the file.

    Filtering, matching `ipsae.py`'s own atom selection:

    - `group_PDB == 'ATOM'` (HETATM, including any ligand or water, is dropped);
    - `label_atom_id` in `('CA', 'CB')`;
    - `label_seq_id` of `'.'` or `'?'`, or any non-integer, is skipped;
    - `pdbx_PDB_model_num == '1'` when that column exists, so an NMR-style
      multi-model file contributes only its first model.

    Args:
        records: From `parse_mmcif_atoms`.
        col_idx: From `parse_mmcif_atoms` (lower-cased, prefix-stripped keys).

    Returns:
        `{chain_id: ChainCoords}`, each chain's residues sorted by `label_seq_id`.

    Raises:
        ValueError: If a required `_atom_site` column is missing.

    Example
    -------
    >>> cif = '''loop_
    ... _atom_site.group_PDB
    ... _atom_site.label_atom_id
    ... _atom_site.label_asym_id
    ... _atom_site.label_seq_id
    ... _atom_site.label_comp_id
    ... _atom_site.Cartn_x
    ... _atom_site.Cartn_y
    ... _atom_site.Cartn_z
    ... _atom_site.B_iso_or_equiv
    ... ATOM CA A 1 GLY 0.0 0.0 0.0 90.0
    ... ATOM CA A 2 ALA 4.0 0.0 0.0 80.0
    ... ATOM CB A 2 ALA 5.0 0.0 0.0 80.0
    ... HETATM CA A 3 HOH 9.0 0.0 0.0 50.0
    ... '''
    >>> chains = extract_chain_coords(*parse_mmcif_atoms(cif))
    >>> sorted(chains)
    ['A']
    >>> chains['A'].res_names.tolist()
    ['GLY', 'ALA']
    >>> chains['A'].coords[:, 0].tolist()
    [0.0, 5.0]
    """
    missing = [name for name in _REQUIRED_ATOM_SITE_COLUMNS if name not in col_idx]
    if missing:
        raise ValueError(
            "mmCIF _atom_site loop is missing required column(s): "
            + ", ".join(missing)
        )

    ix_group = col_idx["group_pdb"]
    ix_atom = col_idx["label_atom_id"]
    ix_chain = col_idx["label_asym_id"]
    ix_seq = col_idx["label_seq_id"]
    ix_comp = col_idx["label_comp_id"]
    ix_x = col_idx["cartn_x"]
    ix_y = col_idx["cartn_y"]
    ix_z = col_idx["cartn_z"]
    ix_plddt = col_idx["b_iso_or_equiv"]
    ix_model = col_idx.get("pdbx_pdb_model_num")

    # chain_id -> res_id -> (atom_id, res_name, x, y, z, plddt)
    per_chain: Dict[str, Dict[int, Tuple[str, str, float, float, float, float]]] = {}

    for tokens in records:
        if tokens[ix_group] != "ATOM":
            continue
        if ix_model is not None and tokens[ix_model] != "1":
            continue
        atom_id = tokens[ix_atom]
        if atom_id not in BACKBONE_ATOMS:
            continue
        seq_token = tokens[ix_seq]
        if seq_token in (".", "?") or not seq_token.lstrip("-").isdigit():
            continue
        try:
            x, y, z = float(tokens[ix_x]), float(tokens[ix_y]), float(tokens[ix_z])
            plddt = float(tokens[ix_plddt])
        except ValueError:
            continue

        residues = per_chain.setdefault(tokens[ix_chain], {})
        res_id = int(seq_token)
        existing = residues.get(res_id)
        # CB wins over CA for the same residue, whichever order they appear in.
        if existing is None or (atom_id == "CB" and existing[0] == "CA"):
            residues[res_id] = (atom_id, tokens[ix_comp], x, y, z, plddt)

    chains: Dict[str, ChainCoords] = {}
    for chain_id, residues in per_chain.items():
        order = sorted(residues)
        chains[chain_id] = ChainCoords(
            chain_id=chain_id,
            coords=np.array([residues[r][2:5] for r in order], dtype=np.float32).reshape(-1, 3),
            res_ids=np.array(order, dtype=np.int32),
            res_names=np.array([residues[r][1] for r in order], dtype="<U3"),
            plddt=np.array([residues[r][5] for r in order], dtype=np.float32),
        )
    return chains


def parse_structure(cif_text: str) -> Dict[str, ChainCoords]:
    """
    `parse_mmcif_atoms` then `extract_chain_coords`, for the common case.

    Args:
        cif_text: The full text of an mmCIF file.

    Returns:
        `{chain_id: ChainCoords}`.
    """
    return extract_chain_coords(*parse_mmcif_atoms(cif_text))


# ---------------------------------------------------------------------------
# PAE / pLDDT parsing
# ---------------------------------------------------------------------------
# Both documents carry a `chains` array giving each chain's span in the
# concatenated complex; both are parsed into the same `ChainSpan` records so the
# two agree by construction. D4's ordered chain pair is made explicit here:
# `PAEMatrix.ordered_pair(x, y)` is the only supported way to get the four
# quadrants, and it always names which chain is rows and which is columns,
# because PAE is asymmetric and `block_xy != block_yx.T`.


@dataclass(frozen=True)
class ChainSpan:
    """
    One chain's residue span within a concatenated complex document.

    Attributes:
        chain_id: `label_asym_id`, e.g. `'A'`.
        start:    `sequenceStart`, 1-based inclusive.
        end:      `sequenceEnd`, inclusive.
        name:     The document's `name` field where present, `''` otherwise.
                  Used for real protein labels instead of "Chain A" (R021).

    Example
    -------
    >>> ChainSpan('A', 1, 172).length
    172
    """

    chain_id: str
    start: int
    end: int
    name: str = ""

    @property
    def length(self) -> int:
        """Number of residues spanned, inclusive of both ends."""
        return self.end - self.start + 1


def _chain_spans(
    entries: Any,
    fallback_lengths: Optional[Mapping[str, int]] = None,
) -> Tuple[ChainSpan, ...]:
    """
    Build sorted `ChainSpan`s from a document's `chains` array.

    Sorting by `label_asym_id` is what makes the slicing deterministic: the
    concatenated matrix is laid out in chain-label order, and the document's own
    array order is not guaranteed.

    Args:
        entries:          The document's `chains` array; may be empty or absent.
        fallback_lengths: `{chain_id: n_residues}`, normally structure-derived,
                          used when the document omits `chains` entirely.

    Returns:
        Spans sorted by chain id.

    Raises:
        ValueError: If neither a `chains` array nor a fallback is available.
    """
    if entries:
        spans = [
            ChainSpan(
                chain_id=str(entry["label_asym_id"]),
                start=int(entry["sequenceStart"]),
                end=int(entry["sequenceEnd"]),
                name=str(entry.get("name", "") or ""),
            )
            for entry in entries
        ]
        return tuple(sorted(spans, key=lambda span: span.chain_id))
    if fallback_lengths:
        return tuple(
            ChainSpan(chain_id=chain_id, start=1, end=int(fallback_lengths[chain_id]))
            for chain_id in sorted(fallback_lengths)
        )
    raise ValueError(
        "Document has no 'chains' array and no fallback chain lengths were given."
    )


@dataclass(frozen=True, eq=False)
class ChainPairPAE:
    """
    The four PAE quadrants of one **ordered** chain pair, `(x, y)`.

    PAE is asymmetric: `PAE[i][j]` is the error on residue `j` when the structure
    is aligned on residue `i`, so `block_xy` and `block_yx` are genuinely
    different measurements rather than transposes of one another. Every score
    that reports a `max` over both directions consumes exactly these two blocks.

    Attributes:
        chain_x:  Chain id used for the rows of `block_xy`.
        chain_y:  Chain id used for the columns of `block_xy`.
        block_xy: `(nx, ny)` inter-chain PAE, aligned on `chain_x`.
        block_yx: `(ny, nx)` inter-chain PAE, aligned on `chain_y`.
        block_xx: `(nx, nx)` intra-chain PAE for `chain_x`.
        block_yy: `(ny, ny)` intra-chain PAE for `chain_y`.
    """

    chain_x: str
    chain_y: str
    block_xy: np.ndarray
    block_yx: np.ndarray
    block_xx: np.ndarray
    block_yy: np.ndarray

    @property
    def nx(self) -> int:
        """Residue count of `chain_x`."""
        return int(self.block_xy.shape[0])

    @property
    def ny(self) -> int:
        """Residue count of `chain_y`."""
        return int(self.block_xy.shape[1])


@dataclass(frozen=True, eq=False)
class PAEMatrix:
    """
    A parsed PAE document.

    Attributes:
        matrix:  `(L, L)` float32, `L` being the total residue count.
        max_pae: `max_predicted_aligned_error`, the natural `vmax` for a heatmap.
        spans:   Per-chain spans, sorted by chain id, in matrix layout order.
    """

    matrix: np.ndarray
    max_pae: float
    spans: Tuple[ChainSpan, ...]

    @property
    def chain_ids(self) -> Tuple[str, ...]:
        """Chain ids in matrix layout order (sorted)."""
        return tuple(span.chain_id for span in self.spans)

    @property
    def chain_lengths(self) -> Dict[str, int]:
        """`{chain_id: n_residues}` as derived from the PAE document."""
        return {span.chain_id: span.length for span in self.spans}

    def chain_length(self, chain_id: str) -> int:
        """Residue count of one chain. Raises `KeyError` if it is unknown."""
        return self.chain_lengths[chain_id]

    def chain_slice(self, chain_id: str) -> slice:
        """
        The half-open row/column range one chain occupies in `matrix`.

        Args:
            chain_id: Chain id.

        Returns:
            A `slice` usable on either axis.

        Raises:
            KeyError: If the chain is not in the document.
        """
        offset = 0
        for span in self.spans:
            if span.chain_id == chain_id:
                return slice(offset, offset + span.length)
            offset += span.length
        raise KeyError(
            f"Chain {chain_id!r} is not in the PAE document; "
            f"available: {', '.join(self.chain_ids)}."
        )

    def block(self, row_chain: str, col_chain: str) -> np.ndarray:
        """
        One sub-block of the matrix, rows aligned on `row_chain`.

        Args:
            row_chain: Chain supplying the rows (the alignment frame).
            col_chain: Chain supplying the columns.

        Returns:
            A `(n_row, n_col)` view into `matrix`.
        """
        return self.matrix[self.chain_slice(row_chain), self.chain_slice(col_chain)]

    def ordered_pair(self, chain_x: str, chain_y: str) -> ChainPairPAE:
        """
        The four quadrants of an ordered chain pair (D4).

        Args:
            chain_x: Chain treated as "first"; rows of `block_xy`.
            chain_y: Chain treated as "second"; columns of `block_xy`.

        Returns:
            A `ChainPairPAE`.

        Raises:
            KeyError:   If either chain is absent.
            ValueError: If the two chain ids are the same, which would make the
                        "inter-chain" blocks intra-chain and silently wrong.
        """
        if chain_x == chain_y:
            raise ValueError(
                f"An ordered chain pair needs two distinct chains, got {chain_x!r} twice."
            )
        return ChainPairPAE(
            chain_x=chain_x,
            chain_y=chain_y,
            block_xy=self.block(chain_x, chain_y),
            block_yx=self.block(chain_y, chain_x),
            block_xx=self.block(chain_x, chain_x),
            block_yy=self.block(chain_y, chain_y),
        )


def parse_pae(
    pae_json: Any,
    fallback_lengths: Optional[Mapping[str, int]] = None,
) -> PAEMatrix:
    """
    Parse an AFDB PAE document.

    The document is a JSON array whose single element holds
    `predicted_aligned_error`, `max_predicted_aligned_error` and `chains`.

    Args:
        pae_json:         The parsed document (array, or the bare object).
        fallback_lengths: `{chain_id: n_residues}` used only when the document
                          omits `chains`; normally the structure-derived lengths.

    Returns:
        A `PAEMatrix`.

    Raises:
        ValueError: If the document is empty, if the matrix is not square, or if
            the chain spans do not account for exactly the matrix's extent --
            any of which would make every quadrant slice silently misaligned.

    Example
    -------
    >>> doc = [{'predicted_aligned_error': [[0.0, 5.0, 9.0],
    ...                                     [5.0, 0.0, 8.0],
    ...                                     [9.0, 8.0, 0.0]],
    ...         'max_predicted_aligned_error': 9.0,
    ...         'chains': [{'label_asym_id': 'B', 'sequenceStart': 1, 'sequenceEnd': 1},
    ...                    {'label_asym_id': 'A', 'sequenceStart': 1, 'sequenceEnd': 2}]}]
    >>> pae = parse_pae(doc)
    >>> pae.chain_ids
    ('A', 'B')
    >>> pae.chain_length('A'), pae.chain_length('B')
    (2, 1)
    >>> pair = pae.ordered_pair('A', 'B')
    >>> pair.block_xy.shape, pair.block_yx.shape
    ((2, 1), (1, 2))
    >>> pair.block_xy.ravel().tolist()
    [9.0, 8.0]
    """
    entry = pae_json[0] if isinstance(pae_json, list) else pae_json
    if not entry:
        raise ValueError("PAE document is empty.")

    matrix = np.array(entry["predicted_aligned_error"], dtype=np.float32)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"PAE matrix must be square, got shape {matrix.shape}.")

    spans = _chain_spans(entry.get("chains", []), fallback_lengths)
    total = sum(span.length for span in spans)
    if total != matrix.shape[0]:
        raise ValueError(
            f"PAE chain spans total {total} residues but the matrix is "
            f"{matrix.shape[0]}x{matrix.shape[0]}."
        )

    return PAEMatrix(
        matrix=matrix,
        max_pae=float(entry.get("max_predicted_aligned_error", 31.75)),
        spans=spans,
    )


@dataclass(frozen=True, eq=False)
class PLDDTScores:
    """
    A parsed AFDB pLDDT document, sliced per chain on demand.

    Attributes:
        scores:          `(L,)` float32, in complex order.
        residue_numbers: `(L,)` int32, the document's `residueNumber` column.
        spans:           Per-chain spans, sorted by chain id, in document order.

    Note:
        `ChainCoords.plddt` carries the same numbers, read from the mmCIF
        B-factor column. This class exists because the JSON document is also
        available when only PAE-side work is being done, and because AFDB's
        `confidenceCategory` bands are attached to it rather than to the mmCIF.
    """

    scores: np.ndarray
    residue_numbers: np.ndarray
    spans: Tuple[ChainSpan, ...]

    @property
    def chain_ids(self) -> Tuple[str, ...]:
        """Chain ids in document layout order (sorted)."""
        return tuple(span.chain_id for span in self.spans)

    def chain_slice(self, chain_id: str) -> slice:
        """The half-open range one chain occupies in `scores`."""
        offset = 0
        for span in self.spans:
            if span.chain_id == chain_id:
                return slice(offset, offset + span.length)
            offset += span.length
        raise KeyError(
            f"Chain {chain_id!r} is not in the pLDDT document; "
            f"available: {', '.join(self.chain_ids)}."
        )

    def for_chain(self, chain_id: str) -> np.ndarray:
        """
        One chain's pLDDT values.

        Args:
            chain_id: Chain id.

        Returns:
            A `(n,)` float32 view into `scores`.
        """
        return self.scores[self.chain_slice(chain_id)]


def parse_plddt(
    plddt_json: Any,
    fallback_lengths: Optional[Mapping[str, int]] = None,
) -> PLDDTScores:
    """
    Parse an AFDB pLDDT document.

    Args:
        plddt_json:       The parsed document object.
        fallback_lengths: `{chain_id: n_residues}` used only when the document
                          omits `chains`.

    Returns:
        A `PLDDTScores`.

    Raises:
        ValueError: If the chain spans do not account for exactly as many
            residues as there are scores.

    Example
    -------
    >>> doc = {'residueNumber': [1, 2, 1],
    ...        'confidenceScore': [90.0, 80.0, 70.0],
    ...        'chains': [{'label_asym_id': 'A', 'sequenceStart': 1, 'sequenceEnd': 2},
    ...                   {'label_asym_id': 'B', 'sequenceStart': 1, 'sequenceEnd': 1}]}
    >>> plddt = parse_plddt(doc)
    >>> plddt.for_chain('A').tolist()
    [90.0, 80.0]
    >>> plddt.for_chain('B').tolist()
    [70.0]
    """
    scores = np.array(plddt_json["confidenceScore"], dtype=np.float32)
    residue_numbers = np.array(plddt_json.get("residueNumber", []), dtype=np.int32)
    spans = _chain_spans(plddt_json.get("chains", []), fallback_lengths)

    total = sum(span.length for span in spans)
    if total != scores.shape[0]:
        raise ValueError(
            f"pLDDT chain spans total {total} residues but the document carries "
            f"{scores.shape[0]} scores."
        )

    return PLDDTScores(scores=scores, residue_numbers=residue_numbers, spans=spans)


def _document_kind(document: "PAEMatrix | PLDDTScores") -> str:
    """`'PAE'` or `'pLDDT'`, so an error message can name the guilty document."""
    return "PAE" if isinstance(document, PAEMatrix) else "pLDDT"


_MAPPING_REFUSAL: str = (
    "Refusing to pair these chains positionally. Neither source carries an order "
    "that could be mapped on: the documents' 'chains' arrays are re-sorted by "
    "chain id here precisely because their arrival order is not authoritative, "
    "and the AFDB prediction endpoint returns its per-chain entries in a "
    "demonstrably non-deterministic order (R020). A positional pairing would "
    "therefore be a guess, and a wrong guess mis-slices every PAE quadrant into "
    "plausible-looking but wrong scores rather than an error. Chain length is no "
    "tie-breaker either: a homodimer's chains are the same length by definition."
)
"""Why a chain-id mismatch fails instead of falling back to positional mapping.

Shared by every mismatch message so the reasoning is stated once (R021)."""

_MISMATCH_REMEDY: str = (
    "Check that the mmCIF, the PAE document and the pLDDT document are all for "
    "the same accession and model version -- mixing files from two predictions, "
    "or re-uploading one of the three in local-file mode, is the usual cause. "
    "Fetch all three from one AFDBPrediction if in doubt."
)
"""What the user should actually do about a mismatch."""


def verify_chain_lengths(
    chains: Mapping[str, ChainCoords],
    document: PAEMatrix | PLDDTScores,
) -> Dict[str, int]:
    """
    Check structure-derived chain lengths against a document's own spans.

    A mismatch means every quadrant slice and every per-chain pLDDT slice is
    misaligned, which produces plausible-looking but wrong scores rather than an
    error, so it is checked rather than assumed (R021).

    Args:
        chains:   `{chain_id: ChainCoords}` from `extract_chain_coords`.
        document: A `PAEMatrix` or `PLDDTScores` to compare against.

    Returns:
        The agreed `{chain_id: n_residues}`.

    Raises:
        ValueError: If the chain sets or any length disagree. The message names
            the two chain sets, or the disagreeing lengths, says why a
            positional fallback is refused, and says what to check.

    Example
    -------
    >>> spans = (ChainSpan('A', 1, 2), ChainSpan('B', 1, 1))
    >>> doc = PLDDTScores(np.zeros(3, dtype=np.float32),
    ...                   np.zeros(3, dtype=np.int32), spans)
    >>> chains = {'A': ChainCoords('A', np.zeros((2, 3), dtype=np.float32),
    ...                            np.array([1, 2]), np.array(['ALA', 'ALA']),
    ...                            np.zeros(2, dtype=np.float32)),
    ...           'B': ChainCoords('B', np.zeros((1, 3), dtype=np.float32),
    ...                            np.array([1]), np.array(['GLY']),
    ...                            np.zeros(1, dtype=np.float32))}
    >>> verify_chain_lengths(chains, doc)
    {'A': 2, 'B': 1}

    A chain the document does not describe is named, not skipped:

    >>> del chains['B']
    >>> chains['C'] = ChainCoords('C', np.zeros((1, 3), dtype=np.float32),
    ...                           np.array([1]), np.array(['GLY']),
    ...                           np.zeros(1, dtype=np.float32))
    >>> verify_chain_lengths(chains, doc)      # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ValueError: Chain identity mismatch between the structure and the pLDDT document.
      structure chains : A, C
      pLDDT chains     : A, B
      only in structure: C
      only in pLDDT    : B
    Refusing to pair these chains positionally...
    """
    kind = _document_kind(document)
    from_structure = {chain_id: chain.n_residues for chain_id, chain in chains.items()}
    from_document = {span.chain_id: span.length for span in document.spans}

    if set(from_structure) != set(from_document):
        only_structure = sorted(set(from_structure) - set(from_document))
        only_document = sorted(set(from_document) - set(from_structure))
        rows = (
            ("structure chains", sorted(from_structure)),
            (f"{kind} chains", sorted(from_document)),
            ("only in structure", only_structure),
            (f"only in {kind}", only_document),
        )
        detail = "\n".join(
            f"  {label:<17}: {', '.join(ids) or '(none)'}" for label, ids in rows
        )
        raise ValueError(
            f"Chain identity mismatch between the structure and the {kind} document.\n"
            f"{detail}\n"
            f"{_MAPPING_REFUSAL}\n{_MISMATCH_REMEDY}"
        )
    disagreeing = {
        chain_id: (from_structure[chain_id], from_document[chain_id])
        for chain_id in from_structure
        if from_structure[chain_id] != from_document[chain_id]
    }
    if disagreeing:
        detail = "\n".join(
            f"  chain {chain_id}: structure {s} residues, {kind} {d} residues"
            for chain_id, (s, d) in sorted(disagreeing.items())
        )
        raise ValueError(
            f"Chain length mismatch between the structure and the {kind} document.\n"
            f"{detail}\n"
            f"Every quadrant slice and every per-chain pLDDT slice is cut at these "
            f"offsets, so continuing would score one chain against the wrong "
            f"residues and report a number rather than an error.\n"
            f"{_MISMATCH_REMEDY}"
        )
    return dict(sorted(from_structure.items()))


def verify_document_agreement(pae: PAEMatrix, plddt: PLDDTScores) -> Dict[str, int]:
    """
    Check the PAE and pLDDT documents describe the same chains at the same lengths.

    `verify_chain_lengths` compares each document against the *structure*, which
    already makes them agree transitively -- but only when both calls are made,
    and only when a structure is present. This is the direct check, so the third
    edge of the triangle is closed explicitly rather than by inference (R021).

    Args:
        pae:   The parsed PAE document.
        plddt: The parsed pLDDT document.

    Returns:
        The agreed `{chain_id: n_residues}`.

    Raises:
        ValueError: If the chain sets or any length disagree.

    Example
    -------
    >>> spans = (ChainSpan('A', 1, 2), ChainSpan('B', 1, 1))
    >>> pae = PAEMatrix(np.zeros((3, 3), dtype=np.float32), 30.0, spans)
    >>> plddt = PLDDTScores(np.zeros(3, dtype=np.float32),
    ...                     np.zeros(3, dtype=np.int32), spans)
    >>> verify_document_agreement(pae, plddt)
    {'A': 2, 'B': 1}
    """
    from_pae = {span.chain_id: span.length for span in pae.spans}
    from_plddt = {span.chain_id: span.length for span in plddt.spans}

    if set(from_pae) != set(from_plddt):
        raise ValueError(
            f"Chain identity mismatch between the PAE and pLDDT documents.\n"
            f"  PAE chains  : {', '.join(sorted(from_pae)) or '(none)'}\n"
            f"  pLDDT chains: {', '.join(sorted(from_plddt)) or '(none)'}\n"
            f"{_MAPPING_REFUSAL}\n{_MISMATCH_REMEDY}"
        )
    disagreeing = sorted(
        chain_id for chain_id in from_pae if from_pae[chain_id] != from_plddt[chain_id]
    )
    if disagreeing:
        detail = "\n".join(
            f"  chain {chain_id}: PAE {from_pae[chain_id]} residues, "
            f"pLDDT {from_plddt[chain_id]} residues"
            for chain_id in disagreeing
        )
        raise ValueError(
            f"Chain length mismatch between the PAE and pLDDT documents.\n"
            f"{detail}\n{_MISMATCH_REMEDY}"
        )
    return dict(sorted(from_pae.items()))


# --- chain labelling (R021) -------------------------------------------------
# "Chain A" identifies a column of a matrix; it does not identify a protein. The
# PAE and pLDDT documents both carry a `name` per chain, and the prediction
# endpoint carries the gene, the UniProt accession and the entry name, so every
# label in the notebook can say *which protein* as well as which chain.
#
# Three forms rather than one, because the contexts differ by an order of
# magnitude in the room they have: `token` for a bare identifier, `short` for an
# axis label or a tick, `full` for a caption, a heading or printed output.

_TRUNCATION_MARKERS: Tuple[str, ...] = ("-", ",", "(", "/", "+")
"""Trailing characters that mean the service cut a protein name short.

The homodimer fixture's `name` arrives as `'3-hydroxydecanoyl-'`, a prefix of
"3-hydroxydecanoyl-[acyl-carrier-protein] dehydratase". The name is still worth
showing -- it is the only human-readable identity in local-file mode -- so it is
marked with an ellipsis rather than discarded."""


def _tidy_name(name: Optional[str]) -> str:
    """
    Normalise a protein name, marking an obviously truncated one.

    Args:
        name: A `name` field from a PAE / pLDDT `chains` entry, or a
              `uniprotDescription`. May be `None` or empty.

    Returns:
        The stripped name, with `'…'` appended when it ends mid-word; `''` when
        there is nothing to show.

    Example
    -------
    >>> _tidy_name('Small ubiquitin-related modifier 1')
    'Small ubiquitin-related modifier 1'
    >>> _tidy_name('3-hydroxydecanoyl-')
    '3-hydroxydecanoyl-…'
    >>> _tidy_name(None), _tidy_name('   ')
    ('', '')
    """
    text = str(name or "").strip()
    if not text:
        return ""
    if text.endswith(_TRUNCATION_MARKERS):
        return text + "…"
    return text


@dataclass(frozen=True)
class ChainLabel:
    """
    How one chain should be named on a figure, in a caption, or in printed output.

    A homodimer gives both chains the same protein, so **every** form carries the
    chain id: dropping it would make the two panels of a per-chain figure
    indistinguishable. Every field is optional and every form degrades to
    `'Chain A'` when nothing is known, which is what local-file mode and a
    metadata-less document produce.

    Attributes:
        chain_id:     `label_asym_id`, e.g. `'A'`. The only required field.
        protein_name: Full protein name, e.g. the PAE document's `name`.
        gene:         Gene name, e.g. `'ISG20'`.
        uniprot:      UniProt accession, e.g. `'Q96AZ6'`.
        entry_name:   UniProt entry name, e.g. `'ISG20_HUMAN'`.

    Example
    -------
    >>> isg20 = ChainLabel('A', 'Interferon-stimulated gene 20 kDa protein',
    ...                    gene='ISG20', uniprot='Q96AZ6', entry_name='ISG20_HUMAN')
    >>> isg20.token
    'ISG20'
    >>> isg20.short
    'ISG20 (A)'
    >>> isg20.full
    'Chain A: Interferon-stimulated gene 20 kDa protein (Q96AZ6)'

    With nothing known it degrades to the chain id, never to an empty string:

    >>> ChainLabel('B').token, ChainLabel('B').short, ChainLabel('B').full
    ('', 'Chain B', 'Chain B')

    A homodimer's two chains stay distinguishable even though the protein is one:

    >>> [ChainLabel(c, gene='fabA').short for c in 'AB']
    ['fabA (A)', 'fabA (B)']
    """

    chain_id: str
    protein_name: str = ""
    gene: str = ""
    uniprot: str = ""
    entry_name: str = ""

    @property
    def token(self) -> str:
        """
        The shortest identifier for the *protein*, with no chain id: `'ISG20'`.

        Gene first, then accession, then entry name -- a gene symbol is the form
        a biologist reads fastest, and it is what fits a tick label. `''` when
        none is known; the protein name is deliberately not used here, because a
        41-character name is not a token.
        """
        return str(self.gene or self.uniprot or self.entry_name or "")

    @property
    def short(self) -> str:
        """
        Axis-label and tick form: `'ISG20 (A)'`, falling back to `'Chain A'`.

        Short enough for an axis label, a tick and a legend entry, and still
        unambiguous for a homodimer because the chain id is always present.
        """
        token = self.token
        return f"{token} ({self.chain_id})" if token else f"Chain {self.chain_id}"

    @property
    def full(self) -> str:
        """
        Caption and printed-output form:
        `'Chain A: Interferon-stimulated gene 20 kDa protein (Q96AZ6)'`.

        Degrades one part at a time: without an accession the parenthesis is
        dropped, without a name the accession takes its place, and with neither
        it is just `'Chain A'`.
        """
        head = f"Chain {self.chain_id}"
        name = self.protein_name
        if name and self.uniprot:
            return f"{head}: {name} ({self.uniprot})"
        if name:
            return f"{head}: {name}"
        if self.token:
            return f"{head}: {self.token}"
        return head

    def __str__(self) -> str:
        """`short`, so a `ChainLabel` can be dropped into an f-string."""
        return self.short


@dataclass(frozen=True)
class ChainIdentity:
    """
    The verified, cross-checked identity of every chain in one complex.

    Produced by `verify_chain_identity`, which is the single place the structure,
    the PAE document, the pLDDT document and (when online) the prediction
    metadata are reconciled. Holding the agreed lengths and the labels together
    is deliberate: a caller cannot get the labels without the check having passed.

    Attributes:
        lengths: The agreed `{chain_id: n_residues}`, sorted by chain id.
        labels:  `{chain_id: ChainLabel}`, same keys as `lengths`.
        notes:   Non-fatal observations -- a chain the metadata never described,
                 a protein name the two documents disagree about. Worth printing;
                 not worth refusing to run over, because none of them can
                 misalign a slice.
        assembly: What this complex is (R023), reconciled from the AFDB
                 declaration and the chains themselves. `None` only when a
                 caller constructs a `ChainIdentity` by hand;
                 `verify_chain_identity` always fills it, and always after
                 `AssemblyDescription.require_dimer` has passed.
    """

    lengths: Dict[str, int]
    labels: Dict[str, ChainLabel]
    notes: Tuple[str, ...] = ()
    assembly: Optional[AssemblyDescription] = None

    @property
    def chain_ids(self) -> Tuple[str, ...]:
        """Chain ids, sorted."""
        return tuple(self.lengths)

    def label(self, chain_id: str) -> ChainLabel:
        """
        One chain's `ChainLabel`.

        Args:
            chain_id: Chain id.

        Returns:
            Its label.

        Raises:
            KeyError: If the chain is not part of this complex.
        """
        try:
            return self.labels[chain_id]
        except KeyError:
            raise KeyError(
                f"Chain {chain_id!r} is not part of this complex; "
                f"it has chains {', '.join(self.chain_ids) or '(none)'}."
            ) from None

    def short(self, chain_id: str) -> str:
        """`ChainLabel.short` for one chain: `'ISG20 (A)'`."""
        return self.label(chain_id).short

    def full(self, chain_id: str) -> str:
        """`ChainLabel.full` for one chain."""
        return self.label(chain_id).full

    def legend(self) -> str:
        """
        One line per chain, mapping the id used on the axes to the protein.

        The companion to the compact forms: figures and tabular output stay
        narrow by using `short`, and this says once, in full, what each id means.

        Example
        -------
        >>> identity = ChainIdentity(
        ...     {'A': 181, 'B': 101},
        ...     {'A': ChainLabel('A', 'Interferon-stimulated gene 20 kDa protein',
        ...                      gene='ISG20', uniprot='Q96AZ6'),
        ...      'B': ChainLabel('B', 'Small ubiquitin-related modifier 1',
        ...                      gene='Sumo1', uniprot='P63166')})
        >>> print(identity.legend())
        Chain A: Interferon-stimulated gene 20 kDa protein (Q96AZ6), 181 residues, labelled "ISG20 (A)"
        Chain B: Small ubiquitin-related modifier 1 (P63166), 101 residues, labelled "Sumo1 (B)"
        """
        return "\n".join(
            f"{self.labels[chain_id].full}, {self.lengths[chain_id]} residues, "
            f'labelled "{self.labels[chain_id].short}"'
            for chain_id in self.chain_ids
        )


def verify_chain_identity(
    chains: Mapping[str, ChainCoords],
    pae: PAEMatrix,
    plddt: Optional[PLDDTScores] = None,
    prediction: Optional[AFDBPrediction] = None,
    require_dimer: bool = True,
) -> ChainIdentity:
    """
    Reconcile the three chain-describing sources, then build a label per chain.

    The structure, the PAE document and the pLDDT document each carry their own
    chain list, and every PAE quadrant and every per-chain pLDDT slice is cut on
    the assumption that all three agree. `CLAUDE.md`'s edge-case table warns that
    asymmetric chain labels can misalign that slicing; this makes the agreement a
    checked precondition instead (R021). All three edges of the triangle are
    tested -- structure/PAE, structure/pLDDT and PAE/pLDDT -- so the check does
    not depend on which arguments a caller happened to pass.

    **Chain ids are never mapped positionally.** A structure labelled `A`/`C`
    against a PAE document labelled `A`/`B` raises. See `_MAPPING_REFUSAL` for
    why: no source here carries a trustworthy order to map on, so a positional
    pairing would be a guess whose failure mode is a plausible wrong number.

    Args:
        chains:     `{chain_id: ChainCoords}` from `parse_structure`.
        pae:        The parsed PAE document.
        plddt:      The parsed pLDDT document, when there is one.
        prediction: The fetched AFDB metadata, when online. Supplies the gene,
                    the UniProt accession and the entry name; without it the
                    labels fall back to the documents' own `name` field, and
                    then to the bare chain id.
        require_dimer: Enforce the notebook's dimer scope (D4) here, rather
                    than leaving it to be discovered eight cells later as an
                    `IndexError` on `chain_ids[1]`. Pass `False` to describe the
                    assembly without restricting it.

    Returns:
        A `ChainIdentity` carrying the agreed lengths, one `ChainLabel` per chain,
        any non-fatal notes, and the reconciled `AssemblyDescription`.

    Raises:
        ValueError: If the chain sets or the chain lengths disagree between any
            two of the three sources.
        UnsupportedAssemblyError: If `require_dimer` and this is not a
            two-chain complex. This is the gate that cannot be skipped: every
            cell downstream consumes this function's return value, so there is
            no path to a score that does not pass through it.

    Example
    -------
    >>> spans = (ChainSpan('A', 1, 2, name='Alpha protein'),
    ...          ChainSpan('B', 1, 1, name='Beta protein'))
    >>> pae = PAEMatrix(np.zeros((3, 3), dtype=np.float32), 30.0, spans)
    >>> chains = {'A': ChainCoords('A', np.zeros((2, 3), dtype=np.float32),
    ...                            np.array([1, 2]), np.array(['ALA', 'ALA']),
    ...                            np.zeros(2, dtype=np.float32)),
    ...           'B': ChainCoords('B', np.zeros((1, 3), dtype=np.float32),
    ...                            np.array([1]), np.array(['GLY']),
    ...                            np.zeros(1, dtype=np.float32))}
    >>> identity = verify_chain_identity(chains, pae)
    >>> identity.lengths
    {'A': 2, 'B': 1}
    >>> identity.short('A'), identity.full('A')
    ('Chain A', 'Chain A: Alpha protein')

    With metadata, the gene becomes the short form:

    >>> pred = AFDBPrediction('AF-1', ({'chainId': 'A', 'gene': 'alp',
    ...                                 'uniprotAccession': 'P1'},
    ...                                {'chainId': 'B', 'gene': 'bet',
    ...                                 'uniprotAccession': 'P2'}))
    >>> identity = verify_chain_identity(chains, pae, prediction=pred)
    >>> identity.short('A'), identity.short('B')
    ('alp (A)', 'bet (B)')

    A metadata chain the structure does not have is a note, not a failure:

    >>> partial = AFDBPrediction('AF-2', ({'chainId': 'A', 'gene': 'alp'},))
    >>> verify_chain_identity(chains, pae, prediction=partial).notes
    ('Chain B: the AFDB metadata describes no such chain (it described: A); labelling it from the PAE/pLDDT documents alone.',)

    The assembly comes back with the identity, reconciled and already checked:

    >>> verify_chain_identity(chains, pae, prediction=pred).assembly.noun
    'heterodimer'

    A single-chain model is refused here, not eight cells later on an index that
    does not exist:

    >>> solo_pae = PAEMatrix(np.zeros((2, 2), dtype=np.float32), 30.0,
    ...                      (ChainSpan('A', 1, 2),))
    >>> try:
    ...     verify_chain_identity({'A': chains['A']}, solo_pae)
    ... except UnsupportedAssemblyError as exc:
    ...     print(str(exc).splitlines()[0])
    This model is a monomer; this notebook analyses two-chain dimers only.
    """
    lengths = verify_chain_lengths(chains, pae)
    if plddt is not None:
        verify_chain_lengths(chains, plddt)
        verify_document_agreement(pae, plddt)

    doc_names: Dict[str, Dict[str, str]] = {}
    for source, document in (("PAE", pae), ("pLDDT", plddt)):
        if document is None:
            continue
        for span in document.spans:
            doc_names.setdefault(span.chain_id, {})[source] = _tidy_name(span.name)

    notes: List[str] = []
    labels: Dict[str, ChainLabel] = {}
    for chain_id in lengths:
        per_source = doc_names.get(chain_id, {})
        distinct = {name for name in per_source.values() if name}
        if len(distinct) > 1:
            notes.append(
                f"Chain {chain_id}: the PAE and pLDDT documents give different "
                f"protein names ({'; '.join(sorted(distinct))}); using the PAE "
                f"document's. Labels only -- no slice depends on this."
            )
        name = per_source.get("PAE") or per_source.get("pLDDT") or ""

        gene = uniprot = entry_name = ""
        if prediction is not None:
            if prediction.describes_chain(chain_id):
                gene = str(prediction.chain_field(
                    chain_id, "geneNames", "gene", default="") or "")
                uniprot = str(prediction.chain_field(
                    chain_id, "uniprotAccession", default="") or "")
                entry_name = str(prediction.chain_field(
                    chain_id, "uniprotId", default="") or "")
                name = _tidy_name(prediction.chain_field(
                    chain_id, "proteinFullName", "uniprotDescription",
                    default="")) or name
            else:
                notes.append(
                    f"Chain {chain_id}: the AFDB metadata describes no such chain "
                    f"(it described: {', '.join(prediction.chain_ids) or '(none)'}); "
                    f"labelling it from the PAE/pLDDT documents alone."
                )

        labels[chain_id] = ChainLabel(
            chain_id=chain_id,
            protein_name=name,
            gene=gene,
            uniprot=uniprot,
            entry_name=entry_name,
        )

    if prediction is not None:
        extra = [c for c in prediction.chain_ids if c not in lengths]
        if extra:
            notes.append(
                f"The AFDB metadata describes chains absent from the structure: "
                f"{', '.join(extra)}. They are not analysed."
            )

    # What is this thing? Answered from the chains just verified *and* from the
    # AFDB declaration, with any disagreement recorded rather than resolved
    # (R023). The scope check goes here, after the length checks and before any
    # caller can reach a score, because this is the one call every downstream
    # cell depends on: a monomer or a three-chain model is refused with an
    # explanation instead of an `IndexError` on `chain_ids[1]`.
    assembly = describe_assembly(prediction, lengths, labels=labels)
    if require_dimer:
        assembly.require_dimer()

    return ChainIdentity(lengths=lengths, labels=labels, notes=tuple(notes),
                         assembly=assembly)


# ---------------------------------------------------------------------------
# Interface detection
# ---------------------------------------------------------------------------
# CB-CB (CA for glycine) contacts at `DIST_CUTOFF`, matching `ipsae_v4.py:385`
# and its inclusive `dist <= cutoff` test (`ipsae_v4.py:652`).
#
# The two counts are named apart deliberately (R015). `interface.py` calls
# `mask_x.sum() + mask_y.sum()` "n_contacts" while the notebook uses the same
# word for `contact_mask.sum()`; they are different numbers, and pDockQ needs the
# **pair** count (`ipsae_v4.py:653`, `npairs`). Here they are
# `n_interface_residues` and `n_contact_pairs`, so neither can be passed where
# the other is meant.
#
# Intentional divergence from AFDB production code, carried over from
# `interface.py`'s module docstring: AFDB's own production `interface.py`
# (Majewski, Apache 2.0) uses **CA-CA** distances via a PyTorch `radius_graph`
# on GPU. This module uses **CB-CB** (CA for glycine) because that is the
# contact definition `ipsae.py` scores against, and torch is a prohibited
# dependency here. Do not "fix" the atom-type selection without also updating
# `compute_pdockq` and `compute_pdockq2`, which consume these contacts.
#
# R015 reconciliation record. `src/insightfold/interface.py` -- the earlier,
# unused module -- was compared symbol by symbol against this section and the
# structure-parsing section above. Everything it does is covered here, and the
# two agree exactly on both fixtures (`AF-0000000065889468`,
# `AF-0000000211034637`): identical chain ids, residue counts, coordinates,
# pLDDT, distance matrices, interface masks and contact pairs. The only thing
# ported was the CA-CA/CB-CB warning immediately above; the remaining
# differences are all improvements already present here (a `ValueError` naming
# the missing `_atom_site` columns rather than a generic message; a `ValueError`
# when there is no `_atom_site` loop at all rather than a silent empty result;
# blank lines inside the atom loop skipped rather than truncating it). Per
# decision D10 `interface.py` stays on disk until R095, so an early deletion
# cannot masquerade as a later bug.


@dataclass(frozen=True, eq=False)
class InterfaceContacts:
    """
    Inter-chain contacts of one ordered chain pair.

    Attributes:
        chain_x:      Chain id along the rows.
        chain_y:      Chain id along the columns.
        dist_matrix:  `(nx, ny)` float32, all pairwise CB/CA distances in Angstrom.
        contact_mask: `(nx, ny)` bool, `dist_matrix <= dist_cutoff`.
        mask_x:       `(nx,)` bool, residues of `chain_x` with any contact.
        mask_y:       `(ny,)` bool, residues of `chain_y` with any contact.
        dist_cutoff:  The cutoff in force, in Angstrom.
    """

    chain_x: str
    chain_y: str
    dist_matrix: np.ndarray
    contact_mask: np.ndarray
    mask_x: np.ndarray
    mask_y: np.ndarray
    dist_cutoff: float

    @property
    def n_contact_pairs(self) -> int:
        """Residue **pairs** within the cutoff. This is pDockQ's `npairs`."""
        return int(self.contact_mask.sum())

    @property
    def n_interface_residues(self) -> int:
        """Interface **residues**, both chains summed. Not `n_contact_pairs`."""
        return int(self.mask_x.sum()) + int(self.mask_y.sum())

    @property
    def n_interface_residues_x(self) -> int:
        """Interface residues in `chain_x` alone."""
        return int(self.mask_x.sum())

    @property
    def n_interface_residues_y(self) -> int:
        """Interface residues in `chain_y` alone."""
        return int(self.mask_y.sum())

    @property
    def contact_pairs(self) -> List[Tuple[int, int]]:
        """`(i, j)` index pairs within the cutoff, positional not residue-numbered."""
        rows, cols = np.where(self.contact_mask)
        return [(int(i), int(j)) for i, j in zip(rows, cols)]


def detect_interface(
    chain_x: ChainCoords,
    chain_y: ChainCoords,
    dist_cutoff: float = DIST_CUTOFF,
) -> InterfaceContacts:
    """
    Detect the interface between an ordered pair of chains by CB-CB distance.

    Glycine has no CB atom; `extract_chain_coords` has already substituted its
    CA, so callers never handle glycine themselves.

    Memory is `O(nx * ny)`: a 1000-residue pair is about 8 MB of float32 distances
    plus the `(nx, ny, 3)` broadcast temporary, which is fine on Colab.

    Args:
        chain_x:     First chain of the ordered pair (rows).
        chain_y:     Second chain of the ordered pair (columns).
        dist_cutoff: Contact cutoff in Angstrom, tested inclusively.

    Returns:
        An `InterfaceContacts` carrying the distance matrix, the pairwise contact
        mask, both per-chain interface masks and both counts.

    Example
    -------
    >>> a = ChainCoords('A', np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], dtype=np.float32),
    ...                 np.array([1, 2], dtype=np.int32), np.array(['ALA', 'ALA']),
    ...                 np.array([90.0, 85.0], dtype=np.float32))
    >>> b = ChainCoords('B', np.array([[3.0, 0.0, 0.0], [20.0, 0.0, 0.0]], dtype=np.float32),
    ...                 np.array([1, 2], dtype=np.int32), np.array(['GLY', 'ALA']),
    ...                 np.array([88.0, 80.0], dtype=np.float32))
    >>> contacts = detect_interface(a, b)
    >>> contacts.n_contact_pairs          # A1-B1 at 3 A and A2-B1 at 7 A
    2
    >>> contacts.n_interface_residues     # both A residues, one B residue
    3
    >>> contacts.contact_pairs
    [(0, 0), (1, 0)]
    """
    diff = chain_x.coords[:, np.newaxis, :] - chain_y.coords[np.newaxis, :, :]
    dist_matrix = np.sqrt((diff ** 2).sum(axis=-1)).astype(np.float32)
    contact_mask = dist_matrix <= dist_cutoff

    return InterfaceContacts(
        chain_x=chain_x.chain_id,
        chain_y=chain_y.chain_id,
        dist_matrix=dist_matrix,
        contact_mask=contact_mask,
        mask_x=contact_mask.any(axis=1),
        mask_y=contact_mask.any(axis=0),
        dist_cutoff=float(dist_cutoff),
    )


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
# Score functions
# ---------------------------------------------------------------------------
# The seven reported values, all from `ipsae.py` v4 as transcribed in
# `specs/homodimer_diagnostic/formula-reference.md`.
#
# Every function takes an **ordered** chain pair (D4) and returns a dataclass
# carrying the headline value together with every intermediate a caller could
# otherwise be tempted to recompute: per-residue arrays, per-direction values,
# the `d0` and `n0` behind each one, contact counts, mean pLDDT, mean ptm, and
# the index of the argmax residue. Nothing downstream needs to re-derive a
# number that was already computed here.
#
# How the two directions combine, per `formula-reference.md`:
#
#   | value                          | per direction     | reported |
#   |--------------------------------|-------------------|----------|
#   | ipTM_d0chn                     | max over residues | max      |
#   | ipSAE_d0res / _d0chn / _d0dom  | max over residues | max      |
#   | pDockQ                         | -- symmetric --   | the value|
#   | pDockQ2                        | pooled over pairs | max      |
#   | LIS                            | pooled over pairs | mean     |
#
# Four corrections to the notebook's inline copies land here (R003, R005, R006,
# R007); each is documented on the function that carries it.

_DIRECTION_FORWARD = "xy"
_DIRECTION_REVERSE = "yx"


@dataclass(frozen=True, eq=False)
class ResidueProfile:
    """
    One direction of a per-residue score: the profile plus what produced it.

    "One direction" means the rows of a single PAE block. `values[i]` is the
    score of residue `i` of `chain_row` measured against the whole of
    `chain_col`, and the direction's reported value is `values.max()` -- the
    score `ipsae.py` quotes is literally one residue's number
    (`ipsae_v4.py:827-844`), which is why `argmax_index` is carried.

    Attributes:
        chain_row:     Chain supplying the rows, i.e. the alignment frame.
        chain_col:     Chain supplying the columns.
        values:        `(n_row,)` per-residue score. `0.0` for a residue with no
                       contributing pair, matching `ipsae_v4.py:747, 796, 800`.
        d0:            `(n_row,)` the `d0` actually used for each row. Constant
                       across rows for the `d0chn` and `d0dom` variants, genuinely
                       per-residue for `d0res`.
        n0:            `(n_row,)` the `L` that produced each `d0`, so a caller can
                       report the `n0chn` / `n0dom` / `n0res` that `ipsae.py`
                       prints alongside the score.
        n_valid_pairs: `(n_row,)` number of partner residues that contributed to
                       each row's mean. Equals `n0` for `d0res`; for the other
                       variants the two are different quantities and are kept
                       apart deliberately.

    Example
    -------
    >>> profile = ResidueProfile('A', 'B', np.array([0.2, 0.8, 0.5]),
    ...                          np.full(3, 6.76), np.full(3, 344),
    ...                          np.array([4, 9, 7]))
    >>> profile.score, profile.argmax_index
    (0.8, 1)
    >>> profile.d0_at_argmax, profile.n0_at_argmax
    (6.76, 344)
    """

    chain_row: str
    chain_col: str
    values: np.ndarray
    d0: np.ndarray
    n0: np.ndarray
    n_valid_pairs: np.ndarray

    @property
    def score(self) -> float:
        """This direction's reported value: the maximum over residues."""
        return float(self.values.max()) if self.values.size else 0.0

    @property
    def argmax_index(self) -> int:
        """Positional index of the residue whose value the direction reports.

        `-1` for an empty chain. Ties go to the lowest index, as `np.argmax` and
        `ipsae.py` both do.
        """
        return int(np.argmax(self.values)) if self.values.size else -1

    @property
    def d0_at_argmax(self) -> float:
        """The `d0` of the argmax residue. This is the `d0` `ipsae.py` prints."""
        return float(self.d0[self.argmax_index]) if self.values.size else 0.0

    @property
    def n0_at_argmax(self) -> int:
        """The `n0` of the argmax residue. This is the `n0` `ipsae.py` prints."""
        return int(self.n0[self.argmax_index]) if self.values.size else 0

    @property
    def n_interface_residues(self) -> int:
        """Rows with at least one contributing pair."""
        return int((self.n_valid_pairs > 0).sum())


@dataclass(frozen=True, eq=False)
class DirectionalPair:
    """
    A per-residue score in both directions, plus the combination `ipsae.py` reports.

    PAE is asymmetric, so the two directions are two measurements rather than one
    measurement seen twice. Collapsing them to the reported `max` early is what
    hides bugs like R003 and R007, so both are kept and `delta` is offered as a
    first-class diagnostic (R008).

    Attributes:
        name:    Key into `THRESHOLDS`, e.g. `'ipsae_d0res'`.
        chain_x: First chain of the ordered pair.
        chain_y: Second chain of the ordered pair.
        forward: The `x -> y` profile.
        reverse: The `y -> x` profile.

    Example
    -------
    >>> f = ResidueProfile('A', 'B', np.array([0.4, 0.9]), np.full(2, 6.0),
    ...                    np.full(2, 300), np.array([5, 7]))
    >>> r = ResidueProfile('B', 'A', np.array([0.7, 0.3]), np.full(2, 6.0),
    ...                    np.full(2, 300), np.array([6, 2]))
    >>> pair = DirectionalPair('ipsae_d0chn', 'A', 'B', f, r)
    >>> pair.score, pair.winner
    (0.9, 'xy')
    >>> round(pair.delta, 6)
    0.2
    """

    name: str
    chain_x: str
    chain_y: str
    forward: ResidueProfile
    reverse: ResidueProfile

    @property
    def forward_score(self) -> float:
        """The `x -> y` value."""
        return self.forward.score

    @property
    def reverse_score(self) -> float:
        """The `y -> x` value."""
        return self.reverse.score

    @property
    def score(self) -> float:
        """The reported value: `max(x -> y, y -> x)` (`ipsae_v4.py:851, 859, 867, 885`)."""
        return max(self.forward_score, self.reverse_score)

    @property
    def winner(self) -> str:
        """`'xy'` or `'yx'`, whichever direction supplied `score`.

        An exact tie resolves to the forward direction. `ipsae.py` resolves a tie
        to its own `chain1 > chain2` direction (`ipsae_v4.py:852-854`); the two
        rules can only disagree about *which* identical number is quoted, never
        about the number.
        """
        return _DIRECTION_FORWARD if self.forward_score >= self.reverse_score else _DIRECTION_REVERSE

    @property
    def winning_profile(self) -> ResidueProfile:
        """The `ResidueProfile` that supplied `score`."""
        return self.forward if self.winner == _DIRECTION_FORWARD else self.reverse

    @property
    def delta(self) -> float:
        """`|x -> y  -  y -> x|`. Zero on a perfectly symmetric pair (R008)."""
        return abs(self.forward_score - self.reverse_score)

    @property
    def d0(self) -> float:
        """The `d0` of the winning direction's argmax residue (R003)."""
        return self.winning_profile.d0_at_argmax

    @property
    def n0(self) -> int:
        """The `n0` of the winning direction's argmax residue (R003)."""
        return self.winning_profile.n0_at_argmax

    @property
    def argmax_index(self) -> int:
        """Positional index, within the winning direction's chain, of the residue
        whose value is the reported score."""
        return self.winning_profile.argmax_index


def _row_means(
    block: np.ndarray,
    mask: np.ndarray,
    d0: np.ndarray | float,
) -> np.ndarray:
    """
    Per-row mean of `ptm(block, d0)` over the entries `mask` selects.

    Rows selecting nothing get `0.0`, which is `ipsae.py`'s
    `... .mean() if valid.any() else 0.0` (`ipsae_v4.py:747, 796, 800`) without
    the Python-level loop.

    Args:
        block: `(n_row, n_col)` PAE block.
        mask:  `(n_row, n_col)` bool, which entries contribute.
        d0:    Scalar, or `(n_row, 1)` for a per-residue `d0`.

    Returns:
        `(n_row,)` float64.
    """
    ptm = ptm_func(block.astype(np.float64), d0)
    counts = mask.sum(axis=1)
    totals = np.where(mask, ptm, 0.0).sum(axis=1)
    means = np.zeros(block.shape[0], dtype=np.float64)
    np.divide(totals, counts, out=means, where=counts > 0)
    return means


def compute_iptm_d0chn(pair: ChainPairPAE) -> DirectionalPair:
    """
    ipTM_d0chn: the TM-transformed PAE averaged over the **whole** partner chain.

    No PAE cutoff at all -- every partner residue contributes
    (`valid_pairs_iptm = (chains == chain2)`, `ipsae_v4.py:736`). That is the only
    thing separating it from `ipSAE_d0chn`, which shares its `d0`.

    `d0chn = d0_scalar(nx + ny)` (`ipsae_v4.py:731-732`, the **scalar** helper).
    Reported value is `max` over the two directions (`ipsae_v4.py:851`).

    Args:
        pair: The ordered chain pair's PAE quadrants.

    Returns:
        A `DirectionalPair` named `'iptm_d0chn'`.

    Note:
        **This is not AlphaFold's own ipTM.** `ipsae.py` calls that `ipTM_af` and
        reads it from the model's summary file (`ipsae_v4.py:944-946`); AFDB does
        not expose it on these endpoints. This value is a reimplementation from
        the PAE matrix and will not equal a number quoted by AlphaFold (R004).

    Example
    -------
    >>> block_xy = np.array([[2.0, 30.0], [1.0, 1.0]], dtype=np.float32)
    >>> block_yx = np.array([[3.0, 2.0], [30.0, 30.0]], dtype=np.float32)
    >>> pair = ChainPairPAE('A', 'B', block_xy, block_yx,
    ...                     np.zeros((2, 2), dtype=np.float32),
    ...                     np.zeros((2, 2), dtype=np.float32))
    >>> result = compute_iptm_d0chn(pair)
    >>> result.n0                      # n0chn = nx + ny, never per residue
    4
    >>> round(result.d0, 6)            # d0_scalar(4) -> the 1.0 floor
    1.0
    >>> round(result.score, 6)         # residue A2: both partners at PAE 1.0
    0.5
    >>> result.winner, result.argmax_index
    ('xy', 1)
    """
    d0chn = d0_scalar(pair.nx + pair.ny)
    n0chn = pair.nx + pair.ny

    def _profile(block: np.ndarray, row_chain: str, col_chain: str) -> ResidueProfile:
        n_row, n_col = block.shape
        # No cutoff: the mask is all-True, i.e. the whole partner chain.
        all_pairs = np.ones(block.shape, dtype=bool)
        return ResidueProfile(
            chain_row=row_chain,
            chain_col=col_chain,
            values=_row_means(block, all_pairs, d0chn),
            d0=np.full(n_row, d0chn, dtype=np.float64),
            n0=np.full(n_row, n0chn, dtype=np.int64),
            n_valid_pairs=np.full(n_row, n_col, dtype=np.int64),
        )

    return DirectionalPair(
        name="iptm_d0chn",
        chain_x=pair.chain_x,
        chain_y=pair.chain_y,
        forward=_profile(pair.block_xy, pair.chain_x, pair.chain_y),
        reverse=_profile(pair.block_yx, pair.chain_y, pair.chain_x),
    )


@dataclass(frozen=True, eq=False)
class IPSAEResult:
    """
    The three ipSAE variants, which differ only in what `L` feeds `d0`.

    All three share one mask -- inter-chain pairs with `pae < pae_cutoff`
    (`ipsae_v4.py:737, 783`) -- so they are computed together and the mask is
    handed back rather than rebuilt by each caller.

    Attributes:
        chain_x:       First chain of the ordered pair.
        chain_y:       Second chain of the ordered pair.
        pae_cutoff:    The cutoff in force, tested strictly (`pae < cutoff`).
        d0res:         Per-residue `d0`, from that residue's own valid-pair count.
        d0chn:         Chain-pair `d0`, from `nx + ny`.
        d0dom:         Domain `d0`, from the interacting residue count **of that
                       direction** (R003).
        valid_xy:      `(nx, ny)` bool mask of contributing pairs, `x -> y`.
        valid_yx:      `(ny, nx)` bool mask of contributing pairs, `y -> x`.
        n0dom_xy:      Interacting-residue count behind `d0dom` for `x -> y`.
        n0dom_yx:      The same for `y -> x`. Differs from `n0dom_xy` in general.
        d0dom_xy:      `d0_scalar(n0dom_xy)`.
        d0dom_yx:      `d0_scalar(n0dom_yx)`.
        n0chn:         `nx + ny`, shared by every direction and variant.
        d0chn_value:   `d0_scalar(n0chn)`.
    """

    chain_x: str
    chain_y: str
    pae_cutoff: float
    d0res: DirectionalPair
    d0chn: DirectionalPair
    d0dom: DirectionalPair
    valid_xy: np.ndarray
    valid_yx: np.ndarray
    n0dom_xy: int
    n0dom_yx: int
    d0dom_xy: float
    d0dom_yx: float
    n0chn: int
    d0chn_value: float

    @property
    def n0dom_delta(self) -> int:
        """`|n0dom(x->y) - n0dom(y->x)|`, the mechanism behind R003.

        Zero on a symmetric pair, non-zero on a real heterodimer, which is why
        the bug was invisible while only homodimers were tested (R008).
        """
        return abs(self.n0dom_xy - self.n0dom_yx)

    @property
    def variants(self) -> Dict[str, DirectionalPair]:
        """`{THRESHOLDS key: DirectionalPair}` for all three variants."""
        return {
            "ipsae_d0res": self.d0res,
            "ipsae_d0chn": self.d0chn,
            "ipsae_d0dom": self.d0dom,
        }


def compute_ipsae(
    pair: ChainPairPAE,
    pae_cutoff: float = PAE_CUTOFF,
) -> IPSAEResult:
    """
    The three ipSAE variants for one ordered chain pair.

    ipSAE is ipTM restricted to inter-chain pairs with `pae < pae_cutoff`, and
    the three variants differ only in the `L` that sets `d0`:

    | variant | `L` | helper | citation |
    |---------|-----|--------|----------|
    | `d0res` | that residue's own count of valid pairs | `d0_array` | `ipsae_v4.py:786-787` |
    | `d0chn` | `nx + ny`, the whole chain pair | `d0_scalar` | `ipsae_v4.py:731-732` |
    | `d0dom` | interacting residues in **this direction** | `d0_scalar` | `ipsae_v4.py:775-778` |

    Args:
        pair:       The ordered chain pair's PAE quadrants.
        pae_cutoff: Strict inter-chain PAE cutoff.

    Returns:
        An `IPSAEResult`.

    Note:
        **R003 -- `d0dom` is per direction.** For `x -> y` with block `P`,
        `n0dom = (P < cutoff).any(axis=1).sum() + (P < cutoff).any(axis=0).sum()`,
        both terms read off *that direction's own block*. The notebook's inline
        copy computed `n0dom` once from the `A -> B` block and reused it for
        `B -> A`, so the reverse direction's whole per-residue array carried the
        wrong `d0`. On the heterodimer fixture the two counts are 250 (`A -> B`)
        and 252 (`B -> A`); on the homodimer both are 342, which is why the bug
        survived homodimer-only testing. Measured effect of the bug on the
        heterodimer: `ipSAE_d0dom(B -> A)` becomes 0.766806 instead of 0.768065,
        an error of 0.0013 and so outside the +/-0.001 tolerance. It does not
        move the *reported* value on either fixture only because `A -> B` wins
        the `max` in both cases. The reported `n0dom` / `d0dom` are the winning
        direction's (`ipsae_v4.py:870-883`).

    Note:
        **R005 -- `d0res` uses `d0_array`, not `d0_scalar`.** The two differ at
        exactly `L == 27` (1.0 vs 1.038891). `d0chn` and `d0dom` keep the scalar
        helper, because `ipsae.py` genuinely uses both.

    Example
    -------
    >>> block_xy = np.array([[2.0, 30.0], [1.0, 1.0]], dtype=np.float32)
    >>> block_yx = np.array([[3.0, 2.0], [30.0, 30.0]], dtype=np.float32)
    >>> pair = ChainPairPAE('A', 'B', block_xy, block_yx,
    ...                     np.zeros((2, 2), dtype=np.float32),
    ...                     np.zeros((2, 2), dtype=np.float32))
    >>> result = compute_ipsae(pair)

    Row A1 keeps one partner, row A2 keeps two; rows B1 / B2 keep two and none:

    >>> result.valid_xy.sum(axis=1).tolist(), result.valid_yx.sum(axis=1).tolist()
    ([1, 2], [2, 0])

    R003 in miniature -- the two directions disagree about `n0dom`:

    >>> result.n0dom_xy, result.n0dom_yx
    (4, 3)
    >>> result.n0dom_delta
    1

    R005 in miniature -- `d0res` is per residue, `d0chn` and `d0dom` are not:

    >>> result.d0res.forward.d0.round(6).tolist()
    [1.0, 1.0]
    >>> sorted({round(float(v), 6) for v in result.d0chn.forward.d0})
    [1.0]
    >>> round(result.d0dom.score, 6) == round(result.d0chn.score, 6)
    True
    """
    n0chn = pair.nx + pair.ny
    d0chn = d0_scalar(n0chn)

    valid_xy = pair.block_xy < pae_cutoff
    valid_yx = pair.block_yx < pae_cutoff

    # R003: rows-with-any plus cols-with-any of *this direction's own* block.
    n0dom_xy = int(valid_xy.any(axis=1).sum()) + int(valid_xy.any(axis=0).sum())
    n0dom_yx = int(valid_yx.any(axis=1).sum()) + int(valid_yx.any(axis=0).sum())
    d0dom_xy = d0_scalar(n0dom_xy)
    d0dom_yx = d0_scalar(n0dom_yx)

    def _fixed_d0_profile(
        block: np.ndarray,
        mask: np.ndarray,
        row_chain: str,
        col_chain: str,
        d0: float,
        n0: int,
    ) -> ResidueProfile:
        n_row = block.shape[0]
        return ResidueProfile(
            chain_row=row_chain,
            chain_col=col_chain,
            values=_row_means(block, mask, d0),
            d0=np.full(n_row, d0, dtype=np.float64),
            n0=np.full(n_row, n0, dtype=np.int64),
            n_valid_pairs=mask.sum(axis=1).astype(np.int64),
        )

    def _per_residue_d0_profile(
        block: np.ndarray,
        mask: np.ndarray,
        row_chain: str,
        col_chain: str,
    ) -> ResidueProfile:
        n0res = mask.sum(axis=1).astype(np.int64)
        d0res = d0_array(n0res)  # R005: the array helper, deliberately
        return ResidueProfile(
            chain_row=row_chain,
            chain_col=col_chain,
            values=_row_means(block, mask, d0res[:, np.newaxis]),
            d0=d0res,
            n0=n0res,
            n_valid_pairs=n0res,
        )

    x, y = pair.chain_x, pair.chain_y
    return IPSAEResult(
        chain_x=x,
        chain_y=y,
        pae_cutoff=float(pae_cutoff),
        d0res=DirectionalPair(
            name="ipsae_d0res",
            chain_x=x,
            chain_y=y,
            forward=_per_residue_d0_profile(pair.block_xy, valid_xy, x, y),
            reverse=_per_residue_d0_profile(pair.block_yx, valid_yx, y, x),
        ),
        d0chn=DirectionalPair(
            name="ipsae_d0chn",
            chain_x=x,
            chain_y=y,
            forward=_fixed_d0_profile(pair.block_xy, valid_xy, x, y, d0chn, n0chn),
            reverse=_fixed_d0_profile(pair.block_yx, valid_yx, y, x, d0chn, n0chn),
        ),
        d0dom=DirectionalPair(
            name="ipsae_d0dom",
            chain_x=x,
            chain_y=y,
            forward=_fixed_d0_profile(pair.block_xy, valid_xy, x, y, d0dom_xy, n0dom_xy),
            reverse=_fixed_d0_profile(pair.block_yx, valid_yx, y, x, d0dom_yx, n0dom_yx),
        ),
        valid_xy=valid_xy,
        valid_yx=valid_yx,
        n0dom_xy=n0dom_xy,
        n0dom_yx=n0dom_yx,
        d0dom_xy=d0dom_xy,
        d0dom_yx=d0dom_yx,
        n0chn=n0chn,
        d0chn_value=d0chn,
    )


@dataclass(frozen=True, eq=False)
class PDockQResult:
    """
    pDockQ for one chain pair. Symmetric, so there is nothing to combine.

    Attributes:
        chain_x:              First chain of the ordered pair.
        chain_y:              Second chain of the ordered pair.
        score:                The reported pDockQ.
        n_contact_pairs:      `npairs`, the count of residue **pairs** within the
                              cutoff. Not the interface residue count.
        n_interface_residues: `nres`, both chains summed. `ipsae.py` computes this
                              (`ipsae_v4.py:662`) and only prints it; it does not
                              enter the score.
        mean_plddt:           Mean pLDDT over the union of interface residues.
        x:                    `mean_plddt * log10(n_contact_pairs)`, the sigmoid's
                              argument.
        mask_x:               `(nx,)` bool, interface residues of `chain_x`.
        mask_y:               `(ny,)` bool, interface residues of `chain_y`.
        dist_cutoff:          Contact cutoff in Angstrom, tested inclusively.
        symmetric:            Always `True`, and checked rather than assumed --
                              see `compute_pdockq`.
    """

    chain_x: str
    chain_y: str
    score: float
    n_contact_pairs: int
    n_interface_residues: int
    mean_plddt: float
    x: float
    mask_x: np.ndarray
    mask_y: np.ndarray
    dist_cutoff: float
    symmetric: bool = True


def _pdockq_from(
    contact_mask: np.ndarray,
    plddt_row: np.ndarray,
    plddt_col: np.ndarray,
) -> Tuple[float, int, int, float, float]:
    """
    pDockQ from one orientation of a contact mask.

    Split out so `compute_pdockq` can run it in both orientations and check the
    symmetry claim by execution instead of asserting it in a comment.

    Args:
        contact_mask: `(n_row, n_col)` bool.
        plddt_row:    `(n_row,)` pLDDT of the row chain.
        plddt_col:    `(n_col,)` pLDDT of the column chain.

    Returns:
        `(score, n_contact_pairs, n_interface_residues, mean_plddt, x)`.
    """
    n_pairs = int(contact_mask.sum())
    mask_row = contact_mask.any(axis=1)
    mask_col = contact_mask.any(axis=0)
    n_residues = int(mask_row.sum()) + int(mask_col.sum())
    if n_pairs == 0:
        # R006: `ipsae_v4.py:669` short-circuits to 0.0, not to the sigmoid's
        # x -> -inf limit of 0.018.
        return 0.0, 0, n_residues, 0.0, 0.0
    interface_plddt = np.concatenate(
        [np.asarray(plddt_row, dtype=np.float64)[mask_row],
         np.asarray(plddt_col, dtype=np.float64)[mask_col]]
    )
    mean_plddt = float(interface_plddt.mean())
    x = mean_plddt * float(np.log10(n_pairs))
    score = 0.724 / (1.0 + np.exp(-0.052 * (x - 152.611))) + 0.018
    return float(score), n_pairs, n_residues, mean_plddt, x


def compute_pdockq(
    contacts: InterfaceContacts,
    plddt_x: np.ndarray,
    plddt_y: np.ndarray,
) -> PDockQResult:
    """
    pDockQ (Bryant 2022) as `ipsae.py` computes it.

    `x = mean_plddt * log10(npairs)`, then
    `0.724 / (1 + exp(-0.052 * (x - 152.611))) + 0.018` (`ipsae_v4.py:663-665`).

    Two things are easy to get wrong and are therefore spelled out:

    - `npairs` is the count of **contact pairs**, not of interface residues
      (`ipsae_v4.py:653`). `CLAUDE.md`'s version uses the residue count, which for
      a typical interface makes `x` several times too small.
    - `mean_plddt` is the unweighted mean over the **union** of interface
      residues from both chains (`ipsae_v4.py:663`), so a residue with forty
      contacts counts exactly once.

    Args:
        contacts: The ordered pair's `InterfaceContacts`.
        plddt_x:  `(nx,)` pLDDT of `contacts.chain_x`, CB-atom values.
        plddt_y:  `(ny,)` pLDDT of `contacts.chain_y`.

    Returns:
        A `PDockQResult`.

    Raises:
        ValueError: If a pLDDT array's length does not match its chain, which
            would silently score the wrong residues.
        AssertionError: If the two orientations disagree -- see the symmetry note.

    Note:
        **R006 -- zero contacts return `0.0`.** `ipsae.py` short-circuits to
        `0.0` (`ipsae_v4.py:669`) rather than to `0.018`, the `x -> -inf` limit of
        Bryant's sigmoid that Bryant's own code returns. This is a real
        disagreement between two defensible conventions, not a transcription
        error, and it is settled in favour of `ipsae.py` because +/-0.001
        agreement with `ipsae.py` is this project's acceptance criterion. A
        reader comparing against Bryant's published implementation on a
        non-interacting pair will see 0.0 here and 0.018 there.

    Note:
        **Symmetry is checked, not assumed.** `npairs` and the interface residue
        set are both invariant under swapping the chains, which is why `ipsae.py`
        prints pDockQ unmaxed (`ipsae_v4.py:989`) while it maxes pDockQ2. Rather
        than take that on trust, the score is computed in both orientations and
        the two are required to agree; R008 then reports pDockQ as symmetric
        instead of showing it an empty directional diff.

    Example
    -------
    >>> a = ChainCoords('A', np.array([[0.0, 0.0, 0.0]], dtype=np.float32),
    ...                 np.array([1], dtype=np.int32), np.array(['ALA']),
    ...                 np.array([90.0], dtype=np.float32))
    >>> b = ChainCoords('B', np.array([[3.0, 0.0, 0.0], [40.0, 0.0, 0.0]], dtype=np.float32),
    ...                 np.array([1, 2], dtype=np.int32), np.array(['ALA', 'ALA']),
    ...                 np.array([80.0, 70.0], dtype=np.float32))
    >>> contacts = detect_interface(a, b)
    >>> result = compute_pdockq(contacts, a.plddt, b.plddt)
    >>> result.n_contact_pairs, result.n_interface_residues
    (1, 2)
    >>> result.mean_plddt          # (90 + 80) / 2, the far B residue excluded
    85.0
    >>> round(result.x, 6)         # log10(1) == 0
    0.0
    >>> round(result.score, 6)     # deep in the sigmoid's lower tail
    0.018259

    R006, no contacts at all:

    >>> far = ChainCoords('B', np.array([[99.0, 0.0, 0.0]], dtype=np.float32),
    ...                   np.array([1], dtype=np.int32), np.array(['ALA']),
    ...                   np.array([80.0], dtype=np.float32))
    >>> compute_pdockq(detect_interface(a, far), a.plddt, far.plddt).score
    0.0
    """
    nx, ny = contacts.contact_mask.shape
    if len(plddt_x) != nx or len(plddt_y) != ny:
        raise ValueError(
            f"pLDDT lengths {len(plddt_x)}/{len(plddt_y)} do not match the contact "
            f"matrix {nx}x{ny} for chains {contacts.chain_x}/{contacts.chain_y}."
        )

    score, n_pairs, n_residues, mean_plddt, x = _pdockq_from(
        contacts.contact_mask, plddt_x, plddt_y
    )
    # Same score, chains swapped. If this ever disagrees, the symmetry claim that
    # lets `ipsae.py` print pDockQ unmaxed is false and every caller must know.
    reverse = _pdockq_from(contacts.contact_mask.T, plddt_y, plddt_x)
    assert abs(reverse[0] - score) < 1e-12 and reverse[1] == n_pairs, (
        f"pDockQ is not symmetric for {contacts.chain_x}/{contacts.chain_y}: "
        f"{score} vs {reverse[0]}."
    )

    return PDockQResult(
        chain_x=contacts.chain_x,
        chain_y=contacts.chain_y,
        score=score,
        n_contact_pairs=n_pairs,
        n_interface_residues=n_residues,
        mean_plddt=mean_plddt,
        x=x,
        mask_x=contacts.mask_x,
        mask_y=contacts.mask_y,
        dist_cutoff=contacts.dist_cutoff,
    )


@dataclass(frozen=True, eq=False)
class PDockQ2Direction:
    """
    One direction of pDockQ2, plus its per-residue decomposition.

    Attributes:
        chain_row:           Chain supplying the rows of the PAE block read.
        chain_col:           Chain supplying the columns.
        score:               This direction's pDockQ2.
        n_contact_pairs:     `npairs`, shared by both directions.
        mean_plddt:          Mean pLDDT over the union of interface residues.
                             Symmetric: identical in both directions.
        mean_ptm:            `sum(ptm(PAE, d0=10)) / npairs` over contact pairs,
                             read from **this direction's** PAE block. This is the
                             only asymmetric ingredient (R007).
        x:                   `mean_plddt * mean_ptm`.
        mean_ptm_by_residue: `(n_row,)` the mean of `ptm(PAE, d0=10)` over each
                             row residue's own contact pairs, `np.nan` for a
                             residue with no contacts. The per-residue quantity
                             pDockQ2 pools away; exposed for R074's contact
                             quality view.
        contact_counts:      `(n_row,)` contacts per row residue.
    """

    chain_row: str
    chain_col: str
    score: float
    n_contact_pairs: int
    mean_plddt: float
    mean_ptm: float
    x: float
    mean_ptm_by_residue: np.ndarray
    contact_counts: np.ndarray

    @property
    def n_interface_residues_row(self) -> int:
        """Row-chain residues with at least one contact."""
        return int((self.contact_counts > 0).sum())


@dataclass(frozen=True, eq=False)
class PDockQ2Result:
    """
    pDockQ2 in both directions, and the `max` `ipsae.py` reports.

    Attributes:
        chain_x: First chain of the ordered pair.
        chain_y: Second chain of the ordered pair.
        forward: The `x -> y` direction.
        reverse: The `y -> x` direction.
        name:    Key into `THRESHOLDS`.
    """

    chain_x: str
    chain_y: str
    forward: PDockQ2Direction
    reverse: PDockQ2Direction
    name: str = "pdockq2"

    @property
    def forward_score(self) -> float:
        """The `x -> y` value."""
        return self.forward.score

    @property
    def reverse_score(self) -> float:
        """The `y -> x` value."""
        return self.reverse.score

    @property
    def score(self) -> float:
        """The reported value: `max(x -> y, y -> x)` (`ipsae_v4.py:977, 990`)."""
        return max(self.forward_score, self.reverse_score)

    @property
    def winner(self) -> str:
        """`'xy'` or `'yx'`, whichever direction supplied `score`; ties go forward."""
        return _DIRECTION_FORWARD if self.forward_score >= self.reverse_score else _DIRECTION_REVERSE

    @property
    def winning_direction(self) -> PDockQ2Direction:
        """The `PDockQ2Direction` that supplied `score`."""
        return self.forward if self.winner == _DIRECTION_FORWARD else self.reverse

    @property
    def delta(self) -> float:
        """`|x -> y  -  y -> x|`. Non-zero whenever the PAE block is asymmetric (R007)."""
        return abs(self.forward_score - self.reverse_score)


def _pdockq2_direction(
    block: np.ndarray,
    contact_mask: np.ndarray,
    row_chain: str,
    col_chain: str,
    mean_plddt: float,
    n_pairs: int,
) -> PDockQ2Direction:
    """
    One direction of pDockQ2.

    Args:
        block:        `(n_row, n_col)` PAE, rows aligned on `row_chain`.
        contact_mask: `(n_row, n_col)` bool, same orientation as `block`.
        row_chain:    Row chain id.
        col_chain:    Column chain id.
        mean_plddt:   Interface mean pLDDT, already computed by pDockQ.
        n_pairs:      Contact pair count.

    Returns:
        A `PDockQ2Direction`.
    """
    counts = contact_mask.sum(axis=1).astype(np.int64)
    ptm = ptm_func(block.astype(np.float64), 10.0)  # fixed d0, `ipsae_v4.py:687`
    row_totals = np.where(contact_mask, ptm, 0.0).sum(axis=1)

    by_residue = np.full(block.shape[0], np.nan, dtype=np.float64)
    np.divide(row_totals, counts, out=by_residue, where=counts > 0)

    if n_pairs == 0:
        # R006: `ipsae_v4.py:700` returns 0.0, not Zhu's 0.005 sigmoid minimum.
        return PDockQ2Direction(
            chain_row=row_chain,
            chain_col=col_chain,
            score=0.0,
            n_contact_pairs=0,
            mean_plddt=0.0,
            mean_ptm=0.0,
            x=0.0,
            mean_ptm_by_residue=by_residue,
            contact_counts=counts,
        )

    mean_ptm = float(row_totals.sum() / n_pairs)
    x = mean_plddt * mean_ptm
    score = float(1.31 / (1.0 + np.exp(-0.075 * (x - 84.733))) + 0.005)
    return PDockQ2Direction(
        chain_row=row_chain,
        chain_col=col_chain,
        score=score,
        n_contact_pairs=n_pairs,
        mean_plddt=mean_plddt,
        mean_ptm=mean_ptm,
        x=x,
        mean_ptm_by_residue=by_residue,
        contact_counts=counts,
    )


def compute_pdockq2(
    contacts: InterfaceContacts,
    pair: ChainPairPAE,
    plddt_x: np.ndarray,
    plddt_y: np.ndarray,
) -> PDockQ2Result:
    """
    pDockQ2 (Zhu 2023) as `ipsae.py` computes it, in both directions.

    `mean_ptm` is the mean of `ptm(PAE, d0=10)` over the contact pairs,
    `x = mean_plddt * mean_ptm`, and
    `1.31 / (1 + exp(-0.075 * (x - 84.733))) + 0.005` (`ipsae_v4.py:686-695`).
    `mean_plddt` is pDockQ's, over the union of interface residues from both
    chains -- an unweighted mean over residues, not over contact pairs.

    Args:
        contacts: The ordered pair's `InterfaceContacts`.
        pair:     The same ordered pair's PAE quadrants. Must name the same two
                  chains in the same order.
        plddt_x:  `(nx,)` pLDDT of `chain_x`, CB-atom values.
        plddt_y:  `(ny,)` pLDDT of `chain_y`.

    Returns:
        A `PDockQ2Result` carrying both directions.

    Raises:
        ValueError: If the contacts and the PAE pair disagree about the chains or
            their lengths, or if a pLDDT array's length is wrong.

    Note:
        **R007 -- pDockQ2 is directional and the reported value is a `max`.**
        `mean_ptm` reads only the `x -> y` PAE block (`ipsae_v4.py:686`) while
        `mean_plddt` is symmetric, so the two directions differ whenever the PAE
        matrix is asymmetric. `ipsae.py` reports
        `max(pDockQ2[A][B], pDockQ2[B][A])` (`ipsae_v4.py:977, 990`). The
        notebook's inline copy is passed a single block, `pae_AB`, and reports
        the `A -> B` value alone. On the heterodimer fixture the two directions
        are 0.705404 (`A -> B`) and 0.685271 (`B -> A`), a gap of 0.020; the
        notebook's number is correct there only because `A -> B` happens to be
        the larger. Reverse the chain order and the same code is wrong by 0.020,
        twenty times the tolerance.

    Note:
        **R006 -- zero contacts return `0.0`**, per `ipsae_v4.py:700`, rather than
        the 0.005 minimum of Zhu's sigmoid that Zhu's own code returns. Same
        reasoning as `compute_pdockq`.

    Note:
        **R074 -- `mean_ptm_by_residue`.** Each direction also carries the mean
        `ptm(PAE, d0=10)` over each row residue's own contact pairs, which is the
        per-residue quantity pDockQ2 pools into one number. `np.nan` marks a
        residue with no contacts, so a colour map can leave non-interface
        residues unpainted rather than painting them a misleading zero.

    Example
    -------
    >>> a = ChainCoords('A', np.array([[0.0, 0.0, 0.0], [40.0, 0.0, 0.0]], dtype=np.float32),
    ...                 np.array([1, 2], dtype=np.int32), np.array(['ALA', 'ALA']),
    ...                 np.array([90.0, 60.0], dtype=np.float32))
    >>> b = ChainCoords('B', np.array([[3.0, 0.0, 0.0], [41.0, 0.0, 0.0]], dtype=np.float32),
    ...                 np.array([1, 2], dtype=np.int32), np.array(['ALA', 'ALA']),
    ...                 np.array([80.0, 70.0], dtype=np.float32))
    >>> contacts = detect_interface(a, b)
    >>> contacts.n_contact_pairs
    2
    >>> block_xy = np.array([[0.0, 20.0], [20.0, 10.0]], dtype=np.float32)
    >>> block_yx = np.array([[10.0, 20.0], [20.0, 0.0]], dtype=np.float32)
    >>> pae = ChainPairPAE('A', 'B', block_xy, block_yx,
    ...                    np.zeros((2, 2), dtype=np.float32),
    ...                    np.zeros((2, 2), dtype=np.float32))
    >>> result = compute_pdockq2(contacts, pae, a.plddt, b.plddt)

    Contacts are A1-B1 and A2-B2. Forward reads PAE 0.0 and 10.0, reverse reads
    10.0 and 0.0 -- the same pairs, different measurements, so the same
    `mean_ptm` here by construction:

    >>> round(result.forward.mean_ptm, 6), round(result.reverse.mean_ptm, 6)
    (0.75, 0.75)
    >>> result.forward.mean_ptm_by_residue.round(6).tolist()
    [1.0, 0.5]
    >>> result.reverse.mean_ptm_by_residue.round(6).tolist()
    [0.5, 1.0]
    >>> result.delta
    0.0

    R006, no contacts at all:

    >>> far = ChainCoords('B', np.array([[99.0, 0.0, 0.0], [99.0, 9.0, 0.0]], dtype=np.float32),
    ...                   np.array([1, 2], dtype=np.int32), np.array(['ALA', 'ALA']),
    ...                   np.array([80.0, 70.0], dtype=np.float32))
    >>> empty = compute_pdockq2(detect_interface(a, far), pae, a.plddt, far.plddt)
    >>> empty.score, empty.forward.score, empty.reverse.score
    (0.0, 0.0, 0.0)
    """
    if (contacts.chain_x, contacts.chain_y) != (pair.chain_x, pair.chain_y):
        raise ValueError(
            f"Contacts are for {contacts.chain_x}->{contacts.chain_y} but the PAE "
            f"pair is {pair.chain_x}->{pair.chain_y}; the blocks would not align."
        )
    if contacts.contact_mask.shape != pair.block_xy.shape:
        raise ValueError(
            f"Contact matrix {contacts.contact_mask.shape} and PAE block "
            f"{pair.block_xy.shape} disagree about the chain lengths."
        )
    if len(plddt_x) != pair.nx or len(plddt_y) != pair.ny:
        raise ValueError(
            f"pLDDT lengths {len(plddt_x)}/{len(plddt_y)} do not match chain "
            f"lengths {pair.nx}/{pair.ny}."
        )

    # Same union-of-interface-residues mean pLDDT pDockQ uses; `ipsae.py` reuses
    # the set built during the pDockQ pass (`ipsae_v4.py:692`).
    _, n_pairs, _, mean_plddt, _ = _pdockq_from(contacts.contact_mask, plddt_x, plddt_y)

    return PDockQ2Result(
        chain_x=pair.chain_x,
        chain_y=pair.chain_y,
        forward=_pdockq2_direction(
            pair.block_xy, contacts.contact_mask,
            pair.chain_x, pair.chain_y, mean_plddt, n_pairs,
        ),
        reverse=_pdockq2_direction(
            pair.block_yx, contacts.contact_mask.T,
            pair.chain_y, pair.chain_x, mean_plddt, n_pairs,
        ),
    )


@dataclass(frozen=True, eq=False)
class LISDirection:
    """
    One direction of LIS.

    Attributes:
        chain_row:     Chain supplying the rows of the PAE block read.
        chain_col:     Chain supplying the columns.
        score:         `mean((cutoff - PAE) / cutoff)` over sub-cutoff pairs.
        n_valid_pairs: Pairs with `PAE < cutoff`.
        n_pairs:       All inter-chain pairs in the block, `n_row * n_col`.
    """

    chain_row: str
    chain_col: str
    score: float
    n_valid_pairs: int
    n_pairs: int

    @property
    def fraction_valid(self) -> float:
        """Share of the block below the cutoff; `0.0` for an empty block."""
        return self.n_valid_pairs / self.n_pairs if self.n_pairs else 0.0


@dataclass(frozen=True, eq=False)
class LISResult:
    """
    LIS in both directions, and the **mean** `ipsae.py` reports.

    Attributes:
        chain_x:    First chain of the ordered pair.
        chain_y:    Second chain of the ordered pair.
        forward:    The `x -> y` direction.
        reverse:    The `y -> x` direction.
        lis_cutoff: The PAE cutoff in force.
        name:       Key into `THRESHOLDS`.
    """

    chain_x: str
    chain_y: str
    forward: LISDirection
    reverse: LISDirection
    lis_cutoff: float
    name: str = "lis"

    @property
    def forward_score(self) -> float:
        """The `x -> y` value."""
        return self.forward.score

    @property
    def reverse_score(self) -> float:
        """The `y -> x` value."""
        return self.reverse.score

    @property
    def score(self) -> float:
        """The reported value: the **mean** of the two directions (`ipsae_v4.py:982`).

        LIS is the one value `ipsae.py` averages rather than maxes. `max >= mean`
        always, so using `max` here -- as `CLAUDE.md` does -- systematically
        overestimates it.
        """
        return (self.forward_score + self.reverse_score) / 2.0

    @property
    def delta(self) -> float:
        """`|x -> y  -  y -> x|` (R008)."""
        return abs(self.forward_score - self.reverse_score)


def compute_lis(
    pair: ChainPairPAE,
    lis_cutoff: float = LIS_CUTOFF,
) -> LISResult:
    """
    LIS (Kim 2024) as `ipsae.py` computes it.

    Per direction, take every inter-chain pair with `PAE < 12`, rescale it to
    `(12 - PAE) / 12`, and average (`ipsae_v4.py:712-718`). The reported value is
    the **mean** of the two directions (`ipsae_v4.py:982`), not the max.

    Args:
        pair:       The ordered chain pair's PAE quadrants.
        lis_cutoff: PAE cutoff, tested strictly. Hardcoded to 12 in `ipsae.py`
                    and independent of `PAE_CUTOFF`: Kim 2024 chose 12 as the
                    value maximising AUC, and the published LIS threshold is only
                    valid at that cutoff. Changing it invalidates the threshold.

    Returns:
        A `LISResult` carrying both directions.

    Note:
        A direction with no sub-cutoff pair scores `0.0` (`ipsae_v4.py:717-720`).
        Unlike pDockQ and pDockQ2 this needs no R006 decision: `ipsae.py` and the
        published implementation agree on `0.0` here.

    Example
    -------
    >>> block_xy = np.array([[0.0, 24.0], [6.0, 24.0]], dtype=np.float32)
    >>> block_yx = np.array([[12.0, 24.0], [24.0, 24.0]], dtype=np.float32)
    >>> pair = ChainPairPAE('A', 'B', block_xy, block_yx,
    ...                     np.zeros((2, 2), dtype=np.float32),
    ...                     np.zeros((2, 2), dtype=np.float32))
    >>> result = compute_lis(pair)

    Forward keeps PAE 0.0 and 6.0 -> 1.0 and 0.5; reverse keeps nothing, because
    the cutoff is strict and 12.0 is not below 12.0:

    >>> result.forward.score, result.forward.n_valid_pairs
    (0.75, 2)
    >>> result.reverse.score, result.reverse.n_valid_pairs
    (0.0, 0)
    >>> result.score            # the MEAN, not the max
    0.375
    >>> result.delta
    0.75
    """

    def _direction(block: np.ndarray, row_chain: str, col_chain: str) -> LISDirection:
        values = block.astype(np.float64).ravel()
        valid = values[values < lis_cutoff]
        score = float(((lis_cutoff - valid) / lis_cutoff).mean()) if valid.size else 0.0
        return LISDirection(
            chain_row=row_chain,
            chain_col=col_chain,
            score=score,
            n_valid_pairs=int(valid.size),
            n_pairs=int(values.size),
        )

    return LISResult(
        chain_x=pair.chain_x,
        chain_y=pair.chain_y,
        forward=_direction(pair.block_xy, pair.chain_x, pair.chain_y),
        reverse=_direction(pair.block_yx, pair.chain_y, pair.chain_x),
        lis_cutoff=float(lis_cutoff),
    )


# --- directional transparency (R008) ---------------------------------------
# Six of the seven values are computed per direction and reported as one number.
# Collapsing them silently is exactly what let R003 and R007 survive; these two
# objects make the collapse visible so a later task can display it.

DIRECTIONAL_DELTA_TOLERANCE: float = 0.05
"""Provisional flag threshold for `|x->y - y->x|`.

A directional difference is normal and is not an error -- PAE is asymmetric by
construction. This is the point past which the difference is large enough that a
reader should be told which direction they are looking at. R008 owns the final
value and the wording that goes with it.
"""


@dataclass(frozen=True)
class DirectionalDelta:
    """
    One score's directional summary: both directions, the reported value, the gap.

    Attributes:
        name:      Key into `THRESHOLDS`.
        directional: `False` only for pDockQ, which is provably symmetric.
        combine:   How the two directions become the reported value: `'max'`,
                   `'mean'`, or `'symmetric'` when there is nothing to combine.
        score:     The reported value.
        forward:   The `x -> y` value, or `None` when not directional.
        reverse:   The `y -> x` value, or `None` when not directional.
        tolerance: The threshold `flagged` compares `delta` against.

    Example
    -------
    >>> row = DirectionalDelta('lis', True, 'mean', 0.6004, 0.6088, 0.5921)
    >>> round(row.delta, 4), row.flagged
    (0.0167, False)
    >>> DirectionalDelta('pdockq', False, 'symmetric', 0.1452).delta is None
    True
    """

    name: str
    directional: bool
    combine: str
    score: float
    forward: Optional[float] = None
    reverse: Optional[float] = None
    tolerance: float = DIRECTIONAL_DELTA_TOLERANCE

    @property
    def delta(self) -> Optional[float]:
        """`|forward - reverse|`, or `None` for a symmetric score.

        `None` rather than `0.0` deliberately: pDockQ has no directional
        difference to show, and printing a zero would suggest it was measured.
        """
        if not self.directional or self.forward is None or self.reverse is None:
            return None
        return abs(self.forward - self.reverse)

    @property
    def flagged(self) -> bool:
        """Whether `delta` exceeds `tolerance`. Always `False` when symmetric."""
        delta = self.delta
        return delta is not None and delta > self.tolerance


def directional_deltas(
    iptm_d0chn: Optional[DirectionalPair] = None,
    ipsae: Optional[IPSAEResult] = None,
    pdockq: Optional[PDockQResult] = None,
    pdockq2: Optional[PDockQ2Result] = None,
    lis: Optional[LISResult] = None,
    tolerance: float = DIRECTIONAL_DELTA_TOLERANCE,
) -> List[DirectionalDelta]:
    """
    Per-score directional summary for every result handed in (R008).

    Args:
        iptm_d0chn: From `compute_iptm_d0chn`.
        ipsae:      From `compute_ipsae`; contributes all three variants.
        pdockq:     From `compute_pdockq`; reported as symmetric, never diffed.
        pdockq2:    From `compute_pdockq2`.
        lis:        From `compute_lis`.
        tolerance:  Passed through to each row's `flagged`.

    Returns:
        One `DirectionalDelta` per supplied score, in the notebook's reporting
        order: ipTM_d0chn, the three ipSAE variants, pDockQ, pDockQ2, LIS.
        Arguments left `None` contribute no row.

    Example
    -------
    >>> f = ResidueProfile('A', 'B', np.array([0.4, 0.9]), np.full(2, 6.0),
    ...                    np.full(2, 300), np.array([5, 7]))
    >>> r = ResidueProfile('B', 'A', np.array([0.7, 0.3]), np.full(2, 6.0),
    ...                    np.full(2, 300), np.array([6, 2]))
    >>> rows = directional_deltas(iptm_d0chn=DirectionalPair('iptm_d0chn', 'A', 'B', f, r))
    >>> [(row.name, row.combine, round(row.delta, 6), row.flagged) for row in rows]
    [('iptm_d0chn', 'max', 0.2, True)]
    """
    rows: List[DirectionalDelta] = []

    def _add_pair(pair: DirectionalPair) -> None:
        rows.append(
            DirectionalDelta(
                name=pair.name,
                directional=True,
                combine="max",
                score=pair.score,
                forward=pair.forward_score,
                reverse=pair.reverse_score,
                tolerance=tolerance,
            )
        )

    if iptm_d0chn is not None:
        _add_pair(iptm_d0chn)
    if ipsae is not None:
        for variant in ("ipsae_d0res", "ipsae_d0chn", "ipsae_d0dom"):
            _add_pair(ipsae.variants[variant])
    if pdockq is not None:
        # Not diffed: `npairs` and the interface residue set are invariant under
        # swapping the chains, so there is no second direction to compare.
        rows.append(
            DirectionalDelta(
                name="pdockq",
                directional=False,
                combine="symmetric",
                score=pdockq.score,
                tolerance=tolerance,
            )
        )
    if pdockq2 is not None:
        rows.append(
            DirectionalDelta(
                name="pdockq2",
                directional=True,
                combine="max",
                score=pdockq2.score,
                forward=pdockq2.forward_score,
                reverse=pdockq2.reverse_score,
                tolerance=tolerance,
            )
        )
    if lis is not None:
        rows.append(
            DirectionalDelta(
                name="lis",
                directional=True,
                combine="mean",
                score=lis.score,
                forward=lis.forward_score,
                reverse=lis.reverse_score,
                tolerance=tolerance,
            )
        )
    return rows


# --- direction selection and the printed report (R008) ----------------------
# The user asked for "A to B by default, with a parameter to view B to A".
# Taken literally for the *reported* score that would break D2: `ipsae.py` does
# not report x -> y, it reports `max(x -> y, y -> x)` for ipTM_d0chn, the three
# ipSAE variants and pDockQ2, and the `mean` for LIS. On the heterodimer fixture
# ipSAE_d0res is 0.5555 one way and 0.7057 the other, so quoting the forward
# direction would miss AFDB's own published value by 0.15 -- 150x the +/-0.001
# tolerance.
#
# So `direction` is a **viewing** control. It selects which direction is
# highlighted in the report and which panel `plot_residue_score_profiles` draws.
# It never reaches a `compute_*` function and never changes a reported number,
# and the report says so on screen rather than only here.

DIRECTION_FORWARD: str = _DIRECTION_FORWARD
"""The `x -> y` direction: rows of `block_xy`, i.e. residues of the first chain."""

DIRECTION_REVERSE: str = _DIRECTION_REVERSE
"""The `y -> x` direction: rows of `block_yx`, i.e. residues of the second chain."""

_DIRECTION_ALIASES: Dict[str, str] = {
    "xy": _DIRECTION_FORWARD,
    "yx": _DIRECTION_REVERSE,
    "ab": _DIRECTION_FORWARD,
    "ba": _DIRECTION_REVERSE,
    "forward": _DIRECTION_FORWARD,
    "reverse": _DIRECTION_REVERSE,
}
"""Accepted spellings of a viewing direction, lower-cased, whitespace- and
arrow-stripped.

`'ab'` / `'ba'` are **positional**, not chain letters: they mean "first chain of
the ordered pair to second" and the reverse, so they stay correct for a pair
whose ids are not `A` and `B`. The real chain ids are accepted too when they are
handed to `resolve_direction`.
"""


def resolve_direction(
    direction: Optional[str],
    chain_x: Optional[str] = None,
    chain_y: Optional[str] = None,
) -> Optional[str]:
    """
    Normalise a user-supplied viewing direction to `'xy'`, `'yx'` or `None`.

    `None` means "no direction selected": show the combination `ipsae.py`
    reports, which is what every default in this module does.

    Args:
        direction: `None`, or any spelling in `_DIRECTION_ALIASES`, or the two
                   chain ids in either order (`'AB'`, `'B->A'`, `'B → A'`).
        chain_x:   First chain of the ordered pair, if the ids should be accepted.
        chain_y:   Second chain of the ordered pair.

    Returns:
        `'xy'`, `'yx'`, or `None`.

    Raises:
        ValueError: On anything else. A typo must fail loudly: silently falling
            back to a direction would mean the reader is told they are looking at
            one measurement while shown the other.

    Example
    -------
    >>> resolve_direction(None) is None
    True
    >>> resolve_direction('BA', 'A', 'B'), resolve_direction('reverse')
    ('yx', 'yx')
    >>> resolve_direction('C -> A', 'C', 'A')
    'xy'
    >>> resolve_direction('sideways')
    Traceback (most recent call last):
        ...
    ValueError: Unknown direction 'sideways'. Use None for the reported combination, 'xy' for first -> second, or 'yx' for second -> first.
    """
    if direction is None:
        return None
    token = str(direction).strip().lower()
    for junk in (" ", "\t", "→", "-", ">", "<", "_"):
        token = token.replace(junk, "")
    if chain_x and chain_y:
        low_x, low_y = str(chain_x).strip().lower(), str(chain_y).strip().lower()
        if token == low_x + low_y:
            return _DIRECTION_FORWARD
        if token == low_y + low_x:
            return _DIRECTION_REVERSE
    if token in _DIRECTION_ALIASES:
        return _DIRECTION_ALIASES[token]
    raise ValueError(
        f"Unknown direction {direction!r}. Use None for the reported "
        f"combination, 'xy' for first -> second, or 'yx' for second -> first."
    )


def describe_direction(
    direction: Optional[str],
    chain_x: str,
    chain_y: str,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
) -> str:
    """
    Human-readable name of a viewing direction, e.g. `'ISG20 (A) → CALM1 (B)'`.

    Args:
        direction: `'xy'`, `'yx'`, or `None` for the reported combination.
        chain_x:   First chain of the ordered pair.
        chain_y:   Second chain of the ordered pair.
        label_x:   Display name for `chain_x`.
        label_y:   Display name for `chain_y`.

    Returns:
        The arrow form, or a description of the default when `direction` is
        `None`.

    Example
    -------
    >>> describe_direction('xy', 'A', 'B', 'ISG20 (A)', 'CALM1 (B)')
    'ISG20 (A) → CALM1 (B)'
    >>> describe_direction('yx', 'A', 'B')
    'Chain B → Chain A'
    >>> describe_direction(None, 'A', 'B')
    'both directions, combined as ipsae.py reports'
    """
    resolved = resolve_direction(direction, chain_x, chain_y)
    name_x = _chain_label(chain_x, label_x)
    name_y = _chain_label(chain_y, label_y)
    if resolved is None:
        return "both directions, combined as ipsae.py reports"
    if resolved == _DIRECTION_FORWARD:
        return f"{name_x} → {name_y}"
    return f"{name_y} → {name_x}"


_SYMMETRIC_NOTE: str = (
    "pDockQ has no second measurement to show. Swapping the chains leaves both "
    "the contact-pair count and the interface residue set unchanged, so it is "
    "symmetric by construction rather than equal by coincidence, and ipsae.py "
    "prints it without a max (ipsae_v4.py:989). It is marked symmetric with no "
    "delta rather than given a zero, because a zero would suggest a difference "
    "had been measured."
)


def format_directional_report(
    rows: Sequence[DirectionalDelta],
    chain_x: str,
    chain_y: str,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
    direction: Optional[str] = None,
    width: int = 78,
) -> str:
    """
    The directional breakdown as printed text: both directions, the gap, the
    reported value, and which direction supplied it (R008).

    The headline number is never displaced. Every row still carries the value
    `ipsae.py` reports, in its own column, and the two directional columns sit
    beside it as the evidence behind it. `direction` changes only which column is
    marked for inspection; no score is recomputed and none can change.

    Args:
        rows:      From `directional_deltas`.
        chain_x:   First chain of the ordered pair.
        chain_y:   Second chain of the ordered pair.
        label_x:   Display name for `chain_x`.
        label_y:   Display name for `chain_y`.
        direction: Viewing direction, per `resolve_direction`. `None` shows the
                   reported combination and marks only the winning direction.
        width:     Wrap width for the prose paragraphs.

    Returns:
        A block of text with no trailing newline, ready to `print`.

    Example
    -------
    The two directional columns flank the reported value rather than replacing
    it, `◄` marks the direction the reported value came from, and a symmetric
    score is labelled instead of diffed:

    >>> rows = [DirectionalDelta('ipsae_d0res', True, 'max', 0.705718,
    ...                          0.555450, 0.705718),
    ...         DirectionalDelta('lis', True, 'mean', 0.600423, 0.608776, 0.592070),
    ...         DirectionalDelta('pdockq', False, 'symmetric', 0.145200)]
    >>> report = format_directional_report(rows, 'A', 'B', 'ISG20 (A)', 'CALM1 (B)')
    >>> print('\\n'.join(report.splitlines()[9:15]))
      Score              ISG20 (A) →      CALM1 (B) →          |Δ|    Reported  Combined
                           CALM1 (B)        ISG20 (A)
      ────────────────────────────────────────────────────────────────────────────────
      ipSAE_d0res             0.5554        ◄  0.7057       0.1503      0.7057  max        [!]
      LIS                     0.6088           0.5921       0.0167      0.6004  mean
      pDockQ               symmetric        symmetric          n/a      0.1452  n/a
    """
    resolved = resolve_direction(direction, chain_x, chain_y)
    name_x = _chain_label(chain_x, label_x)
    name_y = _chain_label(chain_y, label_y)
    head_xy, head_yx = f"{name_x} →", f"{name_y} →"
    col = max(len(head_xy), len(head_yx), len(name_x), len(name_y), 13)
    rule = 2 + 15 + 2 * (col + 2) + 13 + 12 + 10

    out: List[str] = []
    title = "── Directional Breakdown "
    out.append(title + "─" * max(rule - len(title), 3))
    out.append("")
    out.extend(textwrap.wrap(
        "PAE is asymmetric (PAE[i, j] is not PAE[j, i]), so every inter-chain "
        "score is measured twice, once from each chain's frame of reference. "
        "The two are different measurements, not one measurement seen twice, "
        "and the reported score is a combination of them rather than either one "
        "of them. Both are shown here so that a reader can never mistake the "
        "reported number for a property of the complex when it is a property of "
        "one direction.", width))
    out.append("")
    out.append(f"  {'Score':<15s}{head_xy:>{col + 2}s}  {head_yx:>{col + 2}s}"
               f"{'|Δ|':>13s}{'Reported':>12s}  Combined")
    out.append(f"  {'':<15s}{name_y:>{col + 2}s}  {name_x:>{col + 2}s}")
    out.append("  " + "─" * (rule - 2))

    flagged: List[DirectionalDelta] = []
    n_directional = 0
    for row in rows:
        display = SCORE_DISPLAY_NAMES.get(row.name, row.name)
        if not row.directional:
            sym = "symmetric"
            out.append(f"  {display:<15s}{sym:>{col + 2}s}  {sym:>{col + 2}s}"
                       f"{'n/a':>13s}{row.score:>12.4f}  n/a")
            continue
        n_directional += 1
        # `max` has a winning direction; `mean` does not -- both feed the value
        # equally, so marking one would misdescribe how LIS is combined.
        winner = None
        if row.combine == "max":
            winner = (_DIRECTION_FORWARD if row.forward >= row.reverse
                      else _DIRECTION_REVERSE)
        cell_xy = _direction_cell(row.forward, winner == _DIRECTION_FORWARD,
                                  resolved == _DIRECTION_FORWARD)
        cell_yx = _direction_cell(row.reverse, winner == _DIRECTION_REVERSE,
                                  resolved == _DIRECTION_REVERSE)
        flag = "  [!]" if row.flagged else ""
        out.append(
            f"  {display:<15s}{cell_xy:>{col + 2}s}  {cell_yx:>{col + 2}s}"
            f"{row.delta:>13.4f}{row.score:>12.4f}  {row.combine:<9s}{flag}"
        )
        if row.flagged:
            flagged.append(row)

    tolerance = rows[0].tolerance if rows else DIRECTIONAL_DELTA_TOLERANCE
    out.append("")
    out.append("  ◄  the direction the reported value came from")
    if resolved is not None:
        out.append("  »  the direction selected for inspection (DIRECTION)")
    out.append(f"  [!] the two directions differ by more than {tolerance:.4f}")
    out.append("")
    out.extend(textwrap.wrap(_SYMMETRIC_NOTE, width,
                             initial_indent="  ", subsequent_indent="  "))
    out.append("")

    if flagged:
        verb = "differs" if len(flagged) == 1 else "differ"
        out.extend(textwrap.wrap(
            f"{len(flagged)} of the {n_directional} directional scores {verb} "
            f"between the two directions by more than {tolerance:.2f}. For "
            f"that score the reported number is a statement about one chain's "
            f"view of the other, not about the complex as a whole:"
            if len(flagged) == 1 else
            f"{len(flagged)} of the {n_directional} directional scores {verb} "
            f"between the two directions by more than {tolerance:.2f}. For "
            f"those the reported number is a statement about one chain's view "
            f"of the other, not about the complex as a whole:", width))
        for row in flagged:
            out.append("")
            out.extend(_flag_sentence(row, name_x, name_y, width))
    else:
        out.extend(textwrap.wrap(
            f"No directional score differs between the two directions by more "
            f"than {tolerance:.2f}, so on this complex the reported values do "
            f"not depend on which chain is read as the frame of reference. That "
            f"is the expected result when the two chains are copies of one "
            f"protein, but it is measured here, not assumed, and the two "
            f"columns above are where it can be checked.", width))

    out.append("")
    if resolved is None:
        out.extend(textwrap.wrap(
            "Inspecting: both directions (DIRECTION = None, the default). Set "
            "DIRECTION to 'xy' or 'yx' in Section 1 to mark one direction here "
            "and draw only that chain's panel in the per-residue profile figure "
            "below. That is a viewing choice: the Reported column, the summary "
            "table and every traffic light are unaffected by it, because "
            "ipsae.py's combination of the two directions is the definition of "
            "the score.", width))
    else:
        out.extend(textwrap.wrap(
            f"Inspecting: "
            f"{describe_direction(resolved, chain_x, chain_y, label_x, label_y)}"
            f" (DIRECTION = {resolved!r}), marked » above. A viewing choice "
            f"only: the Reported column is unchanged and nothing was "
            f"recomputed. Set DIRECTION = None to go back to both.", width))

    return "\n".join(line.rstrip() for line in out)


def _direction_cell(value: float, is_winner: bool, is_inspected: bool) -> str:
    """One directional value plus its two fixed-width marker slots.

    Example
    -------
    >>> _direction_cell(0.7057, True, False)
    '◄  0.7057'
    >>> _direction_cell(0.5555, False, True)
    ' » 0.5555'
    >>> _direction_cell(0.5555, True, True)
    '◄» 0.5555'
    """
    return f"{'◄' if is_winner else ' '}{'»' if is_inspected else ' '} {value:.4f}"


def _flag_sentence(
    row: DirectionalDelta,
    name_x: str,
    name_y: str,
    width: int = 78,
) -> List[str]:
    """Wrapped prose for one flagged score: which direction won, and what it means."""
    display = SCORE_DISPLAY_NAMES.get(row.name, row.name)
    forward_wins = row.forward >= row.reverse
    high_dir = f"{name_x} → {name_y}" if forward_wins else f"{name_y} → {name_x}"
    low_dir = f"{name_y} → {name_x}" if forward_wins else f"{name_x} → {name_y}"
    low = min(row.forward, row.reverse)
    high = max(row.forward, row.reverse)

    if row.combine == "mean":
        lead = (f"* {display} is reported as {row.score:.4f}, the mean of "
                f"{row.forward:.4f} ({name_x} → {name_y}) and {row.reverse:.4f} "
                f"({name_y} → {name_x}). Neither direction is the score; the "
                f"reported value is one no single measurement produced.")
    else:
        lead = (f"* {display} is reported as {row.score:.4f}, which is the "
                f"{high_dir} measurement. Read {low_dir} the same score is "
                f"{low:.4f}, {row.delta:.4f} lower. The reported value says how "
                f"confidently the model places the first-named chain's residues "
                f"against the second; asked the other way round it is less sure, "
                f"and the max keeps the more optimistic of the two answers.")
    try:
        band_hi = traffic_light(high, row.name)[1]
        band_lo = traffic_light(low, row.name)[1]
        if band_hi != band_lo:
            lead += (f" The gap crosses a band edge: {band_hi} one way, "
                     f"{band_lo} the other.")
    except KeyError:  # pragma: no cover - only for a name outside THRESHOLDS
        pass
    return textwrap.wrap(lead, width, initial_indent="  ",
                         subsequent_indent="    ")


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
# Section 7: one summary object, read by both the table and the prose (R080)
# ---------------------------------------------------------------------------
# The defect this section exists to make impossible: the summary table and the
# plain-text diagnostic used to derive their verdicts separately, so a run could
# print "consistently HIGH confidence across all metrics" underneath a table
# showing an amber badge. Reproduced on FX-008 (AF-0000000211619209, ipSAE_d0res
# 0.6366 amber) and FX-009 (AF-0000000211157965, ipTM_d0chn 0.5769 amber).
#
# The fix is structural rather than a corrected conditional. `summarise_scores`
# builds one `ConfidenceSummary` whose `overall` sentence is *derived from the
# same band counts* the table renders, and both `format_summary_table_html` and
# `format_diagnostic_report` render that one object. There is no second code
# path left that could reach a different conclusion.

BAND_COLOUR_HEX: Dict[Band, str] = {
    "green": "#4CAF50",
    "amber": "#FF9800",
    "red": "#F44336",
}
"""Badge colour per traffic-light band. Status colours, reserved for state: they
are never reused as a series palette. Every badge that carries one also carries
the band's own name, so the verdict is never conveyed by colour alone."""


SCORE_DESCRIPTIONS: Dict[str, str] = {
    "ipsae_d0res": (
        "ipSAE with a per-residue d0, over inter-chain pairs with PAE < 10. "
        "AFDB's ipSAEmax, and one half of its release criterion."
    ),
    "ipsae_d0chn": (
        "Same pairs, but one d0 from the whole chain-pair length, which is the "
        "largest of the three and so the most forgiving variant."
    ),
    "ipsae_d0dom": (
        "Same pairs, d0 from the number of residues on either chain with at "
        "least one sub-cutoff pair. Lies between d0res and d0chn."
    ),
    "iptm_d0chn": (
        "ipSAE_d0chn without the PAE cutoff: every inter-chain cell counts. "
        "A reimplementation from PAE, not AlphaFold's own ipTM."
    ),
    "pdockq": (
        "Sigmoid of mean interface pLDDT x log10(contact pairs). pLDDT is "
        "averaged over the contacting residues of both chains. No PAE."
    ),
    "pdockq2": (
        "Sigmoid of that same mean pLDDT x the mean PAE-to-TM transform "
        "(fixed d0 = 10) over the contact pairs. A product, not a sum."
    ),
    "lis": (
        "Mean of (12 - PAE) / 12 over inter-chain cells below PAE 12: how far "
        "below the cutoff the confident cells sit, not how many there are."
    ),
}
"""One audited line per score, checked against
`specs/homodimer_diagnostic/formula-reference.md` (R081).

Three of the previous seven were wrong and are corrected here:

- **pDockQ** was "Contact count x interface pLDDT". The count is a count of
  contact *pairs*, it enters as its `log10`, the pLDDT is a mean over the union
  of contacting residues in *both* chains, and the product is then passed
  through a sigmoid (formula-reference.md section 5).
- **pDockQ2** was "Contact PAE + interface pLDDT". It is a product, not a sum,
  and the PAE term is a mean of `ptm_func(PAE, 10.0)` over contact pairs rather
  than a mean PAE (formula-reference.md section 6).
- **LIS** was "Density of inter-chain PAE < 12". It is not a density: the mean
  runs over the sub-cutoff cells *only*, so the number of them never enters. A
  block with one cell at PAE 0 scores 1.0 (formula-reference.md section 7).

`ipsae_d0res` also loses "Primary AFDB classifier": the classifier is the
conjunction in `AFDB_JOINT_CRITERION`, of which this score is one side."""


@dataclass(frozen=True)
class SummaryRow:
    """
    One score's row in the Section 7 summary, values and provenance together.

    Attributes:
        name:        `THRESHOLDS` key.
        display:     Display label from `SCORE_DISPLAY_NAMES`.
        value:       The computed score.
        colour:      `'green'` / `'amber'` / `'red'`.
        label:       Band name; AFDB's published wording for `ipsae_d0res`.
        colour_hex:  Badge colour for `colour`.
        provenance:  Provenance of the edge that put the value in this band.
        green:       The green edge, for display.
        amber:       The amber edge, for display.
        description: Audited one-liner from `SCORE_DESCRIPTIONS`.
        source:      Short citation for the threshold, with a page.
        note:        The threshold's caveat, if it carries one.
        independent: Whether this score is one of `AGREEMENT_SCORES`, i.e. one
                     of the values the overall verdict is counted over.
    """

    name: str
    display: str
    value: float
    colour: Band
    label: str
    colour_hex: str
    provenance: Provenance
    green: float
    amber: float
    description: str
    source: str
    note: str
    independent: bool


@dataclass(frozen=True)
class ConfidenceSummary:
    """
    Everything Section 7 says about one complex, derived exactly once.

    `overall` is computed from `n_green` / `n_amber` / `n_red`, which are counted
    from the very `rows` the table renders, so the prose cannot contradict the
    badges (R080). `divergences` are the pairwise readings, each phrased from the
    same band assignments.

    Attributes:
        rows:        Every score supplied, in the order given.
        independent: The subset counted for the overall verdict; by default the
                     five of `AGREEMENT_SCORES`, because the other two ipSAE
                     variants are ordered by a theorem and so cannot supply
                     independent agreement.
        by_name:     `rows` keyed by `THRESHOLDS` key.
        afdb:        The joint release verdict, the headline of the section.
        n_green:     Green count over `independent`.
        n_amber:     Amber count over `independent`.
        n_red:       Red count over `independent`.
        overall:     The OVERALL sentence.
        divergences: Pairwise readings that apply to this complex.
    """

    rows: Tuple[SummaryRow, ...]
    independent: Tuple[SummaryRow, ...]
    by_name: Dict[str, SummaryRow]
    afdb: AFDBHighConfidence
    n_green: int
    n_amber: int
    n_red: int
    overall: str
    divergences: Tuple[str, ...]

    def not_green(self) -> Tuple[SummaryRow, ...]:
        """The independent rows that are not green, worst band first."""
        order = {"red": 0, "amber": 1}
        return tuple(
            sorted(
                (row for row in self.independent if row.colour != "green"),
                key=lambda row: (order[row.colour], row.value),
            )
        )


def _join_names(rows: Sequence[SummaryRow]) -> str:
    """`'a'`, `'a and b'`, `'a, b and c'`."""
    names = [f"{row.display} {row.value:.3f} ({row.label})" for row in rows]
    if len(names) == 1:
        return names[0]
    return ", ".join(names[:-1]) + " and " + names[-1]


def summarise_scores(
    scores: Mapping[str, float],
    independent_names: Optional[Sequence[str]] = None,
) -> ConfidenceSummary:
    """
    Turn a `{THRESHOLDS key: value}` mapping into the one Section 7 summary.

    Every band, colour, count, verdict and sentence Section 7 shows comes from
    the object this returns, so the table and the prose are one derivation with
    two renderings rather than two derivations that happen to usually agree.

    Args:
        scores:            `{THRESHOLDS key: value}`. `ipsae_d0res` and
                           `pdockq2` must be present, since the AFDB criterion
                           is evaluated on them.
        independent_names: Which scores the overall verdict is counted over;
                           defaults to `AGREEMENT_SCORES`.

    Returns:
        A `ConfidenceSummary`.

    Raises:
        KeyError: If a score is not a `THRESHOLDS` key, if `ipsae_d0res` or
            `pdockq2` is missing, or if a requested independent name has no
            value.

    Example
    -------
    >>> s = summarise_scores({'ipsae_d0res': 0.6366, 'ipsae_d0chn': 0.8313,
    ...                       'ipsae_d0dom': 0.7735, 'iptm_d0chn': 0.8183,
    ...                       'pdockq': 0.2418, 'pdockq2': 0.3254,
    ...                       'lis': 0.4783})
    >>> s.n_green, s.n_amber, s.n_red
    (4, 1, 0)
    >>> s.afdb.verdict
    'PASS'

    The overall sentence names the score the table paints amber, so the two
    cannot disagree:

    >>> 'ipSAE_d0res 0.637' in s.overall
    True
    >>> 'consistently' in s.overall
    False
    """
    unknown = [name for name in scores if name not in THRESHOLDS]
    if unknown:
        raise KeyError(
            f"No canonical threshold for: {', '.join(sorted(unknown))}. "
            f"Expected keys from THRESHOLDS: {', '.join(sorted(THRESHOLDS))}."
        )
    for required in ("ipsae_d0res", "pdockq2"):
        if required not in scores:
            raise KeyError(
                f"{required!r} is required: the AFDB joint criterion is "
                "evaluated on ipsae_d0res and pdockq2."
            )

    names = (
        tuple(AGREEMENT_SCORES) if independent_names is None
        else tuple(independent_names)
    )
    missing = [name for name in names if name not in scores]
    if missing:
        raise KeyError(f"No value supplied for: {', '.join(missing)}.")

    rows = []
    for name, value in scores.items():
        threshold = THRESHOLDS[name]
        band = threshold.confidence_band(float(value))
        rows.append(
            SummaryRow(
                name=name,
                display=SCORE_DISPLAY_NAMES.get(name, name),
                value=float(value),
                colour=band.colour,
                label=band.label,
                colour_hex=BAND_COLOUR_HEX[band.colour],
                # The band's own edge provenance, which for a published ladder
                # is per-band rather than the coarse green/amber split.
                provenance=band.provenance,
                green=threshold.green,
                amber=threshold.amber,
                description=SCORE_DESCRIPTIONS.get(name, ""),
                source=band.source or threshold.source,
                note=threshold.note,
                independent=name in names,
            )
        )
    rows = tuple(rows)
    by_name = {row.name: row for row in rows}
    independent = tuple(by_name[name] for name in names)

    n_green = sum(1 for row in independent if row.colour == "green")
    n_amber = sum(1 for row in independent if row.colour == "amber")
    n_red = sum(1 for row in independent if row.colour == "red")
    n_total = len(independent)

    afdb = afdb_high_confidence(
        float(scores["ipsae_d0res"]), float(scores["pdockq2"])
    )

    below = tuple(
        sorted(
            (row for row in independent if row.colour != "green"),
            key=lambda row: ({"red": 0, "amber": 1}[row.colour], row.value),
        )
    )
    if not below:
        overall = (
            f"OVERALL: all {n_total} independent scores are green. "
            "The interface is well resolved, structurally plausible, and the "
            "PAE matrix is confident about the relative chain placement."
        )
    elif n_green == 0:
        overall = (
            f"OVERALL: none of the {n_total} independent scores reaches green "
            f"({_join_names(below)}). Treat any structural conclusion drawn "
            "from this model as unsupported."
        )
    else:
        overall = (
            f"OVERALL: mixed. {n_green} of {n_total} independent scores are "
            f"green; {_join_names(below)} "
            f"{'is' if len(below) == 1 else 'are'} not. The scores disagree, "
            "and the paragraphs below say where."
        )

    divergences = _score_divergences(by_name)

    return ConfidenceSummary(
        rows=rows,
        independent=independent,
        by_name=by_name,
        afdb=afdb,
        n_green=n_green,
        n_amber=n_amber,
        n_red=n_red,
        overall=overall,
        divergences=divergences,
    )


def _score_divergences(by_name: Mapping[str, SummaryRow]) -> Tuple[str, ...]:
    """The pairwise readings, phrased from the same bands the table paints."""
    out: List[str] = []

    ipsae = by_name.get("ipsae_d0res")
    pdockq = by_name.get("pdockq")
    pdockq2 = by_name.get("pdockq2")
    lis = by_name.get("lis")
    iptm = by_name.get("iptm_d0chn")

    if ipsae and pdockq:
        if ipsae.colour == "green" and pdockq.colour != "green":
            out.append(
                f"PAE vs STRUCTURE: AlphaFold is confident about the relative "
                f"chain placement (ipSAE_d0res {ipsae.value:.3f}, "
                f"{ipsae.label}), but the interface itself is thin or poorly "
                f"resolved (pDockQ {pdockq.value:.3f}, {pdockq.label}). This "
                f"happens with a small, tight interface, or when the predicted "
                f"separation is just wide of the {DIST_CUTOFF:.0f} A contact "
                f"cutoff."
            )
        if pdockq.colour == "green" and ipsae.colour != "green":
            out.append(
                f"STRUCTURE vs PAE: there are many contacts between well "
                f"resolved residues (pDockQ {pdockq.value:.3f}, "
                f"{pdockq.label}), but PAE is not confident about how the two "
                f"chains sit relative to one another (ipSAE_d0res "
                f"{ipsae.value:.3f}, {ipsae.label}). Each chain can be "
                f"predicted well while the docking orientation is not. Note "
                f"that Zhu 2023 found pDockQ over-optimistic on exactly this "
                f"case: large, confident, incorrect interfaces."
            )

    if pdockq and pdockq2:
        # Each in units of its own canonical green threshold, so the comparison
        # is between two "fraction of the bar it has to clear" numbers rather
        # than between two incomparable scales.
        pdockq_margin = pdockq.value / pdockq.green
        pdockq2_margin = pdockq2.value / pdockq2.green
        if abs(pdockq_margin - pdockq2_margin) > 0.4:
            lower = pdockq2_margin < pdockq_margin
            out.append(
                f"pDockQ vs pDockQ2 (pDockQ2 "
                f"{'<' if lower else '>'} pDockQ, "
                f"{pdockq2_margin:.2f} vs {pdockq_margin:.2f} times their own "
                f"green thresholds): the two differ by more than they share. "
                f"They use the same contacts and the same mean pLDDT, so the "
                f"gap is entirely the PAE at those contacts, which is "
                f"{'low' if lower else 'high'}. pDockQ2 is the more "
                f"informative of the two and is the one AFDB uses."
            )

    if lis and ipsae and lis.colour == "green" and ipsae.colour != "green":
        out.append(
            f"LIS vs ipSAE: plenty of inter-chain PAE sits below "
            f"{LIS_CUTOFF:.0f} A (LIS {lis.value:.3f}, {lis.label}), but "
            f"tightening the cutoff to {PAE_CUTOFF:.0f} A and applying the "
            f"TM-score transform drops the score (ipSAE_d0res "
            f"{ipsae.value:.3f}, {ipsae.label}). That reads as a broad, "
            f"diffuse contact region rather than a tight, well defined "
            f"interface."
        )

    if iptm and ipsae and iptm.colour != ipsae.colour:
        # The two read the same PAE block and differ in exactly two ways:
        # ipSAE drops every pair at or above the cutoff, which pushes its score
        # up, and it normalises by how many pairs survived rather than by the
        # whole chain pair, which pushes it down. Which of the two wins is the
        # finding, so the sentence has to branch on the sign rather than assume
        # one of them.
        if iptm.value > ipsae.value:
            reading = (
                "Here the normalisation wins: the confident region is small "
                "relative to the length of the chain pair, so the full-length "
                "average that ipTM_d0chn takes is flattering the model, and "
                "ipSAE_d0res is the number to trust."
            )
        else:
            reading = (
                "Here the cutoff wins: there is a genuinely confident patch, "
                "and it was being diluted in ipTM_d0chn by a large surrounding "
                "area of high inter-chain PAE that ipSAE discards. A low "
                "ipTM_d0chn beside a high ipSAE_d0res is the signature of a "
                "well defined interface on chains that are otherwise "
                "uncertain about each other."
            )
        out.append(
            f"ipTM_d0chn vs ipSAE_d0res: the same PAE block read two ways "
            f"disagrees (ipSAE_d0res {ipsae.value:.3f}, {ipsae.label}; "
            f"ipTM_d0chn {iptm.value:.3f}, {iptm.label}). The two differ in "
            f"exactly two things: ipSAE drops every pair with PAE >= "
            f"{PAE_CUTOFF:.0f} A, and it normalises by how many pairs survived "
            f"rather than by the whole chain pair. {reading}"
        )

    return tuple(out)


def format_afdb_verdict(summary: ConfidenceSummary, width: int = 63) -> str:
    """
    The headline verdict: would this model qualify for AFDB high-confidence release?

    This, not the per-score traffic lights, is the classifier AFDB used on ~31
    million candidate complexes. The text states the criterion, the outcome, the
    validation figures, and the two things a reader most often gets wrong: that
    a FAIL is not a claim of "no interaction", and that 0.6 is a conservative
    operating point rather than an optimum.

    Args:
        summary: From `summarise_scores`.
        width:   Wrap width for the prose lines.

    Returns:
        A ready-to-print block of text.

    Example
    -------
    >>> s = summarise_scores({'ipsae_d0res': 0.55, 'pdockq2': 0.64,
    ...                       'iptm_d0chn': 0.8, 'pdockq': 0.4, 'lis': 0.3})
    >>> print(format_afdb_verdict(s).splitlines()[1])
    AFDB HIGH-CONFIDENCE RELEASE CRITERION:  FAIL
    """
    criterion = AFDB_JOINT_CRITERION
    lines = [
        "=" * 65,
        "AFDB HIGH-CONFIDENCE RELEASE CRITERION:  "
        f"{summary.afdb.verdict}",
        "=" * 65,
        "",
        f"  Rule:  ipSAE_d0res >= {criterion.ipsae_d0res_min:.2f}  AND  "
        f"pDockQ2 >= {criterion.pdockq2_min:.2f}",
        f"         ipSAE_d0res {summary.afdb.ipsae_d0res:.4f}  "
        f"{'PASS' if summary.afdb.ipsae_d0res_pass else 'FAIL'}",
        f"         pDockQ2     {summary.afdb.pdockq2:.4f}  "
        f"{'PASS' if summary.afdb.pdockq2_pass else 'FAIL'}",
        "",
    ]
    for paragraph in (
        summary.afdb.reason,
        f"This one conjunction, not the seven traffic lights below, is the "
        f"classifier AlphaFold DB actually applied, to over "
        f"{AFDB_RELEASE_SCALE['candidate_complexes'] / 1e6:.0f} million "
        f"candidate complexes; "
        f"{AFDB_RELEASE_SCALE['homodimers_high_confidence']:,} homodimers "
        f"({AFDB_RELEASE_SCALE['homodimer_high_confidence_fraction']:.1%}) and "
        f"{AFDB_RELEASE_SCALE['heterodimers_high_confidence']:,} heterodimers "
        f"({AFDB_RELEASE_SCALE['heterodimer_high_confidence_fraction']:.1%}) "
        f"passed it. A 9% pass rate is the filter working as designed, not a "
        f"measure of how many predictions are wrong.",
        criterion.note,
        f"A FAIL means the model would not be surfaced as a high-confidence "
        f"AFDB entry. It does not mean the pair does not interact: "
        f"below-threshold dimers are published with their interface scores at "
        f"{AFDB_BELOW_THRESHOLD_FTP}.",
    ):
        lines.append(textwrap.fill(paragraph, width=width,
                                   initial_indent="  ", subsequent_indent="  "))
        lines.append("")
    lines.append(f"  Source: {criterion.source.split(':')[0]}.")
    lines.append("=" * 65)
    return "\n".join(lines)


def format_summary_table_html(
    summary: ConfidenceSummary,
    title: str = "",
    subtitle: str = "",
) -> str:
    """
    Render the summary table, provenance column included.

    Every number, band, colour and citation comes from the `ConfidenceSummary`,
    which is the same object `format_diagnostic_report` renders, so the table
    and the prose are two views of one derivation (R080).

    The provenance column is the point of the rebuild: a coloured badge alone
    gives a reader no way to tell a published operating cutoff from a round
    number somebody picked. `PUBLISHED` / `DERIVED` / `HEURISTIC` are defined in
    `specs/homodimer_diagnostic/threshold-reference.md` and carried on
    `Threshold.green_provenance` / `.amber_provenance` and on each
    `ConfidenceBand`.

    Args:
        summary:  From `summarise_scores`.
        title:    Heading text; omitted when empty. Escaped.
        subtitle: Sub-heading **HTML**, e.g. the two chain descriptions joined
                  by a `<br>`. Interpolated verbatim, so a caller passing a
                  protein name straight from AFDB should escape it first; every
                  string taken from the summary itself is escaped here.

    Returns:
        An HTML string, for `IPython.display.HTML`.

    Example
    -------
    >>> s = summarise_scores({'ipsae_d0res': 0.81, 'pdockq2': 0.64,
    ...                       'iptm_d0chn': 0.8, 'pdockq': 0.4, 'lis': 0.3})
    >>> html = format_summary_table_html(s, title='Demo')
    >>> 'VERY HIGH-CONFIDENCE' in html and 'PUBLISHED' in html
    True
    """
    prov_hex = {"PUBLISHED": "#37474F", "DERIVED": "#607D8B",
                "HEURISTIC": "#B0850F"}

    dagger = ('<span style="font-weight:normal;color:#999;"> '
              '&dagger;</span>')

    # Descriptions and reasons contain '<' ("PAE < 10", "ipSAE_d0res 0.55 <
    # 0.60"), which a browser reads as the start of a tag and silently eats
    # along with everything after it. Escape every interpolated string.
    esc = _html.escape

    rows_html = ""
    for row in summary.rows:
        edges = f"green &ge; {row.green:g}, amber &ge; {row.amber:g}"
        rows_html += (
            '<tr style="border-top:1px solid #eee;">'
            '<td style="font-weight:bold;padding:6px 12px;white-space:nowrap;">'
            f'{esc(row.display)}{"" if row.independent else dagger}'
            '</td>'
            '<td style="padding:6px 12px;text-align:center;font-size:1.15em;'
            f'font-weight:bold;">{row.value:.4f}</td>'
            '<td style="padding:6px 12px;text-align:center;">'
            f'<span style="background:{row.colour_hex};color:white;'
            'padding:3px 10px;border-radius:12px;font-weight:bold;'
            f'font-size:0.85em;white-space:nowrap;">{esc(row.label)}</span></td>'
            '<td style="padding:6px 12px;text-align:center;">'
            f'<span style="color:{prov_hex.get(row.provenance, "#607D8B")};'
            'font-size:0.78em;font-weight:bold;letter-spacing:0.04em;'
            f'white-space:nowrap;">{esc(row.provenance)}</span>'
            '<br><span style="color:#999;font-size:0.72em;white-space:nowrap;">'
            f'{edges}</span></td>'
            '<td style="padding:6px 12px;font-size:0.88em;color:#555;">'
            f'{esc(row.description)}</td>'
            '</tr>'
        )

    heading = f'<h3 style="margin-bottom:4px;">{esc(title)}</h3>' if title else ""
    sub = (
        f'<p style="font-family:sans-serif;color:#555;margin:0 0 10px;">'
        f'{subtitle}</p>' if subtitle else ""
    )
    verdict_hex = "#4CAF50" if summary.afdb.passed else "#F44336"
    banner = (
        '<div style="font-family:sans-serif;margin:0 0 14px;padding:10px 14px;'
        f'border-left:6px solid {verdict_hex};background:#FAFAFA;">'
        '<span style="font-weight:bold;font-size:1.05em;">'
        'AFDB high-confidence release criterion: '
        f'<span style="color:{verdict_hex};">{summary.afdb.verdict}</span>'
        '</span><br>'
        '<span style="color:#555;font-size:0.88em;">'
        f'ipSAE_d0res &ge; {AFDB_IPSAE_D0RES_MIN:.2f} AND pDockQ2 &ge; '
        f'{AFDB_PDOCKQ2_MIN:.2f}: {esc(summary.afdb.reason)}'
        '</span></div>'
    )
    footer = (
        '<p style="font-family:sans-serif;color:#777;font-size:0.8em;'
        'margin:8px 0 0;">'
        'PUBLISHED = the cited paper states or applies this cutoff as an '
        'operating criterion. DERIVED = obtained from a paper\'s own '
        'quantitative statements by one documented inferential step. '
        'HEURISTIC = a judgement call with no literature basis. '
        '&dagger; = not counted in the overall verdict: '
        'ipSAE_d0chn &ge; ipSAE_d0dom &ge; ipSAE_d0res is a theorem, so these '
        'two cannot supply independent agreement.'
        '</p>'
    )
    return (
        heading + sub + banner +
        '<table style="border-collapse:collapse;width:100%;'
        'font-family:sans-serif;">'
        '<thead><tr style="background:#f5f5f5;">'
        '<th style="padding:8px 12px;text-align:left;">Score</th>'
        '<th style="padding:8px 12px;">Value</th>'
        '<th style="padding:8px 12px;">Confidence</th>'
        '<th style="padding:8px 12px;">Threshold provenance</th>'
        '<th style="padding:8px 12px;text-align:left;">What it measures</th>'
        f'</tr></thead><tbody>{rows_html}</tbody></table>' + footer
    )


def format_diagnostic_report(
    summary: ConfidenceSummary,
    width: int = 63,
) -> str:
    """
    The plain-text diagnostic, rendered from the same object as the table.

    The overall sentence is `summary.overall`, which is computed from the band
    counts over `summary.independent` -- the same bands the table paints. There
    is no separate threshold comparison here and no `green / 2` anywhere, so the
    prose cannot label a score differently from the row above it (R080).

    Args:
        summary: From `summarise_scores`.
        width:   Wrap width.

    Returns:
        A ready-to-print block of text.

    Example
    -------
    On FX-008's values the old code printed "consistently HIGH confidence across
    all metrics" while the table showed an amber badge. It now names the amber
    score in the same sentence:

    >>> s = summarise_scores({'ipsae_d0res': 0.6366, 'ipsae_d0chn': 0.8313,
    ...                       'ipsae_d0dom': 0.7735, 'iptm_d0chn': 0.8183,
    ...                       'pdockq': 0.2418, 'pdockq2': 0.3254,
    ...                       'lis': 0.4783})
    >>> report = format_diagnostic_report(s)
    >>> 'consistently HIGH confidence' in report
    False
    >>> 'ipSAE_d0res 0.637 (LOW-CONFIDENCE) is not' in ' '.join(report.split())
    True
    """
    lines = ["=" * 65, "DIAGNOSTIC INTERPRETATION", "=" * 65, ""]
    lines.append(textwrap.fill(summary.overall, width=width))
    for statement in summary.divergences:
        lines.append("")
        lines.append(textwrap.fill(statement, width=width))
    lines.append("")
    lines.append("─" * 65)
    lines.append("Score summary († = not independent, see the table):")
    for row in summary.rows:
        mark = " " if row.independent else "†"
        lines.append(
            f" {mark}{row.display:14s}: {row.value:.4f}  "
            f"[{row.label}]  {row.provenance}"
        )
    lines.append("─" * 65)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
# One `plot_*` function per figure, each returning a `matplotlib` `Figure` and
# never calling `plt.show()`, so a notebook cell is one call plus a title and the
# same figure can be saved head-lessly by a test.
#
# Three deliberate properties of this section:
#
# 1. **Figures are built from `matplotlib.figure.Figure`, not `plt.subplots`.**
#    A pyplot figure is registered in pyplot's global list and the inline backend
#    then draws it at the end of the cell *in addition to* the returned object
#    being rendered, which shows every figure twice. A bare `Figure` is owned by
#    the caller, renders once, and never leaks between cells.
# 2. **Nothing here mutates global state at import time.** The notebook's style
#    block lives in `apply_plot_style()` and is applied only when called, so
#    `import complex_interface_utils` leaves `plt.rcParams` untouched.
# 3. **Chain-pair-generic (D4).** No `'A'` / `'B'` literals, no assumption that
#    `nx == ny`. Every label is derived from the chain ids on the dataclass being
#    plotted, so the heterodimer fixture renders as correctly as the homodimer.
#
# Palette seam (R051)
# -------------------
# `PAE_CMAP` is the single module-level name every PAE figure resolves through,
# at call time rather than at definition time, and every such function also takes
# a `cmap=` override. Changing the palette for the whole notebook is therefore
# one assignment -- `ciu.PAE_CMAP = 'colourblind_safe'` -- and not a rewrite of
# the functions below. `PAE_CMAP_CHOICES` names the alternatives.
#
# Milestone M5, landed
# --------------------
# R013 moved these figures out of the notebook without changing a pixel, and
# left three known defects for M5 to fix. All three are now fixed here:
#
# * **R050, size and aspect.** `figsize=(10, 9)` at `figure.dpi = 150` made the
#   full PAE matrix 1500x1350 px, which overflows a notebook output area and is
#   what put a scrollbar around it. It is now `PAE_MATRIX_FIGSIZE`, 900x810 px.
#   Every PAE and contact image also draws with `aspect='equal'` instead of
#   `'auto'`, so one residue is one square and a 563x69 block reads as 8.16:1
#   rather than being stretched to the axes box. Both the contact map and the
#   score-mask panel size themselves *from* the block's aspect ratio, because an
#   equal-aspect image in a fixed box is either a sliver or a sea of white.
# * **R051, palette.** The default is a sequential green with dark = low PAE =
#   confident. See `PAE_CMAP`.
# * **R052, the score-mask panel.** One shared colour bar across all four axes
#   instead of one stolen from the fourth; a locator strip that says in words and
#   in a diagram which quadrant of the full matrix is being shown; and a hatched
#   slate underlay for "not used by this score" that stays separable from a green
#   ramp whose pale end is nearly white. See `plot_pae_score_masks`.
#
# Ink on a pale ramp
# ------------------
# `RdBu_r` was dark at *both* ends, so white annotation lines and white quadrant
# labels were legible wherever they landed. A sequential ramp is pale at one end,
# so every annotation drawn on top of PAE data now carries a contrasting stroke
# (`_haloed`) rather than relying on the map being dark.


NOTEBOOK_RC_PARAMS: Dict[str, Any] = {
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
}
"""The notebook's own `rcParams` block, moved here verbatim so the notebook does
not have to carry it (R016). Applied by `apply_plot_style()`, never on import."""


# -- colours ----------------------------------------------------------------
# Moved from the notebook's import cell. Named per role, not per chain letter,
# so a chain pair that is not literally A and B still gets the right colours.

COLOUR_A: str = "#009688"
"""Teal. The first chain of the ordered pair."""

COLOUR_B: str = "#FF7043"
"""Coral. The second chain of the ordered pair."""

COLOUR_IF: str = "#FFC107"
"""Amber. Interface residues, in every figure and every 3D view."""

COLOUR_NON_IF: str = "#CCCCCC"
"""Light grey. Non-interface residues in the interface coverage bars."""

COLOUR_NON_IF_HIST: str = "#607D8B"
"""Slate. Non-interface residues in the pLDDT histogram, where the coverage
bars' light grey would be invisible under 55% alpha."""

CHAIN_COLOURS: Tuple[str, str] = (COLOUR_A, COLOUR_B)
"""`(first, second)` chain colours, indexed by position in the ordered pair."""

COLOUR_SIDE_CHAIN_X: str = "#FFD400"
"""Gold. Interface side chains of the **first** chain of the ordered pair, drawn
over that chain's teal cartoon (R071)."""

COLOUR_SIDE_CHAIN_Y: str = "#6495ED"
"""Cornflower. Interface side chains of the **second** chain of the ordered
pair, drawn over that chain's coral cartoon (R071)."""

SIDE_CHAIN_COLOURS: Tuple[str, str] = (COLOUR_SIDE_CHAIN_X, COLOUR_SIDE_CHAIN_Y)
"""`(first, second)` interface side-chain colours, indexed by position in the
ordered pair, parallel to `CHAIN_COLOURS`.

**Why gold and cornflower, and why in that order.** Four colours share one
scene: two cartoons and two sets of side chains. The user proposed cornflower
for the first chain and orange for the second; simulated under deuteranopia
(Machado, Oliveira and Fernandes 2009, severity 1.0) that pairing collapses.
Coral `#FF7043` simulates to a mustard `#BBA83F` and orange `#F5A623` to
`#D4BE29`, CIEDE2000 `dE = 7.3` apart, so the second chain's side chains would
have been invisible against the second chain's own cartoon: the single most
important contrast in the view.

The fix is structural, not a matter of taste. The second chain's cartoon is
already warm, so its side chains have to be cool; the first chain's cartoon is
teal, which desaturates to a neutral grey under both deuteranopia and
protanopia, so its side chains have to be strongly chromatic. Cornflower is
therefore kept but moved to the second chain, and the first chain takes gold.
The worst pair over the four colours plus the viewer's white background is then
`dE = 15.1` (deuteranopia, gold side chains against the coral cartoon of the
*other* chain) and every within-chain contrast is `dE >= 29.9`.

Gold is also within a hair of `COLOUR_IF`, the amber the 2D figures already use
for interface residues, so the first chain's sticks read as "interface" in the
same colour language as the rest of the notebook."""

SCORE_PROFILE_COLOURS: Dict[str, str] = {
    "iptm_d0chn": "#1976D2",
    "ipsae_d0res": "#388E3C",
    "ipsae_d0chn": "#F57C00",
    "ipsae_d0dom": "#7B1FA2",
}
"""Line colour per per-residue score profile, keyed by `THRESHOLDS` key."""

PLDDT_BANDS: Tuple[Tuple[float, str, str], ...] = (
    (90.0, "#1565C0", ">90 (very high)"),
    (70.0, "#42A5F5", "70–90 (confident)"),
    (50.0, "#FFCA28", "50–70 (low)"),
    (float("-inf"), "#EF6C00", "<50 (very low)"),
)
"""AlphaFold's own pLDDT colour ladder, ordered high to low as
`(exclusive lower edge, colour, label)`. A residue takes the first band whose
edge it is strictly above, matching AlphaFold's `> 90` / `> 70` / `> 50` tests.
The en dashes in the labels are numeric ranges, not em dashes, and stay."""

PLDDT_BAND_COLOURS: Tuple[str, ...] = tuple(colour for _, colour, _ in PLDDT_BANDS)
"""Just the four colours of `PLDDT_BANDS`, high band first."""


# -- colormaps --------------------------------------------------------------

PAE_CMAP: str = "Greens_r"
"""Default colormap for every PAE figure. **This is the palette seam (R051).**

Resolved at call time by `resolve_pae_cmap`, so assigning
`complex_interface_utils.PAE_CMAP = 'colourblind_safe'` -- one line, in the
notebook's user-input cell -- retunes every PAE figure at once.

Why `'Greens_r'`
----------------
It is sequential, single-hue, and **dark at PAE 0**: the cells AlphaFold is
confident about are the ones that carry ink, and "no idea" fades toward the page
instead of shouting in red. That orientation is the requirement, not a
preference -- a pale low-PAE end would put the informative half of the figure at
1.0:1 against a white background.

Its lightness decreases monotonically along the whole ramp under normal vision
*and* under simulated deuteranopia and protanopia (Machado-Oliveira-Fernandes
2009, severity 1.0), which is what makes a single-hue ramp safe: the reader is
following lightness, and lightness is the channel that survives colour-vision
deficiency. `'colourblind_safe'` is offered anyway, for readers who would rather
not depend on that.

The trade-off this ramp accepts is at the *pale* end: `Greens_r` reaches
`#f7fcf5`, 1.01:1 against white, so a high-PAE cell has almost no contrast with
the page. For a discrete mark that would be disqualifying. For a contiguous
heatmap it is the domain convention (AlphaFold's own PAE viewers all fade to
near-background at the uncertain end) and it is deliberate: uncertainty should
recede. The figures compensate structurally rather than by truncating the ramp
-- the axes frame stays drawn, so the block's extent is always visible, and the
"not used" underlay is a hatched slate rather than a grey that would collide
with the pale end.
"""

PAE_CMAP_CHOICES: Dict[str, str] = {
    "green": "Greens_r",
    "colourblind_safe": "viridis",
    "rdbu": "RdBu_r",
}
"""Named shorthands accepted anywhere a PAE colormap is taken.

`'green'`
    The default. Sequential, dark = low PAE = confident.
`'colourblind_safe'`
    `viridis`, perceptually uniform and multi-hue, with the largest
    step-to-step separation of the three under both simulated deuteranopia and
    protanopia. Note the name is `viridis`, **not** `viridis_r`: `viridis_r`
    puts bright yellow at PAE 0, which inverts the dark = confident rule this
    module depends on, and is the orientation bug this entry used to carry.
`'rdbu'`
    The pre-M5 default, kept so an older figure can be reproduced. It is
    *diverging*, and PAE is not a diverging quantity: its white midpoint falls
    at an arbitrary PAE (half of `max_pae`), and its lightness is non-monotone,
    so PAE 3 and PAE 28 render at nearly the same lightness -- which under
    deuteranopia or protanopia, where the red and blue arms converge, makes a
    confident cell and a hopeless one genuinely hard to tell apart. Not
    recommended; retained for continuity only.
"""

PAE_UNUSED_COLOUR: str = "#90A4AE"
"""Fill for cells a score does not read, in `plot_pae_score_masks`.

A slate that has to stay separable from *every* step of the PAE ramp, which is
harder than it sounds: a single-hue ramp spans the full lightness range, so a
neutral of any lightness collides with some step of it. This one is chosen by
search to maximise the worst-case OKLab separation against `Greens_r` under
simulated colour-vision deficiency -- minimum dE 14.0 normal, 9.1 deuteranopia,
11.0 protanopia, all above the 8.0 target -- while keeping OKLab chroma at 0.027,
low enough that it still reads as grey rather than as a fifth series colour.

Against the non-default ramps the colour alone is weaker (viridis: 6.0
deuteranopia; RdBu_r: 4.8 protanopia), which is why `PAE_UNUSED_HATCH` exists:
texture is the one channel that does not depend on which palette is selected.
"""

PAE_UNUSED_HATCH: str = "///"
"""Hatch drawn over `PAE_UNUSED_COLOUR`, so "not used" is carried by texture as
well as by colour and survives a palette switch, a monochrome print, and
tritanopia (where the slate's separation from `Greens_r` drops to 6.5)."""

PAE_UNUSED_HATCH_COLOUR: str = "#FFFFFF"
"""Hatch line colour for `PAE_UNUSED_HATCH`."""

DIST_CMAP: str = "viridis_r"
"""Colormap for the contact map, which shows *distance*, not PAE, and so is not
governed by `PAE_CMAP`: dark = close = a tighter contact."""

AGREEMENT_CMAP: str = "RdYlGn_r"
"""Colormap for the score agreement matrix: green = the two scores agree."""

SCORE_DISPLAY_NAMES: Dict[str, str] = {
    "ipsae_d0res": "ipSAE_d0res",
    "ipsae_d0chn": "ipSAE_d0chn",
    "ipsae_d0dom": "ipSAE_d0dom",
    "iptm_d0chn": "ipTM_d0chn",
    "pdockq": "pDockQ",
    "pdockq2": "pDockQ2",
    "lis": "LIS",
}
"""Display label per `THRESHOLDS` key. `iptm_d0chn` is labelled with its full
name rather than a bare "ipTM" because it is a PAE-derived reimplementation and
not AlphaFold's own ipTM, which AFDB does not expose (R004)."""

AGREEMENT_SCORES: Tuple[str, ...] = (
    "ipsae_d0res",
    "iptm_d0chn",
    "pdockq",
    "pdockq2",
    "lis",
)
"""The five genuinely independent values: one ipSAE variant plus the four other
scores. The other two ipSAE variants are excluded because
`d0chn >= d0dom >= d0res` is a theorem, so including them would show three
guaranteed agreements as if they were three confirmations.

Read by `plot_threshold_margins`, by `plot_score_agreement`, and by
`summarise_scores`, which counts its overall verdict over exactly this set."""


# -- figure geometry (R050) -------------------------------------------------
# Every figure here is drawn at `figure.dpi = 150`, so an inch of figure is 150
# rendered pixels and the sizes below are chosen in pixels first. The target is
# that a figure fits a notebook output area -- roughly 900 px of usable width in
# Colab and in a default JupyterLab window -- without the browser having to
# scroll it or shrink it so far that the axis labels stop being readable.
#
# The two block figures cannot take a fixed size, because with `aspect='equal'`
# the *data* dictates the shape of the image: the four project fixtures span
# 1:1, 1.79:1, 8.16:1 and 1:7.11. A fixed box would give the extreme cases a
# sliver of image in a field of white. So they size themselves from the block's
# aspect ratio via `_block_panel_size`, and the score-mask panel additionally
# picks its grid from it via `_mask_panel_grid`.

PAE_MATRIX_FIGSIZE: Tuple[float, float] = (6.0, 5.4)
"""Figure size for `plot_pae_matrix`, in inches: 900x810 px at
`figure.dpi = 150`. The PAE matrix is always square, so this one *can* be fixed.
Was `(10, 9)` -- 1500x1350 px -- which is what made the notebook scroll it."""

PLDDT_FIGSIZE: Tuple[float, float] = (6.0, 5.4)
"""Figure size for `plot_plddt_distribution`, in inches: 900x810 px at
`figure.dpi = 150`.

Was `(14, 5)` with the two panels side by side -- 2100x750 px, which the output
area scales down to about 43%, so each panel arrived roughly 450 px wide and the
per-residue profile had one pixel per two residues. The panels are now stacked,
so each one gets the *whole* column, and the figure is sized to the same
900x810 px budget as `PAE_MATRIX_FIGSIZE`: 900 px is the width a notebook column
renders 1:1, and 810 px is the height above which the output area starts to
scroll (M5)."""

PLDDT_PANEL_RATIOS: Tuple[float, float] = (1.0, 1.25)
"""Height split between the two stacked pLDDT panels, histogram then profile.

The profile is the panel that gains from the full width, and it is also the one
carrying a six-entry legend inside the axes, so it takes the larger share. The
histogram is two overlaid distributions and reads fine in the smaller one."""

CONTACT_MAP_PANEL_IN: float = 3.4
"""Longest side, in inches, of the contact-map image in
`plot_interface_contact_map`. The short side follows from the block's aspect."""

CONTACT_MAP_COVERAGE_IN: float = 4.4
"""Width, in inches, of the interface-coverage panel beside the contact map.
Fixed: its content is one residue axis, whose readable width does not depend on
the contact block's shape."""

MASK_PANEL_MAX_IN: float = 2.0
"""Longest side, in inches, of one score-mask panel in the 2x2 layout."""

MASK_PANEL_STRIP_IN: float = 3.0
"""Longest side, in inches, of one score-mask panel in the 1x4 and 4x1 layouts.

Larger than `MASK_PANEL_MAX_IN` because a single row or column spends the space
in one direction only, so the panels can be longer without the figure growing in
the direction that would put a scrollbar back."""

MASK_STRIP_MIN_WIDTH_IN: float = 5.6
"""Floor on the score-mask figure's width, in inches.

The provenance strip has a fixed amount to say, and a 4x1 layout of very wide
blocks would otherwise produce a figure too narrow to say it in. The strip text
is wrapped to whatever width results."""

MASK_LOCATOR_STRIP_IN: float = 1.25
"""Height, in inches, of the provenance strip above the score-mask panels: the
quadrant locator diagram and the sentence that says which block is drawn."""

MASK_LAYOUT_TALL_RATIO: float = 2.0
"""At or above this rows:cols ratio the four score-mask panels are laid out in a
single row. Four tall-narrow blocks side by side use the space a 2x2 would
waste, and keep the figure short enough not to need scrolling."""

MASK_LAYOUT_WIDE_RATIO: float = 0.5
"""At or below this rows:cols ratio the four panels are stacked in a single
column -- the user's own fallback, and the right answer for a wide-short block
such as the 108x768 fixture."""


def _block_panel_size(
    n_rows: int,
    n_cols: int,
    longest: float,
    min_width: float,
    min_height: float,
) -> Tuple[float, float]:
    """
    Inches for an axes box that shows an `n_rows` x `n_cols` block at true aspect.

    The longer side of the block gets `longest` inches and the shorter side
    follows from the ratio, so the box is the shape of the data. The two minima
    are not aspect corrections -- an `aspect='equal'` image never fills a box
    wider than itself -- they only reserve room for tick labels and an axis
    title, which a 0.2 inch wide panel would have nowhere to put.

    Args:
        n_rows:     Block rows (the vertical extent).
        n_cols:     Block columns (the horizontal extent).
        longest:    Inches for the longer side.
        min_width:  Floor on the returned width.
        min_height: Floor on the returned height.

    Returns:
        `(width_in, height_in)`.

    Example
    -------
    >>> _block_panel_size(172, 172, 3.0, 0.1, 0.1)
    (3.0, 3.0)
    >>> w, h = _block_panel_size(563, 69, 3.0, 0.1, 0.1)
    >>> round(h / w, 2), h
    (8.16, 3.0)
    >>> w, h = _block_panel_size(108, 768, 3.0, 0.1, 0.1)
    >>> round(w / h, 2), w
    (7.11, 3.0)
    >>> _block_panel_size(563, 69, 3.0, 1.0, 1.0)   # floors reserve label room
    (1.0, 3.0)
    """
    ratio = n_rows / n_cols
    if ratio >= 1.0:
        height = longest
        width = longest / ratio
    else:
        width = longest
        height = longest * ratio
    return max(width, min_width), max(height, min_height)


def _mask_panel_grid(n_rows: int, n_cols: int) -> Tuple[int, int]:
    """
    `(nrows, ncols)` of the score-mask grid for a block of this shape.

    A 2x2 is right for a squarish block and wrong for an extreme one: four
    8.16:1 slivers in a 2x2 leave two columns of white, and four 1:7.11 ribbons
    leave two rows of it. So the grid follows the block.

    Example
    -------
    >>> _mask_panel_grid(172, 172), _mask_panel_grid(181, 101)
    ((2, 2), (2, 2))
    >>> _mask_panel_grid(563, 69), _mask_panel_grid(108, 768)
    ((1, 4), (4, 1))
    """
    ratio = n_rows / n_cols
    if ratio >= MASK_LAYOUT_TALL_RATIO:
        return 1, 4
    if ratio <= MASK_LAYOUT_WIDE_RATIO:
        return 4, 1
    return 2, 2


def _haloed(foreground: str = "white", linewidth: float = 2.4) -> List[Any]:
    """
    Path effects that keep an annotation readable on any step of a PAE ramp.

    `RdBu_r` was dark at both ends, so plain white annotations worked. A
    sequential ramp is pale at one end and dark at the other, so anything drawn
    on top of PAE data needs a contrasting outline instead of a single colour
    chosen for one end of the map.
    """
    return [mpatheffects.withStroke(linewidth=linewidth, foreground=foreground)]


def apply_plot_style(overrides: Optional[Mapping[str, Any]] = None) -> None:
    """
    Apply the notebook's seaborn style and `rcParams` to the global pyplot state.

    Called explicitly, never on import: a module that restyled every figure in
    the host process merely by being imported would be a side effect nobody
    asked for, and would make this module unusable from a script that has its own
    style.

    Args:
        overrides: Extra `rcParams` applied after `NOTEBOOK_RC_PARAMS`.

    Note:
        `figure.dpi` is 150 here, so one figure inch is 150 rendered pixels.
        R050 sized every figure against that, rather than lowering it: dropping
        the dpi would have shrunk the text along with the figure, and the point
        was to keep the labels readable while the image fits the output area.
    """
    sns.set_style("white")
    plt.rcParams.update(NOTEBOOK_RC_PARAMS)
    if overrides:
        plt.rcParams.update(dict(overrides))


def resolve_pae_cmap(cmap: Optional[str | Colormap] = None) -> Colormap:
    """
    Resolve a PAE colormap, honouring the module default.

    The lookup order is: the explicit argument, then the module-level `PAE_CMAP`
    read *now* rather than captured at definition time. That late read is the
    whole point of the seam -- it is what lets the notebook's user-input cell
    assign `ciu.PAE_CMAP` once and retune every PAE figure drawn afterwards.

    Args:
        cmap: A `Colormap`, a key of `PAE_CMAP_CHOICES`, a matplotlib colormap
              name, or `None` for the module default.

    Returns:
        A `Colormap`.

    Raises:
        KeyError: If the name is neither a shorthand nor a registered colormap.

    Example
    -------
    >>> resolve_pae_cmap().name
    'Greens_r'
    >>> resolve_pae_cmap('colourblind_safe').name
    'viridis'
    >>> resolve_pae_cmap('rdbu').name
    'RdBu_r'
    >>> resolve_pae_cmap('magma').name
    'magma'
    """
    requested = PAE_CMAP if cmap is None else cmap
    if isinstance(requested, Colormap):
        return requested
    return matplotlib.colormaps[PAE_CMAP_CHOICES.get(requested, requested)]


def _chain_label(chain_id: str, label: "Optional[str | ChainLabel]" = None) -> str:
    """
    `label` if given, else `'Chain <id>'`.

    Accepts a `ChainLabel` as well as a plain string, so a caller can pass the
    record straight through and get its `short` form (R021); `ChainLabel.__str__`
    is `short` for exactly this reason.
    """
    if label is None or label == "":
        return f"Chain {chain_id}"
    return str(label)


def plddt_band(
    value: float,
    bands: Sequence[Tuple[float, str, str]] = PLDDT_BANDS,
) -> Tuple[str, str]:
    """
    AlphaFold's `(colour, label)` for one pLDDT value.

    **The single pLDDT ladder in this module (R075).** Both the 2D figures and
    the 3D View 2 route through here, so the legend and the colours it explains
    cannot disagree. A band is the first one the value is *strictly above*,
    which is how AlphaFold states its own ladder ("Very high (pLDDT > 90)",
    "Confident (90 > pLDDT > 70)", "Low (70 > pLDDT > 50)", "Very low
    (pLDDT < 50)"). The lowest band's edge is `-inf`, so every finite value
    lands in exactly one band and nothing falls through undrawn.

    Args:
        value: A pLDDT score, 0-100.
        bands: `(exclusive lower edge, colour, label)` ordered high to low.

    Returns:
        `(colour, label)` of the band `value` falls in.

    Example
    -------
    >>> plddt_band(95.0)
    ('#1565C0', '>90 (very high)')

    The two documented edges. 90.0 is *not* above 90, so it is confident, not
    very high; 100.0 is above 90 and is drawn rather than falling through:

    >>> plddt_band(90.0)[1]
    '70\u201390 (confident)'
    >>> plddt_band(100.0)[1]
    '>90 (very high)'
    >>> plddt_band(0.0)[1]
    '<50 (very low)'
    """
    for minimum, colour, label in bands:
        if value > minimum:
            return colour, label
    return bands[-1][1], bands[-1][2]


def plddt_band_colour(value: float) -> str:
    """
    AlphaFold's colour for one pLDDT value.

    Args:
        value: A pLDDT score, 0-100.

    Returns:
        A hex colour from `PLDDT_BANDS`.

    Example
    -------
    >>> plddt_band_colour(95.0), plddt_band_colour(90.0)
    ('#1565C0', '#42A5F5')
    >>> plddt_band_colour(60.0), plddt_band_colour(10.0)
    ('#FFCA28', '#EF6C00')
    """
    return plddt_band(value)[0]


def score_masks(
    pair: ChainPairPAE,
    contacts: InterfaceContacts,
    pae_cutoff: float = PAE_CUTOFF,
    lis_cutoff: float = LIS_CUTOFF,
) -> Dict[str, np.ndarray]:
    """
    Which cells of the `x -> y` inter-chain PAE block each score actually reads.

    One boolean mask per score, all on `pair.block_xy`, so the panel figure and
    any printed summary agree by construction instead of rebuilding the cutoff
    logic twice.

    Args:
        pair:       The ordered chain pair's PAE quadrants.
        contacts:   Interface contacts of the **same** ordered pair; supplies the
                    pDockQ2 mask.
        pae_cutoff: ipSAE's inter-chain PAE cutoff, tested strictly.
        lis_cutoff: LIS's PAE cutoff, tested strictly and independent of
                    `pae_cutoff`.

    Returns:
        `{THRESHOLDS key: (nx, ny) bool mask}`, insertion-ordered as the panels
        are drawn: ipTM, ipSAE, LIS, pDockQ2.

    Raises:
        ValueError: If `contacts` describes a different ordered pair, or a
            differently shaped block, than `pair` -- which would silently mask
            the PAE block with somebody else's contacts.

    Example
    -------
    >>> block = np.array([[1.0, 11.0], [13.0, 2.0]], dtype=np.float32)
    >>> pair = ChainPairPAE('A', 'B', block, block.T,
    ...                     np.zeros((2, 2), dtype=np.float32),
    ...                     np.zeros((2, 2), dtype=np.float32))
    >>> contacts = InterfaceContacts('A', 'B', np.full((2, 2), 4.0, dtype=np.float32),
    ...                              np.eye(2, dtype=bool), np.ones(2, dtype=bool),
    ...                              np.ones(2, dtype=bool), 8.0)
    >>> {name: int(mask.sum()) for name, mask in score_masks(pair, contacts).items()}
    {'iptm_d0chn': 4, 'ipsae': 2, 'lis': 3, 'pdockq2': 2}
    """
    _check_same_pair(pair, contacts)
    block = pair.block_xy
    return {
        "iptm_d0chn": np.ones(block.shape, dtype=bool),
        "ipsae": block < pae_cutoff,
        "lis": block < lis_cutoff,
        "pdockq2": contacts.contact_mask.copy(),
    }


def _check_same_pair(pair: ChainPairPAE, contacts: InterfaceContacts) -> None:
    """Raise unless `contacts` and `pair` describe the same ordered chain pair.

    Mixing the two up produces a plausible figure rather than an error, so it is
    checked. The shape test also catches the case where the PAE document and the
    structure disagree on a chain length (`verify_chain_lengths`).
    """
    if (pair.chain_x, pair.chain_y) != (contacts.chain_x, contacts.chain_y):
        raise ValueError(
            f"PAE pair is ({pair.chain_x}, {pair.chain_y}) but contacts are "
            f"({contacts.chain_x}, {contacts.chain_y}); both must be the same "
            f"ordered pair."
        )
    if pair.block_xy.shape != contacts.contact_mask.shape:
        raise ValueError(
            f"PAE block is {pair.block_xy.shape} but the contact mask is "
            f"{contacts.contact_mask.shape}."
        )


def plot_interface_contact_map(
    contacts: InterfaceContacts,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
    cmap: str | Colormap = DIST_CMAP,
    figsize: Optional[Tuple[float, float]] = None,
) -> Figure:
    """
    The interface: a distance-coloured contact map beside per-chain coverage bars.

    Left panel: every contact within the cutoff, coloured by CB-CB distance, rows
    `chain_x` and columns `chain_y`. Non-contacts are `nan` and so render as the
    axes background. Right panel: one bar per residue of each chain, amber where
    that residue touches the partner chain, showing *where along the sequence*
    the interface sits -- one contiguous patch reads very differently from a
    scatter of isolated residues.

    The two coverage tracks share one residue axis so their lengths are
    comparable, which for a heterodimer means the shorter chain's track stops
    part-way across. Its terminus is therefore drawn as an end cap and labelled
    with the residue it ends at, and each track's residue count is in its tick
    label, so a shorter chain reads as shorter rather than as truncated (R022).

    The contact panel is drawn at `aspect='equal'`, so one residue is one square
    and the panel is literally the shape of the contact block (R050). That makes
    the figure's own width depend on the data: a 563x69 block is a narrow strip
    and a 108x768 block is a wide ribbon, and a fixed figure size would give one
    of them a sliver of image in a field of white. `figsize=None` therefore
    computes the size from the block; pass a tuple to override.

    Args:
        contacts: Interface contacts of one ordered chain pair.
        label_x:  Display name for `chain_x`; defaults to `'Chain <id>'`.
        label_y:  Display name for `chain_y`.
        cmap:     Colormap for the distance panel. Not `PAE_CMAP`: this panel
                  shows distance, and the two must stay visually distinct.
        figsize:  Figure size in inches, or `None` to size it from the block's
                  aspect ratio.

    Returns:
        The `Figure`. Nothing is shown or saved; the caller decides.
    """
    name_x = _chain_label(contacts.chain_x, label_x)
    name_y = _chain_label(contacts.chain_y, label_y)
    nx, ny = contacts.contact_mask.shape

    panel_w, panel_h = _block_panel_size(nx, ny, CONTACT_MAP_PANEL_IN,
                                         min_width=1.1, min_height=1.0)
    if figsize is None:
        figsize = (panel_w + 1.45 + CONTACT_MAP_COVERAGE_IN + 0.55,
                   max(panel_h, 2.6) + 1.35)

    fig = Figure(figsize=figsize, layout='constrained')
    # The contact panel's share of the width tracks its own shape; the coverage
    # panel's does not, because its content is one residue axis whose readable
    # width has nothing to do with the contact block.
    axes = fig.subplots(1, 2,
                        width_ratios=[panel_w + 0.95, CONTACT_MAP_COVERAGE_IN])

    ax = axes[0]
    contact_distances = np.where(contacts.contact_mask, contacts.dist_matrix, np.nan)
    im = ax.imshow(contact_distances, aspect='equal', origin='lower',
                   cmap=cmap, vmin=0, vmax=contacts.dist_cutoff)
    # A wide-short block leaves its axes box mostly empty, and a colour bar that
    # spans the box rather than the image reads as a scale for whitespace. Both
    # the map and the bar are therefore pinned to the top of the row and the bar
    # is shrunk to roughly the image's own height.
    ax.set_anchor('N')
    true_h = CONTACT_MAP_PANEL_IN * min(1.0, nx / ny)
    fig.colorbar(im, ax=ax, label='CB-CB distance (Å)', fraction=0.05, pad=0.03,
                 shrink=min(1.0, max(0.3, true_h / max(panel_h, 2.6))),
                 anchor=(0.0, 1.0), panchor=False)
    ax.set_xlabel(f'{name_y} residue index', fontsize=10)
    ax.set_ylabel(f'{name_x} residue index', fontsize=10)
    # The panel is as narrow as the block is, so the title is wrapped to the
    # width it actually has rather than being allowed to run off the figure.
    ax.set_title(
        'Interface Contact Map\n' + textwrap.fill(
            f'(contacts ≤ {contacts.dist_cutoff:.0f} Å, coloured by distance; '
            f'{nx}×{ny} block, drawn to scale)',
            width=max(24, int((panel_w + 1.0) * 13))),
        fontsize=9)

    ax2 = axes[1]
    bar_height = 0.35
    tracks = ((contacts.mask_x, 0.6, name_x), (contacts.mask_y, 0.1, name_y))
    for mask, bottom, _name in tracks:
        ax2.bar(np.arange(len(mask)), bar_height, bottom=bottom,
                color=[COLOUR_IF if is_if else COLOUR_NON_IF for is_if in mask],
                width=1.0, linewidth=0)

    # Both tracks share one residue axis, so their lengths are directly
    # comparable -- which is the point of drawing them together, and which means
    # the shorter chain necessarily stops part-way across. Say where it stops
    # (R022): otherwise the empty remainder reads as a truncated bar, i.e. as
    # missing data, rather than as the end of a shorter protein.
    n_longest = max(nx, ny)
    for mask, bottom, name in tracks:
        end = len(mask)
        middle = bottom + bar_height / 2.0
        ax2.plot([end - 0.5, end - 0.5], [bottom, bottom + bar_height],
                 color='#333333', linewidth=1.4, solid_capstyle='butt',
                 zorder=3)
        if end == n_longest:
            continue
        # Annotate into whichever side of the end cap has room, so a chain that
        # ends close to the right edge does not write over the legend.
        near_right = end > 0.75 * n_longest
        ax2.annotate(
            f'{name} ends at residue {end}',
            xy=(end - 0.5, middle),
            xytext=(-6 if near_right else 6, 0), textcoords='offset points',
            ha='right' if near_right else 'left', va='center',
            fontsize=8, color='#333333', zorder=4)

    ax2.set_xlim(-0.5, n_longest - 0.5)
    ax2.set_ylim(0, 1.1)
    ax2.set_yticks([0.275, 0.775])
    ax2.set_yticklabels([f'{name_y}\n{ny} res', f'{name_x}\n{nx} res'])
    ax2.set_xlabel('Residue index (both chains on one scale)')
    ax2.set_title('Interface Coverage\n'
                  '(amber = at interface, grey = non-interface)', fontsize=9)
    ax2.legend(
        handles=[mpatches.Patch(color=COLOUR_IF, label='Interface'),
                 mpatches.Patch(color=COLOUR_NON_IF, label='Non-interface')],
        loc='upper right', fontsize=9)

    return fig


def plot_pae_matrix(
    pae: PAEMatrix,
    chain_x: Optional[str] = None,
    chain_y: Optional[str] = None,
    accession: str = "",
    labels: "Optional[Mapping[str, str | ChainLabel]]" = None,
    cmap: Optional[str | Colormap] = None,
    figsize: Tuple[float, float] = PAE_MATRIX_FIGSIZE,
) -> Figure:
    """
    The full PAE matrix, with the chain boundary and the four quadrants labelled.

    The point of showing the whole matrix rather than the inter-chain block alone
    is that the two intra-chain quadrants are the control: a model can be
    confident about each chain in isolation and have no idea how they sit
    together, and that reads instantly as two dark diagonal blocks with a pale
    off-diagonal.

    Args:
        pae:       The parsed PAE document.
        chain_x:   Chain to label as the first of the pair; defaults to the first
                   chain in matrix layout order.
        chain_y:   The second; defaults to the second chain in layout order.
        accession: Shown in the title when given.
        labels:    `{chain_id: display name}` for **every** chain in the matrix,
                   not just the pair, because the axis names every block. Pass
                   `ChainIdentity.labels` to get real protein names; omit it and
                   each chain is called `'Chain <id>'` (R021).
        cmap:      Colormap override; `None` uses `PAE_CMAP` (the palette seam).
        figsize:   Figure size in inches. The default is `PAE_MATRIX_FIGSIZE`,
                   900x810 px at `figure.dpi = 150`, sized to fit a notebook
                   output area whole (R050).

    Returns:
        The `Figure`.

    Raises:
        ValueError: If the document has fewer than two chains, or if `chain_x`
            and `chain_y` are the same chain.
        KeyError:   If a named chain is not in the document.

    Note:
        `aspect='equal'`, so the matrix renders square -- which it always is,
        being `(nx + ny)` on both sides -- and each quadrant renders at its true
        proportions. Under `'auto'`, which is what this drew before R050, the
        image was stretched to whatever shape the axes box happened to be, and
        a heterodimer's quadrants came out the wrong shape relative to each
        other.
    """
    ids = pae.chain_ids
    if len(ids) < 2:
        raise ValueError(
            f"A PAE matrix needs at least two chains to have quadrants; got {ids}."
        )
    chain_x = ids[0] if chain_x is None else chain_x
    chain_y = ids[1] if chain_y is None else chain_y
    if chain_x == chain_y:
        raise ValueError(f"Need two distinct chains, got {chain_x!r} twice.")

    label_of = {} if labels is None else dict(labels)
    name = {cid: _chain_label(cid, label_of.get(cid)) for cid in ids}

    fig = Figure(figsize=figsize, layout='constrained')
    ax = fig.subplots()

    im = ax.imshow(pae.matrix, aspect='equal', origin='upper',
                   cmap=resolve_pae_cmap(cmap), vmin=0, vmax=pae.max_pae)
    fig.colorbar(im, ax=ax, label='PAE (Å): lower = more confident',
                 fraction=0.046, pad=0.03)

    # One dashed line per internal chain boundary. Two chains give the notebook's
    # single pair of lines at nx - 0.5; more chains give one pair each. The line
    # is dark with a white halo rather than plain white, because the sequential
    # default is pale at the high-PAE end and a white line on a near-white
    # inter-chain block is invisible.
    offset = 0
    for span in pae.spans[:-1]:
        offset += span.length
        for draw in (ax.axhline, ax.axvline):
            draw(offset - 0.5, color='#263238', linewidth=1.6, linestyle='--',
                 path_effects=_haloed('white', 3.6))

    def _centre(chain_id: str) -> float:
        span_slice = pae.chain_slice(chain_id)
        return (span_slice.start + span_slice.stop) / 2.0

    cx, cy = _centre(chain_x), _centre(chain_y)
    name_x, name_y = name[chain_x], name[chain_y]
    for col, row, text in ((cx, cx, f'Intra\n{name_x}'),
                           (cy, cx, f'Inter\n{name_x}→{name_y}'),
                           (cx, cy, f'Inter\n{name_y}→{name_x}'),
                           (cy, cy, f'Intra\n{name_y}')):
        ax.text(col, row, text, ha='center', va='center',
                color='#212121', fontsize=9, fontweight='bold',
                path_effects=_haloed('white', 2.8))

    if len(ids) == 2:
        first, second = ids
        n_first = pae.chain_length(first)
        ax.set_xlabel(f'Residue index ({name[first]}: 0 to {n_first - 1}, '
                      f'{name[second]}: {n_first} to end)', fontsize=10)
    else:
        ax.set_xlabel('Residue index (' + ', then '.join(name[c] for c in ids) + ')',
                      fontsize=10)
    ax.set_ylabel('Residue index', fontsize=10)
    head = f'Full PAE Matrix: {accession}' if accession else 'Full PAE Matrix'
    ax.set_title(f'{head}\n'
                 f'(dashed line = chain boundary between {name_x} and {name_y})',
                 fontsize=11)

    return fig


def _draw_quadrant_locator(
    ax: Any,
    n_rows: int,
    n_cols: int,
    name_x: str,
    name_y: str,
    cmap: Colormap,
) -> None:
    """
    Draw a thumbnail of the full PAE matrix with the plotted quadrant filled in.

    The axis labels on the panels already say which chain is on which axis, but
    they do not say *where in the full matrix* the block came from, and a reader
    who has just looked at the whole matrix has to take that on trust. This says
    it in a picture: the same four blocks in the same arrangement and the same
    proportions as the figure above, with one of them coloured.
    """
    total = n_rows + n_cols
    split = n_rows / total

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(1.0, 0.0)          # origin upper, matching plot_pae_matrix
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    blocks = (
        # (x0, y0, w, h, highlighted)
        (0.0, 0.0, split, split, False),                    # intra x
        (split, 0.0, 1.0 - split, split, True),             # inter x -> y
        (0.0, split, split, 1.0 - split, False),            # inter y -> x
        (split, split, 1.0 - split, 1.0 - split, False),    # intra y
    )
    for x0, y0, width, height, hot in blocks:
        ax.add_patch(mpatches.Rectangle(
            (x0, y0), width, height,
            facecolor=cmap(0.30) if hot else '#ECEFF1',
            edgecolor='#37474F' if hot else '#B0BEC5',
            linewidth=1.6 if hot else 0.8, zorder=2 if hot else 1))

    ax.set_title('Full PAE matrix', fontsize=8)


def plot_pae_score_masks(
    pair: ChainPairPAE,
    contacts: InterfaceContacts,
    max_pae: float,
    pae_cutoff: float = PAE_CUTOFF,
    lis_cutoff: float = LIS_CUTOFF,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
    cmap: Optional[str | Colormap] = None,
    figsize: Optional[Tuple[float, float]] = None,
) -> Figure:
    """
    Four views of the same inter-chain PAE block: what each score actually reads.

    Every panel shows `pair.block_xy` -- the inter-chain quadrant, rows
    `chain_x` and columns `chain_y` -- and they differ only in which cells are
    left coloured. This is the figure that explains why four scores computed
    from one matrix can disagree: they are not weighting the same evidence
    differently, they are reading different subsets of it.

    Three things about how it is drawn (R052):

    * **Provenance is stated, not implied.** A strip above the panels carries a
      thumbnail of the full matrix with the plotted quadrant filled in, and a
      sentence naming the slice (`pae_matrix[:nx, nx:nx+ny]`, the upper-right
      block of the previous figure). The axis labels alone leave the reader to
      infer which quadrant this is.
    * **One colour bar for all four panels.** It is attached to the whole set of
      axes, so every panel gives up the same width. Attaching it to the last
      axis alone -- which is what this did before R052 -- made the fourth panel
      about 20% narrower than the other three, so the four blocks were not
      comparable even though the whole point is comparing them.
    * **Cells a score ignores are hatched slate, not plain grey.** Against the
      old diverging map a flat grey was distinct enough. Against a sequential
      ramp whose high-PAE end is nearly white it is not, so the underlay carries
      a texture as well as a colour; see `PAE_UNUSED_COLOUR`.

    The grid follows the block's shape rather than always being 2x2, and the
    figure sizes itself from it, because at `aspect='equal'` a 8.16:1 block in a
    2x2 grid is four slivers in a field of white. See `_mask_panel_grid`.

    Args:
        pair:       The ordered chain pair's PAE quadrants.
        contacts:   Interface contacts of the **same** ordered pair.
        max_pae:    `vmax`; use `PAEMatrix.max_pae` so every PAE figure of one
                    model shares a scale.
        pae_cutoff: ipSAE's cutoff, named in the panel title rather than
                    hard-coded into it.
        lis_cutoff: LIS's cutoff, likewise. So is pDockQ2's distance cutoff,
                    which is read off `contacts`.
        label_x:    Display name for `chain_x`.
        label_y:    Display name for `chain_y`.
        cmap:       Colormap override; `None` uses `PAE_CMAP`.
        figsize:    Figure size in inches, or `None` to size it from the block.

    Returns:
        The `Figure`. Per-panel cell counts are in the panel titles; the same
        masks are available from `score_masks` if a caller wants the numbers.

    Raises:
        ValueError: If `contacts` and `pair` are not the same ordered pair.
    """
    name_x = _chain_label(contacts.chain_x, label_x)
    name_y = _chain_label(contacts.chain_y, label_y)
    masks = score_masks(pair, contacts, pae_cutoff=pae_cutoff, lis_cutoff=lis_cutoff)
    # Score name, then the rule that selects its cells -- read off the arguments
    # in force, never hard-coded, so a caller who changes a cutoff sees the
    # change in the figure rather than being quietly contradicted by it.
    titles = {
        'iptm_d0chn': ('ipTM_d0chn', 'all inter-chain, no cutoff'),
        'ipsae': ('ipSAE', f'PAE < {pae_cutoff:.0f} Å'),
        'lis': ('LIS', f'PAE < {lis_cutoff:.0f} Å'),
        'pdockq2': ('pDockQ2', f'CB-CB ≤ {contacts.dist_cutoff:.0f} Å contacts'),
    }
    block = pair.block_xy
    n_rows, n_cols = block.shape
    n_cells_total = block.size

    nrows, ncols = _mask_panel_grid(n_rows, n_cols)
    longest = MASK_PANEL_MAX_IN if (nrows, ncols) == (2, 2) else MASK_PANEL_STRIP_IN
    panel_w, panel_h = _block_panel_size(n_rows, n_cols, longest,
                                         min_width=0.9, min_height=0.5)
    fig_w = max(ncols * (panel_w + 0.78) + 1.15, MASK_STRIP_MIN_WIDTH_IN)

    # The provenance sentence is wrapped to the width the figure turned out to
    # have, and the strip is then made tall enough for however many lines that
    # took. Sizing the strip first and hoping the text fits is what puts the
    # last line of it through the first panel's title.
    text_share = 4.6
    text_width = int((fig_w * text_share / (1.0 + text_share) - 0.25) * 14)
    provenance = textwrap.fill(
        f'All four panels below show the SAME block: the inter-chain '
        f'quadrant {name_x} → {name_y}, rows = {name_x} ({n_rows} residues), '
        f'columns = {name_y} ({n_cols} residues). It is '
        f'pae_matrix[:{n_rows}, {n_rows}:{n_rows + n_cols}] (the '
        f'upper-right block of the full PAE matrix in the previous figure, '
        f'shaded at left) and not the whole matrix. The panels differ only '
        f'in which of its {n_cells_total:,} cells each score reads; hatched '
        f'slate = cells that score ignores.',
        width=max(38, text_width))
    strip_in = max(MASK_LOCATOR_STRIP_IN,
                   0.34 + (provenance.count(chr(10)) + 1) * 0.165)

    if figsize is None:
        figsize = (fig_w, strip_in + nrows * (panel_h + 0.86) + 0.25)
    strip_in = min(strip_in, figsize[1] * 0.6)

    fig = Figure(figsize=figsize, layout='constrained')
    strip, grid = fig.subfigures(
        2, 1, height_ratios=[strip_in, max(figsize[1] - strip_in, 1.0)])

    resolved = resolve_pae_cmap(cmap)
    # Unused cells are left as NaN and rendered transparent, so what shows
    # through is the hatched patch drawn underneath -- one patch per axes rather
    # than a second full-size image, which also keeps the hatch in display space
    # so it stays the same weight whatever the block's pixel dimensions are.
    painted = resolved.with_extremes(bad=(0.0, 0.0, 0.0, 0.0))

    loc_ax, text_ax = strip.subplots(1, 2, width_ratios=[1.0, text_share])
    _draw_quadrant_locator(loc_ax, n_rows, n_cols, name_x, name_y, resolved)
    text_ax.axis('off')
    text_ax.text(0.0, 0.97, provenance, transform=text_ax.transAxes,
                 ha='left', va='top', fontsize=9, linespacing=1.3,
                 color='#263238')

    # A gridspec with an explicit colour-bar column, rather than a colour bar
    # that steals space from whichever axes it happens to overlap. This is the
    # R052 fix: the four panel columns are uniform by construction, so the four
    # blocks are the same size and are therefore comparable, which is the entire
    # point of drawing them together.
    spec = grid.add_gridspec(nrows, ncols + 1,
                             width_ratios=[1.0] * ncols + [0.08 * ncols])
    axes = np.array([grid.add_subplot(spec[r, c])
                     for r in range(nrows) for c in range(ncols)])
    cax = grid.add_subplot(spec[:, -1])

    im = None
    for index, (ax, (key, mask)) in enumerate(zip(axes, masks.items())):
        row, col = divmod(index, ncols)
        ax.add_patch(mpatches.Rectangle(
            (0.0, 0.0), 1.0, 1.0, transform=ax.transAxes, zorder=0,
            facecolor=PAE_UNUSED_COLOUR, edgecolor=PAE_UNUSED_HATCH_COLOUR,
            hatch=PAE_UNUSED_HATCH, linewidth=0.0))

        display_pae = block.astype(float).copy()
        display_pae[~mask] = np.nan
        im = ax.imshow(display_pae, aspect='equal', origin='upper',
                       cmap=painted, vmin=0, vmax=max_pae, zorder=1)

        n_cells = int(mask.sum())
        frac = 100.0 * n_cells / n_cells_total if n_cells_total else 0.0
        score_name, rule = titles[key]
        ax.set_title(f'{score_name}\n{rule}\n{n_cells:,} cells ({frac:.1f}%)',
                     fontsize=8.5)
        # Axis names on the edge panels only: every panel shows the same block,
        # so repeating both names four times spends width that a tall-sliver
        # block does not have. Tick labels stay on every panel.
        if row == nrows - 1:
            ax.set_xlabel(f'{name_y} residue', fontsize=8)
        if col == 0:
            ax.set_ylabel(f'{name_x} residue', fontsize=8)
        # A panel can be under an inch wide when the block is a tall sliver, so
        # the tick count is derived from the panel size rather than left to the
        # default, which would overlap its own labels.
        ax.xaxis.set_major_locator(MaxNLocator(nbins=max(2, int(panel_w * 2.2)),
                                               integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=max(2, int(panel_h * 2.2)),
                                               integer=True))
        ax.tick_params(labelsize=7)

    # One colour bar across every panel, so all four give up the same width.
    bar = grid.colorbar(im, cax=cax)
    bar.set_label('PAE (Å): lower = more confident', fontsize=9)
    bar.ax.tick_params(labelsize=8)
    fig.suptitle(textwrap.fill(f'{name_x} → {name_y} inter-chain PAE block: '
                               'cells used by each score',
                               width=max(30, int(fig_w * 11))),
                 fontsize=12)
    return fig


PROFILE_SERIES: Tuple[str, ...] = (
    "iptm_d0chn",
    "ipsae_d0res",
    "ipsae_d0chn",
    "ipsae_d0dom",
)
"""The four values that have a genuine per-residue decomposition (R062).

Each of these is defined residue by residue -- row `i` of a PAE block is
residue `i`'s own measurement against the whole partner chain -- and the value
`ipsae.py` reports for the direction is `values.max()`, i.e. literally one
residue's number. That is what makes a profile meaningful and an argmax marker
worth drawing.

The other three scores are **not** omitted by oversight:

- `pdockq` pools over the interface residue set (one mean pLDDT, one contact
  count) and has no per-residue term at all.
- `lis` pools over every sub-cutoff inter-chain *pair* in the block; a pair is
  not a residue, and nothing in the formula is indexed by row.
- `pdockq2` pools over contact *pairs* too. `PDockQ2Direction.mean_ptm_by_residue`
  does give a per-residue projection of one of its two ingredients, but the
  reported pDockQ2 is a sigmoid of `mean_plddt * mean_ptm` over all pairs, not
  the max over that array, so plotting it here would put a series on an
  argmax-annotated figure whose argmax means nothing. It belongs on the
  structure instead, where the contact set it lives on is visible (R074).
"""

_PROFILE_STYLE: Dict[str, Dict[str, Any]] = {
    "iptm_d0chn": {"linewidth": 1.5, "linestyle": "-"},
    "ipsae_d0res": {"linewidth": 1.5, "linestyle": "-"},
    "ipsae_d0chn": {"linewidth": 1.2, "linestyle": "--"},
    "ipsae_d0dom": {"linewidth": 1.2, "linestyle": ":"},
}

PEAK_MARKER_COLOUR: str = "#37474F"
"""Slate. The legend proxy for the peak-residue star; the stars themselves take
their series colour."""


def _residue_name(index: int, res_ids: Optional[np.ndarray]) -> str:
    """`'res 137'` when residue numbering is known, `'index 136'` when it is not.

    Example
    -------
    >>> _residue_name(2, np.array([11, 12, 13]))
    'res 13'
    >>> _residue_name(2, None)
    'index 2'
    """
    if res_ids is None:
        return f"index {index}"
    return f"res {int(res_ids[index])}"


def _annotate_peaks(
    ax: Any,
    profiles: Mapping[str, ResidueProfile],
    n_res: int,
    res_ids: Optional[np.ndarray],
) -> None:
    """
    Mark each series' argmax residue and label it with residue number and value.

    The reported score of a direction *is* the argmax residue's value, so this is
    the single most consequential point on the panel and it was previously
    invisible (R062). Series peaking on the same residue -- the usual case, since
    the three ipSAE variants differ only in `d0` -- share one label box rather
    than stacking four overlapping annotations.
    """
    groups: Dict[int, List[Tuple[str, float]]] = {}
    for key in PROFILE_SERIES:
        profile = profiles[key]
        index = profile.argmax_index
        if index < 0:
            continue
        value = float(profile.values[index])
        ax.plot([index], [value], marker="*", markersize=12, linestyle="none",
                color=SCORE_PROFILE_COLOURS[key], markeredgecolor="white",
                markeredgewidth=0.8, zorder=6)
        groups.setdefault(index, []).append((key, value))

    # Tallest peak first, so a lower neighbour stacks below an already-placed
    # box rather than through it.
    placed: List[Tuple[int, float, float, int]] = []
    for index, entries in sorted(groups.items(),
                                 key=lambda kv: -max(v for _, v in kv[1])):
        top = max(value for _, value in entries)
        lines = [f"peak · {_residue_name(index, res_ids)}"]
        lines += [f"{SCORE_DISPLAY_NAMES[key]}  {value:.4f}"
                  for key, value in entries]
        # Keep the box inside the axes: left of the star on the right-hand half
        # of the chain, right of it otherwise; below it when the peak sits high.
        on_right_half = index > n_res / 2
        dx = -16 if on_right_half else 16
        ha = "right" if on_right_half else "left"
        high = top > 0.6
        sign = -1 if high else 1
        va = "top" if high else "bottom"
        dy = sign * 10.0
        height = 11.0 * len(lines) + 8.0
        near = max(0.15 * n_res, 5.0)
        for prev_index, prev_dy, prev_height, prev_sign in placed:
            if prev_sign != sign or abs(prev_index - index) > near:
                continue
            candidate = prev_dy + sign * (prev_height + 6.0)
            if abs(candidate) > abs(dy):
                dy = candidate
        placed.append((index, dy, height, sign))
        ax.annotate(
            "\n".join(lines),
            xy=(index, top),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha, va=va, fontsize=8, family="monospace", zorder=7,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                      edgecolor="#B0BEC5", alpha=0.92),
            arrowprops=dict(arrowstyle="-", color="#B0BEC5", linewidth=0.8),
        )


def plot_residue_score_profiles(
    iptm: DirectionalPair,
    ipsae: IPSAEResult,
    contacts: InterfaceContacts,
    plddt_x: np.ndarray,
    plddt_y: np.ndarray,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
    figsize: Tuple[float, float] = (14.0, 10.0),
    direction: Optional[str] = None,
    res_ids_x: Optional[np.ndarray] = None,
    res_ids_y: Optional[np.ndarray] = None,
) -> Figure:
    """
    Per-residue score profiles, one panel per direction of the chain pair.

    Every score reported for a complex is one residue's number -- the maximum
    over the profile -- so this figure is where a headline value stops being a
    verdict and becomes a location: which residues carry the interface, whether
    the three ipSAE variants rank the same residue highest, and whether the peak
    sits on well-predicted backbone or on a low-pLDDT loop. The peak of each
    series is marked with a star and labelled with its residue number and value,
    because that number is the score (R062).

    The top panel is the `x -> y` direction (rows of `block_xy`, so residues of
    `chain_x`); the bottom is `y -> x`. They are genuinely two measurements, not
    one seen twice, because PAE is asymmetric.

    Only the four values in `PROFILE_SERIES` are drawn. pDockQ, pDockQ2 and LIS
    are pooled statistics over the interface residue set or over contact pairs
    and have no per-residue value to plot; see `PROFILE_SERIES` for the full
    reasoning, including why `mean_ptm_by_residue` is left to the 3D views.

    Args:
        iptm:      `compute_iptm_d0chn` result for the ordered pair.
        ipsae:     `compute_ipsae` result for the same ordered pair.
        contacts:  Interface contacts of the same ordered pair; supplies the
                   shaded interface regions.
        plddt_x:   `(nx,)` pLDDT for `chain_x`.
        plddt_y:   `(ny,)` pLDDT for `chain_y`.
        label_x:   Display name for `chain_x`.
        label_y:   Display name for `chain_y`.
        figsize:   Figure size in inches for the two-panel form. A single-panel
                   form takes half the height.
        direction: `None` (default) draws both directions; `'xy'` or `'yx'`, per
                   `resolve_direction`, draws only that one. A viewing choice --
                   it selects a panel, never a score (R008).
        res_ids_x: Optional `(nx,)` residue numbers for `chain_x`, used to label
                   the peak with the residue's real number rather than its
                   positional index.
        res_ids_y: Optional `(ny,)` residue numbers for `chain_y`.

    Returns:
        The `Figure`.

    Raises:
        ValueError: If the three results do not describe the same ordered chain
            pair, if a pLDDT array's length does not match its chain, or if
            `direction` is not a recognised direction.
    """
    pair_ids = (contacts.chain_x, contacts.chain_y)
    for name, result in (('ipTM', iptm), ('ipSAE', ipsae)):
        if (result.chain_x, result.chain_y) != pair_ids:
            raise ValueError(
                f"{name} describes pair ({result.chain_x}, {result.chain_y}) but "
                f"the contacts describe {pair_ids}; all three must agree."
            )

    resolved = resolve_direction(direction, contacts.chain_x, contacts.chain_y)
    name_x = _chain_label(contacts.chain_x, label_x)
    name_y = _chain_label(contacts.chain_y, label_y)

    all_panels = (
        (_DIRECTION_FORWARD, name_x, name_y, 'forward',
         contacts.mask_x, np.asarray(plddt_x), res_ids_x),
        (_DIRECTION_REVERSE, name_y, name_x, 'reverse',
         contacts.mask_y, np.asarray(plddt_y), res_ids_y),
    )
    panels = [panel for panel in all_panels
              if resolved is None or panel[0] == resolved]

    height = figsize[1] if len(panels) == 2 else figsize[1] / 2.0
    fig = Figure(figsize=(figsize[0], height))
    axes = np.atleast_1d(fig.subplots(len(panels), 1, sharex=False))

    for ax, (_, row_label, col_label, attr, if_mask, plddt_chain, res_ids) in zip(axes, panels):
        profiles = {
            'iptm_d0chn': getattr(iptm, attr),
            'ipsae_d0res': getattr(ipsae.d0res, attr),
            'ipsae_d0chn': getattr(ipsae.d0chn, attr),
            'ipsae_d0dom': getattr(ipsae.d0dom, attr),
        }
        n_res = profiles['iptm_d0chn'].values.shape[0]
        if plddt_chain.shape[0] != n_res:
            raise ValueError(
                f"{row_label} has {n_res} residues in its score profile but "
                f"{plddt_chain.shape[0]} pLDDT values."
            )
        if res_ids is not None and np.asarray(res_ids).shape[0] != n_res:
            raise ValueError(
                f"{row_label} has {n_res} residues in its score profile but "
                f"{np.asarray(res_ids).shape[0]} residue numbers."
            )
        res_ids = None if res_ids is None else np.asarray(res_ids)
        x = np.arange(n_res)

        for key in PROFILE_SERIES:
            ax.plot(x, profiles[key].values, label=SCORE_DISPLAY_NAMES[key],
                    color=SCORE_PROFILE_COLOURS[key], **_PROFILE_STYLE[key])

        # Secondary axis: pLDDT, on its own 0-100 scale.
        ax_plddt = ax.twinx()
        ax_plddt.fill_between(x, plddt_chain, alpha=0.15, color='grey', label='pLDDT')
        ax_plddt.set_ylabel('pLDDT (grey fill)', color='grey', fontsize=10)
        ax_plddt.set_ylim(0, 100)
        ax_plddt.tick_params(axis='y', labelcolor='grey')
        ax_plddt.set_zorder(0)
        ax.set_zorder(1)
        ax.patch.set_visible(False)

        ax.fill_between(x, 0, 1, where=if_mask, alpha=0.12, color=COLOUR_IF,
                        label='Interface region')
        _annotate_peaks(ax, profiles, n_res, res_ids)
        ax.plot([], [], linestyle='none', marker='*', markersize=11,
                color=PEAK_MARKER_COLOUR, markeredgecolor='white',
                markeredgewidth=0.8,
                label='peak residue = this direction’s reported score')

        ax.set_ylim(0, 1)
        ax.set_xlim(-0.5, n_res - 0.5)
        ax.set_xlabel(f'Residue index along {row_label}')
        ax.set_ylabel('Per-residue score (0–1)')
        ax.set_title(f'Per-Residue Score Profiles: {row_label} → {col_label}')
        ax.legend(loc='upper left', fontsize=9, ncol=2)

    fig.tight_layout()
    return fig


def plot_plddt_distribution(
    contacts: InterfaceContacts,
    plddt_x: np.ndarray,
    plddt_y: np.ndarray,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
    figsize: Tuple[float, float] = PLDDT_FIGSIZE,
) -> Figure:
    """
    Interface pLDDT against the rest of the model, as a histogram and a profile.

    pDockQ and pDockQ2 both average pLDDT over interface residues, so a low score
    has two very different causes: a globally uncertain protein, or a confident
    protein with an uncertain interface. The upper panel separates them; the
    lower panel says *where* the uncertain residues are, in AlphaFold's own
    colours.

    The two panels are **stacked, not side by side**. They are two independent
    views of the same residues rather than one figure read across, so neither
    needs to sit beside the other, and the lower one is a per-residue track whose
    readability is purely a matter of pixels per residue: side by side it got
    half of a figure that the output area then scaled to 43%, and stacked it gets
    the full column. See `PLDDT_FIGSIZE` for the sizing.

    Args:
        contacts: Interface contacts of one ordered chain pair.
        plddt_x:  `(nx,)` pLDDT for `chain_x`.
        plddt_y:  `(ny,)` pLDDT for `chain_y`.
        label_x:  Display name for `chain_x`.
        label_y:  Display name for `chain_y`.
        figsize:  Figure size in inches.

    Returns:
        The `Figure`.

    Raises:
        ValueError: If a pLDDT array's length does not match its chain's, which
            would silently pair each residue with somebody else's confidence.
    """
    plddt_x = np.asarray(plddt_x)
    plddt_y = np.asarray(plddt_y)
    nx, ny = contacts.contact_mask.shape
    if plddt_x.shape[0] != nx or plddt_y.shape[0] != ny:
        raise ValueError(
            f"Chains are ({nx}, {ny}) residues but pLDDT arrays are "
            f"({plddt_x.shape[0]}, {plddt_y.shape[0]})."
        )

    # These two strings read as sequence landmarks ("A, then B") rather than as
    # panel headings, so they want the compact form; `ChainLabel.short` already
    # is that form, and the bare chain id is its fallback.
    name_x = _chain_label(contacts.chain_x, label_x)
    name_y = _chain_label(contacts.chain_y, label_y)

    if_plddt = np.concatenate([plddt_x[contacts.mask_x], plddt_y[contacts.mask_y]])
    ni_plddt = np.concatenate([plddt_x[~contacts.mask_x], plddt_y[~contacts.mask_y]])

    fig = Figure(figsize=figsize, layout='constrained')
    axes = fig.subplots(2, 1, height_ratios=list(PLDDT_PANEL_RATIOS))

    ax = axes[0]
    bins = np.linspace(0, 100, 26)
    ax.hist(ni_plddt, bins=bins, alpha=0.55, color=COLOUR_NON_IF_HIST,
            label=f'Non-interface (n={len(ni_plddt)})', density=True)
    ax.hist(if_plddt, bins=bins, alpha=0.65, color=COLOUR_IF,
            label=f'Interface (n={len(if_plddt)})', density=True)
    ax.axvline(70, color='grey', linestyle='--', linewidth=1,
               label='pLDDT = 70 (threshold)')
    ax.set_xlabel('pLDDT score')
    ax.set_ylabel('Density')
    ax.set_title('pLDDT distribution: interface vs non-interface')
    ax.legend(fontsize=9)

    ax2 = axes[1]
    full_plddt = np.concatenate([plddt_x, plddt_y])
    full_x = np.arange(len(full_plddt))
    full_if = np.concatenate([contacts.mask_x, contacts.mask_y])

    ax2.bar(full_x, full_plddt,
            color=[plddt_band_colour(v) for v in full_plddt],
            width=1.0, linewidth=0)
    ax2.fill_between(full_x, 0, 100, where=full_if, alpha=0.18, color=COLOUR_IF,
                     label='Interface residues')
    ax2.axvline(nx - 0.5, color='black', linewidth=1.5, linestyle='--',
                label=f'{name_x} | {name_y} boundary')
    ax2.axhline(70, color='grey', linestyle=':', linewidth=1)
    ax2.set_xlabel(f'Residue index ({name_x}, then {name_y})')
    ax2.set_ylabel('pLDDT')
    ax2.set_title('Per-residue pLDDT profile (AlphaFold colour scheme)')

    legend_patches = [
        mpatches.Patch(color=colour, label=label)
        for _, colour, label in PLDDT_BANDS
    ]
    legend_patches.append(
        mpatches.Patch(color=COLOUR_IF, alpha=0.5, label='Interface region'))
    ax2.legend(handles=legend_patches, fontsize=8, loc='lower right', ncol=2)

    return fig


def plot_score_agreement(
    scores: Mapping[str, float],
    score_names: Optional[Tuple[str, ...]] = None,
    labels: Optional[Mapping[str, str]] = None,
    cap: float = 1.5,
    cmap: str | Colormap = AGREEMENT_CMAP,
    figsize: Tuple[float, float] = (7.0, 6.0),
) -> Figure:
    """
    How far apart the scores are once each is expressed in units of its own cutoff.

    **Superseded by `plot_threshold_margins` (R082), which is what Section 7 now
    draws.** Kept because it is a working public function and a caller may still
    want the matrix form, but read that function's docstring before choosing it:
    over scalar scores `|a_i - a_j|` is a distance matrix over points on a line,
    so every off-diagonal cell here is determined by the handful of normalised
    values, and `cap` makes two scores that are both far past their thresholds
    register as agreeing perfectly.

    The scores are on incomparable scales, so `0.30` means something different for
    LIS than for pDockQ. Dividing each by its own green threshold puts them all in
    "fraction of the bar it has to clear", after which `|a - b|` is a meaningful
    disagreement. A dark red cell is the interesting case: two methods reading the
    same model and reaching opposite verdicts.

    Args:
        scores:      `{THRESHOLDS key: value}`. Extra keys are ignored.
        score_names: Which values to compare, in display order; defaults to
                     `AGREEMENT_SCORES`.
        labels:      Display names; defaults to `SCORE_DISPLAY_NAMES`.
        cap:         Upper clamp on the normalised value, so one score far past
                     its threshold cannot flatten the rest of the scale.
        cmap:        Colormap; green = agreement.
        figsize:     Figure size in inches.

    Returns:
        The `Figure`.

    Raises:
        KeyError: If a requested score is missing from `scores`, or is not a
            `THRESHOLDS` key -- there would be no cutoff to normalise it by.

    Note:
        Normalisation reads `THRESHOLDS[name].green`, the canonical table, so the
        figure inherits whatever uncertainty those cutoffs carry -- four of the
        seven are DERIVED or HEURISTIC. Unlike `plot_threshold_margins` this
        figure has nowhere to say so, which is the second half of why R082
        replaced it.
    """
    names = AGREEMENT_SCORES if score_names is None else tuple(score_names)
    display = dict(SCORE_DISPLAY_NAMES if labels is None else labels)

    missing = [name for name in names if name not in scores]
    if missing:
        raise KeyError(f"No value supplied for: {', '.join(missing)}.")
    unknown = [name for name in names if name not in THRESHOLDS]
    if unknown:
        raise KeyError(
            f"No canonical threshold for: {', '.join(unknown)}. "
            f"Expected keys from THRESHOLDS: {', '.join(sorted(THRESHOLDS))}."
        )

    normalised = np.array(
        [min(float(scores[name]) / THRESHOLDS[name].green, cap) for name in names]
    )
    agreement = np.abs(normalised[:, np.newaxis] - normalised[np.newaxis, :])
    n_s = len(names)

    fig = Figure(figsize=figsize)
    ax = fig.subplots()
    im = ax.imshow(agreement, cmap=cmap, vmin=0, vmax=1)
    fig.colorbar(im, ax=ax,
                 label='Normalised disagreement (0=agree, 1=strongly disagree)')

    tick_labels = [display.get(name, name) for name in names]
    ax.set_xticks(range(n_s))
    ax.set_yticks(range(n_s))
    ax.set_xticklabels(tick_labels, rotation=30, ha='right')
    ax.set_yticklabels(tick_labels)
    for i in range(n_s):
        for j in range(n_s):
            ax.text(j, i, f'{agreement[i, j]:.2f}', ha='center', va='center',
                    fontsize=9, color='black')
    ax.set_title('Pairwise Score Agreement Matrix\n'
                 '(dark green = agree, dark red = strongly disagree)')

    fig.tight_layout()
    return fig


THRESHOLD_MARGIN_CAP: float = 2.0
"""Right-hand limit of the `plot_threshold_margins` axis, in units of the green
threshold. Unlike the agreement matrix's `cap`, this clips only the *drawing* of
a bar: the raw score and its exact margin are printed on every row, so a value
past the edge is still readable as a number."""


def plot_threshold_margins(
    scores: Mapping[str, float],
    score_names: Optional[Sequence[str]] = None,
    labels: Optional[Mapping[str, str]] = None,
    axis_cap: float = THRESHOLD_MARGIN_CAP,
    figsize: Tuple[float, float] = (7.4, 3.6),
) -> Figure:
    """
    Each score's distance from its own green threshold, on one shared axis.

    This replaces the pairwise agreement matrix (R082). The matrix plotted
    `|a_i - a_j|` over scalars, and a distance matrix over scalars is a distance
    matrix over points on a line: all 10 of its off-diagonal numbers are
    determined by the 5 normalised values, so it spent 25 cells re-encoding 5.
    It also capped the normalised values at 1.5, which made any two scores far
    past their thresholds register as *perfect* agreement -- precisely the case
    on an AFDB accession, which has already passed the release filter. The
    matrix was therefore closest to blank exactly where the notebook is most
    often used.

    Plotting the 5 normalised values directly keeps every disagreement the
    matrix could show -- the gap between any two rows *is* the pairwise
    `|a_i - a_j|`, read straight off the axis -- and restores three things the
    matrix destroyed: which side of its threshold each score falls on, how far
    past it a score goes, and what each score's "1.0" is actually worth. The
    green line at 1.0 is 0.70 for ipSAE and 0.23 for pDockQ, and one of those is
    published while the other is derived; the matrix hid both facts inside a
    single number.

    The normalisation is still `value / THRESHOLDS[name].green`, so the figure
    still depends on the thresholds being right. The difference is that it now
    says so: the axis is labelled with what the division is, and every row
    prints its own green edge and that edge's provenance.

    Args:
        scores:      `{THRESHOLDS key: value}`. Extra keys are ignored.
        score_names: Which values to show, in display order; defaults to
                     `AGREEMENT_SCORES`.
        labels:      Display names; defaults to `SCORE_DISPLAY_NAMES`.
        axis_cap:    Right-hand axis limit, in units of the green threshold.
                     Bars are drawn clipped to it; the printed numbers are not.
        figsize:     Figure size in inches.

    Returns:
        The `Figure`.

    Raises:
        KeyError: If a requested score is missing from `scores`, or is not a
            `THRESHOLDS` key -- there would be no cutoff to normalise it by.

    Example
    -------
    >>> fig = plot_threshold_margins({'ipsae_d0res': 0.64, 'iptm_d0chn': 0.82,
    ...                               'pdockq': 0.24, 'pdockq2': 0.33,
    ...                               'lis': 0.48})
    >>> len(fig.axes[0].get_yticklabels())
    5
    """
    names = tuple(AGREEMENT_SCORES) if score_names is None else tuple(score_names)
    display = dict(SCORE_DISPLAY_NAMES if labels is None else labels)

    missing = [name for name in names if name not in scores]
    if missing:
        raise KeyError(f"No value supplied for: {', '.join(missing)}.")
    unknown = [name for name in names if name not in THRESHOLDS]
    if unknown:
        raise KeyError(
            f"No canonical threshold for: {', '.join(unknown)}. "
            f"Expected keys from THRESHOLDS: {', '.join(sorted(THRESHOLDS))}."
        )

    fig = Figure(figsize=figsize)
    ax = fig.subplots()

    # Top row at the top: matplotlib counts y upwards, so the display order is
    # reversed here rather than in the caller.
    order = list(reversed(names))
    y_pos = np.arange(len(order), dtype=float)

    for y, name in zip(y_pos, order):
        threshold = THRESHOLDS[name]
        value = float(scores[name])
        margin = value / threshold.green
        band = threshold.confidence_band(value)
        colour = BAND_COLOUR_HEX[band.colour]

        clipped = margin > axis_cap
        drawn = min(margin, axis_cap)
        # A thin bar plus an end mark: the bar carries the magnitude, the mark
        # says whether that end is the value or the edge of the axis. A round
        # dot means the value is where it looks; an arrowhead means the score
        # runs past the axis and the printed number is the one to read.
        ax.plot([0, drawn], [y, y], color=colour, linewidth=3.0,
                solid_capstyle='butt' if clipped else 'round', zorder=3)
        ax.plot([drawn], [y], marker='>' if clipped else 'o',
                markersize=10 if clipped else 9, color=colour,
                markeredgecolor='white', markeredgewidth=1.4, zorder=4)

        # The band name travels with the colour, so the verdict is never
        # carried by colour alone. A label that would run off the right-hand
        # edge is set inside the bar instead, stroked in white so it stays
        # legible over the fill.
        label_text = f'{value:.3f}  {band.label}  ({margin:.2f}x)'
        # A soft knockout behind the text, so neither the grid nor the
        # threshold rule shows through the glyphs.
        knockout = dict(facecolor='white', edgecolor='none', alpha=0.8,
                        boxstyle='square,pad=0.15')
        if drawn + 0.03 + 0.011 * len(label_text) < axis_cap:
            ax.text(drawn + 0.05, y, label_text, va='center', ha='left',
                    fontsize=8, color='#333', zorder=6, bbox=knockout)
        else:
            # Sit the label above the bar rather than on it: a white stroke
            # over a coloured fill is legible but ugly, and there is a clear
            # row of space here because the rows are a whole unit apart.
            ax.text(drawn - 0.03, y + 0.24, label_text, va='bottom',
                    ha='right', fontsize=8, color='#333', zorder=6,
                    bbox=knockout)

        # The amber edge, in the same units, so the three bands are visible per
        # row rather than only the one the score landed in.
        ax.plot([threshold.amber / threshold.green], [y], marker='|',
                markersize=11, color='#9E9E9E', markeredgewidth=1.4, zorder=2)

    ax.axvline(1.0, color='#424242', linewidth=1.2, zorder=1)
    ax.text(1.03, -0.64, 'green threshold', fontsize=8, color='#424242',
            ha='left', va='bottom')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [
            f'{display.get(name, name)}\n'
            f'green {THRESHOLDS[name].green:g} · '
            f'{THRESHOLDS[name].green_provenance.lower()}'
            for name in order
        ],
        fontsize=8,
    )
    ax.set_xlim(0, axis_cap)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel('score ÷ its own green threshold\n'
                  '(1.0 = exactly at the cutoff; grey tick = the amber edge)',
                  fontsize=9)
    ax.set_title('Distance from threshold, all five on one scale', fontsize=11)
    ax.grid(axis='x', linestyle=':', linewidth=0.6, alpha=0.5)
    ax.set_axisbelow(True)
    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# MolViewSpec views
# ---------------------------------------------------------------------------
# One builder function per 3D view, plus the shared `show_mol_view` helper.
#
# `molviewspec` MUST be imported lazily inside each builder, never at module top
# level, so that importing this module stays free and so the notebook degrades
# gracefully with a clear message when the package is absent. `molviewspec_available()`
# is the capability check a caller uses to decide whether to draw the section at
# all; `_require_molviewspec()` is what every builder calls on entry.
#
# Three deliberate properties of this section, mirroring the plotting section:
#
# 1. **Builders return a `State`, they do not render.** `show_mol_view` is the
#    only function that touches `IPython.display`, so a builder can be exercised
#    head-lessly by a test and a caller decides when and how big to draw.
# 2. **Chain-pair-generic (D4).** No `'A'` / `'B'` literals. Every builder takes
#    the ordered pair from the `InterfaceContacts` it is handed and looks the
#    chain up in the `{chain_id: ChainCoords}` mapping `parse_structure` returns,
#    so a complex whose chains are labelled `A` and `C` renders correctly.
# 3. **Per-residue colouring goes through one primitive.** `colour_runs` plus
#    `add_residue_colours` is the only place a residue is turned into a
#    component, and `value_colours` is the only place a number is turned into a
#    colour. That pairing -- an arbitrary per-residue value array plus a
#    colormap -- is exactly what View 3 does today and exactly what R072 and
#    R074 need, so the two new pDockQ2 views are calls to
#    `build_interface_value_view`, not new colouring code.
#
# Run collapsing (R075)
# ---------------------
# `ComponentExpression` has always taken `beg_label_seq_id` / `end_label_seq_id`;
# the notebook simply always passed `beg == end`, so View 2 emitted one component
# per residue -- 282 of them on the heterodimer fixture. `colour_runs` merges
# residues that are *both* consecutive in `label_seq_id` *and* the same colour
# into one ranged expression. Consecutiveness is required, not just equal colour,
# because a ranged selector covers everything between the two endpoints: merging
# across a gap in the numbering would paint residues that were never in the
# input. A residue with no colour therefore also breaks a run, which is what
# keeps the pLDDT view's uncoloured residues uncoloured.
#
# What this section deliberately does NOT do
# ------------------------------------------
# R014 is a behaviour-preserving move, on the same reasoning as R013: a 3D view
# has no objective oracle, so a visual change landed here would be
# indistinguishable from a regression. Preserved exactly as the notebook draws
# them today, and owned by M6:
#
# - side chains on the first chain only in Views 1 and 4 (R071, R073);
# - the bare "red=low, yellow=mid, green=high" colour scale of View 3, with no
#   legend and no stated numbers (R072);
# - View 4's magic `0.5` cutoff, unrelated to any threshold in `THRESHOLDS`, its
#   colour-key category names, and its missing legend (R073);
# - no pDockQ2 view at all (R074);
# - no supporting prose anywhere (R070).
#
# The seams those tasks change are named on each function: `side_chain_colours`,
# `bands`, `cmap` / `vmin` / `vmax`, `threshold`, `categories` and
# `MVS_VIEW_LABELS`.
#
# Four implementation defects were fixed in R014, because each was a defect
# rather than a design choice and none of them changed a pixel: the deprecated
# `matplotlib.cm.get_cmap`, the one-component-per-residue payload, the hard-coded
# chain letters, and the structure URL and format being resolved from two
# independent expressions that disagree when `bcifUrl` is present but empty.
#
# The fifth, the one that did change pixels, is fixed in R075: the module used
# to carry two disagreeing pLDDT ladders, one for the figures and one for the
# 3D view. There is now one, `PLDDT_BANDS`, and `MVS_PLDDT_BANDS` is an alias of
# it. Every view also carries a legend built from the constants it is coloured
# from (R073, R074), which is what made reconciling the ladders compulsory
# rather than tidy.


# -- view constants ---------------------------------------------------------
# Every colour, cutoff and label the four views use, named once. Section 6 of the
# notebook currently spells all of these inline; naming them here is what makes
# R070-R074 edits to *this list* rather than edits inside the builders.

MVS_VIEW_WIDTH: str | int = "100%"
"""Default viewer width: a CSS length, or an `int` read as pixels.

`'100%'` so the viewer fills whatever column it is rendered into, rather than
being clipped on a narrow screen and stranded in white space on a wide one. Was
`950` px, which is a guess at one particular window. Pass `width=` to
`show_mol_view` to override per call."""

MVS_VIEW_HEIGHT: str | int = 600
"""Default viewer height: a CSS length, or an `int` read as pixels.

**Stays a pixel value.** An iframe in normal document flow has no height to take
a percentage *of*, so `'100%'` collapses it to nothing. Only the width can be
made responsive."""

MVS_CONTEXT_COLOUR: str = "#BDBDBD"
"""Mid grey. The whole-complex cartoon behind View 3's coloured interface, dark
enough to read as structure and light enough not to compete with the colours."""

MVS_FAINT_COLOUR: str = "#EEEEEE"
"""Near-white. The whole-complex cartoon behind the two categorical views,
fainter than `MVS_CONTEXT_COLOUR` because those views' categories are their
entire message."""

MVS_VALUE_CMAP: str = "RdYlGn"
"""Colormap for any per-residue score painted onto a structure: red = low,
green = high. Resolved through `matplotlib.colormaps`, never through the
deprecated `matplotlib.cm.get_cmap` the notebook calls (R075).

`value_ramp_legend` samples this same ramp for the legend under Views 3 and 5,
so the numbers behind "low" and "high" are on screen. This constant is the seam
if the ramp itself changes."""

MVS_PLDDT_BANDS: Tuple[Tuple[float, str, str], ...] = PLDDT_BANDS
"""pLDDT bands for the 3D view. **Now literally `PLDDT_BANDS` (R075).**

There used to be two ladders. This one tested `low <= value < high` over
half-open intervals, so a residue at exactly 100.00 matched no band and was not
drawn at all, and 90.0 landed in the top band; `plddt_band_colour` tested
`value > 90` and put 90.0 in the second band. R014 preserved both rather than
silently reconciling them, because a 3D view has no oracle and changing either
shifts colours in a figure nobody had reviewed.

They are reconciled here, in favour of the `>` ladder, because that is the one
AlphaFold publishes: the AlphaFold DB entry page and the AlphaFold FAQ both
state the bands as "Very high (pLDDT > 90)", "Confident (90 > pLDDT > 70)",
"Low (70 > pLDDT > 50)" and "Very low (pLDDT < 50)". The half-open version was
an artefact of the notebook spelling the ladder out a second time, not a
different reading of AlphaFold. Two consequences, both wanted: a residue at
100.00 is now drawn, and 90.0 is confident rather than very high.

The name is kept because it is exported and because `build_plddt_view` still
takes a `bands` argument; it is an alias, not a copy, so the two cannot drift
apart again."""

MVS_DISAGREEMENT_THRESHOLD: float = THRESHOLDS["ipsae_d0res"].amber
"""Per-residue ipSAE_d0res at or above which View 4 calls PAE "confident".

**Sourced from `THRESHOLDS`, not invented (R073).** It used to be a hard-coded
0.5 unrelated to any published cutoff. It is now
`THRESHOLDS['ipsae_d0res'].amber`, currently **0.60**: AFDB's release cutoff,
the same edge below which the Section 7 traffic light turns red, and the ipSAE
side of the AFDB joint high-confidence criterion. Sourcing it means the view and
the summary can never disagree about what "confident" means, and that raising
the cutoff in one place moves both.

The test is `>=`, matching `traffic_light`, not the `>` the notebook used."""

MVS_CONTACT_PTM_THRESHOLD: float = 0.5
"""Per-residue mean `ptm(PAE, d0=10)` at or above which View 6 calls a residue's
own contacts well-placed (R074).

Not a fourth magic number: `ptm(PAE, d0=10) = 1 / (1 + (PAE/10)**2)` is exactly
0.5 when `PAE == 10`, and 10 A is `PAE_CUTOFF`, the same inter-chain PAE cutoff
ipSAE uses to decide a pair is informative at all. So this cutoff says "the mean
PAE over this residue's contacts is better than ipSAE's own admission
threshold", which is the comparison View 6 is built to draw. Verified by
construction: `ptm_func(PAE_CUTOFF, 10.0) == 0.5`."""

MVS_DISAGREEMENT_CATEGORIES: Tuple[Tuple[str, str, str], ...] = (
    ("confirmed", "#009E73",
     "Confirmed contact: CB within {d:.0f} A AND ipSAE_d0res >= {t:.2f}"),
    ("predicted_not_touching", "#0072B2",
     "Predicted but not touching: ipSAE_d0res >= {t:.2f}, no CB within {d:.0f} A"),
    ("touching_not_trusted", "#A02020",
     "Touching but not trusted: CB within {d:.0f} A, ipSAE_d0res < {t:.2f}"),
)
"""View 4's categories, as `(key, colour, label template)`.

**Named for what the reader should do about them, not for their colour (R073).**
The notebook's caption read "green=PAE+contact, blue=PAE confident/no contact,
red=contact/low PAE", which names the colours and the inputs and never says
whether any of it is good news. Each label here is a sentence stating both
conditions and the numbers behind them; `{t}` is filled with the ipSAE cutoff in
force and `{d}` with the contact distance cutoff, by `format_category_labels`.

**Colours changed, with a measurement.** The old triple was Material green
`#4CAF50`, blue `#2196F3` and red `#F44336`. Under deuteranopia the green and
the red simulate to `#A59857` and `#A4932E`, CIEDE2000 `dE = 6.0` apart, so the
two categories a reader most needs to tell apart were the two that vanished
into each other. These three are Okabe-Ito's bluish green and blue plus a dark
red chosen for luminance separation; the worst pair, including the faint base
cartoon, is `dE = 21.5`.

A residue that is neither confident nor in contact matches no category and is
left undrawn. `MVS_UNLIT_LABELS` carries the sentence that says so, because a
legend that lists only what is painted invites the reader to assume the rest was
not analysed."""

MVS_PDOCKQ2_AGREEMENT_CATEGORIES: Tuple[Tuple[str, str, str], ...] = (
    ("both", "#009E73",
     "Both agree, good: contact PAE better than {p:.0f} A AND "
     "ipSAE_d0res >= {t:.2f}"),
    ("contacts_only", "#0072B2",
     "pDockQ2 only: this residue's own contacts are well placed "
     "(mean contact PAE better than {p:.0f} A) but ipSAE_d0res < {t:.2f}, "
     "so it is not confident about the partner chain as a whole"),
    ("ipsae_only", "#F0E442",
     "ipSAE only: confident about the partner chain (ipSAE_d0res >= {t:.2f}) "
     "but its own contacts are among the worse-placed "
     "(mean contact PAE worse than {p:.0f} A)"),
    ("neither", "#A02020",
     "Both agree, poor: mean contact PAE worse than {p:.0f} A AND "
     "ipSAE_d0res < {t:.2f}"),
)
"""View 6's categories, as `(key, colour, label template)` (R074).

Two per-residue signals answer the same question from different distances.
`ipSAE_d0res` measures a residue against the **whole** partner chain, over every
pair that clears `PAE_CUTOFF`. pDockQ2's `mean_ptm_by_residue` measures it
against **only the partner residues it actually touches**. The middle two
categories are the residues where those two answers differ, which is the "why do
the scores disagree" question the notebook exists to answer, localised onto the
structure.

Only interface residues are classified: a residue with no contacts has no
pDockQ2 value at all (`np.nan`), so there is nothing to compare and it is left
undrawn. `{t}` is the ipSAE cutoff, `{p}` the PAE in Angstrom equivalent to
`MVS_CONTACT_PTM_THRESHOLD`, both filled by `format_category_labels`.

Colours are Okabe-Ito plus the same dark red as View 4, with `confirmed` /
`both` deliberately sharing `#009E73` so the two agreement views read as one
family. Worst pair including the base cartoon, over normal, deuteranopia and
protanopia vision: `dE = 21.5`."""

MVS_UNLIT_LABELS: Dict[str, str] = {
    "disagreement": ("Not drawn: no CB within {d:.0f} A and "
                     "ipSAE_d0res < {t:.2f}, so neither signal fires"),
    "pdockq2_agreement": ("Not drawn: no CB within {d:.0f} A, so this residue "
                          "has no pDockQ2 value to compare"),
}
"""What the residues a category view leaves unpainted actually mean.

Rendered as the last legend entry, with no swatch. Both category views cover
only part of the structure, and the reason differs between them; saying so is
the difference between "not analysed" and "analysed, and the answer was no"."""

MVS_VIEW_LABELS: Dict[str, str] = {
    "chain_overview": (
        "View 1: Chain overview. Cartoons: {x} teal, {y} coral. "
        "Interface side chains: {x} gold, {y} cornflower."),
    "plddt": (
        "View 2: pLDDT per residue, {x} and {y}, in AlphaFold's own bands."),
    "interface_value": (
        "View 3: ipSAE d0res per interface residue, {x} and {y}. "
        "Red 0.00 to green 1.00."),
    "disagreement": (
        "View 4: ipSAE d0res against physical contact, {x} and {y}. "
        "Where the two signals for 'this residue is at the interface' differ."),
    "pdockq2_value": (
        "View 5: pDockQ2 contact quality per interface residue, {x} and {y}. "
        "Mean ptm of the PAE at each residue's own contacts, red 0.00 to "
        "green 1.00."),
    "pdockq2_agreement": (
        "View 6: pDockQ2 contact quality against ipSAE d0res, {x} and {y}. "
        "Where a residue's own contacts and its view of the whole partner "
        "chain disagree."),
}
"""The caption `show_mol_view` draws above each view.

Templates, not finished strings: `{x}` and `{y}` are filled by
`format_view_label` with the two chains' display names, so a caption says which
*proteins* are teal and coral rather than the hard-coded `'A'` and `'B'` the
notebook used to spell out (R021). Every template names both chains, because
every view now draws both (R071, R073).

These are one-line captions, not explanations. The legend under each view
carries the colour key with the cutoffs in force, and R070 owns the supporting
prose above each view."""

def format_view_label(
    key: str,
    label_x: "str | ChainLabel",
    label_y: "str | ChainLabel",
    labels: Mapping[str, str] = MVS_VIEW_LABELS,
) -> str:
    """
    Fill one `MVS_VIEW_LABELS` template with the chain pair's display names.

    Args:
        key:     A key of `labels`, e.g. `'chain_overview'`.
        label_x: Display name for the first chain of the ordered pair. A
                 `ChainLabel` renders as its `short` form.
        label_y: Display name for the second.
        labels:  Template table; defaults to `MVS_VIEW_LABELS`.

    Returns:
        The caption, ready for `show_mol_view`.

    Raises:
        KeyError: If `key` is not in `labels`.

    Example
    -------
    >>> format_view_label('chain_overview',
    ...                   ChainLabel('A', gene='ISG20'),
    ...                   ChainLabel('B', gene='Sumo1'))
    'View 1: Chain overview. Cartoons: ISG20 (A) teal, Sumo1 (B) coral. Interface side chains: ISG20 (A) gold, Sumo1 (B) cornflower.'

    Every view names both chains, because every view draws both:

    >>> format_view_label('disagreement', 'ISG20 (A)', 'Sumo1 (B)')[:64]
    'View 4: ipSAE d0res against physical contact, ISG20 (A) and Sumo'
    """
    try:
        template = labels[key]
    except KeyError:
        raise KeyError(
            f"No view label for {key!r}; known views: {', '.join(sorted(labels))}."
        ) from None
    return template.format(x=str(label_x), y=str(label_y))

# -- legends ----------------------------------------------------------------
# Every 3D view now carries a legend, and every legend is built from the same
# constants the view is coloured from, so the two cannot disagree. That is the
# whole reason the pLDDT ladders had to be reconciled first (R075): a legend
# that contradicts the colours it explains is worse than a silent
# inconsistency, because it looks authoritative.
#
# `LegendEntry` is deliberately dumber than a matplotlib colour bar. These
# legends sit outside the Mol* iframe as ordinary HTML, so they survive the
# viewer failing to load, they can be read by a screen reader, and they can be
# printed as text when `IPython` is not available.


@dataclass(frozen=True)
class LegendEntry:
    """
    One row of a view's legend.

    Attributes:
        colour: `'#RRGGBB'` swatch, or `''` for a row that describes something
                not painted at all (the residues a category view leaves out).
        label:  The sentence shown next to the swatch.
        counts: `((chain display name, n), ...)`, empty when a row has no
                population. Shown as `(ISG20 (A): 41, Sumo1 (B): 45)`.
    """

    colour: str
    label: str
    counts: Tuple[Tuple[str, int], ...] = ()

    @property
    def total(self) -> int:
        """Residues in this row across every chain."""
        return sum(n for _name, n in self.counts)


def contact_ptm_to_pae(ptm: float, d0: float = 10.0) -> float:
    """
    The PAE, in Angstrom, that a given `ptm(PAE, d0)` corresponds to.

    The inverse of `ptm_func`, used only to state `MVS_CONTACT_PTM_THRESHOLD` in
    the units a reader thinks in. `ptm = 1 / (1 + (PAE/d0)**2)`, so
    `PAE = d0 * sqrt(1/ptm - 1)`.

    Args:
        ptm: A ptm value in `(0, 1]`.
        d0:  The normalisation the ptm was computed with. pDockQ2 fixes it at 10.

    Returns:
        The PAE in Angstrom.

    Raises:
        ValueError: If `ptm` is outside `(0, 1]`, which has no finite inverse.

    Example
    -------
    >>> contact_ptm_to_pae(0.5)
    10.0
    >>> round(contact_ptm_to_pae(0.8), 3)
    5.0
    """
    if not 0.0 < ptm <= 1.0:
        raise ValueError(f"ptm must be in (0, 1]; got {ptm}.")
    return float(d0 * math.sqrt(1.0 / ptm - 1.0))


def format_category_labels(
    categories: Sequence[Tuple[str, str, str]],
    threshold: float,
    dist_cutoff: float = DIST_CUTOFF,
    contact_ptm: Optional[float] = None,
) -> Tuple[Tuple[str, str, str], ...]:
    """
    Fill a category table's label templates with the cutoffs actually in force.

    The category tables carry `{t}` for the ipSAE cutoff, `{d}` for the contact
    distance cutoff and `{p}` for the contact-PAE cutoff. Filling them here,
    from the same values the classifier is called with, is what stops a legend
    quoting a number the view did not use.

    Args:
        categories:  `(key, colour, label template)` triples.
        threshold:   Fills `{t}`. The ipSAE_d0res cutoff in force.
        dist_cutoff: Fills `{d}`. The CB-CB contact cutoff.
        contact_ptm: The `ptm` cutoff whose PAE equivalent fills `{p}`. Required
                     only for tables that use `{p}`.

    Returns:
        `(key, colour, finished label)` triples, in input order.

    Raises:
        KeyError: If a template uses a placeholder that was not supplied, which
            is a table/caller mismatch rather than a formatting nicety.

    Example
    -------
    >>> format_category_labels(MVS_DISAGREEMENT_CATEGORIES, 0.6)[0][2]
    'Confirmed contact: CB within 8 A AND ipSAE_d0res >= 0.60'
    >>> format_category_labels(MVS_PDOCKQ2_AGREEMENT_CATEGORIES, 0.6,
    ...                        contact_ptm=0.5)[0][2]
    'Both agree, good: contact PAE better than 10 A AND ipSAE_d0res >= 0.60'
    """
    pae = None if contact_ptm is None else contact_ptm_to_pae(contact_ptm)
    filled = []
    for key, colour, template in categories:
        filled.append((key, colour, template.format(t=threshold, d=dist_cutoff, p=pae)))
    return tuple(filled)


def category_legend(
    categories: Sequence[Tuple[str, str, str]],
    assigned: Mapping[str, np.ndarray],
    threshold: float,
    dist_cutoff: float = DIST_CUTOFF,
    contact_ptm: Optional[float] = None,
    unlit_label: Optional[str] = None,
) -> List[LegendEntry]:
    """
    Legend rows for a category view, with the per-chain population of each row.

    Args:
        categories:  `(key, colour, label template)` triples.
        assigned:    `{chain display name: (n,) array of category keys}`, in the
                     order the chains should be listed. `''` marks a residue in
                     no category.
        threshold:   The ipSAE cutoff in force; fills `{t}`.
        dist_cutoff: The contact cutoff; fills `{d}`.
        contact_ptm: The contact-ptm cutoff; its PAE equivalent fills `{p}`.
        unlit_label: Template for the trailing "not drawn" row. Omit to leave it
                     out; the row is added even when its population is zero, so
                     the legend accounts for every residue either way.

    Returns:
        One `LegendEntry` per category, in table order, then the unlit row.

    Example
    -------
    >>> import numpy as np
    >>> rows = category_legend(
    ...     MVS_DISAGREEMENT_CATEGORIES,
    ...     {'A': np.array(['confirmed', 'confirmed', '']),
    ...      'B': np.array(['touching_not_trusted', '', ''])},
    ...     threshold=0.6, unlit_label=MVS_UNLIT_LABELS['disagreement'])
    >>> [(e.colour, e.total) for e in rows]
    [('#009E73', 2), ('#0072B2', 0), ('#A02020', 1), ('', 3)]
    >>> rows[0].counts
    (('A', 2), ('B', 0))
    """
    labelled = format_category_labels(categories, threshold, dist_cutoff, contact_ptm)
    arrays = {name: np.asarray(keys, dtype=object) for name, keys in assigned.items()}

    entries: List[LegendEntry] = []
    for key, colour, label in labelled:
        counts = tuple((name, int((arr == key).sum())) for name, arr in arrays.items())
        entries.append(LegendEntry(colour=colour, label=label, counts=counts))

    if unlit_label is not None:
        known = {key for key, _colour, _label in categories}
        counts = tuple(
            (name, int(sum(1 for k in arr if k not in known)))
            for name, arr in arrays.items()
        )
        entries.append(LegendEntry(
            colour='',
            label=unlit_label.format(t=threshold, d=dist_cutoff,
                                     p=None if contact_ptm is None
                                     else contact_ptm_to_pae(contact_ptm)),
            counts=counts,
        ))
    return entries


def plddt_legend(
    values: Optional[Mapping[str, "Sequence[float] | np.ndarray"]] = None,
    bands: Sequence[Tuple[float, str, str]] = PLDDT_BANDS,
) -> List[LegendEntry]:
    """
    Legend rows for View 2, from the one pLDDT ladder (R075).

    Reads `PLDDT_BANDS`, exactly as `build_plddt_view` does, so the legend and
    the colours are the same four rows.

    Args:
        values: `{chain display name: (n,) pLDDT array}` to populate the counts.
                Omit for a legend with no populations.
        bands:  The ladder; defaults to the module's single ladder.

    Returns:
        One `LegendEntry` per band, highest first.

    Example
    -------
    >>> rows = plddt_legend({'A': [95.0, 90.0, 100.0, 40.0]})
    >>> [(e.label, e.total) for e in rows][0]
    ('>90 (very high)', 2)
    >>> [e.total for e in rows]
    [2, 1, 0, 1]
    """
    arrays = {} if values is None else {
        name: np.asarray(vals, dtype=float) for name, vals in values.items()
    }
    entries: List[LegendEntry] = []
    for minimum, colour, label in bands:
        counts = tuple(
            (name, int(sum(1 for v in arr if plddt_band(float(v), bands)[1] == label)))
            for name, arr in arrays.items()
        )
        entries.append(LegendEntry(colour=colour, label=label, counts=counts))
    return entries


def value_ramp_legend(
    title: str,
    vmin: float = 0.0,
    vmax: float = 1.0,
    cmap: "str | Colormap" = MVS_VALUE_CMAP,
    n_steps: int = 5,
) -> List[LegendEntry]:
    """
    Legend rows sampling a continuous colormap, for the two value views.

    A discrete key rather than a colour bar, because it has to sit outside the
    Mol* iframe as plain HTML. Colours come from `value_colours`, the same
    function the view paints with, so the swatches are the ramp itself.

    Args:
        title:   Row prefix naming the quantity, e.g. `'ipSAE d0res'`.
        vmin:    Low end of the ramp.
        vmax:    High end.
        cmap:    Colormap; defaults to `MVS_VALUE_CMAP`.
        n_steps: Number of sampled swatches, at least 2.

    Returns:
        `n_steps` `LegendEntry` rows, high value first.

    Raises:
        ValueError: If `n_steps < 2`.

    Example
    -------
    >>> [(e.colour, e.label) for e in value_ramp_legend('ipSAE d0res', n_steps=3)]
    [('#006837', 'ipSAE d0res 1.00'), ('#FEFEBD', 'ipSAE d0res 0.50'), ('#A50026', 'ipSAE d0res 0.00')]
    """
    if n_steps < 2:
        raise ValueError(f"value_ramp_legend needs n_steps >= 2; got {n_steps}.")
    steps = [vmin + (vmax - vmin) * i / (n_steps - 1) for i in range(n_steps)][::-1]
    colours = value_colours(steps, cmap=cmap, vmin=vmin, vmax=vmax)
    return [LegendEntry(colour=colour, label=f'{title} {value:.2f}')
            for value, colour in zip(steps, colours)]


def chain_overview_legend(
    contacts: InterfaceContacts,
    label_x: str,
    label_y: str,
    chain_colours: Optional[Mapping[str, str]] = None,
    side_chain_colours: Optional[Mapping[str, str]] = None,
) -> List[LegendEntry]:
    """
    Legend rows for View 1: two cartoons and two sets of interface side chains.

    Args:
        contacts:           The ordered pair being drawn; supplies the interface
                            residue counts.
        label_x, label_y:   Display names for the two chains.
        chain_colours:      `{chain_id: colour}` as passed to the builder.
        side_chain_colours: `{chain_id: colour}` as passed to the builder.

    Returns:
        Four `LegendEntry` rows: cartoon x, cartoon y, side chains x, side
        chains y. The side-chain rows carry the interface residue counts, which
        is the number a reader wants and the proof that both chains are drawn.
    """
    chain_x, chain_y = contacts.chain_x, contacts.chain_y
    cartoon = (dict(zip((chain_x, chain_y), CHAIN_COLOURS))
               if chain_colours is None else dict(chain_colours))
    sides = (dict(zip((chain_x, chain_y), SIDE_CHAIN_COLOURS))
             if side_chain_colours is None else dict(side_chain_colours))
    names = {chain_x: label_x, chain_y: label_y}
    n_if = {chain_x: contacts.n_interface_residues_x,
            chain_y: contacts.n_interface_residues_y}

    entries = [
        LegendEntry(colour=cartoon[cid], label=f'{names[cid]} cartoon')
        for cid in (chain_x, chain_y) if cid in cartoon
    ]
    entries += [
        LegendEntry(colour=sides[cid],
                    label=f'{names[cid]} interface side chains',
                    counts=((names[cid], n_if[cid]),))
        for cid in (chain_x, chain_y) if cid in sides
    ]
    return entries


def _format_counts(counts: Sequence[Tuple[str, int]]) -> str:
    """`(A: 41, B: 45)`, or `''` when there are no counts."""
    if not counts:
        return ''
    return ' (' + ', '.join(f'{name}: {n}' for name, n in counts) + ')'


def legend_text(entries: Sequence[LegendEntry], title: Optional[str] = None) -> str:
    """
    A view's legend as plain text, for a terminal or a failed viewer.

    Args:
        entries: The rows.
        title:   Optional heading line.

    Returns:
        The legend, newline separated, with no trailing newline.

    Example
    -------
    >>> print(legend_text([LegendEntry('#009E73', 'Confirmed', (('A', 2),))],
    ...                   title='View 4'))
    View 4
      [#009E73] Confirmed (A: 2)
    """
    lines = [] if title is None else [title]
    for entry in entries:
        swatch = f'[{entry.colour}]' if entry.colour else '[not drawn]'
        lines.append(f'  {swatch} {entry.label}{_format_counts(entry.counts)}')
    return '\n'.join(lines)


def legend_html(entries: Sequence[LegendEntry], title: Optional[str] = None) -> str:
    """
    A view's legend as an HTML fragment, to `display(HTML(...))` under the view.

    Args:
        entries: The rows.
        title:   Optional heading.

    Returns:
        HTML. Inline styles only, so it renders identically in Jupyter, PyCharm
        and Colab, none of which share a stylesheet.

    Example
    -------
    >>> html = legend_html([LegendEntry('#009E73', 'Confirmed', (('A', 2),))])
    >>> '#009E73' in html and 'Confirmed (A: 2)' in html
    True
    """
    parts = ['<div style="margin:2px 0 14px; font-size:13px; line-height:1.6;">']
    if title is not None:
        parts.append(f'<div style="font-weight:bold; margin-bottom:3px;">{title}</div>')
    for entry in entries:
        if entry.colour:
            swatch = (f'<span style="display:inline-block; width:13px; height:13px; '
                      f'margin-right:7px; vertical-align:-2px; border:1px solid #999; '
                      f'background:{entry.colour};"></span>')
        else:
            swatch = ('<span style="display:inline-block; width:13px; height:13px; '
                      'margin-right:7px; vertical-align:-2px; border:1px dashed #999; '
                      'background:transparent;"></span>')
        parts.append(f'<div>{swatch}{entry.label}{_format_counts(entry.counts)}</div>')
    parts.append('</div>')
    return ''.join(parts)


MOLVIEWSPEC_MISSING_MESSAGE: str = (
    "molviewspec is not installed, so the 3D views are unavailable. "
    "Install it with `pip install molviewspec`; every other section of this "
    "module works without it."
)
"""What to print when `molviewspec_available()` is `False`. The 3D section is
optional by design -- the dependency policy keeps it out of the module's
top-level imports -- so a caller reports this and carries on."""


def molviewspec_available() -> bool:
    """
    Whether the optional `molviewspec` dependency can be imported.

    The capability check a caller uses to decide whether to draw the 3D section.
    It imports the package (and so pays for it once) but never raises.

    Returns:
        `True` if `import molviewspec` succeeds.

    Example
    -------
    >>> isinstance(molviewspec_available(), bool)
    True
    """
    try:
        import molviewspec  # noqa: F401  (lazy by policy; see the module docstring)
    except Exception:
        return False
    return True


def _require_molviewspec() -> Any:
    """
    Import and return the `molviewspec` module, or raise with a usable message.

    Every builder calls this on entry. The import is deliberately *inside* the
    function: the module must stay importable, and every non-3D function must
    stay usable, on a machine where `molviewspec` is absent.

    Returns:
        The imported `molviewspec` module.

    Raises:
        ImportError: With `MOLVIEWSPEC_MISSING_MESSAGE`.
    """
    try:
        import molviewspec as mvs
    except Exception as exc:  # pragma: no cover - exercised only without the package
        raise ImportError(MOLVIEWSPEC_MISSING_MESSAGE) from exc
    return mvs


# -- structure source -------------------------------------------------------

_STRUCTURE_URL_FIELDS: Tuple[Tuple[str, str], ...] = (
    ("bcifUrl", "bcif"),
    ("cifUrl", "mmcif"),
)
"""`(metadata field, MolViewSpec parse format)` in preference order. BinaryCIF
first: it is several times smaller over the wire and Mol* parses it natively."""


@dataclass(frozen=True)
class StructureSource:
    """
    A URL for Mol* to download, together with the format it should be parsed as.

    The two travel together because deriving them separately is a bug (R075).
    The notebook picks the URL with `meta.get('bcifUrl', meta.get('cifUrl', ...))`
    -- which falls through only on a *missing* key -- and the format with
    `'bcif' if 'bcifUrl' in meta else 'mmcif'` -- which tests only for presence.
    The `USE_LOCAL_FILE` path builds a metadata dict whose `bcifUrl` is present
    and empty, so the two expressions disagree: the URL is `''` while the format
    says `bcif`. `resolve_structure_source` makes one decision instead of two.

    Attributes:
        url:    Absolute URL Mol* downloads.
        format: `'bcif'` or `'mmcif'`, matching `url`.
    """

    url: str
    format: str


def resolve_structure_source(
    metadata: "AFDBPrediction | Mapping[str, Any]",
    prefer_binary: bool = True,
) -> StructureSource:
    """
    Resolve the structure URL and its parse format together, from one decision.

    A field that is present but empty is treated as absent, which is the whole
    point: the `USE_LOCAL_FILE` metadata dict carries `bcifUrl: ''`.

    Args:
        metadata:      An `AFDBPrediction`, or a plain metadata mapping such as
                       the notebook's `meta` dict.
        prefer_binary: `True` to prefer BinaryCIF over mmCIF. Set `False` only to
                       debug against a human-readable file.

    Returns:
        A `StructureSource` whose `url` is non-empty and whose `format` matches it.

    Raises:
        ValueError: If neither field carries a usable URL, which is the offline
            `USE_LOCAL_FILE` case. Mol* downloads the structure itself and cannot
            be handed the already-parsed text, so the 3D section genuinely cannot
            run and says so rather than emitting a viewer that silently shows
            nothing.

    Example
    -------
    >>> resolve_structure_source({'bcifUrl': 'https://x/y.bcif'})
    StructureSource(url='https://x/y.bcif', format='bcif')
    >>> resolve_structure_source({'bcifUrl': '', 'cifUrl': 'https://x/y.cif'})
    StructureSource(url='https://x/y.cif', format='mmcif')
    >>> resolve_structure_source({'bcifUrl': '', 'cifUrl': ''})
    Traceback (most recent call last):
        ...
    ValueError: No usable structure URL: 'bcifUrl' and 'cifUrl' are both absent or empty.
    """
    fields = _STRUCTURE_URL_FIELDS if prefer_binary else tuple(reversed(_STRUCTURE_URL_FIELDS))

    for field, fmt in fields:
        if isinstance(metadata, AFDBPrediction):
            try:
                url = metadata.document_url(field)
            except (KeyError, ValueError):
                continue
        else:
            url = str(metadata.get(field) or "")
        if url:
            return StructureSource(url=url, format=fmt)

    names = " and ".join(repr(field) for field, _ in _STRUCTURE_URL_FIELDS)
    raise ValueError(f"No usable structure URL: {names} are both absent or empty.")


# -- per-residue colouring --------------------------------------------------

@dataclass(frozen=True)
class ColourRun:
    """
    A contiguous stretch of residues that all take the same colour.

    One `ColourRun` becomes one `ComponentExpression` with
    `beg_label_seq_id=beg, end_label_seq_id=end`. A single residue is a run with
    `beg == end`, which is what the notebook emitted for *every* residue.

    Attributes:
        beg:    First `label_seq_id` of the run, inclusive.
        end:    Last `label_seq_id` of the run, inclusive.
        colour: Hex colour applied to the whole run.
    """

    beg: int
    end: int
    colour: str

    @property
    def n_residues(self) -> int:
        """Residues covered by this run."""
        return self.end - self.beg + 1


def colour_runs(
    res_ids: "Sequence[int] | np.ndarray",
    colours: Sequence[str],
) -> List[ColourRun]:
    """
    Collapse a per-residue colour list into contiguous same-colour runs (R075).

    A residue joins the previous run only when it is the same colour **and** its
    `label_seq_id` is exactly one more than the previous residue's. Requiring
    both is what makes the collapse safe: a ranged `ComponentExpression` covers
    every residue between its endpoints, so merging across a gap in the numbering
    would colour residues that were never in the input.

    Args:
        res_ids: `label_seq_id` per residue, ascending. Positional indices are
                 *not* accepted; Mol* selects on the file's own numbering.
        colours: One colour per entry of `res_ids`.

    Returns:
        The runs, in input order.

    Raises:
        ValueError: If the two sequences have different lengths.

    Example
    -------
    >>> colour_runs([1, 2, 3], ['#FF0000', '#FF0000', '#00FF00'])
    [ColourRun(beg=1, end=2, colour='#FF0000'), ColourRun(beg=3, end=3, colour='#00FF00')]
    >>> colour_runs([1, 2, 5, 6], ['#FF0000'] * 4)          # 3 and 4 are absent
    [ColourRun(beg=1, end=2, colour='#FF0000'), ColourRun(beg=5, end=6, colour='#FF0000')]
    >>> [run.n_residues for run in colour_runs(range(1, 11), ['#FF0000'] * 10)]
    [10]
    """
    ids = [int(r) for r in res_ids]
    cols = list(colours)
    if len(ids) != len(cols):
        raise ValueError(
            f"colour_runs got {len(ids)} residue ids but {len(cols)} colours; "
            "they must be parallel."
        )

    runs: List[ColourRun] = []
    for res_id, colour in zip(ids, cols):
        if runs and runs[-1].colour == colour and res_id == runs[-1].end + 1:
            runs[-1] = ColourRun(beg=runs[-1].beg, end=res_id, colour=colour)
        else:
            runs.append(ColourRun(beg=res_id, end=res_id, colour=colour))
    return runs


def _rgba_to_hex(rgba: Sequence[float]) -> str:
    """
    An RGBA tuple from a matplotlib colormap as an upper-case `#RRGGBB` string.

    Truncating rather than rounding, and dropping alpha, reproduces the
    notebook's own conversion exactly; `matplotlib.colors.to_hex` rounds instead
    and so differs by one unit on some channels.

    Args:
        rgba: `(r, g, b, a)` floats in `[0, 1]`; alpha is ignored.

    Returns:
        `'#RRGGBB'`.

    Example
    -------
    >>> _rgba_to_hex((1.0, 0.0, 0.0, 1.0))
    '#FF0000'
    >>> _rgba_to_hex((0.5, 0.5, 0.5, 1.0))
    '#7F7F7F'
    """
    return '#{:02X}{:02X}{:02X}'.format(
        int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255))


def value_colours(
    values: "Sequence[float] | np.ndarray",
    cmap: str | Colormap = MVS_VALUE_CMAP,
    vmin: float = 0.0,
    vmax: float = 1.0,
) -> List[str]:
    """
    Map an arbitrary per-residue value array through a colormap to hex colours.

    This is the seam R072 and R074 use. Any per-residue quantity -- ipSAE
    `d0res` today, pDockQ2's per-residue mean `ptm` next -- becomes a list of
    colours here, and `add_residue_colours` turns that into components. Neither
    function knows which score it is painting.

    `matplotlib.cm.get_cmap`, which the notebook calls, is deprecated and slated
    for removal; the lookup here goes through `matplotlib.colormaps` (R075).

    Args:
        values: Per-residue numbers.
        cmap:   A `Colormap` or a registered colormap name.
        vmin:   Value mapped to the low end of the colormap.
        vmax:   Value mapped to the high end. Values outside `[vmin, vmax]` are
                clipped, as a matplotlib `Normalize` would clip them.

    Returns:
        One `'#RRGGBB'` string per value.

    Raises:
        ValueError: If `vmax <= vmin`, which has no meaningful normalisation.

    Example
    -------
    >>> value_colours([0.0, 1.0], cmap='RdYlGn')
    ['#A50026', '#006837']
    >>> value_colours([-3.0, 4.0], cmap='RdYlGn')          # clipped to 0 and 1
    ['#A50026', '#006837']
    >>> value_colours([0.0, 50.0, 100.0], cmap='RdYlGn', vmin=0.0, vmax=100.0)
    ['#A50026', '#FEFEBD', '#006837']
    """
    if vmax <= vmin:
        raise ValueError(f"value_colours needs vmax > vmin, got vmin={vmin}, vmax={vmax}.")
    colormap = matplotlib.colormaps[cmap] if isinstance(cmap, str) else cmap
    scaled = np.clip((np.asarray(values, dtype=float) - vmin) / (vmax - vmin), 0.0, 1.0)
    return [_rgba_to_hex(colormap(float(value))) for value in scaled]


def add_residue_colours(
    structure: Any,
    chain_id: str,
    res_ids: "Sequence[int] | np.ndarray",
    colours: Sequence[str],
    representation: str = "cartoon",
) -> List[ColourRun]:
    """
    Colour named residues of one chain, one component per contiguous run (R075).

    The only place in this module where a residue becomes a MolViewSpec
    component. Residues absent from `res_ids` get no component at all and so keep
    whatever the surrounding representation gave them.

    Args:
        structure:      A MolViewSpec structure node.
        chain_id:       `label_asym_id` of the chain being coloured.
        res_ids:        `label_seq_id` per residue to colour, ascending.
        colours:        One colour per entry of `res_ids`.
        representation: `'cartoon'`, `'ball_and_stick'`, `'surface'`, ...

    Returns:
        The `ColourRun`s that were emitted, so a caller can report the component
        count without walking the state tree.

    Raises:
        ImportError: If `molviewspec` is not installed.
        ValueError:  If `res_ids` and `colours` differ in length.
    """
    mvs = _require_molviewspec()
    runs = colour_runs(res_ids, colours)
    for run in runs:
        (structure
         .component(selector=mvs.ComponentExpression(
             label_asym_id=chain_id,
             beg_label_seq_id=run.beg,
             end_label_seq_id=run.end))
         .representation(type=representation)
         .color(color=run.colour))
    return runs


def _ordered_chain_ids(contacts: InterfaceContacts) -> Tuple[str, str]:
    """The ordered pair's chain ids. The one place chain identity is resolved."""
    return contacts.chain_x, contacts.chain_y


def _chain_res_ids(chains: Mapping[str, ChainCoords], chain_id: str) -> np.ndarray:
    """One chain's `label_seq_id` array, with a message naming what is available."""
    try:
        return chains[chain_id].res_ids
    except KeyError:
        raise KeyError(
            f"Chain {chain_id!r} is not in the parsed structure; "
            f"available: {', '.join(sorted(chains)) or '(none)'}."
        ) from None


def _new_structure(source: StructureSource) -> Tuple[Any, Any]:
    """`(builder, structure)` for one view, downloading and parsing `source`."""
    mvs = _require_molviewspec()
    builder = mvs.create_builder()
    structure = (
        builder
        .download(url=source.url)
        .parse(format=source.format)
        .model_structure()
    )
    return builder, structure


# -- the display helper -----------------------------------------------------

def _css_length(value: "str | int | float") -> str:
    """
    A CSS length from either a number of pixels or an already-written length.

    Example
    -------
    >>> _css_length(600), _css_length('100%'), _css_length('40vh')
    ('600px', '100%', '40vh')
    """
    if isinstance(value, str):
        return value
    return f"{value:g}px"


def mol_view_html(
    state: Any,
    label: str,
    width: "str | int" = MVS_VIEW_WIDTH,
    height: "str | int" = MVS_VIEW_HEIGHT,
) -> str:
    """
    The label-plus-iframe HTML `show_mol_view` displays. Returned, not shown.

    Split out from `show_mol_view` so the markup can be asserted on without a
    running kernel, and so a caller that is composing its own HTML can embed a
    viewer rather than display one.

    The iframe carries its size in an inline `style` rather than in the `width` /
    `height` attributes. Those attributes are only presentational hints, which
    any host stylesheet rule outranks; an inline style outranks the stylesheet,
    so `width:100%` survives whatever the notebook front end does to the output
    area. `display:block` drops the inline-element baseline gap under the frame.

    Args:
        state:  A MolViewSpec `State`, from any `build_*_view` function.
        label:  Caption drawn above the viewer.
        width:  CSS length, or an `int` read as pixels. Defaults to
                `MVS_VIEW_WIDTH`, i.e. the full width of the column.
        height: CSS length, or an `int` read as pixels. Keep this a pixel value:
                see `MVS_VIEW_HEIGHT`.

    Returns:
        One `<div>` caption followed by one `<iframe>`.

    Example
    -------
    >>> class _S:
    ...     def molstar_html(self): return '<b>hi</b>'
    >>> markup = mol_view_html(_S(), 'View 1')
    >>> 'style="width:100%; height:600px;' in markup
    True
    >>> 'src="data:text/html;base64,PGI+aGk8L2I+"' in markup
    True
    >>> 'width:950px' in mol_view_html(_S(), 'View 1', width=950)
    True
    """
    encoded = base64.b64encode(state.molstar_html().encode()).decode()
    return (
        f'<div style="margin:10px 0 4px; font-weight:bold;">{label}</div>'
        f'<iframe src="data:text/html;base64,{encoded}" '
        f'style="width:{_css_length(width)}; height:{_css_length(height)}; '
        f'display:block; border:0;" allowfullscreen></iframe>'
    )


def show_mol_view(
    state: Any,
    label: str,
    width: "str | int" = MVS_VIEW_WIDTH,
    height: "str | int" = MVS_VIEW_HEIGHT,
) -> None:
    """
    Render a MolViewSpec `State` inline, above a bold label.

    The viewer HTML is inlined as a base64 `data:` URI rather than written to a
    file and served, which is what makes it work in PyCharm and in Colab as well
    as in classic Jupyter: none of the three agree on how a notebook-relative
    file URL resolves, and all three render a `data:` iframe.

    The markup comes from `mol_view_html`, which is where the responsive sizing
    is explained. `IPython.display.IFrame` is deliberately not used: it writes
    the size into the `width` / `height` *attributes*, which a host stylesheet
    can override, and it offers no way to add an inline style.

    The only function in this section that touches `IPython`, so every builder
    stays usable head-lessly. `IPython` is imported lazily here for the same
    reason `molviewspec` is: the module must import outside a notebook.

    Args:
        state:  A MolViewSpec `State`, from any `build_*_view` function.
        label:  Caption drawn above the viewer.
        width:  CSS length, or an `int` read as pixels; `'100%'` by default.
        height: CSS length, or an `int` read as pixels.

    Returns:
        `None`. Displays as a side effect; this is the one function here that does.
    """
    from IPython.display import HTML, display

    display(HTML(mol_view_html(state, label, width=width, height=height)))


# -- View 1: chain overview -------------------------------------------------

def build_chain_overview_view(
    source: StructureSource,
    chains: Mapping[str, ChainCoords],
    contacts: InterfaceContacts,
    chain_colours: Optional[Mapping[str, str]] = None,
    side_chain_colours: Optional[Mapping[str, str]] = None,
) -> Any:
    """
    View 1: both chains as cartoons, with interface side chains picked out.

    The orientation shot. Each chain gets its own cartoon colour, and the
    residues that actually touch the partner chain are drawn as ball-and-stick on
    top, so the interface is visible as a patch of sticks rather than having to
    be inferred from the contact map.

    **Both chains get side chains (R071).** The notebook drew them for the first
    chain only, which on a homodimer made the interface look one-sided: the
    contact is symmetric, and half of it was missing from the picture. Each
    chain's side chains take their own colour, so a stick can be attributed to a
    chain without tracing it back to its cartoon; see `SIDE_CHAIN_COLOURS` for
    why those two colours and not the two first proposed.

    Args:
        source:             Structure URL and format, from `resolve_structure_source`.
        chains:             `{chain_id: ChainCoords}`, from `parse_structure`.
        contacts:           Interface contacts naming the ordered pair to draw.
        chain_colours:      `{chain_id: colour}` for the cartoons. Defaults to
                            `CHAIN_COLOURS` assigned in ordered-pair order.
        side_chain_colours: `{chain_id: colour}` for the interface side chains.
                            Defaults to `SIDE_CHAIN_COLOURS` assigned in
                            ordered-pair order, i.e. **both** chains. Pass a
                            one-entry mapping to draw only one.

    Returns:
        A MolViewSpec `State`. Nothing is rendered; pass it to `show_mol_view`.

    Raises:
        ImportError: If `molviewspec` is not installed.
        KeyError:    If a chain named by `contacts` is absent from `chains`.
    """
    chain_x, chain_y = _ordered_chain_ids(contacts)
    cartoon = (dict(zip((chain_x, chain_y), CHAIN_COLOURS))
               if chain_colours is None else dict(chain_colours))
    side_chains = (dict(zip((chain_x, chain_y), SIDE_CHAIN_COLOURS))
                   if side_chain_colours is None else dict(side_chain_colours))
    masks = {chain_x: contacts.mask_x, chain_y: contacts.mask_y}

    mvs = _require_molviewspec()
    builder, structure = _new_structure(source)

    for chain_id in (chain_x, chain_y):
        (structure
         .component(selector=mvs.ComponentExpression(label_asym_id=chain_id))
         .representation(type='cartoon')
         .color(color=cartoon[chain_id]))

    for chain_id, colour in side_chains.items():
        res_ids = _chain_res_ids(chains, chain_id)
        selected = [int(res_ids[i]) for i in np.where(masks[chain_id])[0]]
        add_residue_colours(structure, chain_id, selected, [colour] * len(selected),
                            representation='ball_and_stick')

    return builder.get_state()


# -- View 2: pLDDT -----------------------------------------------------------

def build_plddt_view(
    source: StructureSource,
    chains: Mapping[str, ChainCoords],
    contacts: InterfaceContacts,
    plddt: Optional[Mapping[str, np.ndarray]] = None,
    bands: Sequence[Tuple[float, str, str]] = MVS_PLDDT_BANDS,
) -> Any:
    """
    View 2: every residue of both chains coloured by its own pLDDT band.

    The control for every other view. A low-confidence interface and a
    high-confidence interface can produce the same contact map, and this is where
    the difference shows.

    **Every residue is drawn (R075).** Band lookup goes through `plddt_band`,
    the module's single ladder, whose lowest edge is `-inf`; the view used to
    carry its own half-open copy under which a residue at exactly 100.00 matched
    no band and silently vanished. `plddt_legend` reads the same table.

    Args:
        source:   Structure URL and format.
        chains:   `{chain_id: ChainCoords}`, from `parse_structure`.
        contacts: Interface contacts naming the ordered pair to draw.
        plddt:    `{chain_id: (n,) array}` overriding the values carried on
                  `ChainCoords`. Pass `PLDDTScores.for_chain(...)` to colour from
                  the pLDDT JSON document instead of the mmCIF B-factor column;
                  the two carry the same numbers.
        bands:    `(exclusive lower edge, colour, label)` per band, ordered
                  high to low, as `PLDDT_BANDS`. The seam for a band-scheme
                  change; changing it changes the legend too, by construction.

    Returns:
        A MolViewSpec `State`.

    Raises:
        ImportError: If `molviewspec` is not installed.
        ValueError:  If a supplied pLDDT array's length does not match its chain.
    """
    builder, structure = _new_structure(source)

    for chain_id in _ordered_chain_ids(contacts):
        res_ids = _chain_res_ids(chains, chain_id)
        values = np.asarray(chains[chain_id].plddt if plddt is None else plddt[chain_id],
                            dtype=float)
        if values.shape[0] != res_ids.shape[0]:
            raise ValueError(
                f"Chain {chain_id!r} has {res_ids.shape[0]} residues but "
                f"{values.shape[0]} pLDDT values."
            )

        selected = [int(r) for r in res_ids]
        selected_colours = [plddt_band(float(value), bands)[0] for value in values]

        add_residue_colours(structure, chain_id, selected, selected_colours,
                            representation='cartoon')

    return builder.get_state()


# -- View 3: an arbitrary per-residue value on the interface -----------------

def build_interface_value_view(
    source: StructureSource,
    chains: Mapping[str, ChainCoords],
    contacts: InterfaceContacts,
    values_x: "Sequence[float] | np.ndarray",
    values_y: "Sequence[float] | np.ndarray",
    cmap: str | Colormap = MVS_VALUE_CMAP,
    vmin: float = 0.0,
    vmax: float = 1.0,
    base_colour: str = MVS_CONTEXT_COLOUR,
    representation: str = "ball_and_stick",
) -> Any:
    """
    View 3: interface residues of both chains, coloured by a per-residue value.

    The whole complex is drawn as a neutral grey cartoon for context, then every
    interface residue is drawn as ball-and-stick coloured by its own number,
    which turns a single headline score back into a location along the chain.

    **This function is score-agnostic, and both value views are calls to it.**
    View 3 is this called with ipSAE `d0res` per-residue values
    (`ipsae.d0res.forward.values` and `.reverse.values`); View 5 is the same call
    with `PDockQ2Direction.mean_ptm_by_residue` (R074). Nothing here knows which
    score it is painting, so a new view is a new call, not new colouring code.

    **Non-finite values are left unpainted.** `mean_ptm_by_residue` is `np.nan`
    for a residue with no contacts, and a colormap turns `nan` into its "bad"
    colour, which is opaque black at the default settings: a residue with *no
    measurement* would have been painted the most emphatic colour in the scene.
    Such residues are dropped from the selection instead and keep the base
    cartoon, so "not painted" means "no value" in both value views.

    Args:
        source:         Structure URL and format.
        chains:         `{chain_id: ChainCoords}`, from `parse_structure`.
        contacts:       Interface contacts naming the ordered pair; its two masks
                        choose which residues are drawn.
        values_x:       `(nx,)` per-residue value for `contacts.chain_x`. Only
                        the interface entries are read.
        values_y:       `(ny,)` per-residue value for `contacts.chain_y`.
        cmap:           Colormap. Default `MVS_VALUE_CMAP`, red-yellow-green.
        vmin, vmax:     Value range mapped onto the colormap; values outside are
                        clipped. `0..1` suits any of the ptm-derived scores.
        base_colour:    Cartoon colour for the rest of the complex.
        representation: Representation for the coloured residues.

    Returns:
        A MolViewSpec `State`.

    Raises:
        ImportError: If `molviewspec` is not installed.
        ValueError:  If a value array's length does not match its chain.
    """
    chain_x, chain_y = _ordered_chain_ids(contacts)
    builder, structure = _new_structure(source)

    (structure
     .component()
     .representation(type='cartoon')
     .color(color=base_colour))

    for chain_id, mask, values in ((chain_x, contacts.mask_x, values_x),
                                   (chain_y, contacts.mask_y, values_y)):
        res_ids = _chain_res_ids(chains, chain_id)
        value_array = np.asarray(values, dtype=float)
        if value_array.shape[0] != res_ids.shape[0]:
            raise ValueError(
                f"Chain {chain_id!r} has {res_ids.shape[0]} residues but "
                f"{value_array.shape[0]} values."
            )
        indices = np.where(np.asarray(mask, dtype=bool) & np.isfinite(value_array))[0]
        add_residue_colours(
            structure,
            chain_id,
            [int(res_ids[i]) for i in indices],
            value_colours(value_array[indices], cmap=cmap, vmin=vmin, vmax=vmax),
            representation=representation,
        )

    return builder.get_state()


# -- Views 4 and 6: categorical agreement views ------------------------------
# Two views ask "do two signals about this residue agree?", of two different
# signal pairs. They share one builder and one legend path, so the second view
# is a classifier plus a colour table rather than a second copy of the drawing
# code -- the same discipline `build_interface_value_view` applies to the two
# continuous views.


def build_category_view(
    source: StructureSource,
    chains: Mapping[str, ChainCoords],
    contacts: InterfaceContacts,
    keys_x: Sequence[str],
    keys_y: Sequence[str],
    categories: Sequence[Tuple[str, str, str]],
    base_colour: str = MVS_FAINT_COLOUR,
    representation: str = "ball_and_stick",
) -> Any:
    """
    Paint both chains from a per-residue category key, over a faint cartoon.

    The shared engine behind View 4 and View 6. A residue whose key is not in
    `categories` -- including the empty string the classifiers use for "no
    category applies" -- gets no component and keeps the base cartoon.

    Args:
        source:         Structure URL and format.
        chains:         `{chain_id: ChainCoords}`, from `parse_structure`.
        contacts:       Interface contacts naming the ordered pair.
        keys_x:         `(nx,)` category key per residue of `contacts.chain_x`.
        keys_y:         `(ny,)` category key per residue of `contacts.chain_y`.
        categories:     `(key, colour, label)` triples; only `key` and `colour`
                        are read here.
        base_colour:    Cartoon colour for the whole complex.
        representation: Representation for the categorised residues.

    Returns:
        A MolViewSpec `State`.

    Raises:
        ImportError: If `molviewspec` is not installed.
        ValueError:  If a key array's length does not match its chain.
    """
    _require_molviewspec()
    chain_x, chain_y = _ordered_chain_ids(contacts)
    builder, structure = _new_structure(source)

    (structure
     .component()
     .representation(type='cartoon')
     .color(color=base_colour))

    colour_of = {key: colour for key, colour, _label in categories}

    for chain_id, keys in ((chain_x, keys_x), (chain_y, keys_y)):
        res_ids = _chain_res_ids(chains, chain_id)
        key_array = np.asarray(keys, dtype=object)
        if key_array.shape[0] != res_ids.shape[0]:
            raise ValueError(
                f"Chain {chain_id!r} has {res_ids.shape[0]} residues but "
                f"{key_array.shape[0]} category keys."
            )
        selected: List[int] = []
        selected_colours: List[str] = []
        for index, key in enumerate(key_array):
            if key in colour_of:
                selected.append(int(res_ids[index]))
                selected_colours.append(colour_of[key])
        add_residue_colours(structure, chain_id, selected, selected_colours,
                            representation=representation)

    return builder.get_state()


# -- View 4: ipSAE confidence against physical contact ------------------------

def disagreement_categories(
    values: "Sequence[float] | np.ndarray",
    interface_mask: np.ndarray,
    threshold: float = MVS_DISAGREEMENT_THRESHOLD,
) -> np.ndarray:
    """
    Classify each residue by whether PAE confidence and physical contact agree.

    Two independent signals claim a residue is at the interface: a CB atom
    within `DIST_CUTOFF` of the partner chain, and a per-residue ipSAE_d0res at
    or above `threshold`. Three of the four combinations are categories; the
    fourth -- neither signal -- is left uncategorised and undrawn, because
    painting "nothing to report" over most of a structure hides the report.

    **The cutoff is `THRESHOLDS['ipsae_d0res'].amber`, not a magic 0.5 (R073),**
    and the test is `>=`, matching `traffic_light`. The notebook used `> 0.5`,
    a number with no relation to any published or AFDB cutoff, so a residue
    could be "confident" here and red in the Section 7 summary.

    Args:
        values:         `(n,)` per-residue ipSAE_d0res, `0..1`.
        interface_mask: `(n,)` bool, `True` where the residue touches the partner.
        threshold:      Score at or above which PAE is called confident.

    Returns:
        `(n,)` array of `MVS_DISAGREEMENT_CATEGORIES` keys, `''` where none applies.

    Raises:
        ValueError: If the two arrays have different lengths.

    Example
    -------
    >>> import numpy as np
    >>> mask = np.array([True, False, True, False])
    >>> disagreement_categories([0.9, 0.9, 0.1, 0.1], mask).tolist()
    ['confirmed', 'predicted_not_touching', 'touching_not_trusted', '']

    The edge is inclusive, and it is the one `THRESHOLDS` publishes:

    >>> MVS_DISAGREEMENT_THRESHOLD == THRESHOLDS['ipsae_d0res'].amber == 0.6
    True
    >>> disagreement_categories([0.6, 0.599999], np.array([True, True])).tolist()
    ['confirmed', 'touching_not_trusted']
    """
    value_array = np.asarray(values, dtype=float)
    mask = np.asarray(interface_mask, dtype=bool)
    if value_array.shape[0] != mask.shape[0]:
        raise ValueError(
            f"disagreement_categories got {value_array.shape[0]} values but "
            f"{mask.shape[0]} mask entries."
        )
    confident = value_array >= threshold
    categories = np.full(value_array.shape[0], '', dtype=object)
    categories[confident & mask] = 'confirmed'
    categories[confident & ~mask] = 'predicted_not_touching'
    categories[~confident & mask] = 'touching_not_trusted'
    return categories


def build_disagreement_view(
    source: StructureSource,
    chains: Mapping[str, ChainCoords],
    contacts: InterfaceContacts,
    values_x: "Sequence[float] | np.ndarray",
    values_y: "Sequence[float] | np.ndarray",
    threshold: float = MVS_DISAGREEMENT_THRESHOLD,
    categories: Sequence[Tuple[str, str, str]] = MVS_DISAGREEMENT_CATEGORIES,
    base_colour: str = MVS_FAINT_COLOUR,
    representation: str = "ball_and_stick",
) -> Any:
    """
    View 4: where ipSAE confidence and physical contact disagree, on both chains.

    Agreement is the common case; the disagreements are what the notebook exists
    to explain. Pair with `disagreement_legend`, which states the cutoff in force
    and the population of every category.

    **Both chains are drawn (R073).** The notebook coloured `chain_x` only, so on
    a homodimer the view showed one copy of a symmetric answer and on a
    heterodimer it silently omitted a whole protein.

    Args:
        source:         Structure URL and format.
        chains:         `{chain_id: ChainCoords}`, from `parse_structure`.
        contacts:       Interface contacts naming the ordered pair.
        values_x:       `(nx,)` per-residue ipSAE_d0res for `contacts.chain_x`,
                        i.e. `ipsae.d0res.forward.values`.
        values_y:       `(ny,)` for `contacts.chain_y`, i.e.
                        `ipsae.d0res.reverse.values`.
        threshold:      Confidence cutoff, defaulting to the `THRESHOLDS` value.
        categories:     `(key, colour, label template)` per category.
        base_colour:    Cartoon colour for the rest of the complex.
        representation: Representation for the categorised residues.

    Returns:
        A MolViewSpec `State`.

    Raises:
        ImportError: If `molviewspec` is not installed.
        ValueError:  If a value array does not match its chain's length.
    """
    return build_category_view(
        source, chains, contacts,
        disagreement_categories(values_x, contacts.mask_x, threshold=threshold),
        disagreement_categories(values_y, contacts.mask_y, threshold=threshold),
        categories=categories,
        base_colour=base_colour,
        representation=representation,
    )


def disagreement_legend(
    contacts: InterfaceContacts,
    values_x: "Sequence[float] | np.ndarray",
    values_y: "Sequence[float] | np.ndarray",
    label_x: str,
    label_y: str,
    threshold: float = MVS_DISAGREEMENT_THRESHOLD,
    dist_cutoff: float = DIST_CUTOFF,
    categories: Sequence[Tuple[str, str, str]] = MVS_DISAGREEMENT_CATEGORIES,
) -> List[LegendEntry]:
    """
    View 4's legend: the categories, the cutoffs in force, and the populations.

    The legend the notebook computed and threw away (R073). It is built from the
    same classifier call the view is painted from, so a row cannot claim a colour
    the view does not use or a count the view does not draw.

    Args:
        contacts:         The ordered pair drawn.
        values_x:         Per-residue ipSAE_d0res for `contacts.chain_x`.
        values_y:         Per-residue ipSAE_d0res for `contacts.chain_y`.
        label_x, label_y: Display names for the two chains.
        threshold:        The cutoff in force; appears in every label.
        dist_cutoff:      The contact cutoff; appears in every label.
        categories:       `(key, colour, label template)` per category.

    Returns:
        One `LegendEntry` per category, then the "not drawn" row.
    """
    return category_legend(
        categories,
        {label_x: disagreement_categories(values_x, contacts.mask_x, threshold=threshold),
         label_y: disagreement_categories(values_y, contacts.mask_y, threshold=threshold)},
        threshold=threshold,
        dist_cutoff=dist_cutoff,
        unlit_label=MVS_UNLIT_LABELS['disagreement'],
    )


# -- Views 5 and 6: pDockQ2, per residue (R074) -------------------------------

def pdockq2_ipsae_categories(
    mean_ptm: "Sequence[float] | np.ndarray",
    ipsae_values: "Sequence[float] | np.ndarray",
    interface_mask: np.ndarray,
    contact_ptm_threshold: float = MVS_CONTACT_PTM_THRESHOLD,
    ipsae_threshold: float = MVS_DISAGREEMENT_THRESHOLD,
) -> np.ndarray:
    """
    Classify each interface residue by whether pDockQ2 and ipSAE agree about it.

    The two scores read the same PAE matrix at different scopes. `ipSAE_d0res`
    scores a residue against **every** partner residue whose PAE clears
    `PAE_CUTOFF`; `PDockQ2Direction.mean_ptm_by_residue` scores it against
    **only the partner residues it physically touches**. Where those disagree is
    where the two headline scores disagree, and this puts that on the structure.

    Non-interface residues have no pDockQ2 value (`np.nan`) and are left
    uncategorised: there is nothing to compare, which is a different statement
    from "the comparison came out badly".

    Args:
        mean_ptm:              `(n,)` mean `ptm(PAE, d0=10)` over each residue's
                               own contacts; `np.nan` where there are none.
        ipsae_values:          `(n,)` per-residue ipSAE_d0res for the same chain
                               and the same direction.
        interface_mask:        `(n,)` bool, `True` where the residue has contacts.
        contact_ptm_threshold: Contact quality at or above which the contacts are
                               called well placed. Defaults to
                               `MVS_CONTACT_PTM_THRESHOLD`, i.e. mean contact PAE
                               better than `PAE_CUTOFF`.
        ipsae_threshold:       ipSAE_d0res cutoff, from `THRESHOLDS`.

    Returns:
        `(n,)` array of `MVS_PDOCKQ2_AGREEMENT_CATEGORIES` keys, `''` off the
        interface.

    Raises:
        ValueError: If the three arrays have different lengths.

    Example
    -------
    >>> import numpy as np
    >>> mask = np.array([True, True, True, True, False])
    >>> pdockq2_ipsae_categories([0.9, 0.9, 0.1, 0.1, np.nan],
    ...                          [0.9, 0.1, 0.9, 0.1, 0.9], mask).tolist()
    ['both', 'contacts_only', 'ipsae_only', 'neither', '']
    """
    ptm_array = np.asarray(mean_ptm, dtype=float)
    ipsae_array = np.asarray(ipsae_values, dtype=float)
    mask = np.asarray(interface_mask, dtype=bool)
    if not (ptm_array.shape[0] == ipsae_array.shape[0] == mask.shape[0]):
        raise ValueError(
            f"pdockq2_ipsae_categories got {ptm_array.shape[0]} ptm values, "
            f"{ipsae_array.shape[0]} ipSAE values and {mask.shape[0]} mask "
            "entries; all three must be parallel."
        )
    scored = mask & np.isfinite(ptm_array)
    good_contacts = scored & (ptm_array >= contact_ptm_threshold)
    good_ipsae = ipsae_array >= ipsae_threshold

    categories = np.full(ptm_array.shape[0], '', dtype=object)
    categories[scored & good_contacts & good_ipsae] = 'both'
    categories[scored & good_contacts & ~good_ipsae] = 'contacts_only'
    categories[scored & ~good_contacts & good_ipsae] = 'ipsae_only'
    categories[scored & ~good_contacts & ~good_ipsae] = 'neither'
    return categories


def build_pdockq2_agreement_view(
    source: StructureSource,
    chains: Mapping[str, ChainCoords],
    contacts: InterfaceContacts,
    mean_ptm_x: "Sequence[float] | np.ndarray",
    mean_ptm_y: "Sequence[float] | np.ndarray",
    ipsae_x: "Sequence[float] | np.ndarray",
    ipsae_y: "Sequence[float] | np.ndarray",
    contact_ptm_threshold: float = MVS_CONTACT_PTM_THRESHOLD,
    ipsae_threshold: float = MVS_DISAGREEMENT_THRESHOLD,
    categories: Sequence[Tuple[str, str, str]] = MVS_PDOCKQ2_AGREEMENT_CATEGORIES,
    base_colour: str = MVS_FAINT_COLOUR,
    representation: str = "ball_and_stick",
) -> Any:
    """
    View 6: where pDockQ2's contact quality and ipSAE's chain-wide confidence
    reach different conclusions about the same residue (R074).

    Args:
        source:                Structure URL and format.
        chains:                `{chain_id: ChainCoords}`.
        contacts:              Interface contacts naming the ordered pair.
        mean_ptm_x:            `(nx,)` `pdockq2.forward.mean_ptm_by_residue`.
        mean_ptm_y:            `(ny,)` `pdockq2.reverse.mean_ptm_by_residue`.
        ipsae_x:               `(nx,)` `ipsae.d0res.forward.values`.
        ipsae_y:               `(ny,)` `ipsae.d0res.reverse.values`.
        contact_ptm_threshold: Contact-quality cutoff.
        ipsae_threshold:       ipSAE_d0res cutoff, from `THRESHOLDS`.
        categories:            `(key, colour, label template)` per category.
        base_colour:           Cartoon colour for the rest of the complex.
        representation:        Representation for the categorised residues.

    Returns:
        A MolViewSpec `State`.

    Raises:
        ImportError: If `molviewspec` is not installed.
        ValueError:  If an array does not match its chain's length.
    """
    return build_category_view(
        source, chains, contacts,
        pdockq2_ipsae_categories(mean_ptm_x, ipsae_x, contacts.mask_x,
                                 contact_ptm_threshold, ipsae_threshold),
        pdockq2_ipsae_categories(mean_ptm_y, ipsae_y, contacts.mask_y,
                                 contact_ptm_threshold, ipsae_threshold),
        categories=categories,
        base_colour=base_colour,
        representation=representation,
    )


def pdockq2_agreement_legend(
    contacts: InterfaceContacts,
    mean_ptm_x: "Sequence[float] | np.ndarray",
    mean_ptm_y: "Sequence[float] | np.ndarray",
    ipsae_x: "Sequence[float] | np.ndarray",
    ipsae_y: "Sequence[float] | np.ndarray",
    label_x: str,
    label_y: str,
    contact_ptm_threshold: float = MVS_CONTACT_PTM_THRESHOLD,
    ipsae_threshold: float = MVS_DISAGREEMENT_THRESHOLD,
    dist_cutoff: float = DIST_CUTOFF,
    categories: Sequence[Tuple[str, str, str]] = MVS_PDOCKQ2_AGREEMENT_CATEGORIES,
) -> List[LegendEntry]:
    """
    View 6's legend: the four categories, both cutoffs, and the populations.

    Args:
        contacts:              The ordered pair drawn.
        mean_ptm_x:            `pdockq2.forward.mean_ptm_by_residue`.
        mean_ptm_y:            `pdockq2.reverse.mean_ptm_by_residue`.
        ipsae_x:               `ipsae.d0res.forward.values`.
        ipsae_y:               `ipsae.d0res.reverse.values`.
        label_x, label_y:      Display names for the two chains.
        contact_ptm_threshold: Contact-quality cutoff, stated as a PAE.
        ipsae_threshold:       ipSAE_d0res cutoff.
        dist_cutoff:           Contact cutoff, for the "not drawn" row.
        categories:            `(key, colour, label template)` per category.

    Returns:
        One `LegendEntry` per category, then the "not drawn" row.
    """
    return category_legend(
        categories,
        {label_x: pdockq2_ipsae_categories(mean_ptm_x, ipsae_x, contacts.mask_x,
                                           contact_ptm_threshold, ipsae_threshold),
         label_y: pdockq2_ipsae_categories(mean_ptm_y, ipsae_y, contacts.mask_y,
                                           contact_ptm_threshold, ipsae_threshold)},
        threshold=ipsae_threshold,
        dist_cutoff=dist_cutoff,
        contact_ptm=contact_ptm_threshold,
        unlit_label=MVS_UNLIT_LABELS['pdockq2_agreement'],
    )


# ---------------------------------------------------------------------------
# Notebook orchestration (R100)
# ---------------------------------------------------------------------------
# Everything below is plumbing the notebook used to spell out inline: the
# environment banner, the upload widgets, the download narration, and one
# `format_*` per printed block. None of it computes a score. It lives here for
# the same reason the figures do -- a reader of the notebook should see *which*
# function is called with *what* data, and not the twenty lines of `print` that
# lay the answer out.
#
# The rule that decides what moves. A cell keeps the call and the arguments,
# because that is the teaching content: `compute_pdockq2(contacts, pair,
# plddt_x, plddt_y)` says exactly which four things pDockQ2 is made of. A cell
# gives up its formatting, because "%.4f in a column of width 12" teaches
# nothing. No print was deleted in the move; every number the notebook used to
# explain, it still explains, from here.
#
# One thing deliberately did *not* move: the notebook's bootstrap cell. Its job
# is to make this module importable, so it cannot call this module to do it.
# What it can hand over is everything that happens *after* the import succeeds,
# which is `prepare_environment` below.


# -- the bootstrap's second half --------------------------------------------

def describe_checkout(root: "str | Path") -> str:
    """
    `'branch @ sha'` for a git checkout, or a plain note when it is not one.

    Args:
        root: Directory to describe.

    Returns:
        `'<branch> @ <short sha>'`, or `'unknown (not a git checkout)'` when
        `git` is absent, fails, or `root` is not a work tree.

    Example
    -------
    >>> describe_checkout('/definitely/not/a/checkout')
    'unknown (not a git checkout)'
    """
    try:
        rev = subprocess.run(['git', '-C', str(root), 'rev-parse', '--abbrev-ref', 'HEAD'],
                             capture_output=True, text=True)
        sha = subprocess.run(['git', '-C', str(root), 'rev-parse', '--short', 'HEAD'],
                             capture_output=True, text=True)
        if rev.returncode == 0 and sha.returncode == 0:
            return f'{rev.stdout.strip()} @ {sha.stdout.strip()}'
    except (OSError, subprocess.SubprocessError):
        pass
    return 'unknown (not a git checkout)'


def prepare_environment(
    repo_root: "str | Path",
    branch: str = "",
    colab: bool = False,
    install_missing: Optional[bool] = None,
) -> str:
    """
    Finish the bootstrap: ensure `molviewspec`, then say where everything came from.

    Called by the notebook's first cell immediately after the import succeeds.
    Splitting it here is what keeps that cell down to the part that genuinely
    cannot use this module: finding or creating a checkout.

    `molviewspec` powers the optional 3D views. It is installed only when it is
    genuinely missing *and* only on Colab: locally it comes from the environment,
    and a `pip` call on every run is pure latency. A failed install is not fatal
    -- Section 6 says it is unavailable and every score is unaffected.

    Args:
        repo_root:       The checkout in use, for the banner.
        branch:          Branch the notebook expects, for the banner. Blank omits
                         the line.
        colab:           Whether this is Google Colab. Gates the install.
        install_missing: Overrides that gate. `None` means "install iff `colab`".

    Returns:
        The revision string, as `describe_checkout` reports it.

    Note:
        Prints the banner as a side effect; that is the point of calling it.
    """
    if install_missing is None:
        install_missing = colab

    if not molviewspec_available() and install_missing:
        print('Installing molviewspec ...')
        subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', 'molviewspec'],
                       check=False)
        importlib.invalidate_caches()

    revision = describe_checkout(repo_root)
    print(f'Environment: {"Colab" if colab else "local"}')
    print(f'Repo root:   {repo_root}')
    print(f'Revision:    {revision}')
    if branch:
        print(f'Branch expected: {branch}')
    print('molviewspec: ' + ('available' if molviewspec_available()
                             else 'not installed (Section 6 will be skipped)'))
    return revision


# -- Section 1: where the three documents come from --------------------------

@dataclass
class UploadPanel:
    """
    The three `USE_LOCAL_FILE` upload slots, or an empty stand-in when online.

    All three slots are required, not one plus two optionals (R025): six of the
    seven values are read off the PAE matrix and the seventh, pDockQ, needs
    per-residue pLDDT, so a partial upload leaves the traffic light with nothing
    to colour. `load_documents` refuses a partial upload by name.

    Attributes:
        enabled: `True` in local-file mode. `False` leaves all three slots `None`.
        cif:     `ipywidgets.FileUpload` for the mmCIF, or `None`.
        pae:     `FileUpload` for the PAE JSON, or `None`.
        plddt:   `FileUpload` for the pLDDT JSON, or `None`.
    """

    enabled: bool
    cif: Any = None
    pae: Any = None
    plddt: Any = None

    @property
    def slots(self) -> Tuple[Tuple[str, Any], ...]:
        """`(label, widget)` for the three slots, in the order they are printed."""
        return (('mmCIF', self.cif), ('PAE  ', self.pae), ('pLDDT', self.plddt))

    @staticmethod
    def content(widget: Any) -> Optional[bytes]:
        """
        Bytes from a `FileUpload` widget, or `None` if nothing was uploaded.

        Example
        -------
        >>> UploadPanel.content(None) is None
        True
        """
        if widget is None or not widget.value:
            return None
        return bytes(widget.value[0]['content'])

    @staticmethod
    def filename(widget: Any) -> str:
        """
        The uploaded file's name, or a marker saying the slot is empty.

        Example
        -------
        >>> UploadPanel.filename(None)
        '(not uploaded)'
        """
        if widget is None or not widget.value:
            return '(not uploaded)'
        return widget.value[0]['name']


UPLOAD_INSTRUCTIONS: str = (
    '<b>Upload all three files, then run the next cell.</b><br>'
    'All three are required: the PAE matrix supplies six of the seven '
    'scores, and per-residue pLDDT supplies the seventh.<br>'
    'For an AFDB model they are the <code>cifUrl</code>, '
    '<code>paeDocUrl</code> and <code>plddtDocUrl</code> downloads: '
    '<code>&hellip;-model_v1.cif</code>, '
    '<code>&hellip;-predicted_aligned_error_v1.json</code> and '
    '<code>&hellip;-confidence_v1.json</code>.'
)
"""What sits above the three upload slots, naming the AFDB file each one wants."""


def local_upload_panel(enabled: bool) -> UploadPanel:
    """
    Build and display the three upload slots, or say that we are fetching online.

    Args:
        enabled: The notebook's `USE_LOCAL_FILE`.

    Returns:
        An `UploadPanel`. When `enabled` is `False` it holds no widgets, and
        `load_documents` will download instead.

    Note:
        Displays as a side effect. `ipywidgets` and `IPython` are imported lazily,
        so importing this module outside a notebook stays free.

    Example
    -------
    >>> panel = local_upload_panel(False)
    Online mode: files will be downloaded from AFDB.
    >>> panel.enabled, panel.cif
    (False, None)
    """
    if not enabled:
        print('Online mode: files will be downloaded from AFDB.')
        return UploadPanel(enabled=False)

    import ipywidgets as widgets
    from IPython.display import display

    panel = UploadPanel(
        enabled=True,
        cif=widgets.FileUpload(accept='.cif,.mmcif', multiple=False,
                               description='mmCIF (required)'),
        pae=widgets.FileUpload(accept='.json', multiple=False,
                               description='PAE JSON (required)'),
        plddt=widgets.FileUpload(accept='.json', multiple=False,
                                 description='pLDDT JSON (required)'),
    )
    display(widgets.VBox([widgets.HTML(UPLOAD_INSTRUCTIONS),
                          panel.cif, panel.pae, panel.plddt]))
    return panel


LOCAL_MODE_GATE_NOTE: str = (
    'No metadata to check, so the declared-assembly gate is skipped here. '
    'The\nstructural gate inside verify_chain_identity still runs, and is '
    'what refuses a\nmonomer or a model with more than two chains.'
)
"""Printed in local-file mode in place of the declared-assembly gate.

With no metadata there is nothing for that gate to read, so local-file mode
passes only one of the two. The one it passes is the structural gate inside
`verify_chain_identity`, which is the authoritative one: it reads the chains
themselves, and a monomer or a three-chain upload is refused there either way."""


def fetch_prediction(
    accession: str,
    online: bool = True,
    timeout: float = DEFAULT_TIMEOUT,
) -> Optional[AFDBPrediction]:
    """
    Fetch the AFDB metadata for `accession`, narrating what it did.

    The URL is printed *before* the request, so a refused accession shows what
    was asked for as well as why it failed.

    Args:
        accession: AFDB accession, e.g. `'AF-0000000065889468'`.
        online:    `False` (the notebook's `USE_LOCAL_FILE`) skips the API
                   entirely and returns `None`.
        timeout:   Per-request timeout, in seconds.

    Returns:
        An `AFDBPrediction`, or `None` in local-file mode.

    Raises:
        AccessionLookupError: Propagated from `fetch_afdb_metadata`.

    Example
    -------
    >>> fetch_prediction('AF-x', online=False) is None
    Local file mode: skipping AFDB API for AF-x
    No metadata to check, so the declared-assembly gate is skipped here. The
    structural gate inside verify_chain_identity still runs, and is what refuses a
    monomer or a model with more than two chains.
    True
    """
    if not online:
        print(f'Local file mode: skipping AFDB API for {accession}')
        print(LOCAL_MODE_GATE_NOTE)
        return None

    print(f'Fetching: {AFDB_PREDICTION_URL.format(accession=accession)}')
    prediction = fetch_afdb_metadata(accession, timeout=timeout)
    # One entry per chain, and the endpoint's order is non-deterministic: the
    # same accession answers ['A', 'B'] on one call and ['B', 'A'] on the next.
    # `AFDBPrediction` sorts the entries by chain id at construction and every
    # field is read by chain id, so nothing indexes an entry by position.
    print(f'Chains described: {list(prediction.chain_ids)}')
    print('Available fields:', sorted({field_name
                                       for chain in prediction.chain_ids
                                       for field_name in prediction.entry_for_chain(chain)}))
    return prediction


class SourceDocuments(NamedTuple):
    """
    The three documents every later cell reads, however they were obtained.

    Attributes:
        cif_text: mmCIF text, for `parse_structure`.
        pae:      Parsed PAE JSON, for `parse_pae`.
        plddt:    Parsed pLDDT JSON, for `parse_plddt`.
    """

    cif_text: str
    pae: Any
    plddt: Any


def load_documents(
    prediction: Optional[AFDBPrediction],
    uploads: Optional[UploadPanel] = None,
    accession: str = "",
    timeout: float = DEFAULT_TIMEOUT,
) -> SourceDocuments:
    """
    Get the mmCIF, PAE and pLDDT documents, by download or from the upload slots.

    Which path is taken is decided by `uploads.enabled`, i.e. by the notebook's
    `USE_LOCAL_FILE`, not by whether `prediction` happens to be `None`.

    Each download URL is printed before its request, so a partial failure names
    the document that failed rather than only the exception.

    Args:
        prediction: Fetched metadata. Required unless `uploads.enabled`.
        uploads:    The `UploadPanel` from `local_upload_panel`.
        accession:  Accession, used only in the missing-document message.
        timeout:    Per-request timeout, in seconds.

    Returns:
        A `SourceDocuments`, unpackable as `cif_text, pae_raw, plddt_raw`.

    Raises:
        MissingLocalDocumentError: If local-file mode is on and any of the three
            slots is empty. Every missing file is named at once, while the upload
            widget is still on screen.
        ValueError: If online mode was asked for without metadata to download from.
    """
    if uploads is not None and uploads.enabled:
        cif_bytes = UploadPanel.content(uploads.cif)
        pae_bytes = UploadPanel.content(uploads.pae)
        plddt_bytes = UploadPanel.content(uploads.plddt)
        for label, widget in uploads.slots:
            print(f'{label} : {UploadPanel.filename(widget)}')
        # R025. This used to print "PAE file not uploaded -- PAE-dependent
        # analyses will be skipped" and then skip nothing. There is no useful
        # partial run to skip *to*, so the honest answer is to refuse here.
        require_local_documents(cif_bytes, pae_bytes, plddt_bytes, accession=accession)
        documents = SourceDocuments(
            cif_text=cif_bytes.decode('utf-8', errors='replace'),
            pae=json.loads(pae_bytes),
            plddt=json.loads(plddt_bytes),
        )
        print('All three local files loaded.')
        return documents

    if prediction is None:
        raise ValueError(
            'load_documents needs either fetched metadata or an enabled '
            'UploadPanel; it was given neither.'
        )

    print(f'Downloading mmCIF: {prediction.cif_url}')
    cif_text = download_structure(prediction, timeout=timeout)
    print(f'Downloading PAE:   {prediction.pae_url}')
    pae_raw = download_pae(prediction, timeout=timeout)
    print(f'Downloading pLDDT: {prediction.plddt_url}')
    plddt_raw = download_plddt(prediction, timeout=timeout)
    print('All downloads complete.')
    return SourceDocuments(cif_text=cif_text, pae=pae_raw, plddt=plddt_raw)


# -- Section 1: what the documents turned out to contain ---------------------

def format_chain_lengths(chains: Mapping[str, ChainCoords]) -> str:
    """
    The chains the mmCIF parse found, and how long each one is.

    Args:
        chains: `parse_structure` output.

    Returns:
        One header line plus one indented line per chain, sorted by chain id.

    Example
    -------
    >>> coords = np.zeros((3, 3))
    >>> ids = np.array([1, 2, 3])
    >>> names = np.array(['ALA', 'GLY', 'SER'])
    >>> chains = {'B': ChainCoords('B', coords, ids, names, None),
    ...           'A': ChainCoords('A', coords[:2], ids[:2], names[:2], None)}
    >>> print(format_chain_lengths(chains))
    Chains found: ['A', 'B']
      Chain A: 2 residues
      Chain B: 3 residues
    """
    chain_ids = sorted(chains)
    lines = [f'Chains found: {chain_ids}']
    lines += [f'  Chain {cid}: {chains[cid].n_residues} residues' for cid in chain_ids]
    return '\n'.join(lines)


def format_chain_report(
    identity: ChainIdentity,
    pair: ChainPairPAE,
    plddt_x: np.ndarray,
    plddt_y: np.ndarray,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
    pae_shape: "Optional[Tuple[int, ...]]" = None,
) -> str:
    """
    What the three documents agreed the chains are, and how the PAE splits up.

    Three sources describe the chains -- the mmCIF, the PAE document and the
    pLDDT document -- and every quadrant slice assumes all three agree. This is
    the report of the check that established it, followed by the four quadrant
    shapes that check makes safe to take.

    Args:
        identity:  `verify_chain_identity` result.
        pair:      The ordered chain pair the notebook goes on to score.
        plddt_x:   `(nx,)` pLDDT for `pair.chain_x`.
        plddt_y:   `(ny,)` pLDDT for `pair.chain_y`.
        label_x:   Display name for `chain_x`.
        label_y:   Display name for `chain_y`.
        pae_shape: Shape of the full PAE matrix. Omit it to derive the square
                   `(nx + ny, nx + ny)` from the pair.

    Returns:
        The multi-line report, ending with the two mean pLDDT values.
    """
    name_x = _chain_label(pair.chain_x, label_x)
    name_y = _chain_label(pair.chain_y, label_y)
    total = pair.nx + pair.ny
    shape = tuple(pae_shape) if pae_shape is not None else (total, total)

    lines = [identity.assembly.headline, '',
             'Chains verified across structure, PAE and pLDDT:',
             identity.legend()]
    lines += [f'NOTE: {note}' for note in identity.notes]
    lines += ['',
              f'PAE matrix shape: {shape}',
              'PAE quadrants:',
              f'  intra {name_x}: {pair.block_xx.shape}',
              f'  inter {name_x} → {name_y}: {pair.block_xy.shape}',
              f'  inter {name_y} → {name_x}: {pair.block_yx.shape}',
              f'  intra {name_y}: {pair.block_yy.shape}',
              f'pLDDT mean: {name_x} {np.asarray(plddt_x).mean():.1f}, '
              f'{name_y} {np.asarray(plddt_y).mean():.1f}']
    return '\n'.join(lines)


# -- Section 2: the interface ------------------------------------------------

def format_interface_report(
    contacts: InterfaceContacts,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
    res_ids_x: Optional[np.ndarray] = None,
    res_ids_y: Optional[np.ndarray] = None,
) -> str:
    """
    Contact cutoff, contact count, and how much of each chain is at the interface.

    Both chains, not just the first: for a heterodimer the two coverage figures
    are genuinely different numbers, and reporting one of them was only ever
    harmless while the chains were copies of each other.

    Args:
        contacts:  `detect_interface` result.
        label_x:   Display name for `chain_x`.
        label_y:   Display name for `chain_y`.
        res_ids_x: Optional `(nx,)` residue numbers for `chain_x`. Supplying them
                   adds the interface residue range, in the model's own numbering
                   rather than in positional indices.
        res_ids_y: Optional `(ny,)` residue numbers for `chain_y`.

    Returns:
        Four lines, plus one range line per chain that has residue numbers and at
        least one interface residue.
    """
    name_x = _chain_label(contacts.chain_x, label_x)
    name_y = _chain_label(contacts.chain_y, label_y)
    nx, ny = contacts.contact_mask.shape

    lines = [f'Contact cutoff         : {contacts.dist_cutoff} Å (CB-CB; CA for GLY)',
             f'Number of contact pairs: {contacts.n_contact_pairs}']
    for name, n_if, n_total in ((name_x, contacts.n_interface_residues_x, nx),
                                (name_y, contacts.n_interface_residues_y, ny)):
        lines.append(f'Interface residues, {name}: {n_if} / {n_total} '
                     f'({100 * n_if / n_total:.1f}%)')

    for name, res_ids, mask in ((name_x, res_ids_x, contacts.mask_x),
                                (name_y, res_ids_y, contacts.mask_y)):
        if res_ids is None:
            continue
        if_res = np.asarray(res_ids)[mask]
        if if_res.size:
            lines.append(f'Interface residue range, {name}: {if_res[0]} – {if_res[-1]}')
    return '\n'.join(lines)


# -- Section 3: which PAE cells each score reads -----------------------------

def format_score_mask_usage(
    pair: ChainPairPAE,
    contacts: InterfaceContacts,
    pae_cutoff: float = PAE_CUTOFF,
    lis_cutoff: float = LIS_CUTOFF,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
) -> str:
    """
    How many cells of the inter-chain PAE block each score actually reads.

    The same `score_masks` the four panels of `plot_pae_score_masks` are drawn
    from, counted rather than drawn, so the figure and the numbers beside it
    cannot disagree.

    Args:
        pair:       The ordered chain pair.
        contacts:   Interface contacts of the same pair.
        pae_cutoff: ipSAE's PAE cutoff.
        lis_cutoff: LIS's PAE cutoff.
        label_x:    Display name for `chain_x`.
        label_y:    Display name for `chain_y`.

    Returns:
        A header naming the block and its slice, then one line per score.
    """
    name_x = _chain_label(pair.chain_x, label_x)
    name_y = _chain_label(pair.chain_y, label_y)
    masks = score_masks(pair, contacts, pae_cutoff=pae_cutoff, lis_cutoff=lis_cutoff)
    block_cells = pair.block_xy.size
    nx, ny = pair.nx, pair.ny

    lines = [f'Cells of the {name_x} → {name_y} inter-chain block '
             f'(pae_matrix[:{nx}, {nx}:{nx + ny}]) used by each score:']
    for name, mask in masks.items():
        n_cells = int(mask.sum())
        lines.append(f'  {name:12s}: {n_cells:6d} cells '
                     f'({100 * n_cells / block_cells:.1f}%)')
    return '\n'.join(lines)


# -- Section 4: one report per score -----------------------------------------
# Each of these is the print block that used to sit under the `compute_*` call
# in its own notebook cell, moved verbatim. The cell keeps the call, because the
# arguments are what the section is teaching; it gives up the formatting.

def _threshold_line(label: str, score_name: str, both_provenances: bool = False) -> str:
    """One `threshold : green >= g, amber >= a (provenance)` line, aligned to `label`."""
    threshold = THRESHOLDS[score_name]
    if both_provenances:
        return (f'  {label}: green ≥ {threshold.green} '
                f'({threshold.green_provenance}), '
                f'amber ≥ {threshold.amber} ({threshold.amber_provenance})')
    return (f'  {label}: green ≥ {threshold.green}, '
            f'amber ≥ {threshold.amber} ({threshold.green_provenance})')


def format_iptm_report(
    result: DirectionalPair,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
    full_x: str = "",
    full_y: str = "",
    res_ids_x: Optional[np.ndarray] = None,
    res_ids_y: Optional[np.ndarray] = None,
) -> str:
    """
    Section 4.1: ipTM_d0chn, its `d0`, its peak residue in each direction.

    Opens with the ordered pair itself, because 4.1 is the first subsection and
    every "x -> y" label in Section 4 reads against that order.

    Args:
        result:    `compute_iptm_d0chn` output.
        label_x:   Compact display name for `chain_x`.
        label_y:   Compact display name for `chain_y`.
        full_x:    Full name for `chain_x`, for the ordered-pair header. Blank
                   omits the header.
        full_y:    Full name for `chain_y`.
        res_ids_x: Optional `(nx,)` residue numbers, so the peak is named by the
                   model's own residue number and not only by an index.
        res_ids_y: Optional `(ny,)` residue numbers.

    Returns:
        The multi-line report.
    """
    name_x = _chain_label(result.chain_x, label_x)
    name_y = _chain_label(result.chain_y, label_y)
    nx = int(result.forward.values.shape[0])
    ny = int(result.reverse.values.shape[0])
    _, band = traffic_light(result.score, 'iptm_d0chn')

    lines: List[str] = []
    if full_x or full_y:
        lines += [f'Ordered pair: {name_x} → {name_y}',
                  f'  {name_x} = {full_x}',
                  f'  {name_y} = {full_y}',
                  '']
    lines += [f'ipTM_d0chn : {result.score:.4f}   [{band}]',
              f'  d0chn                 : {result.d0:.4f} '
              f'(from n0chn = {result.n0} = {nx} + {ny} residues)']
    # The reported value is one residue's number, so name that residue.
    for row, col, profile, res_ids in ((name_x, name_y, result.forward, res_ids_x),
                                       (name_y, name_x, result.reverse, res_ids_y)):
        index = profile.argmax_index
        number = index if res_ids is None else np.asarray(res_ids)[index]
        lines.append(f'  peak residue, {row} → {col}: {number} (index {index})')
    lines.append(_threshold_line('threshold             ', 'iptm_d0chn'))
    return '\n'.join(lines)


def format_ipsae_report(result: IPSAEResult) -> str:
    """
    Section 4.2: the three ipSAE variants, their three `d0` values, and the ordering.

    `d0chn >= d0dom >= d0res` is a theorem, not a measurement, so it is checked
    here rather than left for the reader to notice. It is *reported* rather than
    asserted: a violation would mean the computation is wrong, not the model, and
    that is worth seeing next to the numbers rather than as a traceback.

    Args:
        result: `compute_ipsae` output.

    Returns:
        The multi-line report.
    """
    lines = [f'PAE cutoff: {result.pae_cutoff:.0f} Å (strict <)', '']
    for key, variant in result.variants.items():
        _, band = traffic_light(variant.score, key)
        lines.append(f'  {SCORE_DISPLAY_NAMES[key]:12s}: {variant.score:.4f}   [{band}]')

    # The d0 values are the point of the exercise: same PAE cells, three scales.
    lines += ['',
              f'  d0chn : {result.d0chn_value:.4f}  (n0chn = {result.n0chn})',
              f'  d0dom : {result.d0dom.d0:.4f}  (n0dom = {result.d0dom.n0}, '
              'of the direction that supplied the reported value)',
              f'  d0res : {result.d0res.d0:.4f}  (n0res = {result.d0res.n0}, '
              'of the peak residue alone)']

    ordered = (result.d0chn.score >= result.d0dom.score - 1e-6
               and result.d0dom.score >= result.d0res.score - 1e-6)
    lines += ['',
              '  ordering d0chn ≥ d0dom ≥ d0res : '
              + ('holds' if ordered else 'VIOLATED (this should be impossible)'),
              '  d0chn - d0res spread          : '
              f'{result.d0chn.score - result.d0res.score:+.4f}',
              '  A wide spread means the confidently predicted interface is small',
              '  relative to the two chains. Quote ipSAE_d0res.']
    return '\n'.join(lines)


def format_pdockq_report(result: PDockQResult) -> str:
    """
    Section 4.3: pDockQ, its two ingredients, and the symmetry check.

    Args:
        result: `compute_pdockq` output.

    Returns:
        The multi-line report.
    """
    _, band = traffic_light(result.score, 'pdockq')
    return '\n'.join([
        f'pDockQ : {result.score:.4f}   [{band}]',
        f'  contact pairs (npairs) : {result.n_contact_pairs}   '
        f'(CB-CB ≤ {result.dist_cutoff:.0f} Å)',
        f'  interface residues     : {result.n_interface_residues}   '
        '(counted once each; reported by ipsae.py but not used in the score)',
        f'  mean interface pLDDT   : {result.mean_plddt:.2f}',
        f'  x = mean_pLDDT × log10(npairs) : {result.x:.4f}',
        f'  symmetric              : {result.symmetric} '
        '(checked by computing both orientations, not assumed)',
        _threshold_line('threshold              ', 'pdockq')
        + '; sigmoid bounded on [0.018, 0.742]',
    ])


def format_pdockq2_report(
    result: PDockQ2Result,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
) -> str:
    """
    Section 4.4: pDockQ2, both directions, and the ingredient pDockQ does not have.

    Args:
        result:  `compute_pdockq2` output.
        label_x: Display name for `chain_x`.
        label_y: Display name for `chain_y`.

    Returns:
        The multi-line report.
    """
    name_x = _chain_label(result.chain_x, label_x)
    name_y = _chain_label(result.chain_y, label_y)
    winner = result.winning_direction
    _, band = traffic_light(result.score, 'pdockq2')
    return '\n'.join([
        f'pDockQ2 : {result.score:.4f}   [{band}]',
        f'  reported direction  : {winner.chain_row} → {winner.chain_col} '
        '(the larger of the two)',
        f'  {name_x} → {name_y} : {result.forward_score:.4f}      '
        f'{name_y} → {name_x} : {result.reverse_score:.4f}      '
        f'|Δ| {result.delta:.4f}',
        f'  contact pairs       : {winner.n_contact_pairs}',
        f'  mean_ptm (d0 = 10)  : {winner.mean_ptm:.4f}   '
        '← the ingredient pDockQ does not have',
        f'  mean interface pLDDT: {winner.mean_plddt:.2f}',
        f'  x = mean_pLDDT × mean_ptm : {winner.x:.4f}',
        _threshold_line('threshold           ', 'pdockq2', both_provenances=True)
        + '; bounded on [0.005, 1.315], not a probability',
    ])


def format_lis_report(
    result: LISResult,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
) -> str:
    """
    Section 4.5: LIS, both directions, and the size of the local interaction area.

    The LIA is reported alongside the score because Kim 2024 p. 8 gives high LIS
    over a very small LIA as a false-positive signature.

    Args:
        result:  `compute_lis` output.
        label_x: Display name for `chain_x`.
        label_y: Display name for `chain_y`.

    Returns:
        The multi-line report.
    """
    name_x = _chain_label(result.chain_x, label_x)
    name_y = _chain_label(result.chain_y, label_y)
    names = {result.chain_x: name_x, result.chain_y: name_y}
    _, band = traffic_light(result.score, 'lis')

    lines = [f'LIS : {result.score:.4f}   [{band}]   '
             '(the MEAN of the two directions, not the max)',
             f'  {name_x} → {name_y} : {result.forward_score:.4f}      '
             f'{name_y} → {name_x} : {result.reverse_score:.4f}      '
             f'|Δ| {result.delta:.4f}']
    for direction in (result.forward, result.reverse):
        lines.append(
            f'  LIA, {names[direction.chain_row]} → {names[direction.chain_col]} : '
            f'{direction.n_valid_pairs} of {direction.n_pairs} inter-chain cells '
            f'below {result.lis_cutoff:.0f} Å '
            f'({100 * direction.fraction_valid:.1f}%)')
    lines.append(_threshold_line('threshold        ', 'lis', both_provenances=True))
    return '\n'.join(lines)


def format_score_values(scores: Mapping[str, float]) -> str:
    """
    Section 4.6: the seven collected values, one per line, under a rule.

    Args:
        scores: `{THRESHOLDS key: value}`.

    Returns:
        A leading blank line, the rule, then one line per score.

    Example
    -------
    >>> print(format_score_values({'lis': 0.5, 'pdockq': 0.25}))
    <BLANKLINE>
    ── Score Results ──────────────────────────────
      LIS               : 0.5000
      pDockQ            : 0.2500
    """
    lines = ['', '── Score Results ──────────────────────────────']
    lines += [f'  {SCORE_DISPLAY_NAMES[name]:18s}: {value:.4f}'
              for name, value in scores.items()]
    return '\n'.join(lines)


def format_domain_sizes(
    result: IPSAEResult,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
) -> str:
    """
    The per-direction `n0dom` and `d0dom` behind the directional report.

    n0dom is the mechanism, not a separate finding: the domain size that sets
    d0dom is itself counted per direction (R003), so the two d0dom columns of the
    directional report are not the same measurement under two names.

    Args:
        result:  `compute_ipsae` output.
        label_x: Display name for `chain_x`.
        label_y: Display name for `chain_y`.

    Returns:
        A leading blank line, then the two lines.
    """
    name_x = _chain_label(result.chain_x, label_x)
    name_y = _chain_label(result.chain_y, label_y)
    return '\n'.join([
        '',
        f'  n0dom  {name_x} → {name_y}: {result.n0dom_xy}'
        f'   {name_y} → {name_x}: {result.n0dom_yx}'
        f'   |Δ| {result.n0dom_delta}',
        f'  d0dom  {name_x} → {name_y}: {result.d0dom_xy:.4f}'
        f'   {name_y} → {name_x}: {result.d0dom_yx:.4f}',
    ])


# -- Section 5: interface pLDDT ----------------------------------------------

@dataclass(frozen=True)
class InterfacePLDDT:
    """
    Interface pLDDT against the rest of the model, as numbers rather than a figure.

    pDockQ and pDockQ2 both average pLDDT over interface residues, so a low score
    has two very different causes: a globally uncertain protein, or a confident
    protein with an uncertain interface. These are the numbers that separate them,
    and the same values `plot_plddt_distribution` draws.

    Attributes:
        interface:     pLDDT of every interface residue, both chains pooled.
        non_interface: pLDDT of every other residue.
        low_cutoff:    The pLDDT below which a residue counts as low.
    """

    interface: np.ndarray
    non_interface: np.ndarray
    low_cutoff: float = 70.0

    @property
    def n_interface(self) -> int:
        """Interface residues, both chains."""
        return int(self.interface.shape[0])

    @property
    def n_low_interface(self) -> int:
        """Interface residues below `low_cutoff`."""
        return int((self.interface < self.low_cutoff).sum())

    @property
    def mean_interface(self) -> float:
        """Mean interface pLDDT. This is pDockQ's `mean_plddt`."""
        return float(self.interface.mean())


def interface_plddt_stats(
    contacts: InterfaceContacts,
    plddt_x: np.ndarray,
    plddt_y: np.ndarray,
    low_cutoff: float = 70.0,
) -> InterfacePLDDT:
    """
    Split both chains' pLDDT into interface and non-interface.

    Args:
        contacts:   `detect_interface` result.
        plddt_x:    `(nx,)` pLDDT for `chain_x`.
        plddt_y:    `(ny,)` pLDDT for `chain_y`.
        low_cutoff: The "low pLDDT" line, 70 by AlphaFold's own convention.

    Returns:
        An `InterfacePLDDT`.

    Example
    -------
    >>> contacts = detect_interface(
    ...     ChainCoords('A', np.array([[0., 0., 0.], [50., 0., 0.]]),
    ...                 np.array([1, 2]), np.array(['ALA', 'ALA']), None),
    ...     ChainCoords('B', np.array([[1., 0., 0.], [60., 0., 0.]]),
    ...                 np.array([1, 2]), np.array(['ALA', 'ALA']), None))
    >>> stats = interface_plddt_stats(contacts, np.array([90.0, 40.0]),
    ...                               np.array([50.0, 30.0]))
    >>> stats.n_interface, stats.n_low_interface, round(stats.mean_interface, 1)
    (2, 1, 70.0)
    """
    interface = np.concatenate([np.asarray(plddt_x)[contacts.mask_x],
                                np.asarray(plddt_y)[contacts.mask_y]])
    non_interface = np.concatenate([np.asarray(plddt_x)[~contacts.mask_x],
                                    np.asarray(plddt_y)[~contacts.mask_y]])
    return InterfacePLDDT(interface=interface, non_interface=non_interface,
                          low_cutoff=low_cutoff)


def format_plddt_stats(stats: InterfacePLDDT) -> str:
    """
    Interface against non-interface pLDDT, as the three lines above the figure.

    Args:
        stats: `interface_plddt_stats` result.

    Returns:
        Three lines.
    """
    n_if = max(stats.n_interface, 1)
    return '\n'.join([
        f'Interface pLDDT   mean={stats.interface.mean():.1f}, '
        f'median={np.median(stats.interface):.1f}',
        f'Non-interface     mean={stats.non_interface.mean():.1f}, '
        f'median={np.median(stats.non_interface):.1f}',
        f'Low-pLDDT (<{stats.low_cutoff:.0f}) interface residues: '
        f'{stats.n_low_interface} / {stats.n_interface} '
        f'({100 * stats.n_low_interface / n_if:.1f}%)',
    ])


# -- Sections 2 and 6: drawing the 3D views ----------------------------------

MOLVIEWSPEC_NO_URL_MESSAGE: str = (
    'Skipping the 3D views: Mol* downloads the structure itself, so local file '
    'mode has no URL to hand it.\nThe contact map above, and every score in '
    'this notebook, are unaffected.'
)
"""Why local-file mode draws no 3D views. Not a failure: Mol* fetches the
structure from a URL and an already-parsed mmCIF string cannot be handed to it."""


@dataclass(frozen=True)
class ViewRenderer:
    """
    One resolved structure source, reused by every 3D view in the notebook.

    Section 2 resolves the source once for View 1, and Section 6's Views 2 to 6
    reuse the same object (R030), so the two sections cannot end up pointing Mol*
    at different files. Each view is built lazily and its failures are caught per
    view, so one broken viewer does not cost you the others.

    Attributes:
        source:  The structure Mol* downloads, or `None` when the 3D section
                 cannot run. `view_renderer` has already said why.
        label_x: Compact display name for the first chain, for view captions.
        label_y: Compact display name for the second chain.
    """

    source: Optional[StructureSource]
    label_x: str = ""
    label_y: str = ""

    @property
    def available(self) -> bool:
        """Whether any view can be drawn at all."""
        return self.source is not None

    def render(self, key: str, build: Any, legend: Optional[Any] = None) -> None:
        """
        Draw one view with its caption and legend, or say why it was skipped.

        Args:
            key:    A `MVS_VIEW_LABELS` key, which supplies the caption.
            build:  `source -> State`. Called only if a source was resolved, so
                    the builder never runs in an environment that cannot show it.
            legend: `() -> Sequence[LegendEntry]`, rendered under the viewer.

        Returns:
            `None`. Displays as a side effect.
        """
        from IPython.display import HTML, display

        label = format_view_label(key, self.label_x, self.label_y)
        if self.source is None:
            print(f'{label}\n  skipped: no structure source resolved above.')
            return
        try:
            show_mol_view(build(self.source), label)
            if legend is not None:
                display(HTML(legend_html(legend())))
        except Exception as exc:  # noqa: BLE001  (one bad view must not stop the rest)
            print(f'{label} failed: {exc}')


def view_renderer(
    prediction: "Optional[AFDBPrediction | Mapping[str, Any]]",
    label_x: str = "",
    label_y: str = "",
) -> ViewRenderer:
    """
    Resolve the structure Mol* will download, once, and say so if it cannot.

    Three things can stop the 3D section, and each gets its own sentence rather
    than a shared "unavailable": `molviewspec` is not installed, local-file mode
    left no URL to hand Mol*, or the metadata carries no usable structure URL.

    Args:
        prediction: Fetched metadata, or `None` in local-file mode.
        label_x:    Compact display name for the first chain.
        label_y:    Compact display name for the second chain.

    Returns:
        A `ViewRenderer`, whose `.source` is `None` if any of the three applies.

    Note:
        Prints the reason as a side effect when there is one.

    Example
    -------
    >>> renderer = view_renderer({'bcifUrl': 'https://x/y.bcif'}, 'A', 'B')
    >>> renderer.available and renderer.source.format
    'bcif'
    """
    if not molviewspec_available():
        print(MOLVIEWSPEC_MISSING_MESSAGE)
        return ViewRenderer(None, label_x, label_y)
    if prediction is None:
        print(MOLVIEWSPEC_NO_URL_MESSAGE)
        return ViewRenderer(None, label_x, label_y)
    try:
        source = resolve_structure_source(prediction)
    except ValueError as exc:
        print(f'Skipping the 3D views: {exc}')
        return ViewRenderer(None, label_x, label_y)
    return ViewRenderer(source, label_x, label_y)


def format_view_conventions(full_x: str = "", full_y: str = "") -> str:
    """
    The three cutoffs Section 6's views are painted against, stated once.

    Views 3 to 6 all divide residues by the same two lines -- a contact distance
    and an ipSAE_d0res level -- plus, for View 6, a contact-quality level. Saying
    them once above the views keeps six legends from each restating them.

    Args:
        full_x: Full name of the first chain.
        full_y: Full name of the second chain.

    Returns:
        Four lines, ASCII only: these sit above a Mol* iframe, and the notebook
        keeps that block plain so it reads the same in every front end.

    Example
    -------
    >>> print(format_view_conventions('P1 (A)', 'P1 (B)'))
    Chains in these views: P1 (A); P1 (B)
    Contact, Views 4, 5, 6:     CB within 8.0 A (CA for glycine)
    Confident PAE, Views 4, 6:  ipSAE_d0res >= 0.60 (THRESHOLDS['ipsae_d0res'].amber, AFDB's release edge)
    Well-placed contacts, V6:   mean contact ptm >= 0.50, i.e. mean contact PAE better than 10 A
    """
    contact_pae = contact_ptm_to_pae(MVS_CONTACT_PTM_THRESHOLD)
    return '\n'.join([
        f'Chains in these views: {full_x}; {full_y}',
        f'Contact, Views 4, 5, 6:     CB within {DIST_CUTOFF:.1f} A '
        '(CA for glycine)',
        f'Confident PAE, Views 4, 6:  ipSAE_d0res >= '
        f'{MVS_DISAGREEMENT_THRESHOLD:.2f} '
        "(THRESHOLDS['ipsae_d0res'].amber, AFDB's release edge)",
        f'Well-placed contacts, V6:   mean contact ptm >= '
        f'{MVS_CONTACT_PTM_THRESHOLD:.2f}, i.e. mean contact PAE better '
        f'than {contact_pae:.0f} A',
    ])


# -- Section 7: what closes the notebook -------------------------------------

def format_interface_statistics(
    contacts: InterfaceContacts,
    plddt: InterfacePLDDT,
    label_x: "Optional[str | ChainLabel]" = None,
    label_y: "Optional[str | ChainLabel]" = None,
    rule_width: int = 65,
) -> str:
    """
    The physical facts about the interface, under the diagnostic prose.

    Every number here was computed earlier in the notebook and is read back, not
    recomputed: this block cannot disagree with the sections above it.

    Args:
        contacts:   `detect_interface` result.
        plddt:      `interface_plddt_stats` result.
        label_x:    Display name for `chain_x`.
        label_y:    Display name for `chain_y`.
        rule_width: Width of the closing rule.

    Returns:
        Six lines and the closing rule.
    """
    name_x = _chain_label(contacts.chain_x, label_x)
    name_y = _chain_label(contacts.chain_y, label_y)
    return '\n'.join([
        'Interface statistics:',
        f'  Contact pairs     : {contacts.n_contact_pairs}',
        f'  Interface res, {name_x}: {contacts.n_interface_residues_x}',
        f'  Interface res, {name_y}: {contacts.n_interface_residues_y}',
        f'  Mean pLDDT (if)   : {plddt.mean_interface:.1f}',
        f'  Low pLDDT (<{plddt.low_cutoff:.0f}) if: {plddt.n_low_interface} / '
        f'{plddt.n_interface}',
        '=' * rule_width,
    ])


def format_references() -> str:
    """
    Every source the notebook's numbers rest on, metric and threshold separately.

    The citation policy is in the module docstring: a metric belongs to whoever
    published it, and the ipSAE 0.6 cutoff and the four confidence bands are
    AlphaFold DB's, not Dunbrack's. This block is that policy rendered.

    Returns:
        A leading blank line, then the reference block.

    Example
    -------
    >>> print(format_references())      # doctest: +ELLIPSIS
    <BLANKLINE>
    References. The metric and its threshold are credited separately:
      ipSAE metric    : Dunbrack (2025) biorxiv 2025.02.10.637595
      ...
      Background      : https://...
    """
    return '\n'.join([
        '',
        'References. The metric and its threshold are credited separately:',
        '  ipSAE metric    : Dunbrack (2025) biorxiv 2025.02.10.637595',
        '  pDockQ          : Bryant et al. (2022) Nat Commun s41467-022-28865-w',
        '  pDockQ2 metric  : Zhu et al. (2023) Bioinformatics btad424',
        '  LIS             : Kim et al. (2024) biorxiv 2024.02.19.580970',
        '  ipSAE / pDockQ2 thresholds, bands and the joint release criterion:',
        '                    Han, Tsenkov, Venanzi et al. (2026)',
        '                    biorxiv 10.64898/2026.03.27.714458v2',
        '  ipsae.py v4     : github.com/DunbrackLab/IPSAE',
        f'  Background      : {AFDB_NEWS_URL}',
    ])
