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

The plotting section is a behaviour-preserving move of the notebook's inline
figures, not a redesign: the sizing (R050), the palette (R051) and the
score-mask panel (R052) are corrected later, and `PAE_CMAP` is the seam R051
changes.

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
from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Mapping, Optional, Sequence, Tuple

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import requests
import seaborn as sns
from matplotlib.colors import Colormap
from matplotlib.figure import Figure

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
    "SCORE_PROFILE_COLOURS",
    "PLDDT_BANDS",
    "PLDDT_BAND_COLOURS",
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
    "MVS_VIEW_LABELS",
    # MolViewSpec: structure source
    "StructureSource",
    "resolve_structure_source",
    # MolViewSpec: per-residue colouring
    "ColourRun",
    "colour_runs",
    "value_colours",
    "add_residue_colours",
    # MolViewSpec: display and builders
    "show_mol_view",
    "build_chain_overview_view",
    "build_plddt_view",
    "build_interface_value_view",
    "disagreement_categories",
    "build_disagreement_view",
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
        ValueError:         If the endpoint returns an empty array, i.e. the
                            accession does not exist.
        requests.HTTPError: On a non-2xx status.
    """
    payload = download_json(AFDB_PREDICTION_URL.format(accession=accession), timeout=timeout)
    entries = payload if isinstance(payload, list) else [payload]
    if not entries:
        raise ValueError(f"No AFDB prediction found for {accession!r}.")
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

    if prediction is None:
        assembly = "not fetched (local file mode)"
        version = "not fetched (local file mode)"
    else:
        kind = prediction.shared_field("assemblyType", default=None)
        if kind is None:
            # Derivable from the identity grouping when the endpoint omits it --
            # but only when every chain was actually described.
            if undescribed or not groups:
                kind = "Unknown"
            else:
                kind = "Homo" if len(groups) == 1 else "Hetero"
        state = prediction.shared_field("oligomericState", default=None)
        distinct = len(groups)
        described = f", {len(chain_ids) - len(undescribed)} described" if undescribed else ""
        assembly = (
            f"{kind}{' ' + str(state) if state else ''} "
            f"({len(chain_ids)} chains{described}, {distinct} distinct "
            f"protein{'' if distinct == 1 else 's'})"
        )
        version = str(prediction.shared_field("modelVersion", "latestVersion",
                                              default="N/A"))

    lines: List[str] = ["=" * width, "COMPLEX METADATA REPORT", "=" * width]
    lines.append(row("Accession", str(
        accession or (prediction.accession if prediction is not None else "N/A"))))
    lines.append(row("Assembly", assembly))
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
        ValueError: If the chain sets or any length disagree.

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
    """
    from_structure = {chain_id: chain.n_residues for chain_id, chain in chains.items()}
    from_document = {span.chain_id: span.length for span in document.spans}

    if set(from_structure) != set(from_document):
        raise ValueError(
            f"Chain sets disagree: structure has {sorted(from_structure)}, "
            f"document has {sorted(from_document)}."
        )
    disagreeing = {
        chain_id: (from_structure[chain_id], from_document[chain_id])
        for chain_id in from_structure
        if from_structure[chain_id] != from_document[chain_id]
    }
    if disagreeing:
        detail = ", ".join(
            f"{chain_id}: structure={s}, document={d}"
            for chain_id, (s, d) in sorted(disagreeing.items())
        )
        raise ValueError(f"Chain lengths disagree ({detail}).")
    return dict(sorted(from_structure.items()))


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
# a `cmap=` override. R051 -- "green default palette with a switchable
# alternative" -- is therefore a change of *this default*, not a rewrite of the
# functions below. The default is `'RdBu_r'` today only because that is what the
# notebook draws; R051 flips it to a sequential green with dark = low PAE =
# confident, and `PAE_CMAP_CHOICES` is where the candidates are named.
#
# What this section deliberately does NOT do
# ------------------------------------------
# R013 is a behaviour-preserving move. The sizing (R050), the palette (R051) and
# the score-mask panel's aspect ratio and shared colour bar (R052) are all known
# defects, and all three are left exactly as the notebook draws them today: the
# scoring code could be corrected early because `ipsae.py` is an objective
# oracle, but a figure has no oracle, so a visual change landed here would be
# indistinguishable from a regression. Specifically preserved on purpose:
# `aspect='auto'` everywhere, `figsize=(10, 9)` for the full PAE matrix, and
# `fig.colorbar(im, ax=axes[-1])` in the 2x2 panel, which steals width from the
# fourth axis alone and is why that panel renders narrower than the other three.


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

PAE_CMAP: str = "RdBu_r"
"""Default colormap for every PAE figure. **This is the R051 seam.**

Resolved at call time by `resolve_pae_cmap`, so reassigning
`complex_interface_utils.PAE_CMAP` retunes every PAE figure at once, and R051
becomes a one-line change of this default plus the accompanying prose.

`'RdBu_r'` is a diverging map (0 = dark blue, mid = white, max = dark red) and is
kept only because it is what the notebook draws today. R051 replaces it with a
sequential green in which *dark* = low PAE = confident, so that confident regions
stay legible against a white background; see `PAE_CMAP_CHOICES`.
"""

PAE_CMAP_CHOICES: Dict[str, str] = {
    "rdbu": "RdBu_r",
    "green": "Greens_r",
    "colourblind_safe": "viridis_r",
}
"""Named shorthands accepted anywhere a PAE colormap is taken.

Candidates for R051, which owns the final choice: `'green'` is the AFDB-style
sequential green with dark = low PAE, and `'colourblind_safe'` is the
perceptually uniform alternative for readers who cannot separate the green ramp.
`'rdbu'` is the current default, retained for continuity.
"""

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
"""The five values the agreement matrix compares: one ipSAE variant plus the four
independent scores. The other two ipSAE variants are excluded because
`d0chn >= d0dom >= d0res` is a theorem, so including them would show three
guaranteed agreements as if they were three confirmations."""


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
        `figure.dpi` is 150 here, which is what makes the 10x9 inch PAE figure
        1500x1350 px and triggers the notebook's in-cell scrollbar. R050 fixes
        that by resizing the figure, not by lowering this.
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
    whole point of the seam -- it is what lets R051 change one name and retune
    every PAE figure.

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
    'RdBu_r'
    >>> resolve_pae_cmap('green').name
    'Greens_r'
    >>> resolve_pae_cmap('viridis').name
    'viridis'
    """
    requested = PAE_CMAP if cmap is None else cmap
    if isinstance(requested, Colormap):
        return requested
    return matplotlib.colormaps[PAE_CMAP_CHOICES.get(requested, requested)]


def _chain_label(chain_id: str, label: Optional[str] = None) -> str:
    """`label` if given, else `'Chain <id>'`. R021 supplies real protein names."""
    return label if label else f"Chain {chain_id}"


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
    for minimum, colour, _ in PLDDT_BANDS:
        if value > minimum:
            return colour
    return PLDDT_BANDS[-1][1]


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
    label_x: Optional[str] = None,
    label_y: Optional[str] = None,
    cmap: str | Colormap = DIST_CMAP,
    figsize: Tuple[float, float] = (14.0, 6.0),
) -> Figure:
    """
    The interface: a distance-coloured contact map beside per-chain coverage bars.

    Left panel: every contact within the cutoff, coloured by CB-CB distance, rows
    `chain_x` and columns `chain_y`. Non-contacts are `nan` and so render as the
    axes background. Right panel: one bar per residue of each chain, amber where
    that residue touches the partner chain, showing *where along the sequence*
    the interface sits -- one contiguous patch reads very differently from a
    scatter of isolated residues.

    Args:
        contacts: Interface contacts of one ordered chain pair.
        label_x:  Display name for `chain_x`; defaults to `'Chain <id>'`.
        label_y:  Display name for `chain_y`.
        cmap:     Colormap for the distance panel. Not `PAE_CMAP`: this panel
                  shows distance, and the two must stay visually distinct.
        figsize:  Figure size in inches.

    Returns:
        The `Figure`. Nothing is shown or saved; the caller decides.
    """
    name_x = _chain_label(contacts.chain_x, label_x)
    name_y = _chain_label(contacts.chain_y, label_y)
    nx, ny = contacts.contact_mask.shape

    fig = Figure(figsize=figsize)
    axes = fig.subplots(1, 2)

    ax = axes[0]
    contact_distances = np.where(contacts.contact_mask, contacts.dist_matrix, np.nan)
    im = ax.imshow(contact_distances, aspect='auto', origin='lower',
                   cmap=cmap, vmin=0, vmax=contacts.dist_cutoff)
    fig.colorbar(im, ax=ax, label='CB-CB distance (Å)')
    ax.set_xlabel(f'{name_y} residue index')
    ax.set_ylabel(f'{name_x} residue index')
    ax.set_title('Interface Contact Map\n'
                 f'(contacts ≤ {contacts.dist_cutoff:.0f} Å, coloured by distance)')

    ax2 = axes[1]
    bar_height = 0.35
    for mask, bottom in ((contacts.mask_x, 0.6), (contacts.mask_y, 0.1)):
        ax2.bar(np.arange(len(mask)), bar_height, bottom=bottom,
                color=[COLOUR_IF if is_if else COLOUR_NON_IF for is_if in mask],
                width=1.0, linewidth=0)

    ax2.set_xlim(0, max(nx, ny))
    ax2.set_ylim(0, 1.1)
    ax2.set_yticks([0.275, 0.775])
    ax2.set_yticklabels([name_y, name_x])
    ax2.set_xlabel('Residue index')
    ax2.set_title('Interface Coverage\n'
                  '(amber = at interface, grey = non-interface)')
    ax2.legend(
        handles=[mpatches.Patch(color=COLOUR_IF, label='Interface'),
                 mpatches.Patch(color=COLOUR_NON_IF, label='Non-interface')],
        loc='upper right', fontsize=9)

    fig.tight_layout()
    return fig


def plot_pae_matrix(
    pae: PAEMatrix,
    chain_x: Optional[str] = None,
    chain_y: Optional[str] = None,
    accession: str = "",
    cmap: Optional[str | Colormap] = None,
    figsize: Tuple[float, float] = (10.0, 9.0),
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
        cmap:      Colormap override; `None` uses `PAE_CMAP` (the R051 seam).
        figsize:   Figure size in inches. `(10, 9)` at `figure.dpi = 150` is
                   1500x1350 px, which overflows the notebook output area --
                   preserved deliberately here and fixed by R050.

    Returns:
        The `Figure`.

    Raises:
        ValueError: If the document has fewer than two chains, or if `chain_x`
            and `chain_y` are the same chain.
        KeyError:   If a named chain is not in the document.

    Note:
        `aspect='auto'` is what the notebook uses, and it stretches the matrix to
        the axes box, so a heterodimer's square matrix renders non-square. R050
        switches it to `'equal'`.
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

    fig = Figure(figsize=figsize)
    ax = fig.subplots()

    im = ax.imshow(pae.matrix, aspect='auto', origin='upper',
                   cmap=resolve_pae_cmap(cmap), vmin=0, vmax=pae.max_pae)
    fig.colorbar(im, ax=ax, label='PAE (Å) — lower = more confident')

    # One dashed line per internal chain boundary. Two chains give the notebook's
    # single pair of lines at nx - 0.5; more chains give one pair each.
    offset = 0
    for span in pae.spans[:-1]:
        offset += span.length
        ax.axhline(offset - 0.5, color='white', linewidth=2, linestyle='--')
        ax.axvline(offset - 0.5, color='white', linewidth=2, linestyle='--')

    def _centre(chain_id: str) -> float:
        span_slice = pae.chain_slice(chain_id)
        return (span_slice.start + span_slice.stop) / 2.0

    cx, cy = _centre(chain_x), _centre(chain_y)
    for col, row, text in ((cx, cx, f'Intra {chain_x}'),
                           (cy, cx, f'Inter\n{chain_x}→{chain_y}'),
                           (cx, cy, f'Inter\n{chain_y}→{chain_x}'),
                           (cy, cy, f'Intra {chain_y}')):
        ax.text(col, row, text, ha='center', va='center',
                color='white', fontsize=11, fontweight='bold', alpha=0.8)

    if len(ids) == 2:
        first, second = ids
        ax.set_xlabel(f'Residue index (chain {first}: 0 to n{first}-1, '
                      f'chain {second}: n{first} to end)')
    else:
        ax.set_xlabel('Residue index')
    ax.set_ylabel('Residue index')
    head = f'Full PAE Matrix — {accession}' if accession else 'Full PAE Matrix'
    ax.set_title(f'{head}\n'
                 f'(dashed line = chain boundary between {chain_x} and {chain_y})')

    fig.tight_layout()
    return fig


def plot_pae_score_masks(
    pair: ChainPairPAE,
    contacts: InterfaceContacts,
    max_pae: float,
    pae_cutoff: float = PAE_CUTOFF,
    lis_cutoff: float = LIS_CUTOFF,
    label_x: Optional[str] = None,
    label_y: Optional[str] = None,
    cmap: Optional[str | Colormap] = None,
    figsize: Tuple[float, float] = (14.0, 12.0),
) -> Figure:
    """
    Four views of the same inter-chain PAE block: what each score actually reads.

    Every panel shows `pair.block_xy`; they differ only in which cells are left
    coloured. This is the figure that explains why four scores computed from one
    matrix can disagree -- they are not weighting the same evidence differently,
    they are reading different subsets of it.

    Args:
        pair:       The ordered chain pair's PAE quadrants.
        contacts:   Interface contacts of the **same** ordered pair.
        max_pae:    `vmax`; use `PAEMatrix.max_pae` so every PAE figure of one
                    model shares a scale.
        pae_cutoff: ipSAE's cutoff, named in the panel title rather than
                    hard-coded into it.
        lis_cutoff: LIS's cutoff, likewise.
        label_x:    Display name for `chain_x`.
        label_y:    Display name for `chain_y`.
        cmap:       Colormap override; `None` uses `PAE_CMAP`.
        figsize:    Figure size in inches.

    Returns:
        The `Figure`. Per-panel cell counts are in the panel titles; the same
        masks are available from `score_masks` if a caller wants the numbers.

    Raises:
        ValueError: If `contacts` and `pair` are not the same ordered pair.

    Note:
        Two known defects are preserved here on purpose and are R052's to fix.
        The colour bar is attached to the last axis alone, so the fourth panel
        renders narrower than the other three; and `aspect='auto'` distorts every
        block whenever `nx != ny`.
    """
    name_x = _chain_label(contacts.chain_x, label_x)
    name_y = _chain_label(contacts.chain_y, label_y)
    masks = score_masks(pair, contacts, pae_cutoff=pae_cutoff, lis_cutoff=lis_cutoff)
    titles = {
        'iptm_d0chn': 'ipTM  (all inter-chain, no cutoff)',
        'ipsae': f'ipSAE (PAE < {pae_cutoff:.0f} Å)',
        'lis': f'LIS   (PAE < {lis_cutoff:.0f} Å)',
        'pdockq2': f'pDockQ2 (CB-CB ≤ {contacts.dist_cutoff:.0f} Å contacts)',
    }
    block = pair.block_xy
    n_cells_total = block.size

    fig = Figure(figsize=figsize)
    axes = fig.subplots(2, 2).ravel()

    im = None
    for ax, (key, mask) in zip(axes, masks.items()):
        display_pae = block.astype(float).copy()
        display_pae[~mask] = np.nan
        grey_bg = np.ones(block.shape) * 35   # out-of-range sentinel

        ax.imshow(grey_bg, aspect='auto', origin='upper',
                  cmap='Greys', vmin=0, vmax=40, alpha=0.3)
        im = ax.imshow(display_pae, aspect='auto', origin='upper',
                       cmap=resolve_pae_cmap(cmap), vmin=0, vmax=max_pae)
        n_cells = int(mask.sum())
        frac = 100.0 * n_cells / n_cells_total if n_cells_total else 0.0
        ax.set_title(f'{titles[key]}\n({n_cells} cells used, '
                     f'{frac:.1f}% of inter-chain block)', fontsize=10)
        ax.set_xlabel(f'{name_y} residue')
        ax.set_ylabel(f'{name_x} residue')

    # Preserved defect (R052): one colour bar on the last axis only.
    fig.colorbar(im, ax=axes[-1], label='PAE (Å)')
    fig.suptitle(f'{contacts.chain_x}→{contacts.chain_y} Inter-chain PAE Block: '
                 'cells used by each score\n(grey = not used by this score)',
                 fontsize=12, y=1.01)
    fig.tight_layout()
    return fig


def plot_residue_score_profiles(
    iptm: DirectionalPair,
    ipsae: IPSAEResult,
    contacts: InterfaceContacts,
    plddt_x: np.ndarray,
    plddt_y: np.ndarray,
    label_x: Optional[str] = None,
    label_y: Optional[str] = None,
    figsize: Tuple[float, float] = (14.0, 10.0),
) -> Figure:
    """
    Per-residue score profiles, one panel per direction of the chain pair.

    Every score reported for a complex is one residue's number -- the maximum
    over the profile -- so this figure is where a headline value stops being a
    verdict and becomes a location: which residues carry the interface, whether
    the three ipSAE variants rank the same residue highest, and whether the peak
    sits on well-predicted backbone or on a low-pLDDT loop.

    The top panel is the `x -> y` direction (rows of `block_xy`, so residues of
    `chain_x`); the bottom is `y -> x`. They are genuinely two measurements, not
    one seen twice, because PAE is asymmetric.

    Args:
        iptm:     `compute_iptm_d0chn` result for the ordered pair.
        ipsae:    `compute_ipsae` result for the same ordered pair.
        contacts: Interface contacts of the same ordered pair; supplies the
                  shaded interface regions.
        plddt_x:  `(nx,)` pLDDT for `chain_x`.
        plddt_y:  `(ny,)` pLDDT for `chain_y`.
        label_x:  Display name for `chain_x`.
        label_y:  Display name for `chain_y`.
        figsize:  Figure size in inches.

    Returns:
        The `Figure`.

    Raises:
        ValueError: If the three results do not describe the same ordered chain
            pair, or if a pLDDT array's length does not match its chain.
    """
    pair_ids = (contacts.chain_x, contacts.chain_y)
    for name, result in (('ipTM', iptm), ('ipSAE', ipsae)):
        if (result.chain_x, result.chain_y) != pair_ids:
            raise ValueError(
                f"{name} describes pair ({result.chain_x}, {result.chain_y}) but "
                f"the contacts describe {pair_ids}; all three must agree."
            )

    panels = (
        (_chain_label(contacts.chain_x, label_x), 'forward',
         contacts.mask_x, np.asarray(plddt_x)),
        (_chain_label(contacts.chain_y, label_y), 'reverse',
         contacts.mask_y, np.asarray(plddt_y)),
    )

    fig = Figure(figsize=figsize)
    axes = fig.subplots(2, 1, sharex=False)

    for ax, (chain_label, direction, if_mask, plddt_chain) in zip(axes, panels):
        profiles = {
            'iptm_d0chn': getattr(iptm, direction),
            'ipsae_d0res': getattr(ipsae.d0res, direction),
            'ipsae_d0chn': getattr(ipsae.d0chn, direction),
            'ipsae_d0dom': getattr(ipsae.d0dom, direction),
        }
        n_res = profiles['iptm_d0chn'].values.shape[0]
        if plddt_chain.shape[0] != n_res:
            raise ValueError(
                f"{chain_label} has {n_res} residues in its score profile but "
                f"{plddt_chain.shape[0]} pLDDT values."
            )
        x = np.arange(n_res)

        ax.plot(x, profiles['iptm_d0chn'].values, label='ipTM',
                color=SCORE_PROFILE_COLOURS['iptm_d0chn'], linewidth=1.5)
        ax.plot(x, profiles['ipsae_d0res'].values, label='ipSAE d0res',
                color=SCORE_PROFILE_COLOURS['ipsae_d0res'], linewidth=1.5)
        ax.plot(x, profiles['ipsae_d0chn'].values, label='ipSAE d0chn',
                color=SCORE_PROFILE_COLOURS['ipsae_d0chn'], linewidth=1.2,
                linestyle='--')
        ax.plot(x, profiles['ipsae_d0dom'].values, label='ipSAE d0dom',
                color=SCORE_PROFILE_COLOURS['ipsae_d0dom'], linewidth=1.2,
                linestyle=':')

        # Secondary axis: pLDDT, on its own 0-100 scale.
        ax_plddt = ax.twinx()
        ax_plddt.fill_between(x, plddt_chain, alpha=0.15, color='grey', label='pLDDT')
        ax_plddt.set_ylabel('pLDDT (grey fill)', color='grey', fontsize=10)
        ax_plddt.set_ylim(0, 100)
        ax_plddt.tick_params(axis='y', labelcolor='grey')

        ax.fill_between(x, 0, 1, where=if_mask, alpha=0.12, color=COLOUR_IF,
                        label='Interface region')
        ax.set_ylim(0, 1)
        ax.set_xlabel('Residue index')
        ax.set_ylabel('Per-residue score (0–1)')
        ax.set_title(f'Per-Residue Score Profiles — {chain_label}')
        ax.legend(loc='upper left', fontsize=9)

    fig.tight_layout()
    return fig


def plot_plddt_distribution(
    contacts: InterfaceContacts,
    plddt_x: np.ndarray,
    plddt_y: np.ndarray,
    label_x: Optional[str] = None,
    label_y: Optional[str] = None,
    figsize: Tuple[float, float] = (14.0, 5.0),
) -> Figure:
    """
    Interface pLDDT against the rest of the model, as a histogram and a profile.

    pDockQ and pDockQ2 both average pLDDT over interface residues, so a low score
    has two very different causes: a globally uncertain protein, or a confident
    protein with an uncertain interface. The left panel separates them; the right
    panel says *where* the uncertain residues are, in AlphaFold's own colours.

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

    # Short labels here: these two strings read as sequence landmarks
    # ("A then B"), not as panel headings, so they take the bare chain id.
    name_x = label_x if label_x else contacts.chain_x
    name_y = label_y if label_y else contacts.chain_y

    if_plddt = np.concatenate([plddt_x[contacts.mask_x], plddt_y[contacts.mask_y]])
    ni_plddt = np.concatenate([plddt_x[~contacts.mask_x], plddt_y[~contacts.mask_y]])

    fig = Figure(figsize=figsize)
    axes = fig.subplots(1, 2)

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
                label=f'Chain {name_x}/{name_y} boundary')
    ax2.axhline(70, color='grey', linestyle=':', linewidth=1)
    ax2.set_xlabel(f'Residue index ({name_x} then {name_y})')
    ax2.set_ylabel('pLDDT')
    ax2.set_title('Per-residue pLDDT profile (AlphaFold colour scheme)')

    legend_patches = [
        mpatches.Patch(color=colour, label=label)
        for _, colour, label in PLDDT_BANDS
    ]
    legend_patches.append(
        mpatches.Patch(color=COLOUR_IF, alpha=0.5, label='Interface region'))
    ax2.legend(handles=legend_patches, fontsize=8, loc='lower right')

    fig.tight_layout()
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
        Normalisation reads `THRESHOLDS[name].green`, the canonical table, which
        is not the ad-hoc threshold dict the notebook still carries inline; the
        numbers here therefore differ from the notebook's until R016 switches its
        summary table over too.
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
# Four implementation defects ARE fixed here, because each is a defect rather
# than a design choice and none of them changes a pixel (R075): the deprecated
# `matplotlib.cm.get_cmap`, the one-component-per-residue payload, the hard-coded
# chain letters, and the structure URL and format being resolved from two
# independent expressions that disagree when `bcifUrl` is present but empty.


# -- view constants ---------------------------------------------------------
# Every colour, cutoff and label the four views use, named once. Section 6 of the
# notebook currently spells all of these inline; naming them here is what makes
# R070-R074 edits to *this list* rather than edits inside the builders.

MVS_VIEW_WIDTH: int = 950
"""Default viewer width in pixels. Fits a notebook cell at the usual zoom."""

MVS_VIEW_HEIGHT: int = 600
"""Default viewer height in pixels."""

MVS_CONTEXT_COLOUR: str = "#BDBDBD"
"""Mid grey. The whole-complex cartoon behind View 3's coloured interface, dark
enough to read as structure and light enough not to compete with the colours."""

MVS_FAINT_COLOUR: str = "#EEEEEE"
"""Near-white. The whole-complex cartoon behind View 4, fainter than
`MVS_CONTEXT_COLOUR` because View 4's three categories are its entire message."""

MVS_VALUE_CMAP: str = "RdYlGn"
"""Colormap for any per-residue score painted onto a structure: red = low,
green = high. Resolved through `matplotlib.colormaps`, never through the
deprecated `matplotlib.cm.get_cmap` the notebook calls (R075).

R072 owns making the mapping legible -- the numbers behind "low" and "high", and
a colour bar or legend. This constant is the seam if the ramp itself changes."""

MVS_PLDDT_BANDS: Tuple[Tuple[float, float, str, str], ...] = (
    (90.0, 100.0, "#1565C0", ">90 (very high)"),
    (70.0, 90.0, "#42A5F5", "70–90 (confident)"),
    (50.0, 70.0, "#FFCA28", "50–70 (low)"),
    (0.0, 50.0, "#EF6C00", "<50 (very low)"),
)
"""pLDDT bands for the 3D view, as `(low, high, colour, label)` tested
`low <= value < high`.

Two differences from `PLDDT_BANDS`, both preserved from the notebook on purpose:

- **The bands are half-open intervals, not one-sided tests.** A residue at
  exactly 100.0 or below 0.0 matches no band and is therefore not drawn at all.
  AFDB writes pLDDT to two decimals in the B-factor column and a residue does
  occasionally reach 100.00.
- **The edges are inclusive-below, not exclusive.** 90.0 lands in the top band
  here and in the second band under `plddt_band_colour`, which tests
  `value > 90`.

Both are visible defects and neither is R014's to fix: this task is a
behaviour-preserving move and a 3D view has no oracle. Reconciling the two
ladders belongs with the M6 pass that gives the view a legend."""

MVS_DISAGREEMENT_THRESHOLD: float = 0.5
"""Per-residue score above which View 4 calls PAE "confident".

**A magic number, and known to be one.** It is unrelated to every cutoff in
`THRESHOLDS`, and R073 replaces it with the R002 ipSAE threshold. Kept at 0.5
here only so that R014 reproduces what the notebook draws today."""

MVS_DISAGREEMENT_CATEGORIES: Tuple[Tuple[str, str, str], ...] = (
    ("agree", "#4CAF50", "PAE+contact agree"),
    ("pae_only", "#2196F3", "PAE confident, no contact"),
    ("contact_only", "#F44336", "Contact, low PAE confidence"),
)
"""View 4's three categories, as `(key, colour, label)`.

The notebook builds exactly this table and then discards the label element with
`_`, so the legend was clearly intended and never drawn. The labels are carried
here so R073 can draw it, and so the category names can be rewritten in terms a
reader can act on rather than as a colour key."""

MVS_VIEW_LABELS: Dict[str, str] = {
    "chain_overview": "View 1: Chain Overview (teal=A, coral=B, amber=interface)",
    "plddt": ("View 2: pLDDT Mapping "
              "(dark blue>90, light blue 70–90, yellow 50–70, orange<50)"),
    "interface_value": ("View 3: Interface ipSAE d0res score "
                        "(red=low, yellow=mid, green=high)"),
    "disagreement": ("View 4: Disagreement (green=PAE+contact, "
                     "blue=PAE confident/no contact, red=contact/low PAE)"),
}
"""The caption `show_mol_view` draws above each view, verbatim from the notebook.

Every one of these is a colour key rather than an explanation, and the chain
letters in View 1's are hard-coded. R070 replaces them with real supporting
text; this dict is where that lands."""

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

def show_mol_view(
    state: Any,
    label: str,
    width: int = MVS_VIEW_WIDTH,
    height: int = MVS_VIEW_HEIGHT,
) -> None:
    """
    Render a MolViewSpec `State` inline, above a bold label.

    The viewer HTML is inlined as a base64 `data:` URI rather than written to a
    file and served, which is what makes it work in PyCharm and in Colab as well
    as in classic Jupyter: none of the three agree on how a notebook-relative
    file URL resolves, and all three render a `data:` iframe.

    The only function in this section that touches `IPython`, so every builder
    stays usable head-lessly. `IPython` is imported lazily here for the same
    reason `molviewspec` is: the module must import outside a notebook.

    Args:
        state:  A MolViewSpec `State`, from any `build_*_view` function.
        label:  Caption drawn above the viewer.
        width:  Iframe width in pixels.
        height: Iframe height in pixels.

    Returns:
        `None`. Displays as a side effect; this is the one function here that does.
    """
    from IPython.display import HTML, IFrame, display

    html = state.molstar_html()
    encoded = base64.b64encode(html.encode()).decode()
    display(HTML(f'<div style="margin:10px 0 4px; font-weight:bold;">{label}</div>'))
    display(IFrame(src=f'data:text/html;base64,{encoded}', width=width, height=height))


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

    Args:
        source:             Structure URL and format, from `resolve_structure_source`.
        chains:             `{chain_id: ChainCoords}`, from `parse_structure`.
        contacts:           Interface contacts naming the ordered pair to draw.
        chain_colours:      `{chain_id: colour}` for the cartoons. Defaults to
                            `CHAIN_COLOURS` assigned in ordered-pair order.
        side_chain_colours: `{chain_id: colour}` for the interface side chains.
                            **This is the R071 seam.** Defaults to
                            `{chain_x: COLOUR_IF}`, i.e. the first chain only,
                            which is what the notebook draws today; R071 passes
                            both chains with a colour each.

    Returns:
        A MolViewSpec `State`. Nothing is rendered; pass it to `show_mol_view`.

    Raises:
        ImportError: If `molviewspec` is not installed.
        KeyError:    If a chain named by `contacts` is absent from `chains`.
    """
    chain_x, chain_y = _ordered_chain_ids(contacts)
    cartoon = (dict(zip((chain_x, chain_y), CHAIN_COLOURS))
               if chain_colours is None else dict(chain_colours))
    side_chains = ({chain_x: COLOUR_IF}
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
    bands: Sequence[Tuple[float, float, str, str]] = MVS_PLDDT_BANDS,
) -> Any:
    """
    View 2: every residue of both chains coloured by its own pLDDT band.

    The control for every other view. A low-confidence interface and a
    high-confidence interface can produce the same contact map, and this is where
    the difference shows.

    A residue in no band gets no component and so is not drawn at all. That is
    the notebook's behaviour and is preserved deliberately: `MVS_PLDDT_BANDS`'
    top band is half-open at 100, so a residue at exactly 100.0 falls through
    every band. See that constant for why this differs from `PLDDT_BANDS`.

    Args:
        source:   Structure URL and format.
        chains:   `{chain_id: ChainCoords}`, from `parse_structure`.
        contacts: Interface contacts naming the ordered pair to draw.
        plddt:    `{chain_id: (n,) array}` overriding the values carried on
                  `ChainCoords`. Pass `PLDDTScores.for_chain(...)` to colour from
                  the pLDDT JSON document instead of the mmCIF B-factor column;
                  the two carry the same numbers.
        bands:    `(low, high, colour, label)` per band, tested as
                  `low <= value < high`. The seam for a band-scheme change.

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

        selected: List[int] = []
        selected_colours: List[str] = []
        for index, value in enumerate(values):
            for low, high, colour, _label in bands:
                if low <= value < high:
                    selected.append(int(res_ids[index]))
                    selected_colours.append(colour)
                    break

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

    **This function is the R072 / R074 seam, and is deliberately score-agnostic.**
    Section 6's View 3 is this called with ipSAE `d0res` per-residue values
    (`ipsae.d0res.forward.values` and `.reverse.values`); R074's two pDockQ2 views
    are the same call with the per-residue mean `ptm` that `compute_pdockq2`
    exposes. Nothing here knows which score it is painting, so a new view is a
    new call, not new colouring code.

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
        indices = np.where(mask)[0]
        add_residue_colours(
            structure,
            chain_id,
            [int(res_ids[i]) for i in indices],
            value_colours(value_array[indices], cmap=cmap, vmin=vmin, vmax=vmax),
            representation=representation,
        )

    return builder.get_state()


# -- View 4: PAE / contact disagreement --------------------------------------

def disagreement_categories(
    values: "Sequence[float] | np.ndarray",
    interface_mask: np.ndarray,
    threshold: float = MVS_DISAGREEMENT_THRESHOLD,
) -> np.ndarray:
    """
    Classify each residue by whether PAE confidence and physical contact agree.

    The three categories are mutually exclusive and do not cover every residue:
    a residue that is neither confident nor in contact falls in none of them and
    is left uncoloured, which is the notebook's behaviour.

    Args:
        values:         `(n,)` per-residue PAE-derived score, `0..1`.
        interface_mask: `(n,)` bool, `True` where the residue touches the partner.
        threshold:      Score above which PAE is called confident. **R073 sources
                        this from `THRESHOLDS` instead of the notebook's magic
                        `0.5`, which is unrelated to any published cutoff.**

    Returns:
        `(n,)` array of `MVS_DISAGREEMENT_CATEGORIES` keys, `''` where none applies.

    Raises:
        ValueError: If the two arrays have different lengths.

    Example
    -------
    >>> import numpy as np
    >>> mask = np.array([True, False, True, False])
    >>> disagreement_categories([0.9, 0.9, 0.1, 0.1], mask).tolist()
    ['agree', 'pae_only', 'contact_only', '']
    """
    value_array = np.asarray(values, dtype=float)
    mask = np.asarray(interface_mask, dtype=bool)
    if value_array.shape[0] != mask.shape[0]:
        raise ValueError(
            f"disagreement_categories got {value_array.shape[0]} values but "
            f"{mask.shape[0]} mask entries."
        )
    confident = value_array > threshold
    categories = np.full(value_array.shape[0], '', dtype=object)
    categories[confident & mask] = 'agree'
    categories[confident & ~mask] = 'pae_only'
    categories[~confident & mask] = 'contact_only'
    return categories


def build_disagreement_view(
    source: StructureSource,
    chains: Mapping[str, ChainCoords],
    contacts: InterfaceContacts,
    values_x: "Sequence[float] | np.ndarray",
    threshold: float = MVS_DISAGREEMENT_THRESHOLD,
    categories: Sequence[Tuple[str, str, str]] = MVS_DISAGREEMENT_CATEGORIES,
    base_colour: str = MVS_FAINT_COLOUR,
    representation: str = "ball_and_stick",
) -> Any:
    """
    View 4: where PAE confidence and physical contact disagree.

    Two independent signals say whether a residue is at the interface -- a CB
    atom within the distance cutoff, and a confident inter-chain PAE -- and this
    view paints the residues where they differ. Agreement is the common case;
    the disagreements are what the notebook exists to explain.

    Note:
        Only `contacts.chain_x` is coloured, which is what the notebook draws.
        That is a defect, and it is **R073's** to fix along with the category
        naming, the legend the notebook computed and threw away, and sourcing
        `threshold` from `THRESHOLDS`. It is left alone here because R014 is a
        behaviour-preserving move.

    Args:
        source:         Structure URL and format.
        chains:         `{chain_id: ChainCoords}`, from `parse_structure`.
        contacts:       Interface contacts naming the ordered pair.
        values_x:       `(nx,)` per-residue score for `contacts.chain_x`,
                        typically `ipsae.d0res.forward.values`.
        threshold:      Confidence cutoff. The R073 seam.
        categories:     `(key, colour, label)` per category, keys matching
                        `disagreement_categories`. The R073 legend seam: the
                        labels are carried here rather than discarded.
        base_colour:    Cartoon colour for the rest of the complex.
        representation: Representation for the coloured residues.

    Returns:
        A MolViewSpec `State`.

    Raises:
        ImportError: If `molviewspec` is not installed.
        ValueError:  If `values_x` does not match `contacts.chain_x`'s length.
    """
    chain_x = contacts.chain_x
    builder, structure = _new_structure(source)

    (structure
     .component()
     .representation(type='cartoon')
     .color(color=base_colour))

    res_ids = _chain_res_ids(chains, chain_x)
    value_array = np.asarray(values_x, dtype=float)
    if value_array.shape[0] != res_ids.shape[0]:
        raise ValueError(
            f"Chain {chain_x!r} has {res_ids.shape[0]} residues but "
            f"{value_array.shape[0]} values."
        )

    assigned = disagreement_categories(value_array, contacts.mask_x, threshold=threshold)
    colour_of = {key: colour for key, colour, _label in categories}

    selected: List[int] = []
    selected_colours: List[str] = []
    for index, key in enumerate(assigned):
        if key in colour_of:
            selected.append(int(res_ids[index]))
            selected_colours.append(colour_of[key])

    add_residue_colours(structure, chain_x, selected, selected_colours,
                        representation=representation)

    return builder.get_state()
