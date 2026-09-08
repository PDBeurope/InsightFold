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
This file is being filled in over several tasks; the sections still empty carry
a note naming the task that populates them. Implemented so far: AFDB access,
structure parsing, PAE / pLDDT parsing and interface detection (R011), the shared
scoring primitives (the two `d0` helpers and `ptm_func`), and the threshold
table with its traffic light and the AFDB joint criterion. Still to come: the
score functions (R012), the plots (R013) and the MolViewSpec views (R014).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Mapping, Optional, Tuple

import numpy as np
import requests

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


@dataclass(frozen=True, eq=False)
class AFDBPrediction:
    """
    The whole prediction-endpoint response, with every per-chain entry preserved.

    The endpoint returns **one entry per chain**, and for the fixtures checked so
    far `entries[0]` is chain *B*, not chain A. Collapsing to `entries[0]` is
    therefore a real bug for heterodimers: it reports one chain's sequence, gene
    and UniProt accession as if they were the complex's. This carrier exists so
    that the fix is a change of *which entry a caller asks for*, not a rewrite of
    the fetch path.

    Attributes:
        accession: The `AF-...` accession that was requested.
        entries:   Every entry the endpoint returned, in API order.
    """

    accession: str
    entries: Tuple[Dict[str, Any], ...]

    @property
    def chain_ids(self) -> Tuple[str, ...]:
        """`chainId` of each entry, in API order (so typically `('B', 'A')`)."""
        return tuple(str(entry.get("chainId", "")) for entry in self.entries)

    def entry_for_chain(self, chain_id: str) -> Dict[str, Any]:
        """
        The metadata entry describing one chain.

        Args:
            chain_id: Structure chain label, e.g. `'A'`.

        Returns:
            That chain's entry.

        Raises:
            KeyError: If no entry carries that `chainId`.
        """
        for entry in self.entries:
            if str(entry.get("chainId", "")) == chain_id:
                return entry
        raise KeyError(
            f"No entry for chain {chain_id!r} in {self.accession}; "
            f"available: {', '.join(self.chain_ids) or '(none)'}."
        )

    @property
    def primary_entry(self) -> Dict[str, Any]:
        """
        The entry used for whole-complex display fields.

        Note:
            Currently the first entry, which reproduces the notebook's existing
            behaviour exactly. That is wrong for heterodimers; see the TODO below.
        """
        # TODO(R020): per-chain selection goes here. Replace this with a lookup
        # keyed on the chain being described -- `entry_for_chain(chain_id)` --
        # and report sequence / geneNames / proteinFullName / uniprotAccession /
        # monomer length once per chain instead of once per complex. The
        # document URLs are unaffected: they are identical across entries.
        if not self.entries:
            raise ValueError(f"No prediction entries returned for {self.accession}.")
        return self.entries[0]

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
        values = {entry[field] for entry in self.entries if entry.get(field)}
        if not values:
            raise KeyError(f"{field!r} is absent from every entry of {self.accession}.")
        if len(values) > 1:
            raise ValueError(
                f"{field!r} differs between chains of {self.accession}: {sorted(values)}. "
                "Whole-complex documents were assumed identical across entries."
            )
        return str(values.pop())

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
