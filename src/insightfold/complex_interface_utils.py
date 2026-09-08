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
scoring primitives (the two `d0` helpers and `ptm_func`), the seven score
functions (R012, landing R003 / R005 / R006 / R007), and the threshold table with
its traffic light and the AFDB joint criterion. Still to come: the plots (R013)
and the MolViewSpec views (R014).

The notebook still carries its own inline copies of the scoring code and still
runs off them; R016 switches the call sites over.
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
