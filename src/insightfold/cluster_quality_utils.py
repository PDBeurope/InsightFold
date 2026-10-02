"""AFDB cluster quality diagnostics.

Implements the spec pack at ``specs/cluster_quality_diagnostic/``. The notebook
``notebooks/cluster_quality_diagnostic.ipynb`` orchestrates; all logic lives here so it can be
tested against pinned fixtures without importing a notebook cell.

Requirement IDs in docstrings refer to ``specs/cluster_quality_diagnostic/requirements.md``.
"""

from __future__ import annotations

import gzip
import http.client
import json
import random
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

# --------------------------------------------------------------------------------------
# Constants (NFR-008: each base URL defined exactly once)
# --------------------------------------------------------------------------------------

CLUSTER_API_BASE = "https://alphafold.ebi.ac.uk/api/workbench/cluster-family"
PREDICTION_API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"

CLUSTER_FLAGS = {"sequence": "AFDB50/MMseqs2", "structure": "AFDB/Foldseek"}

#: Minimum gap between cluster-endpoint requests. Measured 2026-09-21: 0.8 s spacing produced
#: 35 spurious HTTP 500s across 40 accessions; 2.5 s produced zero. The prediction endpoint has
#: no such requirement (45 back-to-back requests, zero failures), so the two policies are
#: deliberately separate and must not be merged (NFR-003).
CLUSTER_REQUEST_SPACING_S = 2.5
PREDICTION_REQUEST_SPACING_S = 0.0

#: AFDB confidence bands (Tunyasuvunakool et al. 2021, Nature 596:590).
PLDDT_BANDS = ((-np.inf, 50.0, "<50"), (50.0, 70.0, "50-70"), (70.0, 90.0, "70-90"), (90.0, np.inf, ">90"))
PLDDT_BAND_COLOURS = {"<50": "#EE8453", "50-70": "#F9DC4D", "70-90": "#7FC9EF", ">90": "#2152CE"}

#: A query outside these percentiles of its own cluster is called an extreme (REQ-004).
EXTREME_LOW_PERCENTILE = 5.0
EXTREME_HIGH_PERCENTILE = 95.0

#: Below this many members a percentile and an IQR are unstable, so the caveat fires (REQ-004).
LOW_N_THRESHOLD = 30

#: The notebook needs at least this many members in the sequence cluster to run (REQ-014).
MIN_CLUSTER_MEMBERS = 3

#: Timeouts. Small JSON responses should fail fast; the streamed files must not, because a
#: tight timeout on a slow link is itself a cause of truncation rather than a guard against it.
API_TIMEOUT_S = 30.0
FILE_TIMEOUT_S = 180.0

#: Transient-failure retry policy. AFDB's file server sends no `Accept-Ranges`, so a partial
#: download cannot be resumed with a Range request and every attempt restarts from zero.
DOWNLOAD_ATTEMPTS = 4
RETRY_BACKOFF_S = 1.5

_ACCESSION_RE = re.compile(r"^[A-Z0-9]{6,10}(-\d+)?$")
_AF_FORM_RE = re.compile(r"^AF-(?P<acc>.+?)-F\d+$")


# --------------------------------------------------------------------------------------
# Errors (REQ-014, REQ-017: every refusal names its condition)
# --------------------------------------------------------------------------------------


class ClusterQualityError(RuntimeError):
    """Base class so a notebook can catch every refusal from this module."""


class AccessionFormatError(ClusterQualityError):
    """The accession is not a bare UniProt accession and could not be normalised."""


class AccessionNotFoundError(ClusterQualityError):
    """AFDB returned 404. See REQ-017 for the required user-facing explanation."""


class ClusterTooSmallError(ClusterQualityError):
    """Fewer members than the analysis needs."""


class ServiceError(ClusterQualityError):
    """A 5xx that survived retries."""


class NetworkError(ClusterQualityError):
    """A transient network failure that survived retries.

    Raised rather than letting `IncompleteRead(163840 bytes read)` reach the user, which names
    neither the file nor what to do about it.
    """


class TruncatedDownloadError(ClusterQualityError):
    """A file arrived short of what it should contain.

    AFDB's file server streams with `Transfer-Encoding: chunked` and sends **no
    `Content-Length`**, so a cut connection cannot be detected from the headers. A truncated
    mmCIF parses perfectly happily and yields a partial structure, so completeness is checked
    against the expected residue count instead.
    """


# --------------------------------------------------------------------------------------
# Accession handling (REQ-001, REQ-001a)
# --------------------------------------------------------------------------------------


def normalise_accession(raw: str) -> tuple[str, str | None]:
    """Return ``(bare_accession, note)``.

    The ``AF-{ACC}-F1`` form is **normalised, not rejected** (REQ-001), and ``note`` carries the
    interpretation so the notebook can print what it did. Isoforms such as ``P04637-2`` are
    accepted and returned unchanged.
    """
    text = (raw or "").strip()
    if not text:
        raise AccessionFormatError("No accession supplied. Set ACCESSION to a bare UniProt accession, for example P69905.")

    note = None
    match = _AF_FORM_RE.match(text)
    if match:
        text = match.group("acc")
        note = f"Interpreted {raw.strip()} as {text}."

    text = text.upper()
    if not _ACCESSION_RE.match(text):
        raise AccessionFormatError(
            f"{raw.strip()!r} is not a recognisable UniProt accession. "
            "Expected a bare accession such as P69905, Q9I1F6 or an isoform such as P04637-2."
        )
    return text, note


def is_isoform(accession: str) -> bool:
    return "-" in accession


def afdb_id(accession: str) -> str:
    """AFDB's own identifier form for a UniProt accession."""
    return f"AF-{accession}-F1"


# --------------------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class QueryMeta:
    accession: str
    description: str
    organism: str
    length: int
    plddt: float
    cif_url: str
    msa_url: str | None
    pae_url: str | None
    latest_version: int | None
    sequence: str


@dataclass
class ClusterMembers:
    """One AFDB cluster. ``flag`` is 'sequence' or 'structure'."""

    accession: str
    flag: str
    cluster_total: int
    table: pd.DataFrame
    sort_was_descending: bool

    def __len__(self) -> int:
        return len(self.table)

    @property
    def plddt(self) -> np.ndarray:
        return self.table["average_plddt"].to_numpy(float)

    def index_of(self, accession: str) -> int | None:
        """Row index of ``accession``, or ``None``.

        Absence is a normal path, not an error: on the structure clustering the query is
        reliably absent because Foldseek clusters AFDB50 representatives (REQ-006), and an
        isoform is absent from the parent's cluster (REQ-001a).
        """
        hits = np.flatnonzero(self.table["accession"].to_numpy() == afdb_id(accession))
        return int(hits[0]) if hits.size else None


@dataclass(frozen=True)
class FamilySummary:
    n: int
    mean: float
    median: float
    sd: float | None
    iqr_width: float
    bands: dict[str, float]
    length_median: float
    length_cv: float | None
    low_n: bool


@dataclass(frozen=True)
class Position:
    present: bool
    index: int | None
    percentile: float | None
    category: str
    iqr_width: float
    low_n: bool


@dataclass(frozen=True)
class ViewTarget:
    role: str  # 'query' | 'best' | 'worst'
    accession: str
    afdb_id: str
    description: str
    plddt: float


@dataclass(frozen=True)
class AlignmentResult:
    target: str
    pct_identity: float
    aligned_columns: int
    query_coverage: float
    query_start: int
    query_end: int
    query_length: int
    target_length: int
    aligned_query: str = field(repr=False, default="")
    aligned_target: str = field(repr=False, default="")

    def summary_line(self) -> str:
        """The mandatory text output (REQ-011a). Must print even when the figure fails."""
        return (
            f"{self.pct_identity:.2f}% identity over {self.aligned_columns} aligned columns, "
            f"covering {self.query_coverage:.1f}% of the query "
            f"(residues {self.query_start}-{self.query_end} of {self.query_length})"
        )


# --------------------------------------------------------------------------------------
# HTTP (NFR-003: the two endpoints have separate rate policies)
# --------------------------------------------------------------------------------------

_last_cluster_request_at = 0.0


def _get(url: str, *, timeout: float = API_TIMEOUT_S, attempts: int = DOWNLOAD_ATTEMPTS,
         pace=None, what: str | None = None) -> bytes:
    """GET with retries on transient failures.

    Retries `IncompleteRead` (a cut chunked stream, the common failure on a slow or unstable
    link), socket timeouts, connection resets, and HTTP 429 and 5xx. Does **not** retry 4xx
    other than 429, which are permanent.

    `pace` is called before every attempt, including retries, so an endpoint's rate policy is
    honoured by the backoff rather than bypassed by it.
    """
    label = what or url
    last: Exception | None = None
    for attempt in range(attempts):
        if pace is not None:
            pace()
        try:
            request = urllib.request.Request(url, headers={"Accept": "*/*"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            # HTTPError subclasses URLError subclasses OSError, so it must be caught first.
            if exc.code != 429 and exc.code < 500:
                raise
            last = exc
        except (http.client.IncompleteRead, http.client.HTTPException, OSError) as exc:
            last = exc
        if attempt < attempts - 1:
            delay = RETRY_BACKOFF_S * (2 ** attempt)
            time.sleep(delay + random.uniform(0.0, 0.4 * delay))
    detail = f"{type(last).__name__}: {last}" if last else "unknown error"
    raise NetworkError(
        f"Could not download {label} after {attempts} attempts ({detail}).\n"
        "This is usually a slow or unstable connection rather than a problem with the data. "
        "The file is streamed without a length header and cannot be resumed, so each attempt "
        "restarts. Re-running the cell will try again."
    )


def _get_cluster(url: str, *, attempts: int = DOWNLOAD_ATTEMPTS) -> bytes:
    """Cluster-endpoint GET. Spacing is enforced before **every** attempt, retries included.

    Cluster responses reach 8.8 MB and are streamed chunked, so they are the most exposed to
    truncation of anything the notebook fetches.
    """
    def pace() -> None:
        global _last_cluster_request_at
        wait = CLUSTER_REQUEST_SPACING_S - (time.monotonic() - _last_cluster_request_at)
        if wait > 0:
            time.sleep(wait)
        _last_cluster_request_at = time.monotonic()

    return _get(url, timeout=FILE_TIMEOUT_S, attempts=attempts, pace=pace,
                what=f"cluster members from {url.rsplit('/', 1)[-1].split('?')[0]}")


def _not_found_message(accession: str) -> str:
    """REQ-017: the entry-point failure must be actionable, not a bare 404."""
    return (
        f"AFDB has no cluster data for {accession}.\n"
        "Two common reasons:\n"
        "  1. The protein is not catalogued in AFDB. Proteins split into multiple fragments, for\n"
        "     example very long ones such as titin, are absent from clustering entirely.\n"
        "  2. The accession has been deleted from UniProtKB.\n"
        "Try an equivalent or closely related protein, for example a well-characterised homologue\n"
        "or the reviewed (SwissProt) entry for the same gene."
    )


# --------------------------------------------------------------------------------------
# Fetching
# --------------------------------------------------------------------------------------


def fetch_prediction(accession: str) -> QueryMeta:
    """Fetch AFDB prediction metadata and sequence. No rate spacing required (NFR-003)."""
    try:
        payload = json.loads(_get(f"{PREDICTION_API_BASE}/{accession}"))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise AccessionNotFoundError(_not_found_message(accession)) from exc
        raise ServiceError(f"AFDB prediction service error for {accession}: {exc}") from exc

    if not payload:
        raise AccessionNotFoundError(_not_found_message(accession))
    entry = payload[0]
    sequence = entry.get("uniprotSequence") or entry.get("sequence") or ""
    bcif = entry.get("bcifUrl") or ""  # may be present but empty; treat empty as absent
    return QueryMeta(
        accession=accession,
        description=entry.get("uniprotDescription", "(no description)"),
        organism=entry.get("organismScientificName", "(unknown organism)"),
        length=int(entry.get("uniprotEnd") or len(sequence) or 0),
        plddt=float(entry.get("globalMetricValue", float("nan"))),
        cif_url=entry.get("cifUrl") or bcif,
        msa_url=entry.get("msaUrl") or None,
        pae_url=entry.get("paeDocUrl") or None,
        latest_version=entry.get("latestVersion"),
        sequence=sequence,
    )


def _cluster_from_payload(accession: str, flag: str, payload: dict) -> ClusterMembers:
    members = payload.get("clusterMembers", {})
    columns = {
        "accession": members.get("afdbAccessions", []),
        "description": members.get("uniprotDescriptions", []),
        "species": members.get("speciesNames", []),
        "sequence_length": members.get("sequenceLength", []),
        "average_plddt": members.get("averagePlddt", []),
        "reviewed": members.get("reviewedStatus", []),
        "reference": members.get("referenceLabel", []),
    }
    lengths = {len(v) for v in columns.values() if v}
    if len(lengths) > 1:
        raise ClusterQualityError(
            f"Malformed cluster response for {accession}: member arrays have differing lengths {sorted(lengths)}."
        )
    table = pd.DataFrame(columns)
    plddt = table["average_plddt"].to_numpy(float) if len(table) else np.array([])
    descending = bool(plddt.size <= 1 or np.all(np.diff(plddt) <= 1e-9))
    if not descending:
        # The API has always returned descending order, but assert rather than assume (REQ-009).
        table = table.sort_values("average_plddt", ascending=False).reset_index(drop=True)
    return ClusterMembers(
        accession=accession,
        flag=flag,
        cluster_total=int(payload.get("clusterTotal", len(table))),
        table=table,
        sort_was_descending=descending,
    )


def fetch_cluster(accession: str, flag: str = "sequence") -> ClusterMembers:
    """Fetch one clustering. ``flag`` is 'sequence' (AFDB50/MMseqs2) or 'structure' (AFDB/Foldseek)."""
    if flag not in CLUSTER_FLAGS:
        raise ValueError(f"flag must be one of {sorted(CLUSTER_FLAGS)}, got {flag!r}")
    # The flag values contain a slash. It is legal unencoded in a query value and the service
    # accepts both, but percent-encoding is what the published example uses and what any proxy
    # in front of the service is guaranteed to pass through unchanged.
    url = (f"{CLUSTER_API_BASE}/{urllib.parse.quote(accession, safe='')}"
           f"?cluster_flag={urllib.parse.quote(CLUSTER_FLAGS[flag], safe='')}")
    try:
        payload = json.loads(_get_cluster(url))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise AccessionNotFoundError(_not_found_message(accession)) from exc
        raise ServiceError(f"AFDB cluster service error for {accession} ({flag}): {exc}") from exc
    return _cluster_from_payload(accession, flag, payload)


def load_pinned_cluster(path, accession: str, flag: str = "sequence") -> ClusterMembers:
    """Load a pinned fixture response. For the offline test harness only; the notebook never uses it."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as handle:
        return _cluster_from_payload(accession, flag, json.load(handle))


# --------------------------------------------------------------------------------------
# Statistics (REQ-003, REQ-004)
# --------------------------------------------------------------------------------------


def band_composition(plddt: Sequence[float]) -> dict[str, float]:
    values = np.asarray(plddt, float)
    if values.size == 0:
        return {name: 0.0 for _, _, name in PLDDT_BANDS}
    return {name: float(((values >= low) & (values < high)).mean()) for low, high, name in PLDDT_BANDS}


def summarise_family(cluster: ClusterMembers) -> FamilySummary:
    """REQ-003. Statistics describe the cluster the API returned, not a subset."""
    plddt = cluster.plddt
    lengths = cluster.table["sequence_length"].to_numpy(float)
    n = int(plddt.size)
    return FamilySummary(
        n=n,
        mean=float(plddt.mean()),
        median=float(np.median(plddt)),
        sd=float(plddt.std(ddof=1)) if n > 1 else None,
        iqr_width=float(np.percentile(plddt, 75) - np.percentile(plddt, 25)),
        bands=band_composition(plddt),
        length_median=float(np.median(lengths)),
        length_cv=float(lengths.std(ddof=1) / lengths.mean()) if n > 1 and lengths.mean() else None,
        low_n=n < LOW_N_THRESHOLD,
    )


def locate_query(cluster: ClusterMembers, accession: str) -> Position:
    """REQ-004. Position is reported as a **category**; the percentile is secondary.

    The category is relative to this cluster's own distribution, never a global cutoff, because
    the median interquartile width across clusters is about 4.2 pLDDT points and a rank within
    that band carries almost no information.
    """
    index = cluster.index_of(accession)
    plddt = cluster.plddt
    iqr_width = float(np.percentile(plddt, 75) - np.percentile(plddt, 25)) if plddt.size else 0.0
    low_n = len(cluster) < LOW_N_THRESHOLD

    if index is None:
        return Position(False, None, None, "not a member of this cluster", iqr_width, low_n)

    value = float(plddt[index])
    percentile = float((plddt < value).mean() * 100.0)
    low_cut = float(np.percentile(plddt, EXTREME_LOW_PERCENTILE))
    high_cut = float(np.percentile(plddt, EXTREME_HIGH_PERCENTILE))
    if value <= low_cut:
        category = "low extreme"
    elif value >= high_cut:
        category = "high extreme"
    else:
        category = "typical"
    return Position(True, index, percentile, category, iqr_width, low_n)


def length_deviation(cluster: ClusterMembers) -> pd.Series:
    """Absolute deviation from the cluster's median length, as a fraction."""
    lengths = cluster.table["sequence_length"].to_numpy(float)
    median = float(np.median(lengths))
    return pd.Series(np.abs(lengths - median) / median, index=cluster.table.index, name="length_deviation")


# --------------------------------------------------------------------------------------
# Extremes and the collapse rule (REQ-009, REQ-010)
# --------------------------------------------------------------------------------------


def identify_view_targets(cluster: ClusterMembers, query: QueryMeta) -> list[ViewTarget]:
    """Return 3 targets, or 2 when the query is itself an extreme (REQ-010).

    Computed **once** here and passed forward, so the 3D views and the alignments can never
    disagree about which models they are showing.
    """
    if len(cluster) < MIN_CLUSTER_MEMBERS:
        raise ClusterTooSmallError(
            f"The sequence cluster for {query.accession} has {len(cluster)} member(s). "
            f"At least {MIN_CLUSTER_MEMBERS} are needed to show a distribution and its extremes. "
            "Try a protein with a larger family."
        )

    def target(row_index: int, role: str) -> ViewTarget:
        row = cluster.table.iloc[row_index]
        return ViewTarget(
            role=role,
            accession=re.sub(r"^AF-|-F\d+$", "", str(row["accession"])),
            afdb_id=str(row["accession"]),
            description=str(row["description"]),
            plddt=float(row["average_plddt"]),
        )

    best, worst = target(0, "best"), target(len(cluster) - 1, "worst")
    query_target = ViewTarget("query", query.accession, afdb_id(query.accession), query.description, query.plddt)

    if best.accession == query.accession:
        return [query_target, worst]
    if worst.accession == query.accession:
        return [best, query_target]
    return [query_target, best, worst]


# --------------------------------------------------------------------------------------
# Alignment (REQ-011, REQ-011a)
# --------------------------------------------------------------------------------------


def align_pair(query_sequence: str, target_sequence: str, target_label: str) -> AlignmentResult:
    """Local alignment, reported as identity **with aligned length and query coverage**.

    Local rather than global: global identity penalises length and domain-architecture mismatch,
    which is exactly the situation this notebook exists to surface, and nothing bounds the
    query-to-member length ratio because AFDB50's overlap criterion holds against the cluster
    representative rather than between members.

    Coverage is not decoration. Measured on the fixture set, identity alone ranks a functionally
    divergent member (26.73%) *above* a legitimate distant homologue (25.09%), while coverage
    separates them cleanly (21.3% against 74.6%).
    """
    from Bio import Align  # imported lazily so the module imports without biopython
    from Bio.Align import substitution_matrices

    if not query_sequence or not target_sequence:
        raise ClusterQualityError(f"No sequence available for {target_label}; cannot align.")

    # BLOSUM62 with -11/-1 is the BLASTP local default. Do not tune these to make a case
    # separate; coverage is parameter-sensitive and tuning it would be fitting to the fixtures.
    aligner = Align.PairwiseAligner(
        mode="local",
        open_gap_score=-11,
        extend_gap_score=-1,
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
    )
    alignment = aligner.align(query_sequence, target_sequence)[0]
    aligned_query, aligned_target = str(alignment[0]), str(alignment[1])
    query_blocks = alignment.aligned[0]
    query_start = int(query_blocks[0][0]) + 1
    query_end = int(query_blocks[-1][-1])
    columns = len(aligned_query)
    pairs = [(a, b) for a, b in zip(aligned_query, aligned_target) if a != "-" and b != "-"]
    identities = sum(a == b for a, b in pairs)
    return AlignmentResult(
        target=target_label,
        pct_identity=100.0 * identities / columns if columns else 0.0,
        aligned_columns=columns,
        query_coverage=100.0 * len(pairs) / len(query_sequence),
        query_start=query_start,
        query_end=query_end,
        query_length=len(query_sequence),
        target_length=len(target_sequence),
        aligned_query=aligned_query,
        aligned_target=aligned_target,
    )


# --------------------------------------------------------------------------------------
# Plotting. Every function returns a bare Figure built without pyplot (NFR-005).
# --------------------------------------------------------------------------------------


def _new_figure(width: float, height: float) -> tuple[Figure, "object"]:
    figure = Figure(figsize=(width, height), layout="constrained")
    return figure, figure.add_subplot(111)


def shared_plddt_axis_limits(*clusters, pad: float = 20.0, step: float = 10.0) -> tuple[float, float]:
    """One x-axis range covering every supplied cluster, so the rainclouds are comparable.

    Letting each figure autoscale is actively misleading: two distributions drawn side by side
    on different axes look similar when they are not. The range is computed from the pooled
    minimum and maximum across all clusters, padded by ``pad``, rounded outward to the nearest
    ``step``, then clamped to the 0 to 100 range pLDDT is defined on.

    Rounding is half-up rather than Python's banker's rounding, so 45 gives 50 rather than 40.
    The padding always exceeds half a step, so the rounded bound can never clip a data point:
    the lower bound is at worst ``min - pad + step/2``, which is still below ``min``.
    """
    import math

    pools = [c.plddt for c in clusters if c is not None and len(c)]
    if not pools:
        return 0.0, 100.0
    values = np.concatenate(pools)
    low = math.floor((float(values.min()) - pad) / step + 0.5) * step
    high = math.floor((float(values.max()) + pad) / step + 0.5) * step
    return max(0.0, low), min(100.0, high)


def _half_violin(axis, values, position, width=0.85, colour="#BBBBBB"):
    """Upper half of a horizontal violin at ``position``. The 'cloud' of a raincloud."""
    parts = axis.violinplot([values], positions=[position], vert=False,
                            showextrema=False, showmedians=False, widths=width)
    for body in parts["bodies"]:
        vertices = body.get_paths()[0].vertices
        vertices[:, 1] = np.clip(vertices[:, 1], position, np.inf)
        body.set_facecolor(colour)
        body.set_edgecolor("none")
        body.set_alpha(0.55)


def _rain(axis, values, position, colours, rng, max_points=2000, spread=0.19, offset=0.05):
    """Jittered points below the violin, coloured by pLDDT band. The 'rain'.

    Subsampled above ``max_points``: a 93,793-member cluster draws as a solid block otherwise,
    and the figure takes longer than the API call that produced it.
    """
    values = np.asarray(values, float)
    colours = np.asarray(colours, object)
    n = values.size
    if n > max_points:
        pick = rng.choice(n, size=max_points, replace=False)
        values, colours = values[pick], colours[pick]
    y = position - offset - rng.uniform(0.0, spread, size=values.size)
    axis.scatter(values, y, c=list(colours), s=7, alpha=0.55, linewidths=0)
    return n > max_points


def _band_colours(plddt) -> list[str]:
    return [plddt_colour(float(v)) for v in plddt]


def plot_cluster_raincloud(
    cluster: "ClusterMembers",
    label: str,
    summary: "FamilySummary | None" = None,
    position: "Position | None" = None,
    query_label: str = "your protein",
    seed: int = 0,
    xlim: "tuple[float, float] | None" = None,
) -> Figure:
    """Raincloud of one cluster's pLDDT: half violin, jittered points by band, boxplot.

    One clustering per figure. The two clusterings are **never** drawn on the same Axes: a
    structure-cluster point is one sequence family scored by its best member, so laying the two
    side by side invites a comparison that is not like for like (REQ-005, REQ-006).
    """
    rng = np.random.default_rng(seed)
    figure, axis = _new_figure(9.0, 3.6)
    plddt = cluster.plddt
    pos = 0.0

    _half_violin(axis, plddt, pos)
    subsampled = _rain(axis, plddt, pos, _band_colours(plddt), rng)
    axis.boxplot([plddt], positions=[pos - 0.34], vert=False, widths=0.07, showfliers=False,
                 showcaps=False, medianprops={"linewidth": 1.6, "color": "#333333"},
                 boxprops={"linewidth": 1.2, "color": "#555555"},
                 whiskerprops={"linewidth": 1.2, "color": "#555555"})

    if position is not None and position.present:
        value = float(plddt[position.index])
        axis.axvline(value, color="#C0392B", lw=2.0, zorder=5,
                     label=f"{query_label} {value:.1f} ({position.category})")

    from matplotlib.patches import Patch
    # Ascending, so the legend reads left to right in increasing confidence.
    handles = [Patch(facecolor=PLDDT_BAND_COLOURS[name], label=name)
               for _, _, name in PLDDT_BANDS]
    if position is not None and position.present:
        handles = axis.get_legend_handles_labels()[0] + handles
    # Below the axes: anywhere inside collides with either the violin's tail or the rain.
    axis.legend(handles=handles, frameon=False, fontsize=8, ncol=len(handles),
                loc="upper center", bbox_to_anchor=(0.5, -0.22), handlelength=1.2,
                columnspacing=1.2)

    axis.set_yticks([])
    axis.set_ylim(pos - 0.46, pos + 0.66)
    if xlim is not None:
        axis.set_xlim(*xlim)
    axis.set_xlabel("Average pLDDT")
    title = f"{label}, n = {len(cluster):,}"
    if summary is not None:
        title += f"   (median {summary.median:.1f}, p25 to p75 {summary.iqr_width:.2f} wide)"
    if subsampled:
        title += "\npoints subsampled to 2,000; violin and box use every member"
    axis.set_title(title, fontsize=10)
    return figure


def plot_family_distribution(cluster, summary, position=None, query_label="your protein",
                             xlim=None) -> Figure:
    """The sequence-family raincloud (REQ-V1)."""
    return plot_cluster_raincloud(cluster, "Sequence family", summary, position, query_label,
                                  xlim=xlim)


def plot_structure_cluster_distribution(cluster: "ClusterMembers", xlim=None) -> Figure:
    """The structure-cluster raincloud (REQ-V2).

    No query marker: measured across 8 accessions the structure-cluster members are disjoint
    from the sequence-cluster members, and the query is normally in neither. Marking it would
    place an individual on a distribution it is not part of.
    """
    figure = plot_cluster_raincloud(
        cluster, "Structurally similar proteins", xlim=xlim)
    figure.axes[0].set_xlabel("Average pLDDT of each structurally similar protein")
    return figure


def plot_length_vs_plddt(cluster: "ClusterMembers", label: str = "sequence family",
                         colour_by_band: bool = True, max_points: int = 6000,
                         seed: int = 0) -> Figure:
    """Length against confidence, as a joint density with marginals.

    Points are coloured by AFDB pLDDT band, matching the raincloud and the 3D views. Because
    the y axis *is* pLDDT, that colouring is a function of height and reads as horizontal
    stripes: it carries no information the axis does not already, and is there for continuity
    with the other figures and to make the band boundaries legible without reading the axis.
    Pass ``colour_by_band=False`` for a plain scatter.

    The 2D density underneath is deliberately monochrome. A `mako` density under band-coloured
    points produces two competing colour scales in one panel, and neither reads.

    Built on a bare Figure with gridspec rather than ``seaborn.JointGrid``, which would create
    its own pyplot figure and leak into the registry (NFR-005).
    """
    import seaborn as sns

    rng = np.random.default_rng(seed)
    lengths = cluster.table["sequence_length"].to_numpy(float)
    plddt = cluster.plddt

    figure = Figure(figsize=(8.0, 7.0), layout="constrained")
    if colour_by_band:
        # No colourbar: a 2D density under band-coloured points is invisible, and its
        # colourbar is then dead space. The KDE contours carry the density instead.
        grid = figure.add_gridspec(2, 2, width_ratios=(7, 1.4), height_ratios=(1.4, 7),
                                   wspace=0.05, hspace=0.05)
        cax = None
    else:
        grid = figure.add_gridspec(2, 3, width_ratios=(7, 1.4, 0.35), height_ratios=(1.4, 7),
                                   wspace=0.05, hspace=0.05)
        cax = figure.add_subplot(grid[1, 2])
    main = figure.add_subplot(grid[1, 0])
    top = figure.add_subplot(grid[0, 0], sharex=main)
    right = figure.add_subplot(grid[1, 1], sharey=main)

    if not colour_by_band:
        sns.histplot(x=lengths, y=plddt, bins=60, pthresh=0.05, cmap="mako",
                     ax=main, cbar=True, cbar_ax=cax)

    x, y = lengths, plddt
    if x.size > max_points:
        pick = rng.choice(x.size, size=max_points, replace=False)
        x, y = x[pick], y[pick]
    if colour_by_band:
        main.scatter(x, y, c=[plddt_colour(float(v)) for v in y], s=7, alpha=0.7, linewidths=0)
    else:
        main.scatter(x, y, s=4, color=".15", alpha=0.35, linewidths=0)

    if len(cluster) >= 50:
        try:
            sns.kdeplot(x=lengths, y=plddt, levels=8,
                        color="#B00020" if not colour_by_band else "#222222",
                        linewidths=0.9, ax=main)
        except Exception:
            pass

    sns.histplot(x=lengths, bins=60, element="step", color="#03012d", ax=top)
    if colour_by_band:
        counts, edges = np.histogram(plddt, bins=60)
        centres = (edges[:-1] + edges[1:]) / 2.0
        right.barh(centres, counts, height=np.diff(edges),
                   color=[plddt_colour(float(v)) for v in centres], linewidth=0)
    else:
        sns.histplot(y=plddt, bins=60, element="step", color="#03012d", ax=right)

    top.set_xlabel(""); top.set_ylabel(""); top.tick_params(labelbottom=False)
    right.set_xlabel(""); right.set_ylabel(""); right.tick_params(labelleft=False)
    for spine in ("top", "right"):
        top.spines[spine].set_visible(False); right.spines[spine].set_visible(False)
    main.set_xlabel("Sequence length (aa)", fontsize=12)
    main.set_ylabel("Average pLDDT", fontsize=12)
    if cax is not None:
        cax.set_ylabel("Members per bin", fontsize=9)
    subtitle = f", points subsampled to {max_points:,}" if len(cluster) > max_points else ""
    top.set_title(f"Length against confidence, {label}, n = {len(cluster):,}{subtitle}", fontsize=11)
    if colour_by_band:
        from matplotlib.patches import Patch
        # Descending here, unlike the horizontal legends elsewhere: this one is vertical and
        # sits beside a pLDDT y axis, so reading top to bottom should track the axis.
        main.legend(handles=[Patch(facecolor=PLDDT_BAND_COLOURS[n], label=n)
                             for _, _, n in reversed(PLDDT_BANDS)],
                    title="pLDDT band", frameon=True, framealpha=0.85,
                    loc="lower right", fontsize=8, title_fontsize=8)
    return figure


def plot_plddt_profile(residue_numbers: Sequence[int], plddt: Sequence[float], label: str) -> Figure:
    """REQ-V8. A chain average is largely a function of disorder content, so the profile is what
    distinguishes a uniformly poor model from an ordered domain plus a disordered tail."""
    figure, axis = _new_figure(9.0, 2.6)
    x = np.asarray(residue_numbers, float)
    y = np.asarray(plddt, float)
    axis.plot(x, y, lw=1.1, color="#333333")
    for low, high, name in PLDDT_BANDS:
        axis.axhspan(max(low, 0.0), min(high, 100.0), color=PLDDT_BAND_COLOURS[name], alpha=0.16, zorder=0)
    axis.set_ylim(0, 100)
    axis.set_xlabel("Residue")
    axis.set_ylabel("pLDDT")
    axis.set_title(f"Per-residue confidence, {label} (mean {np.nanmean(y):.1f})")
    return figure


def parse_plddt_from_cif(text: str, expected_residues: int | None = None) -> tuple[list[int], list[float]]:
    """Per-residue pLDDT from an AFDB mmCIF, using CA atoms' B-factor column.

    A blank line inside the ``_atom_site`` loop is **not** a terminator; the loop ends at the
    next data name, ``loop_`` or ``#``.
    """
    lines = text.splitlines()
    columns: list[str] = []
    in_loop = False
    residues: dict[int, float] = {}
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("_atom_site."):
            columns.append(stripped.split(".", 1)[1].split()[0].lower())
            in_loop = True
            continue
        if not in_loop:
            continue
        if stripped.startswith(("_", "loop_", "#")):
            break
        if not stripped:
            continue
        fields = stripped.split()
        if len(fields) != len(columns):
            continue
        row = dict(zip(columns, fields))
        if row.get("group_pdb") != "ATOM" or row.get("label_atom_id") != "CA":
            continue
        seq_id = row.get("label_seq_id", ".")
        if seq_id in (".", "?"):
            continue
        model = row.get("pdbx_pdb_model_num")
        if model is not None and model != "1":
            continue
        try:
            residues[int(seq_id)] = float(row["b_iso_or_equiv"])
        except (KeyError, ValueError):
            continue
    ordered = sorted(residues)
    if expected_residues is not None and len(ordered) < expected_residues:
        raise TruncatedDownloadError(
            f"The structure file parsed to {len(ordered)} residues but the prediction is "
            f"{expected_residues} residues long, so the download was cut short. "
            "The file is streamed without a length header, so this cannot be detected from "
            "the response itself and is checked here instead."
        )
    return ordered, [residues[i] for i in ordered]


# --------------------------------------------------------------------------------------
# MolViewSpec views (REQ-010, REQ-V6)
#
# This section follows the pattern established in `complex_interface_utils.py`. Three rules
# carried over from there, each for a reason:
#
#   1. `molviewspec` is imported **lazily**, inside `_require_molviewspec()`, never at module
#      top level. The module must stay importable and every non-3D function must stay usable
#      on a machine without the package.
#   2. Builders return a MolViewSpec `State` and do **not** render. `show_mol_view` is the
#      only function here that touches IPython, so builders stay usable head-lessly.
#   3. Contiguous same-colour residues collapse into **one** ranged `ComponentExpression`.
#      One component per residue is what the original notebook did; collapsing cut one view
#      from 344 components to 12.
#
# The viewer HTML is inlined as a base64 `data:` URI rather than served from a file, because
# PyCharm, Colab and classic Jupyter disagree about how a notebook-relative file URL resolves
# and all three render a `data:` iframe.
# --------------------------------------------------------------------------------------

MVS_VIEW_WIDTH = "100%"
MVS_VIEW_HEIGHT = 520

MOLVIEWSPEC_MISSING_MESSAGE = (
    "molviewspec is not installed, so the 3D views are unavailable. "
    "Install it with `pip install molviewspec`; every other section of this "
    "notebook runs without it."
)


def molviewspec_available() -> bool:
    """Whether the optional `molviewspec` dependency can be imported."""
    try:
        import molviewspec  # noqa: F401
    except Exception:
        return False
    return True


def _require_molviewspec():
    try:
        import molviewspec as mvs
    except Exception as exc:  # pragma: no cover
        raise ImportError(MOLVIEWSPEC_MISSING_MESSAGE) from exc
    return mvs


@dataclass(frozen=True)
class StructureSource:
    url: str
    format: str


def resolve_structure_source(cif_url: str, bcif_url: str | None = None, prefer_binary: bool = True) -> StructureSource:
    """Decide URL and parse format together.

    A field present but **empty** is treated as absent. An empty `bcifUrl` labelled `bcif` is
    a silently blank viewer, which is worse than no viewer at all.
    """
    candidates = [(bcif_url, "bcif"), (cif_url, "mmcif")]
    if not prefer_binary:
        candidates.reverse()
    for url, fmt in candidates:
        if url:
            return StructureSource(url=url, format=fmt)
    raise ValueError(
        "No usable structure URL. Mol* downloads the structure itself and cannot be handed "
        "already-parsed text, so the 3D view genuinely cannot run for this model."
    )


def plddt_colour(value: float) -> str:
    """The AFDB band colour for one pLDDT value."""
    for low, high, name in PLDDT_BANDS:
        if low <= value < high:
            return PLDDT_BAND_COLOURS[name]
    return PLDDT_BAND_COLOURS[">90"]


@dataclass(frozen=True)
class ColourRun:
    beg: int
    end: int
    colour: str


def colour_runs(res_ids: Sequence[int], colours: Sequence[str]) -> list[ColourRun]:
    """Collapse a per-residue colour list into contiguous same-colour runs.

    A residue joins the previous run only when it is the same colour **and** its
    `label_seq_id` is exactly one more than the previous residue's. Requiring both is what
    makes the collapse safe: a ranged `ComponentExpression` covers every residue between its
    endpoints, so merging across a gap in the numbering would colour residues never in the
    input.
    """
    ids = [int(r) for r in res_ids]
    cols = list(colours)
    if len(ids) != len(cols):
        raise ValueError(f"colour_runs got {len(ids)} residue ids but {len(cols)} colours; they must be parallel.")
    runs: list[ColourRun] = []
    for res_id, colour in zip(ids, cols):
        if runs and colour == runs[-1].colour and res_id == runs[-1].end + 1:
            runs[-1] = ColourRun(runs[-1].beg, res_id, colour)
        else:
            runs.append(ColourRun(res_id, res_id, colour))
    return runs


def build_plddt_view(source: StructureSource, res_ids: Sequence[int], plddt: Sequence[float],
                     chain_id: str = "A", base_colour: str = "#DDDDDD"):
    """One model, cartoon, coloured by AFDB pLDDT confidence band. Returns a `State`.

    Colour is applied as repeated ``.color(color=..., selector=[...])`` calls on a **single**
    cartoon representation, which is the idiomatic MolViewSpec pattern and the reason this
    works at all. The earlier version created one extra `representation` per colour run on top
    of a base representation of the same atoms: Mol* then renders overlapping cartoon geometry
    for every residue, the copies z-fight, and the base colour wins, so the model appears
    uniformly grey. One representation with several scoped colour nodes has no such overlap.

    Runs are grouped by band, so a 142-residue model emits at most five colour nodes rather
    than one per run.
    """
    mvs = _require_molviewspec()
    builder = mvs.create_builder()
    structure = builder.download(url=source.url).parse(format=source.format).model_structure()
    representation = structure.component(selector="polymer").representation(type="cartoon")
    # Base colour first: anything not covered by a run (a gap in the numbering, a ligand)
    # keeps this rather than inheriting whatever Mol* defaults to.
    representation.color(color=base_colour)

    grouped: dict[str, list] = {}
    for run in colour_runs(res_ids, [plddt_colour(float(v)) for v in plddt]):
        grouped.setdefault(run.colour, []).append(
            mvs.ComponentExpression(label_asym_id=chain_id,
                                    beg_label_seq_id=run.beg, end_label_seq_id=run.end))
    for colour, expressions in grouped.items():
        representation.color(color=colour, selector=expressions)
    return builder.get_state()


def plddt_band_legend_html() -> str:
    """A small inline legend, so the 3D colours are readable without cross-referencing.

    Band names are HTML-escaped. `<50` is a literal `<` in interpolated markup: a browser reads
    it as the start of a tag and silently swallows the label and whatever follows. Bands run in
    ascending order so the legend reads left to right in increasing confidence.
    """
    import html as _html

    swatches = "".join(
        f'<span style="display:inline-flex; align-items:center; margin-right:14px;">'
        f'<span style="width:13px; height:13px; background:{PLDDT_BAND_COLOURS[name]}; '
        f'display:inline-block; margin-right:5px; border:1px solid #999;"></span>'
        f'{_html.escape(name)}</span>'
        for _, _, name in PLDDT_BANDS)
    return (f'<div style="font-size:12px; color:#333; margin:2px 0 8px;">'
            f'<b>pLDDT band:</b> {swatches}</div>')


def _viewer_iframe(viewer_html: str, style: str) -> str:
    """Embed a Mol* document in an iframe via `srcdoc` rather than a `data:` URI.

    `complex_interface_utils.py` uses `src="data:text/html;base64,..."`, which works in
    JupyterLab and Colab. PyCharm's renderer treats a `data:` URI as a navigation to a separate
    document and hands it to the **system browser**, so running the notebook there opens one
    Chrome window per viewer on top of rendering it inline.

    `srcdoc` carries the document in an attribute, so there is no navigation for a host to
    delegate and the content renders in place. Scripts and WebGL still run, because a `srcdoc`
    frame without a `sandbox` attribute inherits the parent's origin.

    Only `&` and `"` need escaping for an attribute value, and `&` must go first.
    """
    document = viewer_html.replace("&", "&amp;").replace('"', "&quot;")
    return f'<iframe srcdoc="{document}" style="{style}" allowfullscreen></iframe>'


def _css_length(value) -> str:
    return value if isinstance(value, str) else f"{value:g}px"


def mol_view_html(state, label: str, width=MVS_VIEW_WIDTH, height=MVS_VIEW_HEIGHT) -> str:
    """The label-plus-iframe markup. Returned, not shown, so it can be asserted on."""
    import html as _html

    style = (f"width:{_css_length(width)}; height:{_css_length(height)}; "
             "display:block; border:0;")
    return (
        f'<div style="margin:10px 0 4px; font-weight:bold;">{_html.escape(label)}</div>'
        + _viewer_iframe(state.molstar_html(), style)
    )


def show_mol_view(state, label: str, width=MVS_VIEW_WIDTH, height=MVS_VIEW_HEIGHT) -> None:
    """Render a `State` inline above a bold label. The only function here that touches IPython."""
    from IPython.display import HTML, display

    display(HTML(mol_view_html(state, label, width=width, height=height)))


def fetch_structure_text(url: str, label: str | None = None) -> str:
    """Download an mmCIF as text, for the per-residue pLDDT profile (REQ-010a)."""
    return _get(url, timeout=FILE_TIMEOUT_S,
                what=label or f"structure {url.rsplit('/', 1)[-1]}").decode("utf-8", errors="replace")


# --------------------------------------------------------------------------------------
# PAE (REQ-V10). Reuses `parse_pae` from `complex_interface_utils`, which is this repo's
# canonical PAE document parser. Per D8 that module is now the shared home for PAE parsing:
# this is its second consumer, so the function is used rather than reimplemented.
#
# Monomer PAE documents omit the `chains` array that a dimer document carries, so the
# sequence length is supplied as `fallback_lengths`.
# --------------------------------------------------------------------------------------

PAE_CMAP = "Greens_r"
"""Dark green = low PAE = confident, and legible on white. Matches the repo default.
`viridis` is the colourblind-safe alternative."""


def fetch_pae(meta: QueryMeta | str, sequence_length: int | None = None):
    """Download and parse one model's PAE matrix. Returns a `PAEMatrix`.

    `meta` may be a `QueryMeta` or a bare accession; in the latter case the prediction
    endpoint is called to find `paeDocUrl`.
    """
    from insightfold.complex_interface_utils import parse_pae

    if isinstance(meta, str):
        payload = json.loads(_get(f"{PREDICTION_API_BASE}/{meta}"))
        if not payload:
            raise AccessionNotFoundError(_not_found_message(meta))
        entry = payload[0]
        url = entry.get("paeDocUrl")
        length = sequence_length or len(entry.get("uniprotSequence") or entry.get("sequence") or "")
    else:
        url = getattr(meta, "pae_url", None)
        length = sequence_length or meta.length or len(meta.sequence)
    if not url:
        raise ClusterQualityError("This prediction carries no paeDocUrl, so no PAE matrix is available.")
    raw = _get(url, timeout=FILE_TIMEOUT_S, what=f"PAE matrix {url.rsplit(chr(47), 1)[-1]}")
    return parse_pae(json.loads(raw), fallback_lengths={"A": int(length)})


def plot_pae_matrix(pae, label: str, cmap: str | None = None) -> Figure:
    """PAE heatmap for one model. `aspect='equal'` so the matrix is never rendered skewed."""
    figure = Figure(figsize=(4.6, 4.3), layout="constrained")
    axis = figure.add_subplot(111)
    image = axis.imshow(pae.matrix, cmap=cmap or PAE_CMAP, vmin=0.0, vmax=pae.max_pae,
                        origin="upper", aspect="equal", interpolation="nearest")
    axis.set_xlabel("Scored residue")
    axis.set_ylabel("Aligned residue")
    axis.set_title(f"Predicted aligned error\n{label}", fontsize=9)
    bar = figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    bar.set_label("PAE (Angstrom)", fontsize=8)
    return figure


def _figure_to_data_uri(figure: Figure, dpi: int = 110) -> str:
    import base64
    import io

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight")
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode()


def show_model_panel(state, pae_figure: Figure | None, label: str,
                     height: int = MVS_VIEW_HEIGHT) -> None:
    """One model: the Mol\\* viewer and its PAE matrix, **side by side**.

    The viewer is an iframe and the PAE is a matplotlib Figure, so they are composed into one
    flex row with the Figure inlined as a PNG data URI. Falls back to the viewer alone when no
    PAE is available.
    """
    from IPython.display import HTML, display

    import html as _html

    viewer = _viewer_iframe(state.molstar_html(),
                            f"flex:1 1 60%; height:{height}px; display:block; border:0;")
    if pae_figure is not None:
        panel = (f'<img class="cq-panel-img" src="{_figure_to_data_uri(pae_figure)}" '
                 f'style="flex:0 1 38%; height:{height}px; object-fit:contain;" />')
    else:
        panel = ('<div style="flex:0 1 38%; display:flex; align-items:center; '
                 'justify-content:center; color:#777; font-size:13px;">PAE unavailable</div>')
    display(HTML(
        f'<div style="margin:12px 0 4px; font-weight:bold;">{_html.escape(label)}</div>'
        f'<div style="display:flex; gap:12px; align-items:stretch;">{viewer}{panel}</div>'
    ))


# --------------------------------------------------------------------------------------
# Notebook display style
# --------------------------------------------------------------------------------------


def notebook_display_css(figure_dpi: int = 144) -> str:
    """CSS that makes matplotlib output fill the notebook column and reflow with the window.

    A matplotlib figure is a raster image at a fixed pixel size, so it does not resize on its
    own. Setting `width:100%` with `height:auto` lets the browser scale it to whatever the
    column is, in both directions, and the aspect ratio is preserved. Rendering at a higher dpi
    keeps it from softening when the window is wide.

    The selectors are deliberately narrow. `:not(.cq-panel-img)` exempts the PAE images inside
    the Mol\* panels, which are laid out by their own flex rules and would be stretched out of
    shape by a blanket `width:100%`. Several selector spellings are included because
    JupyterLab, Notebook 7, classic Jupyter and PyCharm each wrap output images differently.
    """
    return f"""<style>
.jp-RenderedImage > img:not(.cq-panel-img),
.jp-OutputArea-output > img:not(.cq-panel-img),
div.output_area .output_png img:not(.cq-panel-img),
div.output_subarea > img:not(.cq-panel-img) {{
    width: 100% !important;
    height: auto !important;
    max-width: none !important;
}}
.jp-Cell-outputWrapper, .jp-OutputArea-output, div.output_subarea {{ max-width: none !important; }}
.jp-RenderedHTMLCommon, .jp-MarkdownOutput {{ max-width: none !important; }}
</style><!-- figure dpi {figure_dpi} -->"""


def apply_notebook_style(figure_dpi: int = 144) -> None:
    """Make figures full width and crisp. Explicit, never applied on import.

    Importing this module must not mutate global matplotlib state, so the dpi change lives
    here rather than at module level.
    """
    import matplotlib
    from IPython.display import HTML, display

    matplotlib.rcParams["figure.dpi"] = figure_dpi
    matplotlib.rcParams["savefig.dpi"] = figure_dpi
    display(HTML(notebook_display_css(figure_dpi)))


def fetch_model_plddt(meta: QueryMeta, attempts: int = 3) -> tuple[list[int], list[float]]:
    """Download a model's mmCIF and return its per-residue pLDDT, verified complete.

    Two layers, because they catch different failures. `_get` retries a **cut connection**,
    which raises. This retries a download that arrived intact-looking but **short**, which
    does not raise on its own: a truncated mmCIF parses happily and yields a partial
    structure. Measured on `AF-O15552-F1`, 160 KiB of a 311 KB file parsed to 150 residues of
    330 with no error at all.
    """
    _, res_ids, plddt = fetch_verified_structure(meta, attempts=attempts)
    return res_ids, plddt


def fetch_verified_structure(meta: QueryMeta, attempts: int = 3
                             ) -> tuple[str, list[int], list[float]]:
    """The mmCIF **text** alongside its per-residue pLDDT, verified complete.

    The superposition needs the text itself, for coordinates and for the secondary structure
    records, and it needs it whole: a short file parses without complaint into a partial chain,
    and TM-align would then fit that partial chain and report a confident transform for it.
    """
    expected = meta.length or (len(meta.sequence) or None)
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            text = fetch_structure_text(meta.cif_url, label=f"structure for {meta.accession}")
            res_ids, plddt = parse_plddt_from_cif(text, expected_residues=expected)
            return text, res_ids, plddt
        except TruncatedDownloadError as exc:
            last = exc
            if attempt < attempts - 1:
                time.sleep(RETRY_BACKOFF_S * (2 ** attempt))
    raise TruncatedDownloadError(
        f"The structure for {meta.accession} arrived incomplete on {attempts} separate "
        f"downloads ({last}). This usually means a slow or unstable connection."
    )


ALIGNMENT_COLOUR_SCHEMES = {
    "clustal": "Clustal", "clustal2": "Clustal", "zappo": "Zappo", "taylor": "Taylor",
    "hydrophobicity": "Hydrophobicity", "helix": "HelixPropensity", "strand": "StrandPropensity",
    "turn": "TurnPropensity", "buried": "BuriedIndex", "flower": "Flower", "blossom": "Blossom",
    "sunset": "Sunset", "ocean": "Ocean",
}
"""Maps the familiar MSA scheme names onto the ones pyMSAviz ships.

The notebook's Colab form offers these as a dropdown. A dropdown's option list has to be a
literal in a `#@param` comment, so it is a second copy of these keys and can drift from them;
`resolve_colour_scheme` exists so that drift is reported rather than silently rendered in the
wrong colours.
"""

#: The dropdown's options, in the order the form shows them. Keep the notebook's `#@param`
#: list in step with this: `ALIGNMENT_SCHEME_CHOICES` is the source, the comment is the copy.
ALIGNMENT_SCHEME_CHOICES = tuple(ALIGNMENT_COLOUR_SCHEMES)


def resolve_colour_scheme(name: str) -> str:
    """The pyMSAviz scheme for a familiar name, naming the fallback rather than taking it silently."""
    key = str(name).lower()
    if key not in ALIGNMENT_COLOUR_SCHEMES:
        print(f"Unknown alignment colour scheme {name!r}; using Clustal. "
              f"Valid choices: {', '.join(ALIGNMENT_SCHEME_CHOICES)}.")
        return "Clustal"
    return ALIGNMENT_COLOUR_SCHEMES[key]


def plot_alignment(result: AlignmentResult, query_label: str, role: str,
                   colour_scheme: str = "clustal2", wrap_length: int = 80) -> Figure:
    """Render one pairwise alignment, titled without overlapping the first sequence row.

    pyMSAviz sizes its figure to exactly fit the alignment rows, leaving no margin, so a plain
    `suptitle` lands on top of the first row. The figure is grown by a fixed strip and every
    axes rescaled into the space below it, which leaves the title somewhere to sit at any
    alignment length.
    """
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from pymsaviz import MsaViz

    alignment = MultipleSeqAlignment([
        SeqRecord(Seq(result.aligned_query), id=f"{query_label} (query)"),
        SeqRecord(Seq(result.aligned_target), id=f"{result.target} ({role})"),
    ])
    scheme = resolve_colour_scheme(colour_scheme)
    viz = MsaViz(alignment, color_scheme=scheme, wrap_length=wrap_length,
                 show_count=True, show_consensus=False)
    figure = viz.plotfig()

    width, height = figure.get_size_inches()
    strip = 0.5                                  # inches reserved for the title
    figure.set_size_inches(width, height + strip)
    scale = height / (height + strip)            # old [0, 1] maps into [0, scale]
    for axes in figure.axes:
        box = axes.get_position()
        axes.set_position([box.x0, box.y0 * scale, box.width, box.height * scale])
    figure.suptitle(f"{query_label} vs {role}: {result.summary_line()}",
                    fontsize=10, y=1.0, va="top")

    # pyMSAviz builds through pyplot, so its figure lands in pyplot's registry. Under
    # `%matplotlib inline` the backend flushes everything registered at the end of a cell, so
    # the caller's own `display()` plus that flush renders each alignment twice. Closing
    # deregisters it while leaving the Figure object fully renderable, which restores the
    # invariant every other plot function here already holds: nothing leaks into the registry.
    import matplotlib.pyplot as plt

    plt.close(figure)
    return figure


# --------------------------------------------------------------------------------------
# Structural superposition (TM-align) and the alignment it produces
#
# Mol* cannot compute an alignment. The MolViewSpec format has 32 node kinds and none of them
# aligns anything; its own superposition example hard-codes a nine-element matrix. What it does
# offer is `transform(rotation=..., translation=...)`, so the matrix is computed here and handed
# to it.
#
# TM-align rather than a sequence-driven Kabsch fit, because it is structure-based and does not
# depend on the sequences being alignable. Measured on the fixture set, that matters exactly
# where this notebook looks: at 25% identity it gave 2.94 A against Kabsch's 6.59 A. Mol* itself
# ships a TM-align implementation for the same reason, and its own UI recommends it "when
# sequence identity is low"; `tmtools` wraps the original Zhang and Skolnick C++ and is a small
# wheel rather than a 91 MB npm tree.
# --------------------------------------------------------------------------------------

#: TM-score is normalised to [0, 1] and, unlike RMSD, is interpretable on its own.
#: Zhang and Skolnick 2004: below 0.30 is the score of a random pair, above 0.50 implies the
#: same fold. RMSD alone cannot say this, because TM-align maximises TM-score rather than
#: minimising RMSD and will happily align two unrelated chains at a large RMSD.
TM_SCORE_SAME_FOLD = 0.50
TM_SCORE_RANDOM = 0.30


@dataclass(frozen=True)
class Superposition:
    """A structural alignment of one member onto the user's protein."""

    target: str
    rotation: tuple[float, ...]      # 9 values, column-major, as MolViewSpec expects
    translation: tuple[float, float, float]
    rmsd: float
    tm_score: float                  # normalised by the query's length
    aligned_length: int
    identity: float
    aligned_query: str               # TM-align's own structure-based alignment
    aligned_target: str
    query_res_ids: tuple[int, ...]   # label_seq_id per ungapped position, in alignment order
    target_res_ids: tuple[int, ...]

    @property
    def fold_verdict(self) -> str:
        if self.tm_score >= TM_SCORE_SAME_FOLD:
            return "same fold"
        if self.tm_score < TM_SCORE_RANDOM:
            return "no shared fold"
        return "uncertain"

    def summary_line(self) -> str:
        return (f"TM-score {self.tm_score:.3f} ({self.fold_verdict}), RMSD {self.rmsd:.2f} A "
                f"over {self.aligned_length} aligned residues, {self.identity:.1%} identity")


def ca_coordinates(cif_text: str) -> tuple["np.ndarray", list[int]]:
    """CA coordinates and their `label_seq_id`s, in file order."""
    columns: list[str] = []
    coords: list[tuple[float, float, float]] = []
    res_ids: list[int] = []
    in_loop = False
    for line in cif_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("_atom_site."):
            columns.append(stripped.split(".", 1)[1].split()[0].lower())
            in_loop = True
            continue
        if not in_loop:
            continue
        if stripped.startswith(("_", "loop_", "#")):
            break
        if not stripped:
            continue
        fields = stripped.split()
        if len(fields) != len(columns):
            continue
        row = dict(zip(columns, fields))
        if row.get("group_pdb") != "ATOM" or row.get("label_atom_id") != "CA":
            continue
        seq_id = row.get("label_seq_id", ".")
        if seq_id in (".", "?"):
            continue
        if row.get("pdbx_pdb_model_num") not in (None, "1"):
            continue
        coords.append((float(row["cartn_x"]), float(row["cartn_y"]), float(row["cartn_z"])))
        res_ids.append(int(seq_id))
    return np.asarray(coords, float), res_ids


def superpose(query_cif: str, target_cif: str, query_sequence: str, target_sequence: str,
              target_label: str) -> Superposition:
    """TM-align the target onto the query. Returns the transform for the **target**.

    tmtools returns the transform that moves chain 1 onto chain 2. The user's protein should
    stay put and the member should move onto it, so the inverse is returned: for a rotation
    matrix the inverse is its transpose, and the translation becomes ``-Rᵀ t``.
    """
    from tmtools import tm_align

    query_ca, query_res_ids = ca_coordinates(query_cif)
    target_ca, target_res_ids = ca_coordinates(target_cif)
    if query_ca.size == 0 or target_ca.size == 0:
        raise ClusterQualityError(f"No CA atoms parsed for {target_label}; cannot superpose.")

    result = tm_align(query_ca, target_ca,
                      query_sequence[:len(query_ca)], target_sequence[:len(target_ca)])
    rotation = result.u.T
    translation = -rotation @ result.t

    pairs = sum(1 for a, b in zip(result.seqxA, result.seqyA) if a != "-" and b != "-")
    identities = sum(1 for a, b in zip(result.seqxA, result.seqyA)
                     if a != "-" and b != "-" and a == b)
    return Superposition(
        target=target_label,
        # MolViewSpec wants the 9 values column-major (j * 3 + i), which is Fortran order.
        # Row-major here would silently render the structure in the wrong orientation.
        rotation=tuple(float(v) for v in rotation.flatten(order="F")),
        translation=tuple(float(v) for v in translation),
        rmsd=float(result.rmsd),
        tm_score=float(result.tm_norm_chain1),
        aligned_length=pairs,
        identity=identities / pairs if pairs else 0.0,
        aligned_query=result.seqxA,
        aligned_target=result.seqyA,
        query_res_ids=tuple(query_res_ids),
        target_res_ids=tuple(target_res_ids),
    )


def build_superposition_view(query_source: StructureSource, target_source: StructureSource,
                             fit: Superposition, query_colour: str | None = None,
                             target_colour: str | None = None):
    """Two structures in one viewer, the target rotated onto the query. Returns a `State`."""
    mvs = _require_molviewspec()
    query_colour = query_colour or SUPERPOSITION_QUERY_COLOUR
    target_colour = target_colour or SUPERPOSITION_TARGET_COLOUR
    builder = mvs.create_builder()
    (builder.download(url=query_source.url).parse(format=query_source.format).model_structure()
            .component(selector="polymer").representation(type="cartoon").color(color=query_colour))
    (builder.download(url=target_source.url).parse(format=target_source.format).model_structure()
            .transform(rotation=list(fit.rotation), translation=list(fit.translation))
            .component(selector="polymer").representation(type="cartoon").color(color=target_colour))
    return builder.get_state()


def parse_secondary_structure(cif_text: str) -> dict[int, str]:
    """`{label_seq_id: 'H' | 'E' | 'T'}` from the mmCIF `_struct_conf` records.

    AFDB files carry DSSP-style assignments, so no external tool is needed. Helix subtypes
    (alpha, 3-10, pi, polyproline) collapse to `H`, strands to `E`, turns and bends to `T`.
    """
    columns: list[str] = []
    rows: list[dict[str, str]] = []
    in_loop = False
    for line in cif_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("_struct_conf."):
            columns.append(stripped.split(".", 1)[1].split()[0])
            in_loop = True
            continue
        if not in_loop:
            continue
        if stripped.startswith(("_", "loop_", "#")) or not stripped:
            break
        fields = stripped.split()
        if len(fields) == len(columns):
            rows.append(dict(zip(columns, fields)))

    assignment: dict[int, str] = {}
    for row in rows:
        kind = row.get("conf_type_id", "")
        if kind.startswith("HELX"):
            code = "H"
        elif kind.startswith("STRN") or kind.startswith("BETA"):
            code = "E"
        else:
            code = "T"
        try:
            begin, end = int(row["beg_label_seq_id"]), int(row["end_label_seq_id"])
        except (KeyError, ValueError):
            continue
        for res_id in range(begin, end + 1):
            assignment[res_id] = code
    return assignment


#: DSSP-style track colours. Deliberately outside the pLDDT palette's orange-to-blue ramp, so a
#: structure row is never mistaken for a confidence row in the pLDDT-coloured view.
SS_TRACK_COLOURS = {"H": "#C994C7", "E": "#7BCCC4", "T": "#DDDDDD"}
SS_TRACK_LABELS = {"H": "helix", "E": "strand", "T": "turn or bend"}


def _project_onto_columns(aligned: str, res_ids: Sequence[int],
                          lookup: dict, default=None) -> list:
    """Spread a per-residue quantity across the alignment's columns.

    Keyed on `label_seq_id` rather than on position, because the pLDDT parser and the
    coordinate parser walk the file independently and a residue missing from one of them would
    otherwise shift every value after it onto the wrong column.
    """
    projected, index = [], 0
    for char in aligned:
        if char == "-":
            projected.append(default)
            continue
        res_id = res_ids[index] if index < len(res_ids) else None
        projected.append(lookup.get(res_id, default) if res_id is not None else default)
        index += 1
    return projected


def _reserve_title_strip(figure: Figure, title: str, strip: float = 0.5) -> None:
    """Grow a pyMSAviz figure upwards and rescale its axes, leaving room for a title.

    pyMSAviz sizes the figure to fit its rows exactly, so a plain `suptitle` lands on top of
    the first sequence.
    """
    width, height = figure.get_size_inches()
    figure.set_size_inches(width, height + strip)
    scale = height / (height + strip)
    for axes in figure.axes:
        box = axes.get_position()
        axes.set_position([box.x0, box.y0 * scale, box.width, box.height * scale])
    figure.suptitle(title, fontsize=10, y=1.0, va="top")


def plot_structural_alignment(fit: Superposition, query_label: str, role: str,
                              query_plddt: dict[int, float] | None = None,
                              target_plddt: dict[int, float] | None = None,
                              query_ss: dict[int, str] | None = None,
                              target_ss: dict[int, str] | None = None,
                              colour_by: str = "scheme",
                              colour_scheme: str = "clustal2",
                              wrap_length: int = 80) -> Figure:
    """Render TM-align's own alignment, with a secondary structure track under each sequence.

    This is the *structural* alignment, the one the superposition above it is made of, not a
    sequence alignment recomputed underneath it. Columns are paired because the residues
    superpose in space, so a column of two dissimilar residues is a real statement: the fold
    puts them in the same place.

    `colour_by='scheme'` uses the requested residue colour scheme; `colour_by='plddt'` paints
    each residue in its own model's confidence band, which turns the same picture into a map of
    where the two models are confident and where they are not.
    """
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from pymsaviz import MsaViz

    query_ss_row = _project_onto_columns(fit.aligned_query, fit.query_res_ids, query_ss or {}, "")
    target_ss_row = _project_onto_columns(fit.aligned_target, fit.target_res_ids, target_ss or {}, "")
    query_band = _project_onto_columns(fit.aligned_query, fit.query_res_ids, query_plddt or {})
    target_band = _project_onto_columns(fit.aligned_target, fit.target_res_ids, target_plddt or {})

    def as_track(codes: list) -> str:
        return "".join(code if code else "-" for code in codes)

    # Two spaces and three spaces: the labels must differ or the two tracks collide, but the
    # difference must not be visible, because both mean the same thing.
    rows = [
        SeqRecord(Seq(fit.aligned_query), id=f"{query_label} (query)"),
        SeqRecord(Seq(as_track(query_ss_row)), id="  SSE"),
        SeqRecord(Seq(fit.aligned_target), id=f"{fit.target} ({role})"),
        SeqRecord(Seq(as_track(target_ss_row)), id="   SSE"),
    ]
    scheme = resolve_colour_scheme(colour_scheme)
    viz = MsaViz(MultipleSeqAlignment(rows), color_scheme=scheme, wrap_length=wrap_length,
                 show_count=False, show_consensus=False)

    bands = {0: query_band, 2: target_band}
    use_plddt = str(colour_by).lower() == "plddt"

    def colour(row: int, col: int, char: str, _msa) -> str | None:
        if row in (1, 3):
            return SS_TRACK_COLOURS.get(char)          # the structure rows, in both views
        if not use_plddt or char == "-":
            return None                                # None falls through to the scheme
        value = bands[row][col] if col < len(bands[row]) else None
        return plddt_colour(float(value)) if value is not None else "#FFFFFF"

    viz.set_custom_color_func(colour)
    figure = viz.plotfig()

    view = "confidence bands" if use_plddt else f"{scheme} colours"
    _reserve_title_strip(
        figure, f"{query_label} vs {role} ({fit.target}), {view}: {fit.summary_line()}")

    import matplotlib.pyplot as plt

    plt.close(figure)                                   # see plot_alignment for why
    return figure


def superposition_legend_html() -> str:
    """A swatch row for the structure track, escaped for the same reason the pLDDT one is."""
    swatches = "".join(_swatch(colour, f"{SS_TRACK_LABELS[code]} ({code})")
                       for code, colour in SS_TRACK_COLOURS.items())
    return f'<div style="font-size:0.9em; margin:6px 0;">SSE track: {swatches}</div>'


#: The two chains in a superposition view. Grey and green rather than anything from the pLDDT
#: ramp, because here colour means "which protein", not "how confident".
SUPERPOSITION_QUERY_COLOUR = "#5A5A5A"
SUPERPOSITION_TARGET_COLOUR = "#1B9E77"


def _swatch(colour: str, text: str) -> str:
    import html as _html

    return (f'<span style="display:inline-block; margin-right:14px;">'
            f'<span style="display:inline-block; width:12px; height:12px; background:{colour}; '
            f'border:1px solid #999; vertical-align:middle; margin-right:5px;"></span>'
            f'{_html.escape(text)}</span>')


def superposition_key_html(query_label: str, target_label: str) -> str:
    """Which chain is which, in the colours `build_superposition_view` uses."""
    return ('<div style="font-size:0.9em; margin:6px 0;">'
            + _swatch(SUPERPOSITION_QUERY_COLOUR, f"{query_label} (your protein, held fixed)")
            + _swatch(SUPERPOSITION_TARGET_COLOUR, f"{target_label} (rotated onto it)")
            + "</div>")


# --------------------------------------------------------------------------------------
# Notebook orchestration
#
# One function per notebook cell. The notebook is a reading surface, not a place to keep
# formatting loops and try/except ladders, so everything below is the cell bodies moved here
# verbatim: same calls, same order, same wording. Each returns whatever the next cell needs, so
# the notebook's data flow stays visible as named variables rather than hidden in a session
# object.
#
# `report_*` prints, `show_*` renders. The split matters: the reports are mandatory output and
# must never sit behind a figure's try/except (REQ-011a).
# --------------------------------------------------------------------------------------

#: Installed by the notebook's setup cell when absent. `tmtools` is deliberately **not** here:
#: it is installed by the superposition section itself, so that section stays deletable whole.
OPTIONAL_PACKAGES = (("molviewspec", "molviewspec"), ("biopython", "Bio"), ("pymsaviz", "pymsaviz"))


def ensure_optional_packages(packages: Sequence[tuple[str, str]] = OPTIONAL_PACKAGES) -> None:
    """Install the optional packages Colab does not ship. A no-op when they are present."""
    import importlib.util
    import subprocess
    import sys as _sys

    missing = [dist for dist, module in packages if importlib.util.find_spec(module) is None]
    if missing:
        subprocess.run([_sys.executable, "-m", "pip", "install", "-q", *missing], check=True)
    print("Optional 3D and alignment packages:", "installed" if missing else "already present")


def resolve_accession(raw: str) -> str:
    """Normalise the accession and print the interpretation when one was needed (REQ-001)."""
    accession, note = normalise_accession(raw)
    if note:
        print(note)
    return accession


def report_query(accession: str) -> QueryMeta:
    """The identity card for the protein being asked about. If it names the wrong one, stop."""
    import datetime

    meta = fetch_prediction(accession)
    print(f"{meta.accession}  {meta.description}")
    print(f"  organism      {meta.organism}")
    print(f"  length        {meta.length} aa")
    print(f"  global pLDDT  {meta.plddt:.2f}")
    print(f"  AFDB version  {meta.latest_version}   run {datetime.datetime.now():%Y-%m-%d %H:%M}")
    if is_isoform(accession):
        print(f"\n  Note: {accession} is an isoform. Its own sequence and structure are used for the\n"
              f"  alignments and 3D views, while the cluster is the one AFDB returns for it.")
    return meta


@dataclass(frozen=True)
class ClusterPair:
    """Both clusterings for one accession, with the size gate already applied."""

    accession: str
    sequence: ClusterMembers
    structure: ClusterMembers | None
    show_structure_section: bool
    plddt_axis: tuple[float, float]


def load_clusters(accession: str) -> ClusterPair:
    """Fetch both clusterings, report their sizes, and apply the refusal gates.

    Degrades **per clustering** (REQ-014a): an unusable structure cluster skips one section,
    while an unusable sequence family stops the run, because everything downstream rests on it.
    """
    sequence = fetch_cluster(accession, "sequence")
    try:
        structure = fetch_cluster(accession, "structure")
    except ClusterQualityError as exc:
        structure = None
        print(f"Structure clustering unavailable: {exc}")

    print(f"Sequence family   {len(sequence):>7,} members  (clusterTotal {sequence.cluster_total:,})")
    if structure is not None:
        print(f"Structure cluster {len(structure):>7,} members  (clusterTotal {structure.cluster_total:,})")
    if len(sequence) != sequence.cluster_total:
        print("  Note: clusterTotal disagrees with the returned array length.")

    if len(sequence) < MIN_CLUSTER_MEMBERS:
        raise ClusterTooSmallError(
            f"The sequence family for {accession} has {len(sequence)} member(s). "
            f"At least {MIN_CLUSTER_MEMBERS} are needed to show a distribution and its extremes. "
            "Try a protein with a larger family.")

    show_structure_section = structure is not None and len(structure) >= MIN_CLUSTER_MEMBERS
    if not show_structure_section:
        n = 0 if structure is None else len(structure)
        print(f"The structure-cluster section will be skipped: it has {n} member(s). "
              "The sequence-family analysis below is unaffected.")

    return ClusterPair(
        accession=accession,
        sequence=sequence,
        structure=structure,
        show_structure_section=show_structure_section,
        plddt_axis=shared_plddt_axis_limits(sequence, structure if show_structure_section else None),
    )


def report_family_summary(clusters: ClusterPair) -> FamilySummary:
    """The cluster summary table (REQ-003)."""
    family = summarise_family(clusters.sequence)
    print(f"Cluster summary: {family.n:,} members\n")
    print(f"  Mean                {family.mean:8.2f}")
    print(f"  Median              {family.median:8.2f}")
    print(f"  Standard deviation  {family.sd:8.2f}" if family.sd is not None
          else "  Standard deviation       n/a  (needs more than one member)")
    print(f"  p25 to p75          {family.iqr_width:8.2f}   pLDDT points wide")
    print("  Confidence bands    " + "   ".join(f"{k} {v:5.1%}" for k, v in family.bands.items()))
    if family.low_n:
        print(f"\n  Only {family.n} members: with a cluster this small the interquartile range above,")
        print("  and the percentile in the next section, are both unstable.")
    return family


def report_query_position(clusters: ClusterPair) -> Position:
    """Where the user's protein sits in its own family (REQ-004)."""
    accession = clusters.accession
    position = locate_query(clusters.sequence, accession)
    if position.present:
        value = float(clusters.sequence.plddt[position.index])
        print(f"{accession}  mean pLDDT {value:.2f}\n")
        print(f"  Descriptor          {position.category}")
        print(f"  Percentile          {position.percentile:5.1f}   rank within this cluster")
        print(f"  Cluster p25 to p75  {position.iqr_width:5.2f}   pLDDT points wide")
        if position.iqr_width < 2.0:
            print("\n  That width is narrow, so the percentile carries little information here.")
    else:
        print(f"{accession} is not itself a member of the cluster that came back, so it has no rank in it.")
        print("  This is a normal outcome; the structure-cluster section below explains why.")
    return position


def show_family_distribution(clusters: ClusterPair, family: FamilySummary, position: Position) -> None:
    """The sequence-family raincloud, on the shared axis (REQ-V1, REQ-V1a)."""
    from IPython.display import display

    display(plot_family_distribution(clusters.sequence, family, position, clusters.accession,
                                     xlim=clusters.plddt_axis))


def show_structure_cluster(clusters: ClusterPair) -> None:
    """The structure-cluster raincloud, or a named skip (REQ-V2, REQ-014a)."""
    from IPython.display import display

    if not clusters.show_structure_section:
        print("Skipped: no usable structure cluster for this protein.")
        return
    structure = clusters.structure
    representative = structure.table.iloc[0]["accession"]
    in_structure_cluster = structure.index_of(clusters.accession) is not None
    print(f"Structure cluster of {len(structure):,} families, representative {representative}")
    if not in_structure_cluster:
        print(f"  {clusters.accession} is not itself in this cluster, which is expected: what is clustered\n"
              f"  structurally is its sequence family's representative, not the protein you asked about.")
    display(plot_structure_cluster_distribution(structure, xlim=clusters.plddt_axis))


def show_length_vs_plddt(clusters: ClusterPair) -> None:
    """The length-versus-pLDDT joint density (REQ-V4)."""
    from IPython.display import display

    try:
        display(plot_length_vs_plddt(clusters.sequence, "sequence family"))
    except Exception as exc:
        print(f"Length figure unavailable: {exc}")


@dataclass
class ModelSelection:
    """The 2 or 3 models the rest of the notebook renders, with their downloads cached.

    The cache is the reason this is an object rather than a tuple. The 3D views, the pLDDT
    profiles and the superposition all need the same mmCIF, and each file is large, served
    without a length header and unresumable. Downloading it once per model rather than once per
    section is the difference between a usable notebook and a slow one on a poor connection.
    """

    accession: str
    query: QueryMeta
    source: str                                  # "sequence" or "structure"
    targets: list[ViewTarget]
    details: dict[str, QueryMeta]
    _structures: dict[str, tuple[str, list[int], list[float]]] = field(default_factory=dict, repr=False)

    def structure(self, role: str) -> tuple[str, list[int], list[float]]:
        """The mmCIF text, residue ids and pLDDT for one role, downloaded at most once."""
        if role not in self._structures:
            self._structures[role] = fetch_verified_structure(self.details[role])
        return self._structures[role]

    @property
    def partners(self) -> list[ViewTarget]:
        """The targets that are not the user's own protein."""
        return [t for t in self.targets if t.role != "query"]


def select_models(clusters: ClusterPair, query: QueryMeta, source: str = "sequence") -> ModelSelection:
    """Pick the best and worst models and report them (REQ-009, REQ-010, REQ-010b).

    The collapse rule runs **once**, here, so the 3D views, the alignments and the superposition
    can never disagree about which models they are showing.
    """
    if source == "structure" and not clusters.show_structure_section:
        print("No usable structure cluster; falling back to the sequence family for the models below.")
        source = "sequence"
    model_cluster = clusters.sequence if source == "sequence" else clusters.structure
    print(f"Models drawn from the {'sequence family' if source == 'sequence' else 'structure cluster'}"
          f" ({len(model_cluster):,} members)\n")

    targets = identify_view_targets(model_cluster, query)
    details = {}
    for target in targets:
        details[target.role] = query if target.role == "query" else fetch_prediction(target.accession)
        print(f"{target.role:6s} {target.afdb_id:20s} pLDDT {target.plddt:6.2f}  {target.description[:60]}")
    if len(targets) == 2:
        print("\nYour protein is itself one of the extremes, so two models are shown rather than three.")
    if source == "structure":
        print("\nNote: these are family representatives, each the best model of its own sequence family,")
        print("not arbitrary individual proteins.")

    return ModelSelection(accession=clusters.accession, query=query, source=source,
                          targets=targets, details=details)


def show_model_panels(models: ModelSelection, show_pae: bool = True) -> None:
    """One Mol\\* view per model, pLDDT-coloured, with its PAE matrix beside it (REQ-010, REQ-V10)."""
    from IPython.display import HTML, display

    if not molviewspec_available():
        print(MOLVIEWSPEC_MISSING_MESSAGE)
        return
    display(HTML(plddt_band_legend_html()))
    for target in models.targets:
        meta = models.details[target.role]
        try:
            _, res_ids, plddt = models.structure(target.role)
            state = build_plddt_view(resolve_structure_source(meta.cif_url), res_ids, plddt)
            pae_figure = None
            if show_pae:
                try:
                    pae_figure = plot_pae_matrix(fetch_pae(meta), target.afdb_id)
                except ClusterQualityError as pae_exc:
                    print(f"  PAE unavailable for {target.afdb_id}: {pae_exc}")
            show_model_panel(
                state, pae_figure,
                f"{target.role}: {target.afdb_id}  mean pLDDT {target.plddt:.1f}  "
                f"({target.description[:55]})")
        except ClusterQualityError as exc:
            print(f"3D view unavailable for {target.afdb_id}:\n  {exc}")
        except Exception as exc:
            print(f"3D view unavailable for {target.afdb_id}: {type(exc).__name__}: {exc}")


def show_plddt_profiles(models: ModelSelection) -> None:
    """Per-residue pLDDT for each rendered model (REQ-010a, REQ-V8)."""
    from IPython.display import display

    for target in models.targets:
        try:
            _, res_ids, plddt = models.structure(target.role)
            display(plot_plddt_profile(res_ids, plddt, f"{target.role}: {target.afdb_id}"))
        except ClusterQualityError as exc:
            print(f"Profile unavailable for {target.afdb_id}:\n  {exc}")


def report_alignments(models: ModelSelection) -> dict[str, AlignmentResult]:
    """Local alignment against each extreme. Mandatory output, so this never guards a figure."""
    alignments = {}
    for target in models.partners:
        result = align_pair(models.query.sequence, models.details[target.role].sequence, target.afdb_id)
        alignments[target.role] = result
        print(f"{models.accession} vs {target.role} ({target.afdb_id})")
        print(f"  {result.summary_line()}")
    return alignments


def show_alignments(models: ModelSelection, alignments: dict[str, AlignmentResult],
                    colour_scheme: str = "clustal2") -> None:
    """The pyMSAviz panels for REQ-V7. Optional: the numbers above are already printed."""
    from IPython.display import display

    try:
        for role, result in alignments.items():
            display(plot_alignment(result, models.accession, role, colour_scheme))
    except Exception as exc:
        print(f"Alignment figures unavailable ({exc}). The numbers above are unaffected.")


def report_superpositions(models: ModelSelection) -> dict[str, Superposition]:
    """TM-align the query onto each extreme and print the result (REQ-017, REQ-017a).

    Installs `tmtools` itself rather than relying on the shared setup cell, so the whole
    superposition section remains deletable as one block.
    """
    import importlib.util
    import subprocess
    import sys as _sys

    if importlib.util.find_spec("tmtools") is None:
        subprocess.run([_sys.executable, "-m", "pip", "install", "-q", "tmtools"], check=True)

    superpositions = {}
    query_text, _, _ = models.structure("query")
    for target in models.partners:
        target_text, _, _ = models.structure(target.role)
        fit = superpose(query_text, target_text, models.query.sequence,
                        models.details[target.role].sequence, target.afdb_id)
        superpositions[target.role] = fit
        print(f"{models.accession} vs {target.role} ({target.afdb_id})")
        print(f"  {fit.summary_line()}")
    return superpositions


def show_superposition_blocks(models: ModelSelection, superpositions: dict[str, Superposition],
                              colour_scheme: str = "clustal2") -> None:
    """One block per model: the superposition, then the same alignment twice (REQ-017b, REQ-V11).

    Grouped by model rather than by output type, so each superposition sits with the alignment
    that explains it.
    """
    from IPython.display import HTML, display

    display(HTML(plddt_band_legend_html() + superposition_legend_html()))
    have_viewer = molviewspec_available()
    if not have_viewer:
        print(MOLVIEWSPEC_MISSING_MESSAGE)

    for role, fit in superpositions.items():
        query_text, query_ids, query_plddt = models.structure("query")
        target_text, target_ids, target_plddt = models.structure(role)

        if have_viewer:
            try:
                state = build_superposition_view(
                    resolve_structure_source(models.query.cif_url),
                    resolve_structure_source(models.details[role].cif_url), fit)
                display(HTML(superposition_key_html(models.accession, fit.target)))
                show_mol_view(state, f"{models.accession} superposed with the {role} model ({fit.target}): "
                                     f"TM-score {fit.tm_score:.3f}, RMSD {fit.rmsd:.2f} A")
            except Exception as exc:
                print(f"Superposition view unavailable for {fit.target}: {type(exc).__name__}: {exc}")

        try:
            shared = dict(
                query_plddt=dict(zip(query_ids, query_plddt)),
                target_plddt=dict(zip(target_ids, target_plddt)),
                query_ss=parse_secondary_structure(query_text),
                target_ss=parse_secondary_structure(target_text))
            for mode in ("scheme", "plddt"):
                display(plot_structural_alignment(
                    fit, models.accession, role, colour_by=mode,
                    colour_scheme=colour_scheme, **shared))
        except Exception as exc:
            print(f"Alignment panels unavailable for {fit.target} ({exc}). "
                  "The numbers above are unaffected.")


def report_msa_status(user_msa_path=None) -> None:
    """REQ-012. The endpoint returns 403, so this path is declared untested rather than guessed at."""
    print("MSA coverage: the AFDB endpoint is under maintenance (HTTP 403) and this path is untested.")
    if user_msa_path:
        print(f"A user-supplied MSA was given at {user_msa_path}; validation against the AFDB sequence "
              "is required before it is used.")


def report_summary(models: ModelSelection, family: FamilySummary, position: Position,
                   alignments: dict[str, AlignmentResult]) -> None:
    """The closing summary and the caveats that bound it (REQ-016, REQ-015)."""
    print(f"{models.accession}  {models.query.description}")
    print(f"  Family        {family.n:,} members, mean pLDDT {family.mean:.1f}, "
          f"p25 to p75 {family.iqr_width:.2f} points wide")
    if position.present:
        print(f"  Your protein  {models.query.plddt:.1f}, {position.category} for this family")
    for role, result in alignments.items():
        print(f"  vs {role:5s}     {result.summary_line()}")
    print("\n  Bounding caveats: the percentile is uninformative in the body of the distribution;")
    print("  a chain-average pLDDT partly measures disorder content; and a cluster member may be a")
    print("  distant homologue with a different biological function.")
