"""Evaluate LLM-written Markdown project specs with a stronger judge model.

The module is intentionally dependency-light so it can be used from a notebook,
plain Python script, or tests without installing an SDK. The default judge uses
OpenAI's Chat Completions API through the Python standard library.
"""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass
from html import escape
from pathlib import Path
from typing import Any, Callable, Optional


CORE_SCORE_KEYS = (
    "correctness",
    "completeness",
    "specificity",
    "coherence",
    "actionability",
)
REFERENCE_SCORE_KEYS = ("reference_alignment",)
ALL_SCORE_KEYS = CORE_SCORE_KEYS + REFERENCE_SCORE_KEYS


@dataclass(frozen=True)
class Candidate:
    project: str
    model: str
    path: Path
    markdown: str


@dataclass(frozen=True)
class ReferenceSpec:
    project: str
    path: Path
    markdown: str


@dataclass
class EvaluationRecord:
    project: str
    model: str
    candidate_path: str
    reference_path: str | None
    used_reference: bool
    status: str
    scores: dict[str, int | None]
    overall_score: float | None
    justification: str
    strengths: list[str]
    weaknesses: list[str]
    suggested_improvements: list[str]
    missing_from_reference: list[str]
    extra_or_conflicting_claims: list[str]
    error: str | None = None
    raw_response: dict[str, Any] | None = None


JudgeCallable = Callable[[Candidate, Optional[ReferenceSpec]], dict[str, Any]]


class EvaluationError(RuntimeError):
    """Raised when an evaluation cannot be completed or normalized."""


def parse_candidate_filename(path: Path, separator: str = "--") -> tuple[str, str]:
    """Parse `<project><separator><model>.md` into `(project, model)`.

    The separator is explicit so both project and model names may safely contain
    hyphens or underscores.
    """
    if path.suffix.lower() != ".md":
        raise ValueError(f"Expected a Markdown file: {path.name}")

    stem = path.stem
    if separator not in stem:
        raise ValueError(
            f"{path.name!r} does not match '<project>{separator}<model>.md'"
        )

    project, model = stem.rsplit(separator, 1)
    if not project or not model:
        raise ValueError(
            f"{path.name!r} must include both project and model names"
        )
    return project, model


def discover_candidates(
    input_dir: str | Path,
    separator: str = "--",
    strict: bool = True,
) -> dict[str, list[Candidate]]:
    """Load candidate Markdown files grouped by project."""
    directory = Path(input_dir)
    if not directory.exists():
        raise FileNotFoundError(f"Input directory does not exist: {directory}")
    if not directory.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {directory}")

    grouped: dict[str, list[Candidate]] = {}
    invalid: list[str] = []

    for path in sorted(directory.glob("*.md")):
        try:
            project, model = parse_candidate_filename(path, separator)
        except ValueError as exc:
            invalid.append(str(exc))
            continue
        grouped.setdefault(project, []).append(
            Candidate(
                project=project,
                model=model,
                path=path,
                markdown=path.read_text(encoding="utf-8"),
            )
        )

    if strict and invalid:
        details = "\n".join(f"- {item}" for item in invalid)
        raise ValueError(f"Invalid candidate filename(s):\n{details}")

    return grouped


def load_reference(
    reference_dir: str | Path | None,
    project: str,
) -> ReferenceSpec | None:
    """Load `<reference_dir>/<project>.md` if a reference directory is set."""
    if reference_dir is None:
        return None

    path = Path(reference_dir) / f"{project}.md"
    if not path.exists():
        return None
    return ReferenceSpec(
        project=project,
        path=path,
        markdown=path.read_text(encoding="utf-8"),
    )


def judge_response_schema() -> dict[str, Any]:
    """JSON schema requested from the judge model."""
    return {
        "type": "object",
        "additionalProperties": False,
        "required": [
            "scores",
            "justification",
            "strengths",
            "weaknesses",
            "suggested_improvements",
            "missing_from_reference",
            "extra_or_conflicting_claims",
        ],
        "properties": {
            "scores": {
                "type": "object",
                "additionalProperties": False,
                "required": list(ALL_SCORE_KEYS),
                "properties": {
                    key: {"type": ["integer", "null"]}
                    for key in ALL_SCORE_KEYS
                },
            },
            "justification": {"type": "string"},
            "strengths": {
                "type": "array",
                "items": {"type": "string"},
            },
            "weaknesses": {
                "type": "array",
                "items": {"type": "string"},
            },
            "suggested_improvements": {
                "type": "array",
                "items": {"type": "string"},
            },
            "missing_from_reference": {
                "type": "array",
                "items": {"type": "string"},
            },
            "extra_or_conflicting_claims": {
                "type": "array",
                "items": {"type": "string"},
            },
        },
    }


def build_judge_messages(
    candidate: Candidate,
    reference: ReferenceSpec | None = None,
) -> list[dict[str, str]]:
    """Build the judging prompt."""
    reference_block = (
        "No reference specification was provided. Run intrinsic evaluation only."
        if reference is None
        else (
            "Reference specification:\n"
            f"<reference_markdown>\n{reference.markdown}\n</reference_markdown>"
        )
    )

    user_prompt = f"""
Evaluate the candidate Markdown project specification below.

Judge only the submitted content. Do not infer unstated intent, missing
requirements, or hidden context. Distinguish factual/technical problems from
missing detail. Use the reference specification only when it is provided.

Rubric:
- correctness: technical accuracy and absence of false claims
- completeness: coverage of required product/spec requirements
- specificity: implementation-level clarity and lack of vague prose
- coherence: organization, consistency, and readability
- actionability: how ready the spec is for an engineer to implement
- reference_alignment: agreement with the reference, or null if no reference

Scores are integers from 1 to 5, where 5 is excellent.

Project: {candidate.project}
Model under evaluation: {candidate.model}

{reference_block}

Candidate specification:
<candidate_markdown>
{candidate.markdown}
</candidate_markdown>
""".strip()

    return [
        {
            "role": "system",
            "content": (
                "You are a rigorous senior software/product specification "
                "reviewer. Return only JSON that matches the provided schema."
            ),
        },
        {"role": "user", "content": user_prompt},
    ]


def call_openai_chat_json(
    messages: list[dict[str, str]],
    model: str,
    api_key: str | None = None,
    timeout_seconds: int = 120,
) -> dict[str, Any]:
    """Call OpenAI Chat Completions and parse a structured JSON response."""
    key = api_key or os.environ.get("OPENAI_API_KEY")
    if not key:
        raise EvaluationError("OPENAI_API_KEY is not set")

    body = {
        "model": model,
        "messages": messages,
        "response_format": {
            "type": "json_schema",
            "json_schema": {
                "name": "spec_evaluation",
                "strict": True,
                "schema": judge_response_schema(),
            },
        },
    }

    request = urllib.request.Request(
        "https://api.openai.com/v1/chat/completions",
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {key}",
            "Content-Type": "application/json",
        },
        method="POST",
    )

    try:
        with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise EvaluationError(f"OpenAI API error {exc.code}: {detail}") from exc
    except urllib.error.URLError as exc:
        raise EvaluationError(f"OpenAI API request failed: {exc}") from exc

    try:
        content = payload["choices"][0]["message"]["content"]
        return json.loads(content)
    except (KeyError, IndexError, TypeError, json.JSONDecodeError) as exc:
        raise EvaluationError(f"Could not parse judge response: {payload}") from exc


def make_openai_judge(model: str, api_key: str | None = None) -> JudgeCallable:
    """Create a judge callable backed by OpenAI."""

    def judge(candidate: Candidate, reference: ReferenceSpec | None) -> dict[str, Any]:
        messages = build_judge_messages(candidate, reference)
        try:
            return call_openai_chat_json(messages, model=model, api_key=api_key)
        except EvaluationError:
            repair_messages = messages + [
                {
                    "role": "user",
                    "content": (
                        "Retry once. Return valid JSON only, exactly matching "
                        "the requested schema."
                    ),
                }
            ]
            return call_openai_chat_json(
                repair_messages,
                model=model,
                api_key=api_key,
            )

    return judge


def normalize_judge_response(
    response: dict[str, Any],
    *,
    used_reference: bool,
) -> tuple[dict[str, int | None], float, dict[str, Any]]:
    """Validate model output and compute a deterministic weighted score."""
    if not isinstance(response, dict):
        raise EvaluationError("Judge response must be a JSON object")

    raw_scores = response.get("scores")
    if not isinstance(raw_scores, dict):
        raise EvaluationError("Judge response is missing scores object")

    scores: dict[str, int | None] = {}
    for key in ALL_SCORE_KEYS:
        value = raw_scores.get(key)
        if key in REFERENCE_SCORE_KEYS and not used_reference:
            scores[key] = None
            continue
        if value is None:
            raise EvaluationError(f"Score {key!r} is required")
        if not isinstance(value, int) or not 1 <= value <= 5:
            raise EvaluationError(f"Score {key!r} must be an integer from 1 to 5")
        scores[key] = value

    weights = {key: 1.0 for key in CORE_SCORE_KEYS}
    if used_reference:
        weights["reference_alignment"] = 1.5

    numerator = sum(float(scores[key] or 0) * weight for key, weight in weights.items())
    denominator = sum(weights.values())
    overall_score = round(numerator / denominator, 2)

    normalized = {
        "justification": _require_string(response, "justification"),
        "strengths": _require_string_list(response, "strengths"),
        "weaknesses": _require_string_list(response, "weaknesses"),
        "suggested_improvements": _require_string_list(
            response,
            "suggested_improvements",
        ),
        "missing_from_reference": _require_string_list(
            response,
            "missing_from_reference",
        ),
        "extra_or_conflicting_claims": _require_string_list(
            response,
            "extra_or_conflicting_claims",
        ),
    }
    return scores, overall_score, normalized


def evaluate_candidate(
    candidate: Candidate,
    reference: ReferenceSpec | None,
    judge: JudgeCallable,
) -> EvaluationRecord:
    """Evaluate a single candidate and return a success or error record."""
    try:
        raw_response = judge(candidate, reference)
        scores, overall_score, normalized = normalize_judge_response(
            raw_response,
            used_reference=reference is not None,
        )
        return EvaluationRecord(
            project=candidate.project,
            model=candidate.model,
            candidate_path=str(candidate.path),
            reference_path=str(reference.path) if reference else None,
            used_reference=reference is not None,
            status="success",
            scores=scores,
            overall_score=overall_score,
            raw_response=raw_response,
            error=None,
            **normalized,
        )
    except Exception as exc:  # Keep batch evaluation moving.
        return EvaluationRecord(
            project=candidate.project,
            model=candidate.model,
            candidate_path=str(candidate.path),
            reference_path=str(reference.path) if reference else None,
            used_reference=reference is not None,
            status="error",
            scores={key: None for key in ALL_SCORE_KEYS},
            overall_score=None,
            justification="",
            strengths=[],
            weaknesses=[],
            suggested_improvements=[],
            missing_from_reference=[],
            extra_or_conflicting_claims=[],
            error=str(exc),
            raw_response=None,
        )


def evaluate_directory(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    reference_dir: str | Path | None = None,
    model: str = "gpt-5",
    separator: str = "--",
    judge: JudgeCallable | None = None,
) -> dict[str, Any]:
    """Evaluate all candidate specs and save JSON plus Markdown outputs."""
    grouped = discover_candidates(input_dir, separator=separator, strict=True)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    active_judge = judge or make_openai_judge(model)
    records: list[EvaluationRecord] = []

    for project, candidates in grouped.items():
        reference = load_reference(reference_dir, project)
        for candidate in candidates:
            records.append(evaluate_candidate(candidate, reference, active_judge))

    result = {
        "metadata": {
            "created_at_unix": int(time.time()),
            "input_dir": str(Path(input_dir)),
            "reference_dir": str(Path(reference_dir)) if reference_dir else None,
            "model": model,
            "separator": separator,
        },
        "evaluations": [asdict(record) for record in records],
        "rankings": build_rankings(records),
    }
    result["plots"] = write_project_plots(result, output_path)

    (output_path / "evaluations.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (output_path / "summary.md").write_text(
        render_markdown_summary(result),
        encoding="utf-8",
    )
    return result


def write_project_plots(
    result: dict[str, Any],
    output_dir: str | Path,
) -> dict[str, dict[str, str]]:
    """Write dependency-free SVG comparison plots and return their paths."""
    output_path = Path(output_dir)
    plots_dir = output_path / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    plots: dict[str, dict[str, str]] = {}
    for project in result["rankings"]:
        records = [
            item
            for item in result["evaluations"]
            if item["project"] == project and item["status"] == "success"
        ]
        if not records:
            continue

        categories = score_categories_for_records(records)
        if not categories:
            continue

        bar_path = plots_dir / f"{safe_filename(project)}_category_bars.svg"
        heatmap_path = plots_dir / f"{safe_filename(project)}_score_heatmap.svg"
        bar_path.write_text(
            render_grouped_bar_svg(project, records, categories),
            encoding="utf-8",
        )
        heatmap_path.write_text(
            render_score_heatmap_svg(project, records, categories),
            encoding="utf-8",
        )
        plots[project] = {
            "category_bars": str(bar_path.relative_to(output_path)),
            "score_heatmap": str(heatmap_path.relative_to(output_path)),
        }

    return plots


def score_categories_for_records(records: list[dict[str, Any]]) -> list[str]:
    """Return score categories worth plotting for a project's records."""
    categories = list(CORE_SCORE_KEYS)
    has_reference = any(
        record.get("scores", {}).get("reference_alignment") is not None
        for record in records
    )
    if has_reference:
        categories.append("reference_alignment")
    return categories


def render_grouped_bar_svg(
    project: str,
    records: list[dict[str, Any]],
    categories: list[str],
) -> str:
    """Render a grouped bar chart comparing models by score category."""
    records = sorted(
        records,
        key=lambda item: item["overall_score"] if item["overall_score"] is not None else -1,
        reverse=True,
    )
    margin_left = 86
    margin_top = 52
    plot_width = max(620, len(categories) * max(72, len(records) * 28))
    plot_height = 260
    legend_height = max(42, ((len(records) + 2) // 3) * 24)
    width = margin_left + plot_width + 28
    height = margin_top + plot_height + legend_height + 74
    baseline = margin_top + plot_height
    group_width = plot_width / len(categories)
    bar_gap = 4
    bar_width = max(10, min(24, (group_width - 22) / max(1, len(records)) - bar_gap))
    colors = chart_palette(len(records))

    parts = [svg_header(width, height)]
    parts.append(
        f'<text x="{margin_left}" y="28" class="title">Category Scores: {escape(project)}</text>'
    )
    parts.extend(render_y_axis(margin_left, margin_top, plot_height, plot_width))

    for category_index, category in enumerate(categories):
        group_start = margin_left + category_index * group_width
        bars_width = len(records) * (bar_width + bar_gap) - bar_gap
        bars_start = group_start + (group_width - bars_width) / 2
        for record_index, record in enumerate(records):
            value = record["scores"].get(category)
            if value is None:
                continue
            bar_height = plot_height * (float(value) / 5.0)
            x = bars_start + record_index * (bar_width + bar_gap)
            y = baseline - bar_height
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" '
                f'height="{bar_height:.1f}" fill="{colors[record_index]}" rx="2">'
                f'<title>{escape(record["model"])} {escape(category)}: {value}/5</title>'
                "</rect>"
            )
        label_x = group_start + group_width / 2
        parts.append(
            f'<text x="{label_x:.1f}" y="{baseline + 22}" class="axis-label" '
            f'text-anchor="middle">{escape(short_category_label(category))}</text>'
        )

    parts.extend(render_legend(records, colors, margin_left, baseline + 50, width))
    parts.append("</svg>")
    return "\n".join(parts)


def render_score_heatmap_svg(
    project: str,
    records: list[dict[str, Any]],
    categories: list[str],
) -> str:
    """Render a heatmap of score categories by model."""
    records = sorted(
        records,
        key=lambda item: item["overall_score"] if item["overall_score"] is not None else -1,
        reverse=True,
    )
    cell_width = 112
    cell_height = 34
    label_width = 150
    top = 58
    left = label_width + 20
    width = left + len(categories) * cell_width + 24
    height = top + len(records) * cell_height + 62

    parts = [svg_header(width, height)]
    parts.append(
        f'<text x="20" y="28" class="title">Score Heatmap: {escape(project)}</text>'
    )
    for index, category in enumerate(categories):
        x = left + index * cell_width + cell_width / 2
        parts.append(
            f'<text x="{x:.1f}" y="48" class="axis-label" text-anchor="middle">'
            f'{escape(short_category_label(category))}</text>'
        )

    for row_index, record in enumerate(records):
        y = top + row_index * cell_height
        parts.append(
            f'<text x="20" y="{y + 22}" class="row-label">{escape(record["model"])}</text>'
        )
        for column_index, category in enumerate(categories):
            value = record["scores"].get(category)
            x = left + column_index * cell_width
            fill = score_color(value)
            display = "" if value is None else str(value)
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell_width - 4}" '
                f'height="{cell_height - 4}" fill="{fill}" rx="4">'
                f'<title>{escape(record["model"])} {escape(category)}: {display or "n/a"}</title>'
                "</rect>"
            )
            parts.append(
                f'<text x="{x + (cell_width - 4) / 2:.1f}" y="{y + 21}" '
                f'class="cell-value" text-anchor="middle">{display}</text>'
            )

    parts.append("</svg>")
    return "\n".join(parts)


def build_rankings(records: list[EvaluationRecord]) -> dict[str, list[dict[str, Any]]]:
    """Rank successful evaluations by project."""
    grouped: dict[str, list[EvaluationRecord]] = {}
    for record in records:
        grouped.setdefault(record.project, []).append(record)

    rankings: dict[str, list[dict[str, Any]]] = {}
    for project, project_records in grouped.items():
        ranked = sorted(
            project_records,
            key=lambda item: item.overall_score if item.overall_score is not None else -1,
            reverse=True,
        )
        rankings[project] = [
            {
                "rank": index + 1,
                "model": record.model,
                "status": record.status,
                "overall_score": record.overall_score,
            }
            for index, record in enumerate(ranked)
        ]
    return rankings


def safe_filename(value: str) -> str:
    """Return a conservative filename stem for generated artifacts."""
    safe = "".join(char if char.isalnum() or char in "-_" else "_" for char in value)
    return safe.strip("_") or "project"


def short_category_label(category: str) -> str:
    labels = {
        "correctness": "Correct",
        "completeness": "Complete",
        "specificity": "Specific",
        "coherence": "Coherent",
        "actionability": "Actionable",
        "reference_alignment": "Reference",
    }
    return labels.get(category, category.replace("_", " ").title())


def chart_palette(count: int) -> list[str]:
    palette = [
        "#2563eb",
        "#dc2626",
        "#16a34a",
        "#9333ea",
        "#ea580c",
        "#0891b2",
        "#be123c",
        "#4f46e5",
    ]
    return [palette[index % len(palette)] for index in range(count)]


def score_color(value: int | None) -> str:
    colors = {
        1: "#fee2e2",
        2: "#fed7aa",
        3: "#fef08a",
        4: "#bbf7d0",
        5: "#86efac",
    }
    return colors.get(value, "#e5e7eb")


def svg_header(width: int | float, height: int | float) -> str:
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">
<style>
  .title {{ font: 700 18px Arial, sans-serif; fill: #111827; }}
  .axis-label {{ font: 12px Arial, sans-serif; fill: #374151; }}
  .row-label {{ font: 13px Arial, sans-serif; fill: #111827; }}
  .cell-value {{ font: 700 13px Arial, sans-serif; fill: #111827; }}
  .tick {{ font: 11px Arial, sans-serif; fill: #6b7280; }}
  .grid {{ stroke: #e5e7eb; stroke-width: 1; }}
</style>
<rect width="100%" height="100%" fill="#ffffff"/>"""


def render_y_axis(
    margin_left: int,
    margin_top: int,
    plot_height: int,
    plot_width: int,
) -> list[str]:
    parts: list[str] = []
    for tick in range(1, 6):
        y = margin_top + plot_height - (plot_height * tick / 5)
        parts.append(
            f'<line x1="{margin_left}" y1="{y:.1f}" '
            f'x2="{margin_left + plot_width}" y2="{y:.1f}" class="grid"/>'
        )
        parts.append(
            f'<text x="{margin_left - 10}" y="{y + 4:.1f}" '
            f'class="tick" text-anchor="end">{tick}</text>'
        )
    parts.append(
        f'<line x1="{margin_left}" y1="{margin_top}" '
        f'x2="{margin_left}" y2="{margin_top + plot_height}" stroke="#9ca3af"/>'
    )
    parts.append(
        f'<line x1="{margin_left}" y1="{margin_top + plot_height}" '
        f'x2="{margin_left + plot_width}" y2="{margin_top + plot_height}" '
        'stroke="#9ca3af"/>'
    )
    return parts


def render_legend(
    records: list[dict[str, Any]],
    colors: list[str],
    x_start: int,
    y_start: int,
    width: int,
) -> list[str]:
    parts: list[str] = []
    x = x_start
    y = y_start
    max_x = width - 180
    for index, record in enumerate(records):
        label = record["model"]
        item_width = max(130, min(220, 42 + len(label) * 7))
        if x + item_width > max_x and x != x_start:
            x = x_start
            y += 24
        parts.append(
            f'<rect x="{x}" y="{y - 12}" width="12" height="12" '
            f'fill="{colors[index]}" rx="2"/>'
        )
        parts.append(
            f'<text x="{x + 18}" y="{y - 2}" class="axis-label">'
            f'{escape(label)}</text>'
        )
        x += item_width
    return parts


def render_markdown_summary(result: dict[str, Any]) -> str:
    """Render a compact human-readable report."""
    lines = ["# LLM Spec Evaluation Summary", ""]

    for project, ranking in result["rankings"].items():
        lines.extend([f"## {project}", "", "| Rank | Model | Status | Overall |"])
        lines.append("| ---: | --- | --- | ---: |")
        for item in ranking:
            overall = "" if item["overall_score"] is None else f"{item['overall_score']:.2f}"
            lines.append(
                f"| {item['rank']} | {item['model']} | {item['status']} | {overall} |"
            )
        lines.append("")
        plot_paths = result.get("plots", {}).get(project, {})
        if plot_paths:
            if "category_bars" in plot_paths:
                lines.extend(
                    [
                        f"![Category score bars]({plot_paths['category_bars']})",
                        "",
                    ]
                )
            if "score_heatmap" in plot_paths:
                lines.extend(
                    [
                        f"![Score heatmap]({plot_paths['score_heatmap']})",
                        "",
                    ]
                )

        records = [
            item
            for item in result["evaluations"]
            if item["project"] == project
        ]
        for record in sorted(
            records,
            key=lambda item: item["overall_score"] if item["overall_score"] is not None else -1,
            reverse=True,
        ):
            lines.extend([f"### {record['model']}", ""])
            if record["status"] != "success":
                lines.extend([f"Error: {record['error']}", ""])
                continue
            lines.append(f"Overall score: {record['overall_score']:.2f}")
            lines.append("")
            lines.append(record["justification"])
            lines.append("")
            lines.append("Strengths:")
            lines.extend(f"- {item}" for item in record["strengths"])
            lines.append("")
            lines.append("Weaknesses:")
            lines.extend(f"- {item}" for item in record["weaknesses"])
            lines.append("")
            lines.append("Suggested improvements:")
            lines.extend(f"- {item}" for item in record["suggested_improvements"])
            lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def _require_string(response: dict[str, Any], key: str) -> str:
    value = response.get(key)
    if not isinstance(value, str):
        raise EvaluationError(f"{key!r} must be a string")
    return value


def _require_string_list(response: dict[str, Any], key: str) -> list[str]:
    value = response.get(key)
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise EvaluationError(f"{key!r} must be a list of strings")
    return value
