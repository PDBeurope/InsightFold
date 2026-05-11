#!/usr/bin/env python3
"""Check that an InsightFold notebook spec pack has the expected files."""

from __future__ import annotations

import argparse
from pathlib import Path


REQUIRED = [
    "spec-pack-overview.md",
    "requirements.md",
    "notebook-ux-contract.md",
    "notebook-design.md",
    "cell-blueprint.md",
    "traceability-matrix.md",
    "data-contracts.md",
    "fixture-manifest.md",
    "validation.md",
    "tasks.md",
    "docs-plan.md",
]

REQUIRED_TEXT = {
    "requirements.md": ["REQ-"],
    "notebook-ux-contract.md": ["First Runnable User Input Cell"],
    "notebook-design.md": ["Dependency Policy", "Visualization Plan"],
    "cell-blueprint.md": ["Hidden-State Risk"],
    "traceability-matrix.md": ["Requirement", "Validation Check"],
    "fixture-manifest.md": ["happy-path"],
    "validation.md": ["Restart-and-run-all"],
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("spec_dir")
    args = parser.parse_args()

    spec_dir = Path(args.spec_dir)
    errors = []

    for name in REQUIRED:
        path = spec_dir / name
        if not path.exists():
            errors.append(f"missing {name}")
            continue
        text = path.read_text(encoding="utf-8")
        for needle in REQUIRED_TEXT.get(name, []):
            if needle not in text:
                errors.append(f"{name} missing marker: {needle}")

    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1

    print(f"PASS: {spec_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
