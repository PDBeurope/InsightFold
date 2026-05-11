#!/usr/bin/env python3
"""Create an InsightFold notebook spec-pack skeleton."""

from __future__ import annotations

import argparse
from pathlib import Path


FILES = {
    "spec-pack-overview.md": "# {title} Spec Pack Overview\n\n| Field | Value |\n|---|---|\n| Source PRD | {prd} |\n| Status | Draft |\n| Target Notebook | notebooks/{slug}.ipynb |\n\n## Summary\n\n## Blocking Questions\n\n## Assumptions Accepted For Build\n\n## Required Advisory Reviews\n",
    "requirements.md": "# Requirements\n\n## User Story\n\n## Functional Requirements\n\n| ID | Requirement | Acceptance Criteria | Priority |\n|---|---|---|---|\n| REQ-001 | | | Must |\n\n## Non-Functional Requirements\n\n## Non-Goals\n\n## Open Questions\n",
    "notebook-ux-contract.md": "# Notebook UX Contract\n\n## Primary User Question\n\n## Target User And Expertise\n\n## First Runnable User Input Cell\n\n## Happy-Path Default Example\n\n## User Flow\n\n## Fixtures Versus User Flow\n\n## Trust, Caveats, And Limitations To Show In Notebook\n",
    "notebook-design.md": "# Notebook Design\n\n## Section Outline\n\n## Data Flow\n\n## Variable And Artifact Handoff\n\n## Dependency Policy\n\n## Visualization Plan\n\nFor 3D molecular views, use MolViewSpec/Mol* through `$molviewspec-rendering`.\n\n## Design Decisions\n",
    "cell-blueprint.md": "# Jupyter Cell Blueprint\n\n| Cell ID | Type | Purpose | Inputs | Outputs | User Editable? | Hidden-State Risk | Validation Hook | Requirements |\n|---|---|---|---|---|---|---|---|---|\n| C001 | code | | | | yes / no | | | REQ-001 |\n",
    "traceability-matrix.md": "# Requirement Traceability Matrix\n\n| Requirement | PRD Goal | Notebook Section / Cell | Fixture | Validation Check | Review Evidence |\n|---|---|---|---|---|---|\n| REQ-001 | G-001 | | | | |\n",
    "data-contracts.md": "# Data Contracts\n\n| Source / Artifact | Required Fields | Optional Fields | Validation | Failure Behavior | Provenance |\n|---|---|---|---|---|---|\n| | | | | | |\n",
    "fixture-manifest.md": "# Fixture Manifest\n\n| Fixture ID | Source | Why Chosen | Expected Outputs | Cache / Version | Notes |\n|---|---|---|---|---|---|\n| happy-path | | | | | |\n| edge-case | | | | | |\n",
    "validation.md": "# Validation Plan\n\n| Check | Method | Pass Criteria | Fixture / Artifact | Requirement |\n|---|---|---|---|---|\n| Restart-and-run-all | Restart kernel and run all | No errors | Notebook | NFR-001 |\n",
    "tasks.md": "# Tasks\n\n## Phase A: Spec Preparation\n\n- [ ] T001 Resolve blocking open questions\n\n## Phase B: Data And Fixtures\n\n- [ ] T010 Create fixture manifest\n- [ ] T011 Document data contracts\n\n## Phase C: Notebook Architecture\n\n- [ ] T020 Create section skeleton\n\n## Phase D: Implementation\n\n- [ ] T030 Implement first executable section\n\n## Phase E: Validation\n\n- [ ] T040 Run notebook top-to-bottom on pinned fixture\n\n## Phase F: Documentation\n\n- [ ] T050 Add user-facing run instructions\n",
    "docs-plan.md": "# Documentation Plan\n\n| Documentation Type | Artifact | Audience | Notes |\n|---|---|---|---|\n| Tutorial | Notebook intro / README | New users | |\n| How-to | Adaptation section | Returning users | |\n| Reference | Inputs, outputs, parameters, APIs | Developers / power users | |\n| Explanation | Scientific background and limitations | All users | |\n",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("slug", help="Feature slug, e.g. protein_model_chem")
    parser.add_argument("--prd", default="", help="Source PRD path")
    parser.add_argument("--root", default="specs", help="Spec root directory")
    parser.add_argument("--force", action="store_true", help="Overwrite existing files")
    args = parser.parse_args()

    target = Path(args.root) / args.slug
    target.mkdir(parents=True, exist_ok=True)
    title = args.slug.replace("_", " ").replace("-", " ").title()

    for name, content in FILES.items():
        path = target / name
        if path.exists() and not args.force:
            continue
        path.write_text(content.format(slug=args.slug, title=title, prd=args.prd), encoding="utf-8")

    print(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

