"""Shared notebook bootstrap: the part that can run *after* `src/` is importable.

Every analysis notebook here opens with the same problem. It must work unchanged
from a local checkout and from a bare Colab VM, and the cell that solves it is
the one cell that cannot use this package, because its whole job is to make the
package importable. That cell therefore keeps only what must run **before** the
import, which is finding a checkout or making one, and hands everything that can
run after it to `setup()` below.

Splitting it this way is what keeps the notebooks' first cell short. Before this
module the same logic was written out three times, at 120, 61 and 27 lines, and
had already drifted: one copy still carried a placeholder organisation in its
clone URL and another still pinned a branch that had long since merged.

**The per-module hook.** `setup()` imports the analysis module and then calls its
`prepare_environment(repo_root, branch=..., colab=...)` if it defines one. That is
where a module installs its own optional packages, applies its own plot style and
prints its own banner. Keeping it there rather than here means this module needs
to know nothing about any particular notebook, and each notebook's output is
decided by the module it already depends on.

Deliberately **not** `pip install`-ing InsightFold itself (D9): that resolves
`pyproject.toml` and drags in biopython, gemmi, scipy and plotly, breaking both
the dependency rules and the 60 s Colab install budget. The repo is public, so
there is no token, no auth header and no `getpass`, which would block forever in
a Run-all notebook.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType

REPO_URL = "https://github.com/PDBeurope/InsightFold.git"
COLAB_CLONE_DIR = Path("/content/InsightFold")


def in_colab() -> bool:
    return "google.colab" in sys.modules


def find_repo_root(start=None) -> Path | None:
    """Walk up from `start` for a directory holding both `pyproject.toml` and `src/`.

    Duplicated as a few inline lines in each notebook's first cell, because it has
    to run before this module can be imported. Kept here as well so that anything
    running after the bootstrap has one definition to call.
    """
    here = (Path(start) if start is not None else Path.cwd()).expanduser().resolve()
    return next((d for d in (here, *here.parents)
                 if (d / "pyproject.toml").is_file() and (d / "src").is_dir()), None)


def add_to_path(repo_root) -> str:
    """Put `repo_root/src` on `sys.path`, once.

    Guarded against duplicates: running the notebook twice in one kernel must not
    grow `sys.path`, which NFR-002 checks for.
    """
    src = str(Path(repo_root) / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    return src


def setup(module: str, repo_root, branch: str = "", *, inline: bool = True) -> ModuleType:
    """Finish the bootstrap and hand back the analysis module.

    Args:
        module:    Name under `insightfold`, for example `"cluster_quality_utils"`.
        repo_root: The checkout in use, as the notebook's first cell resolved it.
        branch:    Branch the notebook expects, for the banner.
        inline:    Select the inline matplotlib backend. A notebook needs this for
                   the figure formatters to be registered, and a module can do it
                   through `run_line_magic` just as a `%matplotlib inline` cell can,
                   which is what lets the notebook drop that line.

    Returns:
        The imported module, so the first cell can bind it in one statement.
    """
    add_to_path(repo_root)

    if inline:
        try:
            from IPython import get_ipython

            shell = get_ipython()
            if shell is not None:
                shell.run_line_magic("matplotlib", "inline")
        except Exception:
            pass        # not under IPython, or no display: the import still matters

    try:
        analysis = importlib.import_module(f"insightfold.{module}")
    except ImportError as exc:
        raise ImportError(
            f"Could not import 'insightfold.{module}': {exc}\n"
            f"  repo root in use: {repo_root}\n"
            f"  branch expected:  {branch or '(unspecified)'}\n"
            "The likely cause is a checkout that predates this module or sits on a "
            "different branch.\n"
            f"    git -C {repo_root} switch {branch or 'main'} && git -C {repo_root} pull\n"
            f"or, on Colab, `!rm -rf {COLAB_CLONE_DIR}` and re-run the first cell."
        ) from exc

    finish = getattr(analysis, "prepare_environment", None)
    if callable(finish):
        finish(repo_root, branch=branch, colab=in_colab())
    return analysis
