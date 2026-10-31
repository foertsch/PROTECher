"""The scripts that call the package still import names it exports.

The slurm jobs, the examples and `msa_cache` are not imported by any other test, and
the slurm jobs carry their Python in a heredoc, so nothing here notices when a name
leaves `design_pipeline/__init__.py`: the job fails with an ImportError on the GPU
node after the queue wait. That is how `predict_structure` broke six jobs on 30 Oct.
What a laptop can check is every `from design_pipeline import ...` in those files
against what the package offers.
"""

from __future__ import annotations

import ast
import importlib.util
import re
from pathlib import Path

import pytest

import design_pipeline

REPO = Path(__file__).resolve().parents[1]

# Code that calls the package without being the package or a test.
SCRIPT_GLOBS = (
    "slurm/*.sh",
    "examples/*.py",
    "msa_cache/*.py",
    "cloud/*.py",
    "run_report/*.py",
    "antifold_score/*.py",
)

# `from design_pipeline import a, b`, or the parenthesised form over several lines.
# Matching the statement instead of parsing the file is what lets one pattern read a
# heredoc in a shell script and a lazy import inside a function alike.
TOP_LEVEL_IMPORT = re.compile(
    r"^[ \t]*from design_pipeline import (\([^)]*\)|[^\n]*)", re.MULTILINE
)


def top_level_imports(script: Path) -> list[tuple[int, str]]:
    """(line, name) for each name the script takes from the package top level."""
    text = script.read_text(encoding="utf-8")
    found = []
    for match in TOP_LEVEL_IMPORT.finditer(text):
        line = text.count("\n", 0, match.start()) + 1
        statement = ast.parse("from design_pipeline import " + match.group(1))
        found.extend((line, alias.name) for alias in statement.body[0].names)
    return found


def importable(name: str) -> bool:
    """What `from design_pipeline import name` needs: an attribute, or a submodule
    to import. `find_spec` does not run the submodule, so nothing heavy loads."""
    return (
        hasattr(design_pipeline, name)
        or importlib.util.find_spec(f"design_pipeline.{name}") is not None
    )


SCRIPTS = sorted(
    script
    for pattern in SCRIPT_GLOBS
    for script in REPO.glob(pattern)
    if top_level_imports(script)
)


def test_the_scan_reaches_the_scripts_it_is_for() -> None:
    """A glob that moved or a pattern that broke would match nothing and pass."""
    scanned = {script.relative_to(REPO).as_posix() for script in SCRIPTS}
    assert {
        "slurm/example_gpu_job.sh",  # a heredoc
        "msa_cache/cache.py",  # an import inside a function
        "examples/maturation.py",  # the parenthesised form
    } <= scanned


def test_importable_tells_an_export_and_a_submodule_from_a_missing_name() -> None:
    assert importable("Design") and importable("structures")
    assert not importable("no_such_name")


@pytest.mark.parametrize(
    "script", SCRIPTS, ids=lambda script: script.relative_to(REPO).as_posix()
)
def test_every_name_taken_from_the_package_exists(script: Path) -> None:
    missing = [
        f"line {line}: {name}"
        for line, name in top_level_imports(script)
        if not importable(name)
    ]
    assert not missing, (
        f"{script.relative_to(REPO)} imports from design_pipeline what it does not "
        f"export ({'; '.join(missing)}). Import it from the module that defines it, "
        "or export it again on purpose."
    )
