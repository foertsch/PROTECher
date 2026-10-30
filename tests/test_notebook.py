"""The validation notebook executes offline, with no GPU and no network.

Running it is the only way to know it still works, and jupyter is not a
dependency of this package - so the code cells are extracted and executed in
order, which is what a kernel does. `display` is the one Jupyter builtin they
use, so the test supplies it.

The parameters cell is overridden after it runs: `FORCE_DUMMY=True` keeps the
test off RCSB and SAbDab even on a machine that has a GPU and a network.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import pytest

NOTEBOOK = Path(__file__).resolve().parents[1] / "notebooks" / (
    "01_structure_pipeline_validation.ipynb"
)

PARAMETERS = (
    "PREDICTOR",
    "NUM_SEEDS",
    "LIMIT",
    "K",
    "OUTPUT_DIR",
    "SABDAB",
    "TRAINING_CUTOFF",
    "FORCE_DUMMY",
)


def code_cells() -> list[str]:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return [
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    ]


def test_the_notebook_is_valid_and_documented() -> None:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))

    assert notebook["nbformat"] == 4
    markdown = [c for c in notebook["cells"] if c["cell_type"] == "markdown"]
    assert len(markdown) >= 5  # it is meant to be read, not only run
    prose = "".join("".join(c["source"]) for c in markdown).lower()
    # The distinction the whole notebook turns on.
    assert "memorized" in prose
    assert "honest" in prose
    assert "epitope" in prose


def test_the_parameters_are_all_in_one_cell() -> None:
    """`TRAINING_CUTOFF` is a parameter, not a constant, and the rest have to be
    reachable without reading the whole notebook."""
    first = code_cells()[0]

    for name in PARAMETERS:
        assert f"{name} =" in first, name
    assert "2021-09-30" in first  # AlphaFold3's cutoff, cited where it is set


def test_the_notebook_runs_end_to_end_on_fixtures(tmp_path) -> None:
    cells = code_cells()
    environment: dict[str, object] = {
        "display": lambda *args, **kwargs: None,
        "__name__": "__main__",
    }

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for index, source in enumerate(cells):
            exec(compile(source, f"<notebook cell {index}>", "exec"), environment)
            if index == 0:
                # Keep the run offline and small whatever the machine has.
                environment["FORCE_DUMMY"] = True
                environment["OUTPUT_DIR"] = str(tmp_path / "runs")
                environment["NUM_SEEDS"] = 4

    assert environment["DUMMY"] is True
    results = environment["results"]
    assert len(results) == 3  # three synthetic targets
    assert set(results["tier"]) == {"fixture"}
    for column in (
        "epitope_recall",
        "interface_rmsd",
        "confidence_iptm",
        "confidence_dominant_cluster_fraction",
        "confidence_decoy_gap",
    ):
        assert column in results.columns, column


def test_a_fixture_run_refuses_to_give_a_verdict(tmp_path, capsys) -> None:
    """The degraded mode exists so this executes in CI, not so it produces a
    result, and the notebook has to say so rather than print a number."""
    cells = code_cells()
    environment: dict[str, object] = {
        "display": lambda *args, **kwargs: None,
        "__name__": "__main__",
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for index, source in enumerate(cells):
            exec(compile(source, f"<notebook cell {index}>", "exec"), environment)
            if index == 0:
                environment["FORCE_DUMMY"] = True
                environment["OUTPUT_DIR"] = str(tmp_path / "runs")
                environment["NUM_SEEDS"] = 4

    printed = capsys.readouterr().out
    assert "meaningless" in printed
    assert "FIXTURE RUN - no verdict" in printed


def test_the_honest_tier_refuses_an_unknown_training_cutoff(tmp_path) -> None:
    """Anything deposited before the cutoff may be training data, so a missing
    cutoff stops the honest tier instead of quietly selecting everything."""
    cells = code_cells()
    environment: dict[str, object] = {
        "display": lambda *args, **kwargs: None,
        "__name__": "__main__",
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exec(compile(cells[0], "<cell 0>", "exec"), environment)
        # A predictor with no recorded cutoff, and pretend we have a GPU.
        environment["FORCE_DUMMY"] = False
        environment["PREDICTOR"] = "boltz2"
        environment["TRAINING_CUTOFF"] = None
        environment["OUTPUT_DIR"] = str(tmp_path / "runs")
        exec(compile(cells[1], "<cell 1>", "exec"), environment)
        assert environment["cutoff"] is None
        exec(compile(cells[2], "<cell 2>", "exec"), environment)  # smoke tier, skipped

        with pytest.raises(RuntimeError, match="no training cutoff known"):
            exec(compile(cells[3], "<cell 3>", "exec"), environment)
