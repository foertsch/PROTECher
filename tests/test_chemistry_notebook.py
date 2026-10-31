"""The chemistry-scoring notebook stays valid and in step with the package.

`notebooks/05_chemistry_scoring.ipynb` needs the GPU box's Boltz folds to run, so there is no offline execution here. What a
laptop checks is what rots: the notebook is valid JSON, every code cell compiles, and the parameters cell defines every name the
notebook relies on. A rename in `scoring.py` that breaks the import surfaces here, not in the middle of a pitch.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

NOTEBOOK = Path(__file__).resolve().parents[1] / "notebooks" / "05_chemistry_scoring.ipynb"


def code_cells() -> list[str]:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return ["".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code"]


def test_notebook_is_valid_json_with_cells() -> None:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    assert notebook["nbformat"] == 4
    assert [c["cell_type"] for c in notebook["cells"]].count("markdown") >= 4
    assert len(code_cells()) >= 4


@pytest.mark.parametrize("index", range(len(code_cells())))
def test_every_cell_compiles(index: int) -> None:
    compile(code_cells()[index], f"cell{index}", "exec")


def test_parameters_cell_defines_every_parameter() -> None:
    namespace: dict = {}
    exec(code_cells()[0], namespace)
    for name in ("RUN", "PARENT", "SEEDS", "CRYSTAL"):
        assert name in namespace, f"the parameters cell does not define {name}"


def test_it_imports_the_metric_the_package_exports() -> None:
    """The notebook calls charge_contact_report; a rename would strand it."""
    from design_pipeline.scoring import charge_contact_report  # noqa: F401
