"""The end-to-end notebook still matches the API it calls.

`notebooks/03_end_to_end.ipynb` needs a GPU and the model environments, so there
is no offline mode and nothing here runs it - running it is what the GPU box is
for. What a laptop can check is the thing that actually rots: the notebook is
valid JSON, every cell compiles, and every name and keyword argument it reaches
for still exists in the package. A rename in `api.py` fails here rather than
forty minutes into a campaign.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path

import pytest

NOTEBOOK = Path(__file__).resolve().parents[1] / "notebooks" / "03_end_to_end.ipynb"


def code_cells() -> list[str]:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return [
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    ]


def test_notebook_is_valid_json_with_cells() -> None:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    assert notebook["nbformat"] == 4
    assert [c["cell_type"] for c in notebook["cells"]].count("markdown") >= 10
    assert len(code_cells()) >= 10


@pytest.mark.parametrize("index", range(len(code_cells())))
def test_every_cell_compiles(index: int) -> None:
    compile(code_cells()[index], f"cell{index}", "exec")


def test_parameters_cell_defines_every_parameter() -> None:
    """The first cell is the only one meant to be edited, so it is the only one
    allowed to introduce a parameter. A name used below but missing here is a
    parameter someone has to go hunting for."""
    namespace: dict = {}
    exec(code_cells()[0], namespace)
    for name in (
        "ANTIBODY", "ANTIGEN", "ANTIGEN_REGION", "EPITOPE", "PAD", "PREDICTOR",
        "NUM_SEEDS", "DECOY_ANTIGEN", "SCAN_RECORD", "POSITIONS", "EXCLUDE_LINKER",
        "EXCLUDE_CYSTEINE", "SCORERS", "RANK_BY", "MAX_PER_VARIANT", "GATE_CHECKS",
        "MAX_PAIRS", "MIN_DISTANCE", "PAIR_SCORERS", "TOP_SINGLES", "VETO_SEEDS",
        "VETO_SD_MULTIPLE", "K", "MAX_SHARED_POSITIONS", "SUBMIT_LINKER",
        "LENGTH_LIMIT", "RUN",
    ):
        assert name in namespace, f"the parameters cell does not define {name}"


def test_named_models_are_registered() -> None:
    """A typo in SCORERS would otherwise surface as an UnknownModelError on the
    GPU box, after the parent fold has already been paid for."""
    from design_pipeline.registry import REGISTRY

    namespace: dict = {}
    exec(code_cells()[0], namespace)
    for name in [
        namespace["PREDICTOR"],
        *namespace["SCORERS"],
        *namespace["PAIR_SCORERS"],
        *namespace["GATE_CHECKS"],
    ]:
        assert name in REGISTRY, f"{name} is not a registered adapter"


def test_the_functions_it_calls_take_the_arguments_it_passes() -> None:
    from design_pipeline import api, variants

    expected = {
        api.get_structure: (
            "model_name", "antibody_seq", "antigen_seq", "num_seeds",
            "decoy_antigen_seq",
        ),
        api.score_variants: ("output_dir", "models", "max_variants"),
        api.gate: ("checks", "output_dir"),
        variants.enumerate_variants: ("chain", "positions", "output_dir"),
        variants.select_pairs: ("max_pairs", "min_distance", "output_dir", "by"),
        variants.mutant_sequence: ("structure", "mutations", "chain"),
    }
    for function, names in expected.items():
        parameters = inspect.signature(function).parameters
        for name in names:
            assert name in parameters, f"{function.__name__} has no {name}="


def test_the_chain_cell_passes_only_flags_the_scorers_accept() -> None:
    """The kwargs the notebook builds become command-line flags verbatim, so a
    scorer whose CLI changes turns them into an argparse error on the GPU box.
    This is the check that catches it here instead: the chain cell is run against
    the committed 1N8Z and every flag it produces is looked up in the parser of
    the script that scorer actually starts.
    """
    import re

    repo = Path(__file__).resolve().parents[1]
    pdb = repo / "data" / "structures" / "1N8Z.pdb"
    if not pdb.is_file():
        pytest.skip(f"{pdb} is not here; it is the only structure this needs")

    from design_pipeline import Complex

    from design_pipeline import structures

    # 1N8Z: heavy B, light A, antigen C - design_pipeline.benchmark says so, and
    # the cell under test is the one that works that out from the sequences.
    sequences = {c: structures.chain_sequence(pdb, c) for c in ("A", "B", "C")}

    cells = code_cells()
    namespace: dict = {"display": lambda obj: None, "structures": structures}
    exec(cells[0], namespace)  # the parameters cell, for ANTIFOLD_SRC
    namespace.update(
        {
            "complex_": Complex(paths=(pdb,), output_dir=repo / "runs" / "unused"),
            "antibody_records": [("VH", sequences["B"]), ("VL", sequences["A"])],
            "antigen_name": "HER2",
            "antigen": sequences["C"],
            "SCAN_RECORD": 0,
        }
    )

    exec(cells[4], namespace)  # the chain-roles cell
    assert namespace["ab_chains"] == ["B", "A"]
    assert namespace["antigen_chain"] == "C"

    # Which script each scorer's flags are read by. `antifold` is the one that
    # moved: its adapter runs antifold_score, not `python -m antifold.main`.
    parsers = {
        "antifold": repo / "antifold_score" / "scorer.py",
        "antifold_seq": repo / "design_pipeline" / "adapters" / "_antifold_seq.py",
    }
    for model, script in parsers.items():
        flags = set(re.findall(r'add_argument\(\s*"(--[A-Za-z0-9_]+)"', script.read_text()))
        assert flags, f"no flags found in {script}"
        for key in namespace["SCORER_KWARGS"].get(model, {}):
            assert f"--{key}" in flags, (
                f"the notebook passes {key}= to {model}, which starts {script.name}; "
                f"that script takes {', '.join(sorted(flags))}"
            )


def test_the_helpers_it_borrows_still_exist() -> None:
    """The notebook loads examples/track_a_fold.py by path and reuses its
    sequence reading and antigen trimming rather than writing them again."""
    import importlib.util

    script = Path(__file__).resolve().parents[1] / "examples" / "track_a_fold.py"
    spec = importlib.util.spec_from_file_location("track_a_fold", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for name in (
        "read_sequences", "parse_region", "parse_residues", "trim_antigen",
        "warnings_for", "STREP_MOTIF",
    ):
        assert hasattr(module, name), f"track_a_fold.{name} is gone"

    from design_pipeline import confidence, preflight, structures, variants

    for module_, name in (
        (preflight, "rows"), (preflight, "render"), (preflight, "notes"),
        (confidence, "warn"), (structures, "chain_ids"), (structures, "chain_sequence"),
        (variants, "mutant_sequence"), (variants, "parse_id"), (variants, "read_table"),
    ):
        assert hasattr(module_, name), f"{module_.__name__}.{name} is gone"


def test_the_linker_exclusion_finds_a_g4s_linker() -> None:
    """EXCLUDE_LINKER drops the (G4S)n span of an scFv, which is not a binding
    position. The regex is in the scan cell; this is the case it is for."""
    import re

    sequence = "EVQL" + "GGGGS" * 3 + "DIQM"
    found = {
        i
        for match in re.finditer(r"(?:G{3,4}S){2,}", sequence)
        for i in range(match.start() + 1, match.end() + 1)
    }
    assert found == set(range(5, 20))
    assert "".join(sequence[i - 1] for i in sorted(found)) == "GGGGS" * 3
