"""examples/maturation.py, the offline run: all six steps, including the scoring of the pairs.

The real run needs AntiFold and a structure; this is the `--dummy` path, which goes through the same six calls with
the fixture scorers, so a change to how step 6 is called is exercised without installing a model.
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "maturation.py"


def test_the_dummy_run_scores_the_pairs(tmp_path, monkeypatch, capsys) -> None:
    monkeypatch.chdir(tmp_path)  # the example writes runs/maturation_dummy under the working directory
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--dummy"])

    runpy.run_path(str(SCRIPT), run_name="__main__")
    printed = capsys.readouterr().out

    assert "pairs:" in printed
    assert "pairs by score_dummy_scorer:" in printed
    assert (tmp_path / "runs" / "maturation_dummy").is_dir()
