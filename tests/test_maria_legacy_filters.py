"""Keep Maria's existing single-substitution reporting contract alongside the new path."""
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("inspect_scores", ROOT / "filters/inspect_scores.py")
inspect = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inspect)


def test_mutation_rows_collapse_without_summing_candidate_scores():
    data = pd.DataFrame([
        {"variant_id": "pair", "chain": "B", "position": 1, "wt": "A", "mut": "C", "score_seq": -5},
        {"variant_id": "pair", "chain": "A", "position": 2, "wt": "D", "mut": "E", "score_seq": -5},
    ])
    frame, _ = inspect.make_one_row_per_variant(data)
    assert len(frame) == 1
    assert frame.iloc[0]["score_seq"] == -5
    assert frame.iloc[0]["mutation_count"] == 2


def test_conflicting_repeated_scores_are_rejected():
    data = pd.DataFrame([
        {"variant_id": "pair", "chain": "B", "position": 1, "wt": "A", "mut": "C", "score_seq": -5},
        {"variant_id": "pair", "chain": "A", "position": 2, "wt": "D", "mut": "E", "score_seq": -3},
    ])
    with pytest.raises(ValueError, match="conflicting"):
        inspect.make_one_row_per_variant(data)


def test_single_mutation_filter_keeps_missing_and_uncalibrated_distinct(tmp_path):
    source, config, out = tmp_path / "scores.csv", tmp_path / "config.toml", tmp_path / "report"
    source.write_text("variant_id,mutation_count,score_seq,score_inv\na,1,1,1\nb,1,,1\nc,1,-1,1\n")
    config.write_text('''[filter]
minimum_sequence_positive = 1
minimum_inverse_positive = 1
require_top_at_site = false
[[sequence_models]]
name = "sequence"
score_column = "score_seq"
auroc = "not_measured"
[[inverse_folding_models]]
name = "inverse"
score_column = "score_inv"
auroc = "not_measured"
''')
    result = subprocess.run([sys.executable, str(ROOT / "filters/filter_variants.py"),
                             "--input", str(source), "--config", str(config), "--output-dir", str(out)],
                            text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    decisions = pd.read_csv(out / "filter_decisions.csv").set_index("variant_id")
    assert decisions["filter_status"].to_dict() == {"a": "PROVISIONAL_PASS", "b": "INCOMPLETE", "c": "PROVISIONAL_FAIL"}
