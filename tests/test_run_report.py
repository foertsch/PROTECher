"""run_report: loading manifests, grading, report files, comparison across runs, leaderboard."""

import json
from pathlib import Path

import pandas as pd
import pytest

from run_report.__main__ import DEFAULT_THRESHOLDS, main
from run_report.config import Metric, load_config
from run_report.leaderboard import build_leaderboard, load_summaries
from run_report.load import load_designs
from run_report.scorecard import evaluate, find_parent, grade


def write_run(root: Path, name: str, rows: list[dict], producer: str = "boltz") -> Path:
    """A run directory shaped like design_pipeline's: runs/<name>/predictions/<key>/DONE."""
    run = root / name
    step = run / "predictions" / "abc123"
    step.mkdir(parents=True)
    records = [
        {
            "id": r["id"],
            "path": f"{r['id']}.pdb",
            "parent_id": None,
            "producer": producer,
            "scores": {k: v for k, v in r.items() if k not in ("id", "meta")},
            "meta": {"seed": 0, "ancestors": ["x"], **r.get("meta", {})},
        }
        for r in rows
    ]
    (step / "DONE").write_text(json.dumps(records))
    return run


ROWS = [
    {"id": "wt", "iptm": 0.90, "plddt": 90.0},
    {"id": "v1", "iptm": 0.85, "plddt": 88.0, "meta": {"group": "stability", "mutations": "VH:K30D+VL:S52D"}},
    {"id": "v2", "iptm": 0.65, "plddt": 82.0, "meta": {"group": "affinity"}},
    {"id": "v3", "iptm": 0.40, "plddt": 91.0, "meta": {"group": "affinity"}},
    {"id": "v4", "iptm": 0.82, "plddt": 70.0, "meta": {"group": "stability"}},
]


def test_grade_higher_is_better():
    m = Metric("iptm", "ipTM", "higher", target=0.8, gate=0.5)
    assert [grade(v, m) for v in (0.9, 0.8, 0.6, 0.5, 0.4)] == ["ok", "ok", "warn", "warn", "fail"]


def test_grade_lower_is_better_and_missing():
    m = Metric("pae", "PAE", "lower", target=6, gate=10)
    assert [grade(v, m) for v in (5, 6, 8, 10, 12)] == ["ok", "ok", "warn", "warn", "fail"]
    assert grade(float("nan"), m) == "n/a"
    assert grade(None, m) == "n/a"


def test_grade_gate_only_and_target_only():
    assert grade(0.7, Metric("a", "A", gate=0.5)) == "ok"
    assert grade(0.4, Metric("a", "A", gate=0.5)) == "fail"
    assert grade(0.7, Metric("a", "A", target=0.8)) == "warn"


def test_metric_rejects_bad_direction():
    with pytest.raises(ValueError):
        Metric("a", "A", better="up")


def test_load_designs_flattens_scores_and_meta(tmp_path):
    run = write_run(tmp_path, "r1", ROWS)
    df = load_designs([run])
    assert set(df["id"]) == {"wt", "v1", "v2", "v3", "v4"}
    row = df.set_index("id").loc["v1"]
    assert row["iptm"] == 0.85 and row["meta_group"] == "stability" and row["step"] == "predictions"
    assert "meta_ancestors" not in df.columns  # lists are not scalar


def test_designs_meta_sidecar_adds_mutations_and_group(tmp_path):
    run = write_run(tmp_path, "r1", [{"id": "p_model_0", "iptm": 0.9, "plddt": 90.0}, {"id": "q_model_0", "iptm": 0.8, "plddt": 85.0}])
    manifest = run / "predictions" / "abc123" / "DONE"
    records = json.loads(manifest.read_text())
    records[0]["parent_id"], records[1]["parent_id"] = "p", "q"
    manifest.write_text(json.dumps(records))
    (run / "designs_meta.csv").write_text("design,mutations,group\np,VH:K30D+VL:S52D,stability\nq,VH:A1T,affinity\n")
    df = load_designs([run]).set_index("id")
    assert df.loc["p_model_0", "meta_mutations"] == "VH:K30D+VL:S52D"
    assert df.loc["q_model_0", "meta_group"] == "affinity"


def test_designs_meta_sidecar_without_design_column_is_ignored(tmp_path):
    run = write_run(tmp_path, "r1", [{"id": "wt", "iptm": 0.9, "plddt": 90.0}])
    (run / "designs_meta.csv").write_text("name,group\nwt,x\n")
    assert "meta_group" not in load_designs([run]).columns


def test_load_designs_ignores_broken_manifests(tmp_path):
    run = write_run(tmp_path, "r1", ROWS)
    (run / "predictions" / "abc123" / "DONE").write_text("{not json")
    assert load_designs([run]).empty


def test_evaluate_pass_and_warn_counts(tmp_path):
    config = load_config(DEFAULT_THRESHOLDS)
    df = evaluate(load_designs([write_run(tmp_path, "r1", ROWS)]), config).set_index("id")
    assert bool(df.loc["wt", "pass"]) and df.loc["wt", "n_warn"] == 0
    assert not bool(df.loc["v3", "pass"])  # ipTM 0.40 is below the 0.5 gate
    assert not bool(df.loc["v4", "pass"])  # pLDDT 70 is below the 80 gate
    assert bool(df.loc["v2", "pass"]) and df.loc["v2", "n_warn"] == 2


def test_find_parent_by_substring(tmp_path):
    df = load_designs([write_run(tmp_path, "r1", ROWS)])
    assert find_parent(df, "wt")["id"] == "wt"
    assert find_parent(df, "nope") is None and find_parent(df, None) is None


def test_parent_with_several_seeds_is_averaged(tmp_path):
    rows = [
        {"id": "wt_model_a", "iptm": 0.90, "plddt": 90.0},
        {"id": "wt_model_b", "iptm": 0.80, "plddt": 88.0},
        {"id": "v1", "iptm": 0.85, "plddt": 89.0},
    ]
    df = load_designs([write_run(tmp_path, "r1", rows)])
    parent = find_parent(df, "wt_model")
    assert parent["iptm"] == pytest.approx(0.85) and parent["plddt"] == pytest.approx(89.0)
    assert parent["id"] == "wt_model_a" and parent["matched_ids"] == "wt_model_a|wt_model_b"
    out = tmp_path / "out"
    assert main(["build", str(tmp_path / "r1"), "--run-id", "r1", "--out", str(out), "--parent", "wt_model"]) == 0
    text = (out / "report.md").read_text()
    assert "parent 0.85" in text and text.count("(parent)") == 2


def test_default_thresholds_file_parses():
    config = load_config(DEFAULT_THRESHOLDS)
    assert config.primary.key == "iptm" and config.by_key("plddt").gate == 80


def test_build_writes_every_file(tmp_path):
    run = write_run(tmp_path, "r1", ROWS)
    out = tmp_path / "reports" / "r1"
    code = main(["build", str(run), "--run-id", "r1", "--out", str(out), "--parent", "wt"])
    assert code == 0
    for name in ("report.md", "designs.csv", "summary.json", "scatter.png", "compare.png"):
        assert (out / name).is_file(), name
    summary = json.loads((out / "summary.json").read_text())
    assert summary["n_designs"] == 5 and summary["n_pass"] == 3
    assert summary["metrics"]["iptm"]["parent"] == 0.90
    text = (out / "report.md").read_text()
    assert "3/5 designs pass all gates" in text and "parent 0.90" in text
    assert "VH:K30D+VL:S52D" in text and "stability" in text
    assert "nan" not in text.split("## Scorecard")[1]  # blank group or mutations print as "-"
    assert "wt (parent)" in text and "Δ ipTM" in text and "+0.00" in text
    assert pd.read_csv(out / "designs.csv").shape[0] == 5


def test_build_fails_clearly_on_empty_or_unknown_scores(tmp_path, capsys):
    assert main(["build", str(tmp_path), "--run-id", "x", "--out", str(tmp_path / "o")]) == 2
    run = write_run(tmp_path, "r2", [{"id": "a", "weird_score": 1.0}])
    assert main(["build", str(run), "--run-id", "x", "--out", str(tmp_path / "o")]) == 2
    assert "weird_score" in capsys.readouterr().err


def test_second_run_is_compared_with_the_first(tmp_path):
    reports = tmp_path / "reports"
    r1 = write_run(tmp_path, "r1", ROWS)
    main(["build", str(r1), "--run-id", "r1", "--out", str(reports / "r1"), "--history", str(reports)])
    better = [{**r, "iptm": min(1.0, r["iptm"] + 0.1)} for r in ROWS]
    r2 = write_run(tmp_path, "r2", better, producer="chai")
    main(["build", str(r2), "--run-id", "r2", "--out", str(reports / "r2"), "--history", str(reports)])
    text = (reports / "r2" / "report.md").read_text()
    assert "earlier best run r1" in text and "r1" in text and "(this run)" in text
    assert [s["run_id"] for s in load_summaries(reports)] == ["r1", "r2"]
    assert [s["run_id"] for s in load_summaries(reports, exclude="r2")] == ["r1"]


def test_rebuilding_a_run_does_not_compare_it_with_itself(tmp_path):
    reports = tmp_path / "reports"
    r1 = write_run(tmp_path, "r1", ROWS)
    for _ in range(2):
        main(["build", str(r1), "--run-id", "r1", "--out", str(reports / "r1"), "--history", str(reports)])
    assert "earlier best run" not in (reports / "r1" / "report.md").read_text()


def test_leaderboard_bolds_best_and_lists_top_designs(tmp_path):
    reports = tmp_path / "reports"
    r1 = write_run(tmp_path, "r1", ROWS)
    r2 = write_run(tmp_path, "r2", [{**r, "iptm": min(1.0, r["iptm"] + 0.05)} for r in ROWS], producer="chai")
    for rid, run in (("r1", r1), ("r2", r2)):
        main(["build", str(run), "--run-id", rid, "--out", str(reports / rid), "--history", str(reports)])
    out = tmp_path / "LEADERBOARD.md"
    assert main(["leaderboard", str(reports), "--out", str(out)]) == 0
    text = out.read_text()
    assert "[r1](reports/r1/report.md)" in text and "[r2](reports/r2/report.md)" in text
    assert "**0.87**" in text and "**0.82**" not in text  # r2's ipTM median beats r1's
    assert "Best passing designs across all runs" in text


def test_leaderboard_with_no_reports(tmp_path):
    assert "No reports yet" in build_leaderboard(tmp_path)
