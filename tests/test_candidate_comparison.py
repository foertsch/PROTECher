"""Maria's whole-candidate reporting: identity, missingness, aggregation and decisions."""
import copy
import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("candidate_comparison", ROOT / "filters/compare_candidates.py")
cmp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cmp)


def write_csv(path, rows, fields=None):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        w.writeheader()
        w.writerows(rows)


@pytest.fixture
def candidates():
    return {"PARENT": "ACDE", "a": "CCDE", "b": "ADDE", "c": "ACDF"}


@pytest.fixture
def config():
    return {"schema_version": 1, "parent_id": "PARENT", "require_developability": False,
            "sources": [{"name": "model", "path": "scores.csv", "sequence_column": "sequence",
                         "protocol": "test_protocol", "replicate_column": "seed", "expected_replicates": 2,
                         "metrics": [{"name": "score", "column": "score", "better": "higher", "required": True}]}],
            "gates": []}


def observation(name, value, metric="model.score", seed="0"):
    return {"design": name, "value": value, "metric": metric, "replicate": seed}


def test_median_counts_parent_deltas_and_no_parent_rank(candidates, config):
    obs = [observation("a", 0), observation("a", 10, seed="1"),
           observation("PARENT", 2), observation("PARENT", 4, seed="1")]
    frame, _ = cmp.compare(candidates, config, obs, {})
    frame = frame.set_index("design")
    assert frame.at["a", "model.score"] == 5
    assert frame.at["a", "model.score.delta_parent"] == 2
    assert frame.at["a", "model.score.n"] == 2
    assert frame.at["a", "decision"] == "ELIGIBLE"
    assert frame.at["b", "decision"] == "INCOMPLETE"
    assert frame.at["PARENT", "decision"] == "REFERENCE"
    assert cmp.pd.isna(frame.at["PARENT", "model.score.rank"])


def test_partial_replicates_do_not_pass_gate(candidates, config):
    config["gates"] = [{"metric": "model.score", "operator": ">=", "threshold": 0, "action": "reject"}]
    frame, _ = cmp.compare(candidates, config, [observation("a", 10)], {})
    row = frame.set_index("design").loc["a"]
    assert row["model.score.status"] == "partial"
    assert row["decision"] == "INCOMPLETE"
    assert cmp.pd.isna(row["model.score.rank"])


@pytest.mark.parametrize("action,expected", [("reject", "REJECT"), ("review", "REVIEW")])
def test_explicit_gate_action(candidates, config, action, expected):
    config["gates"] = [{"metric": "model.score", "operator": ">=", "threshold": 3, "action": action}]
    frame, _ = cmp.compare(candidates, config, [observation("a", 1), observation("a", 2, seed="1")], {})
    assert frame.set_index("design").at["a", "decision"] == expected


def test_lower_is_better_rank_and_descriptive_unranked(candidates, config):
    config["sources"][0]["expected_replicates"] = 1
    config["sources"][0]["metrics"][0]["better"] = "lower"
    obs = [observation("a", 1), observation("b", 2)]
    frame, _ = cmp.compare(candidates, config, obs, {})
    assert frame.set_index("design").at["a", "model.score.rank"] == 1
    config["sources"][0]["metrics"][0]["better"] = "descriptive"
    frame, _ = cmp.compare(candidates, config, obs, {})
    assert frame["model.score.rank"].isna().all()


@pytest.mark.parametrize("bad", ["nan", "inf", "-inf", "failed", "N/A"])
def test_bad_score_is_not_silently_zero_or_missing(bad):
    with pytest.raises(ValueError):
        cmp.number(bad, "test")
    assert cmp.number("", "test") is None
    assert cmp.number("0", "test") == 0


def test_unknown_stale_or_duplicate_predictions_fail(tmp_path, candidates, config):
    source = config["sources"][0]
    good = {"design": "a", "sequence": candidates["a"], "seed": "0", "score": "1"}
    for rows in ([{**good, "design": "unknown"}], [{**good, "sequence": "AAAA"}], [good, good]):
        write_csv(tmp_path / "scores.csv", rows)
        with pytest.raises(ValueError):
            cmp.load_source(source, candidates, tmp_path)


def test_score_hash_validation(tmp_path, candidates, config):
    source = config["sources"][0]
    source.pop("sequence_column")
    source["hash_column"] = "sequence_sha256"
    row = {"design": "a", "sequence_sha256": cmp.digest(candidates["a"]), "seed": "0", "score": "1"}
    write_csv(tmp_path / "scores.csv", [row])
    assert len(cmp.load_source(source, candidates, tmp_path)[0]) == 1
    row["sequence_sha256"] = "0" * 64
    write_csv(tmp_path / "scores.csv", [row])
    with pytest.raises(ValueError, match="identity mismatch"):
        cmp.load_source(source, candidates, tmp_path)


def test_id_only_results_require_frozen_run_mapping(tmp_path, candidates, config):
    source = config["sources"][0]
    source.pop("sequence_column")
    source["identity_file"] = "identities.csv"
    source["id_column"] = "variant_id"
    write_csv(tmp_path / "identities.csv", [{"source_id": "B:A1C", "design": "a", "sequence": candidates["a"]}])
    write_csv(tmp_path / "scores.csv", [{"variant_id": "B:A1C", "seed": "1", "score": "2"}])
    obs, _ = cmp.load_source(source, candidates, tmp_path)
    assert obs[0]["design"] == "a"
    write_csv(tmp_path / "identities.csv", [{"source_id": "B:A1C", "design": "a", "sequence": "AAAA"}])
    with pytest.raises(ValueError):
        cmp.load_source(source, candidates, tmp_path)


def test_mixed_protocols_are_not_averaged(tmp_path, candidates, config):
    source = config["sources"][0]
    source["protocol_column"] = "protocol"
    write_csv(tmp_path / "scores.csv", [{"design": "a", "sequence": candidates["a"], "seed": "0", "score": "1", "protocol": "different_target"}])
    with pytest.raises(ValueError, match="mixed protocols"):
        cmp.load_source(source, candidates, tmp_path)


def test_missing_source_preserves_entire_manifest(tmp_path, candidates, config):
    obs, source = cmp.load_source(config["sources"][0], candidates, tmp_path)
    assert source["status"] == "file_missing"
    frame, _ = cmp.compare(candidates, config, obs, {})
    assert len(frame) == len(candidates)
    assert frame["model.score"].isna().all()
    assert set(frame[~frame["is_reference"]]["decision"]) == {"INCOMPLETE"}


def test_developability_warnings_are_not_hard_failures(candidates, config):
    config["sources"] = []
    config["require_developability"] = True
    evidence = {"input_status": "valid", "rules_status": "pass", "flags": "new_liability", "data_warnings": ""}
    dev = {"a": evidence, "b": {**evidence, "rules_status": "fail", "rule_issues": "linker_changed"}}
    frame, _ = cmp.compare(candidates, config, [], dev)
    f = frame.set_index("design")
    assert f.at["a", "decision"] == "REVIEW"
    assert f.at["b", "decision"] == "REJECT"
    assert f.at["c", "decision"] == "INCOMPLETE"


def test_developability_identity_and_reference_alias(tmp_path, candidates):
    candidates["wt"] = candidates.pop("PARENT")
    row = {"design": "PARENT", "sequence": candidates["wt"], "is_reference": "True",
           "input_status": "valid", "rules_status": "reference", "flags": "", "data_warnings": ""}
    path = tmp_path / "dev.csv"
    write_csv(path, [row])
    assert "wt" in cmp.load_developability(path, candidates, "wt")
    write_csv(path, [{**row, "design": "a", "sequence": "AAAA", "is_reference": "False"}])
    with pytest.raises(ValueError, match="identity mismatch"):
        cmp.load_developability(path, candidates, "wt")


def test_duplicate_sequences_retained_and_shortlist_deduplicates(candidates, config):
    candidates["copy"] = candidates["a"]
    config["sources"][0]["expected_replicates"] = 1
    config["shortlist"] = {"size": 5, "sort_by": ["model.score"], "include_review": True}
    obs = [observation(name, value) for name, value in [("a", 5), ("copy", 5), ("b", 4), ("c", 3)]]
    frame, metrics = cmp.compare(candidates, config, obs, {})
    result, message = cmp.select_shortlist(frame, config, metrics)
    assert len(frame) == 5
    assert list(result["design"]) == ["a", "b", "c"]
    assert "3/5" in message
    assert not result["is_reference"].any()


def test_shortlist_waits_for_priority_decision(candidates, config):
    frame, metrics = cmp.compare(candidates, config, [], {})
    selected, note = cmp.select_shortlist(frame, config, metrics)
    assert selected.empty
    assert "team agreement" in note


def test_ipsae_is_not_primary_shortlist_metric(config):
    config["sources"][0]["metrics"][0]["name"] = "ipsae"
    config["shortlist"] = {"size": 15, "sort_by": ["model.ipsae"]}
    with pytest.raises(ValueError, match="ipSAE"):
        cmp.validate_config(config)


def test_agreement_has_overlap_counts_and_direction(candidates, config):
    source = config["sources"][0]
    source["expected_replicates"] = 1
    other = copy.deepcopy(source)
    other["name"] = "energy"
    other["metrics"][0]["better"] = "lower"
    config["sources"].append(other)
    obs = []
    for n, v in zip(("a", "b", "c"), (1, 2, 3)):
        obs += [observation(n, v), observation(n, -v, "energy.score")]
    frame, metrics = cmp.compare(candidates, config, obs, {})
    agreement = cmp.correlations(frame, metrics)
    assert agreement.iloc[0]["n_common"] == 3
    assert agreement.iloc[0]["spearman"] == pytest.approx(1)
    assert len(agreement) == 2  # overall and the one-mutation stratum


@pytest.mark.parametrize("contents", ["name,sequence\na,ACD\na,ACE\n", "name,sequence\na,ACX\n", "name,sequence\na,\n"])
def test_bad_manifest_is_rejected(tmp_path, contents):
    path = tmp_path / "input.csv"
    path.write_text(contents)
    with pytest.raises(ValueError):
        cmp.read_candidates(path)


def test_wrapped_fasta_keeps_exact_ids(tmp_path):
    path = tmp_path / "input.fasta"
    path.write_text(">candidate_01\nacd\nefg\n")
    assert cmp.read_candidates(path) == {"candidate_01": "ACDEFG"}


def test_demo_end_to_end(tmp_path):
    demo = ROOT / "filters/demo"
    frame, summary = cmp.run(demo / "candidates.csv", demo / "comparison.json", tmp_path, demo / "developability.csv")
    assert summary["counts"] == {"ELIGIBLE": 1, "REVIEW": 1, "INCOMPLETE": 1, "REJECT": 1}
    assert len(frame) == 5
    assert len(cmp.pd.read_csv(tmp_path / "shortlist.csv")) == 1
    assert len(cmp.pd.read_csv(tmp_path / "observations.csv")) == 23
    assert "not a ranking" in (tmp_path / "report.md").read_text()
    meta = json.loads((tmp_path / "summary.json").read_text())
    assert all(record["sha256"] for record in meta["inputs"])


def test_input_cannot_be_overwritten(tmp_path):
    demo = ROOT / "filters/demo"
    candidate = tmp_path / "combined.csv"
    candidate.write_text((demo / "candidates.csv").read_text())
    original = candidate.read_bytes()
    with pytest.raises(ValueError, match="overwrite"):
        cmp.run(candidate, demo / "comparison.json", tmp_path, demo / "developability.csv")
    assert candidate.read_bytes() == original


def test_template_runs_without_gpu_outputs(tmp_path):
    # The shipped config uses paths that may exist in a user's checkout, so move it.
    config_path = tmp_path / "config.json"
    config_path.write_text((ROOT / "filters/candidate_comparison.json").read_text())
    frame, _ = cmp.run(ROOT / "filters/demo/candidates.csv", config_path, tmp_path / "out")
    assert set(frame[~frame["is_reference"]]["decision"]) == {"INCOMPLETE"}
