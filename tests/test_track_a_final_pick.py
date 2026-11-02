"""examples/track_a_final_pick.py: the pick rule on synthetic rows, and the files it writes."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "track_a_final_pick.py"
_spec = importlib.util.spec_from_file_location("track_a_final_pick", SCRIPT)
final = importlib.util.module_from_spec(_spec)
sys.modules["track_a_final_pick"] = final
_spec.loader.exec_module(final)

FIELDS = ["name", "seed", "batch", "eligible", "contact", "coulomb", "epitope_overlap", "n_added_cation", "paratope_net_charge",
          "iptm", "plddt", "ipsae_max", "ipsae_min"]


def row(name, seed, ipsae, contact=0.0, coulomb=0.0, batch="b0"):
    return {"name": name, "seed": str(seed), "batch": batch, "eligible": "1.0", "contact": str(contact), "coulomb": str(coulomb),
            "paratope_net_charge": "0.0", "iptm": "0.8", "ipsae_max": str(ipsae), "ipsae_min": str(ipsae - 0.1)}


PARENT = [row("parent", s, v) for s, v in enumerate((0.80, 0.70, 0.75))]          # mean 0.75, SD 0.05, floor 0.65


def candidate(name, ipsae, contact=0.0, coulomb=0.0):
    return [row(name, s, ipsae, contact, coulomb) for s in (0, 1)]


def test_the_top_by_score_come_first_and_the_chemistry_picks_follow_from_inside_the_band():
    rows = PARENT + candidate("hi", 0.90) + candidate("mid", 0.80) + candidate("chem", 0.70, 5, 5) + candidate("flat", 0.70) + candidate("collapsed", 0.20, 9, 9)
    variants = {"hi": "A:D28W+A:Y49W+A:Y92W", "mid": "B:Y52R+B:D102W+B:Y105W", "chem": "A:S50R+B:Y57W+B:R59F", "flat": "A:T94S+B:Y33W+B:D31R",
                "collapsed": "A:F53Y+B:G101F+B:M107R"}
    picks, band = final.pick(rows, variants, chem_weight=0.15, top=2, chem_picks=2)
    assert [p["name"] for p in picks] == ["hi", "mid", "chem", "flat"]
    assert [p["group"] for p in picks] == [final.GROUP_SCORE] * 2 + [final.GROUP_CHEM] * 2
    assert band["floor"] == pytest.approx(0.65)                                         # "collapsed" is below it, whatever its chemistry


def test_a_chemistry_pick_must_differ_from_every_earlier_pick_by_two_substitutions():
    rows = PARENT + candidate("top", 0.90) + candidate("twin", 0.70, 5, 5) + candidate("other", 0.70, 1, 1)
    variants = {"top": "A:D28W+A:Y49W+A:Y92W", "twin": "A:D28W+A:Y49W+B:Y52R", "other": "A:S50R+B:Y57W+B:R59F"}
    picks, _ = final.pick(rows, variants, top=1, chem_picks=1)
    assert [p["name"] for p in picks] == ["top", "other"]                               # "twin" shares two of its three with "top"


def test_every_parent_fold_counts_but_a_candidate_fold_seen_twice_counts_once(tmp_path):
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    for path, rows in ((a, PARENT[:2] + [row("x", 0, 0.5)]), (b, PARENT[:2] + [row("x", 0, 0.9), row("x", 1, 0.5)])):
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows({k: r.get(k, "") for k in FIELDS} for r in rows)
    rows = final.read_rows([str(tmp_path / "*.csv")])
    assert sum(r["name"] == "parent" for r in rows) == 4
    assert sorted((r["name"], r["seed"], r["ipsae_max"]) for r in rows if r["name"] == "x") == [("x", "0", "0.5"), ("x", "1", "0.5")]


def test_main_writes_the_picks_and_the_parents_band(tmp_path):
    rows = PARENT + candidate("a", 0.90) + candidate("b", 0.80)
    rows_csv, cands, out = tmp_path / "rows.csv", tmp_path / "candidates.csv", tmp_path / "final.csv"
    with open(rows_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows({k: r.get(k, "") for k in FIELDS} for r in rows)
    cands.write_text("name,variant_id\na,A:D28W\nb,B:Y52R\n")
    final.main(["--rows", str(rows_csv), "--candidates", str(cands), "--out", str(out), "--top", "2", "--chem-picks", "0"])
    with open(out, newline="") as handle:
        written = list(csv.DictReader(handle))
    assert [(r["design"], r["name"], r["variant_id"]) for r in written] == [("design_01", "a", "A:D28W"), ("design_02", "b", "B:Y52R")]
    assert json.loads(out.with_suffix(".parent.json").read_text())["parent_folds"] == 3
