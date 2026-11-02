"""examples/track_a_run.py: the evidence-led funnel. The stages that need no weights: evidence modules, exclusions, assembly,
the null and the veto rule, on small synthetic evidence tables."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "track_a_run.py"
_spec = importlib.util.spec_from_file_location("track_a_run", SCRIPT)
run = importlib.util.module_from_spec(_spec)
sys.modules["track_a_run"] = run          # a dataclass looks its module up while it is being defined
_spec.loader.exec_module(run)

VH = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"


def row(study, variant, subs, metric, value, outcome, parent_value="", **extra):
    return {"study": study, "variant": variant, "substitutions_mature_domain_1based": subs, "metric": metric, "value": str(value),
            "parent_value": str(parent_value), "outcome": outcome, **extra}


EVIDENCE = [
    row("adapt2017", "H:D31K", "H:D31K", "affinity_fold_improvement", 3.2, "improved"),
    row("adapt2017", "L:D28K", "L:D28K", "affinity_fold_improvement", 2.4, "improved"),
    row("adapt2017", "H:D102F", "H:D102F", "affinity_fold_improvement", 4.2, "improved"),
    row("adapt2017", "L:S50F", "L:S50F", "affinity_fold_improvement", 0.3, "worsened"),
    row("moon2016", "AH06", "H:D102W;H:M107F;H:D108A;H:Y109L", "KD", 0.06, "improved affinity", parent_value=0.48),
    row("moon2016", "AL07", "L:T93I;L:P95N;L:P96I", "KD", 0.79, "weaker affinity", parent_value=0.48),
    row("moon2016", "A058", "H:D102W;H:M107F;H:D108A;H:Y109L;L:T93Q;L:P96A;L:T97S", "KD", 1.05, "worsened", parent_value=0.06),
    row("ng2025", "VL-Q89H", "L:Q89H", "cell_binding", "maintained", "unchanged"),
    row("ng2025", "VL-Q89H", "L:Q89H", "ADCC", "maintained", "unchanged"),
    row("ng2025", "VH-R50S/R59N", "H:R50S;H:R59N", "cell_binding", "reduced", "worsened"),
]


def test_substitutions_use_the_crystal_chain_ids():
    assert run.parse_substitutions("H:D31K;L:D28K") == (("A", 28, "D", "K"), ("B", 31, "D", "K"))
    with pytest.raises(ValueError):
        run.parse_substitutions("H:D31")


def test_modules_carry_the_measured_fold_and_the_loss_list_keeps_only_what_never_helped():
    modules, losses = run.load_modules(EVIDENCE)
    by_name = {m.name: m for m in modules}
    assert by_name["ADAPT H:D31K"].fold == 3.2 and by_name["ADAPT H:D31K"].site == "H1"
    assert by_name["ADAPT L:D28K"].site == "L1"
    assert by_name["ADAPT H:D102F"].site == "H3"
    assert by_name["Moon AH06"].fold == pytest.approx(8.0)
    assert by_name["Moon AH06"].points == pytest.approx(math.log10(8.0))
    assert "Ng VL-Q89H" in by_name and by_name["Ng VL-Q89H"].site == "FW"
    assert not any("S50F" in name for name in by_name), "a variant that lost binding is not a module"
    assert ("A", 50, "S", "F") in losses
    assert ("A", 93, "T", "Q") in losses and ("A", 96, "P", "A") in losses     # only seen in the worsened variant
    assert ("B", 102, "D", "W") not in losses, "D102W also appears in a gain, so it is not excluded"
    assert ("B", 50, "R", "S") in losses
    assert not any("R50S" in name for name in by_name)


def test_assembly_respects_size_hotspots_and_losses():
    modules, losses = run.load_modules(EVIDENCE)
    sets = run.assemble(modules, losses, min_size=5, max_size=10)
    assert sets, "the synthetic evidence should give at least one set of 5"
    for s in sets:
        assert 5 <= len(s["subs"]) <= 10
        positions = [(c, p) for c, p, _w, _m in s["subs"]]
        assert len(set(positions)) == len(positions)
        assert not any(pos in run.HOTSPOTS for pos in positions)
        assert not any(sub in losses for sub in s["subs"])
        assert not any(c == "A" and p in run.LIGHT_CDR3_WINDOW for c, p, _w, _m in s["subs"])
    best = max(sets, key=lambda s: s["priority"])
    assert "Moon AH06" in best["modules"], "the largest measured gain should lead the priority"


def test_assembly_never_returns_sets_below_the_minimum():
    modules, losses = run.load_modules(EVIDENCE)
    assert all(len(s["subs"]) >= 5 for s in run.assemble(modules, losses))
    assert run.assemble(modules, losses, min_size=11, max_size=12) == []


def test_liabilities_flag_new_cysteines_sequons_and_ng_dg():
    parent = VH
    assert not any(run.liabilities(parent, parent).values())
    new_cys = parent[:10] + "C" + parent[11:]
    assert run.liabilities(parent, new_cys)["new_cysteine"]
    i = parent.index("TISADTSKNTAY")                      # make N-X-T: the second residue after an inserted N
    sequon = parent[: i + 1] + "N" + "A" + "T" + parent[i + 4:]
    assert run.liabilities(parent, sequon)["new_nglyc_sequon"]
    assert run.liabilities("AAAAA", "AANGA")["new_NG"]
    assert run.liabilities("AAAAA", "AADGA")["new_DG"]


def test_random_sets_are_valid_substitutions_over_the_variable_domains():
    sets = run.random_sets(VH, VL, sizes=[5, 7], per_size=20, seed=1)
    assert len(sets) == 40 and sets == run.random_sets(VH, VL, sizes=[5, 7], per_size=20, seed=1)
    for subs in sets:
        assert len({(c, p) for c, p, _w, _m in subs}) == len(subs)
        for chain, position, wt, mut in subs:
            domain = VH if chain == "B" else VL
            assert domain[position - 1] == wt and mut != wt and 1 <= position <= len(domain)
    assert {len(s) for s in sets} == {5, 7}


def test_percentile_and_veto_rule():
    null = list(range(100))
    assert run.percentile_in_null(-1, null) == 0.0
    assert run.percentile_in_null(49, null) == pytest.approx(0.5)
    assert run.percentile_in_null(1000, null) == 1.0
    assert run.vetoed([0.01, 0.02], 0.05)
    assert not run.vetoed([0.01, 0.50], 0.05), "one model alone does not veto"
    assert not run.vetoed([0.05, 0.01], 0.05), "the cutoff itself is not inside the worst fraction"


def test_variant_id_matches_the_pipelines_canonical_form():
    assert run.variant_id([("B", 102, "D", "F"), ("A", 28, "D", "K")]) == "A:D28K+B:D102F"
