"""The verdict table, with one state per test instead of one verdict per scorer.

Re-labels what Part 5 of `notebooks/03_negative_controls.ipynb` already measured. It
scores nothing: every AUROC and every interval below is read out of that notebook's
output CSVs under `runs/negative_controls/`, and the one value computed here is the
ipSAE direction AUROC, from the panel file the notebook writes but never reads back.

The bug it fixes. Part 5 asked one question per test - "is the lower end of the CI
above 0.5?" - and called everything else a failure. That collapses two different
results into one word: a scorer whose interval sits entirely below 0.5 has been shown
to be wrong, and a scorer whose interval spans 0.5 has been shown nothing. Five
scorers were labelled FAIL on the direction panel while their intervals ran to 0.967,
1.000 and 0.600 - a panel of 13 cannot resolve them, and the old table reported that
inability as evidence.

So each test gets one of four states:

    PASS       the 95% interval lies entirely above 0.5
    FAIL       it lies entirely below 0.5
    UNPROVEN   it spans 0.5 - no claim either way
    N/T        not testable: the scorer has no route, or the panel has one label

and the floor checks keep their rates and get no state at all, because a share
against a fixed threshold is not an interval against 0.5.

Writes `runs/negative_controls/verdict_states.csv` beside the notebook's own
`verdict.csv`, which is left alone: that file is the record of what that run printed.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

# The notebook's Part 2 settings, so the one AUROC computed here is computed its way.
BOOTSTRAP_SEED = 20261103
BOOTSTRAP = 10_000
CI = 0.95

PLAUSIBLE_BENEFICIAL = 0.25     # diagnostics.PLAUSIBLE_BENEFICIAL
EXPECTED_WT_ARGMAX = 0.60       # diagnostics.EXPECTED_WT_ARGMAX

ORDER = ["foldx", "antifold_seq", "antifold", "antiberty", "esm1v", "esm2_650m",
         "ablang2", "blosum62", "ipsae_max", "consensus", "gearbind"]

TESTS = ("within_position", "direction", "kinetics")

# Why a row has no state, where "not testable" needs qualifying. N/T is one word for
# three different situations and the difference matters to anyone deciding what to
# re-run, so every N/T in the table names its own.
WHY = {
    "blosum62": "within position: interval computed in the 2026-11-22 run and printed, "
                "not written to conservativeness.csv - recovered here from the executed "
                "notebook's own output. Re-running Part 4 now persists it. Not tested "
                "on direction or kinetics: it ran only as the conservativeness reference",
    "ipsae_max": "within position and kinetics: never scored on either panel. Direction "
                 "is the only test it has, and it is UNPROVEN. Its exclusion from the "
                 "ranking rests on ICC 0.298 over 3 seeds (116 designs), not on this row",
    "consensus": "no route on any box: needs an OAS frequency table that does not exist "
                 "in this repository. Never run, so nothing is shown either way",
    "gearbind": "no route on this box: container image not provisioned. Never run",
}

# The matrix scorers' kinetics N/T has one shared cause, so it is stated once.
ONE_LABEL = ("kinetics: the Chinery panel returns 14 fitted wells and 0 unfittable for "
             "the matrix scorers - one label, so there is no AUROC to have")

# blosum62's within-position interval, as printed by Part 4 of the 2026-11-22 run:
#   "BLOSUM62 AS A SCORER on the position-matched test: 0.536 (0.433 to 0.639) over
#    126 pairs"
# - `runs/negative_controls/03_negative_controls.executed.ipynb`, cell 30 output.
# It was computed by `stratified_ci` and printed, and the write to
# conservativeness.csv dropped it, so the file carries the point estimate alone. The
# notebook now persists it; until that run happens this is the recorded value, and it
# is used ONLY when the file has no interval of its own. It is a recovered measurement,
# not an assumption: blosum62 is UNPROVEN on this test, never untested.
BLOSUM62_RECORDED_CI = (0.433, 0.639)


def state(ci_low: float, ci_high: float) -> str:
    """One of PASS, FAIL, UNPROVEN, N/T for a 95% interval against 0.5."""
    if ci_low != ci_low or ci_high != ci_high:
        return "N/T"
    if ci_low > 0.5:
        return "PASS"
    if ci_high < 0.5:
        return "FAIL"
    return "UNPROVEN"


def auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    """AUROC from average ranks, so ties count half. The notebook's definition."""
    n_pos, n_neg = int((labels == 1).sum()), int((labels == 0).sum())
    ranks = rankdata(scores)
    u = float(ranks[labels == 1].sum()) - n_pos * (n_pos + 1) / 2.0
    return u / float(n_pos * n_neg)


def auroc_ci(scores: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    """Percentile CI from a bootstrap stratified on the label, the notebook's again."""
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    pos, neg = np.flatnonzero(labels == 1), np.flatnonzero(labels == 0)
    draws = np.empty(BOOTSTRAP, dtype=float)
    for i in range(BOOTSTRAP):
        take = np.concatenate([rng.choice(pos, len(pos), True),
                               rng.choice(neg, len(neg), True)])
        draws[i] = auroc(scores[take], labels[take])
    tail = (1.0 - CI) / 2.0
    return float(np.quantile(draws, tail)), float(np.quantile(draws, 1.0 - tail))


def ipsae_direction(path: Path) -> dict:
    """The ipSAE row of the direction test, from the panel `ipsae_panel.py` writes.

    The notebook never read this file back, so the AUROC quoted for ipSAE had no
    source. It does now, and it is UNPROVEN: 0.476 on 21 against 8, interval 0.23 to
    0.73. ipSAE's exclusion from the ranking rests on its ICC of 0.298 over three
    seeds, which is a precision argument and is not in this table.
    """
    panel = pd.read_csv(path)
    scores = panel["score_ipsae"].to_numpy(dtype=float)
    labels = panel["label"].to_numpy(dtype=int)
    low, high = auroc_ci(scores, labels)
    return {"auroc": auroc(scores, labels), "ci_low": low, "ci_high": high,
            "n": len(panel), "n_pos": int((labels == 1).sum()),
            "n_neg": int((labels == 0).sum())}


def build(runs: Path) -> pd.DataFrame:
    matched = pd.read_csv(
        runs / "manufactured" / "position_matched_auroc.csv").set_index("model")
    conservativeness = pd.read_csv(
        runs / "manufactured" / "conservativeness.csv").set_index("model")
    direction = pd.read_csv(runs / "direction" / "auroc_overall.csv").set_index("model")
    kinetics = pd.read_csv(runs / "chinery" / "auroc_clean.csv").set_index("model")
    floors = pd.read_csv(runs / "verdict.csv").set_index("model")
    ipsae = ipsae_direction(runs / "ipsae" / "direction_ipsae.csv")

    rows = []
    for model in ORDER:
        row = {"model": model}

        # Test 1 - within position. The only test anything fails.
        if model in matched.index:
            row.update({"within_position_auroc": matched.loc[model, "auroc"],
                        "within_position_ci_low": matched.loc[model, "ci_low"],
                        "within_position_ci_high": matched.loc[model, "ci_high"],
                        "within_position_n": int(matched.loc[model, "pairs"])})
        elif model in conservativeness.index:
            # blosum62 ran as the conservativeness reference. Its interval exists; the
            # question is only whether this copy of the file carries it yet.
            low = conservativeness.loc[model].get("ci_low", float("nan"))
            high = conservativeness.loc[model].get("ci_high", float("nan"))
            if low != low or high != high:
                low, high = BLOSUM62_RECORDED_CI
            row.update({"within_position_auroc": conservativeness.loc[model, "matched_auroc"],
                        "within_position_ci_low": low,
                        "within_position_ci_high": high,
                        "within_position_n": int(matched["pairs"].iloc[0])})
        else:
            row.update({"within_position_auroc": float("nan"),
                        "within_position_ci_low": float("nan"),
                        "within_position_ci_high": float("nan"),
                        "within_position_n": 0})

        # Test 2 - direction. n is 13 for the matrix scorers and 29 for the rest.
        if model in direction.index:
            row.update({"direction_auroc": direction.loc[model, "auroc"],
                        "direction_ci_low": direction.loc[model, "ci_low"],
                        "direction_ci_high": direction.loc[model, "ci_high"],
                        "direction_n": int(direction.loc[model, "n_pos"]
                                           + direction.loc[model, "n_neg"]),
                        "direction_n_pos": int(direction.loc[model, "n_pos"]),
                        "direction_n_neg": int(direction.loc[model, "n_neg"])})
        elif model == "ipsae_max":
            row.update({"direction_auroc": ipsae["auroc"],
                        "direction_ci_low": ipsae["ci_low"],
                        "direction_ci_high": ipsae["ci_high"],
                        "direction_n": ipsae["n"],
                        "direction_n_pos": ipsae["n_pos"],
                        "direction_n_neg": ipsae["n_neg"]})
        else:
            row.update({"direction_auroc": float("nan"), "direction_ci_low": float("nan"),
                        "direction_ci_high": float("nan"), "direction_n": 0,
                        "direction_n_pos": 0, "direction_n_neg": 0})

        # Test 3 - clean kinetics. One label only for the matrix scorers: 14 fitted
        # wells, 0 unfittable, so there is no AUROC to have.
        if model in kinetics.index:
            row.update({"kinetics_auroc": kinetics.loc[model, "auroc"],
                        "kinetics_ci_low": kinetics.loc[model, "ci_low"],
                        "kinetics_ci_high": kinetics.loc[model, "ci_high"],
                        "kinetics_n": int(kinetics.loc[model, "n_pos"]
                                          + kinetics.loc[model, "n_neg"]),
                        "kinetics_n_pos": int(kinetics.loc[model, "n_pos"]),
                        "kinetics_n_neg": int(kinetics.loc[model, "n_neg"])})
        else:
            row.update({"kinetics_auroc": float("nan"), "kinetics_ci_low": float("nan"),
                        "kinetics_ci_high": float("nan"), "kinetics_n": 0,
                        "kinetics_n_pos": 0, "kinetics_n_neg": 0})

        for test in TESTS:
            row[f"{test}_state"] = state(row[f"{test}_ci_low"], row[f"{test}_ci_high"])

        # Floor checks. Rates against fixed thresholds, deliberately stateless.
        rate = floors.loc[model, "failure_rate"] if model in floors.index else float("nan")
        argmax = floors.loc[model, "wt_argmax"] if model in floors.index else float("nan")
        row["manufactured_failure_rate"] = rate
        row["manufactured_over_threshold"] = (
            "" if rate != rate else str(bool(rate > PLAUSIBLE_BENEFICIAL)))
        row["wt_argmax_rate"] = argmax
        row["wt_argmax_under_threshold"] = (
            "" if argmax != argmax else str(bool(argmax < EXPECTED_WT_ARGMAX)))

        why = [WHY[model]] if model in WHY else []
        if row["kinetics_state"] == "N/T" and row["kinetics_n"] == 14:
            why.append(ONE_LABEL)
        row["why"] = ". ".join(why)
        rows.append(row)

    return pd.DataFrame(rows).set_index("model")


def cell(row: pd.Series, test: str) -> str:
    value = row[f"{test}_auroc"]
    if value != value:
        return "-"
    low, high = row[f"{test}_ci_low"], row[f"{test}_ci_high"]
    if low != low:
        return f"{value:.3f} (no CI)"
    return f"{value:.3f} ({low:.2f}-{high:.2f})"


def report(table: pd.DataFrame) -> None:
    head = (f"{'scorer':<14}{'WITHIN POSITION':>24}{'n':>5}{'  state':<11}"
            f"{'DIRECTION':>24}{'n':>4}{'  state':<11}"
            f"{'CLEAN KINETICS':>24}{'n':>5}{'  state':<11}{'fail':>6}{'wt':>6}")
    print(head)
    print("-" * len(head))
    for model, row in table.iterrows():
        rate, argmax = row["manufactured_failure_rate"], row["wt_argmax_rate"]
        print(f"{model:<14}"
              f"{cell(row, 'within_position'):>24}{(row['within_position_n'] or '-'):>5}"
              f"  {row['within_position_state']:<9}"
              f"{cell(row, 'direction'):>24}{(row['direction_n'] or '-'):>4}"
              f"  {row['direction_state']:<9}"
              f"{cell(row, 'kinetics'):>24}{(row['kinetics_n'] or '-'):>5}"
              f"  {row['kinetics_state']:<9}"
              f"{('-' if rate != rate else f'{rate:.0%}'):>6}"
              f"{('-' if argmax != argmax else f'{argmax:.0%}'):>6}")
    print("-" * len(head))

    for test in TESTS:
        failed = table.index[table[f"{test}_state"] == "FAIL"].tolist()
        print(f"FAIL on {test:<16}{', '.join(failed) if failed else 'none'}")
    print()
    print("Every FAIL in the project is a within-position FAIL. That test holds the")
    print("position, the residue, the neighbours and the burial fixed and varies only")
    print("which amino acid goes in, over 126 pairs - enough for an interval that")
    print("clears 0.5 in one direction or the other.")
    print()
    matrix = table[table["direction_n"] == 13]
    print(f"The direction test cannot resolve the {len(matrix)} matrix scorers: n = 13 "
          f"({int(matrix['direction_n_pos'].iloc[0])} improvements, "
          f"{int(matrix['direction_n_neg'].iloc[0])} losses), because they")
    print("score singles only. Four of their intervals touch 0.000 and ablang2's runs")
    print("0.000 to 1.000. The smallest p this design can reach is 0.0035")
    print("(direction/auroc_overall.csv, p_floor), so a point estimate from it carries")
    print("no information: none of these five is shown to work or shown to be wrong.")
    print()
    over = table.index[table["manufactured_over_threshold"] == "True"].tolist()
    under = table.index[table["wt_argmax_under_threshold"] == "True"].tolist()
    print(f"The floor checks carry no state. A failure rate above "
          f"{PLAUSIBLE_BENEFICIAL:.0%} ({', '.join(over)})")
    print(f"or a wt_argmax below {EXPECTED_WT_ARGMAX:.0%} ({', '.join(under)}) says the "
          f"column is unusable,")
    print("not that the scorer was measured wrong: there is no interval and no 0.5 for")
    print("it to be above or below.")
    print()
    print("Why a row has no state, where N/T needs qualifying:")
    for model, row in table.iterrows():
        if row["why"]:
            print(f"  {model:<14}{row['why']}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--runs", type=Path, default=Path("runs/negative_controls"))
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)

    table = build(args.runs)
    report(table)
    out = args.out or args.runs / "verdict_states.csv"
    table.to_csv(out)
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
