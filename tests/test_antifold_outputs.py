"""The headline AntiFold numbers for trastuzumab, recomputed from the result files under tests/outputs/.

Nothing is run here: the scores were produced on the cluster (tests/outputs/README.md says how). The tests read the
committed files and check the numbers quoted in docs/track_a_filtering.md:

* AbBiBench trastuzumab variants: Spearman of `score` against the measured `binding_score`, overall and by mutation
  count, and the baseline of minus the mutation count with no model.
* SKEMPI 1N8Z rows scored with `antifold_seq`: Spearman of the score against minus the measured ddG.

Rank statistics use pandas ranks and numpy, so the suite needs no scipy.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

OUT = Path(__file__).parent / "outputs"
TOL = 0.005
"""Half a unit of the coarser rounding the numbers were quoted at."""


def spearman(x, y) -> float:
    """Pearson correlation of the average ranks, which is what scipy.stats.spearmanr returns."""
    return float(np.corrcoef(pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy())[0, 1])


def abbibench(tag: str) -> pd.DataFrame:
    """One row per scored variant: complex-route score, measured binding_score, number of mutations."""
    meta = pd.read_csv(OUT / "abbibench_trastuzumab" / f"meta_{tag}.csv")
    scores = pd.read_csv(OUT / "abbibench_trastuzumab" / f"scores_{tag}_complex.csv")
    frame = meta.merge(scores, on="variant_id", how="inner")
    frame["n_mut"] = frame["n_heavy"] + frame["n_light"]
    return frame


# file tag: (scored variants, rho overall, rho of minus the mutation count, [(low, high, n, rho) by mutation count])
ABB_EXPECTED = {
    "n8z": (419, 0.194, 0.257, [(1, 3, 23, 0.695), (4, 6, 86, -0.094), (7, 13, 310, 0.152)]),
    "d5": (2079, 0.314, 0.232, [(4, 6, 645, 0.626), (7, 9, 1429, 0.075)]),
}


@pytest.mark.parametrize("tag", ABB_EXPECTED)
def test_abbibench_spearman_overall_and_baseline(tag):
    n, rho, baseline, _ = ABB_EXPECTED[tag]
    frame = abbibench(tag)
    assert len(frame) == n
    assert spearman(frame["score"], frame["binding_score"]) == pytest.approx(rho, abs=TOL)
    assert spearman(-frame["n_mut"], frame["binding_score"]) == pytest.approx(baseline, abs=TOL)


@pytest.mark.parametrize("tag", ABB_EXPECTED)
def test_abbibench_spearman_by_mutation_count(tag):
    frame = abbibench(tag)
    for low, high, n, rho in ABB_EXPECTED[tag][3]:
        part = frame[frame["n_mut"].between(low, high)]
        assert len(part) == n, f"{low} to {high} mutations"
        assert spearman(part["score"], part["binding_score"]) == pytest.approx(rho, abs=TOL), f"{low} to {high} mutations"


def test_skempi_antifold_seq_spearman_with_minus_ddg():
    frame = pd.read_csv(OUT / "skempi_1n8z" / "skempi_1n8z_antifold_seq.csv")
    assert len(frame) == 8
    assert spearman(frame["score_antifold_seq"], -frame["ddg_exp_kcal_mol"]) == pytest.approx(0.60, abs=TOL)
    single = frame[~frame["label"].str.contains(r"\+|mutations")]
    assert len(single) == 6
    assert spearman(single["score_antifold_seq"], -single["ddg_exp_kcal_mol"]) == pytest.approx(0.71, abs=TOL)
