"""AntiFold's per-residue table for the trastuzumab parent 1N8Z, as committed under tests/outputs/antifold_tables/.

The table is AntiFold's own output (tests/outputs/README.md says how it was made), so these tests pin what is true of
the file and what the rest of tests/outputs relies on: its shape, that each row is a probability distribution over
the 20 amino acids, and that its row order is the position numbering of the AbBiBench files.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

OUT = Path(__file__).parent / "outputs"
AMINO = list("ACDEFGHIKLMNPQRSTVWY")
HEAVY_START = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRY"
LIGHT_START = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPS"


@pytest.fixture(scope="module")
def table() -> pd.DataFrame:
    return pd.read_csv(OUT / "antifold_tables" / "1N8Z_BA.csv", dtype={"pdb_posins": str, "assumed_region": str}, keep_default_na=False)


def test_columns_and_shape(table):
    assert list(table.columns) == ["pdb_posins", "pdb_chain", "pdb_res", "top_res", "pdb_pos", "perplexity", "assumed_region", *AMINO]
    assert table["pdb_chain"].value_counts().to_dict() == {"B": 220, "A": 214}
    assert not table.duplicated(["pdb_chain", "pdb_posins"]).any()
    assert np.isfinite(table[AMINO].to_numpy()).all()


def test_chains_are_the_trastuzumab_heavy_and_light_chains(table):
    heavy = "".join(table.loc[table["pdb_chain"] == "B", "pdb_res"])
    light = "".join(table.loc[table["pdb_chain"] == "A", "pdb_res"])
    assert heavy.startswith(HEAVY_START)
    assert light.startswith(LIGHT_START)


def test_each_row_is_a_distribution_over_the_20_amino_acids(table):
    log_p = table[AMINO].to_numpy()
    assert np.abs(np.log(np.exp(log_p).sum(axis=1))).max() < 1e-3
    assert (log_p <= 0).all()
    assert (np.array(AMINO)[log_p.argmax(axis=1)] == table["top_res"].to_numpy()).all()


def test_perplexity_is_exp_of_the_entropy(table):
    p = np.exp(table[AMINO].to_numpy())
    entropy = -(p * np.log(p)).sum(axis=1)
    assert np.exp(entropy) == pytest.approx(table["perplexity"].to_numpy(), abs=1e-3)


def test_pdb_pos_drifts_from_the_row_number_so_position_means_row(table):
    """Why the README warns: positions in the other files are rows within a chain, not the numbering in `pdb_pos`
    (IMGT, as AntiFold renumbers: the heavy variable domain ends at 128 on row 120, the light at 127 on row 107)."""
    for chain, first_row_that_differs, variable_domain_rows, last_variable_number in (("B", 10, 120, 128), ("A", 30, 107, 127)):
        rows = table[table["pdb_chain"] == chain].reset_index(drop=True)
        differs = np.flatnonzero(rows["pdb_pos"].to_numpy() != np.arange(1, len(rows) + 1))
        assert differs[0] + 1 == first_row_that_differs
        assert rows.loc[variable_domain_rows - 1, "pdb_pos"] == last_variable_number


def test_rows_line_up_with_the_single_mutations_scored_through_the_table_route(table):
    """Each id in scores_singles_table.csv is <chain><position><wt>><mut>; the wild type must be row <position> of that chain."""
    singles = pd.read_csv(OUT / "abbibench_trastuzumab" / "scores_singles_table.csv")
    wild_type = {chain: rows["pdb_res"].tolist() for chain, rows in table.groupby("pdb_chain")}
    assert len(singles) == 236
    for variant_id in singles["variant_id"]:
        chain, position, wt, _mut = re.fullmatch(r"([AB])(\d+)([A-Z])>([A-Z])", variant_id).groups()
        assert wild_type[chain][int(position) - 1] == wt, variant_id
