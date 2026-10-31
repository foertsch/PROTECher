"""The wild-type round trip over an AntiFold table, on the real 1N8Z chains.

The table AntiFold wrote for notebook 03 lives on the cluster node, so what is
exercised here is the check itself: a faithful table built from
`data/structures/1N8Z.pdb` passes, and the three ways a table can be misaligned
each raise. Without this, "the check passed" on the node would only mean the
check cannot fail.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from antifold_score import check_wt
from design_pipeline.matrices import IndexMismatchError
from design_pipeline.structures import AMINO_ACIDS, chain_sequence

PDB = Path("data/structures/1N8Z.pdb")
HEAVY, LIGHT = "B", "A"

pytestmark = pytest.mark.skipif(not PDB.is_file(), reason=f"{PDB} is not in the checkout")


def rows_for(chain: str, sequence: str, offset: int = 0) -> pd.DataFrame:
    """An AntiFold-shaped table for one chain, peaked on the residue it declares.

    `offset` shifts the declared wild type against the structure, which is the
    failure the round trip exists to catch: every row is still a real
    distribution, just the wrong one.
    """
    generator = np.random.default_rng(len(sequence))
    declared = sequence[offset:] + sequence[:offset]
    values = -generator.uniform(4.0, 12.0, size=(len(sequence), len(AMINO_ACIDS)))
    for i, residue in enumerate(declared):
        values[i, AMINO_ACIDS.index(residue)] = -0.01     # the declared residue is the argmax
    frame = pd.DataFrame(values, columns=list(AMINO_ACIDS))
    frame.insert(0, "pdb_posins", [str(i + 1) for i in range(len(sequence))])
    frame.insert(1, "pdb_chain", chain)
    frame.insert(2, "pdb_res", list(declared))
    frame.insert(3, "top_res", [AMINO_ACIDS[i] for i in values.argmax(axis=1)])
    frame.insert(4, "pdb_pos", range(1, len(sequence) + 1))
    return frame


@pytest.fixture
def table(tmp_path):
    def write(frames, name="1N8Z_BA.csv"):
        path = tmp_path / name
        pd.concat(frames, ignore_index=True).to_csv(path, index=False)
        return path

    return write


@pytest.fixture
def sequences():
    return {chain: chain_sequence(PDB, chain) for chain in (HEAVY, LIGHT)}


def test_a_faithful_table_passes_and_says_so_per_chain(table, sequences, capsys):
    path = table([rows_for(chain, seq) for chain, seq in sequences.items()])

    assert check_wt.main(["--table", str(path), "--pdb", str(PDB),
                          "--heavy", HEAVY, "--light", LIGHT]) == 0

    printed = capsys.readouterr().out
    assert "every check passed" in printed
    for chain, seq in sequences.items():
        assert f"chain {chain}: {len(seq)} table rows, {len(seq)} residues" in printed
    # the shift profile is printed, and only offset 0 fits
    assert "+0:100%" in printed and "+1:100%" not in printed and "-1:100%" not in printed


def test_one_position_of_shift_raises_and_names_the_first_difference(table, sequences):
    """The whole point: a shifted table is 20 real distributions per row."""
    path = table([
        rows_for(HEAVY, sequences[HEAVY], offset=1),
        rows_for(LIGHT, sequences[LIGHT]),
    ])

    with pytest.raises(IndexMismatchError, match=r"chain B: matrix rows do not match"):
        check_wt.check(path, PDB, [HEAVY, LIGHT])


def test_a_label_the_chain_does_not_have_is_a_failure(table, sequences):
    frames = [rows_for(chain, seq) for chain, seq in sequences.items()]
    frames[0].loc[7, "pdb_posins"] = "7A"       # an insertion code the input numbering has not got

    with pytest.raises(check_wt.CheckFailed, match=r"label\(s\) AntiFold reports that the chain"):
        check_wt.check(table(frames), PDB, [HEAVY, LIGHT])


def test_two_residue_columns_swapped_is_a_failure_even_though_every_row_is_a_distribution(
    table, sequences
):
    frames = [rows_for(chain, seq) for chain, seq in sequences.items()]
    frames[0][["W", "Y"]] = frames[0][["Y", "W"]].to_numpy()

    with pytest.raises(check_wt.CheckFailed, match=r"highest-probability column is not 'top_res'"):
        check_wt.check(table(frames), PDB, [HEAVY, LIGHT])


def test_the_command_line_exits_1_and_prints_one_line_when_the_table_is_misaligned(
    table, sequences, capsys
):
    path = table([
        rows_for(HEAVY, sequences[HEAVY], offset=2),
        rows_for(LIGHT, sequences[LIGHT]),
    ])

    assert check_wt.main(["--table", str(path), "--pdb", str(PDB),
                          "--heavy", HEAVY, "--light", LIGHT]) == 1
    assert capsys.readouterr().out.startswith("MISALIGNED:")
