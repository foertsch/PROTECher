"""Does AntiFold's per-residue table line up with the structure it was run on?

The check notebook 03 cannot make of itself. The notebook builds each
`LogProbMatrix` with `from_position_mapping(chain, per_position, sequence)`,
where `sequence` is the chain read from the PDB - so the matrix carries the
PDB's residues by construction and `assert_matches` on it is trivially true.
This script builds the matrix with the sequence AntiFold itself reports
(`pdb_res`, row by row) and then runs the same round-trip assertion against the
PDB. The two sequences come from different files, so the assertion can fail,
which is the point: a matrix shifted by one position still yields real-looking
log-odds, for the wrong residue.

Four things are checked per antibody chain, and all four have to hold:

1. **the round trip** - the table's own wild-type column, in file order, is the
   PDB chain's sequence (`LogProbMatrix.assert_matches`);
2. **the shift profile** - no other offset in -3..+3 fits as well, so the
   agreement at offset 0 is not an accident of composition;
3. **the labels** - every `pdb_posins` maps to a position of the chain, and the
   residue there is the one the table claims, which is the mapping the adapter
   indexes with;
4. **the columns** - the residue with the highest log probability in each row
   is the one AntiFold names in `top_res`, so the 20 columns are the residues
   this package thinks they are.

`top_res` disagreeing with `pdb_res` is *not* a failure: that is the model
preferring a different residue to the one the crystal has, which is a result,
not a misalignment. It is counted and printed.

    python -m antifold_score.check_wt --pdb data/structures/1N8Z.pdb --heavy B --light A
        --table runs/adapt_calibration/antifold/1N8Z_BA.csv

Exits 0 if every check passes, 1 otherwise.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from design_pipeline.matrices import IndexMismatchError, LogProbMatrix
from design_pipeline.structures import AMINO_ACIDS, chain_residue_ids, chain_sequence

CHAIN_COLUMNS = ("pdb_chain", "chain")
POSITION_COLUMNS = ("pdb_pos", "position", "pdb_res_num", "resnum")
LABEL_COLUMNS = ("pdb_posins", "pdb_pos_ins", "posins")
WT_COLUMNS = ("pdb_res", "aa_orig", "wt", "aa_wt")
TOP_COLUMNS = ("top_res", "top_aa")
SHIFTS = range(-3, 4)


class CheckFailed(ValueError):
    """A check did not hold. The message names the chain and the position."""


def _column(frame, names, what, required=True):
    for name in names:
        if name in frame.columns:
            return name
    if required:
        raise CheckFailed(f"no {what} column among {names} in {list(frame.columns)}")
    return None


def _matrix_from_table(chain: str, rows: pd.DataFrame, wt_column: str) -> LogProbMatrix:
    """The matrix as AntiFold reports it: its own wild-type column as the sequence.

    This is what makes the round trip real. `from_position_mapping` is not used
    here on purpose - it takes the chain sequence as an argument, which is the
    answer the check is supposed to arrive at.
    """
    declared = "".join(str(v).strip().upper() for v in rows[wt_column])
    values = rows[list(AMINO_ACIDS)].to_numpy(dtype=float)
    return LogProbMatrix(chain=chain, values=values, sequence=declared)


def _shift_profile(table_sequence: str, structure_sequence: str) -> dict[int, float]:
    """Identity of the table's sequence against the chain at each offset."""
    profile = {}
    for shift in SHIFTS:
        start = max(0, shift)
        window = structure_sequence[start : start + len(table_sequence)]
        compared = table_sequence[max(0, -shift) :][: len(window)]
        window = window[: len(compared)]
        profile[shift] = (
            sum(a == b for a, b in zip(window, compared)) / len(compared) if compared else 0.0
        )
    return profile


def check_chain(
    chain: str,
    rows: pd.DataFrame,
    sequence: str,
    label_ids,
    columns: dict,
) -> list[str]:
    """Every check for one chain. Returns the lines to print; raises on failure."""
    out = [f"chain {chain}: {len(rows)} table rows, {len(sequence)} residues in the structure"]

    # 1. the round trip, on a matrix carrying AntiFold's own wild-type column
    matrix = _matrix_from_table(chain, rows, columns["wt"])
    if len(matrix.sequence) != len(sequence):
        raise CheckFailed(
            f"chain {chain}: AntiFold reports {len(matrix.sequence)} residues, the structure "
            f"has {len(sequence)}. The table does not describe this chain."
        )
    matrix.assert_matches(sequence)  # raises IndexMismatchError naming the first difference
    out.append(
        f"  round trip:  the table's {columns['wt']!r} column is the chain, all "
        f"{len(sequence)} positions, offset 0"
    )

    # 2. the shift profile: offset 0 has to be the only one that fits
    profile = _shift_profile(matrix.sequence, sequence)
    neighbours = {s: v for s, v in profile.items() if s != 0}
    worst = max(neighbours, key=neighbours.get)
    if neighbours[worst] >= 0.95:
        raise CheckFailed(
            f"chain {chain}: the table also matches the chain at offset {worst:+d} "
            f"({neighbours[worst]:.1%}); offset 0 is not distinguishable, so the alignment "
            "is not evidence of anything."
        )
    out.append("  shifts:      " + "  ".join(f"{s:+d}:{v:.0%}" for s, v in profile.items()))

    # 3. the labels, which is what the adapter indexes with
    index_of = {label: i + 1 for i, label in enumerate(label_ids)}
    label_column = columns["label"] or columns["position"]
    unmapped, wrong = [], []
    for _, row in rows.iterrows():
        label = str(row[label_column]).strip()
        position = index_of.get(label)
        if position is None:
            unmapped.append(label)
            continue
        declared = str(row[columns["wt"]]).strip().upper()
        if declared != sequence[position - 1]:
            wrong.append((label, position, declared, sequence[position - 1]))
    if unmapped:
        raise CheckFailed(
            f"chain {chain}: {len(unmapped)} label(s) AntiFold reports that the chain does not "
            f"have, e.g. {unmapped[:8]}"
        )
    if wrong:
        raise CheckFailed(
            f"chain {chain}: {len(wrong)} label(s) map to the wrong residue, e.g. "
            + "; ".join(
                f"label {l} -> position {p}: table says {d}, chain has {h}"
                for l, p, d, h in wrong[:5]
            )
        )
    out.append(
        f"  labels:      every {label_column} maps to a position holding the residue the "
        "table claims"
    )

    # 4. the columns, against AntiFold's own argmax
    top_column = columns.get("top")
    if top_column:
        values = rows[list(AMINO_ACIDS)].to_numpy(dtype=float)
        implied = [AMINO_ACIDS[i] for i in np.argmax(values, axis=1)]
        reported = [str(v).strip().upper() for v in rows[top_column]]
        disagree = [(i + 1, a, b) for i, (a, b) in enumerate(zip(implied, reported)) if a != b]
        if disagree:
            raise CheckFailed(
                f"chain {chain}: the highest-probability column is not {top_column!r} at "
                f"{len(disagree)} row(s), e.g. "
                + "; ".join(
                    f"row {i}: columns say {a}, table says {b}" for i, a, b in disagree[:5]
                )
                + f". The 20 columns are not {AMINO_ACIDS} in that order, so every score off "
                "this matrix would be a log-odds for the wrong residue."
            )
        out.append(
            f"  columns:     the argmax of the 20 columns is {top_column!r} at every row, so "
            f"they are {AMINO_ACIDS}"
        )
        recovered = sum(a == b for a, b in zip(implied, matrix.sequence))
        out.append(
            f"  (not a check) AntiFold's top residue is the crystal's at "
            f"{recovered}/{len(implied)} positions ({recovered / len(implied):.0%})"
        )
    return out


def check(table: Path, pdb: Path, chains, imgt_pdb: Path | None = None) -> list[str]:
    """Every check for every chain. Returns the report; raises on the first failure."""
    frame = pd.read_csv(table)
    columns = {
        "chain": _column(frame, CHAIN_COLUMNS, "chain"),
        "position": _column(frame, POSITION_COLUMNS, "position"),
        "label": _column(frame, LABEL_COLUMNS, "pdb label", required=False),
        "wt": _column(frame, WT_COLUMNS, "wild-type residue"),
        "top": _column(frame, TOP_COLUMNS, "top residue", required=False),
    }
    absent = [r for r in AMINO_ACIDS if r not in frame.columns]
    if absent:
        raise CheckFailed(
            f"{table}: no column for residue(s) {''.join(absent)}; the 20 log-probability "
            f"columns this package expects are {AMINO_ACIDS}"
        )

    if imgt_pdb is None:
        candidate = table.parent / "imgt_numbered" / pdb.name
        imgt_pdb = candidate if candidate.is_file() else None
    lines = [
        f"table:     {table}  ({len(frame)} rows)",
        f"structure: {pdb}",
        f"labels:    {imgt_pdb or pdb}"
        + ("" if imgt_pdb else "  (no renumbered copy beside the table)"),
        "",
    ]

    for chain in chains:
        rows = frame[frame[columns["chain"]].astype(str) == chain]
        if rows.empty:
            raise CheckFailed(f"the table has no rows for chain {chain!r}")
        sequence = chain_sequence(pdb, chain)
        if imgt_pdb is not None:
            if chain_sequence(imgt_pdb, chain) != sequence:
                raise CheckFailed(
                    f"chain {chain} of {imgt_pdb} is not the input chain's sequence, so its "
                    "labels cannot be mapped onto the input positions by order"
                )
            label_ids = chain_residue_ids(imgt_pdb, chain)
        else:
            label_ids = chain_residue_ids(pdb, chain)
        lines.extend(check_chain(chain, rows, sequence, label_ids, columns))
        lines.append("")
    return lines


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--table", type=Path, required=True,
        help="AntiFold's per-residue CSV, e.g. runs/adapt_calibration/antifold/1N8Z_BA.csv",
    )
    parser.add_argument("--pdb", type=Path, required=True,
                        help="the structure AntiFold was run on")
    parser.add_argument("--heavy", required=True, help="heavy chain id in that structure")
    parser.add_argument("--light", required=True, help="light chain id")
    parser.add_argument(
        "--imgt-pdb", type=Path, default=None,
        help="AntiFold's renumbered copy; the default looks for imgt_numbered/<pdb> "
             "beside the table",
    )
    args = parser.parse_args(argv)

    try:
        lines = check(args.table, args.pdb, [args.heavy, args.light], args.imgt_pdb)
    except (CheckFailed, IndexMismatchError) as exc:
        print(f"MISALIGNED: {exc}")
        return 1
    print("\n".join(lines))
    print("every check passed: the matrix describes the positions it claims to.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
