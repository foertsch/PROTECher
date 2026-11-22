"""`wt_argmax_rate` split by framework against CDR, for every matrix scorer.

Why this exists. `03_negative_controls.ipynb` reports one `wt_argmax_rate` per
scorer over the paratope - 21 positions, every one of them in or beside a CDR -
and two of those rates sit at chance: ablang2 4.8% and esm1v 14.3% against a
1-in-20 baseline of 5%. A rate at chance has two very different explanations and
the paratope-only number cannot tell them apart:

* the adapter is wrong (bad logit rows, a vocabulary that does not cover the 20,
  a chain offset), in which case every number from that column is meaningless;
* the model is reading the molecule correctly and simply does not prefer the
  matured residue at an affinity-selected position, which is a real finding and
  the thing a paratope is for.

The framework separates them. Framework residues are germline and heavily
conserved, so any masked antibody model that works at all recovers them at a
high rate. Hence the rule this script prints before it computes anything:

    FRAMEWORK at chance          -> the column is broken, whatever CDR says
    FRAMEWORK high, CDR at chance -> the column works and the CDR rate is real

Regions come from IMGT numbering via ANARCI (`design_pipeline.numbering`), not
from a hand-written residue range, so the split transfers to another antibody.
ANARCI labels only the variable domain, which conveniently also restricts the
scan to the positions these models are meant to see.

Run it in the environment that has the scorers:

    python examples/wt_argmax_by_region.py --out runs/negative_controls/wt_argmax_region
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from design_pipeline import numbering
from design_pipeline.api import score_variants
from design_pipeline.errors import PipelineError
from design_pipeline.structures import AMINO_ACIDS, Complex, chain_sequence
from design_pipeline.variants import enumerate_variants, make_id

# IMGT's own region boundaries, inclusive, in IMGT position numbers. The one
# place a region is defined; everything else reads these.
IMGT_CDR = {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)}
IMGT_FR = {"FR1": (1, 26), "FR2": (39, 55), "FR3": (66, 104), "FR4": (118, 128)}

CHANCE = 1.0 / len(AMINO_ACIDS)
BROKEN_AT = 2 * CHANCE      # at or below this, the column is not about this antibody
WORKING_ABOVE = 0.50        # framework above this, the column reads the molecule


def imgt_region(label: str) -> str:
    """One of FR1..FR4, CDR1..CDR3, or empty for a position with no IMGT label."""
    digits = "".join(c for c in label if c.isdigit())
    if not digits:
        return ""
    position = int(digits)
    for name, (low, high) in {**IMGT_CDR, **IMGT_FR}.items():
        if low <= position <= high:
            return name
    return ""


def wt_argmax_by(frame: pd.DataFrame, column: str, by: str) -> pd.DataFrame:
    """`diagnostics.wt_argmax_rate`, one row per value of `by`.

    The same rule as the original: a position counts only if the wild type is
    among its scored rows, the winner is the single highest score at it, and a
    position with no usable rows is skipped rather than guessed at.
    """
    out = []
    for group_value, group in frame.groupby(by, sort=True):
        wins = scored = skipped = 0
        for _, position_rows in group.groupby(["chain", "position"], sort=True):
            usable = position_rows[
                np.isfinite(np.asarray(position_rows[column], dtype=float))
            ]
            wild = set(usable["wt"].astype(str))
            if not len(usable) or not wild & set(usable["mut"].astype(str)):
                skipped += 1
                continue
            scored += 1
            best = usable.loc[usable[column].idxmax()]
            wins += str(best["mut"]) == str(best["wt"])
        out.append(
            {
                by: group_value,
                "positions": scored + skipped,
                "scored": scored,
                "skipped": skipped,
                "wins": wins,
                "rate": wins / float(scored) if scored else float("nan"),
            }
        )
    return pd.DataFrame(out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", default="data/structures/1N8Z.pdb")
    parser.add_argument("--out", default="runs/negative_controls/wt_argmax_region")
    parser.add_argument("--heavy", default="B")
    parser.add_argument("--light", default="A")
    parser.add_argument(
        "--models",
        default="antifold,ablang2,esm1v,esm2_650m,antiberty",
        help="comma-separated matrix scorers",
    )
    parser.add_argument("--device", default=None, help="torch device for the PLMs")
    parser.add_argument("--antifold-src", default=None)
    parser.add_argument("--antifold-python", default=None)
    args = parser.parse_args()

    pdb = Path(args.pdb)
    if not pdb.is_file():
        raise SystemExit(f"no structure at {pdb}")
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    models = [m.strip() for m in args.models.split(",") if m.strip()]

    print("The decision rule, stated before the numbers exist:")
    print(f"  random argmax over {len(AMINO_ACIDS)} residues is {CHANCE:.1%}.")
    print(f"  FRAMEWORK at or below {BROKEN_AT:.0%} -> the adapter is broken and its")
    print("    numbers must not be published, whatever the CDR rate says.")
    print(f"  FRAMEWORK above {WORKING_ABOVE:.0%} with CDR near chance -> the column")
    print("    works and the low paratope rate is a real property of the model.")
    print()

    if not numbering.available():
        raise SystemExit(
            "ANARCI is not installed, so there are no IMGT labels and no region "
            "split. Install it with `uv pip install 'design-pipeline[numbering]'` "
            "(it needs HMMER on PATH). Nothing here guesses a CDR boundary."
        )

    # --- the positions, and the region each one is in ------------------------
    sequences = {
        chain: chain_sequence(pdb, chain) for chain in (args.heavy, args.light)
    }
    regions: dict[tuple[str, int], str] = {}
    labels: dict[tuple[str, int], str] = {}
    for role, chain in (("heavy", args.heavy), ("light", args.light)):
        chain_labels = numbering.imgt_labels(sequences[chain])
        kept = 0
        for index, label in enumerate(chain_labels):
            region = imgt_region(label)
            if not region:
                continue
            regions[(chain, index + 1)] = region
            labels[(chain, index + 1)] = label
            kept += 1
        print(
            f"chain {chain} ({role}): {len(sequences[chain])} residues, "
            f"{kept} with an IMGT region label"
        )

    counts = pd.Series(list(regions.values())).value_counts().sort_index()
    print("\npositions per region: " + ", ".join(f"{k} {v}" for k, v in counts.items()))

    ids = []
    for chain, position in sorted(regions):
        wt = sequences[chain][position - 1]
        if wt not in AMINO_ACIDS:
            continue
        for mut in AMINO_ACIDS:
            if mut != wt:
                ids.append(make_id([(chain, position, wt, mut)]))
    print(
        f"{len(regions)} labelled positions x {len(AMINO_ACIDS) - 1} substitutions "
        f"= {len(ids)} candidates"
    )

    complex_ = Complex(
        paths=(pdb,),
        output_dir=out_dir,
        config={"model_name": "local_file", "source": str(pdb)},
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="anarci is not installed")
        table = enumerate_variants(
            complex_, args.heavy, out_dir / "variants", only=ids, structure=pdb
        )
    print(f"{table.rows} rows, {table.candidates} candidates")

    # --- score, one call per model so one failure costs one column ----------
    plm = {"device": args.device, "antibody_chains": [args.heavy, args.light]}
    kwargs_for = {
        "antifold": {
            "executable": args.antifold_python,
            "heavy": args.heavy,
            "light": args.light,
            "antifold_src": args.antifold_src,
        },
        "ablang2": dict(plm),
        "esm1v": dict(plm),
        "esm2_650m": dict(plm),
        "antiberty": dict(plm),
    }
    ran = []
    print("\nscoring:")
    for model in models:
        kwargs = {k: v for k, v in kwargs_for.get(model, {}).items() if v is not None}
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="anarci is not installed")
                table = score_variants(
                    complex_,
                    table,
                    out_dir / "scores",
                    models=[model],
                    structure=pdb,
                    max_variants=len(ids) + 10,
                    cache=True,
                    **kwargs,
                )
        except PipelineError as exc:
            head = str(exc).splitlines()[0]
            print(f"  {model:<12} SKIPPED  {type(exc).__name__}: {head}")
            continue
        column = f"score_{model}"
        one = table.frame.drop_duplicates("variant_id")
        answered = int(one[column].notna().sum()) if column in one.columns else 0
        if answered == 0:
            print(f"  {model:<12} ran and scored nothing, so treated as skipped")
            continue
        ran.append(model)
        print(f"  {model:<12} {answered:>6} of {table.candidates} candidates scored")
    if not ran:
        raise SystemExit("no scorer ran, so there is no table to report")

    # --- the wild-type rows, at the zero every matrix scorer means ----------
    wt_rows = pd.DataFrame(
        [
            {
                "variant_id": make_id(
                    [(chain, position, wt, wt)]
                ),
                "chain": chain,
                "position": position,
                "wt": wt,
                "mut": wt,
            }
            for chain, position in sorted(regions)
            for wt in [sequences[chain][position - 1]]
            if wt in AMINO_ACIDS
        ]
    )
    for model in ran:
        wt_rows[f"score_{model}"] = 0.0

    columns = ["variant_id", "chain", "position", "wt", "mut"] + [
        f"score_{m}" for m in ran
    ]
    with_wt = pd.concat(
        [
            table.frame.reindex(
                columns=[c for c in columns if c in table.frame.columns]
            ),
            wt_rows.reindex(columns=columns),
        ],
        ignore_index=True,
    )
    with_wt["region"] = [
        regions.get((chain, int(position)), "")
        for chain, position in zip(with_wt["chain"], with_wt["position"])
    ]
    with_wt["imgt"] = [
        labels.get((chain, int(position)), "")
        for chain, position in zip(with_wt["chain"], with_wt["position"])
    ]
    with_wt["part"] = np.where(
        with_wt["region"].str.startswith("CDR"), "CDR", "FRAMEWORK"
    )
    with_wt = with_wt[with_wt["region"] != ""]
    with_wt.to_csv(out_dir / "matrix_with_wt.csv", index=False)

    # --- the 5x2 table -------------------------------------------------------
    rows = []
    for model in ran:
        by_part = wt_argmax_by(with_wt, f"score_{model}", "part").set_index("part")
        row = {"model": model}
        for part in ("FRAMEWORK", "CDR"):
            if part in by_part.index:
                row[f"{part.lower()}_rate"] = by_part.loc[part, "rate"]
                row[f"{part.lower()}_scored"] = int(by_part.loc[part, "scored"])
                row[f"{part.lower()}_wins"] = int(by_part.loc[part, "wins"])
        framework = row.get("framework_rate", float("nan"))
        cdr = row.get("cdr_rate", float("nan"))
        if framework != framework:
            row["verdict"] = "not measured"
        elif framework <= BROKEN_AT:
            row["verdict"] = "ADAPTER BROKEN"
        elif framework > WORKING_ABOVE:
            row["verdict"] = "column works"
        else:
            row["verdict"] = "inconclusive"
        row["cdr_vs_chance"] = (
            "" if cdr != cdr else "at chance" if cdr <= BROKEN_AT else "above chance"
        )
        rows.append(row)
    table_5x2 = pd.DataFrame(rows).set_index("model")
    table_5x2.to_csv(out_dir / "wt_argmax_by_region.csv")

    print(f"\n{'model':<12}{'FRAMEWORK':>22}{'CDR':>20}{'verdict':>18}")
    print("-" * 72)
    for model, row in table_5x2.iterrows():
        fr = row.get("framework_rate", float("nan"))
        cd = row.get("cdr_rate", float("nan"))
        fr_text = (
            "-"
            if fr != fr
            else f"{fr:.1%} ({int(row['framework_wins'])}/{int(row['framework_scored'])})"
        )
        cd_text = (
            "-"
            if cd != cd
            else f"{cd:.1%} ({int(row['cdr_wins'])}/{int(row['cdr_scored'])})"
        )
        print(f"{model:<12}{fr_text:>22}{cd_text:>20}{row['verdict']:>18}")
    print("-" * 72)
    print(f"chance is {CHANCE:.1%}.  per-IMGT-region detail:\n")

    detail = []
    for model in ran:
        one = wt_argmax_by(with_wt, f"score_{model}", "region")
        one.insert(0, "model", model)
        detail.append(one)
    detail = pd.concat(detail, ignore_index=True)
    detail.to_csv(out_dir / "wt_argmax_by_imgt_region.csv", index=False)
    pivot = detail.pivot(index="model", columns="region", values="rate")
    print(pivot.to_string(float_format=lambda v: f"{v:.1%}"))
    print(f"\nwrote {out_dir / 'wt_argmax_by_region.csv'}")
    print(f"      {out_dir / 'wt_argmax_by_imgt_region.csv'}")
    print(f"      {out_dir / 'matrix_with_wt.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
