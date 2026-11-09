"""Score a Boltz-2 candidate screen with ipSAE and the chemistry terms of PR 24, then rank it.

Input is what `examples/track_a_screen.py fold` (cluster) or `cloud/modal_fold.py::fold_dir` (Modal) writes: batches of
folds in which the parent was folded as the record `parent`, so every candidate is compared with the parent of the same batch
and seed. A batch is one `boltz_results_inputs` folder; its seed is the `seed<k>` folder above it (0 when there is none).

    python examples/track_a_chem_rank.py score --root screen/modal_s1 --out rows_modal_s1.csv --workers 8
    python examples/track_a_chem_rank.py rank --rows rows_*.csv --candidates candidates.csv --out ranking.csv
    python examples/track_a_chem_rank.py audit --ranking ranking.csv

`score` is the CPU-heavy part and runs one process per fold. `rank` follows docs/chem_score_spec.md: candidates whose mean
ipSAE (max direction) is at least the parent's mean minus two SD are ranked by chem = z(contact) + z(coulomb), over the seeds
whose pose passes the epitope filter; the rest follow, ordered by ipSAE. A candidate with no on-site seed has no chem score
and is ranked by ipSAE. Nothing is vetoed.

With `--chem-weight W`, `rank` instead gives every candidate the one number the pick uses:

    score = z(ipsae_max) + W * chem_z,   chem_z = z(contact) + z(coulomb)

and nothing else. The other columns of ranking.csv - ipsae_min, iptm, in_band, n_on_site, para_net, chem - are reported and
do not enter the score; `ipsae_min` in particular is kept as a column and deliberately left out of the number, because it is
very nearly a multiple of `ipsae_max` and would weight the co-folding signal twice. `SCORE_INPUTS` names the three columns
the score does depend on, `audit` prints their correlation matrix, and `rank` warns on any pair that reaches 0.9.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import re
import statistics
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROW_FIELDS = ["name", "seed", "batch", "eligible", "contact", "coulomb", "epitope_overlap", "n_added_cation", "paratope_net_charge",
              "iptm", "plddt", "ipsae_max", "ipsae_min"]
PARENT = "parent"


def _fold_files(pae_path: Path) -> tuple[str, Path, Path]:
    name = pae_path.name[len("pae_"):-len("_model_0.npz")]
    return name, pae_path.parent / f"{name}_model_0.cif", pae_path.parent / f"confidence_{name}_model_0.json"


def _score_one(job: tuple) -> dict | None:
    pae_path, root, parent_cif = job
    import numpy as np

    from design_pipeline.scoring import charge_contact_report, ipsae
    from track_a_screen import cif_chains

    name, cif, confidence = _fold_files(Path(pae_path))
    if not cif.is_file() or not confidence.is_file():
        return None
    rel = str(Path(pae_path).relative_to(root))
    seed = re.search(r"seed(\d+)", rel)
    batch = str(Path(pae_path).parents[2].relative_to(root))
    summary = json.loads(confidence.read_text())
    scores = ipsae(np.load(pae_path)["pae"], cif_chains(cif))
    row = {"name": name, "seed": int(seed.group(1)) if seed else 0, "batch": batch, "iptm": round(summary["iptm"], 4),
           "plddt": round(100.0 * summary["complex_plddt"], 2), "ipsae_max": round(scores["max"], 4), "ipsae_min": round(scores["min"], 4)}
    if name != PARENT and parent_cif is not None:
        chem = charge_contact_report(cif, Path(parent_cif))
        if chem is not None:
            row.update({k: round(float(chem[k]), 4) for k in ("eligible", "contact", "coulomb", "epitope_overlap", "n_added_cation", "paratope_net_charge")})
    return row


def score(args: argparse.Namespace) -> int:
    root = Path(args.root).resolve()
    jobs = []
    for pae in sorted(root.rglob("pae_*_model_0.npz")):
        batch_predictions = pae.parents[1]
        parent_cif = batch_predictions / PARENT / f"{PARENT}_model_0.cif"
        jobs.append((str(pae), str(root), str(parent_cif) if parent_cif.is_file() else None))
    print(f"{len(jobs)} folds under {root}", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = [r for r in pool.map(_score_one, jobs, chunksize=2) if r is not None]
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROW_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} rows -> {args.out}")
    return 0


def _mean(values: Sequence[float]) -> float:
    values = [v for v in values if v is not None and not math.isnan(v)]
    return statistics.mean(values) if values else math.nan


def _z(values: dict[str, float]) -> dict[str, float]:
    mean, sd = statistics.mean(values.values()), statistics.pstdev(values.values())
    return {k: (v - mean) / (sd if sd else 1.0) for k, v in values.items()}


def rank_rows(rows: list[dict], band_sd: float = 2.0) -> tuple[list[dict], dict]:
    """One dict per candidate, best first, and the parent's mean and SD of ipSAE (max) that set the band."""
    def num(row, key):
        try:
            return float(row[key])
        except (KeyError, TypeError, ValueError):
            return math.nan

    parent = [num(r, "ipsae_max") for r in rows if r["name"] == PARENT]
    if len(parent) < 2:
        raise SystemExit("need the parent folded at two seeds or more to set the ipSAE band")
    parent_mean, parent_sd = statistics.mean(parent), statistics.stdev(parent)
    floor = parent_mean - band_sd * parent_sd
    by: dict[str, list[dict]] = {}
    for r in rows:
        if r["name"] != PARENT:
            by.setdefault(r["name"], []).append(r)
    out = []
    for name, rs in by.items():
        elig = [r for r in rs if num(r, "eligible") >= 1.0]
        out.append({"name": name, "seeds": len(rs), "ipsae_max": _mean([num(r, "ipsae_max") for r in rs]),
                    "ipsae_min": _mean([num(r, "ipsae_min") for r in rs]), "iptm": _mean([num(r, "iptm") for r in rs]),
                    "n_on_site": len(elig), "contact": _mean([num(r, "contact") for r in elig]),
                    "coulomb": _mean([num(r, "coulomb") for r in elig]), "para_net": _mean([num(r, "paratope_net_charge") for r in rs])})
    for o in out:
        o["in_band"] = o["ipsae_max"] >= floor
    scorable = {o["name"]: o for o in out if o["in_band"] and o["n_on_site"] > 0 and not math.isnan(o["contact"])}
    if len(scorable) >= 2:
        zc, zq = _z({n: o["contact"] for n, o in scorable.items()}), _z({n: o["coulomb"] for n, o in scorable.items()})
        for n, o in scorable.items():
            o["chem"] = zc[n] + zq[n]
    for o in out:
        o.setdefault("chem", math.nan)
    out.sort(key=lambda o: (0 if not math.isnan(o["chem"]) else 1 if o["in_band"] else 2,
                            -o["chem"] if not math.isnan(o["chem"]) else -o["ipsae_max"]))
    for i, o in enumerate(out, 1):
        o["rank"] = i
    return out, {"parent_mean": parent_mean, "parent_sd": parent_sd, "floor": floor, "parent_folds": len(parent)}


def rank_weighted(rows: list[dict], chem_weight: float, band_sd: float = 1.0) -> tuple[list[dict], dict]:
    """Rank by score = z(ipSAE) + chem_weight * chem, over all candidates, with no band cut.

    ipSAE is the mean over seeds (max direction). chem is z(contact) + z(coulomb) over the candidates that have an on-site
    seed; one without any gets 0 (neutral, never a penalty). A small weight keeps chemistry in the ranking without
    letting it override ipSAE. `in_band` is still reported against the parent's mean minus `band_sd` SD."""
    ranked, band = rank_rows(rows, band_sd)
    ip = _z({o["name"]: o["ipsae_max"] for o in ranked})
    have = {o["name"]: o for o in ranked if o["n_on_site"] > 0 and not math.isnan(o["contact"])}
    chem = {}
    if len(have) >= 2:
        zc, zq = _z({n: o["contact"] for n, o in have.items()}), _z({n: o["coulomb"] for n, o in have.items()})
        chem = {n: zc[n] + zq[n] for n in have}
    for o in ranked:
        o["chem_z"] = chem.get(o["name"], 0.0)
        o["score"] = ip[o["name"]] + chem_weight * o["chem_z"]
    ranked.sort(key=lambda o: -o["score"])
    for i, o in enumerate(ranked, 1):
        o["rank"] = i
    return ranked, band


# ------------------------------------------------------------------ what feeds the score, and the correlation audit
# `rank_weighted` is the whole score:   score = z(ipsae_max) + chem_weight * chem_z,   chem_z = z(contact) + z(coulomb)
# so these three columns are its only independent inputs. Everything else ranking.csv carries is reported, never scored:
# ipsae_min, iptm, in_band, n_on_site, para_net and chem do not enter the number. ipsae_min is excluded deliberately - it is
# about 0.79 * ipsae_max with little left over, so adding it would weight the co-folding signal twice against the chemistry
# terms. Keep this tuple beside the formula: the audit below reads it, so a column that joins the score without joining this
# tuple is a column that stops being checked.
SCORE_INPUTS = ("ipsae_max", "contact", "coulomb")

# Every per-candidate column worth auditing, the score's inputs first. `chem` and `chem_z` are the same formula over
# different populations (`chem` over the in-band on-site candidates, `chem_z` over all on-site ones and 0 for a candidate
# with none), so a near-1 correlation between those two is expected and is not evidence that either belongs in the score.
AUDIT_COLUMNS = ("ipsae_max", "contact", "coulomb", "chem_z", "chem", "ipsae_min", "iptm", "in_band", "n_on_site", "para_net")
DERIVED_FROM = {"chem": ("contact", "coulomb"), "chem_z": ("contact", "coulomb"), "in_band": ("ipsae_max",),
                "score": ("ipsae_max", "contact", "coulomb", "chem_z")}
CORRELATION_LIMIT = 0.9


def score_weights(chem_weight: float = 0.15) -> dict[str, float]:
    """The nominal weight each score input carries on its own z-scale. chem_z = z(contact) + z(coulomb), so each chemistry
    term carries the full `chem_weight`, not half of it."""
    return {"ipsae_max": 1.0, "contact": chem_weight, "coulomb": chem_weight}


def _numeric(value) -> float:
    """A column's value as a float. Bools and their CSV spellings count as 1 and 0; anything unparsable is nan."""
    if isinstance(value, bool):
        return 1.0 if value else 0.0
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        text = value.strip()
        if text.lower() in ("true", "false"):
            return 1.0 if text.lower() == "true" else 0.0
        try:
            return float(text)
        except ValueError:
            return math.nan
    return math.nan


def _both(records: Sequence[dict], a: str, b: str) -> list[tuple[float, float]]:
    """The rows where both columns are numbers. Pairwise-complete, so a candidate with no chemistry still counts towards
    the pairs that do not involve a chemistry column."""
    out = []
    for r in records:
        x, y = _numeric(r.get(a)), _numeric(r.get(b))
        if not (math.isnan(x) or math.isnan(y)):
            out.append((x, y))
    return out


def pearson(points: Sequence[tuple[float, float]]) -> float:
    """Pearson r, or nan for fewer than three points or a column that does not vary."""
    if len(points) < 3:
        return math.nan
    xs, ys = [p[0] for p in points], [p[1] for p in points]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sxx, syy = sum((x - mx) ** 2 for x in xs), sum((y - my) ** 2 for y in ys)
    if sxx <= 0.0 or syy <= 0.0:
        return math.nan
    return sum((x - mx) * (y - my) for x, y in points) / math.sqrt(sxx * syy)


def _ranks(values: Sequence[float]) -> list[float]:
    """1-based ranks, ties sharing their mean rank."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = shared
        i = j + 1
    return out


def spearman(points: Sequence[tuple[float, float]]) -> float:
    """Spearman rho: Pearson r of the ranks. This is the number that says whether two columns order the candidates alike."""
    if len(points) < 3:
        return math.nan
    return pearson(list(zip(_ranks([p[0] for p in points]), _ranks([p[1] for p in points]))))


def correlation_matrix(records: Sequence[dict], columns: Sequence[str] = SCORE_INPUTS,
                       method: str = "pearson") -> dict[tuple[str, str], float]:
    """{(a, b): r} for every pair of `columns`, each pair in `columns` order, over the rows where both are numbers.

    `method` is "pearson" or "spearman". A value is nan where a column is constant, absent, or the pair has fewer than
    three usable rows. `records` are the dicts `rank_weighted`/`rank_rows` return or the rows of ranking.csv."""
    correlate = spearman if method == "spearman" else pearson
    present = [c for c in columns if any(not math.isnan(_numeric(r.get(c))) for r in records)]
    return {(a, b): correlate(_both(records, a, b)) for i, a in enumerate(present) for b in present[i + 1:]}


def check_score_correlations(records: Sequence[dict], threshold: float = CORRELATION_LIMIT,
                             columns: Sequence[str] = SCORE_INPUTS, stream=None) -> list[str]:
    """Warn on every pair of score inputs that reaches `threshold` in |r| or |rho|.

    Such a pair is one quantity entering `score` twice, which silently doubles its weight against the other terms - the
    bug that kept ipsae_min out of the score. Returns the warning lines, empty when the inputs are clean, and writes them
    to `stream` (stderr by default). `rank` calls this on every run so a redundant column cannot reappear unnoticed."""
    pear = correlation_matrix(records, columns, "pearson")
    spear = correlation_matrix(records, columns, "spearman")
    lines = []
    for pair in sorted(pear, key=lambda p: -max(abs(v) for v in (pear[p], spear.get(p, math.nan), 0.0) if not math.isnan(v))):
        a, b = pair
        r, rho = pear[pair], spear.get(pair, math.nan)
        worst = max((abs(v) for v in (r, rho) if not math.isnan(v)), default=math.nan)
        if math.isnan(worst) or worst < threshold:
            continue
        parents = DERIVED_FROM.get(b, ()) + DERIVED_FROM.get(a, ())
        why = "  (expected: one is built from the other)" if a in parents or b in parents else ""
        lines.append(f"{a} and {b}: r {r:+.3f}, rho {rho:+.3f} over {len(_both(records, a, b))} candidates{why}")
    if lines:
        out = sys.stderr if stream is None else stream
        what = "score inputs" if tuple(columns) == SCORE_INPUTS else "audited columns"
        tail = ("one of each pair is weighting the same signal twice. Drop one, or say in docs/track_a_final_pipeline.md "
                "why both stay:") if what == "score inputs" else ("a pair is a double-weighting only where both members "
                "are score inputs; see docs/track_a_final_pipeline.md:")
        print(f"warning: {len(lines)} pair(s) of {what} at or above |r| {threshold:.2f}; {tail}", file=out)
        for line in lines:
            print(f"  {line}", file=out)
    return lines


def effective_weights(records: Sequence[dict], chem_weight: float = 0.15, columns: Sequence[str] = SCORE_INPUTS,
                      drop: Sequence[str] = ()) -> dict[str, float]:
    """Each input's share of the score's variance, which is what its nominal weight is really worth once the inputs
    correlate: score = sum_i w_i z_i, so Var = sum_ij w_i w_j r_ij and input i carries w_i * sum_j w_j r_ij / Var.
    Pass `drop` to see the shares after removing a column. The shares sum to 1; a missing pair counts as r = 0."""
    cols = [c for c in columns if c not in drop]
    nominal = score_weights(chem_weight)
    w = {c: nominal.get(c, 1.0) for c in cols}
    r = {}
    for i, a in enumerate(cols):
        r[(a, a)] = 1.0
        for b in cols[i + 1:]:
            value = pearson(_both(records, a, b))
            r[(a, b)] = r[(b, a)] = 0.0 if math.isnan(value) else value
    total = sum(w[a] * w[b] * r[(a, b)] for a in cols for b in cols)
    if total <= 0.0:
        return {c: math.nan for c in cols}
    return {a: w[a] * sum(w[b] * r[(a, b)] for b in cols) / total for a in cols}


def format_correlation_matrix(records: Sequence[dict], columns: Sequence[str] = AUDIT_COLUMNS, method: str = "pearson",
                              threshold: float = CORRELATION_LIMIT) -> str:
    """The square matrix as text, with a `*` on every off-diagonal cell at or above `threshold` and `.` for nan."""
    matrix = correlation_matrix(records, columns, method)
    present = sorted({c for pair in matrix for c in pair}, key=lambda c: list(columns).index(c))
    width = max((len(c) for c in present), default=0)
    head = " " * (width + 1) + "".join(f"{c[:8]:>9}" for c in present)
    lines = [head]
    for a in present:
        cells = []
        for b in present:
            if a == b:
                cells.append(f"{'1.000':>9}")
                continue
            value = matrix.get((a, b), matrix.get((b, a), math.nan))
            if math.isnan(value):
                cells.append(f"{'.':>9}")
            else:
                cells.append(f"{value:+8.3f}{'*' if abs(value) >= threshold else ' '}")
        lines.append(f"{a:<{width}} " + "".join(cells))
    return "\n".join(lines)


def audit(args: argparse.Namespace) -> int:
    """Print the correlation matrix over the audited columns of a ranking.csv and flag the redundant pairs."""
    records = []
    for pattern in args.ranking:
        for path in sorted(glob.glob(pattern)):
            with open(path, newline="") as handle:
                records += list(csv.DictReader(handle))
    columns = list(args.columns) if args.columns else list(AUDIT_COLUMNS)
    print(f"{len(records)} candidates from {', '.join(args.ranking)}")
    print(f"score = z(ipsae_max) + {args.chem_weight} * chem_z, chem_z = z(contact) + z(coulomb); score inputs "
          f"{', '.join(SCORE_INPUTS)}\n")
    for method in ("pearson", "spearman"):
        print(f"{method} over {len(columns)} columns ('*' at or above |{args.threshold:.2f}|, '.' not comparable):")
        print(format_correlation_matrix(records, columns, method, args.threshold), end="\n\n")
    flagged = check_score_correlations(records, args.threshold, columns, stream=sys.stdout)
    if not flagged:
        print(f"no pair at or above |r| {args.threshold:.2f}")
    print("\nshare of the score's variance per input:")
    before = effective_weights(records, args.chem_weight)
    for name, share in before.items():
        print(f"  {name:<12}{share:7.1%}  (nominal weight {score_weights(args.chem_weight)[name]:g})")
    for column in args.drop or ():
        if column in SCORE_INPUTS:
            print(f"\nwithout {column}:")
            for name, share in effective_weights(records, args.chem_weight, drop=(column,)).items():
                print(f"  {name:<12}{share:7.1%}")
        else:
            print(f"\n{column} is not a score input ({', '.join(SCORE_INPUTS)}); dropping it changes no weight")
    return 0


def rank(args: argparse.Namespace) -> int:
    rows = []
    for pattern in args.rows:
        for path in glob.glob(pattern):
            with open(path, newline="") as handle:
                rows += list(csv.DictReader(handle))
    meta = {}
    if args.candidates:
        with open(args.candidates, newline="") as handle:
            meta = {r["name"]: r for r in csv.DictReader(handle)}
    ranked, band = rank_weighted(rows, args.chem_weight) if args.chem_weight is not None else rank_rows(rows)
    check_score_correlations(ranked, args.correlation_limit)      # a score input that duplicates another must not pass quietly
    print(f"parent ipSAE (max) {band['parent_mean']:.3f}, SD {band['parent_sd']:.3f} over {band['parent_folds']} folds; band floor {band['floor']:.3f}")
    print(f"{len(ranked)} candidates; in band {sum(o['in_band'] for o in ranked)}; chem-scored {sum(not math.isnan(o['chem']) for o in ranked)}")
    fields = ["rank", "name", "source", "variant_id", "seeds", "ipsae_max", "ipsae_min", "iptm", "in_band", "n_on_site", "contact", "coulomb", "para_net", "chem", "chem_z", "score"]
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for o in ranked:
            m = meta.get(o["name"], {})
            writer.writerow({**{k: (round(v, 4) if isinstance(v, float) else v) for k, v in o.items() if k in fields},
                             "source": m.get("source", ""), "variant_id": m.get("variant_id", "")})
    print(f"{'#':>3} {'name':<6}{'ipSAE':>7}{'seeds':>6}{'on-site':>8}{'contact':>9}{'coulomb':>9}{'chem':>7}  variant")
    for o in ranked[: args.top]:
        m = meta.get(o["name"], {})
        print(f"{o['rank']:>3} {o['name']:<6}{o['ipsae_max']:7.3f}{o['seeds']:6d}{o['n_on_site']:8d}{o['contact']:9.3f}{o['coulomb']:9.3f}{o['chem']:7.2f}  {m.get('variant_id', '')[:70]}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    s = sub.add_parser("score", help="ipSAE and chemistry terms for every fold under --root")
    s.add_argument("--root", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--workers", type=int, default=4)
    s.set_defaults(run=score)
    r = sub.add_parser("rank", help="aggregate rows files over seeds and rank the candidates")
    r.add_argument("--rows", nargs="+", required=True, help="rows CSVs from `score` (globs allowed)")
    r.add_argument("--candidates", help="candidates.csv, to show the variant ids")
    r.add_argument("--out", required=True)
    r.add_argument("--top", type=int, default=25)
    r.add_argument("--chem-weight", type=float, help="rank by z(ipSAE) + W * chemistry over all candidates instead of the band rule (0.15 keeps chemistry in at a small weight)")
    r.add_argument("--correlation-limit", type=float, default=CORRELATION_LIMIT,
                   help=f"warn when two score inputs correlate at or above this (default {CORRELATION_LIMIT})")
    r.set_defaults(run=rank)
    a = sub.add_parser("audit", help="correlation matrix over the audited columns of a ranking.csv, and the redundant pairs")
    a.add_argument("--ranking", nargs="+", required=True, help="ranking.csv from `rank` (globs allowed)")
    a.add_argument("--columns", nargs="+", help=f"columns to audit (default: {' '.join(AUDIT_COLUMNS)})")
    a.add_argument("--drop", nargs="+", help="also report the weights with these columns removed")
    a.add_argument("--chem-weight", type=float, default=0.15, help="the weight the ranking was made with (default 0.15)")
    a.add_argument("--threshold", type=float, default=CORRELATION_LIMIT, help=f"flag pairs at or above this (default {CORRELATION_LIMIT})")
    a.set_defaults(run=audit)
    args = parser.parse_args(argv)
    return args.run(args)


if __name__ == "__main__":
    sys.exit(main())
