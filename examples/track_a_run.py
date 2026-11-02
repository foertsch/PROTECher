"""Track A pipeline v2: an evidence-led funnel, run end to end.

The model scores of the earlier funnel did not rank affinity on measured data (AntiFold, AbLang2 and FoldX on the ADAPT panel and
SKEMPI; see notebooks/03 and 04), so this funnel does not rank with them. It builds designs from substitutions that were measured
to help, excludes what was measured to hurt, and uses the models only as a veto.

    1 evidence    modules of substitutions from the evidence table: measured gains (ADAPT 2017, Moon 2016) and stability changes
                  that kept binding (Ng 2025), each with its measured fold
    2 exclusions  no hotspot position, no substitution that only appears in a variant that lost binding, no light-chain CDR3
                  change, no new cysteine or N-glycosylation sequon
    3 assemble    every combination of modules with 5 to 10 substitutions; priority = summed log10 measured fold (additive, which
                  the ADAPT panel supports within its error and which is untested across studies)
    4 veto        AntiFold `antifold_seq` and AbLang2 pseudo-log-likelihood of each set against a null of random sets of the same
                  size; vetoed only if BOTH fall in the worst `--veto-percentile` of the null
    5 pick        `pick_designs` on the priority, per arm (evidence loop sets, and loop sets with a stability module)
    6 export      scFv sequences and submission.csv through examples/track_a_export.py

Then `fold` runs Boltz-2 on the parent and the designs against an antigen fragment with the parent's cached alignments, and
`summarize` collects the scores. Every stage is timed in `run_report.json`.

    python examples/track_a_run.py design --out runs/track_a_run      # needs the cluster node (AntiFold and AbLang2 on a GPU)
    python examples/track_a_run.py fold --fasta runs/track_a_run/export/scfv_for_folding.fasta --seeds 2 --out runs/track_a_run/fold_a
    python examples/track_a_run.py summarize --out runs/track_a_run

The substitutions are in the crystal's numbering (heavy chain B, light chain A, 1-based, variable domains first).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import re
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Iterable, Sequence

Mutation = tuple[str, int, str, str]

HEAVY, LIGHT = "B", "A"
CHAIN_OF = {"H": HEAVY, "L": LIGHT}
VH_LENGTH, VL_LENGTH = 120, 107

HOTSPOTS = {
    (HEAVY, 33), (HEAVY, 50), (HEAVY, 57), (HEAVY, 99), (HEAVY, 105),     # VH Y33 R50 Y57 W99 Y105
    (LIGHT, 30), (LIGHT, 31), (LIGHT, 53), (LIGHT, 91), (LIGHT, 94),      # VL N30 T31 F53 H91 T94
}
"""The ten FoldX alanine hotspots (docs/track_a_filtering.md): never touched."""
LIGHT_CDR3_WINDOW = range(93, 98)
"""Light-chain CDR3 positions 93 to 97: every change tried there in Moon 2016 weakened binding."""
LOSS_OUTCOMES = ("worsened", "weaker", "numerically weaker", "lost")
SITE_OF = {(HEAVY, 31): "H1", (LIGHT, 28): "L1"}
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


# --------------------------------------------------------------------------------------------------------- evidence
@dataclass(frozen=True)
class Module:
    """A set of substitutions that was measured together, with the measured fold improvement of binding."""

    name: str
    substitutions: tuple[Mutation, ...]
    fold: float
    source: str
    site: str                  # H1, L1, H3 or FW

    @property
    def points(self) -> float:
        return math.log10(self.fold)


def parse_substitutions(text: str) -> tuple[Mutation, ...]:
    """'H:D31K;H:D102F' to ((B, 31, D, K), (B, 102, D, F)): the benchmark's H and L become the crystal's B and A."""
    out = []
    for token in str(text).split(";"):
        match = re.match(r"^([HL]):([A-Z])(\d+)([A-Z])$", token.strip())
        if match is None:
            raise ValueError(f"{text!r}: {token!r} is not a substitution like H:D31K")
        out.append((CHAIN_OF[match[1]], int(match[3]), match[2], match[4]))
    return tuple(sorted(out))


def _site(substitutions: Sequence[Mutation]) -> str:
    if any(c == HEAVY and 99 <= p <= 111 for c, p, _w, _m in substitutions):
        return "H3"
    sites = {SITE_OF.get((c, p)) for c, p, _w, _m in substitutions}
    sites.discard(None)
    return sites.pop() if len(sites) == 1 else "FW"


def read_rows(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_modules(rows: Sequence[dict]) -> tuple[list[Module], set[Mutation]]:
    """The modules the evidence supports, and the substitutions that only appear in variants that lost binding."""
    modules: dict[str, Module] = {}
    gains: set[Mutation] = set()
    losses: set[Mutation] = set()

    for row in rows:
        subs_text = row.get("substitutions_mature_domain_1based", "")
        if not subs_text or row["study"] not in ("adapt2017", "moon2016", "ng2025", "okazaki2023"):
            continue
        subs = parse_substitutions(subs_text)
        outcome = row["outcome"].lower()
        if row["study"] == "adapt2017" and row["metric"] == "affinity_fold_improvement" and len(subs) == 1:
            fold = float(row["value"])
            if fold > 1:
                modules[f"ADAPT {row['variant']}"] = Module(f"ADAPT {row['variant']}", subs, fold, "ADAPT 2017 (Fab, SPR)", _site(subs))
                gains.update(subs)
            else:
                losses.update(subs)
        elif row["study"] == "moon2016" and row["metric"] == "KD" and row["parent_value"]:
            fold = float(row["parent_value"]) / float(row["value"])
            if fold > 1:
                name = f"Moon {row['variant']}"
                modules[name] = Module(name, subs, fold, "Moon 2016 (IgG, SPR)", _site(subs))
                gains.update(subs)
            else:
                losses.update(subs)
        elif any(word in outcome for word in LOSS_OUTCOMES) and row["metric"] != "Tm_change":
            losses.update(subs)
    # Ng 2025 framework variants: kept if every binding readout was unchanged, with a stability gain a bonus
    by_variant: dict[str, list[dict]] = {}
    for row in rows:
        if row["study"] == "ng2025" and row.get("substitutions_mature_domain_1based"):
            by_variant.setdefault(row["variant"], []).append(row)
    for variant, group in by_variant.items():
        subs = parse_substitutions(group[0]["substitutions_mature_domain_1based"])
        binding = [r["outcome"].lower() for r in group if r["metric"] in ("cell_binding", "proliferation_inhibition", "ADCC")]
        if binding and all(o in ("unchanged",) for o in binding):
            gain = any(r["metric"] == "Tm" and r["outcome"].lower().startswith("improved") for r in group)
            name = f"Ng {variant}" + (" (Tm up)" if gain else "")
            modules[name] = Module(name, subs, 1.0, "Ng 2025 (IgG, binding kept)", "FW")
    # a substitution is excluded only if it appears in a loss and in no gain
    return list(modules.values()), losses - gains


# ------------------------------------------------------------------------------------------------------- assembly
def liabilities(parent: str, variant: str) -> dict[str, bool]:
    """New unpaired cysteine or N-glycosylation sequon excludes a design; a new NG or DG motif only flags it."""
    def count(seq: str, pattern: str) -> int:
        return len(re.findall(f"(?={pattern})", seq))

    return {
        "new_cysteine": variant.count("C") != parent.count("C"),
        "new_nglyc_sequon": count(variant, r"N[^P][ST]") > count(parent, r"N[^P][ST]"),
        "new_NG": count(variant, "NG") > count(parent, "NG"),
        "new_DG": count(variant, "DG") > count(parent, "DG"),
    }


def assemble(modules: Sequence[Module], losses: set[Mutation], min_size: int = 5, max_size: int = 10) -> list[dict]:
    """Every combination of one H3 module, optional H1 and L1 modules and optional framework modules."""
    by_site = {s: [m for m in modules if m.site == s] for s in ("H1", "L1", "H3", "FW")}
    h1 = [None] + by_site["H1"]
    l1 = [None] + by_site["L1"]
    fw_pool = by_site["FW"]
    fw_sets: list[tuple[Module, ...]] = [()]
    for k in range(1, len(fw_pool) + 1):
        from itertools import combinations
        fw_sets += list(combinations(fw_pool, k))
    sets = []
    for a, b, c, fw in product(h1, l1, by_site["H3"], fw_sets):
        parts = [m for m in (a, b, c) if m is not None] + list(fw)
        subs = sorted(s for m in parts for s in m.substitutions)
        positions = [(chain, position) for chain, position, _w, _m in subs]
        if len(set(positions)) != len(positions):
            continue                                            # two modules change the same position
        if not min_size <= len(subs) <= max_size:
            continue
        if any((chain, position) in HOTSPOTS for chain, position in positions):
            continue
        if any(s in losses for s in subs) or any(c == LIGHT and p in LIGHT_CDR3_WINDOW for c, p, _w, _m in subs):
            continue
        sets.append({
            "subs": tuple(subs),
            "modules": " + ".join(m.name for m in parts),
            "priority": sum(m.points for m in parts) + (0.1 if any("Tm up" in m.name for m in fw) else 0.0),
            "arm": "stability" if fw else "loop",
        })
    return sets


def variant_id(subs: Iterable[Mutation]) -> str:
    return "+".join(f"{c}:{w}{p}{m}" for c, p, w, m in sorted(subs))


# ----------------------------------------------------------------------------------------------------- the null
def random_sets(parent_vh: str, parent_vl: str, sizes: Sequence[int], per_size: int, seed: int = 0) -> list[tuple[Mutation, ...]]:
    """Random substitution sets of each size over the variable domains, positions distinct, hotspots allowed (a null for the
    model scores, not a design space)."""
    rng = random.Random(seed)
    domains = {HEAVY: parent_vh, LIGHT: parent_vl}
    out = []
    for size in sizes:
        for _ in range(per_size):
            sites = rng.sample([(c, p) for c, seq in domains.items() for p in range(1, len(seq) + 1)], size)
            out.append(tuple(sorted((c, p, domains[c][p - 1], rng.choice([a for a in AMINO_ACIDS if a != domains[c][p - 1]]))
                                    for c, p in sites)))
    return out


def percentile_in_null(value: float, null: Sequence[float]) -> float:
    """Fraction of the null at or below `value`: 0.02 means the design scores worse than 98% of random sets."""
    return sum(v <= value for v in null) / len(null)


def vetoed(percentiles: Sequence[float], cutoff: float) -> bool:
    """Vetoed only if every model puts the design in the worst `cutoff` of the null."""
    return all(p < cutoff for p in percentiles)


# ----------------------------------------------------------------------------------------------------- scorers
def score_antifold_seq(structure: Path, work: Path, sets: dict[str, tuple[Mutation, ...]], weights: Path) -> dict[str, float]:
    import pandas as pd
    from design_pipeline import enumerate_variants, get_structure, score_variants
    from design_pipeline.variants import _rows_for, write_table

    cx = get_structure(str(work), source=str(structure))
    base = enumerate_variants(cx, chain=HEAVY)
    rows = []
    for subs in sets.values():
        rows += _rows_for(list(subs), {})
    frame = pd.DataFrame(rows)
    frame["n_mutations"] = frame.groupby("variant_id")["position"].transform("size")
    table = write_table(frame, work / "antifold_in.parquet", output_dir=work, structure=base.structure, config={})
    scored = score_variants(cx, table, models=["antifold_seq"], heavy=HEAVY, light=LIGHT, weights=str(weights),
                            device="cuda", cache=False, output_dir=str(work / "antifold_out"))
    return scored.frame.drop_duplicates("variant_id").set_index("variant_id")["score_antifold_seq"].to_dict()


def score_ablang2(parent_vh: str, parent_vl: str, sets: dict[str, tuple[Mutation, ...]]) -> dict[str, float]:
    """Pseudo-log-likelihood change against the parent, nats of total log likelihood (notebook 04's definition)."""
    for entry in os.environ.get("ABLANG2_PYTHONPATH", "").split(os.pathsep):
        if entry and entry not in sys.path:
            sys.path.append(entry)          # ablang2 and a torch may live outside this environment
    import numpy as np
    import torch
    import ablang2

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = ablang2.pretrained(model_to_use="ablang2-paired", device=device)
    n_residues = len(parent_vh) + len(parent_vl)

    def pair(subs: Sequence[Mutation]) -> list[str]:
        heavy, light = list(parent_vh), list(parent_vl)
        for chain, position, wt, mut in subs:
            target = heavy if chain == HEAVY else light
            if target[position - 1] != wt:
                raise ValueError(f"{chain}:{position} is {target[position - 1]!r}, not {wt!r}")
            target[position - 1] = mut
        return ["".join(heavy), "".join(light)]

    names = list(sets)
    pairs = [[parent_vh, parent_vl]] + [pair(sets[n]) for n in names]
    means = np.concatenate([np.asarray(model(pairs[s:s + 8], mode="pseudo_log_likelihood")).reshape(-1) for s in range(0, len(pairs), 8)])
    return {n: float((m - means[0]) * n_residues) for n, m in zip(names, means[1:])}


# -------------------------------------------------------------------------------------------------------- stages
class Clock:
    def __init__(self) -> None:
        self.stages: dict[str, float] = {}

    def run(self, name: str):
        clock = self

        class _Timer:
            def __enter__(self_inner):
                self_inner.t = time.perf_counter()
                print(f"[{name}] ...", flush=True)

            def __exit__(self_inner, *exc):
                clock.stages[name] = round(time.perf_counter() - self_inner.t, 1)
                print(f"[{name}] {clock.stages[name]} s", flush=True)

        return _Timer()


def design(args: argparse.Namespace) -> int:
    import pandas as pd
    from design_pipeline import scfv
    from design_pipeline.combine import pick_designs
    from design_pipeline.variants import _rows_for, write_table

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    clock = Clock()
    candidate = json.loads(Path(args.candidate).read_text())
    parent = candidate["scfv"]
    layout = scfv.layout_from_parent(parent, first=HEAVY, second=LIGHT)
    parent_vh, parent_vl = parent[:VH_LENGTH], parent[len(parent) - VL_LENGTH:]
    assert len(parent) == VH_LENGTH + 15 + VL_LENGTH == 242, "the parent is not the 242-residue trastuzumab scFv"

    with clock.run("1 evidence"):
        rows = read_rows(Path(args.evidence))
        modules, losses = load_modules(rows)
        print(f"  {len(modules)} modules from the evidence table; {len(losses)} substitutions only seen in losses")
        for m in sorted(modules, key=lambda m: -m.points):
            print(f"    {m.site:<3} {m.name:<34} x{m.fold:<5.1f} {variant_id(m.substitutions)}")

    with clock.run("2+3 exclusions and assembly"):
        sets = assemble(modules, losses)
        keep = []
        from collections import Counter
        excluded: Counter = Counter()
        for s in sets:
            seq = scfv.apply_id(layout, variant_id(s["subs"]))
            flags = liabilities(parent, seq)
            s["flags"] = ",".join(k for k, v in flags.items() if v)
            s["sequence"] = seq
            if flags["new_cysteine"] or flags["new_nglyc_sequon"]:
                excluded[("new_cysteine" if flags["new_cysteine"] else "new_nglyc_sequon", s["modules"].count("Ng VH-S85N") > 0)] += 1
                continue
            keep.append(s)
        print(f"  {len(sets)} sets of 5 to 10 substitutions assembled, {len(keep)} after the liability exclusions")
        for (reason, has_ng_pair), count in sorted(excluded.items()):
            print(f"    excluded {count:>3} for {reason}" + (" (all carry the Ng S85N/R87T pair)" if has_ng_pair else ""))
        sets = keep
        for s in sets:
            s["id"] = variant_id(s["subs"])

    candidates = {s["id"]: s["subs"] for s in sets}
    sizes = sorted({len(s["subs"]) for s in sets})
    null = random_sets(parent_vh, parent_vl, sizes, args.null_per_size, seed=args.seed)
    null_sets = {f"null{i}": subs for i, subs in enumerate(null)}

    everything = {**candidates, **null_sets}
    cache_file = out / "scores.json"
    cached = json.loads(cache_file.read_text()) if args.reuse_scores and cache_file.exists() else {}
    with clock.run("4a veto: AntiFold antifold_seq"):
        if cached.get("af") and set(cached["af"]) >= set(everything):
            af_by_name = cached["af"]
            print("  reused cached scores")
        else:
            af = score_antifold_seq(Path(args.structure), out / "work", everything, Path(args.antifold_weights))
            af_by_name = {name: af[variant_id(subs)] for name, subs in everything.items()}
    with clock.run("4b veto: AbLang2 pseudo-log-likelihood"):
        if cached.get("ab") and set(cached["ab"]) >= set(everything):
            ab = cached["ab"]
            print("  reused cached scores")
        else:
            ab = score_ablang2(parent_vh, parent_vl, everything)
    cache_file.write_text(json.dumps({"af": af_by_name, "ab": ab}))

    with clock.run("4c veto percentiles"):
        by_size = {size: {"af": [], "ab": []} for size in sizes}
        for name, subs in null_sets.items():
            by_size[len(subs)]["af"].append(af_by_name[name])
            by_size[len(subs)]["ab"].append(ab[name])
        survivors = []
        for s in sets:
            size = len(s["subs"])
            s["af"], s["ab"] = af_by_name[s["id"]], ab[s["id"]]
            s["af_pct"] = percentile_in_null(s["af"], by_size[size]["af"])
            s["ab_pct"] = percentile_in_null(s["ab"], by_size[size]["ab"])
            s["vetoed"] = vetoed([s["af_pct"], s["ab_pct"]], args.veto_percentile)
            if not s["vetoed"]:
                survivors.append(s)
        print(f"  null: {args.null_per_size} random sets per size {sizes}; {len(sets) - len(survivors)} of {len(sets)} sets vetoed")

    picked: list[dict] = []
    with clock.run("5 pick"):
        for arm, n in (("loop", args.n_loop), ("stability", args.n_stability)):
            pool = [s for s in survivors if s["arm"] == arm]
            if not pool:
                print(f"  arm {arm}: no sets left")
                continue
            frame_rows = []
            for s in pool:
                for row in _rows_for(list(s["subs"]), {}):
                    frame_rows.append({**row, "priority": s["priority"], "n_mutations": len(s["subs"])})
            table = write_table(pd.DataFrame(frame_rows), out / f"pool_{arm}.parquet", output_dir=out, structure=None, config={})
            want = min(n, len(pool))
            for needed in range(args.min_differences, 0, -1):
                chosen = pick_designs(table, by="priority", n=want, min_differences=needed, max_uses=args.max_uses, min_uses=1,
                                      max_pi_shift=None, output_dir=out / f"picked_{arm}_md{needed}")
                ids = list(dict.fromkeys(chosen.frame["variant_id"]))
                print(f"  arm {arm}: min_differences {needed} gives {len(ids)} of {want} designs")
                if len(ids) >= want:
                    break
            by_id = {s["id"]: s for s in pool}
            picked += [{**by_id[i], "arm": arm, "min_differences": needed} for i in ids]
            print(f"  arm {arm}: {len(ids)} designs from {len(pool)} sets (spread rule relaxed to {needed} positions)")

    with clock.run("6 export"):
        export_dir = out / "export"
        parent_csv = out / "parent.csv"
        parent_csv.write_text(f"trastuzumab_scfv,{parent}\n")
        ids = [s["id"] for s in picked]
        command = [sys.executable, str(Path(__file__).with_name("track_a_export.py")), "--parent", str(parent_csv), "--ids", *ids,
                   "--out", str(export_dir)]
        if args.structure_pdb:
            command += ["--structure", args.structure_pdb]
        done = subprocess.run(command, capture_output=True, text=True)
        print(done.stdout[-1500:] + done.stderr[-1500:])
        if done.returncode != 0:
            raise SystemExit(f"track_a_export.py failed with {done.returncode}")

    with open(out / "designs_final.csv", "w", newline="") as handle:
        fields = ["design", "arm", "size", "priority", "modules", "variant_id", "antifold_seq", "antifold_pct", "ablang2_pll",
                  "ablang2_pct", "flags", "min_differences", "sequence"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for i, s in enumerate(picked, start=1):
            writer.writerow({"design": f"design_{i:02d}", "arm": s["arm"], "size": len(s["subs"]), "priority": round(s["priority"], 3),
                             "modules": s["modules"], "variant_id": s["id"], "antifold_seq": round(s["af"], 2),
                             "antifold_pct": round(s["af_pct"], 3), "ablang2_pll": round(s["ab"], 2), "ablang2_pct": round(s["ab_pct"], 3),
                             "flags": s["flags"], "min_differences": s["min_differences"], "sequence": s["sequence"]})
    report = {"stages_seconds": clock.stages, "total_seconds": round(sum(clock.stages.values()), 1), "n_designs": len(picked),
              "n_modules": len(modules), "n_sets": len(sets), "n_vetoed": len(sets) - len(survivors), "args": vars(args)}
    (out / "run_report.json").write_text(json.dumps(report, indent=1, default=str))
    print(json.dumps(report["stages_seconds"], indent=1), f"\ntotal {report['total_seconds']} s")
    return 0


# ----------------------------------------------------------------------------------------------------- the fold
def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, name, chunks = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name, chunks = line[1:].strip(), []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def fold(args: argparse.Namespace) -> int:
    from design_pipeline import Design
    from design_pipeline.denovo import predict_structure
    from msa_cache import MsaCache, use_cached_msas

    out = Path(args.out)
    candidate = json.loads(Path(args.candidate).read_text())
    start, end = (int(x) for x in args.antigen_region.split("-"))
    antigen = candidate["antigen"][start - 1:end]
    records = read_fasta(Path(args.fasta))
    first, last = (int(x) for x in args.subset.split("-")) if args.subset else (0, len(records) - 1)
    records = records[first:last + 1]
    inputs = out / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    designs = []
    for name, scfv in records:
        path = inputs / f"{name}.fa"
        path.write_text(f">{name}\n{scfv}:{antigen}\n")
        designs.append(Design(id=name, path=path, parent_id=None, producer="track_a_run"))
    print(f"folding {len(designs)} complexes ({[d.id for d in designs]}), antigen {start}-{end} ({len(antigen)} residues), seeds {args.seeds}", flush=True)
    cache = MsaCache(args.msa_cache)
    rows, seconds = [], {}
    for seed in range(args.seeds):
        began = time.perf_counter()
        with use_cached_msas(cache):
            done = predict_structure(designs, output_dir=out / "predictions" / f"seed{seed}", model_name="boltz", seed=seed,
                                     use_msa_server=False)
        seconds[f"seed{seed}"] = round(time.perf_counter() - began, 1)
        for d in done:
            rows.append({"complex": d.parent_id, "seed": seed, "iptm": d.scores.get("iptm"), "ptm": d.scores.get("ptm"),
                         "plddt": d.scores.get("plddt"), "path": str(d.path)})
        print(f"seed {seed}: {seconds[f'seed{seed}']} s for {len(designs)} folds", flush=True)
        write_summary(rows, out / "fold_summary.csv")
    (out / "fold_timing.json").write_text(json.dumps({"seconds": seconds, "complexes": len(designs)}, indent=1))
    return 0


def write_summary(rows: Sequence[dict], path: Path) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["complex", "seed", "iptm", "ptm", "plddt", "path"])
        writer.writeheader()
        writer.writerows(rows)


def summarize(args: argparse.Namespace) -> int:
    out = Path(args.out)
    rows = []
    for summary in sorted(out.glob("fold_*/fold_summary.csv")):
        rows += list(csv.DictReader(open(summary)))
    by: dict[str, list[dict]] = {}
    for row in rows:
        by.setdefault(row["complex"], []).append(row)
    designs = {r["design"]: r for r in csv.DictReader(open(out / "designs_final.csv"))}
    parent_name = next((c for c in by if c not in designs), None)
    parent = [float(r["iptm"]) for r in by.get(parent_name, []) if r["iptm"]]
    parent_mean = statistics.mean(parent) if parent else float("nan")
    parent_sd = statistics.stdev(parent) if len(parent) > 1 else float("nan")
    print(f"parent ipTM {parent_mean:.3f} (SD {parent_sd:.3f} over {len(parent)} folds)")
    print(f"{'design':<10}{'arm':<10}{'n':>3}{'ipTM mean':>11}{'SD':>7}{'vs parent':>11}{'pLDDT':>8}   modules")
    lines = []
    for name, d in designs.items():
        folds = by.get(name, [])
        iptm = [float(r["iptm"]) for r in folds if r["iptm"]]
        plddt = [float(r["plddt"]) for r in folds if r["plddt"]]
        if not iptm:
            print(f"{name:<10}{d['arm']:<10}{d['size']:>3}   not folded")
            continue
        m = statistics.mean(iptm)
        sd = statistics.stdev(iptm) if len(iptm) > 1 else float("nan")
        lines.append({"design": name, "arm": d["arm"], "size": d["size"], "iptm_mean": round(m, 4), "iptm_sd": round(sd, 4),
                      "iptm_vs_parent": round(m - parent_mean, 4), "plddt_mean": round(statistics.mean(plddt), 2) if plddt else "",
                      "folds": len(iptm), "modules": d["modules"]})
        print(f"{name:<10}{d['arm']:<10}{d['size']:>3}{m:>11.3f}{sd:>7.3f}{m - parent_mean:>+11.3f}{(statistics.mean(plddt) if plddt else float('nan')):>8.1f}   {d['modules'][:60]}")
    with open(out / "fold_results.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(lines[0]))
        writer.writeheader()
        writer.writerows(lines)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    d = sub.add_parser("design", help="stages 1 to 6")
    d.add_argument("--out", required=True)
    d.add_argument("--evidence", default="trastuzumab_mutation_evidence.csv")
    d.add_argument("--candidate", default="/scratch/foertsch/protecher_cache/msa_candidates/trastuzumab_her2/candidate.json")
    d.add_argument("--structure", default="/scratch/foertsch/protecher_cache/structures/1N8Z.cif", help="the crystal AntiFold reads")
    d.add_argument("--structure-pdb", default="data/structures/1N8Z.pdb", help="the crystal the export checks the parent against")
    d.add_argument("--antifold-weights", default="/scratch/foertsch/protecher_tools/antifold/src/models/model.pt")
    d.add_argument("--null-per-size", type=int, default=60)
    d.add_argument("--veto-percentile", type=float, default=0.05)
    d.add_argument("--n-loop", type=int, default=9)
    d.add_argument("--n-stability", type=int, default=6)
    d.add_argument("--min-differences", type=int, default=3)
    d.add_argument("--max-uses", type=int, default=8)
    d.add_argument("--seed", type=int, default=0)
    d.add_argument("--reuse-scores", action="store_true", help="read scores.json from --out instead of rescoring")
    d.set_defaults(func=design)

    f = sub.add_parser("fold", help="Boltz-2 on the parent and the designs with the parent's cached alignments")
    f.add_argument("--fasta", required=True, help="scfv_for_folding.fasta: the parent first, then the designs")
    f.add_argument("--out", required=True)
    f.add_argument("--subset", help="FIRST-LAST record numbers, 0-based and inclusive (default all)")
    f.add_argument("--seeds", type=int, default=2)
    f.add_argument("--antigen-region", default="480-607", help="START-END of the 630-residue HER2 ectodomain, 1-based")
    f.add_argument("--candidate", default="/scratch/foertsch/protecher_cache/msa_candidates/trastuzumab_her2/candidate.json")
    f.add_argument("--msa-cache", default="/scratch/foertsch/protecher_cache/msa_candidates/trastuzumab_her2")
    f.set_defaults(func=fold)

    s = sub.add_parser("summarize", help="collect the fold summaries under --out")
    s.add_argument("--out", required=True)
    s.set_defaults(func=summarize)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
