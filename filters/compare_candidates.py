#!/usr/bin/env python3
"""Join complete-candidate model results and developability evidence, without inference.

All score columns retain their source/protocol and direction. Missing evidence is
not zero or a failed experiment. Thresholds and shortlist priorities are explicit
configuration, not learned weights or assumptions about the team's model choice.
See filters/MARIA_HANDOFF.md; demo inputs are synthetic, not submission designs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from statistics import median, stdev

import pandas as pd

AAS = set("ACDEFGHIKLMNPQRSTVWY")


def digest(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if not fields or len(set(fields)) != len(fields):
            raise ValueError(f"{path}: missing or duplicate column names")
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(f"{path}: malformed CSV row")
    return fields, rows


def require_columns(fields, columns, path):
    missing = set(columns) - set(fields)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")


def clean_sequence(value: str) -> str:
    seq = "".join(value.split()).upper()
    if not seq or set(seq) - AAS:
        raise ValueError("Empty sequence or non-standard amino acids in candidate manifest")
    return seq


def read_candidates(path: Path) -> dict[str, str]:
    """Authoritative whole-scFv sequences; no inferred IDs or parent sequence."""
    if path.suffix.lower() == ".csv":
        fields, rows = read_csv(path)
        require_columns(fields, ["name", "sequence"], path)
        pairs = [(r["name"], r["sequence"]) for r in rows]
    else:
        pairs, name, parts = [], None, []
        for line in path.read_text(encoding="utf-8-sig").splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    pairs.append((name, "".join(parts)))
                name, parts = line[1:].strip(), []
            elif name is None:
                raise ValueError("FASTA sequence before first header")
            else:
                parts.append(line)
        if name is not None:
            pairs.append((name, "".join(parts)))
    result = {}
    for name, seq in pairs:
        name = name.strip()
        if not name or name in result:
            raise ValueError(f"Empty/duplicate candidate ID: {name!r}")
        result[name] = clean_sequence(seq)
    if not result:
        raise ValueError("No candidates in manifest")
    return result


def number(value, context: str) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context}: not numeric: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{context}: non-finite number")
    return parsed


def validate_config(config: dict) -> dict[str, dict]:
    allowed = {"schema_version", "parent_id", "require_developability", "sources", "gates", "shortlist"}
    if set(config) - allowed or config.get("schema_version") != 1:
        raise ValueError("Config requires schema_version=1 and only documented top-level keys")
    if not isinstance(config.get("require_developability", True), bool):
        raise ValueError("require_developability must be a boolean")
    metrics, names = {}, set()
    for source in config.get("sources", []):
        source_keys = {"name", "path", "id_column", "sequence_column", "hash_column", "identity_file",
                       "replicate_column", "expected_replicates", "protocol", "protocol_column", "metrics"}
        if set(source) - source_keys:
            raise ValueError(f"Unknown source config keys: {sorted(set(source) - source_keys)}")
        name = source["name"]
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name) or name in names:
            raise ValueError(f"Invalid/duplicate source name: {name}")
        names.add(name)
        if not source.get("protocol", "").strip() or not source.get("path"):
            raise ValueError(f"{name}: path and protocol description are required")
        identity_options = sum(bool(source.get(k)) for k in ("sequence_column", "hash_column", "identity_file"))
        if identity_options != 1:
            raise ValueError(f"{name}: specify exactly one sequence_column, hash_column or identity_file")
        expected = source.get("expected_replicates", 1)
        if type(expected) is not int or expected < 1:
            raise ValueError(f"{name}: expected_replicates must be a positive integer")
        if expected > 1 and not source.get("replicate_column"):
            raise ValueError(f"{name}: repeated predictions need replicate_column")
        if not source.get("metrics"):
            raise ValueError(f"{name}: no metrics configured")
        for metric in source["metrics"]:
            if set(metric) - {"name", "column", "better", "required", "units"}:
                raise ValueError("Unknown metric config keys")
            key = f"{name}.{metric['name']}"
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", metric["name"]) or key in metrics:
                raise ValueError(f"Duplicate/invalid metric: {key}")
            if metric.get("better") not in ("higher", "lower", "descriptive"):
                raise ValueError(f"{key}: better must be higher, lower or descriptive")
            if not metric.get("column") or not isinstance(metric.get("required", False), bool):
                raise ValueError(f"{key}: column and boolean required flag expected")
            metrics[key] = {**metric, "source": name, "protocol": source["protocol"], "expected": expected}
    for gate in config.get("gates", []):
        if set(gate) != {"metric", "operator", "threshold", "action"}:
            raise ValueError("Gate keys must be metric, operator, threshold and action")
        if gate["metric"] not in metrics or gate.get("operator") not in (">=", "<="):
            raise ValueError("Gate requires a known metric and >= or <= operator")
        if gate.get("action") not in ("reject", "review") or number(gate.get("threshold"), "gate") is None:
            raise ValueError("Gate requires finite threshold and reject/review action")
    selection = config.get("shortlist")
    if selection:
        if set(selection) - {"size", "min_sequence_distance", "sort_by", "include_review"}:
            raise ValueError("Unknown shortlist config keys")
        if type(selection.get("size")) is not int or selection["size"] < 1:
            raise ValueError("shortlist.size must be a positive integer")
        distance = selection.get("min_sequence_distance", 1)
        if type(distance) is not int or distance < 1:
            raise ValueError("shortlist.min_sequence_distance must be >= 1")
        if not isinstance(selection.get("include_review", False), bool):
            raise ValueError("shortlist.include_review must be boolean")
        keys = selection.get("sort_by", [])
        if not keys or len(set(keys)) != len(keys) or any(k not in metrics or metrics[k]["better"] == "descriptive" for k in keys):
            raise ValueError("shortlist.sort_by needs unique known directional metrics")
        if "ipsae" in keys[0].lower():
            raise ValueError("Use ipSAE as a structural gate or tie-breaker, not the primary shortlist ranking")
    return metrics


def verify_identity(name, value, candidates, *, hashed=False):
    if name not in candidates:
        raise ValueError(f"Unknown candidate ID: {name!r}")
    expected = digest(candidates[name]) if hashed else candidates[name]
    found = value.strip().lower() if hashed else clean_sequence(value)
    if found != expected:
        raise ValueError(f"Sequence identity mismatch for {name}; possible stale results")


def load_source(source, candidates, base: Path):
    """Return raw observations after verifying identity, seed uniqueness and protocol."""
    path = base / source["path"]
    if not path.exists():
        return [], {"source": source["name"], "status": "file_missing", "path": str(path)}
    if source["protocol"].startswith("SET_"):
        raise ValueError(f"{source['name']}: replace the template protocol with the actual run method")
    aliases = {}
    if source.get("identity_file"):
        mapping = base / source["identity_file"]
        fields, rows = read_csv(mapping)
        require_columns(fields, ["source_id", "design", "sequence"], mapping)
        for row in rows:
            key = row["source_id"].strip()
            if not key or key in aliases:
                raise ValueError("Duplicate/empty source_id in run identity file")
            verify_identity(row["design"], row["sequence"], candidates)
            aliases[key] = row["design"]
    fields, rows = read_csv(path)
    id_col = source.get("id_column", "design")
    needed = [id_col] + [m["column"] for m in source["metrics"]]
    needed += [source[k] for k in ("sequence_column", "hash_column", "replicate_column", "protocol_column") if source.get(k)]
    require_columns(fields, needed, path)
    seen, observations = set(), []
    for row in rows:
        source_id = row[id_col].strip()
        if source.get("identity_file"):
            if source_id not in aliases:
                raise ValueError(f"{source_id}: no identity mapping recorded for this run")
            name = aliases[source_id]
        else:
            name = source_id
            hashed = bool(source.get("hash_column"))
            verify_identity(name, row[source.get("hash_column") or source["sequence_column"]], candidates, hashed=hashed)
        seed = row[source["replicate_column"]].strip() if source.get("replicate_column") else "single"
        if not seed or (name, seed) in seen:
            raise ValueError(f"{source['name']}: duplicate/empty replicate for {name}: {seed!r}")
        seen.add((name, seed))
        if source.get("protocol_column") and row[source["protocol_column"]] != source["protocol"]:
            raise ValueError(f"{source['name']}: mixed protocols; split them into separate sources")
        for metric in source["metrics"]:
            key = f"{source['name']}.{metric['name']}"
            observations.append({"design": name, "sequence_sha256": digest(candidates[name]),
                                 "metric": key, "replicate": seed,
                                 "value": number(row[metric["column"]], f"{name}/{key}"),
                                 "protocol": source["protocol"]})
    return observations, {"source": source["name"], "status": "loaded", "path": str(path), "rows": len(rows)}


def load_developability(path: Path | None, candidates, parent_id):
    if path is None or not path.exists():
        return {}
    fields, rows = read_csv(path)
    require_columns(fields, ["design", "sequence", "input_status", "rules_status", "flags", "data_warnings"], path)
    result = {}
    for row in rows:
        name = row["design"]
        # Iason's report labels its parent PARENT; match only a declared reference by sequence.
        is_reference = row.get("is_reference", "").lower() == "true"
        if is_reference and parent_id and row["sequence"] == candidates[parent_id]:
            name = parent_id
        verify_identity(name, row["sequence"], candidates)
        if name in result:
            raise ValueError(f"Duplicate developability row: {name}")
        if row["input_status"] not in ("valid", "invalid") or row["rules_status"] not in (
            "pass", "fail", "reference", "review_required", "not_evaluated"
        ):
            raise ValueError(f"{name}: unknown developability status")
        if row["rules_status"] == "reference" and name != parent_id:
            raise ValueError("Developability reference must match configured parent_id")
        result[name] = row
    return result


def compare(candidates, config, observations, dev):
    metrics = validate_config(config)
    parent = config.get("parent_id")
    if parent is not None and parent not in candidates:
        raise ValueError("Configured parent_id is absent from the candidate manifest")
    grouped = defaultdict(list)
    for observation in observations:
        grouped[observation["design"], observation["metric"]].append(observation["value"])
    duplicates = defaultdict(list)
    for name, seq in candidates.items():
        if name != parent:
            duplicates[seq].append(name)
    rows = []
    for name, seq in candidates.items():
        rejected, pending, warnings = [], [], []
        row = {"design": name, "sequence": seq, "sequence_sha256": digest(seq), "is_reference": name == parent,
               "n_mutations": sum(a != b for a, b in zip(seq, candidates[parent]))
               if parent and len(seq) == len(candidates[parent]) else None,
               "duplicate_ids": ";".join(x for x in duplicates[seq] if x != name) if name != parent else ""}
        if row["duplicate_ids"]:
            warnings.append("duplicate sequence: " + row["duplicate_ids"])
        for key, spec in metrics.items():
            values = [x for x in grouped[name, key] if x is not None]
            row[key] = median(values) if values else None
            row[key + ".n"] = len(values)
            row[key + ".min"] = min(values) if values else None
            row[key + ".max"] = max(values) if values else None
            row[key + ".sd"] = stdev(values) if len(values) > 1 else None
            row[key + ".status"] = "missing" if not values else "partial" if len(values) < spec["expected"] else "available"
            if spec.get("required") and row[key + ".status"] != "available":
                pending.append(f"{key}: {len(values)}/{spec['expected']} required observations")
        evidence = dev.get(name)
        row["developability_status"] = "missing" if evidence is None else evidence["rules_status"]
        for column in ("flags", "rule_issues", "data_warnings", "review_notes"):
            row["developability_" + column] = evidence.get(column, "") if evidence else ""
        if evidence is None:
            if config.get("require_developability", True):
                pending.append("developability report not supplied")
        else:
            if evidence["input_status"] == "invalid" or evidence["rules_status"] == "fail":
                rejected.append("sequence-rule violation: " + evidence.get("rule_issues", ""))
            elif evidence["rules_status"] in ("review_required", "not_evaluated"):
                pending.append("sequence constraints not fully evaluated")
            if evidence["flags"]:
                warnings.append("developability warnings: " + evidence["flags"])
            if evidence["data_warnings"]:
                warnings.append("developability evidence incomplete: " + evidence["data_warnings"])
        for gate in config.get("gates", []):
            key, value = gate["metric"], row[gate["metric"]]
            if row[key + ".status"] != "available":
                pending.append(f"cannot evaluate gate: {key}")
                continue
            threshold = float(gate["threshold"])
            passed = value >= threshold if gate["operator"] == ">=" else value <= threshold
            if not passed:
                reason = f"{key}={value:g} does not satisfy {gate['operator']} {threshold:g}"
                (rejected if gate["action"] == "reject" else warnings).append(reason)
        row["decision"] = ("REFERENCE" if name == parent else "REJECT" if rejected else
                           "INCOMPLETE" if pending else "REVIEW" if warnings else "ELIGIBLE")
        row["reject_reasons"], row["pending_reasons"], row["review_reasons"] = map("; ".join, (rejected, pending, warnings))
        rows.append(row)
    frame = pd.DataFrame(rows).set_index("design", drop=False)
    nonref = ~frame["is_reference"]
    for key, spec in metrics.items():
        frame[key] = pd.to_numeric(frame[key], errors="coerce")
        frame[key + ".delta_parent"] = frame[key] - frame.at[parent, key] if parent else float("nan")
        frame[key + ".rank"] = float("nan")
        frame[key + ".rank_within_n_mutations"] = float("nan")
        if spec["better"] != "descriptive":
            available = nonref & frame[key + ".status"].eq("available")
            values = frame.loc[available, key]
            frame.loc[available, key + ".rank"] = values.rank(ascending=spec["better"] == "lower", method="min")
            frame.loc[available, key + ".rank_within_n_mutations"] = frame.loc[available].groupby("n_mutations")[key].rank(
                ascending=spec["better"] == "lower", method="min")
    return frame.reset_index(drop=True), metrics


def correlations(frame, metrics):
    rows = []
    candidates = frame[~frame["is_reference"]]
    groups = [("all_candidates", candidates)] + [(f"n_mutations={n}", f) for n, f in candidates.groupby("n_mutations")]
    directional = [key for key, m in metrics.items() if m["better"] != "descriptive"]
    for label, group in groups:
        for left, right in combinations(directional, 2):
            pair = group.loc[group[left + ".status"].eq("available") & group[right + ".status"].eq("available"), [left, right]].dropna()
            rho = None
            if len(pair) >= 3 and all(pair[c].nunique() > 1 for c in pair):
                ranks = pair.rank()
                for key in pair:
                    if metrics[key]["better"] == "lower":
                        ranks[key] *= -1
                rho = float(ranks[left].corr(ranks[right]))
            rows.append({"group": label, "metric_a": left, "metric_b": right, "n_common": len(pair), "spearman": rho})
    return pd.DataFrame(rows, columns=["group", "metric_a", "metric_b", "n_common", "spearman"])


def select_shortlist(frame, config, metrics):
    """Deterministic, explicitly configured lexicographic selection with diversity."""
    selection = config.get("shortlist")
    if not selection:
        return frame.iloc[0:0].copy(), "Not requested; metric priorities and shortlist size need team agreement."
    statuses = ["ELIGIBLE"] + (["REVIEW"] if selection.get("include_review", False) else [])
    pool = frame[frame["decision"].isin(statuses)]
    keys = selection["sort_by"]
    for key in keys:
        pool = pool[pool[key + ".status"].eq("available")]
    pool = pool.sort_values(keys + ["design"], ascending=[metrics[k]["better"] == "lower" for k in keys] + [True])
    chosen = []
    distance = selection.get("min_sequence_distance", 1)
    for index, row in pool.iterrows():
        if all(len(row["sequence"]) == len(frame.at[i, "sequence"]) and
               sum(a != b for a, b in zip(row["sequence"], frame.at[i, "sequence"])) >= distance for i in chosen):
            chosen.append(index)
            if len(chosen) == selection["size"]:
                break
    result = frame.loc[chosen].copy()
    result.insert(0, "selection_order", range(1, len(result) + 1))
    return result, f"Selected {len(chosen)}/{selection['size']}; greedy selection, minimum Hamming distance {distance}. No quota is silently relaxed."


def provenance(path):
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None}


def run(candidate_path: Path, config_path: Path, out: Path, developability_path: Path | None = None):
    candidates = read_candidates(candidate_path)
    config = json.loads(config_path.read_text())
    validate_config(config)
    parent = config.get("parent_id")
    if parent is not None and parent not in candidates:
        raise ValueError("Configured parent_id is absent from candidate manifest")
    observations, sources = [], []
    inputs = [candidate_path, config_path, Path(__file__)] + ([developability_path] if developability_path else [])
    for source in config.get("sources", []):
        obs, status = load_source(source, candidates, config_path.parent)
        observations.extend(obs)
        sources.append(status)
        inputs.append(config_path.parent / source["path"])
        if source.get("identity_file"):
            inputs.append(config_path.parent / source["identity_file"])
    dev = load_developability(developability_path, candidates, parent)
    frame, metrics = compare(candidates, config, observations, dev)
    shortlist, selection_note = select_shortlist(frame, config, metrics)
    products = {"combined.csv": frame, "eligible.csv": frame[frame["decision"] == "ELIGIBLE"],
                "shortlist.csv": shortlist, "model_agreement.csv": correlations(frame, metrics),
                "observations.csv": pd.DataFrame(observations, columns=["design", "sequence_sha256", "metric", "replicate", "value", "protocol"])}
    destinations = [out / name for name in [*products, "summary.json", "report.md"]]
    if any(p.resolve() == source.resolve() for p in destinations for source in inputs):
        raise ValueError("An output would overwrite an input file")
    summary = {"schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
               "counts": dict(Counter(frame.loc[~frame["is_reference"], "decision"])),
               "sources": sources, "metrics": metrics, "selection": selection_note,
               "config": config, "inputs": [provenance(p) for p in inputs],
               "limitations": ["No weighted combined score; ranks are per metric and cohort, not calibrated affinity.",
                               "Agreement between models is not experimental validation; missing values are never zero.",
                               "Median aggregation is within the declared source protocol only.",
                               "Protocol metadata are declared by the producer, not independently verified.",
                               "Full-candidate scores must not be replaced by sums of single-mutation scores."]}
    lines = ["# Candidate comparison", "", "Local reporting only: no models were run and nothing was submitted.", "",
             "## Decision counts", ""]
    lines += [f"- {key}: {value}" for key, value in summary["counts"].items()]
    lines += ["", "ELIGIBLE means no configured gate or supplied warning blocked the candidate; it is not experimental validation.",
              "", "## Shortlist", "", selection_note, "", "## Missing source files", ""]
    lines += [f"- {s['source']}: {s['path']}" for s in sources if s["status"] == "file_missing"] or ["None."]
    lines += ["", "See combined.csv for every candidate, decision reasons, score coverage and parent deltas.",
              "See observations.csv for individual seeds; model_agreement.csv includes overlap counts and mutation-count strata.",
              "", "## Comparison (first 50 manifest rows; not a ranking)", ""]
    columns = ["design", "decision", *metrics, "reject_reasons", "pending_reasons", "review_reasons"]
    lines += ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for _, row in frame.head(50).iterrows():
        cells = ["unavailable" if pd.isna(row[key]) else str(row[key]).replace("|", "\\|").replace("\n", " ") for key in columns]
        lines.append("| " + " | ".join(cells) + " |")
    lines += ["", "## Limitations", ""] + [f"- {s}" for s in summary["limitations"]]
    out.mkdir(parents=True, exist_ok=True)
    for name, table in products.items():
        table.to_csv(out / name, index=False)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    (out / "report.md").write_text("\n".join(lines) + "\n")
    return frame, summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--developability", type=Path)
    parser.add_argument("--out", type=Path, default=Path("runs/candidate_comparison"))
    args = parser.parse_args(argv)
    try:
        frame, summary = run(args.candidates, args.config, args.out, args.developability)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        parser.error(str(exc))
    print(json.dumps(summary["counts"], indent=2))
    print(summary["selection"])
    print(f"Wrote {len(frame)} rows to {args.out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
