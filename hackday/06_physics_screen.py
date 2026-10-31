#!/usr/bin/env python3
"""CPU structure screening. Only explicitly mapped structures receive metrics.

Default: the parent Fv and the available HER2 480-607 crystal coordinates.
Use --input for candidate sequences and --manifest for their actual structures.
See PHYSICS_SCREENING.md. No GPU, inference, atom rebuilding or energy calculation.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

from biotite import InvalidFileError

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from design_pipeline.physics import (EVENT_FIELDS, METRICS, PARAMETERS, analyze_geometry,
                                     load_structure, select_structure, sequence_hash)

spec = importlib.util.spec_from_file_location("developability", Path(__file__).with_name("05_developability.py"))
dev = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dev)

LIMITATIONS = [
    "CPU geometry screen only; no force-field energy, ddG, expression, affinity or melting-temperature prediction.",
    "Thresholds are explicit review windows, not calibrated biological rejection rules.",
    "Heavy-atom overlaps exclude 1-2/1-3/1-4 neighbours and plausible disulfides; not MolProbity clashscore.",
    "Asp/Glu assumed negative and Lys/Arg positive. Histidine, termini, water, ions and glycans are not scored as charged groups.",
    "Polar distance contacts are not hydrogen bonds. Hydrogen-bond geometry requires explicit hydrogens; protonation is not optimized.",
    "SASA uses heavy atoms and single-element radii. Hydrophobic side-chain area is descriptive, not an aggregation probability.",
    "No candidate inherits parent coordinates or parent scores. Missing structures stay unassessed.",
    "Fv-only structures omit the scFv linker and Fab constant domains. Target cuts and missing loops affect interface and surface metrics.",
    "No automatic candidate deletion, combined score or final selection.",
]


def read_manifest(path, sequences):
    manifest = json.loads(path.read_text())
    if not isinstance(manifest, dict) or set(manifest) != {"schema_version", "structures"} or manifest["schema_version"] != 1:
        raise ValueError("Manifest needs schema_version=1 and structures")
    if not isinstance(manifest["structures"], list):
        raise ValueError("structures must be a list")
    records = {}
    required = {"design", "sequence_sha256", "path", "scope", "binder_segments"}
    for record in manifest["structures"]:
        if not isinstance(record, dict) or not required <= set(record) or set(record) - required - {"antigen_segments", "model"}:
            raise ValueError("Structure record fields do not match the documented schema")
        design = record["design"]
        if design not in sequences or design in records:
            raise ValueError(f"Unknown or duplicate structure design ID: {design}")
        if record["sequence_sha256"] != sequence_hash(sequences[design]):
            raise ValueError(f"Manifest/sequence hash mismatch for {design}")
        if not isinstance(record["path"], str) or not record["path"]:
            raise ValueError("Structure path must be nonempty")
        records[design] = dict(record, path=(path.parent / record["path"]).resolve())
    return records


def run(input_path, manifest_path, output_dir, csv_no_header=False):
    sequences = dev.prepare_sequences(dev.read_sequences(input_path, csv_no_header=csv_no_header) if input_path else {}, "PARENT")
    records = read_manifest(manifest_path, sequences)
    if "PARENT" not in records:
        raise ValueError("Include the PARENT structure for the baseline")
    outputs = [output_dir / name for name in ("physics_review.csv", "physics_events.csv", "physics_metadata.json", "physics_report.md")]
    sources = [manifest_path, Path(__file__), ROOT / "design_pipeline/physics.py", ROOT / "hackday/05_developability.py"]
    if input_path:
        sources.append(input_path)
    sources += [r["path"] for r in records.values()]
    if any(o.resolve() == s.resolve() for o in outputs for s in sources):
        raise ValueError("Output paths cannot overwrite an input or implementation file")
    rows, all_events, details = [], [], {}
    for design, sequence in sequences.items():
        row = {"design": design, "sequence": sequence, "sequence_sha256": sequence_hash(sequence),
               "physics_status": "not_evaluated_no_structure", "physics_scope": None,
               "physics_flags": "", "physics_missing_checks": "candidate_structure_not_provided",
               "physics_error": "", "hydrogen_bond_status": "not_evaluated", "parent_comparison_status": "not_evaluated"}
        row.update({m: None for m in METRICS})
        row.update({"d_"+m: None for m in METRICS})
        if design in records:
            record = records[design]
            row["physics_scope"] = record["scope"]
            try:
                atoms = load_structure(record["path"], record.get("model", 1))
                atoms, binder, positions, coverage = select_structure(atoms, record, sequence)
                result = analyze_geometry(atoms, binder, positions)
                missing = list(result["missing_checks"])
                if coverage["antigen_missing_residues"]:
                    missing.append("antigen_coordinates_incomplete")
                row.update(result["metrics"])
                row.update(physics_status="review_required" if result["flags"] else "partial" if missing else "no_geometry_warnings",
                           physics_flags=";".join(result["flags"]), physics_missing_checks=";".join(missing),
                           hydrogen_bond_status=result["hydrogen_bond_status"])
                for event in result.pop("events"):
                    all_events.append({"design": design, **event})
                details[design] = {**coverage, **result}
            except (OSError, ValueError, IndexError, KeyError, TypeError, InvalidFileError) as exc:
                row.update(physics_status="invalid_structure", physics_error=str(exc),
                           physics_missing_checks="structure_validation_failed")
        rows.append(row)
    parent = next(r for r in rows if r["design"] == "PARENT")
    for row in rows:
        design = row["design"]
        if design not in details or "PARENT" not in details:
            continue
        if row["physics_scope"] != parent["physics_scope"]:
            row["parent_comparison_status"] = "unavailable_scope_mismatch"
            continue
        same_target = details[design]["antigen_signature"] == details["PARENT"]["antigen_signature"]
        row["parent_comparison_status"] = "same_scope_and_target_coverage" if same_target else "binder_only_target_mismatch"
        for metric in METRICS:
            is_interface = "interface" in metric or "on_binding" in metric
            incomplete = (details[design]["metrics"]["missing_heavy_atoms_binder"] or
                          details["PARENT"]["metrics"]["missing_heavy_atoms_binder"])
            if is_interface:
                incomplete |= bool(details[design]["missing_heavy_atoms_antigen"] or details["PARENT"]["missing_heavy_atoms_antigen"])
            if (is_interface and not same_target) or (incomplete and metric != "missing_heavy_atoms_binder"):
                continue
            if "hydrogen_bonds" in metric and (row["hydrogen_bond_status"] != "available_explicit_hydrogens"
                                               or parent["hydrogen_bond_status"] != "available_explicit_hydrogens"):
                continue
            if row[metric] is not None and parent[metric] is not None:
                row["d_"+metric] = row[metric]-parent[metric]
    metadata = {"schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
                "parameters": PARAMETERS, "limitations": LIMITATIONS,
                "versions": {name: version(name) for name in ("numpy", "biotite")},
                "inputs": [dev.file_info(p) for p in sources], "structure_details": details,
                "structure_manifest": json.loads(manifest_path.read_text())}
    output_dir.mkdir(parents=True, exist_ok=True)
    for path, data, fields in ((outputs[0], rows, list(rows[0])), (outputs[1], all_events, ["design", *EVENT_FIELDS])):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(data)
    outputs[2].write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n")
    lines = ["# CPU physics and geometry review", "", "Geometry evidence only; no experimental pass/fail or energy prediction.", "",
             "| Design | Status | Scope | Geometry flags | Missing checks |", "|---|---|---|---|---|"]
    for row in rows:
        cells = [row[k] for k in ("design", "physics_status", "physics_scope", "physics_flags", "physics_missing_checks")]
        lines.append("| " + " | ".join(str(v or "—").replace("|", "\\|").replace("\n", " ") for v in cells) + " |")
    lines += ["", "## Parent baseline", ""]
    for metric in METRICS:
        lines.append(f"- {metric}: {parent[metric] if parent[metric] is not None else 'unavailable'}")
    lines += ["", "## Interpretation", "", *["- "+note for note in LIMITATIONS], "",
              "See physics_events.csv for individual atom pairs, distances and reasons; physics_metadata.json records coverage and provenance."]
    outputs[3].write_text("\n".join(lines)+"\n")
    return rows, metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="Candidate FASTA or name,sequence CSV; omitted means parent only")
    parser.add_argument("--csv-no-header", action="store_true")
    parser.add_argument("--manifest", default="hackday/physics_parent.json")
    parser.add_argument("--out-dir", default="runs/physics_parent")
    args = parser.parse_args(argv)
    try:
        rows, _ = run(ROOT / args.input if args.input else None, ROOT / args.manifest, ROOT / args.out_dir, args.csv_no_header)
    except (OSError, ValueError, TypeError, csv.Error) as exc:
        parser.error(str(exc))
    for row in rows:
        print(f"{row['design']}: {row['physics_status']} {row['physics_error']}")
    print(f"Reports: {ROOT / args.out_dir}. Blank metrics are unavailable, not zero.")
    return int(any(row["physics_status"] == "invalid_structure" for row in rows))


if __name__ == "__main__":
    raise SystemExit(main())
