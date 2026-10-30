"""Boltz-2 on Modal: the cloud fallback for examples/track_a_fold.py.

Same inputs, same trimming, same seeds, but each fold runs in a Modal container on an L40S
instead of the cluster. Use it if the cluster is down or full. The weights live in a Modal Volume,
so the multi-GB download happens once and later runs start from the cache. Seeds and complexes run
in parallel, up to the plan's GPU concurrency.

    modal run examples/modal_boltz.py --scfv scfv.fasta --antigen antigen.fasta \\
        --antigen-region 480-607 --seeds 5 --out runs/track_a_modal

    modal run examples/modal_boltz.py --scfv scfv.fasta --antigen antigen.fasta --dry-run

`--dry-run` prints the plan and starts nothing. The scFv is chain A and the antigen chain B, as in
track_a_fold.py. Results land in `<out>/predictions/<name>_s<seed>/` (the mmCIF, Boltz's confidence
JSON and the PAE matrix) with a `summary.csv` of ipTM, pTM and pLDDT per fold.

Only public sequences belong here: the cluster, FoldX and anything internal stay off Modal, and the
sequences reach the public MSA server from the container, as they do from the cluster.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import modal

HERE = Path(__file__).resolve().parent
CACHE = "/cache"
GPU = "L40S"

app = modal.App("track-a-boltz")
weights = modal.Volume.from_name("boltz-weights", create_if_missing=True)
# [cuda] is Boltz's own GPU extra, as in pyproject.toml. Without it the triangular kernels fail with
# "No module named cuequivariance_torch".
image = modal.Image.debian_slim(python_version="3.11").pip_install("boltz[cuda]>=2.2,<3")


def boltz_fasta(chains: list[str]) -> str:
    """One record per chain in Boltz's `chain_id|entity_type` header form: A, B, ..."""
    return "".join(f">{chr(ord('A') + i)}|protein\n{seq}\n" for i, seq in enumerate(chains))


def confidence_row(name: str, seed: int, confidence_json: str) -> dict:
    """ipTM, pTM and pLDDT (0 to 100) from Boltz's confidence file, as the cluster adapter reports them."""
    payload = json.loads(confidence_json)
    plddt = payload.get("complex_plddt")
    return {
        "complex": name,
        "seed": seed,
        "iptm": payload.get("iptm"),
        "ptm": payload.get("ptm"),
        "plddt": None if plddt is None else 100.0 * float(plddt),
    }


@app.function(image=image, gpu=GPU, volumes={CACHE: weights}, timeout=1800)
def fold(name: str, chains: list[str], seed: int, msa_server: bool = True) -> dict:
    """Fold one complex with one seed and return the files as text and bytes."""
    import subprocess
    import tempfile

    work = Path(tempfile.mkdtemp())
    fasta = work / f"{name}.fasta"
    fasta.write_text(boltz_fasta(chains))
    command = [
        "boltz", "predict", str(fasta), "--out_dir", str(work / "out"), "--cache", CACHE,
        "--output_format", "mmcif", "--write_full_pae", "--seed", str(seed),
    ]
    if msa_server:
        command.append("--use_msa_server")
    try:
        subprocess.run(command, check=True)
    finally:
        weights.commit()  # keep the downloaded weights even if the fold failed

    def first(pattern: str) -> Path:
        found = sorted((work / "out").rglob(pattern))
        if not found:
            raise RuntimeError(f"boltz wrote nothing matching {pattern}")
        return found[0]

    pae = list((work / "out").rglob("pae_*_model_0.npz"))
    return {
        "name": name,
        "seed": seed,
        "cif": first("*_model_0.cif").read_text(),
        "confidence": first("confidence_*_model_0.json").read_text(),
        "pae": pae[0].read_bytes() if pae else b"",
    }


@app.local_entrypoint()
def main(
    scfv: str,
    antigen: str,
    antigen_region: str = "",
    epitope: str = "",
    pad: int = 20,
    seeds: int = 5,
    out: str = "runs/track_a_modal",
    dry_run: bool = False,
):
    # Imported here: the container only needs boltz, and does not have the examples folder.
    sys.path.insert(0, str(HERE))
    import track_a_fold as inputs

    scfvs = inputs.read_sequences(scfv)
    antigens = inputs.read_sequences(antigen)
    if len(antigens) != 1:
        raise SystemExit(f"--antigen holds {len(antigens)} records; give exactly one")
    if seeds < 1:
        raise SystemExit("--seeds must be at least 1")
    region = inputs.parse_region(antigen_region) if antigen_region else None
    residues = inputs.parse_residues(epitope) if epitope else None
    trimmed, start, end = inputs.trim_antigen(antigens[0][1], region, residues, pad)
    for message in inputs.warnings_for(scfvs, len(trimmed)):
        print(f"warning: {message}", file=sys.stderr)

    jobs = [(name, [seq, trimmed], seed) for name, seq in scfvs for seed in range(seeds)]
    print(
        f"{len(scfvs)} complex(es) x {seeds} seed(s) = {len(jobs)} fold(s), antigen residues "
        f"{start}-{end} of {len(antigens[0][1])} ({len(trimmed)} kept), {GPU} on Modal"
    )
    if dry_run:
        print("dry run: nothing started")
        return

    target = Path(out)
    rows = []
    for result in fold.starmap(jobs):
        folder = target / "predictions" / f"{result['name']}_s{result['seed']}"
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "model_0.cif").write_text(result["cif"])
        (folder / "confidence.json").write_text(result["confidence"])
        if result["pae"]:
            (folder / "pae.npz").write_bytes(result["pae"])
        row = confidence_row(result["name"], result["seed"], result["confidence"])
        row["cif"] = str(folder / "model_0.cif")
        rows.append(row)

    with (target / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["complex", "seed", "iptm", "ptm", "plddt", "cif"])
        writer.writeheader()
        writer.writerows(rows)
    for name, n, mean, sd in inputs.seed_spread(rows):
        spread = f"{sd:.3f}" if sd is not None else "n/a (one seed)"
        print(f"{name}: ipTM mean {mean:.3f}, SD {spread} over {n} seed(s)")
    print(f"wrote {target / 'summary.csv'}")
