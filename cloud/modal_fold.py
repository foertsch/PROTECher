"""Fold a batch of scFv variants with Boltz-2 on Modal, using cached alignments. A fallback for the cluster GPUs.

Same Boltz (2.2.1) and torch (2.15.0) as the FGCZ node, one process for the whole batch (msa_cache.fold_batch),
no alignment server. Run everything from the repository root. cloud/README.md has the checklist.

    modal run cloud/modal_fold.py::prepare                  # once: builds the image, downloads the weights (CPU only)
    modal run cloud/modal_fold.py::smoke                    # once: five sequences, compared with the cluster's numbers
    modal run cloud/modal_fold.py::fold_dir --designs DIR --cache CACHE_DIR --out OUT [--seed 0] [--order a,b,c]
    modal run cloud/modal_fold.py::fetch --fasta PARENT.fa --out CACHE_DIR   # alignments for a new antibody, CPU only

DIR holds one FASTA per design (`>id` then `chain1:chain2`), CACHE_DIR is a folder from `msa_cache` (for example
data/msa_candidates/trastuzumab_her2). OUT receives what the cluster path writes: a DONE manifest, PDBs and scores,
so run_report and scoring read it unchanged. PROTECHER_MODAL_GPU picks the card (default L40S, the one tested).
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import time
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))  # `modal run` does not put the repository root on the path

import modal  # noqa: E402

from msa_cache.remote import localize_manifest, read_chains, stage  # noqa: E402

GPU = os.environ.get("PROTECHER_MODAL_GPU", "L40S")
BOLTZ_VERSION = "2.2.1"  # the versions installed on the cluster GPU node, so results stay comparable
TORCH_VERSION = "2.15.0"

app = modal.App("protecher-fold")
weights = modal.Volume.from_name("protecher-boltz-weights", create_if_missing=True)
work = modal.Volume.from_name("protecher-work", create_if_missing=True)

# Two environments, as on the cluster: Boltz needs numpy 1.x, the pipeline's biotite needs numpy 2.
# The pipeline runs in the image's own Python and starts `boltz` from /opt/boltz.
image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("build-essential")  # triton compiles small C stubs when the first kernel runs
    .pip_install("numpy>=2.1", "biotite>=2.0")
    .run_commands(
        "python -m venv /opt/boltz",
        f"/opt/boltz/bin/pip install --no-cache-dir 'boltz[cuda]=={BOLTZ_VERSION}' 'torch=={TORCH_VERSION}'",
    )
    .add_local_python_source("design_pipeline", "msa_cache")
)

CLUSTER_SEED0 = {  # run msa_batch_01 on the cluster GPU node: batch of five in this order, seed 0, ipTM and pLDDT
    "wt": (0.937, 91.88),
    "v1": (0.929, 92.20),
    "v2": (0.854, 92.12),
    "v3": (0.924, 91.33),
    "v4": (0.927, 91.84),
}


@app.function(image=image, volumes={"/cache": weights}, timeout=1800)
def download_weights() -> str:
    """Boltz-2's weights and molecule data into a volume, about 7.6 GB. CPU only."""
    import subprocess

    code = (
        "from pathlib import Path; from boltz.main import download_boltz2; "
        "p = Path('/cache/boltz'); p.mkdir(parents=True, exist_ok=True); download_boltz2(p)"
    )
    subprocess.run(["/opt/boltz/bin/python", "-c", code], check=True)
    weights.commit()
    return subprocess.run(["du", "-sh", "/cache/boltz"], capture_output=True, text=True, check=True).stdout.strip()


@app.function(image=image, gpu=GPU, cpu=4, memory=32768, volumes={"/cache": weights, "/work": work}, timeout=3600)
def fold(run_id: str, seed: int = 0) -> dict:
    """Fold /work/<run_id>/designs with the cache in /work/<run_id>/cache, writing /work/<run_id>/out."""
    import subprocess

    os.environ["PATH"] += ":/opt/boltz/bin"
    os.environ["BOLTZ_CACHE"] = "/cache/boltz"
    weights.reload()
    work.reload()
    from msa_cache.remote import run_staged

    gpu_name = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"], capture_output=True, text=True).stdout.strip()
    started = time.perf_counter()
    folded = run_staged(Path("/work") / run_id, seed=seed)
    seconds = time.perf_counter() - started
    work.commit()
    return {
        "gpu": gpu_name,
        "seconds": round(seconds, 1),
        "designs": [{"id": d.parent_id, "scores": dict(d.scores)} for d in folded],
    }


@app.function(image=image, volumes={"/work": work}, cpu=2, memory=8192, timeout=1800)
def fetch_remote(run_id: str, chains: list[str]) -> dict:
    """Ask the public alignment server for a complex's alignments from a CPU container and keep them as an MsaCache."""
    from msa_cache import MsaCache, fetch_alignments

    cache = MsaCache(Path("/work") / run_id / "fetched")
    fetch_alignments(cache, chains, Path("/work") / run_id / "fetch_work", python="/opt/boltz/bin/python")
    work.commit()
    return cache.info()


def _download(prefix: str, out_dir: Path) -> int:
    """Copy every file under `prefix` in the work volume into out_dir, keeping the structure below it."""
    count = 0
    for entry in work.iterdir(prefix, recursive=True):
        if entry.type != modal.volume.FileEntryType.FILE:
            continue
        target = out_dir / PurePosixPath(entry.path.strip("/")).relative_to(prefix)
        target.parent.mkdir(parents=True, exist_ok=True)
        with open(target, "wb") as fh:
            for chunk in work.read_file(entry.path):
                fh.write(chunk)
        count += 1
    return count


def _run(designs: Path, cache: Path, out: Path, seed: int, order: list[str] | None) -> dict:
    run_id = time.strftime("run_%Y%m%d_%H%M%S")
    with tempfile.TemporaryDirectory() as tmp:
        ids = stage(designs, cache, Path(tmp) / run_id, order=order)
        print(f"staged {len(ids)} designs: {', '.join(ids)}")
        with work.batch_upload(force=True) as batch:
            batch.put_directory(str(Path(tmp) / run_id), f"/{run_id}")
    started = time.perf_counter()
    result = fold.remote(run_id, seed)
    elapsed = time.perf_counter() - started
    files = _download(f"{run_id}/out", out)
    localize_manifest(out, f"/work/{run_id}/out")
    print(f"folded on {result['gpu']}: {result['seconds']} s of folding, {elapsed:.0f} s including start-up; {files} files in {out}")
    result["elapsed"] = elapsed
    return result


@app.local_entrypoint()
def prepare() -> None:
    """Build the image and fill the weights volume. Run once, a few minutes of CPU, no GPU."""
    print("weights volume now holds:", download_weights.remote())


@app.local_entrypoint()
def fetch(fasta: str, out: str) -> None:
    """Fetch alignments for a new parent (one record, chains separated by ':') and write them as an msa_cache folder."""
    chains = read_chains(Path(fasta))
    run_id = time.strftime("fetch_%Y%m%d_%H%M%S")
    info = fetch_remote.remote(run_id, chains)
    files = _download(f"{run_id}/fetched", Path(out))
    print(f"{files} files in {out}; chain lengths {info['chain_lengths']}, alignment rows {info['alignment_rows']}")


@app.local_entrypoint()
def fold_dir(designs: str, cache: str, out: str, seed: int = 0, order: str = "") -> None:
    """Fold every *.fa in `designs` against the cache in `cache`, results into `out`."""
    result = _run(Path(designs), Path(cache), Path(out), seed, [i for i in order.split(",") if i] or None)
    for d in result["designs"]:
        print(f"  {d['id']:>12s}  ipTM {d['scores'].get('iptm', float('nan')):.3f}  pLDDT {d['scores'].get('plddt', float('nan')):.2f}")


@app.local_entrypoint()
def smoke() -> None:
    """Fold the five trastuzumab sequences of the cluster test and print them next to the cluster's numbers."""
    data = ROOT / "data" / "modal_smoke"
    if not (data / "cache" / "parent.json").is_file():
        raise SystemExit(f"{data} is missing; it holds the five test sequences and their alignment cache (5 MB, copied from the cluster)")
    out = Path(tempfile.mkdtemp(prefix="modal_smoke_"))
    result = _run(data / "designs", data / "cache", out, seed=0, order=list(CLUSTER_SEED0))
    print("\n  design   Modal ipTM  cluster ipTM   diff   Modal pLDDT  cluster pLDDT")
    diffs = []
    for d in result["designs"]:
        ref_iptm, ref_plddt = CLUSTER_SEED0[d["id"]]
        iptm, plddt = d["scores"]["iptm"], d["scores"]["plddt"]
        diffs.append(abs(iptm - ref_iptm))
        print(f"  {d['id']:6s}   {iptm:9.3f}  {ref_iptm:12.3f}  {iptm - ref_iptm:+6.3f}   {plddt:10.2f}  {ref_plddt:13.2f}")
    print(f"\n  mean |ipTM difference| {sum(diffs) / len(diffs):.3f} (the cluster's own seed-to-seed spread is about 0.01 to 0.02 per design)")
    print(f"  results in {out}")
