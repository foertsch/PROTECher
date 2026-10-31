"""Fold a batch of scFv variants with Boltz-2 on Modal, using cached alignments. A fallback for the cluster GPUs.

Same Boltz (2.2.1) and torch (2.15.0) as the FGCZ node, one process for the whole batch (msa_cache.fold_batch),
no alignment server. This is the one Modal script of the repository. Run everything from the repository root.
cloud/README.md has the checklist.

    modal run cloud/modal_fold.py::prepare                  # once: builds the image, fills the weights volume if empty (CPU only)
    modal run cloud/modal_fold.py::smoke                    # once: five sequences, compared with the cluster's numbers
    modal run cloud/modal_fold.py::fold_dir --designs DIR --cache CACHE_DIR --out OUT [--seeds 1] [--order a,b,c]
    modal run cloud/modal_fold.py::fetch --fasta PARENT.fa --out CACHE_DIR   # alignments for a new antibody, CPU only

DIR holds one FASTA per design (`>id` then `chain1:chain2`), CACHE_DIR is a folder from `msa_cache` (for example
data/msa_candidates/trastuzumab_her2). OUT receives what the cluster path writes: a DONE manifest, PDBs and scores,
so run_report and scoring read it unchanged. PROTECHER_MODAL_GPU picks the card (default L40S, the one tested).

`--seeds N` folds seeds 0 to N-1, each in its own GPU container at the same time, so five seeds take about as long
as one. The results go into `OUT/seed<k>/` (one folder with `--seeds 1` goes straight into OUT), and the ipTM mean
and SD per design are printed: the seed noise the structure veto of the Track A plan is measured against.

Track A's veto folds, with the antigen trimmed (examples/track_a_fold.py writes the inputs and stops with --dry-run):

    python examples/track_a_fold.py --scfv scfv.fasta --antigen antigen.fasta --antigen-region 480-607 --dry-run --out runs/veto
    modal run cloud/modal_fold.py::fetch --fasta runs/veto/inputs/PARENT.fa --out runs/veto/cache
    modal run cloud/modal_fold.py::fold_dir --designs runs/veto/inputs --cache runs/veto/cache --out runs/veto/modal --seeds 5
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import time
import uuid
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))  # `modal run` does not put the repository root on the path

import modal  # noqa: E402

from msa_cache.remote import localize_manifest, read_chains, seed_spread, stage  # noqa: E402

GPU = os.environ.get("PROTECHER_MODAL_GPU", "L40S")
BOLTZ_VERSION = "2.2.1"  # the versions installed on the cluster GPU node, so results stay comparable
TORCH_VERSION = "2.15.0"

app = modal.App("protecher-fold")
weights = modal.Volume.from_name("boltz-weights", create_if_missing=True)  # the name Track A's earlier script already filled
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
    """Boltz-2's weights and molecule data into the volume's root, about 7.6 GB, skipping what is already there. CPU only."""
    import subprocess

    code = (
        "from pathlib import Path; from boltz.main import download_boltz2; "
        "p = Path('/cache'); p.mkdir(parents=True, exist_ok=True); download_boltz2(p)"
    )
    subprocess.run(["/opt/boltz/bin/python", "-c", code], check=True)
    weights.commit()
    return subprocess.run(["du", "-sh", "/cache"], capture_output=True, text=True, check=True).stdout.strip()


@app.function(image=image, gpu=GPU, cpu=4, memory=32768, volumes={"/cache": weights, "/work": work}, timeout=3600)
def fold(run_id: str, seed: int = 0) -> dict:
    """Fold /work/<run_id>/designs with the cache in /work/<run_id>/cache, writing /work/<run_id>/out_s<seed>.

    One folder per seed, because several seeds run in separate containers at the same time and must not share one.
    """
    import subprocess

    os.environ["PATH"] += ":/opt/boltz/bin"
    os.environ["BOLTZ_CACHE"] = "/cache"
    weights.reload()
    work.reload()
    from msa_cache.remote import run_staged

    gpu_name = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"], capture_output=True, text=True).stdout.strip()
    started = time.perf_counter()
    folded = run_staged(Path("/work") / run_id, seed=seed, out_name=f"out_s{seed}")
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


def _run(designs: Path, cache: Path, out: Path, seeds: list[int], order: list[str] | None) -> list[dict]:
    """Stage once, fold every seed in its own container at the same time, download each into out (one seed) or out/seed<k>."""
    run_id = f"{time.strftime('run_%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:4]}"  # two runs in one second must not share a folder
    with tempfile.TemporaryDirectory() as tmp:
        ids = stage(designs, cache, Path(tmp) / run_id, order=order)
        print(f"staged {len(ids)} designs: {', '.join(ids)}")
        with work.batch_upload(force=True) as batch:
            batch.put_directory(str(Path(tmp) / run_id), f"/{run_id}")
    started = time.perf_counter()
    results = list(fold.starmap([(run_id, seed) for seed in seeds]))
    elapsed = time.perf_counter() - started
    for seed, result in zip(seeds, results):
        target = out if len(seeds) == 1 else out / f"seed{seed}"
        files = _download(f"{run_id}/out_s{seed}", target)
        localize_manifest(target, f"/work/{run_id}/out_s{seed}")
        print(f"seed {seed}: folded on {result['gpu']}, {result['seconds']} s of folding; {files} files in {target}")
        result["seed"] = seed
    print(f"{len(seeds)} seed(s) done in {elapsed:.0f} s including start-up")
    return results


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
def fold_dir(designs: str, cache: str, out: str, seeds: int = 1, order: str = "") -> None:
    """Fold every *.fa in `designs` against the cache in `cache` for seeds 0 to seeds-1, results into `out`."""
    if seeds < 1:
        raise SystemExit("--seeds must be at least 1")
    results = _run(Path(designs), Path(cache), Path(out), list(range(seeds)), [i for i in order.split(",") if i] or None)
    if seeds == 1:
        for d in results[0]["designs"]:
            print(f"  {d['id']:>12s}  ipTM {d['scores'].get('iptm', float('nan')):.3f}  pLDDT {d['scores'].get('plddt', float('nan')):.2f}")
        return
    print("\n        design  seeds  ipTM mean     SD")
    for name, n, mean, sd in seed_spread(results):
        print(f"  {name:>12s}  {n:5d}  {mean:9.3f}  {sd:.3f}" if sd is not None else f"  {name:>12s}  {n:5d}  {mean:9.3f}    n/a")


@app.local_entrypoint()
def smoke() -> None:
    """Fold the five trastuzumab sequences of the cluster test and print them next to the cluster's numbers."""
    data = ROOT / "data" / "modal_smoke"
    if not (data / "cache" / "parent.json").is_file():
        raise SystemExit(f"{data} is missing; it holds the five test sequences and their alignment cache (5 MB, copied from the cluster)")
    out = Path(tempfile.mkdtemp(prefix="modal_smoke_"))
    result = _run(data / "designs", data / "cache", out, seeds=[0], order=list(CLUSTER_SEED0))[0]
    print("\n  design   Modal ipTM  cluster ipTM   diff   Modal pLDDT  cluster pLDDT")
    diffs = []
    for d in result["designs"]:
        ref_iptm, ref_plddt = CLUSTER_SEED0[d["id"]]
        iptm, plddt = d["scores"]["iptm"], d["scores"]["plddt"]
        diffs.append(abs(iptm - ref_iptm))
        print(f"  {d['id']:6s}   {iptm:9.3f}  {ref_iptm:12.3f}  {iptm - ref_iptm:+6.3f}   {plddt:10.2f}  {ref_plddt:13.2f}")
    print(f"\n  mean |ipTM difference| {sum(diffs) / len(diffs):.3f} (the cluster's own seed-to-seed spread is about 0.01 to 0.02 per design)")
    print(f"  results in {out}")
