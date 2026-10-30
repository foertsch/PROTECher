# Modal fallback for Boltz-2 folds

For the case that the cluster GPUs are busy. It folds a batch of scFv variants with Boltz-2 on a Modal GPU, using the
cached alignments, and writes the same output folder as the cluster path (a `DONE` manifest, PDBs, scores).

**Status.** Not run on Modal yet, because that needs your account. The staging and result-handling code is tested
(`tests/test_msa_remote.py`, 172 tests pass on the cluster) and the script loads and registers its functions. The image
build, the weights download and the GPU run are untested. The steps below find out tonight.

## Tonight

Run everything from the repository root. You sign up and log in; nothing here creates an account.

```bash
uv tool install modal
modal setup
```

```bash
modal run cloud/modal_fold.py::prepare
```
Builds the image (Boltz 2.2.1 and torch 2.15.0, the versions on the cluster) and downloads 7.6 GB of weights into the
volume `protecher-boltz-weights`. CPU only. The first build takes several minutes.

```bash
modal run cloud/modal_fold.py::smoke
```
Folds the five trastuzumab test sequences, in the cluster's order and with seed 0, and prints them next to the cluster's
numbers (run `msa_batch_01`). It prints which GPU ran. The first design should land close to the cluster's 0.937. Folding
in a different place changes the random draws slightly, so expect small differences. The cluster's own seed-to-seed
spread is about 0.01 to 0.02 per design, and its v2 value of 0.854 was a single low draw.

```bash
modal run cloud/modal_fold.py::fetch --fasta data/modal_smoke/designs/wt.fa --out /tmp/fetched_check
```
Optional. Fetches alignments for a parent from a CPU container, which is the way to get them for a new antibody if the
cluster is not reachable. For this parent it should report chain lengths `[242, 128]` and alignment rows near
`[16383, 7777]`.

## Tomorrow, if the cluster GPUs are busy

This is the only Modal script in the repository, for every fold that has to leave the cluster.

```bash
modal run cloud/modal_fold.py::fold_dir --designs DESIGNS_DIR --cache data/msa_candidates/trastuzumab_her2 --out runs/modal_01
```

- `DESIGNS_DIR` holds one FASTA per design: `>id`, then `chain1:chain2`, as the pipeline writes them.
- `--cache` is a folder from `msa_cache`. The six candidates are in `data/msa_candidates/`. For an antibody that is not
  one of them, run `fetch` on its parent first.
- A design that is too far from the parent's chains is refused before anything is uploaded, with its name.
- The order of the designs changes each one's random draw. Pass `--order a,b,c` to fix it, otherwise it is alphabetical.
- `--seeds N` folds seeds 0 to N-1, each in its own GPU container at the same time, so five seeds take about as long as
  one and cost five times as much. Results go into `OUT/seed<k>/`, and the ipTM mean and SD per design are printed.
  With the default `--seeds 1` the results go straight into `OUT`.
- `PROTECHER_MODAL_GPU=H100` (or `A100-80GB`) picks another card. The default, L40S, is the one tested on the cluster.

### Track A: the veto folds, several seeds, trimmed antigen

The plan's structure veto folds the parent and the finalists against a trimmed antigen over at least five seeds.
`examples/track_a_fold.py --dry-run` trims the antigen and writes one two-chain FASTA per scFv, which is the form
`fold_dir` reads. The alignments come from the parent, once:

```bash
python examples/track_a_fold.py --scfv scfv.fasta --antigen antigen.fasta --antigen-region 480-607 --dry-run --out runs/veto
modal run cloud/modal_fold.py::fetch --fasta runs/veto/inputs/PARENT.fa --out runs/veto/cache
modal run cloud/modal_fold.py::fold_dir --designs runs/veto/inputs --cache runs/veto/cache --out runs/veto/modal --seeds 5
```

Replace `PARENT.fa` with the parent's file in `runs/veto/inputs/`. A variant whose scFv has changed too much from the
parent's, or a different antigen fragment, is refused before the upload.

## Costs and cleanup

You pay for GPU seconds while a fold runs. A third-party price list gives about $1.95 per hour for an L40S; check
`modal.com/pricing`. The smoke test is a few minutes. The two volumes keep their files until deleted:

```bash
modal volume delete protecher-work
modal volume delete protecher-boltz-weights
```

## Limits

- Boltz-2 only. Any model Kira swaps in for the antibody work needs its own image.
- Results are not bit-identical to the cluster's, because of different hardware and drivers.
- Sequences and alignments are uploaded to Modal. The candidate antibodies and antigens are public; decide separately
  about the team's own designs.
