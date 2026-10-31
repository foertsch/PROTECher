"""Candidate variants scored with AntiFold's own functions, in the scorer contract of design_pipeline.

Two routes, both AntiFold's own code, neither with scoring maths of ours beyond turning its averages into totals:

* `complex` (default). AntiFold's `score_sequence_in_complex`: the likelihood of a chain's sequence given the
  backbone of the whole antibody, one forward pass per mutated chain. The mutant is scored as a whole, so several
  mutations in one chain are not the sum of the singles, which is what design_pipeline asks of a scorer.
* `table`. AntiFold's per-residue table (the log-likelihood of every amino acid at every residue) and its
  `get_sequence_sampled_global_score`. One distribution per position, so the lookup is exact for a single mutation
  only; a candidate with several mutations gets no score rather than a sum.

It runs inside the AntiFold environment, which has torch and the `antifold` package:

    python -m antifold_score --pdb 1N8Z.pdb --heavy B --light A --out_dir out --mutations_csv candidates.csv \\
        --antifold_src $PROTECHER_SCRATCH/protecher_tools/antifold/src

`--pdb`, `--out_dir` and `--mutations_csv` are the flags the container scorers take, and the output is the same
`scores.csv`, so the command can go into an image later without a change to the pipeline.

What is scored: the heavy and the light chain only, no antigen. A candidate with mutations in both chains gets the
sum of the two chain scores, because the chains are scored one after the other from the backbone alone.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path
from typing import Callable, Mapping, Protocol, Sequence

import pandas as pd

AMINO = "ACDEFGHIKLMNPQRSTVWY"
SCORES = "scores.csv"
"""The name design_pipeline's scorers read back."""

Mutation = tuple[str, int, str, str]
"""(chain, position, wt, mut); position is 1-based into the chain as given, as in design_pipeline.variants."""
Row = tuple[str, float | None, float | None, str]
"""(variant_id, score, raw, note). An empty score always comes with a note saying why."""
Scorer = Callable[[str, pd.DataFrame], tuple[float, float]]


class ScoreError(RuntimeError):
    """The inputs cannot be scored at all (as opposed to one candidate that does not fit)."""


class ChainScorer(Protocol):
    """What `score_with_complex` needs: the chains of the structure and a likelihood for a chain sequence."""

    sequences: Mapping[str, str]

    def log_likelihood(self, chain: str, sequence: str) -> float: ...


def read_candidates(path: Path) -> dict[str, list[Mutation]]:
    """variant_id -> its mutations, in file order. Columns: variant_id,chain,position,wt,mut."""
    grouped: dict[str, list[Mutation]] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            grouped.setdefault(row["variant_id"], []).append((row["chain"], int(row["position"]), row["wt"], row["mut"]))
    return grouped


def check_mutations(sequences: Mapping[str, str], mutations: Sequence[Mutation]) -> str:
    """Why a candidate does not fit the structure, or an empty string when it does.

    The wild-type check is the one that catches a numbering that does not line up with the structure.
    """
    seen: set[tuple[str, int]] = set()
    for chain, position, wt, mut in mutations:
        sequence = sequences.get(chain)
        if sequence is None:
            return f"chain {chain!r} is not in the structure (chains: {sorted(sequences)})"
        if not 0 < position <= len(sequence):
            return f"position {position} is outside chain {chain}, which has {len(sequence)} residues"
        if sequence[position - 1] != wt:
            return f"wild type {wt} at {chain}{position}, the structure has {sequence[position - 1]}"
        if mut not in AMINO:
            return f"{mut!r} is not one of the 20 amino acids"
        if (chain, position) in seen:
            return f"{chain}{position} is mutated twice in one candidate"
        seen.add((chain, position))
    return ""


def _mutant_chains(sequences: Mapping[str, str], mutations: Sequence[Mutation]) -> dict[str, str]:
    """chain -> its sequence with this candidate's mutations applied, for the chains it touches."""
    edited = {chain: list(sequences[chain]) for chain, *_ in mutations}
    for chain, position, _wt, mut in mutations:
        edited[chain][position - 1] = mut
    return {chain: "".join(residues) for chain, residues in edited.items()}


class ComplexScorer:
    """AntiFold's `score_sequence_in_complex`, bound to one structure and one set of weights.

    Chains are given in the order the model reads them: heavy first, then light. `sequences` holds each chain as
    AntiFold reads it from the file, which is what a candidate's 1-based positions index.
    """

    def __init__(self, pdb: Path, chains: Sequence[str], device: str = "cpu", weights: Path | None = None):
        try:
            import antifold.esm
            from antifold.antiscripts import load_model
            from antifold.esm_multichain_util_custom import load_complex_coords, score_sequence_in_complex
        except ModuleNotFoundError as exc:
            raise ScoreError(
                "the antifold package is not importable here. Run this inside the AntiFold environment, "
                "for example $PROTECHER_SCRATCH/protecher_tools/antifold/venv/bin/python -m antifold_score"
            ) from exc
        self._score = score_sequence_in_complex
        self.model = load_model(**({"checkpoint_path": str(weights)} if weights else {}), device=device)
        self.alphabet = antifold.esm.data.Alphabet.from_architecture("invariant_gvp")
        self.coords, sequences, *_ = load_complex_coords(str(pdb), list(chains))
        self.sequences: dict[str, str] = {str(chain): str(sequence) for chain, sequence in sequences.items()}

    def log_likelihood(self, chain: str, sequence: str) -> float:
        """Average log-likelihood per residue of `sequence` as the chain, given the backbone of the complex."""
        full_chain, _without_missing_coordinates = self._score(self.model, self.alphabet, self.coords, chain, sequence)
        return float(full_chain)


def score_with_complex(scorer: ChainScorer, candidates: Mapping[str, Sequence[Mutation]]) -> list[Row]:
    """(variant_id, score, raw, note) for every candidate, from AntiFold's likelihood of each mutated chain.

    `score` is the change in total log-likelihood in nats, summed over the chains the candidate touches: for each,
    the change in AntiFold's per-residue average times the chain length. Higher is better, the convention of
    design_pipeline's scorers. `raw` is the log-likelihood of the whole mutant in nats, summed over all chains.

    A candidate that does not fit the structure gets no score and a note (see `check_mutations`).
    """
    sequences = scorer.sequences
    wild: dict[str, float] = {chain: scorer.log_likelihood(chain, sequence) for chain, sequence in sequences.items()}
    wild_total = sum(len(sequences[chain]) * ll for chain, ll in wild.items())
    seen: dict[tuple[str, str], float] = {}
    scored: list[Row] = []
    for variant_id, mutations in candidates.items():
        problem = check_mutations(sequences, mutations)
        if problem:
            scored.append((variant_id, None, None, problem))
            continue
        change = 0.0
        for chain, mutant in _mutant_chains(sequences, mutations).items():
            if (chain, mutant) not in seen:
                seen[(chain, mutant)] = scorer.log_likelihood(chain, mutant)
            change += (seen[(chain, mutant)] - wild[chain]) * len(mutant)
        scored.append((variant_id, change, wild_total + change, ""))
    return scored


def default_scorer() -> Scorer:
    """AntiFold's own `get_sequence_sampled_global_score`. Imported late so that the rest works without torch."""
    try:
        from antifold.antiscripts import get_sequence_sampled_global_score
    except ModuleNotFoundError as exc:
        raise ScoreError(
            "the antifold package is not importable here. Run this inside the AntiFold environment, "
            "for example $PROTECHER_SCRATCH/protecher_tools/antifold/venv/bin/python -m antifold_score"
        ) from exc
    return get_sequence_sampled_global_score


def score_with_table(
    table: pd.DataFrame, candidates: Mapping[str, Sequence[Mutation]], scorer: Scorer | None = None
) -> list[Row]:
    """(variant_id, score, raw, note) for every candidate, from AntiFold's table for the wild-type structure.

    `raw` is AntiFold's own number for the mutant sequence, a mean negative log-likelihood over all rows (lower is
    a better fit). `score` is the change from the wild type multiplied by the number of rows, so it is the total
    log-likelihood ratio in nats and higher is better. The factor only undoes AntiFold's averaging.

    The table has one distribution per position, so a candidate with several mutations gets no score and a note:
    a lookup would add the single scores, which design_pipeline's scorers never do.
    """
    scorer = scorer or default_scorer()
    missing = [a for a in AMINO if a not in table.columns]
    if missing or not {"pdb_chain", "pdb_res"} <= set(table.columns):
        raise ScoreError("the table needs columns pdb_chain, pdb_res and one per amino acid (" + AMINO + ")")
    wild_type = [str(r) for r in table["pdb_res"]]
    if any(r not in AMINO for r in wild_type):
        odd = sorted({r for r in wild_type if r not in AMINO})
        raise ScoreError(f"the table has residues AntiFold's scoring function cannot take: {odd}")
    sequence = "".join(wild_type)
    rows_of: dict[str, list[int]] = {}
    for index, chain in enumerate(table["pdb_chain"]):
        rows_of.setdefault(str(chain), []).append(index)
    sequences = {chain: "".join(wild_type[i] for i in rows) for chain, rows in rows_of.items()}

    _, wild_mean = scorer(sequence, table)
    scored: list[Row] = []
    for variant_id, mutations in candidates.items():
        if len(mutations) != 1:
            scored.append((variant_id, None, None, f"{len(mutations)} mutations: a lookup would add single scores, so no score"))
            continue
        problem = check_mutations(sequences, mutations)
        if problem:
            scored.append((variant_id, None, None, problem))
            continue
        chain, position, _wt, mut = mutations[0]
        row = rows_of[chain][position - 1]
        _, mean = scorer(sequence[:row] + mut + sequence[row + 1 :], table)
        scored.append((variant_id, (wild_mean - mean) * len(sequence), mean, ""))
    return scored


def write_scores(path: Path, rows: Sequence[Row]) -> Path:
    """scores.csv as design_pipeline reads it (variant_id, score, raw), plus a note column for the empty cells."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["variant_id", "score", "raw", "note"])
        for variant_id, score, raw, note in rows:
            writer.writerow([variant_id, "" if score is None else repr(float(score)), "" if raw is None else repr(float(raw)), note])
    return path


def run_antifold(
    pdb: Path, heavy: str | None, light: str | None, out_dir: Path, src: Path, executable: str | None = None
) -> Path:
    """AntiFold's own command line on one structure; returns the per-residue table it writes (the `table` route).

    Started from its source folder, where it keeps the weights, and with CUDA hidden: AntiFold 0.4.0 moves each
    batch to the GPU inside a forked data-loader worker whenever one is visible, which fails even with
    `--device cpu`. One structure takes seconds on a CPU.

    Without chain ids AntiFold finds the heavy and the light chain itself (on 1N8Z it writes the same table as with
    `--heavy_chain B --light_chain A`), and the table is the one csv it leaves in `out_dir`.
    """
    if (heavy is None) != (light is None):
        raise ScoreError("give both chain ids or neither")
    out_dir = Path(out_dir).resolve()
    chains = ["--heavy_chain", heavy, "--light_chain", light] if heavy else []
    command = [
        executable or sys.executable, "antifold/main.py", "--pdb_file", str(Path(pdb).resolve()), *chains,
        "--out_dir", str(out_dir), "--device", "cpu",
    ]
    done = subprocess.run(command, cwd=src, env={**os.environ, "CUDA_VISIBLE_DEVICES": ""}, capture_output=True, text=True)
    tables = [out_dir / f"{Path(pdb).stem}_{heavy}{light}.csv"] if heavy else sorted(out_dir.glob("*.csv"))
    if done.returncode != 0 or len(tables) != 1 or not tables[0].is_file():
        found = f"; tables found: {[t.name for t in tables]}" if done.returncode == 0 else ""
        raise ScoreError(
            f"AntiFold failed (exit {done.returncode}{found}); last output:\n" + "\n".join((done.stdout + done.stderr).splitlines()[-12:])
        )
    return tables[0]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="antifold_score", description=__doc__.split("\n\n")[0])
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--out_dir", type=Path, required=True)
    parser.add_argument("--mutations_csv", type=Path, required=True, help="variant_id,chain,position,wt,mut, one row per mutation")
    parser.add_argument("--heavy", help="heavy chain id in the structure (required by --route complex; --route table finds it itself)")
    parser.add_argument("--light", help="light chain id in the structure (required by --route complex; --route table finds it itself)")
    parser.add_argument("--route", choices=["complex", "table"], default="complex", help="how to score; see the module docstring")
    parser.add_argument("--device", default="cpu", help="complex route only; with cpu the GPUs are hidden from AntiFold (default cpu)")
    parser.add_argument("--weights", type=Path, help="AntiFold weights; default <antifold_src>/models/model.pt, else the package's own path")
    parser.add_argument("--table", type=Path, help="table route only: an AntiFold table already written for this structure")
    parser.add_argument("--antifold_src", type=Path, default=os.environ.get("ANTIFOLD_SRC"), help="AntiFold source folder (default $ANTIFOLD_SRC)")
    args = parser.parse_args(argv)

    if args.route == "table":
        table_path = args.table
        if table_path is None and not args.antifold_src:
            parser.error("--route table needs --table, or --antifold_src (or ANTIFOLD_SRC)")
        if bool(args.heavy) != bool(args.light):
            parser.error("--heavy and --light go together")
        candidates = read_candidates(args.mutations_csv)
        if table_path is None:
            table_path = run_antifold(args.pdb, args.heavy, args.light, args.out_dir / "antifold", Path(args.antifold_src))
        rows = score_with_table(pd.read_csv(table_path), candidates)
    else:
        if args.table:
            parser.error("--table belongs to --route table")
        if not (args.heavy and args.light):
            parser.error("--route complex needs --heavy and --light")
        candidates = read_candidates(args.mutations_csv)  # before the model loads, so a bad file fails at once
        weights = args.weights
        if weights is None and args.antifold_src and (Path(args.antifold_src) / "models" / "model.pt").is_file():
            weights = Path(args.antifold_src) / "models" / "model.pt"
        if args.device == "cpu":
            os.environ["CUDA_VISIBLE_DEVICES"] = ""  # before AntiFold imports torch: a CPU run must not touch a GPU
        rows = score_with_complex(ComplexScorer(args.pdb, [args.heavy, args.light], device=args.device, weights=weights), candidates)

    path = write_scores(args.out_dir / SCORES, rows)
    empty = [r for r in rows if r[1] is None]
    print(f"{len(rows) - len(empty)} of {len(rows)} candidates scored ({args.route} route) -> {path}")
    for variant_id, _, _, note in empty[:10]:
        print(f"  no score for {variant_id}: {note}")
    return 0
