"""AntiFold's whole-sequence likelihood as a scorer: `antifold_seq`, the counterpart of the matrix scorer `antifold`.

`antifold` reads one distribution per position off AntiFold's table, so it is exact for a single mutation and has no
answer for several. This scorer asks AntiFold's own `score_sequence_in_complex` how likely the mutant chain is, given
the backbone of the antibody, with one forward pass per mutated chain. A candidate with several mutations in a chain
is scored as a whole and never as the sum of its singles. A single mutation scores differently here than in the matrix
lookup (on 20 trastuzumab singles the two correlate at 0.81), so the two are separate columns and neither replaces the
other.

It runs in AntiFold's environment, by path. This file imports nothing from design_pipeline, only the standard library
and AntiFold, so the adapter can start it with AntiFold's own interpreter:

    python _antifold_seq.py --pdb 1N8Z.pdb --heavy B --light A --out_dir out --mutations_csv candidates.csv

`--pdb`, `--out_dir` and `--mutations_csv` are the flags the other scorers take. The result is `antifold_seq.csv` in
`out_dir`, with columns variant_id, score, raw and note. The name differs from `scores.csv` on purpose: the scorer
that wraps this tool writes its own normalised `scores.csv` into the same folder.

What is scored: the heavy and the light chain only, no antigen. A candidate with mutations in both chains gets the sum
of the two chain scores, because each chain is scored from the backbone alone and does not see the other chain's
sequence. A candidate that does not fit the structure gets no score and a note saying why.

Needs AntiFold 0.4.0 or later (Python 3.13 or later), the version this was run against: it provides
`score_sequence_in_complex` and installs next to biotite 2.0. PyPI's 0.3.1 pins biotite==0.38.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Mapping, Protocol, Sequence

AMINO = "ACDEFGHIKLMNPQRSTVWY"
OUTPUT = "antifold_seq.csv"

Mutation = tuple[str, int, str, str]
"""(chain, position, wt, mut); position is 1-based into the chain as given, as in design_pipeline.variants."""
Row = tuple[str, float | None, float | None, str]
"""(variant_id, score, raw, note). An empty score always comes with a note saying why."""


class ScoreError(RuntimeError):
    """The inputs cannot be scored at all, as opposed to one candidate that does not fit."""


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
                "AntiFold 0.4.0 or later is not importable by this interpreter (it provides score_sequence_in_complex). "
                "Run this with the interpreter of an environment that has it: DESIGN_PIPELINE_ANTIFOLD, or executable="
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
    the change in AntiFold's per-residue average times the chain length. Higher is better, the convention of the
    other scorers. `raw` is the log-likelihood of the whole mutant in nats, summed over all chains.

    A candidate that does not fit the structure gets no score and a note (see `check_mutations`).
    """
    sequences = scorer.sequences
    wild = {chain: scorer.log_likelihood(chain, sequence) for chain, sequence in sequences.items()}
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


def write_scores(path: Path, rows: Sequence[Row]) -> Path:
    """variant_id, score, raw, note. An empty score is 'no score for this candidate', and the note says why."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["variant_id", "score", "raw", "note"])
        for variant_id, score, raw, note in rows:
            writer.writerow([variant_id, "" if score is None else repr(float(score)), "" if raw is None else repr(float(raw)), note])
    return path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="_antifold_seq.py", description=__doc__.split("\n\n")[0])
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--out_dir", type=Path, required=True)
    parser.add_argument("--mutations_csv", type=Path, required=True, help="variant_id,chain,position,wt,mut, one row per mutation")
    parser.add_argument("--heavy", required=True, help="heavy chain id in the structure")
    parser.add_argument("--light", required=True, help="light chain id in the structure")
    parser.add_argument("--device", default="cpu", help="with cpu, the GPUs are hidden from AntiFold (default cpu)")
    parser.add_argument("--weights", type=Path, help="AntiFold weights; default is AntiFold's own path, which downloads them if absent")
    args = parser.parse_args(argv)

    candidates = read_candidates(args.mutations_csv)  # before the model loads, so a bad file fails at once
    if args.device == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = ""  # before AntiFold imports torch: a CPU run must not touch a GPU
    rows = score_with_complex(ComplexScorer(args.pdb, [args.heavy, args.light], device=args.device, weights=args.weights), candidates)

    path = write_scores(args.out_dir / OUTPUT, rows)
    empty = [r for r in rows if r[1] is None]
    print(f"{len(rows) - len(empty)} of {len(rows)} candidates scored -> {path}")
    for variant_id, _, _, note in empty[:10]:
        print(f"  no score for {variant_id}: {note}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ScoreError as exc:
        sys.exit(f"antifold_seq: {exc}")
