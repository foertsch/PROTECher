"""Step 1 by download: the experimental structure, when there is one.

Experimental coordinates beat a prediction every time, so this is the default
and the predictors are the fallback. The job is small - fetch by accession,
keep a copy, check the chains are actually there - and the check matters: an
asymmetric unit that is missing the light chain, or that numbers the antibody
differently from the paper, is the kind of thing that otherwise shows up as a
wild-type mismatch 4,000 variants later.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import urlopen

from ..containers import make_id
from ..designs import Design
from ..errors import RunFailedError
from ..structures import AA3TO1, STRUCTURE_SUFFIXES, chain_ids
from ._dummy import _jitter

RCSB = "https://files.rcsb.org/download"
TIMEOUT = 60

RELABELLED_SUFFIXES = (".txt", ".dat", ".ent1", ".pdb1")
"""Suffixes `file_format` is allowed to rename, because they either say nothing
about the format or say it in a spelling the tools downstream do not take: a
biological assembly arrives as .pdb1, and a stripped copy someone mailed you
arrives as whatever their editor saved it as."""

AA1TO3: dict[str, str] = {}
"""One-letter to three-letter, for writing a structure rather than reading one.
Built by inverting AA3TO1 first-wins, so each letter keeps its canonical name
(M is MET, not the selenomethionine that also reads as M)."""
for _three, _one in AA3TO1.items():
    AA1TO3.setdefault(_one, _three)


def _cache_dir(override: Path | str | None) -> Path:
    """Where downloads are kept between campaigns.

    The run directory is content-addressed, so a repeat of the same campaign
    never downloads twice; this is for the next campaign on the same target.
    """
    if override is not None:
        return Path(override)
    configured = os.environ.get("DESIGN_PIPELINE_PDB_CACHE")
    if configured:
        return Path(configured)
    return Path.home() / ".cache" / "design_pipeline" / "pdb"


@dataclass(frozen=True)
class PDBDownload:
    """Fetch by accession from RCSB. In-process: it is one HTTP GET."""

    name: str = "pdb_download"

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        pdb_id: str | None = None,
        chains: Sequence[str] | None = None,
        source: str | Path | None = None,
        file_format: str = "pdb",
        cache_dir: Path | str | None = None,
        extra_args: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        """`source` overrides the download with a local file or an explicit URL.

        That is the escape hatch for a structure that is not in the PDB - a
        collaborator's unreleased coordinates, a stripped and renumbered copy -
        and it is how the test suite and the shipped harness fixtures run with
        no network.

        `file_format` chooses which file to ask RCSB for. It does not rename a
        `source`: that one names its own format and keeps it.
        """
        if pdb_id is None and source is None:
            raise ValueError(
                "pdb_download needs pdb_id='1N8Z', or source=<path or URL>"
            )
        local_name = _source_path(source)
        name = (pdb_id or local_name.stem).lower()
        # One resolved suffix for the copy, the cache entry and the URL: the
        # between-campaign cache is keyed on the file name alone, so a .cif
        # filed under .pdb would be handed to the next campaign that asks for
        # that accession as a PDB.
        suffix = _suffix(local_name, file_format)
        target = output_dir / f"{name}{suffix}"
        target.parent.mkdir(parents=True, exist_ok=True)

        if source is not None and not _is_url(str(source)):
            local = Path(source)
            if not local.is_file():
                raise FileNotFoundError(f"source structure does not exist: {local}")
            shutil.copyfile(local, target)
        else:
            url = str(source) if source is not None else f"{RCSB}/{name}{suffix}"
            cached = _cache_dir(cache_dir) / f"{name}{suffix}"
            if not cached.is_file():
                cached.parent.mkdir(parents=True, exist_ok=True)
                cached.write_bytes(_fetch(url))
            shutil.copyfile(cached, target)

        if chains:
            present = chain_ids(target)
            missing = [c for c in chains if c not in present]
            if missing:
                raise RunFailedError(
                    f"{target.name} does not have chain(s) {', '.join(missing)}; "
                    f"it has {', '.join(present)}. Antibody structures often "
                    f"label the heavy and light chains differently from the "
                    f"paper, and a biological assembly may hold a second copy "
                    f"under other ids."
                )
        return [
            Design(
                id=make_id(target, output_dir),
                path=target,
                parent_id=None,
                producer=self.name,
                meta={"pdb_id": name},
            )
        ]


def _is_url(source: str) -> bool:
    return urlparse(source).scheme in ("http", "https", "ftp")


def _source_path(source: str | Path | None) -> Path:
    """The path part of `source`, so a URL's query string stays out of the
    name we derive from it."""
    if source is None:
        return Path()
    text = str(source)
    return Path(urlparse(text).path if _is_url(text) else text)


def _suffix(source: Path, file_format: str) -> str:
    """The extension to write the copy under.

    `file_format` only says what the bytes are when they come from RCSB, so a
    `source` of its own format wins over it. Writing a .cif as .pdb used to be
    silent and then fatal: `structures.load` dispatches on the suffix, so the
    mmCIF came back out through the PDB reader several steps later.

    A suffix that is neither - .gz above all, since the RCSB download links
    everyone pastes are gzipped - is refused rather than relabelled, because
    renaming it is exactly what hides the problem.
    """
    if source.suffix.lower() in STRUCTURE_SUFFIXES:
        return source.suffix.lower()
    if source.suffix and source.suffix.lower() not in RELABELLED_SUFFIXES:
        raise ValueError(
            f"cannot read {source.name}: this package reads "
            f"{', '.join(STRUCTURE_SUFFIXES)}. Decompress or convert it "
            f"first, or rename it if the suffix is wrong."
        )
    # No suffix, or an uninformative one on a file the caller has told us the
    # format of: `file_format` is then the only thing that knows.
    resolved = f".{file_format.strip().lstrip('.').lower()}"
    if resolved not in STRUCTURE_SUFFIXES:
        raise ValueError(
            f"file_format={file_format!r} is not one this package reads: "
            f"{', '.join(s.lstrip('.') for s in STRUCTURE_SUFFIXES)}"
        )
    return resolved


def _fetch(url: str) -> bytes:
    try:
        with urlopen(url, timeout=TIMEOUT) as response:  # noqa: S310 - fixed scheme
            return response.read()
    except HTTPError as exc:
        if exc.code == 404:
            raise RunFailedError(
                f"{url} returned 404. Check the accession, and note that very "
                f"large structures are mmCIF only - try file_format='cif'."
            ) from exc
        raise RunFailedError(f"{url} returned HTTP {exc.code}") from exc
    except URLError as exc:
        raise RunFailedError(
            f"could not reach {url}: {exc.reason}. Pass source=<local path> to "
            f"work from a file instead."
        ) from exc


def write_ca_pdb(
    path: Path,
    chains: Mapping[str, str],
    *,
    offset: float = 16.0,
    origins: Mapping[str, tuple[float, float, float]] | None = None,
) -> Path:
    """A CA-only PDB from one-letter sequences. For fixtures, not for science.

    Each chain is laid out as an alpha helix - 1.5 A rise, 100 degrees per
    residue, 2.3 A radius - and the chains are offset from one another, so
    CA-CA distances across a chain and between chains are in a plausible range
    and `select_pairs` has something real to filter on.

    `origins` places chains deliberately instead, which is what a fixture
    *complex* needs: at the default spacing no chain is within contact distance
    of another, so there is no interface for `metrics.interface_residues` to
    find and nothing for the confidence metrics to cluster on.
    """
    import math

    lines: list[str] = []
    serial = 1
    for index, (chain, sequence) in enumerate(chains.items()):
        base = (origins or {}).get(chain)
        shift = base[0] if base is not None else index * offset
        lift = base[1] if base is not None else 0.0
        depth = base[2] if base is not None else 0.0
        for position, residue in enumerate(sequence, start=1):
            angle = math.radians(100.0 * position)
            x = 2.3 * math.cos(angle) + shift
            y = 2.3 * math.sin(angle) + lift
            z = 1.5 * position + depth
            name = AA1TO3.get(residue.upper(), "GLY")
            lines.append(
                "ATOM  "
                + f"{serial:5d}"
                + " "
                + " CA "
                + " "
                + f"{name:>3}"
                + " "
                + chain
                + f"{position:4d}"
                + "    "
                + f"{x:8.3f}{y:8.3f}{z:8.3f}"
                + f"{1.0:6.2f}{70.0:6.2f}"
                + "          C"
            )
            serial += 1
        lines.append("TER")
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


# A heavy chain, a light chain and an antigen, with enough sequence variety for
# the liability regexes to have something to find. Short on purpose: the whole
# offline suite runs in under a second.
FIXTURE_CHAINS: Mapping[str, str] = {
    "H": "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG",
    "L": "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRF",
    "A": "TQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSF",
}


@dataclass(frozen=True)
class DummyComplex:
    """Stands in for step 1: a three-chain fixture complex, no network.

    `chains` is overridable so the harness fixtures and a test can pin their
    own sequences, which is what makes the shipped systems reproducible.
    """

    name: str = "dummy_complex"

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        chains: Mapping[str, str] | Sequence[str] | None = None,
        extra_args: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        # A sequence of chain ids means "the fixture chains, these ones only" -
        # which is also what the chain-verification kwarg of pdb_download
        # looks like, so the two adapters stay swappable in a test.
        if chains is None:
            selected = dict(FIXTURE_CHAINS)
        elif isinstance(chains, Mapping):
            selected = dict(chains)
        else:
            selected = {c: FIXTURE_CHAINS[c] for c in chains}
        path = write_ca_pdb(output_dir / "complex.pdb", selected)
        return [
            Design(
                id=make_id(path, output_dir),
                path=path,
                parent_id=None,
                producer=self.name,
                meta={"chains": ",".join(selected)},
            )
        ]


pdb_download = PDBDownload()
dummy_complex = DummyComplex()


# Where a fixture complex puts its chains. The antibody chains sit side by
# side and the antigen against them; `MISPLACED` is the same antigen on the
# opposite face, which is a genuinely wrong pose - different epitope, interface
# RMSD above the clustering threshold - rather than a jittered copy of the
# right one.
DOCKED = (3.0, 7.0, 0.0)
MISPLACED = (3.0, -7.0, 0.0)
ANTIBODY_SPACING = 6.5
AGREEING_FRACTION = 0.7
"""Share of seeds that land on the docked pose. Not 1.0 on purpose: an
ensemble that always agrees would leave `seed_agreement` untested, and
disagreement is the signal the confidence code exists to report."""

MAX_PAE_RESIDUES = 200
"""Above this, the fixture skips the PAE matrix. It is n^2 floats of JSON per
seed and a fixture does not need to scale; the metric that depends on it
reports absent, which is a case worth exercising anyway."""


@dataclass(frozen=True)
class DummyPredictor:
    """Stands in for a co-folding predictor: an ensemble, with confidence.

    `dummy_complex` is enough for a campaign that only needs *a* structure.
    This one exists for the validation notebook, which needs the parts a real
    predictor has and a fixture usually does not: it honours `chain_ids` and
    the input sequences, writes one structure per seed, puts most of them on
    one pose and the rest on another, and reports ipTM, pTM and a PAE matrix.

    So seed agreement, the decoy control and the whole notebook run end to end
    with no GPU, no network and no model weights. The numbers mean nothing; the
    plumbing they exercise is the real thing.
    """

    name: str = "dummy_predictor"

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        chain_ids: Sequence[str] | None = None,
        num_seeds: int | None = 1,
        seed: int = 0,
        extra_args: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        from ._folding import fasta_chains, write_scores_json

        if not inputs:
            raise ValueError(
                f"{self.name} needs the sequences to fold; get_structure writes "
                f"them to a FASTA from antibody_seq and antigen_seq"
            )
        record, sequences = fasta_chains(inputs[0].path)
        ids = list(chain_ids or ("H", "L", "A", "B", "C", "D"))[: len(sequences)]
        if len(ids) < len(sequences):
            raise ValueError(
                f"{len(sequences)} chains to fold and only {len(ids)} chain ids"
            )
        chains = dict(zip(ids, sequences))
        # The FASTA's last record is the antigen, which is the order
        # get_structure writes antibody_seq then antigen_seq in.
        antigen = ids[-1]
        lengths = [len(s) for s in sequences]

        designs: list[Design] = []
        for index in range(max(1, int(num_seeds or 1))):
            this_seed = seed + index
            token = f"{record}:{this_seed}:{':'.join(sequences)}"
            agrees = _jitter(f"{token}:pose") < AGREEING_FRACTION
            origins = {
                chain: (i * ANTIBODY_SPACING, 0.0, 0.0)
                for i, chain in enumerate(ids[:-1])
            }
            origins[antigen] = DOCKED if agrees else MISPLACED

            work = output_dir / f"seed_{this_seed}"
            pdb = write_ca_pdb(work / f"{record}_seed{this_seed}.pdb", chains, origins=origins)
            # Keyed on the antigen sequence, so a decoy antigen gets its own
            # ipTM and the decoy control has something to compare. Nudged up
            # for the dominant pose, as a real predictor's would be.
            base = _jitter(f"{sequences[-1]}:iptm")
            iptm = min(0.95, 0.35 + 0.45 * base + (0.10 if agrees else 0.0))
            ptm = min(0.95, 0.45 + 0.40 * base)
            plddt = 70.0 + 25.0 * _jitter(f"{token}:plddt")
            scores_json = write_scores_json(
                work / f"{pdb.stem}_scores.json",
                plddt=[plddt] * sum(lengths),
                pae=_block_pae(lengths, token),
                ptm=ptm,
                iptm=iptm,
            )
            designs.append(
                Design(
                    id=make_id(pdb, output_dir),
                    path=pdb,
                    parent_id=inputs[0].id,
                    producer=self.name,
                    scores={"plddt": plddt, "ptm": ptm, "iptm": iptm},
                    meta={
                        "scores_json": str(scores_json),
                        "seed": this_seed,
                        "pose": "docked" if agrees else "misplaced",
                    },
                )
            )
        # Ranked by the tool's own confidence, so Complex.path is its top pick
        # and Complex.paths is the ensemble behind it - the same shape a real
        # predictor hands back, and the same reason not to trust the ranking.
        designs.sort(key=lambda d: (-d.scores["iptm"], d.id))
        return designs


def _block_pae(lengths: Sequence[int], token: str) -> list[list[float]] | None:
    """A PAE matrix with worse inter-chain than intra-chain error.

    Rounded and skipped above MAX_PAE_RESIDUES: it is the shape
    `scoring.pae_interaction` reads, not a prediction, and n^2 of it is the
    largest thing a fixture writes.
    """
    total = sum(lengths)
    if total > MAX_PAE_RESIDUES:
        return None
    chain_of: list[int] = []
    for index, length in enumerate(lengths):
        chain_of += [index] * length
    noise = 4.0 * _jitter(f"{token}:pae")
    return [
        [
            round((3.0 if chain_of[i] == chain_of[j] else 12.0) + noise, 1)
            for j in range(total)
        ]
        for i in range(total)
    ]


dummy_predictor = DummyPredictor()
