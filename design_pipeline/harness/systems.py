"""Systems: a named bundle of a structure, a chain, and measured affinities.

The harness exists because on the day of a real campaign there is no feedback.
Nothing tells you which scorer to trust, so the only way to know is to have
measured it beforehand on systems where the answers are already published.

A `System` is deliberately small - an accession, which chain the variants live
on, and a CSV of `variant_id,measured_affinity` - because anything larger would
turn the harness into a second pipeline instead of the ordinary one in a loop.

Ground truth comes from:

* FLAb        https://github.com/Graylab/FLAb   (also s3://graylab-flab)
* AbBiBench   antibody binding benchmark, HuggingFace
* AB-Bind     https://github.com/sarahsirin/AB-Bind-Database

None of it is redistributable here, so the real systems point at files you
download once into the ground-truth directory, and two fixture systems ship
with the package so the harness runs in CI against the dummy scorer.

The one conversion you have to do yourself: variant ids here are
`H:T30Y`, 1-based into the chain as given in the structure file. Published sets
are numbered in IMGT, in Kabat, or in the PDB's own numbering, and a
mis-numbered row is the one error this harness cannot detect for you - except
that `enumerate_variants(only=...)` checks every wild-type residue against the
structure and refuses the ones that disagree.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from ..errors import MissingDependencyError

FIXTURES = Path(__file__).parent / "fixtures"
GROUND_TRUTH_ENV = "DESIGN_PIPELINE_GROUND_TRUTH"
DEFAULT_GROUND_TRUTH = Path("data") / "ground_truth"


@dataclass(frozen=True)
class System:
    """One benchmark: where the structure comes from, and what was measured."""

    name: str
    chain: str
    """The chain the variants are on, as labelled in the structure file."""
    truth: str
    """CSV of variant_id,measured_affinity. Resolved against the fixtures
    directory for a shipped system and the ground-truth directory otherwise."""
    pdb_id: str | None = None
    source: str | None = None
    """A local file instead of a download. Set for the shipped fixtures."""
    model_name: str = "pdb_download"
    chains: tuple[str, ...] = ()
    """Chains that must be present; checked before anything is scored."""
    higher_is_better: bool = False
    """Which way `measured_affinity` points. False for a KD or a ddG, where the
    improvement is the smaller number - the default, because that is what the
    published sets report."""
    citation: str = ""
    fixture: bool = False

    def truth_path(self) -> Path:
        """Where the measurements are, with an error that names the download."""
        if self.fixture:
            return FIXTURES / self.truth
        root = Path(os.environ.get(GROUND_TRUTH_ENV) or DEFAULT_GROUND_TRUTH)
        path = root / self.truth
        if not path.is_file():
            raise MissingDependencyError(
                f"system {self.name!r} needs measured affinities at {path}, "
                f"which does not exist. Download them from {self.citation or 'FLAb'} "
                f"and save a CSV of variant_id,measured_affinity there, or set "
                f"{GROUND_TRUTH_ENV} to point somewhere else. Variant ids are "
                f"1-based into chain {self.chain} as given in the structure. "
                f"The shipped fixture systems "
                f"({', '.join(n for n, s in SYSTEMS.items() if s.fixture)}) "
                f"run with no download."
            )
        return path

    def source_path(self) -> str | None:
        if self.source is None:
            return None
        return str(FIXTURES / self.source) if self.fixture else self.source


SYSTEMS: dict[str, System] = {
    # Shipped fixtures. Not science: a CA-only helix bundle and invented
    # affinities, so that `run_systems` and `report` are exercised end to end
    # offline, in CI, against the dummy scorer. See fixtures/README.md.
    "fixture_a": System(
        name="fixture_a",
        chain="H",
        truth="fixture_a.csv",
        source="fixture_a.pdb",
        model_name="pdb_download",
        chains=("H", "L", "A"),
        higher_is_better=True,
        fixture=True,
    ),
    "fixture_b": System(
        name="fixture_b",
        chain="H",
        truth="fixture_b.csv",
        source="fixture_b.pdb",
        model_name="pdb_download",
        chains=("H", "A"),
        higher_is_better=True,
        fixture=True,
    ),
    # Real systems. The structure downloads; the measurements you supply once.
    "1n8z_her2": System(
        name="1n8z_her2",
        chain="B",
        truth="1n8z_her2.csv",
        pdb_id="1N8Z",
        # 1N8Z is trastuzumab Fab (chains A light, B heavy) with HER2 domain IV
        # as chain C. The heavy chain is 'B', not 'H': this is exactly the
        # labelling trap `chains` is checked for.
        chains=("A", "B", "C"),
        citation="FLAb (github.com/Graylab/FLAb)",
    ),
    "4fqi_flu": System(
        name="4fqi_flu",
        chain="H",
        truth="4fqi_flu.csv",
        pdb_id="4FQI",
        chains=("H", "L", "A"),
        citation="AB-Bind (github.com/sarahsirin/AB-Bind-Database)",
    ),
}

FIXTURE_SYSTEMS: tuple[str, ...] = tuple(
    name for name, system in SYSTEMS.items() if system.fixture
)


def get_system(name: str) -> System:
    try:
        return SYSTEMS[name]
    except KeyError:
        known = ", ".join(sorted(SYSTEMS))
        raise KeyError(
            f"no system {name!r}. Known: {known}. Add your own with "
            f"design_pipeline.harness.SYSTEMS['mine'] = System(...)."
        ) from None
