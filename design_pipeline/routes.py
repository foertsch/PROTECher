"""How each model's tool is obtained, and whether we have checked that route.

This table exists because the previous arrangement did not have one. Ten
adapters named ten container images; eight of those images do not exist, and
the names turned out to be GitHub `org/repo` slugs used as if they were
registry paths. A default that 404s is worse than no default: it fails inside
`docker pull` with a registry error instead of telling you what to install.

So there are no image defaults anywhere in this package except RFdiffusion,
which really is published and is pinned here by digest. Everything else
resolves through this table, and an adapter whose route is not satisfied raises
`MissingRuntimeError` quoting the route verbatim.

`verified` means someone ran the check and recorded the result on VERIFIED_ON -
not that it is likely, nor that the tool's README implies it. An unchecked
route says so, and an adapter with no route at all says that too rather than
guessing one. The target is two environments and two containers, not ten
images:

    Env A   the default campaign: this package, the sequence-only scorers
            (ablang2, antiberty, esm1v, esm2_650m), fair-esm and anarci.
            Build and test this first.
    Env B   whatever will not resolve into A - AntiFold, ThermoMPNN, and the
            co-folding models, which pin torch against each other. Reached by
            path (executable= / DESIGN_PIPELINE_*), not by import.
    GearBind    a container because its TorchDrug stack pins Python <=3.10 and
            an old torch: the one conflict that a second environment cannot
            fix, since nothing in it can share a process with us.
    RFdiffusion a container because it already had a published one, which is
            pinned by digest below.

AntiFold is in Env B rather than A despite being a default scorer: it pins
biotite==0.38 against this package's biotite>=2.0, and pip calls that
combination ResolutionImpossible. Whether AntiFold and ThermoMPNN can share
one Env B is NOT established - their pins may force two - but that is a
configuration choice, since each is reached by its own path.

Nothing here invokes a tool. This module is data plus two probes, so importing
it stays free and `preflight` can read it with no GPU and no network.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from .errors import MissingRuntimeError

VERIFIED_ON = "2026-10-30"
"""The day the verified routes below were last checked."""

ENV_A = "A"
ENV_B = "B"
CONTAINER = "container"
MANUAL = "manual"
BUILTIN = "-"

STRUCTURE = "structure"
SCORE = "score"
GATE = "gate"
DENOVO = "denovo"


@dataclass(frozen=True)
class Route:
    """One tool: where it comes from, and how to tell whether it is here.

    `env_var` is checked first and on its own: a tool in another environment is
    not importable from this one, so asking whether the path has been
    configured is the only probe that can succeed from here. After that, at
    most one of the rest is set. `module` is an import name, `executable` a
    name on PATH, `image` a container reference; none of the four means there
    is nothing to probe, because the route is manual or not established, and
    the adapter will raise rather than assemble a command.
    """

    step: str
    env: str
    how: str
    """The install route, verbatim, as it goes into the error message."""
    verified: bool = False
    module: str | None = None
    executable: str | None = None
    image: str | None = None
    env_var: str | None = None
    """The variable that points at this tool's interpreter or binary, for an
    Env B route. The same variable the adapter reads, so preflight and the run
    agree about where the tool is."""
    needs_env: tuple[str, ...] = ()
    """Further variables the adapter's command reads. They do not decide
    `resolve`, because a keyword argument can stand in for each, so preflight
    only reports the ones that are unset."""
    note: str = ""


RFDIFFUSION_IMAGE = (
    "rosettacommons/rfdiffusion"
    "@sha256:b88dc5d7ccd166557fb9b922bbaeb15bdf8d0568907c50405a8bf282e67254fe"
)
"""The one published image in the set, pinned to the digest that `latest`
pointed at on 2026-10-30 (8.1 GB, pushed 2025-06-25). Pinned rather than
tracked: a tag that moves under a campaign changes the scores without changing
the config, and the cache key would not notice."""

PLM_EXTRA = "pip install 'design-pipeline[plm]'"


ROUTES: Mapping[str, Route] = {
    # ------------------------------------------------------------- Env A, step 1
    "pdb_download": Route(
        step=STRUCTURE, env=BUILTIN, how="bundled: one HTTP GET", verified=True
    ),
    # ------------------------------------------------------------- Env B, step 3
    # Env B, not A, despite being a default scorer: `pip install antifold` pins
    # biotite==0.38 against this package's biotite>=2.0 and pip calls the
    # combination ResolutionImpossible. It is reached by path, so the probe is
    # the variable the adapter reads - `module` would ask this interpreter
    # about a package that is deliberately not in it, and always answer no.
    "antifold": Route(
        step=SCORE,
        env=ENV_B,
        how=(
            "pip install antifold into its own environment, then point the "
            "adapter at it with executable= or DESIGN_PIPELINE_ANTIFOLD"
        ),
        verified=True,
        module="antifold",
        env_var="DESIGN_PIPELINE_ANTIFOLD",
        needs_env=("ANTIFOLD_SRC",),
        note=(
            "PyPI 0.3.1, BSD-3. Scored through the repo's antifold_score module, "
            "which also needs ANTIFOLD_SRC: the AntiFold source folder that "
            "holds models/model.pt."
        ),
    ),
    # The whole-sequence counterpart of `antifold`, run from the same
    # environment, so it reads the same variable. It needs 0.4.0 or later, which
    # is not on PyPI: that release is the GitHub source, and unlike the PyPI
    # ones it does not pin biotite, so it installs next to biotite>=2.0.
    "antifold_seq": Route(
        step=SCORE,
        env=ENV_B,
        how=(
            "install AntiFold 0.4.0 or later from github.com/oxpig/AntiFold "
            "(Python 3.13 or later) into the environment DESIGN_PIPELINE_ANTIFOLD "
            "points at, with executable= as the alternative"
        ),
        verified=True,
        module="antifold",
        env_var="DESIGN_PIPELINE_ANTIFOLD",
        note=(
            "Run on 1N8Z on 2026-10-30 in a Python 3.13 environment holding "
            "biotite 2.0.0 and numpy 2.5.4. Weights, 567 MB, download on first "
            "run into the package folder; pass weights= to reuse a download."
        ),
    ),
    # ------------------------------------------------------------- Env A, step 3
    "esm_if": Route(
        step=SCORE,
        env=ENV_A,
        how="pip install fair-esm, plus torch-geometric, torch-scatter, biotite",
        verified=True,
        module="esm",
        note=(
            "torch-scatter ships no wheel for current torch and builds from "
            "source; INSTALL.md records what the resolver actually said."
        ),
    ),
    "consensus": Route(
        step=SCORE,
        env=BUILTIN,
        how="bundled, but needs an OAS frequency table (DESIGN_PIPELINE_OAS_TABLE)",
        verified=True,
    ),
    # Sequence-only scorers. They never open the structure, so they are the
    # columns that survive a wrong predicted complex - which is why they belong
    # in Env A rather than behind a container.
    "ablang2": Route(
        step=SCORE,
        env=ENV_A,
        how="pip install ablang2",
        verified=True,
        module="ablang2",
        note="Antibody-specific. Also reachable through the plm extra.",
    ),
    "esm1v": Route(
        step=SCORE, env=ENV_A, how=PLM_EXTRA, verified=True, module="transformers"
    ),
    "esm2_650m": Route(
        step=SCORE, env=ENV_A, how=PLM_EXTRA, verified=True, module="transformers"
    ),
    "antiberty": Route(
        step=SCORE,
        env=ENV_A,
        how="pip install antiberty",
        module="antiberty",
        note="Not re-checked on VERIFIED_ON. Skipping it costs one column.",
    ),
    # ------------------------------------------------------------- Env B, step 4
    "thermompnn": Route(
        step=GATE,
        env=ENV_B,
        how=(
            "git clone https://github.com/Kuhlman-Lab/ThermoMPNN and build its "
            "environment.yaml with mamba, then point the adapter at that "
            "environment with executable= or DESIGN_PIPELINE_THERMOMPNN"
        ),
        verified=True,
        executable="thermompnn",
        env_var="DESIGN_PIPELINE_THERMOMPNN",
        note=(
            "Its README warns the environment.yaml can resolve a CPU-only "
            "torch build. Check torch.cuda.is_available() inside Env B."
        ),
    ),
    "liabilities": Route(
        step=GATE, env=BUILTIN, how="bundled: sequence regex", verified=True
    ),
    # ---------------------------------------------------------- The one container
    "gearbind": Route(
        step=SCORE,
        env=CONTAINER,
        how=(
            "build an image from https://github.com/DeepGraphLearning/GearBind "
            "and set it with REGISTRY['gearbind'] - no published image exists"
        ),
        verified=True,
        note=(
            "TorchDrug pins Python 3.7-3.10 and an old torch: the one genuine "
            "conflict in this set, hence a container rather than Env B."
        ),
    ),
    "rfdiffusion": Route(
        step=DENOVO,
        env=CONTAINER,
        how=f"docker pull {RFDIFFUSION_IMAGE}",
        verified=True,
        image=RFDIFFUSION_IMAGE,
        note="The only published image of the ten. Weights are inside it.",
    ),
    # ------------------------------------- Licence-gated, manual, not automated
    "rosetta_ddg": Route(
        step=SCORE,
        env=MANUAL,
        how="licence-gated: obtain Rosetta yourself, then set REGISTRY['rosetta_ddg']",
        verified=True,
        note="An official image does exist; see INSTALL.md. Not configured here.",
    ),
    "foldx": Route(
        step=SCORE,
        env=MANUAL,
        how="licence-gated: your own FoldX binary. Set FOLDX_BINARY, DESIGN_PIPELINE_FOLDX_PDB and DESIGN_PIPELINE_FOLDX (the Python that runs adapters/_foldx_score.py)",
        verified=True,
        env_var="DESIGN_PIPELINE_FOLDX",
        needs_env=("FOLDX_BINARY",),
    ),
    "alphafold3": Route(
        step=STRUCTURE,
        env=MANUAL,
        how=(
            "licence-gated: weights come from DeepMind on request and the image "
            "is your own build. Set REGISTRY['alphafold3']"
        ),
        verified=True,
    ),
    # -------------------------------------------------- No route established yet
    # These four carried invented image names. Nothing replaces them, because
    # we have not checked how they install. Saying so is the point.
    "aggrescan3d": Route(
        step=GATE,
        env=MANUAL,
        how="",
        note="Carried an invented image name. Install route not checked.",
    ),
    "netsolp": Route(
        step=GATE,
        env=MANUAL,
        how="",
        note="Carried an invented image name. Install route not checked.",
    ),
    "abmpnn": Route(
        step=SCORE,
        env=MANUAL,
        how="",
        note=(
            "Carried an invented image name, and AbMPNN has no public "
            "repository: the weights are request-gated. Route not checked."
        ),
    ),
    "proteinmpnn_score": Route(
        step=SCORE,
        env=MANUAL,
        how="",
        note=(
            "Carried an invented image name. The design-side MPNN adapters "
            "install from the mpnn extra; scoring was not checked."
        ),
    ),
    # ------------------------------- Pip extras that predate this table. Their
    # routes come from pyproject and were not re-checked on VERIFIED_ON.
    "chai": Route(
        step=STRUCTURE, env=ENV_B, how="pip install 'design-pipeline[chai]'",
        module="chai_lab",
    ),
    "chai1": Route(
        step=STRUCTURE, env=ENV_B, how="pip install 'design-pipeline[chai]'",
        module="chai_lab",
    ),
    "boltz": Route(
        step=STRUCTURE, env=ENV_B, how="pip install 'design-pipeline[boltz]'",
        module="boltz",
    ),
    "boltz2": Route(
        step=STRUCTURE, env=ENV_B, how="pip install 'design-pipeline[boltz]'",
        module="boltz",
    ),
    "protenix": Route(
        step=STRUCTURE, env=ENV_B, how="pip install 'design-pipeline[protenix]'",
        module="protenix",
    ),
    "colabfold": Route(
        step=STRUCTURE, env=ENV_B, how="pip install 'design-pipeline[colabfold]'",
        executable="colabfold_batch",
    ),
    "protein_mpnn": Route(
        step=DENOVO, env=ENV_B, how="pip install 'design-pipeline[mpnn]'",
        executable="ligandmpnn",
    ),
    "ligand_mpnn": Route(
        step=DENOVO, env=ENV_B, how="pip install 'design-pipeline[mpnn]'",
        executable="ligandmpnn",
    ),
    "soluble_mpnn": Route(
        step=DENOVO, env=ENV_B, how="pip install 'design-pipeline[mpnn]'",
        executable="ligandmpnn",
    ),
}

EXTRAS: Mapping[str, Route] = {
    # In Env A but not an adapter: `numbering` imports ANARCI directly. Listed
    # so that building Env A from this file does not quietly drop it.
    "anarci": Route(
        step="numbering",
        env=ENV_A,
        how="conda install -c bioconda anarci",
        verified=True,
        module="anarci",
        note="Fills the imgt column. Without it that column is empty.",
    ),
}

DUMMY = Route(
    step="fixture",
    env=BUILTIN,
    how="bundled: offline fixture, no GPU and no container",
    verified=True,
)
"""Every `dummy_*` adapter. They are what makes the test suite and the shipped
harness systems run on a laptop, so they always resolve."""


def route_for(name: str) -> Route | None:
    if name.startswith("dummy_"):
        return DUMMY
    return ROUTES.get(name) or EXTRAS.get(name)


def resolve(name: str) -> tuple[bool, str]:
    """(is it here, how we can tell). The probe only; nothing is run.

    An Env B route resolves when its variable is set, which is a claim about
    configuration rather than about the tool: this interpreter cannot see into
    that environment, so nothing here checks that the path exists or that the
    package is in it. The run is where that is found out.

    A container route counts as resolved once the image is pinned in this file,
    because checking it properly means pulling 8 GB. `preflight` prints that
    distinction rather than implying the pull has happened.
    """
    route = route_for(name)
    if route is None:
        return False, "no route recorded"
    if route.env_var is not None:
        # Asked first, and answered from here: a tool in another environment
        # cannot be imported into this one, so a configured path is the only
        # evidence this interpreter can have that the tool exists. Falls
        # through when the variable is unset, which is what lets preflight run
        # inside Env B itself and still find the import.
        configured = configured_path(route.env_var)
        if configured is not None:
            return True, f"${route.env_var}={configured}"
    if route.module is not None:
        found = importlib.util.find_spec(route.module) is not None
        return found, f"import {route.module}" + ("" if found else " failed")
    if route.executable is not None:
        found = shutil.which(route.executable)
        return bool(found), found or f"{route.executable} not on PATH"
    if route.image is not None:
        return True, "image pinned, not pulled"
    if route.env == BUILTIN:
        return True, "bundled"
    return False, "not configured"


def require(name: str) -> None:
    """Raise unless `name`'s tool is here. The error carries the route.

    Called by the adapters that used to carry an invented image, so the failure
    arrives before a runtime is even detected and names something actionable.
    """
    ok, _how = resolve(name)
    if ok:
        return
    route = route_for(name)
    if route is None or not route.how:
        detail = (
            f"no install route has been verified for {name}, so this package "
            f"does not configure it. Find the tool's own instructions, then "
            f"set REGISTRY[{name!r}] to an adapter that runs it"
        )
    else:
        detail = (
            f"{name} is not installed. Verified route (checked "
            f"{VERIFIED_ON}), environment {route.env}: {route.how}"
        )
    note = f". Note: {route.note}" if route is not None and route.note else ""
    raise MissingRuntimeError(
        f"{detail}{note}. `python -m design_pipeline.preflight` prints the "
        f"whole table."
    )


def configured_path(variable: str) -> Path | None:
    """An executable an environment variable points at, or None.

    Env B lives in its own interpreter, so the adapters that run there take
    their path from the environment rather than from PATH.
    """
    value = os.environ.get(variable, "").strip()
    return Path(value) if value else None
