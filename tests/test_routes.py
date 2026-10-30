"""Install routes: no invented images, and a preflight that says what is here.

The registry used to name ten container images, eight of which did not exist on
Docker Hub - they were GitHub `org/repo` slugs used as registry paths. These
tests are what stops that coming back.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

from design_pipeline import preflight, routes
from design_pipeline.adapters import _gates, _scorers
from design_pipeline.errors import MissingRuntimeError
from design_pipeline.registry import REGISTRY

PACKAGE = Path(routes.__file__).parent

IMAGE_REFERENCE = re.compile(
    r"^[a-z0-9]+(?:[._-][a-z0-9]+)*(?:/[a-z0-9]+(?:[._-][a-z0-9]+)*)+"
    r"(?::[A-Za-z0-9_.-]+)?(?:@sha256:[0-9a-f]{64})?$"
)
"""`namespace/repo`, optionally `:tag` or `@sha256:...`: a container image
reference as docker and apptainer accept one."""

DATA_SUFFIXES = (
    ".py", ".csv", ".json", ".fasta", ".fa", ".pdb", ".cif", ".parquet",
    ".md", ".yaml", ".yml", ".txt", ".pt", ".trb", ".npz", ".a3m",
)

ALLOWED_IMAGES = frozenset({routes.RFDIFFUSION_IMAGE})
"""The only image this package may name: RFdiffusion's, which is published and
pinned by digest. Adding to this set means you have pulled the image and are
prepared to say so in INSTALL.md."""


def image_like_strings() -> list[tuple[str, int, str]]:
    """Every string constant in the package that reads as an image reference.

    A path or a filename reads like one too - `inputs/complex.fasta` matches
    the same shape - so anything with a data suffix, a scheme, a brace or
    whitespace is dropped first.
    """
    found: list[tuple[str, int, str]] = []
    for path in sorted(PACKAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Constant) and isinstance(node.value, str)):
                continue
            value = node.value.strip()
            if not value or len(value) > 200:
                continue
            if any(c in value for c in " \t\n{}\\") or "://" in value:
                continue
            if value.startswith("/") or value.endswith("/"):
                continue
            if value.lower().endswith(DATA_SUFFIXES):
                continue
            if IMAGE_REFERENCE.match(value):
                found.append((path.name, node.lineno, value))
    return found


def test_the_package_names_no_container_image_except_the_pinned_one() -> None:
    """The regression test for the whole episode: an image name that nobody has
    pulled must not reach the registry as a default."""
    unexpected = [
        (name, line, value)
        for name, line, value in image_like_strings()
        if value not in ALLOWED_IMAGES
    ]

    assert not unexpected, (
        "these look like container image references and are not in "
        f"ALLOWED_IMAGES: {unexpected}. Verify the image exists and pin it by "
        "digest, or take it out and give the adapter a route instead."
    )


def test_rfdiffusion_is_pinned_by_digest_not_by_tag() -> None:
    """A tag that moves under a campaign changes the scores without changing
    the config, and the cache key would not notice."""
    assert "@sha256:" in routes.RFDIFFUSION_IMAGE
    assert _scorers.ContainerScorer  # imported for the registry side effect
    from design_pipeline.adapters import rfdiffusion

    assert rfdiffusion.image == routes.RFDIFFUSION_IMAGE


def test_every_registered_adapter_has_a_route() -> None:
    """Preflight reports one row per adapter, so a new adapter without a route
    would show up as a hole in the table."""
    missing = [name for name in REGISTRY if routes.route_for(name) is None]

    assert not missing, f"no route recorded for {missing}; add one to routes.ROUTES"


def test_an_unverified_tool_says_there_is_no_route(tmp_path) -> None:
    """aggrescan3d had an invented image. Nothing replaced it, and the error
    says so rather than naming a route nobody checked."""
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("variant_id,chain,position,wt,mut\n")

    with pytest.raises(MissingRuntimeError, match="no install route has been verified"):
        _gates.aggrescan3d.run([], tmp_path, candidates=candidates)


def test_a_verified_route_is_quoted_in_the_error(tmp_path) -> None:
    """ThermoMPNN has a route, so the error is the route."""
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("variant_id,chain,position,wt,mut\n")

    with pytest.raises(MissingRuntimeError, match="Kuhlman-Lab/ThermoMPNN"):
        _gates.thermompnn.run([], tmp_path, candidates=candidates)


def test_the_one_container_names_no_image_of_its_own(tmp_path) -> None:
    """GearBind is the only genuine environment conflict in the set, so it is
    the only scorer that should be a container - built by you, not named here."""
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("variant_id,chain,position,wt,mut\n")

    assert _scorers.gearbind.tool is None
    with pytest.raises(MissingRuntimeError, match="DeepGraphLearning/GearBind"):
        _scorers.gearbind.run([], tmp_path, candidates=candidates)


def test_a_licence_gated_scorer_says_so(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("variant_id,chain,position,wt,mut\n")

    with pytest.raises(MissingRuntimeError, match="licence-gated"):
        _scorers.foldx.run([], tmp_path, candidates=candidates)


def test_preflight_lists_every_adapter_and_fails_on_a_short_campaign(capsys) -> None:
    """Exit 1 when a default-campaign model does not resolve, so CI can gate on
    it. In this environment antifold and esm_if are absent by design - the test
    suite installs no model."""
    code = preflight.main([])
    printed = capsys.readouterr().out

    for name in REGISTRY:
        if not name.startswith("dummy_"):
            assert name in printed
    assert code == 1
    assert "default campaign NOT ready here" in printed


def test_preflight_probes_nothing_it_cannot_probe_cheaply() -> None:
    """A container route resolves on the pin, not on a pull: checking it
    properly means 8 GB, and preflight has to be runnable on a login node."""
    ok, how = routes.resolve("rfdiffusion")

    assert ok and how == "image pinned, not pulled"
