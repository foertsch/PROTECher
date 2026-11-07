"""Every route's declared environment must be able to answer that route's probe.

The failure this exists to catch: a tool that lives in another environment,
probed with `module` or `executable`. `importlib.util.find_spec` and
`shutil.which` both ask *this* interpreter and *this* PATH, so for a tool that
is deliberately not here they do not say "not installed" - they say "not
installed" when it is installed, every time, and there is no way to tell the two
apart from the outside.

It is not hypothetical. `boltz` ran the whole final screen out of `.venv-boltz`
and `python -m design_pipeline.preflight` in Env A reported `import boltz
failed` throughout. The reasoning was written out for `antifold` when that route
was added and simply not applied to the nine routes that came from pyproject.

So: a route in Env B must declare an `env_var`, setting it must resolve, and the
adapter must read the same variable the table names. `test_the_check_catches_the
_nine_routes_it_was_written_for` feeds the pre-change table back in, so the test
is demonstrably the one that would have caught this rather than a claim that it
is.
"""

from __future__ import annotations

import dataclasses

import pytest

from design_pipeline import routes
from design_pipeline.registry import REGISTRY

# The nine that declared `module` or `executable` only. Kept as data so that the
# reconstruction below is exact and so the list cannot quietly shrink.
PRE_CHANGE_WITHOUT_A_PATH_PROBE = frozenset(
    {
        "chai",
        "chai1",
        "boltz",
        "boltz2",
        "protenix",
        "colabfold",
        "protein_mpnn",
        "ligand_mpnn",
        "soluble_mpnn",
    }
)

BY_PATH = (routes.ENV_B,)
"""Environments whose tools cannot be imported from here, so the only probe that
can succeed is a configured path. MANUAL is not in the list because a manual
route may legitimately have nothing automated at all - `foldx` is checked
separately, by name, because it does."""


def missing_a_path_probe(table: dict[str, routes.Route]) -> list[str]:
    """Routes that live in another environment and have no variable to probe.

    The whole check, as a function over a table, so the same logic can be run
    against the current `ROUTES` and against a reconstruction of the one before
    the fix.
    """
    return sorted(
        name
        for name, route in table.items()
        if route.env in BY_PATH and route.env_var is None
    )


def pre_change_table() -> dict[str, routes.Route]:
    """Today's table with the nine `env_var`s taken back off."""
    return {
        name: (
            dataclasses.replace(route, env_var=None)
            if name in PRE_CHANGE_WITHOUT_A_PATH_PROBE
            else route
        )
        for name, route in routes.ROUTES.items()
    }


def test_every_env_b_route_can_be_probed_from_outside_its_environment() -> None:
    """The drift test. A new Env B route with no variable fails here."""
    missing = missing_a_path_probe(dict(routes.ROUTES))

    assert not missing, (
        f"these routes live in Env B and declare no env_var: {missing}. From "
        f"Env A their probe asks this interpreter about a package that is "
        f"deliberately not in it, so it always answers no - which is exactly "
        f"what happened to `boltz` through the whole final screen. Give each "
        f"one the variable its adapter reads, as `antifold` and `thermompnn` "
        f"have; routes.py:118 has the reasoning."
    )


def test_the_check_catches_the_nine_routes_it_was_written_for() -> None:
    """Run the check against the pre-change table: it has to find all nine.

    Without this, the test above passes on an empty rule just as happily as on
    a correct one.
    """
    assert missing_a_path_probe(pre_change_table()) == sorted(
        PRE_CHANGE_WITHOUT_A_PATH_PROBE
    )


@pytest.mark.parametrize(
    "name", sorted(n for n, r in routes.ROUTES.items() if r.env_var is not None)
)
def test_setting_the_variable_resolves_the_route(name: str, monkeypatch) -> None:
    """A declared variable has to be the one `resolve` actually consults.

    Checked per route rather than once, because the variable is read by name: a
    typo in one entry is invisible to a test that only exercises another.
    """
    route = routes.ROUTES[name]
    assert route.env_var is not None  # the parametrization guarantees it
    monkeypatch.setenv(route.env_var, "/opt/some-environment/bin/tool")

    ok, how = routes.resolve(name)

    assert ok, f"{name} declares ${route.env_var} and setting it did not resolve"
    assert route.env_var in how, (
        f"{name} resolved, but `how` is {how!r} and does not name "
        f"${route.env_var}. preflight prints `how`, and a row that does not say "
        f"which variable it believed is not actionable."
    )


@pytest.mark.parametrize(
    "name", sorted(n for n, r in routes.ROUTES.items() if r.env_var is not None)
)
def test_an_unset_variable_does_not_resolve_by_itself(name: str, monkeypatch) -> None:
    """The fallback has to stay a fallback.

    With the variable unset, an Env B route may only resolve on its own
    `module`/`executable` - which is the case that lets preflight run *inside*
    Env B and still find the tool. In this suite no model is installed, so
    nothing should resolve. A route that resolves here anyway is claiming a tool
    is present on a machine that has none.
    """
    route = routes.ROUTES[name]
    assert route.env_var is not None
    monkeypatch.delenv(route.env_var, raising=False)
    # An empty value is "not configured", not a path. The node's env.sh exports
    # these unconditionally, so an empty one is the realistic failure.
    monkeypatch.setenv(route.env_var, "   ")

    ok, _how = routes.resolve(name)

    assert not ok, (
        f"{name} resolved with ${route.env_var} set to whitespace. An empty "
        f"export is not a configured path."
    )


def test_env_a_and_bundled_routes_do_not_claim_a_path() -> None:
    """The converse of the drift test, and the mistake it invites.

    An Env A route is reachable by import and a bundled one needs nothing, so a
    variable on either would be answering a question nobody asked - and because
    `env_var` is consulted first, it would also override the import probe that
    is actually correct there.
    """
    wrong = sorted(
        name
        for name, route in routes.ROUTES.items()
        if route.env in (routes.ENV_A, routes.BUILTIN) and route.env_var is not None
    )

    assert not wrong, (
        f"{wrong} are in Env A or bundled and declare an env_var. `resolve` asks "
        f"the variable first, so this hides the import probe that is the right "
        f"one for them."
    )


def test_foldx_is_manual_and_still_reached_by_path() -> None:
    """The one MANUAL route with a variable, which is not a contradiction.

    `env` says where a tool comes from - a licence, in this case - and `env_var`
    says how it is reached. FoldX is obtained by hand and then run by path.
    """
    route = routes.ROUTES["foldx"]

    assert route.env == routes.MANUAL
    assert route.env_var == "DESIGN_PIPELINE_FOLDX"


def test_foldx_needs_env_lists_every_variable_the_adapter_reads() -> None:
    """All three, because `_foldx_score.py` reads all three.

    `needs_env` listed only FOLDX_BINARY, so preflight never reported
    DESIGN_PIPELINE_FOLDX_PDB as unset - the one with no sensible default, since
    it names the structure every ddG is computed against.
    """
    needs = set(routes.ROUTES["foldx"].needs_env)

    assert needs == {
        "FOLDX_BINARY",
        "DESIGN_PIPELINE_FOLDX_PDB",
        "DESIGN_PIPELINE_FOLDX",
    }


def adapter_env_vars(adapter: object) -> set[str]:
    """The variables an adapter will actually read.

    Both the adapter itself and, for the scorers and gates that wrap one, its
    `tool` - `ContainerScorer` holds a `ProcessAdapter` and the variable is on
    the inner one.
    """
    found: set[str] = set()
    for candidate in (adapter, getattr(adapter, "tool", None)):
        value = getattr(candidate, "env_var", None)
        if isinstance(value, str) and value:
            found.add(value)
    return found


def test_the_adapter_reads_the_variable_the_table_names() -> None:
    """The invariant in `Route.env_var`, enforced.

    Preflight reporting one variable while the run reads another is worse than
    either being wrong on its own: the table says the tool is configured and the
    run cannot find it, and nothing points at the discrepancy.
    """
    mismatched: list[tuple[str, str | None, set[str]]] = []
    for name, adapter in REGISTRY.items():
        if name.startswith("dummy_"):
            continue
        route = routes.route_for(name)
        if route is None:
            continue
        declared = adapter_env_vars(adapter)
        if not declared:
            continue
        if {route.env_var} != declared:
            mismatched.append((name, route.env_var, declared))

    assert not mismatched, (
        "the table and the adapter disagree about which variable points at the "
        f"tool: {mismatched} (route env_var, then what the adapter reads). "
        "Preflight prints the first and the run uses the second."
    )


def test_an_env_b_adapter_that_shells_out_declares_the_variable() -> None:
    """The other direction: an adapter whose route has a variable should read it.

    Chai is the documented exception and is asserted as one rather than skipped,
    so that wiring it up later makes this test fail and ask to be updated.
    """
    silent: list[str] = []
    for name, adapter in REGISTRY.items():
        route = routes.route_for(name)
        if route is None or route.env_var is None or name.startswith("dummy_"):
            continue
        if not adapter_env_vars(adapter):
            silent.append(name)

    assert sorted(silent) == ["chai", "chai1"], (
        f"these routes name a variable that their adapter never reads: {silent}. "
        f"Only chai and chai1 may: that adapter calls run_inference in process, "
        f"so there is no path to point at and the variable is a preflight probe "
        f"only. Anything else here is a route and an adapter that disagree."
    )
