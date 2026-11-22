"""Environment names in prose must resolve to something in the manifest.

`test_route_drift.py` guards the other direction: that a route's declared
environment can answer that route's probe. It says nothing about the words. So
prose has drifted twice now, in the same way both times - a name outlived the
thing it named, and nothing failed:

1. INSTALL.md told the reader to install "into Env A, B or C". There has never
   been an Env C anywhere in this repository, in the manifest, in `routes.py`
   or in `provision/`. It was a third place that did not exist.
2. `routes.py` opened by stating the target as "two environments and two
   containers" long after Env B had fragmented into five, which
   `provision/manifest.toml` had recorded all along. The docstring a reader
   meets first was describing an arrangement the inventory had replaced.

Both were readable by anyone who looked at the manifest next to the prose, and
nobody did, because looking is a thing a person has to remember and a test does
not. So: `provision/manifest.toml` is the inventory, and these checks ask the
prose to agree with it.

The checks are pure functions over text rather than assertions over files, so
`test_the_checks_catch_the_two_drifts_they_were_written_for` can feed the
pre-change sentences back in. The suite is therefore demonstrably the one that
would have caught both, rather than a claim that it is.

`routes.py`'s `ENV_*` constants are the coarser question of whether a route is
reached by import or by path, not a pointer into the manifest. Their CONSTANT
names keep the letters - `ENV_A`, `ENV_B` - and the letter check reads those,
so `Env A` stays legal prose while the constant exists and `Env C` does not.
Their VALUES print, in preflight's `env` column and in `require`'s error, so
`ENV_A` is checked against the manifest: it names exactly one environment and
must name it correctly. `ENV_B` is not checked that way, because it covers five
environments and no single key would be honest.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import pytest

from design_pipeline import routes

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "provision" / "manifest.toml"

SCANNED = (
    "INSTALL.md",
    "docs/provisioning.md",
    "docs/cluster_gpu_node_rebuild.md",
    "design_pipeline/routes.py",
)
"""The prose that names environments. `routes.py` is in the list because its
module docstring is where drift 2 lived: it is the first thing a reader of the
route table meets, and it is prose whatever file it sits in.

This file is NOT scanned, and must not be added: the reconstructions below
quote the drifted sentences verbatim as test data, so scanning it would trip
every check by design."""

NUMBER_WORDS = {
    "one": 1,
    "two": 2,
    "three": 3,
    "four": 4,
    "five": 5,
    "six": 6,
    "seven": 7,
    "eight": 8,
    "nine": 9,
    "ten": 10,
    "eleven": 11,
    "twelve": 12,
}

TOTAL_CLAIM_MARKERS = ("the target is", "there are", "in total", "all told")
"""A count is only checked when the sentence announces a total.

Without this the check fires on true statements about a subset - "Two
environments install from a git tree and neither commit is known" is correct
and will stay correct while the inventory grows. A count is drift only when it
claims to count everything, and in prose that claim is made in words. The cost
is real: a total stated without one of these phrases is not guarded. That is
the trade this check accepts, because a drift test with false positives gets
deleted and then guards nothing at all."""


# ------------------------------------------------------------- the manifest


def manifest() -> dict:
    with MANIFEST.open("rb") as handle:
        return tomllib.load(handle)


def env_keys() -> frozenset[str]:
    """Every `[env.*]` key. The source of truth every check below resolves
    against."""
    return frozenset(manifest()["env"])


def env_scripts() -> frozenset[str]:
    """The `script` of every environment, as the manifest spells it."""
    return frozenset(
        env["script"] for env in manifest()["env"].values() if "script" in env
    )


def utility_scripts() -> frozenset[str]:
    """Scripts in `provision/` that no environment claims - `_common.sh` and
    the scratch cleaner.

    Derived rather than listed, so adding a utility does not mean remembering
    to edit this file. A script is a utility exactly when it is on disk and no
    environment names it.
    """
    on_disk = {
        f"provision/{path.name}" for path in (REPO / "provision").glob("*.sh")
    }
    return frozenset(on_disk - env_scripts())


def defined_letters() -> frozenset[str]:
    """The environment letters `routes.py` defines, read off the constant
    NAMES - `ENV_A` defines the letter A - rather than restated here, so that
    adding an `ENV_C` makes the prose legal in the same commit that makes it
    real.

    The first version of this read the constant's VALUE and kept the ones a
    single character long, which worked only while `ENV_A` happened to be "A".
    When the values became "core" and "by path" - because they print, in
    preflight's `env` column and in `require`'s error - it silently returned
    nothing, and the check flipped from "a letter must be defined" to "no
    letter may be named", failing on 58 lines of correct prose. The letter
    lives in the name; that is what this reads.
    """
    return frozenset(
        name.removeprefix("ENV_")
        for name, value in vars(routes).items()
        if re.fullmatch(r"ENV_[A-Z]", name) and isinstance(value, str)
    )


# ---------------------------------------------------------------- the scans


def scripts_named_in(text: str) -> frozenset[str]:
    return frozenset(re.findall(r"provision/[A-Za-z0-9_]+\.sh", text))


LETTERS = re.compile(r"\bEnv\s+([A-Z](?:\s*(?:,|\bor\b|\band\b)\s*[A-Z])*)\b")
"""`Env A`, and also the elided list form.

The plain `\\bEnv ([A-Z])\\b` this started as read "Env A, B or C" as a single
mention of Env A, which is exactly the sentence that prompted the check: the
undefined letter is the one the word "Env" is furthest from. The trailing `\\b`
does the work - on "Env A and ThermoMPNN" the repetition tries to swallow the
`T`, fails the boundary, and backtracks to just `A`."""


def letters_named_in(text: str) -> frozenset[str]:
    return frozenset(
        letter
        for run in LETTERS.findall(text)
        for letter in re.findall(r"[A-Z]", run)
    )


def total_claims_in(text: str) -> list[tuple[str, int]]:
    """Sentences that announce how many environments there are, as
    (sentence, claimed count).

    Sentence-splitting on `.` is crude and would mis-split "3.12". It does not
    matter here: the markers and the number only have to land in the same
    window, and a window that is too wide errs toward checking more prose, not
    less.
    """
    found: list[tuple[str, int]] = []
    flat = " ".join(text.split())
    for sentence in re.split(r"(?<=[.:;])\s+", flat):
        low = sentence.lower()
        if not any(marker in low for marker in TOTAL_CLAIM_MARKERS):
            continue
        for match in re.finditer(r"\b([a-z]+|\d+) environments\b", low):
            token = match.group(1)
            count = NUMBER_WORDS.get(token)
            if count is None and token.isdigit():
                count = int(token)
            if count is not None:
                found.append((sentence, count))
    return found


def unknown_scripts(text: str, known: frozenset[str]) -> frozenset[str]:
    return scripts_named_in(text) - known


def undefined_letters(text: str, known: frozenset[str]) -> frozenset[str]:
    return letters_named_in(text) - known


def wrong_counts(text: str, actual: int) -> list[tuple[str, int]]:
    return [(where, n) for where, n in total_claims_in(text) if n != actual]


# ---------------------------------------------------------------- the tests


@pytest.mark.parametrize("relative", SCANNED)
def test_every_provision_script_named_in_prose_resolves(relative: str) -> None:
    """A `provision/*.sh` in prose is either an environment's script or a
    utility, and either way it is on disk.

    This is the check that makes a rename finish. `provision/env_a.sh` was
    named in three files; renaming the script and the manifest key while
    leaving those behind would have left the reader following a path to a file
    that is not there.
    """
    known = env_scripts() | utility_scripts()
    stale = unknown_scripts((REPO / relative).read_text(encoding="utf-8"), known)
    assert not stale, (
        f"{relative} names {sorted(stale)}, which no environment in "
        f"provision/manifest.toml claims and which is not a utility script in "
        f"provision/. Either the script was renamed and the prose was not, or "
        f"the manifest lost the environment that owned it."
    )


@pytest.mark.parametrize("relative", SCANNED)
def test_every_environment_letter_in_prose_is_defined(relative: str) -> None:
    """`Env C` is the drift this catches: a letter no module defines.

    The legal letters come from `routes.py`, not from a list here, so the check
    cannot be satisfied by editing the test.
    """
    known = defined_letters()
    unknown = undefined_letters(
        (REPO / relative).read_text(encoding="utf-8"), known
    )
    assert not unknown, (
        f"{relative} names Env {', Env '.join(sorted(unknown))}, which "
        f"design_pipeline/routes.py does not define - it defines only "
        f"{sorted(known)}. A letter in prose that no ENV_* constant backs is a "
        f"place the reader cannot install into."
    )


@pytest.mark.parametrize("relative", SCANNED)
def test_stated_environment_counts_match_the_manifest(relative: str) -> None:
    """A sentence that announces a total must announce the manifest's total.

    This is the drift that outlived Env B's fragmentation: `routes.py` said
    "two environments" while the manifest listed eight.
    """
    actual = len(env_keys())
    wrong = wrong_counts((REPO / relative).read_text(encoding="utf-8"), actual)
    assert not wrong, (
        f"{relative} states a total that provision/manifest.toml contradicts "
        f"({actual} environments): "
        + "; ".join(f"{n} in {where!r}" for where, n in wrong)
    )


def test_every_manifest_environment_is_named_in_the_runbook() -> None:
    """The reverse direction: an environment the manifest builds must appear in
    the runbook.

    Without this, the checks above are satisfiable by deleting prose. An
    environment that nothing documents is as hard to rebuild as one documented
    under a name it no longer has.
    """
    runbook = (REPO / "docs/provisioning.md").read_text(encoding="utf-8")
    missing = sorted(key for key in env_keys() if f"`{key}`" not in runbook)
    assert not missing, (
        f"provision/manifest.toml builds {missing}, which docs/provisioning.md "
        f"never names. The runbook is what someone rebuilding the node reads; "
        f"an environment it does not mention does not get rebuilt."
    )


def test_the_printed_name_of_the_import_environment_is_a_manifest_key() -> None:
    """`routes.ENV_A` is what preflight prints in its `env` column and what
    `require`'s error calls the environment, and it names exactly one of them.
    So it has to be that one's manifest key, spelled the manifest's way.

    This is the check that would have caught the original defect from the code
    side rather than the prose side: `ENV_A` stayed "A" through the rename, so
    preflight said `A` while every document said `core`. One environment, two
    names, and nothing failing.

    `ENV_B` has no equivalent check on purpose - it covers five environments,
    so it describes how they are reached instead of naming one.
    """
    keys = env_keys()
    assert routes.ENV_A in keys, (
        f"routes.ENV_A is {routes.ENV_A!r}, which is not an environment in "
        f"provision/manifest.toml ({sorted(keys)}). preflight prints this value "
        f"in its `env` column, so it is the name a reader sees next to docs "
        f"that use the manifest's."
    )
    assert routes.ENV_B not in keys, (
        f"routes.ENV_B is {routes.ENV_B!r}, which is now a manifest key. ENV_B "
        f"tags routes in five different environments, so a value that names one "
        f"of them tells the reader the wrong thing about the other four."
    )


def test_every_manifest_script_is_on_disk() -> None:
    """The manifest's own half of the rename. A `script` pointing at a file
    that is not there is the same drift read from the other side."""
    missing = sorted(
        script for script in env_scripts() if not (REPO / script).is_file()
    )
    assert not missing, (
        f"provision/manifest.toml points at {missing}, which is not on disk."
    )


# ------------------------------------------------- the reconstruction

DRIFT_ENV_C = (
    "Two adapters use a container. Everything else installs into Env A, B or C,\n"
    "so you can skip this section entirely for a scoring-and-gating campaign."
)
"""INSTALL.md as it stood, verbatim."""

DRIFT_ENV_COUNT = (
    "guessing one. The target is two environments and two containers, not ten\n"
    "images:"
)
"""design_pipeline/routes.py as it stood, verbatim."""

DRIFT_STALE_SCRIPT = (
    "needs, so the Env A that `provision/env_a.sh` builds includes it with a CPU"
)
"""INSTALL.md as it stood, verbatim - the reference a rename leaves behind."""


def test_the_checks_catch_the_two_drifts_they_were_written_for() -> None:
    """The point of the suite, asserted rather than claimed.

    Each drifted sentence is fed to the same function the tests above call. If
    a later edit makes a check permissive, this fails even when the tree is
    clean, which is the only way to know the suite still has teeth.
    """
    assert undefined_letters(DRIFT_ENV_C, defined_letters()) == {"C"}, (
        "the letter check no longer flags Env C - most likely it stopped "
        "reading 'Env A, B or C' as three letters"
    )

    assert wrong_counts(DRIFT_ENV_COUNT, len(env_keys())) == [
        ("The target is two environments and two containers, not ten images:", 2)
    ], "the count check no longer flags the two-environment claim"

    known = env_scripts() | utility_scripts()
    assert unknown_scripts(DRIFT_STALE_SCRIPT, known) == {
        "provision/env_a.sh"
    }, "the script check no longer flags the pre-rename path"


def test_the_checks_stay_quiet_on_true_subset_claims() -> None:
    """The false positive the count check is scoped to avoid.

    "Two environments install from a git tree" is true, is about a subset, and
    must stay legal however many environments the manifest grows to. A check
    that fails on correct prose is a check someone deletes.
    """
    subset = (
        "Two environments install from a git tree and **neither commit is "
        "known**."
    )
    assert wrong_counts(subset, len(env_keys())) == []
