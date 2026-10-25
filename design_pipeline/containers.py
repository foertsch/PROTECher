"""Container runtime detection and `docker run` / `apptainer run` assembly.

One adapter class covers every containerised model. What differs between models
is data — image, argument style, which template carries the input, how to read
the scores back — not a subclass.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Literal, Mapping, Sequence

from .designs import Design
from .errors import MissingRuntimeError, RunFailedError

Runtime = Literal["apptainer", "docker"]
ArgStyle = Literal["hydra", "flag"]

IN_MOUNT = "/mnt/in"
OUT_MOUNT = "/mnt/out"

Collect = Callable[[Path, Path], list[tuple[Path, dict[str, float]]]]
"""(produced_file, run_dir) -> [(design_path, scores)]. Lets an adapter parse
its own output format, and split one file into several designs."""


def detect_runtime() -> Runtime:
    """Prefer apptainer: on shared clusters it is the one we are allowed to use."""
    for runtime in ("apptainer", "docker"):
        if shutil.which(runtime):
            return runtime  # type: ignore[return-value]
    raise MissingRuntimeError(
        "No container runtime found on PATH. design_pipeline needs either apptainer "
        "(preferred: https://apptainer.org/docs/admin/main/installation.html) "
        "or docker (https://docs.docker.com/get-docker/). On a cluster, try "
        "`module load apptainer`. To develop without one, use the dummy "
        "models (model_name='dummy_backbones' and friends)."
    )


def format_kwargs(
    kwargs: Mapping[str, Any],
    style: ArgStyle,
    aliases: Mapping[str, str],
    dash_flags: bool = False,
) -> list[str]:
    """Map Python kwargs to the model's CLI flags.

    Unknown keys pass straight through, so an unanticipated flag never requires
    editing this package. `aliases` only renames the ergonomic names.
    `dash_flags` is for models whose options are dashed (`--host-url`) rather
    than underscored (`--pdb_path`).
    """
    args: list[str] = []
    for key, value in kwargs.items():
        if value is None:
            continue
        name = aliases.get(key, key)
        if style == "hydra":
            args.append(f"{name}={_scalar(value)}")
            continue
        if dash_flags:
            name = name.replace("_", "-")
        flag = ("--" if len(name) > 1 else "-") + name
        if value is True:
            args.append(flag)  # bare switch
        else:
            args += [flag, _scalar(value)]
    return args


def _scalar(value: Any) -> str:
    if isinstance(value, (list, tuple)):
        # Hydra list syntax; argparse-style models that want repeats use extra_args.
        return "[" + ",".join(str(v) for v in value) + "]"
    return str(value)


def stream(cmd: Sequence[str]) -> None:
    """Run a command, inheriting stdout/stderr so a long job stays visible.

    Deliberately not captured: a tee into a log file is the next thing to add,
    and nothing needs it yet.
    """
    print("+ " + " ".join(str(c) for c in cmd), file=sys.stderr, flush=True)
    code = subprocess.run([str(c) for c in cmd], check=False).returncode
    if code != 0:
        raise RunFailedError(f"exit code {code} from: {' '.join(str(c) for c in cmd)}")


def default_collect(path: Path, run_dir: Path) -> list[tuple[Path, dict[str, float]]]:
    return [(path, {})]


@dataclass(frozen=True)
class ContainerAdapter:
    """Runs a model image. Satisfies the ModelAdapter protocol."""

    image: str
    output_globs: tuple[str, ...]
    input_template: tuple[str, ...] = ()
    """How an input path reaches the model, e.g. ``("--pdb_path", "{path}")``."""
    output_template: tuple[str, ...] = ()
    fixed_args: tuple[str, ...] = ()
    command: tuple[str, ...] = ()
    """Prefix before the model's own args; empty when the image has an ENTRYPOINT."""
    arg_style: ArgStyle = "flag"
    aliases: Mapping[str, str] = field(default_factory=dict)
    dash_flags: bool = False
    gpu: bool = True
    per_input: bool = True
    """One invocation per input design, which is also what makes parent_id
    unambiguous. Set False only for a model that truly batches."""
    collect: Collect = default_collect
    path_kwargs: tuple[str, ...] = ()
    """Kwargs whose value is a host path: bind-mounted and rewritten."""

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        extra_args: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        runtime = detect_runtime()
        # No inputs means a fresh run driven entirely by kwargs (a target PDB
        # passed as input_pdb, say), so there is still one invocation to make.
        batches: list[list[Design]] = (
            [[d] for d in inputs] if self.per_input else [list(inputs)]
        )
        if not inputs:
            batches = [[]]

        designs: list[Design] = []
        for batch in batches:
            before = self._existing(output_dir)
            stream(self._command(runtime, batch, output_dir, extra_args or (), kwargs))
            produced = [p for p in self._existing(output_dir) if p not in before]
            # One input per invocation, so the parent is unambiguous; a fresh
            # run from a target path has no parent design.
            parent = batch[0] if len(batch) == 1 else None
            designs += self._designs(produced, output_dir, parent)

        if not designs:
            raise RunFailedError(
                f"{self.image} exited 0 but wrote nothing matching "
                f"{self.output_globs} under {output_dir}"
            )
        return designs

    def _existing(self, output_dir: Path) -> set[Path]:
        found: set[Path] = set()
        for pattern in self.output_globs:
            found |= set(output_dir.rglob(pattern))
        return found

    def _designs(
        self, produced: Sequence[Path], run_dir: Path, parent: Design | None
    ) -> list[Design]:
        designs: list[Design] = []
        for source in sorted(produced):
            for path, scores in self.collect(source, run_dir):
                designs.append(
                    Design(
                        id=make_id(path, run_dir),
                        path=path,
                        parent_id=None if parent is None else parent.id,
                        producer=self.image.split("/")[-1],
                        scores=scores,
                    )
                )
        return designs

    def _command(
        self,
        runtime: Runtime,
        inputs: Sequence[Design],
        output_dir: Path,
        extra_args: Sequence[str],
        kwargs: Mapping[str, Any],
    ) -> list[str]:
        host_paths = [d.path for d in inputs]
        host_paths += [
            Path(str(kwargs[k])) for k in self.path_kwargs if kwargs.get(k) is not None
        ]
        parents = sorted({p.resolve().parent for p in host_paths})
        mounts = {parent: f"{IN_MOUNT}{i}" for i, parent in enumerate(parents)}
        mounts[output_dir.resolve()] = OUT_MOUNT

        def in_container(path: Path) -> str:
            return f"{mounts[path.resolve().parent]}/{path.name}"

        cmd: list[str] = [runtime, "run"]
        if runtime == "docker":
            cmd.append("--rm")
            if self.gpu:
                cmd += ["--gpus", "all"]
            for host, dest in mounts.items():
                cmd += ["-v", f"{host.as_posix()}:{dest}"]
            cmd.append(self.image)
        else:
            if self.gpu:
                cmd.append("--nv")
            for host, dest in mounts.items():
                cmd += ["--bind", f"{host.as_posix()}:{dest}"]
            # Apptainer resolves a bare name against its own library; docker:// is explicit.
            cmd.append(f"docker://{self.image}")

        cmd += list(self.command)
        for design in inputs:
            cmd += [t.format(path=in_container(design.path)) for t in self.input_template]
        cmd += [t.format(path=OUT_MOUNT) for t in self.output_template]
        cmd += list(self.fixed_args)

        rewritten = dict(kwargs)
        for key in self.path_kwargs:
            if rewritten.get(key) is not None:
                rewritten[key] = in_container(Path(str(rewritten[key])))
        cmd += format_kwargs(rewritten, self.arg_style, self.aliases, self.dash_flags)
        cmd += list(extra_args)
        return cmd


def make_id(path: Path, run_dir: Path) -> str:
    """Stable across a resume (the run dir is the cache key) and unique across
    rounds (a different config or seed is a different key)."""
    return f"{path.stem}-{run_dir.name[:8]}"
