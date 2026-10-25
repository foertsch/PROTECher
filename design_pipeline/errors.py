"""Errors that say which thing is missing and how to get it."""

from __future__ import annotations


class PipelineError(Exception):
    """Base class so callers can catch everything this package raises."""


class UnknownModelError(PipelineError):
    """model_name is not in the registry."""


class MissingRuntimeError(PipelineError):
    """No container runtime (apptainer or docker) on PATH."""


class MissingDependencyError(PipelineError):
    """An optional model dependency is not installed."""


class MissingLigandError(PipelineError):
    """A ligand-conditioned model was asked to run on a structure with no ligand."""


class RunFailedError(PipelineError):
    """The model exited non-zero, or produced no output files."""
