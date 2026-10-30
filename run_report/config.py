"""Thresholds and headline metrics, read from one TOML file the team edits."""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Metric:
    key: str
    label: str
    better: str = "higher"  # "higher" or "lower"
    target: float | None = None  # at or beyond this the design is "ok"
    gate: float | None = None  # worse than this the design "fail"s
    digits: int = 2

    def __post_init__(self) -> None:
        if self.better not in ("higher", "lower"):
            raise ValueError(f"metric {self.key!r}: better must be 'higher' or 'lower'")


@dataclass(frozen=True)
class Config:
    metrics: tuple[Metric, ...]
    x: str
    y: str
    group: str = "meta_group"

    @property
    def primary(self) -> Metric:
        return self.by_key(self.x)

    def by_key(self, key: str) -> Metric:
        for metric in self.metrics:
            if metric.key == key:
                return metric
        raise KeyError(f"no metric {key!r} in the thresholds file")


def load_config(path: Path | str) -> Config:
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)
    metrics = tuple(Metric(**entry) for entry in raw.get("metric", []))
    if not metrics:
        raise ValueError(f"{path} defines no [[metric]] entries")
    run = raw.get("run", {})
    config = Config(
        metrics=metrics,
        x=run.get("x", metrics[0].key),
        y=run.get("y", metrics[min(1, len(metrics) - 1)].key),
        group=run.get("group", "meta_group"),
    )
    config.by_key(config.x), config.by_key(config.y)  # fail early on a typo
    return config
