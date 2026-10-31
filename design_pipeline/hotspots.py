"""FoldX alanine-scan hotspots of the trastuzumab scFv (1N8Z), as importable data.

Source: docs/track_a_filtering.md, "Precomputed today": an alanine scan of the 1N8Z scFv with HER2 residues 480 to 607 in
FoldX 5.1 (`examples/foldx_alanine_scan.py`). This module holds the numbers quoted there and nothing newer.

`ddg_bind` is the change in binding energy when the position is set to alanine, in kcal/mol; positive means the alanine
weakens binding. `runs` is how many FoldX runs the value is a mean of; `None` where the doc does not say. A position is a
hotspot when `ddg_bind` is at least 1.0. Alanine removes the side chain beyond its first carbon, so a hotspot is a position
whose side chain carries binding energy, and the plan keeps designs away from those positions. FoldX agreement with
measurements is weak to moderate, so treat every row as a flag and not as a measurement.

Numbering: `position` counts from the start of the domain (VH 1 to 120, VL 1 to 107), which is also the 1N8Z chain
position, with the heavy chain as chain B and the light chain as chain A. `scfv_position` is the position in the 242-residue
scFv (VH, 15-residue linker, VL): the same as `position` for VH, and `position` + 135 for VL. `design_pipeline.scfv` converts
between the two and `tests/test_hotspots.py` checks every row against it and against the parent sequence.

    from design_pipeline import hotspots

    avoid = hotspots.positions()               # scFv positions of the ten hotspots
    caution = hotspots.positions("near_miss")  # just under 1.0

FoldX is licensed per user, and its licence text as summarised in docs/track_a_filtering.md may need the licensors' consent to
publish its output. These values are FoldX output: read the agreement before this file leaves a private repository.
"""

from __future__ import annotations

from dataclasses import dataclass

STATUSES = ("hotspot", "near_miss", "excluded", "alanine_improves")
"""`hotspot`: alanine costs at least 1.0 kcal/mol of binding. `near_miss`: just under 1.0, check before changing.
`excluded`: second-shell position kept out of the framework arm (VL R66). `alanine_improves`: alanine is predicted to
strengthen binding, a hypothesis only, and a loop-arm candidate."""

_CHAIN = {"VH": "B", "VL": "A"}
_SCFV_OFFSET = {"VH": 0, "VL": 135}


@dataclass(frozen=True)
class Hotspot:
    domain: str
    position: int
    wt: str
    ddg_bind: float | None
    runs: int | None
    status: str

    @property
    def chain(self) -> str:
        """The 1N8Z chain: B for the heavy domain, A for the light domain."""
        return _CHAIN[self.domain]

    @property
    def scfv_position(self) -> int:
        return _SCFV_OFFSET[self.domain] + self.position

    @property
    def crystal_id(self) -> str:
        """Chain, wild type and position as a variant id starts, e.g. `B:Y33`."""
        return f"{self.chain}:{self.wt}{self.position}"


HOTSPOTS: tuple[Hotspot, ...] = (
    # Ten positions that lose 1.0 kcal/mol or more on alanine, all within 3.8 A of HER2. VL N30 and VH W99 stay at 3 runs
    # (spread 0.26 and 0.02); the other eight are 10-run means.
    Hotspot("VH", 33, "Y", 2.04, 10, "hotspot"),
    Hotspot("VH", 50, "R", 1.98, 10, "hotspot"),
    Hotspot("VH", 57, "Y", 1.36, 10, "hotspot"),
    Hotspot("VH", 99, "W", 1.97, 3, "hotspot"),
    Hotspot("VH", 105, "Y", 1.13, 10, "hotspot"),
    Hotspot("VL", 30, "N", 2.30, 3, "hotspot"),
    Hotspot("VL", 31, "T", 1.34, 10, "hotspot"),
    Hotspot("VL", 53, "F", 1.09, 10, "hotspot"),
    Hotspot("VL", 91, "H", 1.29, 10, "hotspot"),
    Hotspot("VL", 94, "T", 1.64, 10, "hotspot"),
    # Just under 1.0.
    Hotspot("VH", 59, "R", 0.92, None, "near_miss"),
    Hotspot("VL", 92, "Y", 0.90, None, "near_miss"),
    # Second shell, excluded from the framework arm.
    Hotspot("VL", 66, "R", 0.79, None, "excluded"),
    # Alanine predicted to improve binding; the doc gives no value.
    Hotspot("VL", 28, "D", None, None, "alanine_improves"),
    Hotspot("VL", 50, "S", None, None, "alanine_improves"),
    Hotspot("VH", 103, "G", None, None, "alanine_improves"),
)


def positions(status: str = "hotspot") -> frozenset[int]:
    """The scFv positions with this status."""
    if status not in STATUSES:
        raise ValueError(f"status must be one of {', '.join(STATUSES)}, got {status!r}")
    return frozenset(h.scfv_position for h in HOTSPOTS if h.status == status)
