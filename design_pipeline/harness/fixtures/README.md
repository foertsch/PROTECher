# Fixture systems

These are **not science**. They exist so that `run_systems` and `report` run
end to end in CI, offline, with no GPU, no container runtime and no network.

* `fixture_a.pdb`, `fixture_b.pdb` — CA-only helix bundles written by
  `design_pipeline.adapters._pdb.write_ca_pdb`. Real chain ids and real residue
  names, so chain lookup, sequence extraction and the CA–CA distance filter all
  work on them; no real fold, no real interface.
* `fixture_a.csv`, `fixture_b.csv` — invented `measured_affinity` values,
  pseudo-random but deterministic, on a pKD-like scale where higher is better
  (hence `higher_is_better=True` on both systems). `fixture_a` includes one
  double mutant, because published sets do.

A Spearman computed against these numbers measures nothing. What it checks is
that the loop, the `only=` validation, the join and the per-scorer breakout
work — and that a scorer with no skill reports a correlation near zero rather
than a crash.

The real systems in `systems.py` download their structure from RCSB and read
their measurements from the ground-truth directory
(`$DESIGN_PIPELINE_GROUND_TRUTH`, default `data/ground_truth/`). Get those from
FLAb, AbBiBench or AB-Bind; they are not redistributed here.
