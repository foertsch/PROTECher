# CPU physics and geometry screening

This adds structure-specific geometry checks to Iason's developability report.
It uses the existing NumPy/Biotite environment, without GPU inference, downloads,
atom rebuilding or energy minimization. It reports evidence and review flags,
not experimentally validated expression, stability or affinity predictions.

## Run now

From the repository root:

```bash
# Parent only. The existing local 1N8Z PDB is required.
.venv/bin/python hackday/06_physics_screen.py

# Include candidate IDs, explicitly recording which lack their own structures.
.venv/bin/python hackday/06_physics_screen.py \
  --input data/designs/variants.fasta --out-dir runs/physics_parent

# Attach the observations to the existing sequence/developability report.
# This writes a separate snapshot; the original review CSV remains unchanged.
.venv/bin/python hackday/05_developability.py \
  --input data/designs/variants.fasta \
  --physics runs/physics_parent/physics_review.csv \
  --out runs/physics_parent/developability_with_physics.csv
```

FASTA or `name,sequence` CSV is accepted. The team's headerless submission CSV
requires `--csv-no-header`. All CLI paths are relative to the repository root
unless absolute; structure paths inside a manifest are relative to that manifest.
The parent sequence is fixed to the organiser's 242-aa scFv.

The parent configuration is `hackday/physics_parent.json`: VH B1-120, VL A1-107,
and available HER2 C480-607 coordinates in 1N8Z. This uses the target interval
already used in the repository's FoldX scan. It does not assert that this is
the team's final cofolding construct. It removes Fab constant domains and does
not invent coordinates for the scFv linker or missing HER2 residues.

## What is implemented

| Check | Implementation | Interpretation |
|---|---|---|
| Structure identity | Exact amino-acid agreement at every mapped scFv position, plus manifest sequence SHA-256 | Rejects stale/wrong structures; never relabels the parent as a mutant |
| Missing atoms/coverage | Standard-residue heavy-atom templates; missing binder residues are invalid, missing side-chain atoms and target gaps are recorded | Measurements from incomplete structures can undercount contacts |
| Steric clashes | Element-specific van der Waals radii; nonbonded overlaps, separated into binder and binder-target pairs | Custom heavy-atom flag, **not** MolProbity clashscore |
| Backbone continuity | C-to-next-N distances within each chain | Flags gross peptide geometry problems |
| Canonical disulfides | S-S distances and CB-S-S-CB torsions at full-scFv pairs 22/96 and 158/223 | Checks both structural bridges, not just cysteine letters |
| Charged contacts | Minimum distance between charged side-chain atoms, one count per residue pair | Opposite-charge salt-bridge candidates and separate like-charge contacts; not energies |
| Polar contacts | Chemically restricted donor/acceptor heavy atoms within a distance window | Distance candidates, **not** confirmed hydrogen bonds |
| Hydrogen-bond geometry | Explicit donor-H-acceptor geometry, if H/D atoms are supplied | Missing/partial hydrogens are labelled; the parent has none |
| Surface exposure | Free-binder SASA, hydrophobic side-chain SASA, and binder area buried by the selected antigen | Descriptive areas in square angstroms; not aggregation or expression probabilities |

Individual atom pairs, chain/residue identities, full-scFv positions, distances,
angles and reasons are saved in `physics_events.csv`. Charge counts and surface
areas are descriptive; there is no rule that more salt bridges means better binding.

## Explicit screening windows

Values are recorded in `physics_metadata.json` and implemented in
`design_pipeline.physics.PARAMETERS`. They are review settings, not new challenge
requirements or calibrated biological rejection thresholds.

| Setting | Value |
|---|---|
| Heavy-atom clash | Adjusted overlap >= 0.6 Å |
| Polar donor/acceptor overlap allowance | 0.4 Å subtracted from overlap |
| Topological exclusions | Atoms separated by up to three covalent bonds (1-2, 1-3, 1-4) |
| Salt-bridge / like-charge contact distance | <= 4.0 Å between selected side-chain atoms |
| Polar heavy-atom contact distance | <= 3.5 Å |
| Explicit H-to-acceptor distance and D-H-A angle | <= 2.5 Å and >= 120 degrees |
| Consecutive peptide C-N distance | 1.1-1.6 Å |
| Canonical S-S distance | 1.8-2.3 Å |
| Absolute CB-S-S-CB torsion | 45-135 degrees |
| SASA | 1.4 Å probe, 1000 points, single-element radii |

Intramolecular bonds come from Biotite's bundled chemical dictionary, peptide
connectivity and geometrically plausible disulfides. No chemical-dictionary
download is performed. Distance-based charged groups assume Asp/Glu negative
and Lys/Arg positive; histidine, termini, ions, water and glycans are excluded
from that analysis. Hydrogen-bond atom typing is deliberately restricted and
does not optimize protonation or infer absent hydrogens. Hydrogen-bond counts
with partial hydrogen coverage are observed counts, not complete counts.

Sources and distinction from this implementation:

- [Biotite van der Waals radii](https://www.biotite-python.org/latest/apidoc/biotite.structure.info.vdw_radius_single.html)
  supplies the element radii.
- [MolProbity](https://molprobity.lbl.gov/help/) performs all-atom contact validation.
  Our heavy-atom approximation and thresholds are not that method.
- [VMD salt bridges](https://www.ks.uiuc.edu/Research/vmd/plugins/saltbr/)
  illustrates that distance definitions are configurable (its default is 3.2 Å).
  Our 4.0 Å window is explicitly a candidate-contact screen.
- [Biotite hydrogen bonds](https://www.biotite-python.org/latest/apidoc/biotite.structure.hbond.html)
  documents the H-acceptor/angle criterion. Our implementation additionally
  restricts donors/acceptors by residue and excludes local bonded neighbours.
- [Biotite disulfide example](https://www.biotite-python.org/latest/examples/gallery/structure/contacts/disulfide_bonds.html)
  illustrates the distance/torsion checks. Our windows are broader review windows.

## Parent observations on 31 October 2026

Input PDB SHA-256:
`ed934d36503335c6764b08880f9e50933d99eb13b9bf145339cfcaee150eb21e`.

| Observation | Result |
|---|---|
| Mapped binder residues | 227, Fv only; linker omitted |
| VH disulfide, scFv 22/96 | 2.039 Å, torsion -84.73 degrees; within screening window |
| VL disulfide, scFv 158/223 | 2.090 Å, torsion 69.48 degrees; within screening window |
| Backbone C-N outliers in binder | 0 |
| Binder heavy-atom clash flags | 2 |
| Binder-target clash flags | 0 among the available selected atoms |
| Salt-bridge candidates | 10 within binder, 3 at selected binder-target interface |
| Like-charge contacts within 4 Å | 0 in binder and interface |
| Polar distance contacts | 308 within binder, 7 at interface; not hydrogen-bond counts |
| Free-binder SASA | About 10,270.65 Å² |
| Hydrophobic side-chain SASA | About 1,594.84 Å² |
| Binder area buried by target | About 696.16 Å² |
| Missing binder heavy atoms | B:30 LYS lacks CG, CD, CE, NZ |
| Missing target coordinates | C581-590 within the requested interval |
| Explicit hydrogen-bond check | Unavailable: the PDB has no hydrogens |

The two close-contact flags are A1 ASP O / A3 GLN NE2 (2.046 Å) and
B99 TRP O / B100 GLY C (2.415 Å). These are inspectable coordinate warnings,
not evidence that the known parent fails to fold. Do not tune thresholds merely
to force this parent to pass. The overall parent status is `review_required`,
with incomplete checks recorded separately. Missing K30 atoms are especially
relevant when investigating K30 mutations.

The initial run also retained all 15 provisional local variants as
`not_evaluated_no_structure`, with blank physical measurements. No candidate
was scored by borrowing parent coordinates. When attached to developability,
unassessed physical checks make the combined assessment `needs_review` while
leaving sequence-rule results unchanged. They are not biological failures.

## Adding candidate structures later

Supply the exact sequences and a manifest following `physics_parent.json`.
Keep a parent record. Each structure record needs:

- `design`: exact ID in the sequence file (`PARENT` for the reference).
- `sequence_sha256`: SHA-256 of the normalized uppercase sequence; available
  in the developability report. A renamed parent PDB will fail sequence checking.
- `path`: PDB, CIF/mmCIF or BCIF; `model` optionally selects a 1-based model (default 1).
- `scope`: `fv_only` (positions 1-120 and 136-242) or `full_scfv` (all 242).
- `binder_segments`: each has `chain`, `start`, `end`, `sequence_start`.
  These map author chain/residue identifiers to 1-based full-scFv positions.
- `antigen_segments`: optional list of `chain`, `start`, `end`. Target and
  binder must have separate chains. Use the actual output numbering, not the
  crystal numbering by assumption.

For example, only if the candidate structure really has the complete scFv
as chain X residues 1-242, its binder segment is
`{"chain":"X","start":1,"end":242,"sequence_start":1}`.
Insertion codes and length-changing designs require explicit remapping or a
future alignment adapter; the current code rejects them instead of guessing.

Run with `--manifest path/to/candidate_structures.json --input candidates.fasta`.
Only one structure per design is accepted in a run. Use a separate manifest and
output directory for each seed; seed aggregation and calibrated vetoes remain
future work. Highest-occupancy alternate conformers are selected explicitly.

Parent deltas require matching binder scope. Interface deltas additionally
require identical selected antigen sequence/coverage; gaps enter the signature.
Missing side-chain atoms suppress affected deltas. Partial-hydrogen bond counts
are not subtracted as if complete. Fv-only crystal measurements must not serve
as a full-scFv/linker baseline. Re-run the predicted parent in the same format
as the candidates for that comparison.

## Outputs and remaining dependencies

`physics_review.csv` contains one row per design; `physics_events.csv` explains
the measured contacts; `physics_report.md` is a readable overview;
`physics_metadata.json` records input/script hashes, versions, selection coverage
and all settings. No structures or candidate sequences are uploaded by these
scripts. Output/input path collisions are rejected. Invalid supplied structures
are retained with blank metrics and cause a nonzero CLI exit.

| Needed from the team | Enables | When |
|---|---|---|
| Stable candidate IDs and exact untagged sequences | Identity checks and sequence screening | Available now, even before GPU runs |
| Actual mutant structures, chain/residue mappings and matching parent | Per-candidate geometry; full-scFv comparison including linker | After side-chain rebuilding or structure prediction; CPU rebuilding is also possible |
| HER2 construct sequence/boundaries, included chains, model/seed IDs and preparation details | Comparable interfaces and later seed summaries | With each structure batch |
| Prepared hydrogens/protonation plus assumed pH; buffer/salt context | Better chemical interpretation and explicit hydrogen-bond checks | When prepared structures exist |
| FoldX access through the team's existing setup, or complete candidate results with parent references | Folding and binding energy differences, separately | After complete mutant side chains and an agreed preparation protocol exist |
| Agreed metrics and calibrated tolerances | Biological filtering decisions beyond coordinate flags | After reviewing parent controls and relevant validation data |

FoldX `BuildModel`/`AnalyseComplex` results, OpenMM minimization, buried
unsatisfied-polar energetics, ion coordination, molecular dynamics, and seed
aggregation are **not implemented or run here**. The existing FoldX parent
alanine scan is not a full-candidate stability calculation. Do not sum its
single-mutation energies to claim a multi-mutant energy.

## Verification

```bash
.venv/bin/python -m pytest -q tests/test_physics.py \
  tests/test_hackday_developability.py tests/test_gates.py \
  tests/test_scfv.py tests/test_track_a_export.py
```

The checks cover chemical/topological exclusions, clashes, charge-pair counting,
hydrogen presence and angles, disulfide geometry, missing atoms, PDB/CIF/BCIF,
sequence mapping, stale IDs/hashes, report integration and overwrite protection.
Synthetic fixtures are code tests only. The optional local-parent test performs
no download and skips when 1N8Z is absent.
