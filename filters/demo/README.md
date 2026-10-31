# Synthetic comparison demonstration

All sequences, scores, developability labels and thresholds here are invented
software fixtures. The ten-residue sequences are not scFvs, proposed designs,
measurements or biological evidence. Do not submit them or use the example
thresholds as scientific recommendations.

- `candidates.csv`: a toy parent and four toy candidates.
- `sequence_scores.csv`: one synthetic model score per candidate.
- `structure_scores.csv`: synthetic repeats, deliberately missing a repeat for toy_c.
- `developability.csv`: synthetic rule/warning labels, not calculated with Iason's script.
- `comparison.json`: demo-only configuration illustrating gates and diversity selection.

Run the command in `../MARIA_HANDOFF.md`. Expected decisions are one each of
ELIGIBLE, REVIEW, INCOMPLETE and REJECT; only one of two requested shortlist
places can be filled. Real model agreement or affinity cannot be inferred.
