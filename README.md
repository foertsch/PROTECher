# PROTECher

Team repo for the Breast Cancer Awareness Hackathon (HerCode, part of Swiss {ai} Weeks), Zürich, 31 Oct 2026, 9:00–19:00. We are aiming at the Adaptyv protein design challenge.

Event page: https://luma.com/4ogriih0

## Team

| Name | Background | Proposed focus |
|---|---|---|
| Arion Förtsch | Bioinformatics, R, cluster work | Filtering logic, repo, submission |
| Kira | ML engineer (Berlin), language and graph models, protein track at the last Future of Health hackathon | Generation pipeline |
| Maria | MSc medical physics and computational science, UZH | Target biology, structure prep, epitope choice |
| Iason | MSc physics, ETH; simulations and computer vision | Physics-based rescoring, compute setup |
| Alixe | Life sciences engineering, EPFL; biosensing and neuroscience | Assay readouts (BLI/SPR), literature on known binders |

Focus areas are proposals, not assignments.

## Layout

```
docs/        challenge notes and filtering strategy
pipeline/    candidate generation (Kira)
filters/     scoring and filtering code
notebooks/   exploration
data/        local only, git-ignored
```

## Working rules

- Branch, open a PR, get one review. Keep `main` runnable.
- No data, structures over a few MB, or model weights in git. Share those through the team chat or a shared drive.
- Update `docs/challenge_notes.md` as soon as the organisers announce the target, metric and limits.

## Status

Challenge details are unknown until the event. Preparation notes are in [docs/challenge_notes.md](docs/challenge_notes.md) and [docs/filtering_strategy.md](docs/filtering_strategy.md).
