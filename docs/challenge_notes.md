# Challenge notes

Status as of 2026-10-26, corrected after the research in `research/`. Everything under "Likely challenge" is inference, not announced. Start with `research/SUMMARY.md`.

## What the event page says

- Breast Cancer Awareness Hackathon, 31 Oct 2026, 9:00–19:00, Zürich.
- Hosted by HerCode, part of Swiss {ai} Weeks. Registration needs host approval, so check that all five of us are approved.
- Includes "a protein design challenge by Adaptyv".
- Not stated: target, judging criteria, prizes, schedule, compute, rules, whether wet-lab testing happens.

Source: https://luma.com/4ogriih0 (re-read 2026-10-24).

## Precedents for the format

Details and sources are in `research/01`, `04` and `05`. The metric changed between events, so we cannot assume one.

- **TREM2 hackathon, 28 Feb 2026 (muni/BioArena x Adaptyv, one day, 16 teams):** at most 10 designs per team, ranked by Boltz-2 ipSAE (ipSAE_d0res, PAE and distance cutoff 15). 141 designs received, the top 100 went to the lab, 89 expressed and 37 bound. This is the closest precedent for a one-day event.
- **EGFR round 1 (2024):** 726 designs submitted, 201 tested (top 100 by AF2-Multimer iPAE plus 101 picked by Adaptyv), 5 binders (2.5%).
- **EGFR round 2 (2024):** 1,131 submitted, 400 tested (100 ranked by a rank average of ESM2 likelihood, iPAE and ipTM, plus 300 hand-picked), 53 binders (13%). Best optimised binder 1.21 nM. Best de novo binder 82 nM in the paper and 91.5 nM in Adaptyv's blog post, so check before quoting.
- **Anthropic x Adaptyv competition (28 Sep to 8 Nov 2026):** single chain, 10 to 250 residues, up to 20 designs per participant. Challenge 5 opens 26 Oct. Whether it relates to this event is unknown.
- **Last term's HerCode event (Future of Health, May 2026, ETH):** the protein track was an FGFR2 track run by Lemna Bio, with no public results, no Adaptyv involvement and no reported wet-lab step. It is not a precedent for the Adaptyv format.
- **Removed claim:** "Adaptyv picked 1 to 5 designs per designer" appears only in a secondary write-up and in no primary source. What the primary sources support is that every round-2 designer had at least one design tested.

## Likely challenge

De novo binder design against a breast-cancer-relevant extracellular receptor. The research ranks the candidates in tiers (see `research/02_target_candidates.md` and the table in `SUMMARY.md`):

- **Likely:** HER2 (ERBB2). Trastuzumab binds domain IV and pertuzumab binds domain II. PDB 1N8Z (trastuzumab Fab complex, 2.52 A) and 1S78 (pertuzumab complex, 3.25 A) were verified on RCSB. HER2 numbering in 1N8Z starts from the mature protein (UniProt minus 22).
- **Possible:** EGFR (Adaptyv's best-characterised target) and TROP2 (9PI9, 1.56 A, an approved triple-negative breast cancer drug target).
- **Long tail:** Nectin-4, HER3, PD-L1 (weak breast cancer link) and FGFR2 (last term's track target).
- Probable constraints: a sequence length limit, a cap on submissions, a target structure supplied by the organisers.

## Questions for the organisers at kickoff

The full list of nine is at the end of `research/SUMMARY.md`. The ones that change what we do in the first hour:

1. Which target, and which structure or construct is supplied?
2. What is the ranking metric, computed with which model and settings, and is a scoring script provided?
3. How many designs per team, and how many go to the lab? Is testing guaranteed?
4. Length limit, and does it include the 64-residue tag Adaptyv appends?
5. Which tools are allowed, and what compute is provided?
6. How is the jury scored (metric, novelty, biology, pitch)?

## Before 31 Oct

- Kira: pipeline running end to end on a public target (HER2 and EGFR), with the ipSAE scoring scripted.
- Everyone: read `research/SUMMARY.md` and `filtering_strategy.md`, and send comments to Arion.
- Check what GPUs the team can use on the day. File 03's compute budget depends on it.
- Download the public calibration data (see `filtering_strategy.md` section 5).
