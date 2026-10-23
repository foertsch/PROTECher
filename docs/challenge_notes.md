# Challenge notes

Status as of 2026-10-23. Everything under "Likely challenge" is inference, not announced.

## What the event page says

- Breast Cancer Awareness Hackathon, 31 Oct 2026, 9:00–19:00, Zürich.
- Hosted by HerCode, part of Swiss {ai} Weeks. Registration needs host approval, so check that all five of us are approved.
- Includes "a protein design challenge by Adaptyv".
- Not stated: judging criteria, prizes, schedule, compute, rules.

Source: https://luma.com/4ogriih0

## Precedents: how Adaptyv scored earlier challenges

The metric changed between events, so we cannot assume one.

- EGFR competition: designs ranked by iPAE (AlphaFold-derived interface error). A later account also lists AF2 pTM and ESM2 pseudo-log-likelihood, with the top 100 going to the wet lab. Adaptyv also picked 1–5 designs per designer based on quality and novelty of the approach. ([write-up](https://scour.ing/@minezone/p/https://blog.booleanbiotech.com/adaptyv-binder-design-competition))
- A later hackathon format: up to 10 binders per team, ranked by ipSAE, 100 sent to the wet lab.
- EGFR round 1: 700+ designs submitted, 200 validated, 5 bound (2.5%). ([Adaptyv](https://adaptyvbio.substack.com/p/designing-better-cancer-drugs-using))
- EGFR rounds 1 and 2 together: 1,800+ designs submitted, 601 characterized by BLI. Best optimized binder 1.21 nM, best de novo binder 82 nM. ([Twist/bioRxiv](https://twistbioscience.com/resources/publication/crowdsourced-protein-design-lessons-adaptyv-egfr-binder-competition))
- Last term's event, the HerCode Future of Health hackathon at ETH (May 2026), had a protein design track too. ([SOPHiA GENETICS](https://careers.sophiagenetics.com/blogs/sg-updates-blog/building-the-future-of-health-one-hack-at-a-time))

## Likely challenge

De novo binder design against a breast-cancer-relevant extracellular receptor.

- HER2 (ERBB2) is the most likely target. HER2-positive disease is roughly 15–20% of breast cancers. Trastuzumab binds domain IV, pertuzumab binds domain II (dimerisation arm). Starting structures: PDB 1N8Z, 1S78 (IDs not yet verified).
- Alternatives: EGFR, HER3, TROP2, Nectin-4, PD-L1.
- Probable constraints: a sequence length limit, a cap on submissions, a target structure supplied by the organisers.

## Questions for the organisers at kickoff

1. Which target, and which structure or construct is supplied?
2. What is the ranking metric (ipSAE, ipTM, iPAE, something else)? Is there a second stage?
3. Sequence length limit and number of submissions per team?
4. Is wet-lab validation guaranteed, and for how many designs?
5. Which tools are allowed (BindCraft, RFdiffusion, ProteinMPNN, Boltz, Chai)? Is compute provided?
6. Does the jury score approach and novelty, or only the metric?

## Before 31 Oct

- Kira: environment running, pipeline tested end to end on EGFR.
- Everyone: read up on HER2 epitopes and known HER2 binders (trastuzumab, pertuzumab, affibodies, DARPins).
- Check whether the Adaptyv EGFR competition data is public (see `filtering_strategy.md`, calibration).
