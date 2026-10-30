"""Per-run reports and a cross-run leaderboard, built from design_pipeline's DONE manifests.

`python -m run_report build RUN_DIR --run-id ID --out reports/ID` writes one report.
`python -m run_report leaderboard reports --out LEADERBOARD.md` compares all reports.
Needs pandas and matplotlib only. It reads the manifests as JSON and imports no model.
"""
