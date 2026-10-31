"""AntiFold as a scorer: candidate variants scored with AntiFold's own functions. See scorer.py."""

from .scorer import (
    SCORES,
    ComplexScorer,
    ScoreError,
    main,
    read_candidates,
    run_antifold,
    score_with_complex,
    score_with_table,
    write_scores,
)

__all__ = [
    "SCORES",
    "ComplexScorer",
    "ScoreError",
    "main",
    "read_candidates",
    "run_antifold",
    "score_with_complex",
    "score_with_table",
    "write_scores",
]
