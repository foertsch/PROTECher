import sys

from .scorer import ScoreError, main

try:
    sys.exit(main())
except ScoreError as exc:
    sys.exit(f"antifold_score: {exc}")
