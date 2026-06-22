"""holinbench — benchmarking universal, topology-specific, and context-aware
signatures for bacteriophage holin prediction.

This package is deliberately conservative: it is built to *test* whether a
universal holin sequence signature exists, not to assume one. See PROJECT_SPEC.md.
"""

__version__ = "0.1.0"

DATASET_CATEGORIES = ("gold", "weak", "hard_negative", "unknown")

__all__ = ["__version__", "DATASET_CATEGORIES"]
