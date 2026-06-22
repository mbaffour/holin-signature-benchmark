"""Input/output: read the four dataset classes and auxiliary tables.

All datasets are normalized to a common long-form frame with at least:
    protein_id, sequence, dataset_category
plus whatever metadata the source provided. Sequences come from the CSV
``sequence`` column, or from ``proteins.faa`` for the unknown candidates.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import Config
from .utils import clean_sequence, read_fasta, log

# Expected columns per the spec (used for friendly validation messages).
GOLD_COLUMNS = ["protein_id", "protein_name", "phage_name", "host", "accession",
                "sequence", "evidence_type", "citation", "holin_type",
                "family_label", "notes"]
WEAK_COLUMNS = ["protein_id", "protein_name", "phage_name", "accession",
                "sequence", "annotation", "source_database", "notes"]
NEG_COLUMNS = ["protein_id", "protein_name", "source", "accession", "sequence",
               "negative_type", "annotation", "notes"]
CONTEXT_COLUMNS = ["contig_id", "gene_id", "start", "end", "strand",
                   "protein_id", "product", "sequence"]
TOPOLOGY_COLUMNS = ["protein_id", "length", "tmd_count", "topology", "n_region",
                    "c_region", "signal_peptide", "sar_like", "tool"]


def _read_csv(path: Path, expected: list[str], label: str) -> pd.DataFrame:
    if not path.exists():
        log.warning("%s file not found: %s (skipping)", label, path)
        return pd.DataFrame(columns=expected)
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [c for c in expected if c not in df.columns]
    if missing:
        log.warning("%s file %s is missing columns: %s", label, path.name, missing)
    return df


def load_gold(cfg: Config) -> pd.DataFrame:
    df = _read_csv(cfg.path("gold_holins"), GOLD_COLUMNS, "gold")
    df["dataset_category"] = "gold"
    return df


def load_weak(cfg: Config) -> pd.DataFrame:
    df = _read_csv(cfg.path("weak_holins"), WEAK_COLUMNS, "weak")
    df["dataset_category"] = "weak"
    return df


def load_hard_negatives(cfg: Config) -> pd.DataFrame:
    df = _read_csv(cfg.path("hard_negatives"), NEG_COLUMNS, "hard_negative")
    df["dataset_category"] = "hard_negative"
    return df


def load_unknown(cfg: Config) -> pd.DataFrame:
    path = cfg.path("proteins_faa")
    if not path.exists():
        log.warning("proteins.faa not found: %s (skipping)", path)
        return pd.DataFrame(columns=["protein_id", "sequence", "dataset_category"])
    recs = read_fasta(path)
    df = pd.DataFrame(recs, columns=["protein_id", "sequence"])
    df["dataset_category"] = "unknown"
    df["protein_name"] = "unknown candidate"
    return df


def load_all_datasets(cfg: Config) -> pd.DataFrame:
    """Concatenate the four classes into one long frame (raw, pre-validation)."""
    frames = [load_gold(cfg), load_weak(cfg), load_hard_negatives(cfg), load_unknown(cfg)]
    df = pd.concat(frames, ignore_index=True, sort=False)
    # Normalize sequence whitespace/case but keep raw for QC; validate.py cleans.
    df["sequence"] = df["sequence"].fillna("").map(clean_sequence)
    return df


def load_topology(cfg: Config) -> pd.DataFrame:
    path = cfg.path("topology_tsv")
    if not path.exists():
        log.warning("topology_predictions.tsv not found: %s", path)
        return pd.DataFrame(columns=TOPOLOGY_COLUMNS)
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def load_context(cfg: Config) -> pd.DataFrame:
    path = cfg.path("context_tsv")
    if not path.exists():
        log.warning("context.tsv not found: %s", path)
        return pd.DataFrame(columns=CONTEXT_COLUMNS)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for col in ("start", "end"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_lysis_terms(cfg: Config) -> dict[str, list[str]]:
    """Parse lysis_context_terms.txt -> {category: [terms]} (lowercased)."""
    path = cfg.path("lysis_terms")
    cats: dict[str, list[str]] = {}
    if not path.exists():
        log.warning("lysis_context_terms.txt not found: %s", path)
        return cats
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "\t" in line:
            cat, term = line.split("\t", 1)
        else:
            cat, term = "lysis", line
        cats.setdefault(cat.strip().lower(), []).append(term.strip().lower())
    return cats


def write_table(df: pd.DataFrame, path: Path, sep: str = "\t") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=sep, index=False)
    return path
