"""Stage 1 — data validation and cleaning.

Reads the four dataset classes, validates sequences, removes invalid records and
exact duplicates, standardizes IDs, computes basic biochemistry, and emits a
clear table of QC flags. Conservative by design: most issues are *flagged*, not
silently dropped. Only clearly invalid sequences are removed.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import io
from .config import Config
from .utils import (charge_counts, hydrophobic_fraction, is_valid_sequence,
                    log, molecular_weight, net_charge, write_fasta)


def _standardize_id(pid: str, idx: int) -> str:
    pid = (pid or "").strip()
    return pid if pid else f"unnamed_{idx:05d}"


def run(cfg: Config) -> dict:
    """Execute Stage 1. Returns a dict with the cleaned frame and flags frame."""
    vcfg = cfg.section("validation")
    min_len = int(vcfg.get("min_length", 30))
    max_len = int(vcfg.get("max_length", 300))
    hard_min = int(vcfg.get("hard_min_length", 10))
    max_amb = float(vcfg.get("max_ambiguous_fraction", 0.10))
    drop_invalid = bool(vcfg.get("drop_invalid", True))
    drop_dups = bool(vcfg.get("drop_exact_duplicates", True))

    raw = io.load_all_datasets(cfg)
    raw["protein_id"] = [_standardize_id(p, i) for i, p in enumerate(raw["protein_id"])]

    flags: list[dict] = []

    def flag(pid, category, flag_name, detail=""):
        flags.append({"protein_id": pid, "dataset_category": category,
                      "flag": flag_name, "detail": detail})

    # -- validity -------------------------------------------------------------
    keep_mask = []
    for _, row in raw.iterrows():
        ok, reason = is_valid_sequence(row["sequence"], max_amb)
        if not ok:
            flag(row["protein_id"], row["dataset_category"], "invalid_sequence", reason)
        if len(row["sequence"]) < hard_min:
            ok = False
            flag(row["protein_id"], row["dataset_category"], "below_hard_min_length",
                 f"len={len(row['sequence'])}")
        keep_mask.append(ok or not drop_invalid)
    clean = raw[pd.Series(keep_mask, index=raw.index)].copy()

    # -- duplicate sequences within / across class ----------------------------
    # Cross-class duplicates are always flagged (data-leakage / labeling risk).
    seq_to_cats: dict[str, set[str]] = {}
    for _, row in clean.iterrows():
        seq_to_cats.setdefault(row["sequence"], set()).add(row["dataset_category"])
    for seq, cats in seq_to_cats.items():
        if len(cats) > 1:
            members = clean[clean["sequence"] == seq]["protein_id"].tolist()
            for pid, cat in zip(members, clean[clean["sequence"] == seq]["dataset_category"]):
                flag(pid, cat, "duplicate_sequence_across_classes",
                     f"classes={sorted(cats)}")

    # Exact duplicates within the same class: keep first, drop/flag the rest.
    if drop_dups:
        before = len(clean)
        dup_mask = clean.duplicated(subset=["dataset_category", "sequence"], keep="first")
        for _, row in clean[dup_mask].iterrows():
            flag(row["protein_id"], row["dataset_category"],
                 "exact_duplicate_within_class_removed", "")
        clean = clean[~dup_mask].copy()
        if before != len(clean):
            log.info("Removed %d exact within-class duplicate(s)", before - len(clean))

    # -- length / ambiguity flags (kept records) ------------------------------
    for _, row in clean.iterrows():
        L = len(row["sequence"])
        if L < min_len:
            flag(row["protein_id"], row["dataset_category"], "short_protein", f"len={L}")
        if L > max_len:
            flag(row["protein_id"], row["dataset_category"], "long_protein", f"len={L}")

    # -- gold-specific evidence QC --------------------------------------------
    gold = clean[clean["dataset_category"] == "gold"]
    for _, row in gold.iterrows():
        if not str(row.get("citation", "")).strip():
            flag(row["protein_id"], "gold", "gold_missing_citation", "")
        if not str(row.get("evidence_type", "")).strip():
            flag(row["protein_id"], "gold", "gold_missing_evidence_type", "")

    # -- per-record biochemistry ----------------------------------------------
    metrics = []
    for _, row in clean.iterrows():
        seq = row["sequence"]
        pos, neg = charge_counts(seq)
        metrics.append({
            "protein_id": row["protein_id"],
            "length": len(seq),
            "molecular_weight": molecular_weight(seq),
            "hydrophobic_fraction": round(hydrophobic_fraction(seq), 4),
            "positive_residues": pos,
            "negative_residues": neg,
            "net_charge": net_charge(seq),
        })
    metrics_df = pd.DataFrame(metrics)
    clean = clean.merge(metrics_df, on="protein_id", how="left")

    flags_df = pd.DataFrame(flags, columns=["protein_id", "dataset_category", "flag", "detail"])

    # -- write outputs --------------------------------------------------------
    tables = cfg.out_dir("tables")
    seqs = cfg.out_dir("sequences")
    clean.to_csv(tables / "cleaned_metadata.tsv", sep="\t", index=False)
    flags_df.to_csv(tables / "qc_flags.tsv", sep="\t", index=False)

    # Cleaned FASTA: one per category + a combined file.
    write_fasta(seqs / "all_clean.faa",
                list(zip(clean["protein_id"], clean["sequence"])))
    for cat, sub in clean.groupby("dataset_category"):
        write_fasta(seqs / f"{cat}.faa", list(zip(sub["protein_id"], sub["sequence"])))

    log.info("Stage 1: %d valid records (%s); %d QC flags",
             len(clean),
             ", ".join(f"{c}={n}" for c, n in clean["dataset_category"].value_counts().items()),
             len(flags_df))

    return {"clean": clean, "flags": flags_df,
            "metadata_path": tables / "cleaned_metadata.tsv",
            "flags_path": tables / "qc_flags.tsv"}
