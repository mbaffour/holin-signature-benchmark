"""Transmembrane topology: import predictions, or fall back to a crude estimate.

The preferred path is to IMPORT predictions from an external, trustworthy tool
(DeepTMHMM / Phobius / TMHMM / TOPCONS) supplied as ``topology_predictions.tsv``.
When no import is available, ``builtin_topology`` provides a rough Kyte-Doolittle
estimate that is always tagged ``tool='builtin_kd_estimate'`` so downstream code
and the report can flag it as low-confidence.
"""
from __future__ import annotations

import pandas as pd

from . import io
from .config import Config
from .utils import POSITIVE_AA, NEGATIVE_AA, estimate_tmd_count, log

TOPOLOGY_FIELDS = ["protein_id", "length", "tmd_count", "topology", "n_region",
                   "c_region", "signal_peptide", "sar_like", "tool"]


def import_topology(cfg: Config) -> pd.DataFrame:
    df = io.load_topology(cfg)
    if df.empty:
        return df
    if "tmd_count" in df.columns:
        df["tmd_count"] = pd.to_numeric(df["tmd_count"], errors="coerce")
    if "length" in df.columns:
        df["length"] = pd.to_numeric(df["length"], errors="coerce")
    return df


def builtin_topology(records: list[tuple[str, str]], cfg: Config) -> pd.DataFrame:
    fcfg = cfg.section("features")
    window = int(fcfg.get("kyte_doolittle_window", 19))
    thr = float(fcfg.get("tmd_hydrophobicity_threshold", 1.6))
    sar_max = int(fcfg.get("sar_max_charge_nterm", 1))
    rows = []
    for pid, seq in records:
        ntmd = estimate_tmd_count(seq, window=window, threshold=thr)
        nterm = seq[:30]
        n_charge = sum(1 for c in nterm if c in POSITIVE_AA or c in NEGATIVE_AA)
        rows.append({
            "protein_id": pid, "length": len(seq), "tmd_count": ntmd,
            "topology": "unknown",
            "n_region": "1-30", "c_region": f"{max(1, len(seq) - 30)}-{len(seq)}",
            "signal_peptide": "unknown",
            "sar_like": "yes" if (ntmd >= 1 and n_charge <= sar_max) else "no",
            "tool": "builtin_kd_estimate",
        })
    return pd.DataFrame(rows, columns=TOPOLOGY_FIELDS)


def get_topology(cfg: Config, records: list[tuple[str, str]]) -> pd.DataFrame:
    """Return a topology frame, preferring import, filling gaps with the estimate."""
    source = cfg.dotted("features.topology_source", "import")
    imported = import_topology(cfg) if source == "import" else pd.DataFrame()
    if imported.empty:
        if source == "import":
            log.warning("No imported topology found; using builtin KD estimate "
                        "(low confidence, flagged tool='builtin_kd_estimate').")
        return builtin_topology(records, cfg)
    # Fill any records missing from the import with the builtin estimate.
    have = set(imported["protein_id"])
    missing = [(pid, seq) for pid, seq in records if pid not in have]
    if missing:
        log.info("Topology import covers %d/%d proteins; estimating %d missing.",
                 len(have), len(records), len(missing))
        imported = pd.concat([imported, builtin_topology(missing, cfg)],
                             ignore_index=True, sort=False)
    return imported
