"""Stage 6 — transmembrane topology and architecture features.

Builds one feature row per protein from the cleaned dataset plus topology
(imported or estimated). These features feed both the architecture sub-score
(scoring.py) and the architecture-only benchmark model (benchmark.py).
"""
from __future__ import annotations

import re

import pandas as pd

from . import topology
from .config import Config
from .utils import (HYDROPHOBIC_AA, NEGATIVE_AA, POSITIVE_AA, charge_counts,
                    hydrophobic_fraction, molecular_weight, net_charge)

# Annotation substrings that indicate a non-holin enzymatic / structural domain.
ENZYMATIC_STRUCTURAL_TERMS = [
    "amidase", "glycosidase", "muramidase", "hydrolase", "lysozyme", "transglycosylase",
    "transporter", "permease", "polymerase", "kinase", "protease", "nuclease",
    "capsid", "tail", "baseplate", "portal", "terminase", "structural",
]


def _terminal_charge(seq: str, n: int) -> dict[str, int]:
    nterm, cterm = seq[:n], seq[-n:] if len(seq) >= n else seq
    np_ = sum(1 for c in nterm if c in POSITIVE_AA)
    nn = sum(1 for c in nterm if c in NEGATIVE_AA)
    cp = sum(1 for c in cterm if c in POSITIVE_AA)
    cn = sum(1 for c in cterm if c in NEGATIVE_AA)
    return {"nterm_pos": np_, "nterm_neg": nn, "cterm_pos": cp, "cterm_neg": cn}


def _low_complexity_fraction(seq: str, k: int = 12) -> float:
    """Fraction of windows dominated (>= 50%) by a single residue. Rough LCR proxy."""
    if len(seq) < k:
        return 0.0
    flagged = 0
    total = len(seq) - k + 1
    for i in range(total):
        w = seq[i:i + k]
        top = max(w.count(c) for c in set(w))
        if top / k >= 0.5:
            flagged += 1
    return round(flagged / total, 4)


def _has_enzymatic_or_structural(annotation: str) -> bool:
    a = (annotation or "").lower()
    return any(t in a for t in ENZYMATIC_STRUCTURAL_TERMS)


def run(cfg: Config, clean: pd.DataFrame) -> pd.DataFrame:
    fcfg = cfg.section("features")
    term_n = int(fcfg.get("terminal_window", 30))

    records = list(zip(clean["protein_id"], clean["sequence"]))
    topo = topology.get_topology(cfg, records)
    topo_idx = topo.set_index("protein_id") if not topo.empty else pd.DataFrame()

    # Gather any annotation text available per protein (weak/negative carry it).
    ann_col = {}
    for _, row in clean.iterrows():
        ann = str(row.get("annotation", "") or row.get("product", "") or "")
        ann_col[row["protein_id"]] = ann

    rows = []
    for _, row in clean.iterrows():
        pid, seq = row["protein_id"], row["sequence"]
        pos, neg = charge_counts(seq)
        tchg = _terminal_charge(seq, term_n)
        tmd_count = None
        topo_tool = ""
        sar_like = "unknown"
        topo_str = "unknown"
        sigp = "unknown"
        if not topo_idx.empty and pid in topo_idx.index:
            t = topo_idx.loc[pid]
            tmd_count = t.get("tmd_count")
            topo_tool = str(t.get("tool", ""))
            sar_like = str(t.get("sar_like", "unknown"))
            topo_str = str(t.get("topology", "unknown"))
            sigp = str(t.get("signal_peptide", "unknown"))
        try:
            tmd_count = int(tmd_count) if pd.notna(tmd_count) else None
        except (TypeError, ValueError):
            tmd_count = None
        ann = ann_col.get(pid, "")
        rows.append({
            "protein_id": pid,
            "dataset_category": row["dataset_category"],
            "length": len(seq),
            "molecular_weight": molecular_weight(seq),
            "hydrophobic_fraction": round(hydrophobic_fraction(seq), 4),
            "positive_residues": pos,
            "negative_residues": neg,
            "net_charge": net_charge(seq),
            "nterm_pos_charge": tchg["nterm_pos"],
            "nterm_neg_charge": tchg["nterm_neg"],
            "cterm_pos_charge": tchg["cterm_pos"],
            "cterm_neg_charge": tchg["cterm_neg"],
            "tmd_count": tmd_count if tmd_count is not None else -1,
            "topology": topo_str,
            "signal_peptide": sigp,
            "sar_like": sar_like,
            "topology_tool": topo_tool,
            "low_complexity_fraction": _low_complexity_fraction(seq),
            "annotation": ann,
            "has_enzymatic_or_structural_domain": _has_enzymatic_or_structural(ann),
        })
    feats = pd.DataFrame(rows)
    out = cfg.out("tables", "feature_table.tsv")
    feats.to_csv(out, sep="\t", index=False)
    return feats
