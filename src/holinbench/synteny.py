"""Stage 11 — lysis-cassette synteny visualization.

Draws gene-arrow maps of the genomic neighborhood around selected candidates,
coloured by functional class (candidate / endolysin / spanin / antiholin /
structural / hypothetical / other). Exports PNG and/or SVG. Useful for visually
inspecting whether a high-scoring candidate really sits in a plausible lysis
cassette — supportive evidence, not proof.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
import pandas as pd

from . import io
from .config import Config
from .utils import log

CLASS_COLORS = {
    "candidate": "#542788",
    "endolysin": "#1b7837",
    "sar_endolysin": "#5aae61",
    "spanin": "#2166ac",
    "antiholin": "#f1a340",
    "holin": "#9970ab",
    "structural": "#888888",
    "hypothetical": "#cccccc",
    "other": "#e0e0e0",
}


def _classify(product: str, terms: dict[str, list[str]], is_candidate: bool) -> str:
    if is_candidate:
        return "candidate"
    p = (product or "").lower()
    for cat in ("endolysin", "sar_endolysin", "spanin", "antiholin", "holin"):
        for w in terms.get(cat, []):
            if w in p:
                return cat
    if any(w in p for w in ["capsid", "tail", "portal", "baseplate", "terminase", "structural"]):
        return "structural"
    if "hypothetical" in p:
        return "hypothetical"
    return "other"


def draw_cassette(genes: pd.DataFrame, candidate_id: str, terms: dict,
                  flank: int, formats: list[str], out_dir: Path) -> list[Path]:
    genes = genes.sort_values("start").reset_index(drop=True)
    idx = genes.index[genes["protein_id"] == candidate_id]
    if len(idx) == 0:
        return []
    i = int(idx[0])
    lo, hi = max(0, i - flank), min(len(genes), i + flank + 1)
    window = genes.iloc[lo:hi]

    fig, ax = plt.subplots(figsize=(max(6, len(window) * 1.1), 2.4))
    xmin = window["start"].min()
    xmax = window["end"].max()
    span = (xmax - xmin) or 1
    for _, g in window.iterrows():
        is_cand = (g["protein_id"] == candidate_id)
        cls = _classify(g["product"], terms, is_cand)
        color = CLASS_COLORS.get(cls, "#e0e0e0")
        strand = 1 if str(g.get("strand", "+")) != "-" else -1
        x0 = g["start"] if strand == 1 else g["end"]
        dx = (g["end"] - g["start"]) * strand
        ax.add_patch(FancyArrow(
            x0, 0, dx, 0, width=0.3, length_includes_head=True,
            head_width=0.5, head_length=min(abs(dx) * 0.4, span * 0.03),
            color=color, ec="black", lw=0.6))
        ax.text((g["start"] + g["end"]) / 2, 0.55,
                g["product"][:18], ha="center", va="bottom", fontsize=6, rotation=20)
        if is_cand:
            ax.text((g["start"] + g["end"]) / 2, -0.6, "★ candidate",
                    ha="center", va="top", fontsize=7, color=CLASS_COLORS["candidate"])
    ax.set_xlim(xmin - span * 0.05, xmax + span * 0.05)
    ax.set_ylim(-1.2, 1.5)
    ax.set_yticks([])
    ax.set_xlabel("genome coordinate (bp)")
    ax.set_title(f"Lysis-cassette context: {candidate_id}")
    fig.tight_layout()

    paths = []
    for fmt in formats:
        out = out_dir / f"synteny_{candidate_id}.{fmt}"
        fig.savefig(out, dpi=150)
        paths.append(out)
    plt.close(fig)
    return paths


def run(cfg: Config, ranking: pd.DataFrame | None = None) -> list[Path]:
    ctx = io.load_context(cfg)
    terms = io.load_lysis_terms(cfg)
    scfg = cfg.section("synteny")
    flank = int(scfg.get("flank_genes", 5))
    formats = scfg.get("formats", ["png"])
    max_maps = int(scfg.get("max_maps", 10))
    fig_dir = cfg.out_dir("figures")

    if ctx.empty:
        log.warning("No genome context; skipping synteny maps.")
        return []

    # Candidates to draw: top-ranked that appear in the context, else all context
    # proteins present in candidate set.
    candidates = []
    if ranking is not None and not ranking.empty:
        ctx_ids = set(ctx["protein_id"])
        candidates = [pid for pid in ranking["protein_id"] if pid in ctx_ids][:max_maps]
    if not candidates:
        candidates = list(dict.fromkeys(ctx["protein_id"]))[:max_maps]

    made = []
    for cand in candidates:
        contig = ctx[ctx["protein_id"] == cand]["contig_id"].iloc[0]
        genes = ctx[ctx["contig_id"] == contig]
        made += draw_cassette(genes, cand, terms, flank, formats, fig_dir)
    log.info("Stage 11: drew %d synteny file(s) for %d candidate(s).",
             len(made), len(candidates))
    return made
