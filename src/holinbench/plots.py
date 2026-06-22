"""Manuscript figures (matplotlib, headless).

Generates the dataset/feature distributions, benchmark comparison, ROC/PR curves,
confusion matrix, and final candidate-score distribution. Each function fails
soft: if an input table is empty it logs and skips, so `run-all` never crashes on
sparse data.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import Config
from .utils import log

CATEGORY_ORDER = ["gold", "weak", "hard_negative", "unknown"]


def _palette(cfg: Config) -> dict:
    return cfg.dotted("plotting.palette_by_category", {
        "gold": "#1b7837", "weak": "#7fbf7b",
        "hard_negative": "#b2182b", "unknown": "#999999"})


def _save(fig, cfg: Config, name: str) -> Path:
    fig_dir = cfg.out_dir("figures")
    dpi = int(cfg.dotted("plotting.dpi", 150))
    out = fig_dir / f"{name}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out


def _by_category_hist(feats, cfg, col, title, xlabel, name, bins=20):
    pal = _palette(cfg)
    fig, ax = plt.subplots(figsize=(6, 4))
    plotted = False
    for cat in CATEGORY_ORDER:
        vals = pd.to_numeric(feats[feats["dataset_category"] == cat][col],
                             errors="coerce").dropna()
        if len(vals):
            ax.hist(vals, bins=bins, alpha=0.55, label=cat, color=pal.get(cat))
            plotted = True
    if not plotted:
        plt.close(fig)
        return None
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("count")
    ax.legend(fontsize=8)
    return _save(fig, cfg, name)


def feature_distributions(cfg: Config, feats: pd.DataFrame) -> list[Path]:
    if feats is None or feats.empty:
        return []
    out = []
    for col, title, xlabel, name in [
        ("length", "Protein length by dataset", "length (aa)", "fig_length_distribution"),
        ("tmd_count", "TMD count by dataset", "predicted TMD count", "fig_tmd_distribution"),
        ("hydrophobic_fraction", "Hydrophobic fraction by dataset",
         "hydrophobic fraction", "fig_hydrophobic_distribution"),
    ]:
        bins = 10 if col == "tmd_count" else 20
        p = _by_category_hist(feats, cfg, col, title, xlabel, name, bins=bins)
        if p:
            out.append(p)
    return out


def cluster_plot(cfg: Config, clusters: pd.DataFrame) -> Path | None:
    if clusters is None or clusters.empty:
        return None
    thr_cols = [c for c in clusters.columns if c.startswith("cluster_")]
    if not thr_cols:
        return None
    fig, ax = plt.subplots(figsize=(6, 4))
    xs, ys = [], []
    for c in sorted(thr_cols, key=lambda x: int(x.split("_")[1])):
        thr = int(c.split("_")[1])
        xs.append(thr)
        ys.append(clusters[c].nunique())
    ax.plot(xs, ys, "o-", color="#1b7837")
    ax.set_xlabel("clustering identity threshold (%)")
    ax.set_ylabel("number of clusters")
    ax.set_title("Gold holin sequence-space fragmentation")
    for x, y in zip(xs, ys):
        ax.annotate(str(y), (x, y), textcoords="offset points", xytext=(0, 6), fontsize=8)
    return _save(fig, cfg, "fig_cluster_fragmentation")


def benchmark_comparison(cfg: Config, naive: pd.DataFrame) -> Path | None:
    if naive is None or naive.empty:
        return None
    fig, ax = plt.subplots(figsize=(7, 4))
    df = naive.copy()
    x = np.arange(len(df))
    roc = pd.to_numeric(df["roc_auc"], errors="coerce")
    pr = pd.to_numeric(df["pr_auc"], errors="coerce")
    ax.bar(x - 0.2, roc, width=0.4, label="ROC-AUC", color="#2166ac")
    ax.bar(x + 0.2, pr, width=0.4, label="PR-AUC", color="#b2182b")
    ax.set_xticks(x)
    ax.set_xticklabels(df["model"], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("AUC")
    ax.set_ylim(0, 1.05)
    ax.axhline(0.5, ls="--", c="gray", lw=0.8)
    ax.set_title("Model comparison (NAIVE regime — HMM models inflated)")
    ax.legend(fontsize=8)
    return _save(fig, cfg, "fig_benchmark_comparison")


def roc_pr_curves(cfg: Config, scores: pd.DataFrame, models: list[str]) -> list[Path]:
    if scores is None or scores.empty or scores["label"].nunique() < 2:
        return []
    from sklearn.metrics import precision_recall_curve, roc_curve
    y = scores["label"].to_numpy()
    out = []
    figr, axr = plt.subplots(figsize=(5, 5))
    figp, axp = plt.subplots(figsize=(5, 5))
    for m in models:
        if m not in scores.columns:
            continue
        s = pd.to_numeric(scores[m], errors="coerce").fillna(0).to_numpy()
        fpr, tpr, _ = roc_curve(y, s)
        axr.plot(fpr, tpr, label=m, lw=1.2)
        prec, rec, _ = precision_recall_curve(y, s)
        axp.plot(rec, prec, label=m, lw=1.2)
    axr.plot([0, 1], [0, 1], ls="--", c="gray", lw=0.8)
    axr.set_xlabel("false positive rate")
    axr.set_ylabel("true positive rate")
    axr.set_title("ROC (naive)")
    axr.legend(fontsize=7)
    out.append(_save(figr, cfg, "fig_roc_curves"))
    axp.set_xlabel("recall")
    axp.set_ylabel("precision")
    axp.set_title("Precision-Recall (naive)")
    axp.legend(fontsize=7)
    out.append(_save(figp, cfg, "fig_pr_curves"))
    return out


def confusion_matrix_plot(cfg: Config, fp_fn: pd.DataFrame, scores: pd.DataFrame) -> Path | None:
    if scores is None or scores.empty or scores["label"].nunique() < 2:
        return None
    from sklearn.metrics import confusion_matrix
    from .benchmark import _best_f1_threshold
    y = scores["label"].tolist()
    s = scores["hmm_arch_context"].tolist()
    thr = _best_f1_threshold(y, s)
    pred = [1 if v >= thr else 0 for v in s]
    cm = confusion_matrix(y, pred, labels=[0, 1])
    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(cm, cmap="Blues")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                    color="black", fontsize=12)
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(["neg", "pos"]); ax.set_yticklabels(["neg", "pos"])
    ax.set_xlabel("predicted"); ax.set_ylabel("true")
    ax.set_title(f"Confusion (full model, thr={thr:.2f})")
    fig.colorbar(im, fraction=0.046)
    return _save(fig, cfg, "fig_confusion_matrix")


def candidate_score_distribution(cfg: Config, ranking: pd.DataFrame) -> Path | None:
    if ranking is None or ranking.empty:
        return None
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(pd.to_numeric(ranking["final_holin_score"], errors="coerce").dropna(),
            bins=20, color="#542788", alpha=0.8)
    for cut, lab in [(0.70, "high"), (0.50, "med"), (0.30, "weak")]:
        ax.axvline(cut, ls="--", c="gray", lw=0.8)
        ax.text(cut, ax.get_ylim()[1] * 0.9, lab, fontsize=7, rotation=90)
    ax.set_xlabel("final holin candidate score")
    ax.set_ylabel("count")
    ax.set_title("Candidate score distribution")
    return _save(fig, cfg, "fig_candidate_scores")


def run(cfg: Config, feats=None, clusters=None, naive=None, scores=None,
        fp_fn=None, ranking=None) -> list[Path]:
    made = []
    made += feature_distributions(cfg, feats)
    p = cluster_plot(cfg, clusters); made += [p] if p else []
    p = benchmark_comparison(cfg, naive); made += [p] if p else []
    models = cfg.dotted("benchmark.models", [])
    made += roc_pr_curves(cfg, scores, models)
    p = confusion_matrix_plot(cfg, fp_fn, scores); made += [p] if p else []
    p = candidate_score_distribution(cfg, ranking); made += [p] if p else []
    made = [m for m in made if m]
    log.info("Plots: wrote %d figure(s) to %s", len(made), cfg.out_dir("figures"))
    return made
