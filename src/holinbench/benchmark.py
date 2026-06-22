"""Stage 9 — benchmarking holin-prediction approaches.

Compares seven approaches (A universal HMM, B topology HMM, C family HMM,
D architecture-only, E context-only, F HMM+architecture, G HMM+architecture+
context) at distinguishing gold positives from HARD negatives.

Two evaluation regimes are reported side by side:
  * NAIVE — scores from models built on ALL gold, evaluated on gold + negatives.
    This is CIRCULAR for the HMM models (they were built from the same gold) and
    is explicitly labeled as optimistic / potentially inflated.
  * LEAVE-ONE-FAMILY-OUT (LOFO) — the universal HMM is rebuilt excluding each
    family in turn and asked to recover the held-out family. This is the honest
    test of whether a single universal model GENERALIZES across diverse holins,
    which is the project's central question.

Small sample sizes are reported, never hidden. ROC/PR-AUC come with bootstrap
confidence intervals where both classes are present.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import alignment as alignment_mod, clustering, hmmer, scoring
from .config import Config
from .scoring import architecture_score, context_norm, hmm_score
from .utils import log

POS_LABEL, NEG_LABEL = "gold", "hard_negative"


# --------------------------------------------------------- metrics -------------
def _safe_auc(y_true, y_score, kind: str) -> float:
    from sklearn.metrics import average_precision_score, roc_auc_score
    if len(set(y_true)) < 2:
        return float("nan")
    try:
        return float(roc_auc_score(y_true, y_score) if kind == "roc"
                     else average_precision_score(y_true, y_score))
    except Exception:
        return float("nan")


def _operating_point(y_true, y_score, threshold) -> dict:
    from sklearn.metrics import confusion_matrix
    y_pred = [1 if s >= threshold else 0 for s in y_score]
    if len(set(y_true)) < 2:
        return {"threshold": threshold, "precision": float("nan"),
                "recall": float("nan"), "specificity": float("nan"),
                "f1": float("nan"), "tp": 0, "fp": 0, "tn": 0, "fn": 0}
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {"threshold": round(float(threshold), 3), "precision": round(prec, 3),
            "recall": round(rec, 3), "specificity": round(spec, 3),
            "f1": round(f1, 3), "tp": int(tp), "fp": int(fp),
            "tn": int(tn), "fn": int(fn)}


def _best_f1_threshold(y_true, y_score) -> float:
    # Iterate thresholds ASCENDING and keep the best F1 with `>=`, so among ties
    # the HIGHER threshold wins. This favors specificity rather than silently
    # biasing the operating point toward calling everything positive.
    cands = sorted(set(y_score))
    if not cands:
        return 0.5
    best_t, best_f1 = cands[0], -1.0
    for t in cands:
        op = _operating_point(y_true, y_score, t)
        if op["f1"] == op["f1"] and op["f1"] >= best_f1:
            best_f1, best_t = op["f1"], t
    return best_t


def _bootstrap_auc_ci(y_true, y_score, kind, n=1000, ci=0.95, seed=1729):
    if len(set(y_true)) < 2:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    idx = np.arange(len(y_true))
    vals = []
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        if len(set(y_true[s])) < 2:
            continue
        vals.append(_safe_auc(y_true[s], y_score[s], kind))
    if not vals:
        return (float("nan"), float("nan"))
    lo = float(np.nanpercentile(vals, 100 * (1 - ci) / 2))
    hi = float(np.nanpercentile(vals, 100 * (1 + ci) / 2))
    return (round(lo, 3), round(hi, 3))


def _topk_recovery(y_true, y_score, k=None) -> float:
    n_pos = sum(y_true)
    if n_pos == 0:
        return float("nan")
    k = k or n_pos
    # Deterministic, PESSIMISTIC tie-break: sort by descending score, and among
    # equal scores place negatives (label 0) ahead of positives (label 1), so a
    # tie between a positive and a negative does not flatter recovery.
    order = sorted(range(len(y_score)), key=lambda i: (-y_score[i], y_true[i]))
    top = order[:k]
    return round(sum(y_true[i] for i in top) / n_pos, 3)


# --------------------------------------------------- per-model score vectors ---
def _hmm_type_scores(hits: pd.DataFrame) -> dict[str, dict[str, float]]:
    """protein_id -> {universal,topology,family: max bitscore}."""
    out: dict[str, dict[str, float]] = {}
    if hits is None or hits.empty:
        return out
    for pid, sub in hits.groupby("protein_id"):
        d = {"universal": 0.0, "topology": 0.0, "family": 0.0}
        for mt, g in sub.groupby("model_type"):
            d[mt] = float(g["bitscore"].max())
        out[pid] = d
    return out


def _build_score_table(cfg: Config, feats: pd.DataFrame, context: pd.DataFrame | None,
                       hits: pd.DataFrame, exclude_ids: set | None = None) -> pd.DataFrame:
    """Per-protein score for each of the 7 models, restricted to gold + negatives.

    `exclude_ids` drops leakage cases (e.g. a hard negative whose sequence is
    identical to a gold positive) so they cannot taint the metrics.
    """
    sub = feats[feats["dataset_category"].isin([POS_LABEL, NEG_LABEL])].copy()
    if exclude_ids:
        sub = sub[~sub["protein_id"].isin(exclude_ids)]
    type_scores = _hmm_type_scores(hits)
    ctx_ids = set()
    ctx_lookup = {}
    if context is not None and not context.empty:
        ctx_lookup = dict(zip(context["protein_id"], context["context_score"]))
        ctx_ids = set(context["protein_id"])
    w = cfg.dotted("scoring.weights", {"hmm": 0.45, "architecture": 0.30, "context": 0.25})

    rows = []
    for _, fr in sub.iterrows():
        pid = fr["protein_id"]
        label = 1 if fr["dataset_category"] == POS_LABEL else 0
        ts = type_scores.get(pid, {"universal": 0.0, "topology": 0.0, "family": 0.0})
        phits = hits[hits["protein_id"] == pid] if (hits is not None and not hits.empty) else None
        hmm_val, _, _ = hmm_score(phits, cfg)
        arch, _ = architecture_score(fr.to_dict(), cfg)
        has_ctx = pid in ctx_ids
        raw_ctx = float(ctx_lookup.get(pid, 0.0))
        ctx_val = context_norm(raw_ctx)
        wh, wa, wc = w.get("hmm", 0.45), w.get("architecture", 0.30), w.get("context", 0.25)
        hmm_arch = (wh * hmm_val + wa * arch) / (wh + wa)
        # Mirror scoring.py: when a protein has no context row, drop the context
        # term and renormalize, rather than crediting a flat 0.5.
        if has_ctx:
            full = wh * hmm_val + wa * arch + wc * ctx_val
        else:
            full = (wh * hmm_val + wa * arch) / (wh + wa)
        rows.append({
            "protein_id": pid, "label": label, "has_context": has_ctx,
            "universal": ts["universal"], "topology": ts["topology"], "family": ts["family"],
            "architecture": arch, "context": ctx_val,
            "hmm_arch": round(hmm_arch, 4), "hmm_arch_context": round(full, 4),
        })
    return pd.DataFrame(rows)


def _evaluate_models(score_df: pd.DataFrame, models: list[str], cfg: Config) -> pd.DataFrame:
    bcfg = cfg.section("benchmark")
    boot = bcfg.get("bootstrap", {})
    do_boot = bool(boot.get("enabled", True))
    n_res = int(boot.get("n_resamples", 1000))
    ci = float(boot.get("ci", 0.95))
    y_true = score_df["label"].tolist()
    n_pos, n_neg = sum(y_true), len(y_true) - sum(y_true)

    rows = []
    for model in models:
        if model not in score_df.columns:
            continue
        y_score = score_df[model].tolist()
        roc = _safe_auc(y_true, y_score, "roc")
        pr = _safe_auc(y_true, y_score, "pr")
        thr = _best_f1_threshold(y_true, y_score)
        op = _operating_point(y_true, y_score, thr)
        roc_ci = _bootstrap_auc_ci(y_true, y_score, "roc", n_res, ci) if do_boot else (np.nan, np.nan)
        rows.append({
            "model": model, "regime": "naive_full_data",
            "n_pos": n_pos, "n_neg": n_neg,
            "roc_auc": round(roc, 3) if roc == roc else float("nan"),
            "roc_auc_ci_low": roc_ci[0], "roc_auc_ci_high": roc_ci[1],
            "pr_auc": round(pr, 3) if pr == pr else float("nan"),
            "topk_recall": _topk_recovery(y_true, y_score),
            **op,
            "circularity_warning": "YES (HMM built on these gold)" if model in
                                   ("universal", "topology", "family", "hmm_arch",
                                    "hmm_arch_context") else "no",
            "note": "",
        })
    return pd.DataFrame(rows)


# --------------------------------------------------- LOFO (honest) -------------
def leave_one_family_out_universal(cfg: Config, clean: pd.DataFrame,
                                   feats: pd.DataFrame,
                                   exclude_ids: set | None = None) -> pd.DataFrame:
    """Rebuild the universal HMM excluding each family; test recovery of the
    held-out family vs. false hits on hard negatives. The honest generalization test.
    """
    gold = clean[clean["dataset_category"] == "gold"].copy()
    neg = clean[clean["dataset_category"] == "hard_negative"].copy()
    if exclude_ids:
        # Drop leakage negatives (identical to a gold sequence) so they don't
        # inflate the false-hit count.
        neg = neg[~neg["protein_id"].isin(exclude_ids)]
    if gold.empty:
        return pd.DataFrame()

    # Family key: family_label if present, else cluster at primary threshold.
    # Robust to empty-string OR NaN labels (don't rely on truthiness of NaN).
    fam_label = gold["family_label"].fillna("").astype(str)
    needs_cluster = fam_label.str.strip() == ""
    if needs_cluster.any():
        clusters = clustering.run(cfg, clean)
        pcol = f"cluster_{int(cfg.dotted('clustering.primary_threshold', 0.30) * 100)}"
        cmap = dict(zip(clusters["protein_id"], clusters[pcol])) if pcol in clusters else {}
        fam_key = [lab if lab.strip() else f"cluster{cmap.get(pid, 'NA')}"
                   for lab, pid in zip(fam_label, gold["protein_id"])]
    else:
        fam_key = list(fam_label)
    gold = gold.assign(_fam=fam_key)

    evalue = float(cfg.dotted("hmmer.evalue_cutoff", 0.01))
    neg_records_df = neg[["protein_id", "sequence", "dataset_category"]].copy()

    rows = []
    families = sorted(gold["_fam"].unique())
    for fam in families:
        test = gold[gold["_fam"] == fam]
        train = gold[gold["_fam"] != fam]
        if len(train) < 2:
            rows.append({"held_out_family": fam, "n_test": len(test),
                         "n_train": len(train), "recovered": np.nan,
                         "recall": np.nan, "neg_false_hits": np.nan,
                         "note": "insufficient training sequences"})
            continue
        train_records = list(zip(train["protein_id"], train["sequence"]))
        model = hmmer.build_single(cfg, f"universal_minus_{fam}", "universal", train_records)
        if model is None:
            rows.append({"held_out_family": fam, "n_test": len(test),
                         "n_train": len(train), "recovered": np.nan, "recall": np.nan,
                         "neg_false_hits": np.nan, "note": "model build failed"})
            continue
        # scan held-out gold + all negatives
        fold_clean = pd.concat([
            test[["protein_id", "sequence", "dataset_category"]], neg_records_df],
            ignore_index=True)
        hits = hmmer.scan(cfg, [model], fold_clean, write_output=False)
        hit_ids = set(hits[hits["evalue"] <= evalue]["protein_id"]) if not hits.empty else set()
        recovered = sum(1 for pid in test["protein_id"] if pid in hit_ids)
        neg_fp = sum(1 for pid in neg["protein_id"] if pid in hit_ids)
        rows.append({
            "held_out_family": fam, "n_test": len(test), "n_train": len(train),
            "recovered": recovered, "recall": round(recovered / len(test), 3),
            "neg_false_hits": neg_fp, "n_neg": len(neg),
            "note": "",
        })
    out = pd.DataFrame(rows)
    return out


# --------------------------------------------------- driver --------------------
def run(cfg: Config, feats: pd.DataFrame, context: pd.DataFrame | None,
        hits: pd.DataFrame, clean: pd.DataFrame) -> dict:
    models = cfg.dotted("benchmark.models",
                        ["universal", "topology", "family", "architecture",
                         "context", "hmm_arch", "hmm_arch_context"])

    # Identify leakage: hard negatives whose exact sequence equals a gold
    # positive's. These are unwinnable false positives that would distort both
    # regimes (a sequence identical to a training positive scored as a negative).
    gold_seqs = set(clean[clean["dataset_category"] == "gold"]["sequence"])
    neg_df = clean[clean["dataset_category"] == "hard_negative"]
    leak_ids = set(neg_df[neg_df["sequence"].isin(gold_seqs)]["protein_id"])
    if leak_ids:
        log.warning("Stage 9: excluding %d hard negative(s) identical to a gold "
                    "sequence (label leakage): %s", len(leak_ids), sorted(leak_ids))
        pd.DataFrame({"protein_id": sorted(leak_ids),
                      "reason": "sequence identical to a gold positive"}).to_csv(
            cfg.out("tables", "benchmark_excluded_leakage.tsv"), sep="\t", index=False)

    score_df = _build_score_table(cfg, feats, context, hits, exclude_ids=leak_ids)
    n_pos = int(score_df["label"].sum())
    n_neg = int(len(score_df) - n_pos)
    warn_n = int(cfg.dotted("benchmark.validation.min_class_size_warn", 5))
    if n_pos < warn_n or n_neg < warn_n:
        log.warning("Stage 9: SMALL SAMPLE — %d positives, %d negatives "
                    "(< %d). Metrics are unstable; interpret with caution.",
                    n_pos, n_neg, warn_n)

    # Context-only model is only meaningful if positives actually have context.
    pos_ctx_cov = 0.0
    if n_pos and "has_context" in score_df.columns:
        pos_ctx_cov = float(score_df[score_df["label"] == 1]["has_context"].mean())

    naive = _evaluate_models(score_df, models, cfg)
    # Guard H2: if positives have ~no genomic context, the context-only AUC is an
    # artifact (it ranks the few annotated negatives). Mark it not-evaluable
    # instead of reporting a misleading sub-random number.
    if pos_ctx_cov < 0.2 and "context" in set(naive["model"]):
        mask = naive["model"] == "context"
        for col in ["roc_auc", "roc_auc_ci_low", "roc_auc_ci_high", "pr_auc",
                    "topk_recall", "precision", "recall", "specificity", "f1"]:
            if col in naive.columns:
                naive.loc[mask, col] = float("nan")
        naive.loc[mask, "note"] = (f"not evaluable: context present for only "
                                   f"{pos_ctx_cov:.0%} of positives")
        log.info("Stage 9: context-only model marked not-evaluable "
                 "(context covers %.0f%% of positives).", pos_ctx_cov * 100)
    naive.to_csv(cfg.out("tables", "benchmark_metrics.tsv"), sep="\t", index=False)
    score_df.to_csv(cfg.out("tables", "benchmark_scores.tsv"), sep="\t", index=False)

    lofo = pd.DataFrame()
    if cfg.dotted("benchmark.validation.do_leave_one_family_out", True):
        lofo = leave_one_family_out_universal(cfg, clean, feats, exclude_ids=leak_ids)
        lofo.to_csv(cfg.out("tables", "benchmark_lofo_universal.tsv"), sep="\t", index=False)
        if not lofo.empty and lofo["recall"].notna().any():
            mean_recall = lofo["recall"].dropna().mean()
            log.info("Stage 9 LOFO: universal HMM mean held-out recall = %.2f "
                     "across %d families (honest generalization).",
                     mean_recall, lofo["recall"].notna().sum())

    # False positives / negatives at the full model's best-F1 operating point.
    fp_fn = _false_pos_neg(score_df, feats, cfg)
    fp_fn.to_csv(cfg.out("tables", "false_pos_neg.tsv"), sep="\t", index=False)

    return {"naive": naive, "lofo": lofo, "scores": score_df, "fp_fn": fp_fn}


def _false_pos_neg(score_df: pd.DataFrame, feats: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    y_true = score_df["label"].tolist()
    y_score = score_df["hmm_arch_context"].tolist()
    thr = _best_f1_threshold(y_true, y_score)
    fmap = feats.set_index("protein_id")
    rows = []
    for _, r in score_df.iterrows():
        pred = 1 if r["hmm_arch_context"] >= thr else 0
        if pred != r["label"]:
            kind = "false_positive" if (pred == 1 and r["label"] == 0) else "false_negative"
            ann = fmap.loc[r["protein_id"]]["annotation"] if r["protein_id"] in fmap.index else ""
            rows.append({"protein_id": r["protein_id"], "error_type": kind,
                         "score": r["hmm_arch_context"], "threshold": round(thr, 3),
                         "annotation": ann})
    return pd.DataFrame(rows, columns=["protein_id", "error_type", "score",
                                       "threshold", "annotation"])
