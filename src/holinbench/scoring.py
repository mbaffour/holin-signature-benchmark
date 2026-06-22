"""Stage 8 — composite holin candidate score.

Combines three normalized sub-scores — HMM evidence, architecture, and genomic
context — into a single ``final_holin_score`` in [0, 1], with a confidence
category and a human-readable explanation per candidate.

This is a RANKING/PRIORITIZATION system, not a classifier and not proof of
function. All weights are configurable; nothing biological is hardcoded.
"""
from __future__ import annotations

import math

import pandas as pd

from .config import Config
from .utils import log


# ----------------------------------------------------------- sub-scores --------
def architecture_score(row: dict, cfg: Config) -> tuple[float, list[str]]:
    a = cfg.dotted("scoring.architecture", {})
    reasons: list[str] = []
    length = int(row.get("length", 0) or 0)
    tmd = int(row.get("tmd_count", -1) if row.get("tmd_count") is not None else -1)
    hydro = float(row.get("hydrophobic_fraction", 0.0) or 0.0)

    pref_lo, pref_hi = a.get("preferred_length", [50, 150])
    acc_lo, acc_hi = a.get("acceptable_length", [30, 200])
    if pref_lo <= length <= pref_hi:
        length_comp = a.get("length_in_preferred", 1.0)
        reasons.append(f"length {length} in preferred {pref_lo}-{pref_hi}")
    elif acc_lo <= length <= acc_hi:
        length_comp = a.get("length_in_acceptable", 0.5)
        reasons.append(f"length {length} acceptable")
    else:
        length_comp = a.get("length_out", 0.0)
        reasons.append(f"length {length} outside holin range")

    tlo, thi = a.get("preferred_tmd", [1, 4])
    if tmd >= 0 and tlo <= tmd <= thi:
        tmd_comp = a.get("tmd_in_range", 1.0)
        reasons.append(f"{tmd} TMD(s) in range {tlo}-{thi}")
    else:
        tmd_comp = a.get("tmd_out_of_range", 0.0)
        reasons.append(f"TMD count {tmd} outside {tlo}-{thi}")

    min_hydro = a.get("min_hydrophobic_fraction", 0.40)
    if hydro >= min_hydro:
        hydro_comp = a.get("hydrophobic_ok", 1.0)
        reasons.append(f"hydrophobic fraction {hydro:.2f} >= {min_hydro}")
    else:
        hydro_comp = 0.0
        reasons.append(f"hydrophobic fraction {hydro:.2f} low")

    positives = (length_comp + tmd_comp + hydro_comp) / 3.0
    penalty = 0.0
    if row.get("has_enzymatic_or_structural_domain"):
        penalty += a.get("has_enzymatic_or_structural_domain", -1.0)
        reasons.append("PENALTY: enzymatic/structural domain annotation")
    if tmd == 0:
        penalty += a.get("no_tmd", -1.0)
        reasons.append("PENALTY: no predicted TMD")
    score = max(0.0, min(1.0, positives + penalty))
    return score, reasons


def hmm_score(hits: pd.DataFrame, cfg: Config) -> tuple[float, dict, list[str]]:
    """hits = rows for ONE protein (cols: model_type, model_id, bitscore, evalue)."""
    h = cfg.dotted("scoring.hmm", {})
    strong = float(h.get("strong_bitscore", 30.0))
    weak = float(h.get("weak_bitscore", 12.0))
    info = {"best_hmm_model": "", "best_hmm_bitscore": float("nan"),
            "best_hmm_evalue": float("nan"),
            "universal_hmm_hit": False, "topology_hmm_hit": False, "family_hmm_hit": False}
    reasons: list[str] = []
    if hits is None or hits.empty:
        reasons.append("no HMM hit")
        return 0.0, info, reasons

    hits = hits.sort_values("bitscore", ascending=False)
    best = hits.iloc[0]
    info["best_hmm_model"] = str(best["model_id"])
    info["best_hmm_bitscore"] = float(best["bitscore"])
    info["best_hmm_evalue"] = float(best["evalue"])

    value = 0.0
    for _, hit in hits.iterrows():
        mt = str(hit["model_type"])
        bits = float(hit["bitscore"])
        if mt == "universal":
            info["universal_hmm_hit"] = True
            value = max(value, h.get("universal_hit", 0.4))
        elif mt == "topology":
            info["topology_hmm_hit"] = True
            value = max(value, h.get("topology_hit", 0.5))
        elif mt == "family":
            info["family_hmm_hit"] = True
            if bits >= strong:
                value = max(value, h.get("family_hit_strong", 1.0))
            elif bits >= weak:
                value = max(value, h.get("family_hit_weak", 0.6))
            else:
                value = max(value, h.get("family_hit_weak", 0.6) * 0.5)
    reasons.append(f"best HMM {info['best_hmm_model']} (bits={info['best_hmm_bitscore']:.1f}, "
                   f"E={info['best_hmm_evalue']:.1g})")
    return value, info, reasons


def context_norm(raw_score: float) -> float:
    """Squash an additive context score to [0,1] (logistic, midpoint at 0)."""
    return 1.0 / (1.0 + math.exp(-raw_score / 2.0))


def confidence_category(score: float, cfg: Config) -> str:
    b = cfg.dotted("scoring.confidence_bins", {})
    if score >= b.get("high_confidence_candidate", 0.70):
        return "high_confidence_candidate"
    if score >= b.get("medium_confidence_candidate", 0.50):
        return "medium_confidence_candidate"
    if score >= b.get("weak_candidate", 0.30):
        return "weak_candidate"
    return "unlikely_holin"


# ----------------------------------------------------------- driver ------------
def run(cfg: Config, features: pd.DataFrame,
        context: pd.DataFrame | None = None,
        hmm_hits: pd.DataFrame | None = None,
        restrict_categories: list[str] | None = None) -> pd.DataFrame:
    weights = cfg.dotted("scoring.weights", {"hmm": 0.45, "architecture": 0.30, "context": 0.25})
    ctx_lookup = {}
    if context is not None and not context.empty:
        for _, r in context.iterrows():
            ctx_lookup[r["protein_id"]] = r

    ranking_cols = ["protein_id", "dataset_category", "best_hmm_model", "best_hmm_evalue",
                    "best_hmm_bitscore", "universal_hmm_hit", "topology_hmm_hit",
                    "family_hmm_hit", "length", "tmd_count", "hydrophobic_fraction",
                    "near_endolysin", "near_spanin", "context_score", "architecture_score",
                    "hmm_score", "final_holin_score", "confidence_category", "explanation"]

    rows = []
    feats = features
    if restrict_categories:
        feats = feats[feats["dataset_category"].isin(restrict_categories)]

    if feats.empty:
        empty = pd.DataFrame(columns=ranking_cols)
        empty.to_csv(cfg.out("tables", "candidate_ranking.tsv"), sep="\t", index=False)
        log.warning("Stage 8: no proteins to score (empty/over-filtered input).")
        return empty

    for _, frow in feats.iterrows():
        pid = frow["protein_id"]
        arch, arch_reasons = architecture_score(frow.to_dict(), cfg)

        phits = None
        if hmm_hits is not None and not hmm_hits.empty:
            phits = hmm_hits[hmm_hits["protein_id"] == pid]
        hmm_val, hmm_info, hmm_reasons = hmm_score(phits, cfg)

        cinfo = ctx_lookup.get(pid)
        has_ctx = cinfo is not None
        raw_ctx = float(cinfo["context_score"]) if has_ctx else 0.0
        ctx_val = context_norm(raw_ctx)
        ctx_reason = (cinfo["context_explanation"] if has_ctx
                      else "no genomic context supplied")

        w_hmm = weights.get("hmm", 0.45)
        w_arch = weights.get("architecture", 0.30)
        w_ctx = weights.get("context", 0.25)
        if has_ctx:
            final = w_hmm * hmm_val + w_arch * arch + w_ctx * ctx_val
            ctx_expl = f"Context[{ctx_val:.2f} from raw {raw_ctx:+.1f}]: {ctx_reason}."
        else:
            # No context row for this protein: drop the context term and
            # renormalize over the remaining evidence, so "no data" neither
            # rewards (a flat 0.5) nor penalizes the candidate relative to one
            # whose context was actually measured.
            denom = (w_hmm + w_arch) or 1.0
            final = (w_hmm * hmm_val + w_arch * arch) / denom
            ctx_expl = f"Context[n/a]: {ctx_reason} (context weight renormalized out)."
        final = round(max(0.0, min(1.0, final)), 4)
        category = confidence_category(final, cfg)

        explanation = (
            f"score={final:.2f} ({category}). "
            f"HMM[{hmm_val:.2f}]: {'; '.join(hmm_reasons)}. "
            f"Architecture[{arch:.2f}]: {'; '.join(arch_reasons)}. "
            f"{ctx_expl}"
        )

        rows.append({
            "protein_id": pid,
            "dataset_category": frow["dataset_category"],
            "best_hmm_model": hmm_info["best_hmm_model"],
            "best_hmm_evalue": hmm_info["best_hmm_evalue"],
            "best_hmm_bitscore": hmm_info["best_hmm_bitscore"],
            "universal_hmm_hit": hmm_info["universal_hmm_hit"],
            "topology_hmm_hit": hmm_info["topology_hmm_hit"],
            "family_hmm_hit": hmm_info["family_hmm_hit"],
            "length": int(frow.get("length", 0) or 0),
            "tmd_count": int(frow.get("tmd_count", -1) if frow.get("tmd_count") is not None else -1),
            "hydrophobic_fraction": float(frow.get("hydrophobic_fraction", 0.0) or 0.0),
            "near_endolysin": bool(cinfo["near_endolysin"]) if cinfo is not None else False,
            "near_spanin": bool(cinfo["near_spanin"]) if cinfo is not None else False,
            "context_score": raw_ctx,
            "architecture_score": round(arch, 4),
            "hmm_score": round(hmm_val, 4),
            "final_holin_score": final,
            "confidence_category": category,
            "explanation": explanation,
        })

    ranking = pd.DataFrame(rows).sort_values("final_holin_score", ascending=False).reset_index(drop=True)
    out = cfg.out("tables", "candidate_ranking.tsv")
    ranking.to_csv(out, sep="\t", index=False)
    log.info("Stage 8: scored %d proteins; %s",
             len(ranking),
             ", ".join(f"{c}={n}" for c, n in
                       ranking["confidence_category"].value_counts().items()))
    return ranking
