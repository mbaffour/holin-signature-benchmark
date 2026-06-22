"""Stage 7 — genomic context / lysis-cassette scoring.

For each gene in the supplied genome context, score how holin-plausible its
neighborhood is: proximity to endolysins/spanins, same-strand co-orientation,
SAR-endolysin + pinholin-topology pairing, isolation from lysis genes, and the
presence of a competing holin candidate nearby. All weights come from config.

Context is SUPPORTIVE evidence, never proof — see docs/interpretation_guide.md.
"""
from __future__ import annotations

import pandas as pd

from . import io
from .config import Config
from .utils import log


def _classify(product: str, terms: dict[str, list[str]]) -> set[str]:
    """Return the set of lysis categories matching a product annotation."""
    p = (product or "").lower()
    cats = set()
    for cat, words in terms.items():
        if any(w in p for w in words):
            cats.add(cat)
    return cats


def _gene_distance(genes: pd.DataFrame, i: int, predicate) -> tuple[int | None, int | None, bool]:
    """Nearest gene (by index) on the same contig matching predicate.

    Returns (distance_in_genes, distance_in_bp, same_strand).
    """
    row = genes.iloc[i]
    best_gene = best_bp = None
    same_strand = False
    for j in range(len(genes)):
        if j == i:
            continue
        other = genes.iloc[j]
        if not predicate(other):
            continue
        dgene = abs(j - i)
        dbp = int(min(abs(other["start"] - row["end"]), abs(row["start"] - other["end"])))
        if best_gene is None or dgene < best_gene:
            best_gene, best_bp = dgene, dbp
            same_strand = (other.get("strand") == row.get("strand"))
    return best_gene, best_bp, same_strand


def run(cfg: Config, features: pd.DataFrame | None = None) -> pd.DataFrame:
    ctx = io.load_context(cfg)
    terms = io.load_lysis_terms(cfg)
    ccfg = cfg.section("context")
    win_genes = int(ccfg.get("window_genes", 5))
    win_bp = int(ccfg.get("window_bp", 5000))
    near_genes = int(ccfg.get("endolysin_near_genes", 3))
    spanin_near_genes = int(ccfg.get("spanin_near_genes", win_genes))
    w = ccfg.get("weights", {})

    if ctx.empty:
        log.warning("No genome context supplied; Stage 7 produces an empty table.")
        empty = pd.DataFrame(columns=["protein_id", "context_score"])
        empty.to_csv(cfg.out("tables", "context_table.tsv"), sep="\t", index=False)
        return empty

    # SAR-like / pinholin topology lookup from features (optional).
    sar_topo = {}
    if features is not None and not features.empty:
        for _, r in features.iterrows():
            sar_topo[r["protein_id"]] = (str(r.get("sar_like", "")).lower() == "yes")

    rows = []
    for contig_id, genes in ctx.groupby("contig_id"):
        genes = genes.sort_values("start").reset_index(drop=True)
        # Pre-classify every gene's product.
        cat_sets = [_classify(p, terms) for p in genes["product"]]
        for i in range(len(genes)):
            g = genes.iloc[i]
            pid = g["protein_id"]

            def within_window(j):
                if abs(j - i) > win_genes:
                    return False
                other = genes.iloc[j]
                dbp = min(abs(other["start"] - g["end"]), abs(g["start"] - other["end"]))
                return dbp <= win_bp

            neigh = [j for j in range(len(genes)) if j != i and within_window(j)]
            neigh_cats = set().union(*[cat_sets[j] for j in neigh]) if neigh else set()

            # distances. `genes` was reset_index'd, so each row's positional
            # index == its label (o.name); use that to look up its category set
            # rather than re-matching on gene_id (which is non-unique in messy
            # GFF/TSV exports and would mis-score neighbours).
            d_endo_g, d_endo_bp, endo_same = _gene_distance(
                genes, i, lambda o: bool({"endolysin", "sar_endolysin"} & cat_sets[int(o.name)]))
            d_span_g, _, _ = _gene_distance(
                genes, i, lambda o: bool({"spanin"} & cat_sets[int(o.name)]))

            near_endolysin = d_endo_g is not None and d_endo_g <= near_genes
            near_spanin = d_span_g is not None and d_span_g <= spanin_near_genes
            near_sar = "sar_endolysin" in neigh_cats
            isolated = not ({"endolysin", "sar_endolysin", "spanin", "antiholin", "holin"} & neigh_cats)

            # competing holin candidate nearby: another gene annotated holin-like
            better_holin_nearby = any("holin" in cat_sets[j] and genes.iloc[j]["protein_id"] != pid
                                      for j in neigh)

            this_cats = cat_sets[i]
            alt_annotation = bool(this_cats) and not ({"holin"} & this_cats) and \
                ({"endolysin", "spanin", "sar_endolysin"} & this_cats)

            # -- additive score ----------------------------------------------
            score = 0.0
            reasons = []
            if near_endolysin:
                score += w.get("near_endolysin_within_3genes", 2.0)
                reasons.append(f"endolysin within {d_endo_g} gene(s)")
                if endo_same:
                    score += w.get("same_strand_as_endolysin", 1.0)
                    reasons.append("same strand as endolysin")
            if near_spanin:
                score += w.get("near_spanin_or_rz", 1.0)
                reasons.append("near spanin/Rz/Rz1")
            if near_sar and sar_topo.get(pid, False):
                score += w.get("near_sar_with_pinholin_topology", 1.0)
                reasons.append("near SAR endolysin + pinholin-like topology")
            if not better_holin_nearby:
                score += w.get("no_better_holin_candidate", 1.0)
                reasons.append("no competing holin candidate nearby")
            if alt_annotation:
                score += w.get("has_alternative_annotation", -2.0)
                reasons.append("has alternative (non-holin) lysis annotation")
            if isolated:
                score += w.get("not_near_any_lysis_gene", -1.0)
                reasons.append("isolated from lysis genes")

            rows.append({
                "protein_id": pid,
                "contig_id": contig_id,
                "product": g["product"],
                "strand": g.get("strand", ""),
                "near_endolysin": near_endolysin,
                "nearest_endolysin_genes": d_endo_g if d_endo_g is not None else -1,
                "nearest_endolysin_bp": d_endo_bp if d_endo_bp is not None else -1,
                "near_spanin": near_spanin,
                "nearest_spanin_genes": d_span_g if d_span_g is not None else -1,
                "near_sar_endolysin": near_sar,
                "same_strand_as_endolysin": bool(endo_same),
                "in_plausible_lysis_cassette": near_endolysin and not alt_annotation,
                "better_holin_candidate_nearby": better_holin_nearby,
                "isolated_from_lysis_genes": isolated,
                "context_score": round(score, 3),
                "context_explanation": "; ".join(reasons) if reasons else "no lysis context",
            })

    out = pd.DataFrame(rows)
    out.to_csv(cfg.out("tables", "context_table.tsv"), sep="\t", index=False)
    log.info("Stage 7: scored context for %d genes across %d contig(s)",
             len(out), ctx["contig_id"].nunique())
    return out
