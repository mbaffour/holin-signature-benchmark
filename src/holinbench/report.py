"""Stage 12 — manuscript-style Markdown report.

Reads the result tables and assembles a cautious, analytical report with all 13
sections from the spec. Conclusion language is chosen by the DATA (e.g. a weak
universal HMM yields the "universal model performs poorly" template). The tone is
deliberately conservative: every claim is hedged and tied to a table or figure.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import Config
from .utils import log


def _read(cfg: Config, name: str) -> pd.DataFrame:
    p = cfg.out_dir("tables") / name
    if p.exists():
        try:
            return pd.read_csv(p, sep="\t")
        except Exception:
            return pd.DataFrame()
    return pd.DataFrame()


def _md_table(df: pd.DataFrame, max_rows: int = 20, cols: list[str] | None = None) -> str:
    if df is None or df.empty:
        return "_(no data)_\n"
    if cols:
        df = df[[c for c in cols if c in df.columns]]
    df = df.head(max_rows)
    header = "| " + " | ".join(str(c) for c in df.columns) + " |"
    sep = "| " + " | ".join("---" for _ in df.columns) + " |"
    lines = [header, sep]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(_fmt(v) for v in r.tolist()) + " |")
    return "\n".join(lines) + "\n"


def _fmt(v) -> str:
    if isinstance(v, float):
        if v != v:
            return "NA"
        return f"{v:.3g}"
    s = str(v).replace("|", "\\|")
    return s if len(s) <= 60 else s[:57] + "..."


def _fig(cfg: Config, name: str, caption: str) -> str:
    p = cfg.out_dir("figures") / name
    if p.exists():
        return f"![{caption}](figures/{name})\n\n*{caption}*\n"
    return ""


def _get(df: pd.DataFrame, where_col, where_val, get_col, default="NA"):
    if df.empty or where_col not in df.columns or get_col not in df.columns:
        return default
    sub = df[df[where_col] == where_val]
    if sub.empty:
        return default
    return sub.iloc[0][get_col]


def run(cfg: Config) -> Path:
    clean = _read(cfg, "cleaned_metadata.tsv")
    flags = _read(cfg, "qc_flags.tsv")
    clusters = _read(cfg, "cluster_summary.tsv")
    aln = _read(cfg, "alignment_quality.tsv")
    meta = _read(cfg, "hmm_model_metadata.tsv")
    bench = _read(cfg, "benchmark_metrics.tsv")
    lofo = _read(cfg, "benchmark_lofo_universal.tsv")
    motif = _read(cfg, "motif_conservation_summary.tsv")
    fp_fn = _read(cfg, "false_pos_neg.tsv")
    ranking = _read(cfg, "candidate_ranking.tsv")

    cat_counts = (clean["dataset_category"].value_counts().to_dict()
                  if not clean.empty else {})
    n_gold = cat_counts.get("gold", 0)

    # ---- data-driven conclusion logic --------------------------------------
    uni_roc = pd.to_numeric(pd.Series([_get(bench, "model", "universal", "roc_auc")]),
                            errors="coerce").iloc[0]
    fam_roc = pd.to_numeric(pd.Series([_get(bench, "model", "family", "roc_auc")]),
                            errors="coerce").iloc[0]
    full_roc = pd.to_numeric(pd.Series([_get(bench, "model", "hmm_arch_context", "roc_auc")]),
                             errors="coerce").iloc[0]
    hmmarch_roc = pd.to_numeric(pd.Series([_get(bench, "model", "hmm_arch", "roc_auc")]),
                                errors="coerce").iloc[0]
    lofo_recall = (pd.to_numeric(lofo["recall"], errors="coerce").dropna().mean()
                   if (not lofo.empty and "recall" in lofo.columns) else float("nan"))
    # Specificity under LOFO: fraction of hard negatives falsely hit per fold.
    lofo_fp_frac = float("nan")
    if not lofo.empty and "neg_false_hits" in lofo.columns and "n_neg" in lofo.columns:
        nfh = pd.to_numeric(lofo["neg_false_hits"], errors="coerce")
        nneg = pd.to_numeric(lofo["n_neg"], errors="coerce")
        frac = (nfh / nneg).dropna()
        if len(frac):
            lofo_fp_frac = float(frac.mean())
    uni_interp = str(_get(motif, "group_id", "universal", "interpretation", ""))

    conclusions = []
    if lofo_recall == lofo_recall and lofo_recall < 0.5:
        conclusions.append(
            "Experimentally characterized holins did not support a robust universal "
            "sequence-level HMM: under leave-one-family-out validation the universal model "
            f"recovered only {lofo_recall:.0%} of held-out families on average, indicating it "
            "primarily captured generic hydrophobic membrane-protein features rather than a "
            "transferable holin signature.")
    if lofo_fp_frac == lofo_fp_frac and lofo_fp_frac >= 0.3:
        conclusions.append(
            f"The universal HMM showed poor specificity: under leave-one-family-out it matched "
            f"{lofo_fp_frac:.0%} of hard-negative small membrane proteins on average per fold, "
            "consistent with the model capturing generic transmembrane character rather than a "
            "holin-specific signal. High recall alone is therefore not evidence of a usable "
            "universal model.")
    if fam_roc == fam_roc and uni_roc == uni_roc and fam_roc > uni_roc:
        conclusions.append(
            "Family-specific HMMs improved recovery of related holins relative to the universal "
            "model but, by construction, cannot identify distant holin families — consistent with "
            "extensive sequence diversity and a fragmented sequence space.")
    if full_roc == full_roc and hmmarch_roc == hmmarch_roc and full_roc >= hmmarch_roc:
        conclusions.append(
            "Combining HMM evidence with transmembrane topology and lysis-cassette context "
            "improved (or matched) candidate prioritization compared with sequence-only models.")
    if "hydrophobic" in uni_interp.lower():
        conclusions.append(
            "No universal linear amino-acid motif was detected across the curated holin set; "
            "conserved features were mostly architecture-level (small size, hydrophobicity, and "
            "predicted transmembrane segments).")
    if not conclusions:
        conclusions.append(
            "Results are reported descriptively; sample sizes and class balance should be "
            "considered before drawing strong conclusions.")

    # ---- assemble markdown --------------------------------------------------
    md = []
    md.append("# Benchmarking universal, topology-specific, and context-aware "
              "signatures for bacteriophage holin prediction\n")
    md.append("> **Caution.** This report is generated automatically from the pipeline outputs. "
              "Scores rank candidates; they do not prove holin function. With the bundled "
              "example data the sequences are SYNTHETIC placeholders — replace `data/example/` "
              "with curated data before any biological interpretation.\n")

    md.append("## 1. Background and rationale\n")
    md.append("Holins are small, hydrophobic phage membrane proteins that schedule host lysis. "
              "They are extremely diverse and may share function without deep sequence homology. "
              "This pipeline tests whether a single universal profile HMM can capture holins, and "
              "compares it against topology-specific HMMs, family-specific HMMs, and "
              "architecture/genomic-context models, while explicitly guarding against circular "
              "annotation and homology leakage.\n")

    md.append("## 2. Dataset summary\n")
    md.append("| class | n |\n| --- | --- |\n" +
              "".join(f"| {k} | {v} |\n" for k, v in cat_counts.items()))
    md.append(f"\nTotal QC flags raised: **{len(flags)}**.\n")
    if not flags.empty:
        md.append(_md_table(flags.groupby("flag").size().reset_index(name="count"),
                            max_rows=30))

    md.append("## 3. Gold positive curation summary\n")
    md.append(f"Curated gold positives: **{n_gold}**. "
              "Gold positives must rest on experimental evidence (see "
              "`docs/data_curation_guidelines.md`). The following gold records were flagged for "
              "missing evidence metadata and should be reviewed:\n")
    gold_flags = flags[flags["flag"].isin(
        ["gold_missing_citation", "gold_missing_evidence_type"])] if not flags.empty else pd.DataFrame()
    md.append(_md_table(gold_flags, max_rows=20) if not gold_flags.empty
              else "_No gold records missing citation/evidence metadata._\n")

    md.append("## 4. Sequence diversity and clustering\n")
    md.append(_md_table(clusters))
    md.append(_fig(cfg, "fig_cluster_fragmentation.png",
                   "Number of gold clusters vs identity threshold."))

    md.append("## 5. Universal HMM performance\n")
    md.append("Naive (potentially inflated) metrics:\n")
    md.append(_md_table(bench[bench["model"] == "universal"] if not bench.empty else bench))
    md.append("\n**Leave-one-family-out (honest generalization):**\n")
    md.append(_md_table(lofo))
    if lofo_recall == lofo_recall:
        md.append(f"\nMean held-out recall of the universal HMM: **{lofo_recall:.2f}**.\n")

    md.append("## 6. Topology-specific HMM performance\n")
    md.append(_md_table(bench[bench["model"] == "topology"] if not bench.empty else bench))

    md.append("## 7. Family-specific HMM performance\n")
    md.append(_md_table(bench[bench["model"] == "family"] if not bench.empty else bench))
    md.append("\nModel metadata:\n")
    md.append(_md_table(meta))

    md.append("## 8. Architecture / context-aware scoring performance\n")
    md.append(_md_table(bench[bench["model"].isin(
        ["architecture", "context", "hmm_arch", "hmm_arch_context"])] if not bench.empty else bench))
    md.append(_fig(cfg, "fig_benchmark_comparison.png", "Model comparison (naive regime)."))
    md.append(_fig(cfg, "fig_roc_curves.png", "ROC curves (naive)."))
    md.append(_fig(cfg, "fig_pr_curves.png", "Precision-recall curves (naive)."))

    md.append("## 9. Motif analysis\n")
    md.append(_md_table(motif))
    md.append(_fig(cfg, "logo_universal.png", "Sequence logo of the universal alignment."))

    md.append("## 10. False-positive and false-negative analysis\n")
    md.append(_md_table(fp_fn))
    md.append(_fig(cfg, "fig_confusion_matrix.png", "Confusion matrix (full model)."))

    md.append("## 11. Candidate ranking\n")
    md.append(_md_table(ranking, max_rows=25, cols=[
        "protein_id", "final_holin_score", "confidence_category", "hmm_score",
        "architecture_score", "context_score", "near_endolysin", "tmd_count", "length"]))
    md.append(_fig(cfg, "fig_candidate_scores.png", "Final candidate-score distribution."))

    md.append("## 12. Limitations\n")
    md.append(
        "- Sample sizes for experimentally validated holins are typically small; metrics "
        "(especially per-family AUCs) are unstable and reported with that caveat.\n"
        "- HMM models built from gold and evaluated on gold are circular; the naive regime is "
        "labeled as such, and the leave-one-family-out regime is the honest test.\n"
        "- TMD counts here may be imported predictions or a crude built-in estimate; the latter "
        "is flagged (`tool=builtin_kd_estimate`) and should be replaced with DeepTMHMM/Phobius.\n"
        "- Genomic context and 'small hydrophobic protein' are SUPPORTIVE signals, not proof.\n"
        "- With example data, sequences are synthetic and results are for pipeline validation only.\n")

    md.append("## 13. Recommended experimental follow-up\n")
    md.append(
        "- For high-confidence candidates, test holin function directly: lysis-timing assays, "
        "membrane depolarization/permeabilization, complementation of a holin-defective phage, "
        "and endolysin-dependent lysis.\n"
        "- Prioritize candidates that combine a family-HMM hit, compatible topology, and a "
        "plausible lysis cassette over those supported by hydrophobicity alone.\n")

    md.append("## Conclusions (data-driven)\n")
    for c in conclusions:
        md.append(f"- {c}\n")

    text = "\n".join(md) + "\n"
    out = cfg.output_dir / "report.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    log.info("Stage 12: wrote report to %s", out)
    return out
