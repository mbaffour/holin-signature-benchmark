"""Stage 0b — evidence extraction, conservative classification, review & export.

Reads the papers from Stage 0a, extracts candidate gene/protein names and
experimental-evidence sentences, and classifies each candidate using DELIBERATELY
CONSERVATIVE rules. The word "holin" never, by itself, qualifies a protein as a
gold positive. Only candidates with direct functional/genetic evidence
(evidence_score >= 4) AND manual verification are exported to gold_holins.csv.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .config import Config
from .utils import log

# Evidence-trigger taxonomy used for scoring (lowercased substrings).
GENETIC = ["amber mutant", "nonsense mutant", "knockout", "knock-out",
           "deletion mutant", "complementation", "complements", "complemented"]
FUNCTIONAL_LYSIS = ["required for lysis", "essential for lysis", "necessary for lysis",
                    "sufficient for lysis", "restored lysis", "failed to lyse",
                    "delayed lysis", "rapid lysis", "premature lysis", "lysis timing",
                    "timing of lysis", "abolished lysis", "no lysis"]
MEMBRANE = ["membrane depolarization", "membrane permeabilization", "depolariz",
            "permeabiliz", "proton motive force", "membrane potential",
            "triggered the membrane"]
ENDOLYSIN_DEP = ["endolysin dependent", "endolysin-dependent",
                 "co-expression with endolysin", "coexpression with endolysin",
                 "chloroform"]
ANNOTATION_ONLY = ["putative holin", "predicted holin", "annotated as holin",
                   "holin-like domain", "pfam", "blast", "homolog of"]


def _split_sentences(text: str) -> list[str]:
    if not text:
        return []
    # naive but adequate sentence splitter
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if s.strip()]


def _matches(text: str, terms: list[str]) -> list[str]:
    t = text.lower()
    return [term for term in terms if term in t]


def _find_candidate_names(text: str, patterns: list[str]) -> list[str]:
    found: dict[str, int] = {}
    for pat in patterns:
        for m in re.finditer(pat, text, flags=re.IGNORECASE):
            name = m.group(0).strip()
            found[name] = found.get(name, 0) + 1
    # rank by frequency, keep most common distinct names
    return [n for n, _ in sorted(found.items(), key=lambda kv: -kv[1])][:6]


def _evidence_score(genetic, functional, membrane, endo, annotation_only) -> int:
    has_g = bool(genetic)
    has_func = bool(functional)
    has_mem = bool(membrane)
    has_endo = bool(endo)
    if has_g and (has_func or has_mem or has_endo):
        return 5
    if has_mem or (has_func and has_endo):
        return 4
    if has_func or has_endo:
        return 3
    if annotation_only:
        return 1
    return 0


def _classify(score: int, text_l: str, annotation_only: bool):
    is_pin = "pinholin" in text_l or "sar endolysin" in text_l or "pin-holin" in text_l
    is_like = "holin-like" in text_l or "lysis protein" in text_l
    if score >= 5:
        cls = "experimentally_validated_pinholin" if is_pin else "experimentally_validated_holin"
        reason = "direct genetic + functional/membrane evidence"
    elif score == 4:
        cls = ("experimentally_validated_pinholin" if is_pin else
               ("experimentally_validated_holin_like_lysis_protein" if is_like else
                "experimentally_validated_holin"))
        reason = "direct functional lysis or membrane evidence"
    elif score == 3:
        cls = "experimentally_validated_holin"
        reason = "strong experimental evidence; sequence mapping/identity to verify"
    elif score == 2:
        cls = "literature_claimed_holin_but_no_direct_evidence"
        reason = "claim present but evidence indirect/unclear"
    elif score == 1:
        cls = "database_annotated_only" if annotation_only else "putative_holin"
        reason = "annotation/prediction/neighborhood support only"
    else:
        cls = "insufficient_evidence"
        reason = "no experimental evidence detected in available text"
    return cls, reason


def _paper_text(cfg: Config, row, fulltext_dir: Path) -> tuple[str, str]:
    """Return (combined_text, section_label_of_best). Includes fulltext if present."""
    parts = [str(row.get("title", "")), str(row.get("abstract", ""))]
    ftpath = fulltext_dir / f"{row['paper_id']}.txt"
    if ftpath.exists():
        parts.append(ftpath.read_text(encoding="utf-8"))
    return " ".join(parts), "abstract+fulltext" if ftpath.exists() else "title+abstract"


def run_extract(cfg: Config) -> pd.DataFrame:
    lit = cfg.section("literature")
    out_dir = cfg.resolve(lit.get("output_dir", "results/literature"))
    fulltext_dir = out_dir / "fulltext"
    search_path = out_dir / "literature_search_results.tsv"
    if not search_path.exists():
        log.warning("No literature_search_results.tsv; run `literature-search` first.")
        return pd.DataFrame()
    papers = pd.read_csv(search_path, sep="\t", dtype=str, keep_default_na=False)
    patterns = lit.get("gene_name_patterns", [r"\bholin\b", r"\bpinholin\b",
                                              r"\bgp\d+\b", r"\bS10[57]\b"])

    sent_rows, cand_rows = [], []
    for _, row in papers.iterrows():
        text, section = _paper_text(cfg, row, fulltext_dir)
        text_l = text.lower()
        if "holin" not in text_l and "lysis" not in text_l:
            continue
        genetic = _matches(text, GENETIC)
        functional = _matches(text, FUNCTIONAL_LYSIS)
        membrane = _matches(text, MEMBRANE)
        endo = _matches(text, ENDOLYSIN_DEP)
        annotation_only = _matches(text, ANNOTATION_ONLY)
        all_triggers = genetic + functional + membrane + endo

        # "holin only in introduction" heuristic: holin in title/abstract but no
        # trigger anywhere -> background mention.
        holin_in_intro_only = ("holin" in (str(row.get("title", "")) + " " +
                                           str(row.get("abstract", ""))).lower()
                               and not all_triggers)

        # evidence sentences
        best_sentence = ""
        for sent in _split_sentences(text):
            hit_terms = _matches(sent, all_triggers)
            if hit_terms:
                sent_rows.append({"paper_id": row["paper_id"], "pmid": row.get("pmid", ""),
                                  "trigger_terms": "; ".join(hit_terms),
                                  "section": section, "evidence_sentence": sent[:400]})
                if not best_sentence:
                    best_sentence = sent[:400]

        score = _evidence_score(genetic, functional, membrane, endo, annotation_only)
        cls, reason = _classify(score, text_l, bool(annotation_only))
        names = _find_candidate_names(text, patterns) or ["(unspecified)"]

        for name in names:
            qc = []
            if str(row.get("is_review", "")).lower() in ("true", "1"):
                qc.append("review_article")
            if holin_in_intro_only:
                qc.append("holin_only_in_background")
            if not all_triggers:
                qc.append("annotation_or_claim_only")
            cand_rows.append({
                "candidate_id": f"{row['paper_id']}__{re.sub(r'[^A-Za-z0-9]', '', name)[:20]}",
                "paper_id": row["paper_id"], "pmid": row.get("pmid", ""),
                "pmcid": row.get("pmcid", ""), "doi": row.get("doi", ""),
                "title": row.get("title", ""), "year": row.get("year", ""),
                "phage_name": "", "host_name": "", "gene_name": name, "protein_name": name,
                "accession": "", "sequence": "",
                "evidence_sentence": best_sentence,
                "evidence_section": section,
                "evidence_type": "; ".join(all_triggers[:5]) or "none",
                "evidence_score": score,
                "confidence_class": cls,
                "reason_for_classification": reason,
                "qc_flags": "; ".join(qc),
                "requires_manual_review": True,
                "manually_verified": False,
                "curator_notes": "",
            })

    sentences = pd.DataFrame(sent_rows, columns=["paper_id", "pmid", "trigger_terms",
                                                 "section", "evidence_sentence"])
    candidates = pd.DataFrame(cand_rows)
    sentences.to_csv(out_dir / "candidate_evidence_sentences.tsv", sep="\t", index=False)
    candidates.to_csv(out_dir / "candidate_holin_literature_table.tsv", sep="\t", index=False)
    n_strong = int((candidates["evidence_score"] >= 4).sum()) if not candidates.empty else 0
    log.info("Stage 0: extracted %d candidate mention(s) from %d paper(s); "
             "%d with strong (>=4) evidence (require manual review).",
             len(candidates), len(papers), n_strong)
    return candidates


def run_prepare_review(cfg: Config) -> Path:
    lit = cfg.section("literature")
    out_dir = cfg.resolve(lit.get("output_dir", "results/literature"))
    cand_path = out_dir / "candidate_holin_literature_table.tsv"
    if not cand_path.exists():
        log.warning("No candidate table; run `extract-evidence` first.")
        return out_dir / "manual_review_template.tsv"
    cand = pd.read_csv(cand_path, sep="\t", dtype=str, keep_default_na=False)

    # incorporate sequence mapping if available
    seqmap_path = out_dir / "sequence_mapping_table.tsv"
    if seqmap_path.exists():
        sm = pd.read_csv(seqmap_path, sep="\t", dtype=str, keep_default_na=False)
        if not sm.empty:
            cand = cand.merge(sm[["candidate_id", "source_accession", "sequence"]],
                              on="candidate_id", how="left", suffixes=("", "_mapped"))
            cand["accession"] = cand.get("source_accession", cand.get("accession", ""))
            cand["sequence"] = cand["sequence_mapped"].fillna(cand.get("sequence", "")) \
                if "sequence_mapped" in cand.columns else cand.get("sequence", "")

    review_cols = ["candidate_id", "title", "pmid", "doi", "phage_name", "gene_name",
                   "evidence_sentence", "evidence_type", "evidence_score",
                   "confidence_class", "accession", "sequence", "qc_flags",
                   "manually_verified", "curator_decision", "curator_notes"]
    for c in ["curator_decision"]:
        if c not in cand.columns:
            cand[c] = ""
    template = cand.reindex(columns=review_cols, fill_value="")
    template.to_csv(out_dir / "manual_review_template.tsv", sep="\t", index=False)

    _write_curation_summary(cfg, out_dir, cand)
    log.info("Stage 0: wrote manual_review_template.tsv (%d rows) + curation_summary.md",
             len(template))
    return out_dir / "manual_review_template.tsv"


def _write_curation_summary(cfg: Config, out_dir: Path, cand: pd.DataFrame) -> None:
    search_path = out_dir / "literature_search_results.tsv"
    papers = pd.read_csv(search_path, sep="\t", dtype=str, keep_default_na=False) \
        if search_path.exists() else pd.DataFrame()
    n_papers = len(papers)
    n_abstract = int((papers["abstract"].str.len() > 0).sum()) if not papers.empty else 0
    n_fulltext = int((papers["has_fulltext"].astype(str).str.lower()
                      .isin(["true", "1"])).sum()) if not papers.empty else 0
    scores = pd.to_numeric(cand["evidence_score"], errors="coerce") if not cand.empty else pd.Series([], dtype=float)
    n_accept = int((scores >= 4).sum())
    n_weak = int(((scores >= 1) & (scores < 4)).sum())
    n_reject = int((scores < 1).sum())
    cls_counts = cand["confidence_class"].value_counts().to_dict() if not cand.empty else {}
    high_value = cand[scores >= 4][["title", "pmid", "gene_name", "evidence_score"]].head(25) \
        if not cand.empty else pd.DataFrame()

    lines = ["# Curation summary (Stage 0)\n",
             "> Automated extraction is a CURATION AID, not an authority. Every candidate "
             "below MUST be manually checked before use as a gold positive.\n",
             f"- Papers retrieved: **{n_papers}**",
             f"- Papers with abstracts: **{n_abstract}**",
             f"- Open-access full texts retrieved: **{n_fulltext}**",
             f"- Candidate holin mentions extracted: **{len(cand)}**",
             f"- Candidates with strong evidence (score >= 4): **{n_accept}**",
             f"- Candidates downgraded to weak (1-3): **{n_weak}**",
             f"- Candidates rejected / insufficient (0): **{n_reject}**\n",
             "## Confidence-class breakdown\n"]
    for k, v in sorted(cls_counts.items()):
        lines.append(f"- {k}: {v}")
    lines.append("\n## High-value papers requiring manual reading\n")
    if high_value.empty:
        lines.append("_None reached score >= 4 automatically; review the full table._")
    else:
        lines.append("| title | pmid | gene | score |")
        lines.append("| --- | --- | --- | --- |")
        for _, r in high_value.iterrows():
            t = str(r["title"]).replace("|", "\\|")[:70]
            lines.append(f"| {t} | {r['pmid']} | {r['gene_name']} | {r['evidence_score']} |")
    (out_dir / "curation_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_export_curated(cfg: Config) -> Path:
    lit = cfg.section("literature")
    out_dir = cfg.resolve(lit.get("output_dir", "results/literature"))
    require_manual = bool(lit.get("require_manual_verification_for_gold", True))
    min_score = int(lit.get("min_evidence_score_for_gold", 4))

    review_path = out_dir / "manual_review_template.tsv"
    cand_path = out_dir / "candidate_holin_literature_table.tsv"
    src = review_path if review_path.exists() else cand_path
    if not src.exists():
        log.warning("Nothing to export; run extract-evidence/prepare-review first.")
        return out_dir / "accepted_gold_holins.tsv"
    df = pd.read_csv(src, sep="\t", dtype=str, keep_default_na=False)
    scores = pd.to_numeric(df.get("evidence_score", pd.Series()), errors="coerce").fillna(0)
    verified = df.get("manually_verified", pd.Series([""] * len(df))).astype(str).str.lower()
    is_verified = verified.isin(["true", "1", "yes"])

    accept_mask = scores >= min_score
    if require_manual:
        accept_mask = accept_mask & is_verified
    accepted = df[accept_mask].copy()
    rejected = df[~accept_mask].copy()
    accepted.to_csv(out_dir / "accepted_gold_holins.tsv", sep="\t", index=False)
    rejected.to_csv(out_dir / "rejected_or_weak_holins.tsv", sep="\t", index=False)

    if accepted.empty:
        log.warning("export-curated: 0 candidates met the gold criteria "
                    "(score >= %d%s). gold_holins.csv NOT modified — this is the "
                    "safe default. Manually set manually_verified=true in "
                    "manual_review_template.tsv for genuine holins, then re-run.",
                    min_score, " AND manually_verified=true" if require_manual else "")
        return out_dir / "accepted_gold_holins.tsv"

    # Build a gold_holins.csv-shaped export (curator still owns the canonical file).
    gold = pd.DataFrame({
        "protein_id": accepted["candidate_id"],
        "protein_name": accepted.get("protein_name", accepted.get("gene_name", "")),
        "phage_name": accepted.get("phage_name", ""),
        "host": accepted.get("host_name", ""),
        "accession": accepted.get("accession", ""),
        "sequence": accepted.get("sequence", ""),
        "evidence_type": accepted.get("evidence_type", ""),
        "citation": accepted.get("pmid", ""),
        "holin_type": "unknown",
        "family_label": "",
        "notes": "exported from Stage 0 (verify before use): " + accepted.get("confidence_class", ""),
    })
    export_path = out_dir / "gold_holins_from_literature.csv"
    gold.to_csv(export_path, index=False)
    log.info("export-curated: wrote %d curated gold candidate(s) to %s. "
             "Review and merge into the canonical data/example/gold_holins.csv.",
             len(gold), export_path)
    return export_path
