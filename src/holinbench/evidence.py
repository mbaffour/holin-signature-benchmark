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

# Evidence-trigger taxonomy used for scoring (lowercased substrings). These are
# DEFAULTS; config/config.yaml -> literature.trigger_terms overrides them so
# scoring stays configurable (PROJECT_SPEC: "scoring must be configurable").
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
                 "co-expression with endolysin", "coexpression with endolysin"]
# Weak-only assay terms: NEVER sufficient alone; only count alongside a stronger
# trigger AND a holin subject (chloroform sensitivity is a generic lysis assay).
WEAK_ASSAY = ["chloroform"]
# Subject terms a trigger must co-occur with IN THE SAME SENTENCE to count.
HOLIN_SUBJECTS = ["holin", "pinholin", "antiholin", "lysis protein", "holin-like",
                  "gene s", "s105", "s107", "sar endolysin"]
ANNOTATION_ONLY = ["putative holin", "predicted holin", "annotated as holin",
                   "holin-like domain", "pfam", "blast", "homolog of"]


def _trigger_config(cfg: Config) -> dict:
    """Resolve the (configurable) evidence-trigger taxonomy.

    Reads literature.trigger_terms / weak_assay_terms / holin_subject_terms /
    exclusion_terms from config, falling back to the module-level defaults when a
    key is absent. Returns a dict of lowercased term lists.
    """
    lit = cfg.section("literature") if cfg is not None else {}
    tt = lit.get("trigger_terms") or {}
    def _lc(terms):
        return [str(t).lower() for t in (terms or [])]
    return {
        "genetic": _lc(tt.get("genetic")) or GENETIC,
        "functional": _lc(tt.get("functional_lysis")) or FUNCTIONAL_LYSIS,
        "membrane": _lc(tt.get("membrane")) or MEMBRANE,
        "endo": _lc(tt.get("endolysin_dependent")) or ENDOLYSIN_DEP,
        "weak_assay": _lc(lit.get("weak_assay_terms")) or WEAK_ASSAY,
        "subjects": _lc(lit.get("holin_subject_terms")) or HOLIN_SUBJECTS,
        "annotation_only": _lc(lit.get("exclusion_terms")) or ANNOTATION_ONLY,
    }


def _split_sentences(text: str) -> list[str]:
    if not text:
        return []
    # naive but adequate sentence splitter
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if s.strip()]


def _matches(text: str, terms: list[str]) -> list[str]:
    t = text.lower()
    return [term for term in terms if term in t]


def _has_subject(sentence_l: str, subjects: list[str]) -> bool:
    return any(subj in sentence_l for subj in subjects)


def _score_per_sentence(text: str, tc: dict) -> dict:
    """Subject-bound, per-sentence evidence aggregation.

    A trigger only counts when it appears IN THE SAME SENTENCE as a holin subject
    term. Weak-assay terms (chloroform) are tracked separately and only counted
    when that same sentence also carries a stronger trigger. Returns the
    document-level category sets used by ``_evidence_score`` plus bookkeeping.
    """
    genetic: set[str] = set()
    functional: set[str] = set()
    membrane: set[str] = set()
    endo: set[str] = set()
    weak_assay: set[str] = set()
    subject_seen = False
    best_sentence = ""
    for sent in _split_sentences(text):
        sl = sent.lower()
        if not _has_subject(sl, tc["subjects"]):
            continue
        subject_seen = True
        g = [t for t in tc["genetic"] if t in sl]
        f = [t for t in tc["functional"] if t in sl]
        m = [t for t in tc["membrane"] if t in sl]
        e = [t for t in tc["endo"] if t in sl]
        w = [t for t in tc["weak_assay"] if t in sl]
        strong = g + f + m + e
        # weak-assay terms only count when corroborated by a stronger trigger in
        # the SAME (subject-bound) sentence.
        if w and strong:
            weak_assay.update(w)
        genetic.update(g)
        functional.update(f)
        membrane.update(m)
        endo.update(e)
        if (strong or (w and strong)) and not best_sentence:
            best_sentence = sent[:400]
    return {
        "genetic": sorted(genetic),
        "functional": sorted(functional),
        "membrane": sorted(membrane),
        "endo": sorted(endo),
        "weak_assay": sorted(weak_assay),
        "subject_seen": subject_seen,
        "best_sentence": best_sentence,
    }


def _find_candidate_names(text: str, patterns: list[str]) -> list[str]:
    found: dict[str, int] = {}
    for pat in patterns:
        for m in re.finditer(pat, text, flags=re.IGNORECASE):
            name = m.group(0).strip()
            found[name] = found.get(name, 0) + 1
    # rank by frequency, keep most common distinct names
    return [n for n, _ in sorted(found.items(), key=lambda kv: -kv[1])][:6]


def _evidence_score(genetic, functional, membrane, endo, annotation_only,
                    subject_seen=None) -> int:
    """Map subject-bound evidence sets to the documented 0-5 ladder.

    5 = genetic + (functional/membrane/endo); 4 = direct functional lysis, or
    membrane corroborated by genetic/functional; 3 = strong but single-line /
    uncorroborated (membrane-only or endolysin-dependent-only); 2 = a holin claim
    (subject present) with no qualifying trigger, but not purely annotation-only;
    1 = annotation/prediction only; 0 = none.

    Weak-assay terms (chloroform) are deliberately NOT an input here: by the time
    sets reach this function they have already been filtered to only count when
    corroborated by a stronger trigger, so they can never lift the score alone.
    Annotation-only ALWAYS caps the score at 1 when there is no trigger evidence.
    """
    has_g = bool(genetic)
    has_func = bool(functional)
    has_mem = bool(membrane)
    has_endo = bool(endo)
    has_trigger = has_g or has_func or has_mem or has_endo

    # (C4) Annotation terms present with NO subject-bound trigger -> cap at 1,
    # regardless of unrelated lysis terms elsewhere.
    if annotation_only and not has_trigger:
        return 1

    if has_g and (has_func or has_mem or has_endo):
        return 5
    # Direct functional lysis evidence, or membrane corroborated by genetic/
    # functional, reaches 4. Membrane ALONE is only a 3.
    if has_func or (has_mem and (has_g or has_func)):
        return 4
    if has_mem or has_endo:
        return 3
    # No qualifying trigger sentence. If a holin subject was nonetheless asserted
    # and this is not purely annotation-only -> a claim with indirect evidence (2).
    if subject_seen and not annotation_only:
        return 2
    if annotation_only:
        return 1
    return 0


def _classify(score: int, text_l: str, annotation_only: bool):
    is_pin = "pinholin" in text_l or "sar endolysin" in text_l or "pin-holin" in text_l
    is_like = "holin-like" in text_l or "lysis protein" in text_l
    # Apply pinholin / holin-like qualifiers consistently across the 5 AND 4 tiers.
    if score >= 4:
        cls = ("experimentally_validated_pinholin" if is_pin else
               ("experimentally_validated_holin_like_lysis_protein" if is_like else
                "experimentally_validated_holin"))
        if score >= 5:
            reason = "direct genetic + functional/membrane evidence"
        else:
            reason = "direct functional lysis or corroborated membrane evidence"
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


def _append_qc(row_flags: str, *new_flags: str) -> str:
    """Append new QC flags to an existing semicolon-joined qc_flags string."""
    existing = [f.strip() for f in str(row_flags or "").split(";") if f.strip()]
    for f in new_flags:
        if f and f not in existing:
            existing.append(f)
    return "; ".join(existing)


def _apply_qc_post_pass(candidates: pd.DataFrame) -> pd.DataFrame:
    """Post-pass QC flags required by PROJECT_SPEC (~line 231).

    Adds, on top of the per-row flags already present:
      - duplicate_across_papers: same gene/protein name (or identical mapped
        sequence) appears in more than one distinct paper.
      - conflicting_name_for_accession: one accession maps to >1 distinct
        gene_name.
      - review_citation_only: a review_article whose only holin mention is
        background (holin_only_in_background also present).
    """
    if candidates.empty:
        return candidates

    qc = candidates["qc_flags"].tolist()

    def _norm(s):
        return str(s or "").strip().lower()

    # (a) duplicate protein/gene name (or identical mapped sequence) across papers.
    for key_col in ("gene_name", "sequence"):
        if key_col not in candidates.columns:
            continue
        groups: dict[str, set] = {}
        for i, r in candidates.iterrows():
            val = _norm(r.get(key_col, ""))
            if not val or val in ("(unspecified)",):
                continue
            groups.setdefault(val, set()).add(_norm(r.get("paper_id", "")))
        dup_vals = {v for v, papers in groups.items() if len(papers) > 1}
        if dup_vals:
            for pos, (i, r) in enumerate(candidates.iterrows()):
                if _norm(r.get(key_col, "")) in dup_vals:
                    qc[pos] = _append_qc(qc[pos], "duplicate_across_papers")

    # (b) conflicting gene_names for the same accession.
    if "accession" in candidates.columns:
        acc_names: dict[str, set] = {}
        for i, r in candidates.iterrows():
            acc = _norm(r.get("accession", ""))
            if not acc:
                continue
            acc_names.setdefault(acc, set()).add(_norm(r.get("gene_name", "")))
        conflict_acc = {a for a, names in acc_names.items() if len(names) > 1}
        if conflict_acc:
            for pos, (i, r) in enumerate(candidates.iterrows()):
                if _norm(r.get("accession", "")) in conflict_acc:
                    qc[pos] = _append_qc(qc[pos], "conflicting_name_for_accession")

    # (c) review-citation-only: review_article AND holin_only_in_background.
    for pos, flags in enumerate(qc):
        present = {f.strip() for f in str(flags).split(";")}
        if "review_article" in present and "holin_only_in_background" in present:
            qc[pos] = _append_qc(qc[pos], "review_citation_only")

    candidates = candidates.copy()
    candidates["qc_flags"] = qc
    return candidates


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
    tc = _trigger_config(cfg)

    sent_rows, cand_rows = [], []
    for _, row in papers.iterrows():
        text, section = _paper_text(cfg, row, fulltext_dir)
        text_l = text.lower()
        if "holin" not in text_l and "lysis" not in text_l:
            continue

        # PER-SENTENCE, subject-bound aggregation. A trigger only counts when it
        # co-occurs with a holin subject term in the SAME sentence; this prevents
        # an unrelated "required for lysis"/"chloroform"/"depolariz" elsewhere in
        # the document from promoting an annotation-only paper.
        ev = _score_per_sentence(text, tc)
        genetic = ev["genetic"]
        functional = ev["functional"]
        membrane = ev["membrane"]
        endo = ev["endo"]
        weak_assay = ev["weak_assay"]
        annotation_only = _matches(text, tc["annotation_only"])
        # Only subject-bound triggers; weak-assay terms are appended for the audit
        # trail but already filtered to corroborated occurrences in ev.
        all_triggers = genetic + functional + membrane + endo + weak_assay
        best_sentence = ev["best_sentence"]

        # "holin only in introduction" heuristic: holin in title/abstract but no
        # subject-bound trigger anywhere -> background mention.
        holin_in_intro_only = ("holin" in (str(row.get("title", "")) + " " +
                                           str(row.get("abstract", ""))).lower()
                               and not all_triggers)

        # Emit per-sentence evidence rows (subject-bound only).
        for sent in _split_sentences(text):
            sl = sent.lower()
            if not _has_subject(sl, tc["subjects"]):
                continue
            hit = ([t for t in tc["genetic"] if t in sl] +
                   [t for t in tc["functional"] if t in sl] +
                   [t for t in tc["membrane"] if t in sl] +
                   [t for t in tc["endo"] if t in sl])
            weak_hit = [t for t in tc["weak_assay"] if t in sl]
            if weak_hit and hit:
                hit += weak_hit
            if hit:
                sent_rows.append({"paper_id": row["paper_id"], "pmid": row.get("pmid", ""),
                                  "trigger_terms": "; ".join(hit),
                                  "section": section, "evidence_sentence": sent[:400]})

        score = _evidence_score(genetic, functional, membrane, endo, annotation_only,
                                subject_seen=ev["subject_seen"])
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
    candidates = _apply_qc_post_pass(candidates)
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
    else:
        # The default is require_manual_verification_for_gold=True (safe). When a
        # user explicitly turns it off, human review is being bypassed; warn LOUDLY.
        log.warning(
            "*** WARNING: require_manual_verification_for_gold is FALSE. "
            "Human verification is being BYPASSED -- candidates with "
            "evidence_score >= %d will be exported as gold WITHOUT manual review. "
            "This is NOT the safe default; results may include annotation-only or "
            "misread proteins. Set require_manual_verification_for_gold: true to "
            "restore the conservative gate. ***", min_score)
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
