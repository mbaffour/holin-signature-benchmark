from pathlib import Path

import pandas as pd

from holinbench import evidence
from holinbench.config import Config


def test_evidence_score_genetic_plus_functional_is_5():
    s = evidence._evidence_score(["amber mutant"], ["delayed lysis"], [], [], [])
    assert s == 5


def test_evidence_score_membrane_only_is_3():
    # Fix #2: membrane-only evidence is uncorroborated -> 3 (NOT the gold 4).
    s = evidence._evidence_score([], [], ["membrane depolarization"], [], [])
    assert s == 3


def test_evidence_score_membrane_plus_genetic_is_5():
    s = evidence._evidence_score(["amber mutant"], [], ["depolariz"], [], [])
    assert s == 5


def test_evidence_score_functional_only_is_4():
    # Fix #2: direct functional lysis evidence reaches 4.
    s = evidence._evidence_score([], ["required for lysis"], [], [], [])
    assert s == 4


def test_evidence_score_endo_only_is_3():
    s = evidence._evidence_score([], [], [], ["endolysin-dependent"], [])
    assert s == 3


def test_evidence_score_annotation_only_is_1():
    s = evidence._evidence_score([], [], [], [], ["putative holin"])
    assert s == 1


def test_evidence_score_none_is_0():
    assert evidence._evidence_score([], [], [], [], []) == 0


def test_classify_pinholin():
    cls, _ = evidence._classify(5, "the pinholin s21 depolarized the membrane", False)
    assert cls == "experimentally_validated_pinholin"


def test_classify_annotation_only():
    cls, _ = evidence._classify(1, "a putative holin annotated by pfam", True)
    assert cls == "database_annotated_only"


def test_evidence_score_claim_no_trigger_is_2():
    # Fix #3: a holin subject asserted (subject_seen) but no qualifying trigger and
    # not annotation-only -> 2 (literature_claimed_holin_but_no_direct_evidence).
    s = evidence._evidence_score([], [], [], [], [], subject_seen=True)
    assert s == 2
    cls, _ = evidence._classify(s, "this holin is involved in lysis", False)
    assert cls == "literature_claimed_holin_but_no_direct_evidence"


def test_evidence_score_annotation_caps_even_with_subject():
    # Fix #4: annotation terms present, no subject-bound trigger -> capped at 1,
    # even when subject_seen is True.
    s = evidence._evidence_score([], [], [], [], ["putative holin"], subject_seen=True)
    assert s == 1


def test_classify_pinholin_qualifier_applies_at_score_5():
    # Fix #5: pinholin qualifier must apply at score 5 too, not only at 4.
    cls5, _ = evidence._classify(5, "the pinholin depolarized the membrane", False)
    assert cls5 == "experimentally_validated_pinholin"
    cls4, _ = evidence._classify(4, "the pinholin depolarized the membrane", False)
    assert cls4 == "experimentally_validated_pinholin"


def test_classify_holin_like_qualifier_applies_at_score_5():
    cls5, _ = evidence._classify(5, "this holin-like lysis protein is required", False)
    assert cls5 == "experimentally_validated_holin_like_lysis_protein"


def test_score_per_sentence_requires_subject_binding():
    tc = evidence._trigger_config(None)
    # Trigger and subject in DIFFERENT sentences -> no qualifying evidence.
    text = ("ProteinX was required for lysis of the host. "
            "Separately, the holin was annotated by Pfam.")
    ev = evidence._score_per_sentence(text, tc)
    assert ev["functional"] == []
    assert ev["genetic"] == []
    # Trigger and subject in the SAME sentence -> counts.
    text2 = "The holin was required for lysis of the host."
    ev2 = evidence._score_per_sentence(text2, tc)
    assert "required for lysis" in ev2["functional"]


def test_score_per_sentence_weak_assay_never_alone():
    tc = evidence._trigger_config(None)
    # chloroform with a holin subject but NO stronger trigger -> not counted.
    text = "The holin showed chloroform sensitivity in the assay."
    ev = evidence._score_per_sentence(text, tc)
    assert ev["weak_assay"] == []
    s = evidence._evidence_score(ev["genetic"], ev["functional"], ev["membrane"],
                                 ev["endo"], [], subject_seen=ev["subject_seen"])
    assert s < 4


def test_find_candidate_names():
    text = "The holin S105 and gene S were required; gp17 also studied. holin holin"
    names = evidence._find_candidate_names(text, [r"\bS105\b", r"\bgp\d+\b", r"\bholin\b"])
    assert "holin" in [n.lower() for n in names]


def test_export_gate_refuses_unverified(tmp_path: Path):
    out = tmp_path / "results" / "literature"
    out.mkdir(parents=True)
    # a candidate with strong evidence but NOT manually verified
    pd.DataFrame([{
        "candidate_id": "c1", "title": "t", "pmid": "1", "doi": "",
        "phage_name": "", "gene_name": "S", "evidence_sentence": "s",
        "evidence_type": "delayed lysis", "evidence_score": 5,
        "confidence_class": "experimentally_validated_holin", "accession": "",
        "sequence": "", "qc_flags": "", "manually_verified": False,
        "curator_decision": "", "curator_notes": "",
    }]).to_csv(out / "manual_review_template.tsv", sep="\t", index=False)

    cfg = Config({"literature": {"output_dir": str(out),
                                 "require_manual_verification_for_gold": True,
                                 "min_evidence_score_for_gold": 4}},
                 tmp_path / "config" / "config.yaml")
    result = evidence.run_export_curated(cfg)
    accepted = pd.read_csv(out / "accepted_gold_holins.tsv", sep="\t")
    assert accepted.empty                      # gate held: unverified -> not exported
    assert not (out / "gold_holins_from_literature.csv").exists()


def test_export_gate_allows_verified(tmp_path: Path):
    out = tmp_path / "results" / "literature"
    out.mkdir(parents=True)
    pd.DataFrame([{
        "candidate_id": "c1", "title": "t", "pmid": "1", "doi": "",
        "phage_name": "Lambda", "gene_name": "S", "protein_name": "S holin",
        "host_name": "Ecoli", "evidence_sentence": "s", "evidence_type": "delayed lysis",
        "evidence_score": 5, "confidence_class": "experimentally_validated_holin",
        "accession": "ACC1", "sequence": "MKLLIV", "qc_flags": "",
        "manually_verified": True, "curator_decision": "accept", "curator_notes": "",
    }]).to_csv(out / "manual_review_template.tsv", sep="\t", index=False)

    cfg = Config({"literature": {"output_dir": str(out),
                                 "require_manual_verification_for_gold": True,
                                 "min_evidence_score_for_gold": 4}},
                 tmp_path / "config" / "config.yaml")
    evidence.run_export_curated(cfg)
    gold = pd.read_csv(out / "gold_holins_from_literature.csv")
    assert len(gold) == 1
    assert gold.iloc[0]["protein_id"] == "c1"


# --------------------------------------------------------------------------- #
# End-to-end run_extract conservatism tests (Stage 0 must be CONSERVATIVE).
# --------------------------------------------------------------------------- #
SEARCH_COLS = ["paper_id", "query_group", "query", "source", "pmid", "pmcid", "doi",
               "title", "authors", "year", "journal", "abstract", "is_review",
               "has_fulltext", "url"]


def _write_search(out: Path, papers: list[dict]) -> None:
    rows = []
    for i, p in enumerate(papers):
        row = {c: "" for c in SEARCH_COLS}
        row["paper_id"] = p.get("paper_id", f"lit_p{i}")
        row["pmid"] = p.get("pmid", str(i))
        row["title"] = p.get("title", "")
        row["abstract"] = p.get("abstract", "")
        row["is_review"] = str(p.get("is_review", False))
        rows.append(row)
    pd.DataFrame(rows, columns=SEARCH_COLS).to_csv(
        out / "literature_search_results.tsv", sep="\t", index=False)


def _extract_cfg(out: Path, tmp_path: Path) -> Config:
    return Config({"literature": {"output_dir": str(out)}},
                  tmp_path / "config" / "config.yaml")


def _run_extract(out: Path, tmp_path: Path) -> pd.DataFrame:
    cfg = _extract_cfg(out, tmp_path)
    return evidence.run_extract(cfg)


def test_c1_annotation_only_with_unrelated_lysis_sentence_scores_at_most_1(tmp_path: Path):
    """C1 regression: a putative-holin (annotation-only) abstract that ALSO
    contains an unrelated 'required for lysis' sentence about a DIFFERENT protein
    must NOT be promoted. Document-level OR matching used to score this 3-4."""
    out = tmp_path / "results" / "literature"
    out.mkdir(parents=True)
    abstract = ("This putative holin was identified by BLAST and annotated as "
                "holin-like by Pfam. In an unrelated system, ProteinX was "
                "required for lysis of the host cell.")
    _write_search(out, [{"title": "A putative holin from a BLAST survey",
                         "abstract": abstract}])
    cand = _run_extract(out, tmp_path)
    assert not cand.empty
    scores = pd.to_numeric(cand["evidence_score"], errors="coerce")
    assert scores.max() <= 1, f"annotation-only paper over-promoted: {scores.tolist()}"


def test_lone_chloroform_with_subject_not_gold(tmp_path: Path):
    out = tmp_path / "results" / "literature"
    out.mkdir(parents=True)
    abstract = "The holin showed chloroform sensitivity in a standard lysis assay."
    _write_search(out, [{"title": "Holin chloroform assay", "abstract": abstract}])
    cand = _run_extract(out, tmp_path)
    scores = pd.to_numeric(cand["evidence_score"], errors="coerce")
    assert scores.max() < 4, f"lone chloroform reached gold: {scores.tolist()}"


def test_lone_depolariz_with_subject_not_gold(tmp_path: Path):
    out = tmp_path / "results" / "literature"
    out.mkdir(parents=True)
    # Membrane-only (depolariz) with a holin subject -> score 3, never >= 4.
    abstract = "Expression of the holin depolarized the inner membrane."
    _write_search(out, [{"title": "Holin membrane study", "abstract": abstract}])
    cand = _run_extract(out, tmp_path)
    scores = pd.to_numeric(cand["evidence_score"], errors="coerce")
    assert scores.max() < 4, f"lone depolariz reached gold: {scores.tolist()}"


def test_per_sentence_positive_amber_mutant_delayed_lysis_scores_5(tmp_path: Path):
    out = tmp_path / "results" / "literature"
    out.mkdir(parents=True)
    abstract = ("An amber mutant of the holin showed delayed lysis of the "
                "infected culture, confirming its role in lysis timing.")
    _write_search(out, [{"title": "Amber mutant of a holin", "abstract": abstract}])
    cand = _run_extract(out, tmp_path)
    scores = pd.to_numeric(cand["evidence_score"], errors="coerce")
    assert scores.max() == 5, f"expected 5, got {scores.tolist()}"


def test_score_2_reachable_claim_without_trigger(tmp_path: Path):
    out = tmp_path / "results" / "literature"
    out.mkdir(parents=True)
    # A holin claim, no annotation terms, no qualifying trigger sentence -> 2.
    abstract = "We report a novel holin from this phage and discuss its biology."
    _write_search(out, [{"title": "A novel holin", "abstract": abstract}])
    cand = _run_extract(out, tmp_path)
    scores = pd.to_numeric(cand["evidence_score"], errors="coerce")
    assert 2 in set(scores.tolist()), f"score 2 not reachable: {scores.tolist()}"
    cls = set(cand[scores == 2]["confidence_class"])
    assert "literature_claimed_holin_but_no_direct_evidence" in cls


def test_qc_duplicate_across_papers(tmp_path: Path):
    out = tmp_path / "results" / "literature"
    out.mkdir(parents=True)
    # Same holin subject + functional evidence in two distinct papers -> the
    # shared gene_name should be flagged duplicate_across_papers.
    ab = "The holin was required for lysis."
    _write_search(out, [
        {"paper_id": "lit_a", "pmid": "10", "title": "Holin paper A", "abstract": ab},
        {"paper_id": "lit_b", "pmid": "11", "title": "Holin paper B", "abstract": ab},
    ])
    cand = _run_extract(out, tmp_path)
    flags = "; ".join(cand["qc_flags"].tolist())
    assert "duplicate_across_papers" in flags


def test_export_warns_when_manual_verification_bypassed(tmp_path: Path, caplog):
    out = tmp_path / "results" / "literature"
    out.mkdir(parents=True)
    pd.DataFrame([{
        "candidate_id": "c1", "title": "t", "pmid": "1", "doi": "",
        "phage_name": "", "gene_name": "S", "protein_name": "S", "host_name": "",
        "evidence_sentence": "s", "evidence_type": "delayed lysis", "evidence_score": 5,
        "confidence_class": "experimentally_validated_holin", "accession": "",
        "sequence": "", "qc_flags": "", "manually_verified": False,
        "curator_decision": "", "curator_notes": "",
    }]).to_csv(out / "manual_review_template.tsv", sep="\t", index=False)
    cfg = Config({"literature": {"output_dir": str(out),
                                 "require_manual_verification_for_gold": False,
                                 "min_evidence_score_for_gold": 4}},
                 tmp_path / "config" / "config.yaml")
    import logging
    with caplog.at_level(logging.WARNING, logger="holinbench"):
        evidence.run_export_curated(cfg)
    msg = " ".join(r.getMessage() for r in caplog.records)
    assert "BYPASSED" in msg or "bypass" in msg.lower()
