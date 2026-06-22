from pathlib import Path

import pandas as pd

from holinbench import evidence
from holinbench.config import Config


def test_evidence_score_genetic_plus_functional_is_5():
    s = evidence._evidence_score(["amber mutant"], ["delayed lysis"], [], [], [])
    assert s == 5


def test_evidence_score_membrane_only_is_4():
    s = evidence._evidence_score([], [], ["membrane depolarization"], [], [])
    assert s == 4


def test_evidence_score_functional_only_is_3():
    s = evidence._evidence_score([], ["required for lysis"], [], [], [])
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
