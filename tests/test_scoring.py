import pandas as pd

from holinbench import context, features, scoring, validate


def test_architecture_score_prefers_holin_like(toy_cfg):
    holin = {"length": 80, "tmd_count": 2, "hydrophobic_fraction": 0.55,
             "has_enzymatic_or_structural_domain": False}
    enzyme = {"length": 250, "tmd_count": 0, "hydrophobic_fraction": 0.30,
              "has_enzymatic_or_structural_domain": True}
    s_holin, _ = scoring.architecture_score(holin, toy_cfg)
    s_enzyme, _ = scoring.architecture_score(enzyme, toy_cfg)
    assert s_holin > s_enzyme
    assert 0.0 <= s_enzyme <= 1.0


def test_hmm_score_family_beats_universal(toy_cfg):
    fam = pd.DataFrame([{"protein_id": "p", "model_id": "famX", "model_type": "family",
                         "bitscore": 50.0, "evalue": 1e-9}])
    uni = pd.DataFrame([{"protein_id": "p", "model_id": "universal", "model_type": "universal",
                         "bitscore": 50.0, "evalue": 1e-9}])
    s_fam, info_fam, _ = scoring.hmm_score(fam, toy_cfg)
    s_uni, _, _ = scoring.hmm_score(uni, toy_cfg)
    assert s_fam > s_uni
    assert info_fam["family_hmm_hit"] is True


def test_hmm_score_empty():
    import holinbench.scoring as sc
    s, info, _ = sc.hmm_score(pd.DataFrame(), _cfg_stub())
    assert s == 0.0
    assert info["best_hmm_model"] == ""


def _cfg_stub():
    from holinbench.config import Config
    from pathlib import Path
    return Config({"scoring": {"hmm": {"family_hit_strong": 1.0, "family_hit_weak": 0.6,
                                       "topology_hit": 0.5, "universal_hit": 0.4,
                                       "strong_bitscore": 30, "weak_bitscore": 12}}},
                  Path("config/config.yaml"))


def test_context_norm_monotonic():
    assert scoring.context_norm(-4) < scoring.context_norm(0) < scoring.context_norm(4)
    assert abs(scoring.context_norm(0) - 0.5) < 1e-9


def test_confidence_categories(toy_cfg):
    assert scoring.confidence_category(0.9, toy_cfg) == "high_confidence_candidate"
    assert scoring.confidence_category(0.55, toy_cfg) == "medium_confidence_candidate"
    assert scoring.confidence_category(0.35, toy_cfg) == "weak_candidate"
    assert scoring.confidence_category(0.1, toy_cfg) == "unlikely_holin"


def test_full_scoring_ranks_cassette_candidate_above_soluble(toy_cfg):
    clean = validate.run(toy_cfg)["clean"]
    feats = features.run(toy_cfg, clean)
    ctx = context.run(toy_cfg, feats)
    ranking = scoring.run(toy_cfg, feats, ctx, hmm_hits=None,
                          restrict_categories=["unknown"])
    assert list(ranking["protein_id"])  # non-empty
    score_c1 = ranking[ranking["protein_id"] == "c1"]["final_holin_score"].iloc[0]
    score_c2 = ranking[ranking["protein_id"] == "c2"]["final_holin_score"].iloc[0]
    assert score_c1 > score_c2
    assert ranking.iloc[0]["explanation"]  # explanation present
