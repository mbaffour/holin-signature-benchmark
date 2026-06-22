from holinbench import context, features, validate


def test_context_scores_lysis_cassette_high(toy_cfg):
    clean = validate.run(toy_cfg)["clean"]
    feats = features.run(toy_cfg, clean)
    ctx = context.run(toy_cfg, feats)
    c1 = ctx[ctx["protein_id"] == "c1"].iloc[0]
    # c1 is adjacent to an endolysin and a spanin -> positive context score
    assert bool(c1["near_endolysin"]) is True
    assert bool(c1["near_spanin"]) is True
    assert c1["context_score"] > 0


def test_context_scores_isolated_candidate_low(toy_cfg):
    clean = validate.run(toy_cfg)["clean"]
    feats = features.run(toy_cfg, clean)
    ctx = context.run(toy_cfg, feats)
    c2 = ctx[ctx["protein_id"] == "c2"].iloc[0]
    # c2 is on a contig with no lysis genes nearby -> isolated, non-positive
    assert bool(c2["isolated_from_lysis_genes"]) is True
    assert c2["context_score"] <= 0


def test_context_endolysin_distance(toy_cfg):
    clean = validate.run(toy_cfg)["clean"]
    feats = features.run(toy_cfg, clean)
    ctx = context.run(toy_cfg, feats)
    c1 = ctx[ctx["protein_id"] == "c1"].iloc[0]
    assert int(c1["nearest_endolysin_genes"]) == 1
