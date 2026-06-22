from holinbench import validate


def test_validate_flags_missing_citation(toy_cfg):
    res = validate.run(toy_cfg)
    flags = res["flags"]
    g2_flags = set(flags[flags["protein_id"] == "g2"]["flag"])
    assert "gold_missing_citation" in g2_flags
    assert "gold_missing_evidence_type" in g2_flags


def test_validate_flags_cross_class_duplicate(toy_cfg):
    res = validate.run(toy_cfg)
    flags = res["flags"]
    dup = flags[flags["flag"] == "duplicate_sequence_across_classes"]
    # g1 (gold) and n1 (hard_negative) share a sequence
    ids = set(dup["protein_id"])
    assert {"g1", "n1"} <= ids


def test_validate_computes_metrics(toy_cfg):
    res = validate.run(toy_cfg)
    clean = res["clean"]
    assert "hydrophobic_fraction" in clean.columns
    assert (clean["length"] > 0).all()
    # holin-like records should be quite hydrophobic
    g1 = clean[clean["protein_id"] == "g1"].iloc[0]
    assert g1["hydrophobic_fraction"] > 0.4


def test_validate_writes_outputs(toy_cfg):
    res = validate.run(toy_cfg)
    assert res["metadata_path"].exists()
    assert res["flags_path"].exists()
