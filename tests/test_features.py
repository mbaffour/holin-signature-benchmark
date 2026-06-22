from holinbench import features, validate
from holinbench.utils import (estimate_tmd_count, hydrophobic_fraction,
                              molecular_weight, net_charge)


def test_hydrophobic_fraction_bounds():
    assert hydrophobic_fraction("") == 0.0
    assert hydrophobic_fraction("LLLL") == 1.0
    assert 0.0 <= hydrophobic_fraction("MKRDESTNQ") <= 1.0


def test_net_charge():
    assert net_charge("KR") == 2
    assert net_charge("DE") == -2
    assert net_charge("KKDE") == 0


def test_molecular_weight_monotonic():
    assert molecular_weight("M") < molecular_weight("MM") < molecular_weight("MMM")


def test_estimate_tmd_count_detects_hydrophobic_stretch():
    soluble = "KRDESTNQGH" * 6
    membrane = "KR" + "LLIVFAMWLLIVFAMWLLIV" + "KR"
    assert estimate_tmd_count(soluble) == 0
    assert estimate_tmd_count(membrane) >= 1


def test_feature_table_has_expected_columns(toy_cfg):
    clean = validate.run(toy_cfg)["clean"]
    feats = features.run(toy_cfg, clean)
    for col in ["length", "tmd_count", "hydrophobic_fraction",
                "has_enzymatic_or_structural_domain", "nterm_pos_charge"]:
        assert col in feats.columns
    # amidase negative carries an enzymatic-domain annotation
    n2 = feats[feats["protein_id"] == "n2"].iloc[0]
    assert bool(n2["has_enzymatic_or_structural_domain"]) is True
    # imported topology used for g1
    g1 = feats[feats["protein_id"] == "g1"].iloc[0]
    assert int(g1["tmd_count"]) == 2
