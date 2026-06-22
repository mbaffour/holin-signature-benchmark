import pandas as pd

from holinbench import benchmark, context, features, validate


def test_operating_point_perfect():
    op = benchmark._operating_point([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], 0.5)
    assert op["recall"] == 1.0
    assert op["specificity"] == 1.0
    assert op["fp"] == 0 and op["fn"] == 0


def test_safe_auc_separable():
    auc = benchmark._safe_auc([0, 0, 1, 1], [0.0, 0.1, 0.9, 1.0], "roc")
    assert auc == 1.0


def test_safe_auc_single_class_is_nan():
    auc = benchmark._safe_auc([1, 1, 1], [0.2, 0.5, 0.9], "roc")
    assert auc != auc  # NaN


def test_topk_recovery():
    # 2 positives; top-2 by score are both positives
    assert benchmark._topk_recovery([1, 0, 1, 0], [0.9, 0.1, 0.8, 0.2]) == 1.0


def test_best_f1_threshold_runs():
    t = benchmark._best_f1_threshold([0, 1, 1], [0.2, 0.6, 0.9])
    assert isinstance(t, float)


def test_build_score_table_labels(toy_cfg):
    clean = validate.run(toy_cfg)["clean"]
    feats = features.run(toy_cfg, clean)
    ctx = context.run(toy_cfg, feats)
    # synthetic HMM hits: g1 gets a strong family hit
    hits = pd.DataFrame([{"protein_id": "g1", "model_id": "famA", "model_type": "family",
                          "bitscore": 60.0, "evalue": 1e-12}])
    score_df = benchmark._build_score_table(toy_cfg, feats, ctx, hits)
    # only gold + hard_negative rows, labeled 1/0
    assert set(score_df["label"]) <= {0, 1}
    g1 = score_df[score_df["protein_id"] == "g1"]
    assert not g1.empty and g1.iloc[0]["label"] == 1
    assert g1.iloc[0]["family"] == 60.0


def test_evaluate_models_columns(toy_cfg):
    score_df = pd.DataFrame({
        "protein_id": ["g1", "g2", "n1", "n2"],
        "label": [1, 1, 0, 0],
        "universal": [50, 40, 5, 0],
        "architecture": [0.9, 0.8, 0.2, 0.1],
        "context": [0.7, 0.6, 0.4, 0.3],
        "hmm_arch": [0.8, 0.7, 0.2, 0.1],
        "hmm_arch_context": [0.85, 0.72, 0.25, 0.12],
        "topology": [0, 0, 0, 0], "family": [0, 0, 0, 0],
    })
    res = benchmark._evaluate_models(score_df, ["universal", "architecture",
                                                "hmm_arch_context"], toy_cfg)
    assert set(["model", "roc_auc", "pr_auc", "f1", "circularity_warning"]) <= set(res.columns)
    uni = res[res["model"] == "universal"].iloc[0]
    assert uni["circularity_warning"].startswith("YES")
