from holinbench import clustering


def test_greedy_cluster_groups_identical_separates_different():
    holin = "MKR" + "LLIVFAMWLLIVFAMWLLIV" + "KKDE"
    same = holin
    diff = "MKR" + "KRDESTNQKRDESTNQKRDE" + "STNQ"  # polar, dissimilar
    assign = clustering.greedy_cluster(["a", "b", "c"], [holin, same, diff], 0.7)
    assert assign["a"] == assign["b"]      # identical -> same cluster
    assert assign["c"] != assign["a"]      # dissimilar -> different cluster


def test_greedy_cluster_threshold_monotonic():
    seqs = ["MKLLIVFAMWLLIVFAMWGG",
            "MKLLIVFAMWLLIVFAMWGA",   # 1 diff
            "MKLLIVFAGGLLIVFAMWGG",   # few diffs
            "MQRPSTNQKRDESTNQGHKR"]   # very different
    ids = ["s1", "s2", "s3", "s4"]
    n_strict = len(set(clustering.greedy_cluster(ids, seqs, 0.95).values()))
    n_loose = len(set(clustering.greedy_cluster(ids, seqs, 0.40).values()))
    # higher identity threshold -> at least as many clusters
    assert n_strict >= n_loose


def test_run_writes_assignments_and_fragmentation(toy_cfg):
    clusters = clustering.run(toy_cfg)
    assert "cluster_30" in clusters.columns
    assert (toy_cfg.out_dir("tables") / "cluster_summary.tsv").exists()
