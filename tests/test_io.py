from holinbench import io
from holinbench.utils import read_fasta, write_fasta


def test_load_all_datasets_categories(toy_cfg):
    df = io.load_all_datasets(toy_cfg)
    cats = set(df["dataset_category"])
    assert cats == {"gold", "weak", "hard_negative", "unknown"}
    # gold has 2, weak 1, neg 2, unknown 2
    counts = df["dataset_category"].value_counts().to_dict()
    assert counts["gold"] == 2
    assert counts["unknown"] == 2


def test_fasta_roundtrip(tmp_path):
    recs = [("a", "MKLV"), ("b", "MMMWWW")]
    p = tmp_path / "x.faa"
    write_fasta(p, recs)
    assert read_fasta(p) == recs


def test_load_lysis_terms(toy_cfg):
    terms = io.load_lysis_terms(toy_cfg)
    assert "endolysin" in terms
    assert "endolysin" in terms["endolysin"]


def test_load_context_numeric(toy_cfg):
    ctx = io.load_context(toy_cfg)
    assert not ctx.empty
    assert ctx["start"].dtype.kind in "if"  # coerced to numeric
