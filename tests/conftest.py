"""Shared pytest fixtures: a hermetic toy config + dataset in a temp dir."""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from holinbench.config import Config

# A small holin-like synthetic sequence (2 hydrophobic stretches).
HOLIN_LIKE = ("MKRDS" + "LLIVFAMWLLIVFAMWLLIV" + "GSTNKD" +
              "FAMWLLIVFAMWLLIVFAMW" + "KKDERSTQ")
SOLUBLE = "M" + "KRDESTNQGHKRDESTNQGHKRDESTNQGH" * 4  # polar, no TMD


def _write(p: Path, text: str) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(textwrap.dedent(text).lstrip("\n"), encoding="utf-8")


@pytest.fixture
def toy_root(tmp_path: Path) -> Path:
    d = tmp_path / "data" / "example"
    # gold: 2 records, second is missing citation + evidence_type (QC flags)
    _write(d / "gold_holins.csv", f"""
        protein_id,protein_name,phage_name,host,accession,sequence,evidence_type,citation,holin_type,family_label,notes
        g1,S holin,Lambda,Ecoli,ACC1,{HOLIN_LIKE},amber_mutant,Smith1990,class_I,lambdaS,note
        g2,holin,PhageX,Ecoli,ACC2,{HOLIN_LIKE}A,,,unknown,,missing-evidence
        """)
    # weak
    _write(d / "weak_annotated_holins.csv", f"""
        protein_id,protein_name,phage_name,accession,sequence,annotation,source_database,notes
        w1,putative holin,PhageW,ACCW,{HOLIN_LIKE}KK,putative holin,RefSeq,note
        """)
    # hard negatives: one duplicates g1's sequence (cross-class dup flag);
    # one is a soluble amidase (enzymatic domain).
    _write(d / "hard_negatives.csv", f"""
        protein_id,protein_name,source,accession,sequence,negative_type,annotation,notes
        n1,dup of gold,synthetic,ACCN1,{HOLIN_LIKE},qc_duplicate,cross-class dup,note
        n2,amidase,synthetic,ACCN2,{SOLUBLE},enzymatic_domain,N-acetylmuramoyl-L-alanine amidase,note
        """)
    # unknown candidates
    _write(d / "proteins.faa", f">c1\n{HOLIN_LIKE}\n>c2\n{SOLUBLE}\n")
    # topology import
    _write(d / "topology_predictions.tsv", f"""
        protein_id\tlength\ttmd_count\ttopology\tn_region\tc_region\tsignal_peptide\tsar_like\ttool
        g1\t{len(HOLIN_LIKE)}\t2\tin\t1-30\t40-70\tno\tno\ttest
        g2\t{len(HOLIN_LIKE)+1}\t2\tin\t1-30\t40-70\tno\tno\ttest
        w1\t{len(HOLIN_LIKE)+2}\t2\tin\t1-30\t40-70\tno\tno\ttest
        n1\t{len(HOLIN_LIKE)}\t2\tin\t1-30\t40-70\tno\tno\ttest
        n2\t{len(SOLUBLE)}\t0\tin\t1-30\t40-70\tno\tno\ttest
        c1\t{len(HOLIN_LIKE)}\t2\tin\t1-30\t40-70\tno\tyes\ttest
        c2\t{len(SOLUBLE)}\t0\tin\t1-30\t40-70\tno\tno\ttest
        """)
    # genomic context: c1 sits next to an endolysin and a spanin (lysis cassette).
    _write(d / "context.tsv", f"""
        contig_id\tgene_id\tstart\tend\tstrand\tprotein_id\tproduct\tsequence
        ctgA\tA1\t100\t400\t+\tc1\thypothetical protein\t{HOLIN_LIKE}
        ctgA\tA2\t420\t900\t+\tendoA\tendolysin\t
        ctgA\tA3\t920\t1200\t+\tspanA\tRz spanin\t
        ctgB\tB1\t100\t700\t+\tc2\thypothetical protein\t{SOLUBLE}
        ctgB\tB2\t720\t1400\t+\ttailB\ttail protein\t
        """)
    _write(d / "lysis_context_terms.txt", """
        endolysin\tendolysin
        endolysin\tlysin
        spanin\tspanin
        spanin\tRz
        holin\tholin
        """)
    return tmp_path


@pytest.fixture
def toy_cfg(toy_root: Path) -> Config:
    data = {
        "project": {"random_seed": 1},
        "paths": {
            "output_dir": "results",
            "gold_holins": "data/example/gold_holins.csv",
            "weak_holins": "data/example/weak_annotated_holins.csv",
            "hard_negatives": "data/example/hard_negatives.csv",
            "proteins_faa": "data/example/proteins.faa",
            "context_tsv": "data/example/context.tsv",
            "topology_tsv": "data/example/topology_predictions.tsv",
            "lysis_terms": "data/example/lysis_context_terms.txt",
        },
        "validation": {"min_length": 30, "max_length": 300, "hard_min_length": 10,
                       "max_ambiguous_fraction": 0.10, "drop_invalid": True,
                       "drop_exact_duplicates": True},
        "features": {"kyte_doolittle_window": 19, "tmd_hydrophobicity_threshold": 1.6,
                     "terminal_window": 30, "sar_max_charge_nterm": 1,
                     "topology_source": "import"},
        "context": {"window_genes": 5, "window_bp": 5000, "endolysin_near_genes": 3,
                    "weights": {"near_endolysin_within_3genes": 2.0,
                                "same_strand_as_endolysin": 1.0,
                                "near_spanin_or_rz": 1.0,
                                "near_sar_with_pinholin_topology": 1.0,
                                "no_better_holin_candidate": 1.0,
                                "has_alternative_annotation": -2.0,
                                "not_near_any_lysis_gene": -1.0,
                                "too_large_or_multidomain": -1.0}},
        "scoring": {
            "weights": {"hmm": 0.45, "architecture": 0.30, "context": 0.25},
            "hmm": {"family_hit_strong": 1.0, "family_hit_weak": 0.6, "topology_hit": 0.5,
                    "universal_hit": 0.4, "strong_bitscore": 30.0, "weak_bitscore": 12.0},
            "architecture": {"preferred_length": [50, 150], "acceptable_length": [30, 200],
                             "preferred_tmd": [1, 4], "min_hydrophobic_fraction": 0.40,
                             "length_in_preferred": 1.0, "length_in_acceptable": 0.5,
                             "length_out": 0.0, "tmd_in_range": 1.0, "tmd_out_of_range": 0.0,
                             "hydrophobic_ok": 1.0,
                             "has_enzymatic_or_structural_domain": -1.0, "no_tmd": -1.0},
            "confidence_bins": {"high_confidence_candidate": 0.70,
                                "medium_confidence_candidate": 0.50,
                                "weak_candidate": 0.30},
        },
    }
    cfg_path = toy_root / "config" / "config.yaml"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text("# toy\n", encoding="utf-8")
    return Config(data, cfg_path)
