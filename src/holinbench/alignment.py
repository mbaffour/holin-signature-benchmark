"""Stage 3 — multiple sequence alignment + alignment-quality metrics.

Builds MSAs for (A) all gold positives, (B) topology subsets (by TMD count), and
(C) family/cluster subsets. Crucially, it QUANTIFIES alignment quality (mean
pairwise identity, % gaps, conserved columns, per-column entropy) rather than
hiding a poor universal alignment — a poor universal MSA is itself a key result.

Backends: pyfamsa (preferred, pure-pip) -> MAFFT/Clustal Omega CLI if present ->
import pre-computed MSAs from a directory. HMMs can still be built from whatever
alignment is produced; quality caveats propagate to the report.
"""
from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from . import clustering, features as features_mod, validate
from .config import Config
from .utils import log, read_fasta, write_fasta

GAP = "-"


# ----------------------------------------------------------- backends ----------
def _align_pyfamsa(records: list[tuple[str, str]]) -> list[tuple[str, str]] | None:
    try:
        from pyfamsa import Aligner, Sequence
    except Exception:
        return None
    seqs = [Sequence(pid.encode(), seq.encode()) for pid, seq in records]
    aligner = Aligner(guide_tree="upgma")
    msa = aligner.align(seqs)
    out = []
    for s in msa:
        sid = s.id.decode() if isinstance(s.id, bytes) else str(s.id)
        seq = s.sequence.decode() if isinstance(s.sequence, bytes) else str(s.sequence)
        out.append((sid, seq))
    return out


def _align_cli(records: list[tuple[str, str]], tool: str) -> list[tuple[str, str]] | None:
    exe = shutil.which("mafft" if tool == "mafft" else "clustalo")
    if exe is None:
        return None
    with tempfile.TemporaryDirectory() as td:
        inp = Path(td) / "in.faa"
        write_fasta(inp, records)
        if tool == "mafft":
            cmd = [exe, "--auto", "--quiet", str(inp)]
            res = subprocess.run(cmd, capture_output=True, text=True)
            out = Path(td) / "out.afa"
            out.write_text(res.stdout)
        else:
            out = Path(td) / "out.afa"
            cmd = [exe, "-i", str(inp), "-o", str(out), "--force", "--outfmt=fasta"]
            subprocess.run(cmd, capture_output=True, text=True)
        return read_fasta(out)


def align_group(records: list[tuple[str, str]], cfg: Config) -> list[tuple[str, str]] | None:
    """Align a group of records; returns aligned (id, gapped_seq) or None."""
    if len(records) < 2:
        return None
    backend = cfg.dotted("alignment.backend", "auto")
    if backend in ("auto", "pyfamsa"):
        aln = _align_pyfamsa(records)
        if aln:
            return aln
    if backend in ("auto", "mafft"):
        aln = _align_cli(records, "mafft")
        if aln:
            return aln
    if backend in ("auto", "clustalo"):
        aln = _align_cli(records, "clustalo")
        if aln:
            return aln
    log.warning("No alignment backend available for a group of %d sequences.", len(records))
    return None


# ----------------------------------------------------------- metrics -----------
def alignment_metrics(aligned: list[tuple[str, str]]) -> dict:
    if not aligned:
        return {}
    seqs = [s for _, s in aligned]
    n = len(seqs)
    L = max(len(s) for s in seqs)
    seqs = [s.ljust(L, GAP) for s in seqs]

    # % gaps
    total = n * L
    gaps = sum(s.count(GAP) for s in seqs)
    pct_gaps = round(100 * gaps / total, 2) if total else 0.0

    # per-column entropy + conserved columns
    entropies, conserved = [], 0
    for col in range(L):
        column = [s[col] for s in seqs if s[col] != GAP]
        if not column:
            continue
        counts: dict[str, int] = {}
        for c in column:
            counts[c] = counts.get(c, 0) + 1
        m = len(column)
        ent = -sum((v / m) * math.log2(v / m) for v in counts.values())
        entropies.append(ent)
        top = max(counts.values()) / m
        if top >= 0.9 and m >= max(2, n // 2):
            conserved += 1
    mean_entropy = round(sum(entropies) / len(entropies), 3) if entropies else 0.0

    # mean pairwise identity over non-gap aligned columns
    idents = []
    for i in range(n):
        for j in range(i + 1, n):
            match = aligned_cols = 0
            for a, b in zip(seqs[i], seqs[j]):
                if a == GAP and b == GAP:
                    continue
                aligned_cols += 1
                if a == b and a != GAP:
                    match += 1
            if aligned_cols:
                idents.append(match / aligned_cols)
    mpi = round(100 * sum(idents) / len(idents), 2) if idents else 0.0

    return {"n_sequences": n, "alignment_length": L, "percent_gaps": pct_gaps,
            "mean_pairwise_identity": mpi, "conserved_columns": conserved,
            "mean_column_entropy": mean_entropy}


# ----------------------------------------------------------- driver ------------
def _topology_subsets(clean: pd.DataFrame, feats: pd.DataFrame) -> dict[str, list]:
    gold_ids = set(clean[clean["dataset_category"] == "gold"]["protein_id"])
    seqmap = dict(zip(clean["protein_id"], clean["sequence"]))
    subsets: dict[str, list] = {}
    for _, r in feats.iterrows():
        if r["protein_id"] not in gold_ids:
            continue
        tmd = int(r.get("tmd_count", -1))
        key = f"{tmd}TMD" if tmd >= 0 else "unknownTMD"
        subsets.setdefault(key, []).append((r["protein_id"], seqmap[r["protein_id"]]))
    return subsets


def run(cfg: Config, clean: pd.DataFrame | None = None,
        feats: pd.DataFrame | None = None,
        clusters: pd.DataFrame | None = None) -> pd.DataFrame:
    if clean is None:
        clean = validate.run(cfg)["clean"]
    if feats is None:
        feats = features_mod.run(cfg, clean)
    if clusters is None:
        clusters = clustering.run(cfg, clean)

    aln_dir = cfg.out_dir("alignments")
    min_seqs = int(cfg.dotted("alignment.min_seqs_per_msa", 3))
    gold = clean[clean["dataset_category"] == "gold"]
    seqmap = dict(zip(clean["protein_id"], clean["sequence"]))

    groups: list[tuple[str, str, list]] = []  # (group_id, group_type, records)

    # (A) universal
    groups.append(("universal", "universal",
                   list(zip(gold["protein_id"], gold["sequence"]))))

    # (B) topology subsets
    for key, recs in _topology_subsets(clean, feats).items():
        groups.append((f"topology_{key}", "topology", recs))

    # (C) family/cluster subsets at the family threshold
    fam_threshold = int(cfg.dotted("hmmer.family_threshold", 0.30) * 100)
    fam_col = f"cluster_{fam_threshold}"
    if not clusters.empty and fam_col in clusters.columns:
        for cid, sub in clusters.groupby(fam_col):
            recs = [(pid, seqmap[pid]) for pid in sub["protein_id"] if pid in seqmap]
            if len(recs) >= min_seqs:
                groups.append((f"family_cluster{int(cid)}", "family", recs))

    rows = []
    for gid, gtype, recs in groups:
        if len(recs) < min_seqs and gtype != "universal":
            continue
        aligned = align_group(recs, cfg) if len(recs) >= 2 else None
        if aligned:
            write_fasta(aln_dir / f"{gid}.afa", aligned)
            m = alignment_metrics(aligned)
        else:
            m = {"n_sequences": len(recs), "alignment_length": 0, "percent_gaps": 0.0,
                 "mean_pairwise_identity": 0.0, "conserved_columns": 0,
                 "mean_column_entropy": 0.0}
        m.update({"group_id": gid, "group_type": gtype,
                  "aligned": bool(aligned)})
        rows.append(m)

    report = pd.DataFrame(rows, columns=["group_id", "group_type", "n_sequences",
                                         "alignment_length", "mean_pairwise_identity",
                                         "percent_gaps", "conserved_columns",
                                         "mean_column_entropy", "aligned"])
    report.to_csv(cfg.out("tables", "alignment_quality.tsv"), sep="\t", index=False)

    uni = report[report["group_id"] == "universal"]
    if not uni.empty and bool(uni.iloc[0]["aligned"]):
        mpi = uni.iloc[0]["mean_pairwise_identity"]
        log.info("Stage 3: universal MSA mean pairwise identity = %.1f%% "
                 "(%s).", mpi,
                 "LOW — universal signal likely weak" if mpi < 25 else "moderate/high")
    return report
