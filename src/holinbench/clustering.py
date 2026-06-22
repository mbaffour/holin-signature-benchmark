"""Stage 2 — redundancy reduction and clustering of gold positives.

Clusters the gold-standard holins at several identity thresholds to answer a
central question: do experimentally characterized holins form ONE connected
sequence space, or many isolated islands? Many islands => a single universal
HMM is unlikely to work, which is exactly what this project is testing.

Backends (auto-selected): mmseqs2 / cd-hit CLI if present, else a dependency-light
built-in greedy clusterer using pairwise alignment identity (Biopython). The
built-in path is deterministic and fine for the small gold sets these datasets
contain; it warns if asked to cluster a very large set.
"""
from __future__ import annotations

import shutil

import pandas as pd

from . import io, validate
from .config import Config
from .utils import log

_BUILTIN_MAX = 800  # all-pairs identity above this is slow; warn the user.


# --------------------------------------------------------- identity ------------
def _pairwise_identity(a: str, b: str, aligner) -> tuple[float, float]:
    """Return (identity, coverage).

    identity = identities / length of shorter sequence.
    coverage = aligned residues (match+mismatch) / length of LONGER sequence.

    The coverage term matters because a short sequence that is a prefix/subset of
    a longer one can score identity 1.0 while covering only a fraction of the
    longer (representative) sequence; without a coverage gate such pairs falsely
    merge and inflate cluster connectivity (the project's headline metric).
    """
    if not a or not b:
        return 0.0, 0.0
    if a == b:
        return 1.0, 1.0
    try:
        counts = aligner.align(a, b)[0].counts()
        ident = counts.identities
        aligned = counts.identities + counts.mismatches
    except Exception:
        # Fallback: k-mer Jaccard as a rough proxy if alignment fails.
        ji = _kmer_identity(a, b)
        return ji, ji
    identity = ident / min(len(a), len(b))
    coverage = aligned / max(len(a), len(b))
    return identity, coverage


def _kmer_identity(a: str, b: str, k: int = 3) -> float:
    ka = {a[i:i + k] for i in range(len(a) - k + 1)}
    kb = {b[i:i + k] for i in range(len(b) - k + 1)}
    if not ka or not kb:
        return 0.0
    return len(ka & kb) / len(ka | kb)


def _make_aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices
    al = PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10
    al.extend_gap_score = -0.5
    return al


def greedy_cluster(ids: list[str], seqs: list[str], threshold: float,
                   coverage: float = 0.0) -> dict[str, int]:
    """Deterministic greedy clustering: longest-first representatives.

    A sequence joins a cluster only if it meets BOTH the identity `threshold`
    and the `coverage` fraction against the cluster representative.
    """
    order = sorted(range(len(ids)), key=lambda i: (-len(seqs[i]), ids[i]))
    aligner = _make_aligner()
    reps: list[int] = []
    assignment: dict[str, int] = {}
    for i in order:
        placed = False
        for cid, rep in enumerate(reps):
            identity, cov = _pairwise_identity(seqs[i], seqs[rep], aligner)
            if identity >= threshold and cov >= coverage:
                assignment[ids[i]] = cid
                placed = True
                break
        if not placed:
            assignment[ids[i]] = len(reps)
            reps.append(i)
    return assignment


# --------------------------------------------------------- backends ------------
def _cli_available(name: str) -> bool:
    return shutil.which(name) is not None


def cluster_records(records: list[tuple[str, str]], threshold: float,
                    backend: str = "auto", cfg: Config | None = None) -> dict[str, int]:
    ids = [r[0] for r in records]
    seqs = [r[1] for r in records]
    if not ids:
        return {}
    chosen = backend
    if backend == "auto":
        if _cli_available("mmseqs"):
            chosen = "mmseqs2"
        elif _cli_available("cd-hit"):
            chosen = "cdhit"
        else:
            chosen = "builtin"
    if chosen in ("mmseqs2", "cdhit"):
        # Native CLI clustering is available on Linux/macOS native envs. For the
        # cross-platform default we always have the builtin path; we keep the CLI
        # branch documented but delegate to builtin unless explicitly wired, to
        # avoid silent dependence on tool-specific output parsing.
        log.info("Native clusterer '%s' detected; using builtin greedy clusterer "
                 "for portability/determinism.", chosen)
        chosen = "builtin"
    if chosen == "builtin" and len(ids) > _BUILTIN_MAX:
        log.warning("Builtin clusterer on %d sequences may be slow (all-pairs).",
                    len(ids))
    coverage = float(cfg.dotted("clustering.coverage", 0.0)) if cfg is not None else 0.0
    return greedy_cluster(ids, seqs, threshold, coverage=coverage)


# --------------------------------------------------------- driver --------------
def run(cfg: Config, clean: pd.DataFrame | None = None) -> pd.DataFrame:
    if clean is None:
        clean = validate.run(cfg)["clean"]
    gold = clean[clean["dataset_category"] == "gold"]
    records = list(zip(gold["protein_id"], gold["sequence"]))
    ccfg = cfg.section("clustering")
    thresholds = ccfg.get("identity_thresholds", [0.90, 0.70, 0.50, 0.30])
    backend = ccfg.get("backend", "auto")

    if not records:
        log.warning("No gold positives to cluster.")
        empty = pd.DataFrame(columns=["protein_id"])
        empty.to_csv(cfg.out("tables", "cluster_assignments.tsv"), sep="\t", index=False)
        return empty

    assign = {pid: {} for pid, _ in records}
    summary_rows = []
    for thr in thresholds:
        clusters = cluster_records(records, thr, backend=backend, cfg=cfg)
        for pid, cid in clusters.items():
            assign[pid][f"cluster_{int(thr*100)}"] = cid
        n_clusters = len(set(clusters.values()))
        sizes = pd.Series(list(clusters.values())).value_counts()
        n_singletons = int((sizes == 1).sum())
        summary_rows.append({
            "threshold": thr, "n_sequences": len(records), "n_clusters": n_clusters,
            "n_singletons": n_singletons, "largest_cluster": int(sizes.max()),
            "mean_cluster_size": round(len(records) / n_clusters, 2),
        })

    rows = []
    for pid, _ in records:
        row = {"protein_id": pid}
        row.update(assign[pid])
        rows.append(row)
    assignments = pd.DataFrame(rows)
    assignments = assignments.merge(
        gold[["protein_id", "family_label", "holin_type"]], on="protein_id", how="left")
    assignments.to_csv(cfg.out("tables", "cluster_assignments.tsv"), sep="\t", index=False)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(cfg.out("tables", "cluster_summary.tsv"), sep="\t", index=False)

    # Headline interpretation at the primary threshold.
    primary = ccfg.get("primary_threshold", 0.30)
    pcol = f"cluster_{int(primary*100)}"
    if pcol in assignments.columns:
        n_fam = assignments[pcol].nunique()
        log.info("Stage 2: %d gold holins form %d cluster(s) at %.0f%% identity "
                 "(%s sequence space).", len(records), n_fam, primary * 100,
                 "connected" if n_fam == 1 else "fragmented")
    return assignments
