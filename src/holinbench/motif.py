"""Stage 10 — motif / conservation analysis (conservative by construction).

For each alignment we quantify per-column information content and identify
"conserved" columns. Crucially, we then ask whether that conservation is just
generic transmembrane HYDROPHOBICITY (architecture-level) rather than a genuine
holin-specific linear motif. A column whose conserved residue is hydrophobic and
whose neighbours are also hydrophobic is reported as architecture-level, NOT as a
universal motif. Sequence logos are produced where logomaker is available.
"""
from __future__ import annotations

import math

import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import pandas as pd

from .config import Config
from .utils import HYDROPHOBIC_AA, log, read_fasta

GAP = "-"
MAX_BITS = math.log2(20)


def _column_stats(seqs: list[str], col: int) -> dict:
    column = [s[col] for s in seqs if col < len(s) and s[col] != GAP]
    n = len(column)
    if n == 0:
        return {"n": 0, "info_bits": 0.0, "top_residue": "-", "top_freq": 0.0}
    counts: dict[str, int] = {}
    for c in column:
        counts[c] = counts.get(c, 0) + 1
    ent = -sum((v / n) * math.log2(v / n) for v in counts.values())
    info = MAX_BITS - ent
    top = max(counts, key=counts.get)
    return {"n": n, "info_bits": round(info, 3), "top_residue": top,
            "top_freq": round(counts[top] / n, 3)}


def _logo(seqs: list[str], out_png) -> bool:
    try:
        import logomaker
    except Exception:
        return False
    L = max(len(s) for s in seqs)
    seqs = [s.ljust(L, GAP) for s in seqs]
    rows = []
    for col in range(L):
        counts = {aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWY"}
        tot = 0
        for s in seqs:
            c = s[col]
            if c in counts:
                counts[c] += 1
                tot += 1
        rows.append({aa: (counts[aa] / tot if tot else 0.0) for aa in counts})
    mat = pd.DataFrame(rows)
    try:
        info_mat = logomaker.transform_matrix(mat, from_type="probability",
                                              to_type="information")
        fig, ax = plt.subplots(figsize=(max(4, L * 0.18), 2.2))
        logomaker.Logo(info_mat, ax=ax)
        ax.set_ylabel("bits")
        ax.set_xlabel("alignment column")
        fig.tight_layout()
        fig.savefig(out_png, dpi=130)
        plt.close(fig)
        return True
    except Exception as exc:  # pragma: no cover - defensive
        log.info("Logo generation failed for %s: %s", out_png.name, exc)
        return False


def run(cfg: Config) -> pd.DataFrame:
    aln_dir = cfg.out_dir("alignments")
    fig_dir = cfg.out_dir("figures")
    mcfg = cfg.section("motif")
    cons_bits = float(mcfg.get("entropy_conserved_bits", 1.0))
    make_logos = bool(mcfg.get("make_logos", True))

    rows = []
    for afa in sorted(aln_dir.glob("*.afa")):
        gid = afa.stem
        aligned = read_fasta(afa)
        if len(aligned) < 2:
            continue
        seqs = [s for _, s in aligned]
        L = max(len(s) for s in seqs)
        conserved_cols = []
        hydrophobic_conserved = 0
        for col in range(L):
            st = _column_stats(seqs, col)
            if st["n"] >= max(2, len(seqs) // 2) and st["info_bits"] >= cons_bits:
                conserved_cols.append((col, st["top_residue"], st["info_bits"]))
                if st["top_residue"] in HYDROPHOBIC_AA:
                    hydrophobic_conserved += 1
        n_cons = len(conserved_cols)
        hydro_frac = (hydrophobic_conserved / n_cons) if n_cons else 0.0
        interpretation = "no conserved columns"
        if n_cons:
            if hydro_frac >= 0.7:
                interpretation = ("conservation dominated by hydrophobic residues "
                                  "= architecture-level (TM) signal, NOT a specific motif")
            elif hydro_frac >= 0.4:
                interpretation = "mixed: partly hydrophobic, partly specific residues"
            else:
                interpretation = "conserved non-hydrophobic residues present (candidate motif)"
        logo_made = False
        if make_logos:
            logo_made = _logo(seqs, fig_dir / f"logo_{gid}.png")
        rows.append({
            "group_id": gid, "n_sequences": len(aligned), "alignment_length": L,
            "n_conserved_columns": n_cons,
            "hydrophobic_fraction_of_conserved": round(hydro_frac, 3),
            "interpretation": interpretation, "logo_png": f"logo_{gid}.png" if logo_made else "",
        })

    summary = pd.DataFrame(rows, columns=["group_id", "n_sequences", "alignment_length",
                                          "n_conserved_columns",
                                          "hydrophobic_fraction_of_conserved",
                                          "interpretation", "logo_png"])
    summary.to_csv(cfg.out("tables", "motif_conservation_summary.tsv"), sep="\t", index=False)
    uni = summary[summary["group_id"] == "universal"]
    if not uni.empty:
        log.info("Stage 10: universal alignment — %d conserved column(s); %s",
                 int(uni.iloc[0]["n_conserved_columns"]), uni.iloc[0]["interpretation"])
    return summary
