"""Shared utilities: logging, FASTA/TSV IO, and sequence biochemistry.

Sequence helpers here are intentionally dependency-light (pure numpy/stdlib) so
they are easy to unit-test and reason about. They are used by validate.py and
features.py.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Iterable, Iterator

# -------------------------------------------------------------------- logging --
def setup_logging(level: int = logging.INFO) -> logging.Logger:
    log = logging.getLogger("holinbench")
    if not log.handlers:
        h = logging.StreamHandler(sys.stderr)
        h.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        log.addHandler(h)
    log.setLevel(level)
    return log


log = logging.getLogger("holinbench")


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


# --------------------------------------------------------------------- FASTA ---
def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read a FASTA file into a list of (id, sequence). Robust to blank lines."""
    records: list[tuple[str, str]] = []
    pid: str | None = None
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if pid is not None:
                    records.append((pid, "".join(chunks)))
                pid = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if pid is not None:
        records.append((pid, "".join(chunks)))
    return records


def write_fasta(path: str | Path, records: Iterable[tuple[str, str]], width: int = 60) -> None:
    ensure_dir(Path(path).parent)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for pid, seq in records:
            fh.write(f">{pid}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    yield from read_fasta(path)


# ----------------------------------------------------------- biochemistry ------
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = set("XBZJUO*")

# Residues counted as "hydrophobic" for hydrophobic-fraction (Kyte-Doolittle
# positive / membrane-favoring set). Documented choice; tweak with care.
HYDROPHOBIC_AA = set("AILMFWVC")

POSITIVE_AA = set("KR")
NEGATIVE_AA = set("DE")

# Kyte & Doolittle (1982) hydropathy index.
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Average residue masses (monomer, Da), water added once for MW.
_AA_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
_WATER = 18.01524


def clean_sequence(seq: str) -> str:
    """Uppercase, strip whitespace and gaps/stops; keep letters only."""
    return "".join(c for c in seq.upper() if c.isalpha())


def is_valid_sequence(seq: str, max_ambiguous_fraction: float = 0.10) -> tuple[bool, str]:
    """Return (ok, reason). Empty / non-protein / too-ambiguous sequences fail."""
    if not seq:
        return False, "empty"
    bad = [c for c in seq if c not in STANDARD_AA and c not in AMBIGUOUS_AA]
    if bad:
        return False, f"non_amino_acid_chars:{''.join(sorted(set(bad)))[:8]}"
    amb = sum(1 for c in seq if c in AMBIGUOUS_AA)
    if len(seq) and amb / len(seq) > max_ambiguous_fraction:
        return False, f"too_many_ambiguous:{amb}/{len(seq)}"
    return True, "ok"


def hydrophobic_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(1 for c in seq if c in HYDROPHOBIC_AA) / len(seq)


def aa_composition(seq: str) -> dict[str, float]:
    comp = {aa: 0 for aa in sorted(STANDARD_AA)}
    for c in seq:
        if c in comp:
            comp[c] += 1
    n = len(seq) or 1
    return {aa: comp[aa] / n for aa in comp}


def charge_counts(seq: str) -> tuple[int, int]:
    """Return (positive, negative) residue counts."""
    pos = sum(1 for c in seq if c in POSITIVE_AA)
    neg = sum(1 for c in seq if c in NEGATIVE_AA)
    return pos, neg


def net_charge(seq: str) -> int:
    pos, neg = charge_counts(seq)
    return pos - neg


def molecular_weight(seq: str) -> float:
    """Approximate average molecular weight in Da."""
    if not seq:
        return 0.0
    total = sum(_AA_MASS.get(c, 0.0) for c in seq if c in _AA_MASS)
    return round(total + _WATER, 2)


def kd_hydropathy_profile(seq: str, window: int = 19) -> list[float]:
    """Sliding-window mean Kyte-Doolittle hydropathy (centered)."""
    if not seq:
        return []
    vals = [KD.get(c, 0.0) for c in seq]
    if len(vals) < window:
        return [sum(vals) / len(vals)] * len(vals)
    half = window // 2
    out: list[float] = []
    for i in range(len(vals)):
        lo, hi = max(0, i - half), min(len(vals), i + half + 1)
        out.append(sum(vals[lo:hi]) / (hi - lo))
    return out


def estimate_tmd_count(seq: str, window: int = 19, threshold: float = 1.6,
                       min_gap: int = 5) -> int:
    """Crude TMD count from KD hydropathy peaks (fallback only).

    This is NOT a substitute for DeepTMHMM/Phobius — it is a rough estimate used
    only when no imported topology is available, and is always flagged as such.
    """
    prof = kd_hydropathy_profile(seq, window=window)
    if not prof:
        return 0
    above = [v >= threshold for v in prof]
    count, i, n = 0, 0, len(above)
    while i < n:
        if above[i]:
            count += 1
            while i < n and above[i]:
                i += 1
            # skip a short gap so adjacent peaks aren't double-counted
            gap = 0
            while i < n and not above[i] and gap < min_gap:
                i += 1
                gap += 1
        else:
            i += 1
    return count
