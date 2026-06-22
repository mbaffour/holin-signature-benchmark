"""Stage 0c — link candidate proteins to accessions / sequences (best effort).

Scans each candidate's source paper for accession-like identifiers and tries to
fetch the protein sequence from NCBI. Name->accession mapping from free text is
inherently heuristic, so every mapping is tagged with a confidence and a clear
"inferred" flag. Candidates with no resolvable accession keep an empty sequence
and are flagged — they must not enter the gold set automatically.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from . import literature
from .config import Config
from .utils import log

# Accession patterns (rough): GenBank protein, RefSeq protein, UniProt.
PROTEIN_ACC = re.compile(
    r"\b("
    r"[A-Z]{3}\d{5}(?:\.\d+)?"           # GenBank protein e.g. AAA12345
    r"|[A-Z]{2}_\d{6,9}(?:\.\d+)?"        # RefSeq e.g. NP_123456 / YP_009123456
    r"|[OPQ][0-9][A-Z0-9]{3}[0-9]"        # UniProt
    r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"  # UniProt
    r")\b")


def _paper_text(row, paper_id: str, fulltext_dir: Path) -> str:
    parts = [str(row.get("title", "")), str(row.get("abstract", ""))]
    ft = fulltext_dir / f"{paper_id}.txt"
    if ft.exists():
        parts.append(ft.read_text(encoding="utf-8"))
    return " ".join(parts)


def _fetch_protein_fasta(cfg, acc, email, api_key, delay) -> str:
    params = {"db": "protein", "id": acc, "rettype": "fasta", "retmode": "text",
              "email": email}
    if api_key:
        params["api_key"] = api_key
    text = literature._http_get(cfg, f"{literature.EUTILS}/efetch.fcgi", params,
                                is_json=False, delay=delay)
    if not text or not text.startswith(">"):
        return ""
    seq = "".join(line.strip() for line in text.splitlines()[1:] if not line.startswith(">"))
    return "".join(c for c in seq.upper() if c.isalpha())


def run_map(cfg: Config) -> pd.DataFrame:
    lit = cfg.section("literature")
    out_dir = cfg.resolve(lit.get("output_dir", "results/literature"))
    fulltext_dir = out_dir / "fulltext"
    cand_path = out_dir / "candidate_holin_literature_table.tsv"
    search_path = out_dir / "literature_search_results.tsv"
    if not cand_path.exists():
        log.warning("No candidate table; run `extract-evidence` first.")
        return pd.DataFrame()
    cand = pd.read_csv(cand_path, sep="\t", dtype=str, keep_default_na=False)
    papers = pd.read_csv(search_path, sep="\t", dtype=str, keep_default_na=False) \
        if search_path.exists() else pd.DataFrame()
    paper_idx = papers.set_index("paper_id") if not papers.empty else pd.DataFrame()

    email = lit.get("ncbi_email") or ""
    api_key = lit.get("ncbi_api_key")
    delay = float(lit.get("request_delay_seconds", 0.4))
    do_fetch = "ncbi_protein" in (lit.get("sequence_mapping_sources") or [])

    rows = []
    # cache accessions per paper to avoid rescanning
    paper_accs: dict[str, list[str]] = {}
    for _, c in cand.iterrows():
        pid = c["paper_id"]
        if pid not in paper_accs:
            text = ""
            if not paper_idx.empty and pid in paper_idx.index:
                prow = paper_idx.loc[pid]
                if isinstance(prow, pd.DataFrame):
                    prow = prow.iloc[0]
                text = _paper_text(prow, pid, fulltext_dir)
            paper_accs[pid] = list(dict.fromkeys(PROTEIN_ACC.findall(text)))[:10]
        accs = paper_accs[pid]
        acc = c.get("accession") or (accs[0] if accs else "")
        seq, conf, note = "", "none", "no accession found in paper text"
        db = ""
        if acc:
            db = "ncbi_protein"
            conf = "low_inferred_from_paper_text"
            note = "accession inferred from paper text (verify it is THIS protein)"
            if do_fetch:
                seq = _fetch_protein_fasta(cfg, acc, email, api_key, delay)
                if seq:
                    note = "sequence fetched from NCBI protein"
        rows.append({
            "candidate_id": c["candidate_id"], "paper_id": pid,
            "gene_name": c.get("gene_name", ""),
            "source_accession": acc, "source_database": db,
            "sequence": seq, "sequence_length": len(seq),
            "stated_or_inferred": "inferred" if acc else "none",
            "mapping_confidence": conf,
            "other_accessions_in_paper": "; ".join(accs[:5]),
            "note": note,
        })

    sm = pd.DataFrame(rows)
    sm.to_csv(out_dir / "sequence_mapping_table.tsv", sep="\t", index=False)
    n_mapped = int((sm["sequence_length"] > 0).sum()) if not sm.empty else 0
    log.info("Stage 0: sequence mapping — %d/%d candidate(s) linked to a sequence "
             "(all low-confidence/inferred; verify before use).",
             n_mapped, len(sm))
    return sm
