"""Stages 4 & 5 — profile HMM construction and scanning.

Builds three families of models from the Stage 3 alignments:
  * universal_holin       (Model A) — all gold positives
  * holin_<k>TMD          (Model B) — topology subsets
  * family_cluster<n>     (Model C) — per-cluster/family subsets

then scans them against every dataset class. Backend: pyhmmer (preferred,
pure-pip, works cross-platform) -> native HMMER CLI if present -> import a
pre-computed domtblout. If no backend is available, the pipeline continues with
an empty hit table and the downstream scores simply omit HMM evidence (logged).

A poor/empty universal model is an expected, reportable outcome — not an error.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from . import alignment as alignment_mod, validate
from .config import Config
from .utils import log, read_fasta

GAP = "-"


def _s(x) -> str:
    """Decode pyhmmer names that may be bytes or str depending on version."""
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


def _pyhmmer():
    try:
        import pyhmmer  # noqa: F401
        return pyhmmer
    except Exception:
        return None


def _model_type_from_id(model_id: str) -> str:
    if model_id.startswith("universal"):
        return "universal"
    if model_id.startswith("topology") or "TMD" in model_id:
        return "topology"
    return "family"


# --------------------------------------------------------- build ---------------
def build(cfg: Config, alignment_report: pd.DataFrame | None = None) -> dict:
    """Build HMMs from alignment .afa files. Returns dict with models + metadata."""
    py = _pyhmmer()
    aln_dir = cfg.out_dir("alignments")
    hmm_dir = cfg.out_dir("hmm")
    hcfg = cfg.section("hmmer")

    afa_files = sorted(aln_dir.glob("*.afa"))
    if not afa_files:
        log.warning("No alignments found in %s; run `align` first. No HMMs built.", aln_dir)

    want_universal = hcfg.get("build_universal", True)
    want_topology = hcfg.get("build_topology", True)
    want_family = hcfg.get("build_family", True)

    metadata_rows = []
    models = []  # list of dicts: {model_id, model_type, hmm, members}

    if py is None and not shutil.which("hmmbuild"):
        log.warning("No HMM backend (pyhmmer / hmmbuild) available; skipping HMM build.")
    else:
        alphabet = py.easel.Alphabet.amino() if py else None
        background = py.plan7.Background(alphabet) if py else None
        builder = py.plan7.Builder(alphabet) if py else None

        for afa in afa_files:
            model_id = afa.stem
            mtype = _model_type_from_id(model_id)
            if mtype == "universal" and not want_universal:
                continue
            if mtype == "topology" and not want_topology:
                continue
            if mtype == "family" and not want_family:
                continue
            aligned = read_fasta(afa)
            members = [pid for pid, _ in aligned]
            mean_len = round(sum(len(s.replace(GAP, "")) for _, s in aligned) / len(aligned), 1)
            if py is not None:
                try:
                    hmm = _build_hmm_pyhmmer(py, alphabet, builder, background, model_id, aligned)
                except Exception as exc:  # pragma: no cover - defensive
                    log.warning("Failed to build HMM %s: %s", model_id, exc)
                    continue
                # NB: pyhmmer's ASCII HMM serialization deadlocks on some
                # Windows builds (0.12.x); binary serialization is reliable.
                # We always persist a binary .h3m; ASCII .hmm is opt-in.
                with (hmm_dir / f"{model_id}.h3m").open("wb") as fh:
                    hmm.write(fh, binary=True)
                if hcfg.get("write_text_hmm", False):
                    try:
                        with (hmm_dir / f"{model_id}.hmm").open("wb") as fh:
                            hmm.write(fh, binary=False)
                    except Exception as exc:  # pragma: no cover
                        log.warning("Text HMM write failed for %s: %s", model_id, exc)
                models.append({"model_id": model_id, "model_type": mtype,
                               "hmm": hmm, "members": members})
                m_length = hmm.M
            else:
                # native hmmbuild path
                _build_hmm_cli(afa, hmm_dir / f"{model_id}.hmm")
                models.append({"model_id": model_id, "model_type": mtype,
                               "hmm": None, "members": members,
                               "hmm_path": hmm_dir / f"{model_id}.hmm"})
                m_length = None

            qrow = {}
            if alignment_report is not None and not alignment_report.empty:
                match = alignment_report[alignment_report["group_id"] == model_id]
                if not match.empty:
                    qrow = match.iloc[0].to_dict()
            metadata_rows.append({
                "model_id": model_id, "model_type": mtype,
                "num_sequences": len(members),
                "model_length": m_length,
                "alignment_length": qrow.get("alignment_length", ""),
                "mean_length": mean_len,
                "mean_pairwise_identity": qrow.get("mean_pairwise_identity", ""),
                "notes": ("LOW universal identity" if (mtype == "universal" and
                          qrow.get("mean_pairwise_identity", 100) < 25) else ""),
            })

    # combined, pressed HMM database (binary .h3{m,i,f,p}) via hmmpress
    if models and py is not None:
        combined = hmm_dir / "holin_models.hmm"
        try:
            py.hmmer.hmmpress([m["hmm"] for m in models], combined)
            log.info("hmmpress: wrote pressed database %s.h3{m,i,f,p}", combined.name)
        except Exception as exc:
            log.info("hmmpress skipped (%s).", type(exc).__name__)

    meta = pd.DataFrame(metadata_rows, columns=["model_id", "model_type", "num_sequences",
                                                "model_length", "alignment_length",
                                                "mean_length", "mean_pairwise_identity", "notes"])
    meta.to_csv(cfg.out("tables", "hmm_model_metadata.tsv"), sep="\t", index=False)
    log.info("Stage 4: built %d HMM(s) (%s).", len(models),
             ", ".join(sorted({m["model_type"] for m in models})) or "none")
    return {"models": models, "metadata": meta}


def _build_hmm_pyhmmer(py, alphabet, builder, background, name, aligned):
    TextMSA = py.easel.TextMSA
    TextSequence = py.easel.TextSequence
    seqs = [TextSequence(name=pid.encode(), sequence=seq) for pid, seq in aligned]
    msa = TextMSA(name=name.encode(), sequences=seqs)
    dmsa = msa.digitize(alphabet)
    hmm, _profile, _opt = builder.build_msa(dmsa, background)
    hmm.name = name.encode()
    return hmm


def _build_hmm_cli(afa: Path, out: Path) -> None:
    import subprocess
    subprocess.run(["hmmbuild", "--amino", str(out), str(afa)],
                   capture_output=True, text=True)


# --------------------------------------------------------- scan ----------------
def build_single(cfg: Config, model_id: str, model_type: str,
                 records: list[tuple[str, str]]) -> dict | None:
    """Align + build a single in-memory HMM (no disk writes).

    Used by benchmark.py for leakage-free cross-validation, where models must be
    rebuilt per fold from the training sequences only.
    """
    py = _pyhmmer()
    if py is None or len(records) < 2:
        return None
    aligned = alignment_mod.align_group(records, cfg)
    if not aligned:
        return None
    alphabet = py.easel.Alphabet.amino()
    background = py.plan7.Background(alphabet)
    builder = py.plan7.Builder(alphabet)
    try:
        hmm = _build_hmm_pyhmmer(py, alphabet, builder, background, model_id, aligned)
    except Exception:
        return None
    return {"model_id": model_id, "model_type": model_type, "hmm": hmm,
            "members": [p for p, _ in aligned]}


def scan(cfg: Config, models: list[dict] | None = None,
         clean: pd.DataFrame | None = None,
         write_output: bool = True) -> pd.DataFrame:
    py = _pyhmmer()
    hcfg = cfg.section("hmmer")
    evalue = float(hcfg.get("evalue_cutoff", 0.01))

    if clean is None:
        clean = validate.run(cfg)["clean"]
    records = list(zip(clean["protein_id"], clean["sequence"]))
    cat = dict(zip(clean["protein_id"], clean["dataset_category"]))
    lengths = {pid: len(seq) for pid, seq in records}

    cols = ["protein_id", "dataset_category", "model_id", "model_type", "bitscore",
            "evalue", "bias", "ali_from", "ali_to", "hmm_from", "hmm_to",
            "model_length", "target_length", "query_coverage", "target_coverage"]

    if py is None or not models:
        if write_output:
            log.warning("HMM scan skipped (no pyhmmer backend or no models). "
                        "Downstream scores will omit HMM evidence.")
        empty = pd.DataFrame(columns=cols)
        if write_output:
            empty.to_csv(cfg.out("tables", "hmm_scan_results.tsv"), sep="\t", index=False)
        return empty

    alphabet = py.easel.Alphabet.amino()
    TextSequence = py.easel.TextSequence
    DigitalSequenceBlock = py.easel.DigitalSequenceBlock
    digital = [TextSequence(name=pid.encode(), sequence=seq).digitize(alphabet)
               for pid, seq in records]
    block = DigitalSequenceBlock(alphabet, digital)

    rows = []
    for m in models:
        hmm = m["hmm"]
        if hmm is None:
            continue
        for top_hits in py.hmmer.hmmsearch([hmm], block, E=evalue):
            for hit in top_hits:
                if hit.evalue > evalue:
                    continue
                target = _s(hit.name)
                dom = hit.best_domain
                aln = dom.alignment if dom is not None else None
                ali_from = getattr(aln, "target_from", 0) if aln else 0
                ali_to = getattr(aln, "target_to", 0) if aln else 0
                hmm_from = getattr(aln, "hmm_from", 0) if aln else 0
                hmm_to = getattr(aln, "hmm_to", 0) if aln else 0
                mlen = hmm.M
                tlen = lengths.get(target, 0) or 1
                qcov = round((hmm_to - hmm_from + 1) / mlen, 3) if mlen else 0.0
                tcov = round((ali_to - ali_from + 1) / tlen, 3) if tlen else 0.0
                rows.append({
                    "protein_id": target, "dataset_category": cat.get(target, ""),
                    "model_id": m["model_id"], "model_type": m["model_type"],
                    "bitscore": round(float(hit.score), 2),
                    "evalue": float(hit.evalue),
                    "bias": round(float(getattr(hit, "bias", 0.0)), 2),
                    "ali_from": ali_from, "ali_to": ali_to,
                    "hmm_from": hmm_from, "hmm_to": hmm_to,
                    "model_length": mlen, "target_length": tlen,
                    "query_coverage": max(0.0, qcov), "target_coverage": max(0.0, tcov),
                })

    hits = pd.DataFrame(rows, columns=cols)
    if write_output:
        hits.to_csv(cfg.out("tables", "hmm_scan_results.tsv"), sep="\t", index=False)
        log.info("Stage 5: %d HMM hit(s) across %d model(s) and %d target(s).",
                 len(hits), len(models), len(records))
    return hits


def best_hits_per_protein(hits: pd.DataFrame) -> pd.DataFrame:
    if hits.empty:
        return hits
    idx = hits.groupby("protein_id")["bitscore"].idxmax()
    return hits.loc[idx].reset_index(drop=True)
