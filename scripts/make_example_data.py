#!/usr/bin/env python
"""Generate the toy/placeholder example dataset under data/example/.

The sequences here are SYNTHETIC. They are constructed to have holin-like
*architecture* (small size, high hydrophobic fraction, a controllable number of
transmembrane-like hydrophobic segments) so the pipeline can be exercised
end-to-end. They are NOT real holins and must not be interpreted scientifically.
Replace data/example/ with real curated data (e.g. via Stage 0 literature mining
+ manual review) before drawing any biological conclusions.

Deterministic: seeded RNG -> identical output every run.
"""
from __future__ import annotations

import csv
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "example"
OUT.mkdir(parents=True, exist_ok=True)

RNG = random.Random(1729)

HYDROPHOBIC = "LIVFAMW"           # transmembrane-segment residues
POLAR = "STNQGYHP"                # loop residues
POS = "KR"                         # positively charged
NEG = "DE"                         # negatively charged


def tm_segment(n: int = 20) -> str:
    return "".join(RNG.choice(HYDROPHOBIC) for _ in range(n))


def loop(n: int = 6, charged: bool = True) -> str:
    pool = POLAR + (POS + NEG if charged else "")
    return "".join(RNG.choice(pool) for _ in range(n))


def make_seq(n_tmd: int, *, nterm_charge: int = 2, ctail: int = 12,
             soluble: bool = False, length_target: int | None = None) -> str:
    """Build a synthetic protein with `n_tmd` hydrophobic segments."""
    parts = ["M"]
    parts.append("".join(RNG.choice(POS + NEG) for _ in range(nterm_charge)))
    parts.append(loop(4))
    if soluble:
        # mostly polar/charged, no long hydrophobic stretch
        body = "".join(RNG.choice(POLAR + POS + NEG) for _ in range(length_target or 120))
        return "M" + body
    for i in range(n_tmd):
        parts.append(tm_segment(RNG.randint(18, 22)))
        if i < n_tmd - 1:
            parts.append(loop(RNG.randint(5, 9)))
    parts.append(loop(ctail, charged=True))
    seq = "".join(parts)
    return seq


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", newline="\n") as fh:
        for pid, seq in records:
            fh.write(f">{pid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


# ---------------------------------------------------------------- gold --------
# (protein_id, name, phage, host, accession, n_tmd, evidence_type, citation,
#  holin_type, family_label, notes)
gold_spec = [
    ("gold_lambdaS_01", "S holin (placeholder)", "Lambda-like", "Escherichia coli",
     "PLACEHOLDER_G1", 2, "amber_mutant_delayed_lysis", "Placeholder et al. 1990",
     "class_I", "lambdaS_like", "SYNTHETIC placeholder sequence"),
    ("gold_phi21S_02", "S pinholin (placeholder)", "Phage21-like", "Escherichia coli",
     "PLACEHOLDER_G2", 2, "membrane_depolarization", "Placeholder et al. 2000",
     "pinholin", "phi21S_like", "SYNTHETIC placeholder sequence"),
    ("gold_T4t_03", "t holin (placeholder)", "T4-like", "Escherichia coli",
     "PLACEHOLDER_G3", 1, "complementation_restored_lysis", "Placeholder et al. 1997",
     "class_III", "T4t_like", "SYNTHETIC placeholder sequence"),
    ("gold_P2Y_04", "Y holin (placeholder)", "P2-like", "Escherichia coli",
     "PLACEHOLDER_G4", 2, "deletion_mutant_lysis_defect", "Placeholder et al. 2003",
     "canonical_holin", "P2Y_like", "SYNTHETIC placeholder sequence"),
    ("gold_mycoB_05", "holin (placeholder)", "Mycobacteriophage-like", "Mycobacterium smegmatis",
     "PLACEHOLDER_G5", 3, "endolysin_dependent_lysis", "Placeholder et al. 2011",
     "canonical_holin", "myco_like", "SYNTHETIC placeholder sequence"),
    ("gold_T7gp17_06", "holin (placeholder)", "T7-like", "Escherichia coli",
     "PLACEHOLDER_G6", 1, "lysis_timing_mutation", "Placeholder et al. 2008",
     "class_II", "T7_like", "SYNTHETIC placeholder sequence"),
    ("gold_class3_07", "holin (placeholder)", "Phage-like-X", "Bacillus subtilis",
     "PLACEHOLDER_G7", 4, "required_for_lysis", "Placeholder et al. 2015",
     "canonical_holin", "bsuX_like", "SYNTHETIC placeholder sequence"),
    # Deliberately missing citation + evidence_type to exercise QC flags.
    ("gold_nocite_08", "holin (placeholder)", "Phage-like-Y", "Escherichia coli",
     "PLACEHOLDER_G8", 2, "", "", "unknown", "", "SYNTHETIC placeholder; missing evidence"),
]

gold_rows, gold_fasta = [], []
for pid, name, phage, host, acc, ntmd, ev, cite, htype, fam, notes in gold_spec:
    seq = make_seq(ntmd, nterm_charge=2, ctail=14)
    gold_rows.append({
        "protein_id": pid, "protein_name": name, "phage_name": phage, "host": host,
        "accession": acc, "sequence": seq, "evidence_type": ev, "citation": cite,
        "holin_type": htype, "family_label": fam, "notes": notes,
    })
    gold_fasta.append((pid, seq))

with (OUT / "gold_holins.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["protein_id", "protein_name", "phage_name", "host",
                                       "accession", "sequence", "evidence_type", "citation",
                                       "holin_type", "family_label", "notes"])
    w.writeheader()
    w.writerows(gold_rows)

# ---------------------------------------------------------------- weak --------
weak_spec = [
    ("weak_put_01", "putative holin", "PhageW1", "PLACEHOLDER_W1", 2, "putative holin", "RefSeq"),
    ("weak_put_02", "phage holin family protein", "PhageW2", "PLACEHOLDER_W2", 1, "holin family", "Pfam"),
    ("weak_put_03", "holin-like protein", "PhageW3", "PLACEHOLDER_W3", 3, "holin-like", "RefSeq"),
    ("weak_put_04", "putative pinholin", "PhageW4", "PLACEHOLDER_W4", 2, "putative pinholin", "RefSeq"),
    ("weak_put_05", "hypothetical protein (holin?)", "PhageW5", "PLACEHOLDER_W5", 1, "hypothetical", "GenBank"),
]
weak_rows = []
for pid, name, phage, acc, ntmd, ann, db in weak_spec:
    seq = make_seq(ntmd, nterm_charge=1, ctail=10)
    weak_rows.append({"protein_id": pid, "protein_name": name, "phage_name": phage,
                      "accession": acc, "sequence": seq, "annotation": ann,
                      "source_database": db, "notes": "SYNTHETIC placeholder"})
with (OUT / "weak_annotated_holins.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["protein_id", "protein_name", "phage_name", "accession",
                                       "sequence", "annotation", "source_database", "notes"])
    w.writeheader()
    w.writerows(weak_rows)

# ------------------------------------------------------------ hard negatives --
# These resemble holins superficially (small, membrane) but are not holins,
# OR have alternative annotations inconsistent with holin function.
neg_spec = [
    ("neg_spanin_i_01", "i-spanin component", "phage_membrane_protein", 1, "inner membrane spanin"),
    ("neg_tail_02", "tail membrane protein", "tail_membrane_protein", 1, "tail assembly protein"),
    ("neg_smemb_03", "small bacterial membrane protein", "small_bacterial_membrane", 2, "DUF membrane protein"),
    ("neg_TAtoxin_04", "toxin-antitoxin membrane toxin", "toxin_antitoxin", 1, "TA system toxin"),
    ("neg_transp_05", "transporter fragment", "transporter_fragment", 4, "MFS transporter"),
    ("neg_hypo_06", "hypothetical protein", "hypothetical_no_lysis_context", 2, "hypothetical protein"),
    ("neg_struct_07", "virion structural protein", "structural_virion", 1, "major capsid protein"),
    ("neg_amidase_08", "N-acetylmuramoyl-L-alanine amidase", "enzymatic_domain", 0, "amidase (enzyme)"),
    ("neg_soluble_09", "cytoplasmic protein", "soluble_protein", 0, "soluble cytoplasmic protein"),
]
neg_rows, neg_fasta = [], []
for pid, name, ntype, ntmd, ann in neg_spec:
    if ntmd == 0:
        seq = make_seq(0, soluble=True, length_target=RNG.randint(150, 240))
    else:
        seq = make_seq(ntmd, nterm_charge=2, ctail=10)
    neg_rows.append({"protein_id": pid, "protein_name": name, "source": "synthetic",
                     "accession": f"PLACEHOLDER_{pid.upper()}", "sequence": seq,
                     "negative_type": ntype, "annotation": ann, "notes": "SYNTHETIC placeholder"})
    neg_fasta.append((pid, seq))

# Add an exact-duplicate of a gold sequence appearing as a negative, to exercise
# the "duplicate sequence across classes" QC flag.
dup_seq = gold_rows[0]["sequence"]
neg_rows.append({"protein_id": "neg_dupgold_10", "protein_name": "duplicate-of-gold (QC test)",
                 "source": "synthetic", "accession": "PLACEHOLDER_DUP",
                 "sequence": dup_seq, "negative_type": "qc_duplicate",
                 "annotation": "intentional cross-class duplicate", "notes": "QC flag test"})
neg_fasta.append(("neg_dupgold_10", dup_seq))

with (OUT / "hard_negatives.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["protein_id", "protein_name", "source", "accession",
                                       "sequence", "negative_type", "annotation", "notes"])
    w.writeheader()
    w.writerows(neg_rows)

# ------------------------------------------------- unknown candidates (.faa) --
# A mix: some holin-like (and near lysis genes in context.tsv), some not.
cand_spec = [
    ("cand_hi_01", 2),     # holin-like, in a lysis cassette -> should rank high
    ("cand_hi_02", 1),     # holin-like pinholin-ish, near SAR endolysin
    ("cand_med_03", 2),    # holin-like but isolated from lysis genes
    ("cand_med_04", 3),    # holin-like topology, weak context
    ("cand_low_05", 4),    # 4 TMD transporter-ish, alternative annotation in context
    ("cand_low_06", 0),    # soluble, unlikely
    ("cand_unl_07", 0),    # large soluble, unlikely
]
cand_fasta = []
for pid, ntmd in cand_spec:
    if ntmd == 0:
        seq = make_seq(0, soluble=True, length_target=RNG.randint(160, 260))
    else:
        seq = make_seq(ntmd, nterm_charge=2, ctail=12)
    cand_fasta.append((pid, seq))

# proteins.faa = unknown candidates to scan. Include all four dataset members too
# so `scan` can be demonstrated against everything; but the canonical "unknowns"
# are the cand_* records. Keep proteins.faa = candidates only (others have CSVs).
write_fasta(OUT / "proteins.faa", cand_fasta)

# ---------------------------------------------------- topology_predictions ----
# One row per protein across ALL datasets (import path for Stage 6).
def topo_rows():
    everything = []
    for r, ntmd in zip(gold_rows, [s[5] for s in gold_spec]):
        everything.append((r["protein_id"], r["sequence"], ntmd))
    for r, spec in zip(weak_rows, weak_spec):
        everything.append((r["protein_id"], r["sequence"], spec[4]))
    # negatives: ntmd from spec, dup uses gold[0] ntmd
    neg_ntmds = [s[3] for s in neg_spec] + [gold_spec[0][5]]
    for r, ntmd in zip(neg_rows, neg_ntmds):
        everything.append((r["protein_id"], r["sequence"], ntmd))
    for (pid, seq), (_, ntmd) in zip(cand_fasta, cand_spec):
        everything.append((pid, seq, ntmd))
    return everything

with (OUT / "topology_predictions.tsv").open("w", newline="") as fh:
    fh.write("protein_id\tlength\ttmd_count\ttopology\tn_region\tc_region\tsignal_peptide\tsar_like\ttool\n")
    for pid, seq, ntmd in topo_rows():
        topo = "in" if ntmd % 2 == 0 else "out"
        sar = "yes" if ("phi21S" in pid or "cand_hi_02" in pid) else "no"
        sp = "no"
        fh.write(f"{pid}\t{len(seq)}\t{ntmd}\t{topo}\t1-30\t{max(1,len(seq)-30)}-{len(seq)}\t{sp}\t{sar}\tplaceholder_topology\n")

# ------------------------------------------------------------- context.tsv ----
# A lysis cassette on contig_A around cand_hi_01 / cand_hi_02, plus an isolated
# candidate on contig_B. Columns: contig_id,gene_id,start,end,strand,protein_id,product,sequence
def seqof(pid):
    for p, s in cand_fasta + gold_fasta + neg_fasta:
        if p == pid:
            return s
    return ""

context = [
    # contig_A: structural ... holin candidate, endolysin, spanin(Rz/Rz1)
    ("contig_A", "A1", 100, 1300, "+", "neg_struct_07", "major capsid protein", seqof("neg_struct_07")),
    ("contig_A", "A2", 1350, 1700, "+", "cand_hi_01", "hypothetical protein", seqof("cand_hi_01")),
    ("contig_A", "A3", 1720, 2200, "+", "endolysin_A", "endolysin", ""),
    ("contig_A", "A4", 2220, 2600, "+", "spanin_A", "Rz spanin", ""),
    ("contig_A", "A5", 2300, 2520, "+", "spanin_A2", "Rz1 spanin", ""),
    ("contig_A", "A6", 2650, 3000, "+", "neg_hypo_06", "hypothetical protein", seqof("neg_hypo_06")),
    # contig_A second cassette with SAR endolysin near a pinholin-like candidate
    ("contig_A", "A7", 3100, 3400, "+", "cand_hi_02", "hypothetical protein", seqof("cand_hi_02")),
    ("contig_A", "A8", 3420, 3950, "+", "sar_endolysin_A", "SAR endolysin", ""),
    # contig_B: isolated candidate, no lysis genes nearby
    ("contig_B", "B1", 100, 700, "+", "neg_tail_02", "tail protein", seqof("neg_tail_02")),
    ("contig_B", "B2", 720, 1050, "+", "cand_med_03", "hypothetical protein", seqof("cand_med_03")),
    ("contig_B", "B3", 1100, 2400, "+", "neg_transp_05", "MFS transporter", seqof("neg_transp_05")),
    ("contig_B", "B4", 2450, 2900, "+", "cand_low_05", "transporter fragment", seqof("cand_low_05")),
]
with (OUT / "context.tsv").open("w", newline="") as fh:
    fh.write("contig_id\tgene_id\tstart\tend\tstrand\tprotein_id\tproduct\tsequence\n")
    for row in context:
        fh.write("\t".join(str(x) for x in row) + "\n")

# ------------------------------------------------------ lysis_context_terms ---
(OUT / "lysis_context_terms.txt").write_text(
    "\n".join([
        "# Lysis-related annotation keywords (case-insensitive substring match).",
        "# category<TAB>term  (category drives context scoring)",
        "endolysin\tendolysin",
        "endolysin\tlysin",
        "endolysin\tlysozyme",
        "endolysin\tamidase",
        "endolysin\tglycosidase",
        "endolysin\tmuramidase",
        "endolysin\tpeptidoglycan hydrolase",
        "sar_endolysin\tSAR endolysin",
        "spanin\tspanin",
        "spanin\tRz",
        "spanin\tRz1",
        "spanin\tunimolecular spanin",
        "antiholin\tantiholin",
        "holin\tholin",
        "holin\tpinholin",
    ]) + "\n",
    encoding="utf-8",
)

# ---------------------------------------------------------- search_queries ----
(OUT / "search_queries.txt").write_text(
    "\n".join([
        "# Literature-mining query groups for Stage 0 (one query per line).",
        "# group<TAB>query",
        "general\tbacteriophage holin characterization lysis",
        "general\tphage holin lysis timing",
        "general\tholin endolysin membrane permeabilization",
        "general\tholin defective mutant phage lysis",
        "pinholin\tpinholin SAR endolysin",
        "pinholin\tbacteriophage pinholin membrane depolarization",
        "model_systems\tlambda S holin lysis S105 S107",
        "model_systems\tphage 21 pinholin",
        "model_systems\tmycobacteriophage holin lysis cassette",
        "evidence\tholin amber mutant delayed lysis",
        "evidence\tholin complementation restored lysis",
        "evidence\tholin membrane depolarization potential",
    ]) + "\n",
    encoding="utf-8",
)

print("Wrote example data to", OUT)
for p in sorted(OUT.iterdir()):
    print("  ", p.name)
