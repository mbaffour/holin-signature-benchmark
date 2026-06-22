# holin-signature-benchmark

**Benchmarking universal, topology-specific, and context-aware signatures for bacteriophage holin prediction.**

This repository is a reproducible, conservative bioinformatics pipeline that tests a single
question:

> *Do experimentally characterized bacteriophage holins contain detectable universal
> sequence-level signatures, and how do universal HMMs compare against topology-specific,
> family-specific, and genomic-context-aware approaches for holin prediction?*

It is **designed to be able to fail**. Holins are small, hydrophobic, 1–4 TMD membrane proteins
that schedule host lysis, but they are extraordinarily diverse and may share function without
deep sequence homology. A universal holin HMM may therefore perform poorly — and if it does, the
pipeline reports that clearly rather than overclaiming. See [`PROJECT_SPEC.md`](PROJECT_SPEC.md)
for the full scientific brief.

## Why a universal holin motif may not exist

- Holins span many unrelated families; characterized holins often occupy *isolated* sequence
  islands rather than one connected sequence space (quantified in Stage 2).
- "Small + hydrophobic + has TMDs" describes a huge set of unrelated membrane proteins, so a
  universal model easily matches **hard negatives** (Stage 9 makes this explicit).
- Conservation seen in a holin alignment is frequently just generic transmembrane hydrophobicity,
  not a holin-specific motif (Stage 10 distinguishes the two).

The pipeline therefore treats the final output as a **ranking / prioritization** system, never a
proof of function.

## What it does (stages)

| Stage | Command | Output |
| --- | --- | --- |
| 0 | `literature-search`, `extract-evidence`, `map-sequences`, `prepare-review`, `export-curated` | Curation aids for building a gold set from the literature |
| 1 | `validate` | Cleaned datasets + QC flags |
| 2 | `cluster` | Gold clustering at 90/70/50/30 % identity; fragmentation summary |
| 3 | `align` | MSAs + alignment-quality metrics (incl. how poor the universal MSA is) |
| 4 | `build-hmms` | Universal / topology / family HMMs (+ pressed DB) |
| 5 | `scan` | HMM hits across all dataset classes |
| 6 | `features` | Length, hydrophobicity, TMD count, charge, SAR-like, domain flags |
| 7 | `context` | Lysis-cassette / genomic-neighborhood scoring |
| 8 | `score` | Composite candidate ranking with per-candidate explanations |
| 9 | `benchmark` | Models A–G, naive **and** leave-one-family-out, ROC/PR-AUC + bootstrap CIs |
| 10 | `motifs` | Conservation/entropy + sequence logos, motif-vs-hydrophobicity call |
| 11 | `synteny` | Lysis-cassette gene-arrow maps for top candidates |
| 12 | `report` | Manuscript-style Markdown report |
| — | `run-all` | Stages 1–12 end-to-end |

## Installation

The default environment is **fully cross-platform and self-contained** (Windows/Linux/macOS): all
bioinformatics backends (HMMER, FAMSA) are pure-pip wheels (`pyhmmer`, `pyfamsa`), so no native
tools are required.

```bash
conda env create -f environment.yml
conda activate holinbench
pip install -e .
```

On Linux/macOS you may optionally add the canonical native tools (the pipeline auto-detects and
prefers them when present):

```bash
conda env create -f environment-native.yml   # hmmer, mafft, mmseqs2, diamond
```

### Required vs optional external tools

| Capability | Default backend (always works) | Optional native tool |
| --- | --- | --- |
| HMM build/search/press | `pyhmmer` | HMMER CLI |
| Multiple sequence alignment | `pyfamsa` | MAFFT / Clustal Omega |
| Clustering | built-in greedy identity clusterer | MMseqs2 / CD-HIT |
| Topology (TMD count) | **import** from `topology_predictions.tsv` (preferred) or crude built-in estimate | DeepTMHMM / Phobius / TMHMM / TOPCONS |

> **Windows note:** `pyhmmer`'s ASCII HMM serialization can deadlock on some builds, so models are
> persisted as binary `.h3m` + a pressed combined database. Text `.hmm` output is opt-in
> (`hmmer.write_text_hmm: true`) for platforms where it is known to work.

## Example run

```bash
# regenerate the synthetic example dataset (optional; already committed)
python scripts/make_example_data.py

# run the whole analytical pipeline on the bundled example data
holinbench run-all -c config/config.yaml

# outputs land in results/: tables/, figures/, alignments/, hmm/, sequences/, report.md
```

> The bundled `data/example/` sequences are **synthetic placeholders** for pipeline validation
> only. They are *not* real holins. Replace them with curated data before any biological
> interpretation.

### Literature mining (Stage 0)

```bash
holinbench literature-search  -c config/config.yaml   # PubMed + Europe PMC (cached)
holinbench extract-evidence   -c config/config.yaml   # conservative evidence classification
holinbench map-sequences      -c config/config.yaml   # best-effort accession -> sequence
holinbench prepare-review     -c config/config.yaml   # manual_review_template.tsv + summary
# ... a human edits manual_review_template.tsv, setting manually_verified=true for real holins ...
holinbench export-curated     -c config/config.yaml   # exports ONLY verified, strong-evidence rows
```

## Input file formats

All paths are configurable in [`config/config.yaml`](config/config.yaml).

- `data/example/gold_holins.csv` — `protein_id,protein_name,phage_name,host,accession,sequence,evidence_type,citation,holin_type,family_label,notes`
- `data/example/weak_annotated_holins.csv` — `protein_id,protein_name,phage_name,accession,sequence,annotation,source_database,notes`
- `data/example/hard_negatives.csv` — `protein_id,protein_name,source,accession,sequence,negative_type,annotation,notes`
- `data/example/proteins.faa` — FASTA of unknown candidates to scan
- `data/example/context.tsv` — `contig_id,gene_id,start,end,strand,protein_id,product,sequence`
- `data/example/topology_predictions.tsv` — `protein_id,length,tmd_count,topology,n_region,c_region,signal_peptide,sar_like,tool`
- `data/example/lysis_context_terms.txt` — `category<TAB>term`

## Output explanation

- `results/tables/` — cleaned metadata, QC flags, clusters, alignment quality, HMM metadata + hits,
  features, context, candidate ranking, benchmark metrics (+ LOFO), false positives/negatives,
  motif/conservation summary.
- `results/figures/` — the manuscript figures (distributions, cluster fragmentation, model
  comparison, ROC/PR, confusion matrix, sequence logos, synteny maps, score distribution).
- `results/report.md` — the manuscript-style report with data-driven conclusions.
- `results/literature/` — Stage 0 curation tables + `curation_summary.md`.

## Interpretation warnings

- A high candidate score **does not prove** holin function — it prioritizes for experiments.
- HMM hits may reflect family similarity, not a universal holin signature.
- TMDs and hydrophobicity are necessary-ish but **not sufficient**.
- Genomic context is **supportive**, not definitive.
- The naive benchmark regime is **circular** (HMMs built on the same gold); trust the
  leave-one-family-out regime for generalization claims. See
  [`docs/interpretation_guide.md`](docs/interpretation_guide.md).

## How to add new experimentally characterized holins

1. Confirm the protein meets the gold-positive bar in
   [`docs/data_curation_guidelines.md`](docs/data_curation_guidelines.md) (direct experimental
   evidence — not annotation).
2. Add a row to `data/example/gold_holins.csv` (or your own curated CSV; update the path in
   `config.yaml`). Include `evidence_type` and `citation`.
3. Re-run `holinbench run-all`. Stage 1 will flag any gold record missing citation/evidence.

## How to run benchmarking

`holinbench benchmark -c config/config.yaml` produces `results/tables/benchmark_metrics.tsv`
(naive) and `results/tables/benchmark_lofo_universal.tsv` (honest, leave-one-family-out). Validation
strategies and bootstrap settings are configured under `benchmark:` in `config.yaml`.

## How to cite / use the output

This is a methods/benchmarking scaffold. Treat its rankings as hypotheses for experimental
validation (lysis-timing, membrane permeabilization, complementation, endolysin-dependent lysis).
The `results/report.md` is structured for direct adaptation into a methods manuscript.

## Tests

```bash
pytest -q
```

## License

MIT.
