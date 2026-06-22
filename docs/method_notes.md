# Method notes

Technical notes on how each stage is implemented and the deliberate design choices that keep the
analysis conservative and reproducible.

## Backends and portability

The default environment is cross-platform and pure-pip. Each external-tool stage is backend-pluggable
and degrades gracefully:

| Stage | Preferred | Fallbacks |
| --- | --- | --- |
| HMM build/scan (`hmmer.py`) | `pyhmmer` | native HMMER CLI → import domtbl → skip (empty hits) |
| MSA (`alignment.py`) | `pyfamsa` | MAFFT / Clustal Omega CLI → import `.afa` |
| Clustering (`clustering.py`) | built-in greedy identity clusterer (BLOSUM62 pairwise via Biopython) | MMseqs2 / CD-HIT CLI |
| Topology (`topology.py`, `features.py`) | **import** `topology_predictions.tsv` | crude Kyte–Doolittle estimate, tagged `builtin_kd_estimate` |

If a backend is missing, the pipeline still completes; downstream scores omit that evidence and the
report notes it.

### pyhmmer Windows caveat

On some `pyhmmer` 0.12.x Windows builds, ASCII HMM serialization (`HMM.write(fh)`) deadlocks, while
binary serialization and `hmmpress` work. We therefore persist models as binary `.h3m` plus a pressed
combined database (`results/hmm/holin_models.hmm.h3{m,i,f,p}`). ASCII `.hmm` text output is opt-in via
`hmmer.write_text_hmm: true`. Scanning uses in-memory HMM objects, so no text round-trip is needed.

## Stage 1 — validation

Sequences are uppercased and stripped to letters. Clearly invalid sequences (non-amino-acid
characters, > `max_ambiguous_fraction` ambiguous residues, shorter than `hard_min_length`) are
dropped; everything else is *flagged*, not removed. Cross-class duplicate sequences and within-class
exact duplicates are detected. Gold records missing `citation`/`evidence_type` are flagged.

## Stage 2 — clustering

Greedy, longest-first clustering at 90/70/50/30 % identity. Identity is alignment-based (Biopython
`PairwiseAligner`, BLOSUM62, global) as `identities / len(shorter)`, with a k-mer Jaccard fallback if
alignment fails. The headline metric is the number of clusters at the primary threshold: many
clusters ⇒ fragmented sequence space ⇒ a universal HMM is unlikely to generalize.

## Stage 3 — alignment quality

For every MSA we compute mean pairwise identity, % gaps, alignment length, conserved columns, and
mean per-column entropy. The point is to **quantify** how poor the universal alignment is rather than
hide it.

## Stages 4–5 — HMMs

Universal, per-topology, and per-cluster HMMs are built from the Stage 3 alignments and scanned
against all dataset classes. Derived metrics include best-hit-per-protein, query/target coverage, and
number of models hit.

## Stage 6 — features / architecture

Length, MW, hydrophobic fraction (set `AILMFWVC`), positive/negative charge, N30/C30 charge
distribution, TMD count (imported or estimated), SAR-like flag, low-complexity proxy, and an
enzymatic/structural-domain flag derived from annotation keywords.

## Stage 7 — genomic context

Within a ±gene / ±bp window, neighbours are classified using the lysis-term lexicon. An **additive,
fully configurable** context score rewards proximity to endolysins/spanins, same-strand
co-orientation, and SAR+pinholin pairing, and penalizes alternative annotations and isolation.

## Stage 8 — composite score

Three normalized sub-scores (HMM, architecture, context) are combined with configurable weights into
`final_holin_score ∈ [0,1]`, then binned into confidence categories. Every candidate carries a
human-readable `explanation` describing why it was ranked as it was. This is a ranking system, not a
classifier.

## Stage 9 — benchmarking

Seven approaches (A–G) are scored on gold vs hard negatives in two regimes:

- **Naive** — full-data models; HMM models are circular here and labeled with `circularity_warning`.
- **Leave-one-family-out** — the universal HMM is rebuilt excluding each family and asked to recover
  the held-out family while we count false hits on hard negatives. This is the honest generalization
  test and the one that answers the project's central question.

Metrics: recall/specificity/precision/F1, ROC-AUC, PR-AUC (with bootstrap CIs), top-k recovery,
confusion matrix, and explicit false-positive / false-negative tables. Small class sizes are warned.

## Stage 10 — motif analysis

Per-column information content over each alignment; a column is "conserved" above
`entropy_conserved_bits`. We then report the fraction of conserved columns whose consensus residue is
hydrophobic — if high, the signal is architecture-level (TM hydrophobicity), **not** a holin motif.
Sequence logos are produced with `logomaker` where available.

## Stage 0 — literature mining

PubMed E-utilities + Europe PMC, with on-disk HTTP caching for reproducibility. Evidence is scored
0–5 from trigger-term taxonomy (genetic, functional-lysis, membrane, endolysin-dependent). Nothing
reaches gold without `evidence_score >= 4` **and** manual verification; the export command enforces
this gate.

## Reproducibility

- A single random seed (`project.random_seed`) governs stochastic steps (bootstrap).
- HTTP responses are cached under `data/litcache/` (gitignored) so re-runs are deterministic.
- The example dataset is regenerated deterministically by `scripts/make_example_data.py`.
