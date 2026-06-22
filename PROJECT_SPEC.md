# Project Specification — Benchmarking universal, topology-specific, and context-aware signatures for bacteriophage holin prediction

> This file is the **canonical specification** for the project, stored verbatim. The pipeline
> implemented in this repository is built to satisfy this brief. It is intentionally
> **scientifically conservative**: it is designed to *test* whether a universal holin signature
> exists, not to assume one. A negative result is a valid and useful outcome.

---

You are helping me build a rigorous, reproducible bioinformatics project to test whether bacteriophage holins have a universal sequence signature, and to benchmark whether universal HMMs, topology-specific HMMs, family-specific HMMs, and context-aware models can identify characterized or annotated holins.

This project must be scientifically conservative. Do not assume there is a universal holin motif. The real goal is to test that hypothesis carefully. The expected outcome may be negative: a universal holin HMM may fail because holins are highly diverse and may share function without deep sequence homology. The pipeline should therefore be designed to answer:

“Do experimentally characterized bacteriophage holins contain detectable universal sequence-level signatures, and how do universal HMMs compare against topology-specific, family-specific, and genomic-context-aware approaches for holin prediction?”

Build the project as a reproducible bioinformatics repository suitable for a small methods/bioinformatics paper.

## Core scientific framing

* Holins are small phage-encoded membrane proteins involved in host lysis.
* They are often small, hydrophobic, and predicted to contain 1–4 transmembrane domains.
* They are frequently encoded near endolysins, SAR endolysins, spanins, antiholins, or other lysis genes.
* However, holins are highly diverse and may not share one universal sequence motif.
* Therefore, the project should test the limits of sequence-only prediction and compare it to architecture/context-aware prediction.
* The paper should not overclaim. The conclusion should be data-driven.

## The repository should include

1. A reproducible command-line pipeline.
2. Clear input data formats.
3. Scripts for building positive, weak-positive, negative, and unknown datasets.
4. HMM construction and benchmarking workflows.
5. Transmembrane topology prediction integration.
6. Genomic neighborhood/context scoring.
7. Validation workflows that avoid circular annotation and homology leakage.
8. Figures and tables suitable for a manuscript.
9. A final report summarizing results.
10. Unit tests or sanity checks for core functions.
11. Documentation explaining how to run everything.

Use Python as the main language. Use Snakemake or Nextflow if appropriate, but do not make the workflow overly complicated. A clean Python CLI plus optional Snakemake workflow is acceptable. Prefer modular code.

## Recommended external tools

* HMMER: hmmbuild, hmmsearch, hmmpress, hmmalign.
* MAFFT or Clustal Omega for multiple sequence alignment.
* MMseqs2 or CD-HIT for clustering and redundancy reduction.
* DeepTMHMM, TMHMM, Phobius, TOPCONS, or other topology predictors. If tools are unavailable locally, design the pipeline so topology predictions can be imported from user-supplied TSV files.
* BLASTp or DIAMOND for similarity searches if needed.
* Biopython, pandas, numpy, scipy, scikit-learn, matplotlib, seaborn optional, pyhmmer optional.
* Use matplotlib for plots unless there is a good reason to use another plotting library.

## Important scientific constraints

* Do not train on generic “annotated holin” labels as if they are experimentally confirmed.
* Separate experimentally characterized holins from database-annotated or putative holins.
* Avoid circularity: if the model is trained on database annotations and tested on database annotations, that is not meaningful.
* Use hard negative controls, not just random soluble proteins.
* Avoid random sequence-level train/test splits that put close homologs in both train and test.
* Use cluster-based validation, family-level holdout, topology-level holdout, or leave-one-family-out testing.
* Report failures clearly. A negative result is scientifically useful.
* Do not claim “universal motif found” unless the evidence is extremely strong. Prefer “universal sequence HMM performs poorly” if that is what the data show.
* Treat “small hydrophobic protein” as a weak feature, not proof of holin function.
* Treat “near endolysin” as supportive context, not proof.
* Make the final model a ranking/scoring system, not a magic binary classifier.

## Project title idea

“Benchmarking universal, topology-specific, and context-aware signatures for bacteriophage holin prediction”

## Main hypothesis to test

A single universal holin profile HMM will perform poorly across diverse characterized holins, while family-specific HMMs combined with transmembrane topology and lysis-cassette context will provide better candidate ranking.

## Alternative hypotheses

1. A universal HMM detects many known holins but also has many false positives among unrelated small membrane proteins.
2. Topology-specific HMMs improve sensitivity or specificity modestly.
3. Family-specific HMMs perform best for known families but fail on novel holin-like proteins.
4. Adding genomic context improves interpretability and candidate prioritization.
5. No universal linear motif exists, but architecture-level features may be enriched.

## Data categories

Create a clear data model for four groups.

### A. Gold positives

Experimentally characterized holins only. These are proteins where literature supports holin function through at least one of:

* Required for phage lysis.
* Complements a holin-defective system.
* Causes timed membrane permeabilization or depolarization.
* Enables endolysin-dependent lysis.
* Mutations alter lysis timing.
* Direct experimental evidence identifies the protein as holin, pinholin, or holin-like lysis protein.

Required fields for gold positives: protein_id, protein_name, phage_name, host, accession, sequence, protein_length, evidence_type, citation, notes, holin_type (canonical_holin, pinholin, class_I, class_II, class_III, unknown), known_TMD_count, endolysin_partner, spanin_partner, family_label, source_database, manually_verified.

### B. Weak positives

Proteins annotated as holin, putative holin, phage holin family, holin-like protein, pinholin, etc., but without direct experimental evidence. These may be used for exploratory searches but not as gold-standard training labels.

### C. Hard negatives

Proteins that resemble holins superficially but are not known holins. Include: small phage proteins with 1–4 TMDs not near lysis genes; phage membrane proteins with alternative annotations; spanins (especially inner membrane spanin components when separable); tail membrane proteins; small bacterial membrane proteins; toxin-antitoxin membrane toxins; transporter fragments; hypothetical membrane proteins from phages with no lysis context; proteins with enzymatic or structural domains inconsistent with holin function.

### D. Unknown candidates

Hypothetical phage proteins with holin-like features that could be ranked by the final model.

## Input file formats

Create example input files in `data/example/`.

1. `gold_holins.csv` — protein_id,protein_name,phage_name,host,accession,sequence,evidence_type,citation,holin_type,family_label,notes
2. `weak_annotated_holins.csv` — protein_id,protein_name,phage_name,accession,sequence,annotation,source_database,notes
3. `hard_negatives.csv` — protein_id,protein_name,source,accession,sequence,negative_type,annotation,notes
4. `proteins.faa` — FASTA file containing proteins to scan.
5. `genome_context.gff` or `context.tsv` — simple TSV: contig_id,gene_id,start,end,strand,protein_id,product,sequence
6. `topology_predictions.tsv` — protein_id,length,tmd_count,topology,n_region,c_region,signal_peptide,sar_like,tool
7. `known_lisis_genes.tsv` or `lysis_context_terms.txt` — keywords/annotations for: endolysin, lysin, amidase, glycosidase, muramidase, peptidoglycan hydrolase, SAR endolysin, spanin, Rz, Rz1, unimolecular spanin, antiholin, holin, pinholin

## Pipeline stages

### Stage 1: Data validation and cleaning

* Read positive, weak-positive, negative, and unknown sequence files.
* Validate protein sequences; remove invalid sequences; remove exact duplicates; standardize IDs.
* Compute sequence length, hydrophobic fraction, amino acid composition, predicted molecular weight if useful.
* Flag suspicious records: very short proteins under 30 aa; very long proteins over 300 aa; sequences with many ambiguous residues; duplicate sequences appearing in multiple classes; gold positives without citation or evidence type.
* Output cleaned FASTA and metadata TSV files.

### Stage 2: Redundancy reduction and clustering

* Cluster all gold positive holins at multiple identity thresholds: 90%, 70%, 50%, 30%.
* Use MMseqs2 or CD-HIT. Output cluster assignments.
* Summarize whether characterized holins form one connected sequence space or many isolated islands.
* Generate plots: cluster size distribution; sequence length distribution by cluster; topology/TMD count by cluster if topology data are available; network or matrix showing sequence similarity among gold positives if feasible.

### Stage 3: Multiple sequence alignment

Build MSAs for: (A) all gold positives together; (B) topology-specific subsets (1, 2, 3, 4 TMD, unknown); (C) family/cluster-specific subsets (clusters with at least 3 or 5 proteins); (D) optional known holin family groups if family labels are provided.

Important: the universal alignment may be poor. Quantify this instead of hiding it. Compute alignment quality metrics if possible: average pairwise identity, percentage gaps, alignment length, conserved columns, entropy per position. Output alignment reports. If universal alignment is obviously poor, still build the HMM, but report this caveat.

### Stage 4: HMM building

* Model A: `universal_holin.hmm` from all gold positives.
* Model B: topology-specific HMMs: `holin_1TMD.hmm`, `holin_2TMD.hmm`, `holin_3TMD.hmm`, `holin_4TMD.hmm`.
* Model C: family/cluster-specific HMMs: one HMM per cluster/family when enough sequences are available.
* Also produce: combined HMM database; hmmpress output; model metadata table with model_id,model_type,num_sequences,alignment_length,mean_length,mean_pairwise_identity,notes.

### Stage 5: HMM scanning

Run hmmsearch against: gold positives, weak positives, hard negatives, unknown candidate proteins, optional phage proteomes. Output full HMMER results (target protein, model, bit score, E-value, bias, alignment coordinates, coverage, sequence length, model length, class label, dataset category). Compute derived metrics: query coverage, target coverage, best hit per protein, number of models hit per protein, whether top model matches known family if known.

### Stage 6: Transmembrane topology and architecture features

For every protein compute or import: length, TMD count, predicted topology, N-terminal orientation, C-terminal orientation, hydrophobic fraction, number of charged residues, positive/negative charge count, charge distribution in N-terminal 30 aa and C-terminal 30 aa, predicted signal peptide, possible SAR-like N-terminal signal anchor, low complexity regions, presence/absence of domains from InterPro/Pfam/HHpred if supplied. If topology tools cannot be run automatically, create an import function that accepts `topology_predictions.tsv`.

### Stage 7: Genomic context scoring

If genome context is supplied, score each candidate based on neighborhood. For each protein: find neighboring genes within a window (default ±5 genes or ±5 kb); detect nearby lysis-related annotations (endolysin, lysin, amidase, glycosidase, peptidoglycan hydrolase, SAR endolysin, spanin, Rz, Rz1, antiholin, holin, pinholin); record nearest endolysin distance (genes/bp), nearest spanin distance (genes), same strand, plausible lysis cassette, alternative stronger holin candidate nearby, isolation from lysis genes.

Create a `context_score`. Example scoring: +2 if within ±3 genes of endolysin/lysin; +1 if same strand as endolysin; +1 if near spanin/Rz/Rz1 in Gram-negative-type cassette; +1 if near SAR endolysin and candidate has pinholin-like topology; +1 if no better holin candidate nearby; -2 if protein has strong alternative non-holin annotation; -1 if not near any lysis gene; -1 if too large or multi-domain. Make this configurable in config.yaml.

### Stage 8: Candidate holin score

Create a composite `holin_candidate_score` that combines: best family-specific HMM hit, universal HMM hit, topology-specific HMM hit, protein length compatibility, TMD count compatibility, hydrophobicity, absence of enzymatic/structural domain, lysis-cassette context, synteny/context conservation, penalty for alternative annotation, penalty for being too large/short/multi-domain.

The score should rank candidates, not prove function. Output confidence categories: high_confidence_candidate, medium_confidence_candidate, weak_candidate, unlikely_holin. Include explanations for each score. For each candidate generate a “why this was ranked this way” field.

Example output columns: protein_id, best_hmm_model, best_hmm_evalue, best_hmm_bitscore, universal_hmm_hit, topology_hmm_hit, family_hmm_hit, length, tmd_count, hydrophobic_fraction, near_endolysin, near_spanin, context_score, architecture_score, hmm_score, final_holin_score, confidence_category, explanation.

### Stage 9: Benchmarking

Benchmark performance using gold positives and hard negatives. Compare: (A) Universal HMM only; (B) Topology-specific HMMs; (C) Family-specific HMMs; (D) Architecture-only features; (E) Genomic context-only features; (F) HMM + architecture; (G) HMM + architecture + genomic context.

Metrics: sensitivity/recall, specificity, precision, F1, ROC-AUC if meaningful, PR-AUC (especially if class imbalance exists), confusion matrix, false positives, false negatives, top-k recovery, leave-one-family-out recovery, leave-one-cluster-out recovery.

Validation strategies: (1) Random split, labeled as naive and potentially inflated; (2) Cluster-based split where sequences above a chosen identity threshold do not cross train/test boundary; (3) Leave-one-family-out; (4) Leave-one-topology-class-out if enough data; (5) Test against hard negatives specifically.

Do not let close homologs leak across train/test splits. Report if sample sizes are too small. Report confidence intervals if possible using bootstrapping. Analyze false positives and false negatives.

### Stage 10: Motif analysis

Test whether any universal motif exists, but be conservative. Use alignment conservation/entropy to identify conserved positions in the universal MSA. Run motif discovery tools if available, but do not overinterpret hydrophobic motifs. Compare discovered motifs against hard negatives. Distinguish true conserved sequence motifs from generic transmembrane hydrophobicity. Report whether motifs are family-specific rather than universal. Produce sequence logos for: universal alignment if meaningful, topology-specific alignments, family-specific alignments. If the universal alignment is poor, say so clearly. If sequence logos mainly show hydrophobic residues in TMD regions, interpret that as architecture-level enrichment, not a universal holin motif.

### Stage 11: Synteny visualization

Create simple lysis cassette visualizations when genome context is supplied. For selected candidates: draw gene arrows around candidate ±5 genes; color genes by annotation (candidate holin-like protein, endolysin/lysin, spanin/Rz/Rz1, antiholin, structural genes, hypothetical, other); show strand and relative position; export SVG and PNG; make a summary figure for high-confidence candidates.

### Stage 12: Manuscript-style report

Generate a Markdown or Quarto report that includes: (1) Background and rationale; (2) Dataset summary; (3) Gold positive curation summary; (4) Sequence diversity and clustering results; (5) Universal HMM performance; (6) Topology-specific HMM performance; (7) Family-specific HMM performance; (8) Architecture/context-aware scoring performance; (9) Motif analysis; (10) False-positive and false-negative analysis; (11) Candidate ranking table; (12) Limitations; (13) Recommended experimental follow-up.

The report must be cautious and analytical. Suggested conclusion templates:

* If universal HMM performs poorly: “Experimentally characterized holins did not support a robust universal sequence-level HMM. The universal model primarily captured generic hydrophobic membrane-protein features and showed limited specificity against hard negative small membrane proteins.”
* If family-specific HMMs perform better: “Family-specific HMMs improved recovery of related holins but had limited ability to identify distant holin families, consistent with extensive sequence diversity.”
* If context improves performance: “Combining HMM evidence with transmembrane topology and lysis-cassette context improved candidate prioritization compared with sequence-only models.”
* If no universal motif appears: “No universal linear amino-acid motif was detected across the curated holin set. Conserved features were mostly architecture-level, including small size, hydrophobicity, and predicted transmembrane segments.”

Potential manuscript title: “Benchmarking universal and context-aware signatures for bacteriophage holin prediction”

Potential abstract skeleton: “Bacteriophage holins are membrane proteins that control host lysis, but their extreme sequence diversity complicates annotation. We curated experimentally characterized holins and benchmarked universal, topology-specific, and family-specific profile HMMs against hard negative sets of small membrane proteins. We found that [RESULT]. Universal HMMs [RESULT], whereas [RESULT]. Integration of transmembrane topology and lysis-cassette context [RESULT]. These findings suggest that holin annotation should rely on family-specific profiles and genomic context rather than a single universal sequence motif.”

## Stage 0: Literature mining for experimentally characterized holins

Add a dedicated literature-mining and evidence-curation stage before building the holin HMM benchmark. This stage actively scans the scientific literature to identify experimentally characterized bacteriophage holins, pinholins, and holin-like lysis proteins, to build a conservative gold-standard dataset of experimentally validated holins, not merely proteins annotated as “holin” in databases.

Distinguish between: (1) experimentally_validated_holin; (2) experimentally_validated_pinholin; (3) experimentally_validated_holin_like_lysis_protein; (4) literature_claimed_holin_but_no_direct_evidence; (5) database_annotated_only; (6) putative_holin; (7) insufficient_evidence; (8) rejected_non_holin.

**Important scientific rule:** Do not treat the word “holin” in a paper as proof of experimental validation. A protein should enter the gold-positive set only if the paper presents direct experimental evidence that the protein functions in lysis, membrane permeabilization, lysis timing, or endolysin-dependent lysis.

Literature sources: PubMed via NCBI E-utilities; PMC via NCBI E-utilities; Europe PMC search API; Europe PMC open-access full-text XML; user-supplied PMID list; user-supplied DOI list; user-supplied PDF/full-text files; optional manually curated seed papers.

Search strategy — multiple query groups (general holin characterization; pinholin/SAR systems; classic model systems such as lambda S / phage 21 pinholin / T4 / T7 / P2 / N4 / mycobacteriophage; evidence-specific searches such as holin mutant delayed lysis, complementation, membrane depolarization). The module retrieves for each paper: title, authors, year, journal, PMID, PMCID, DOI, abstract, full text if OA, source database, URL, query that found it.

Text mining — for each paper, search title, abstract, figure legends, results, methods, discussion for candidate evidence. Extract candidate protein/gene names (holin, pinholin, antiholin, lysis protein, gene S, S105, S107, Rz, Rz1, lysin, endolysin, SAR endolysin, gpN, ORFn) and experimental-evidence sentences via trigger phrases (required/essential/sufficient for lysis, complementation, restored/failed/delayed/rapid lysis, lysis timing, amber/nonsense/knockout/deletion mutant, membrane depolarization/permeabilization, endolysin dependent, chloroform rescue, lysis inhibition / LIN, antiholin activity).

For every candidate output: candidate_id, paper_id, PMID, PMCID, DOI, title, year, phage_name, host_name, gene_name, protein_name, accession, sequence, evidence_sentence, evidence_section, evidence_type, confidence_class, reason_for_classification, requires_manual_review, manually_verified, curator_notes.

Evidence scoring: 5 = direct genetic and functional evidence; 4 = direct functional lysis or membrane evidence; 3 = strong experimental evidence but sequence mapping uncertain; 2 = literature claim with limited/indirect evidence; 1 = annotation/gene-neighborhood/TMD prediction only; 0 = rejected or irrelevant. Only evidence_score 4–5 and manually_verified = true should be exported as gold positives by default.

Outputs: literature_search_results.tsv; candidate_evidence_sentences.tsv; candidate_holin_literature_table.tsv; sequence_mapping_table.tsv; manual_review_template.tsv; accepted_gold_holins.tsv; rejected_or_weak_holins.tsv; curation_summary.md.

CLI: `holinbench literature-search`, `extract-evidence`, `map-sequences`, `prepare-review`, `export-curated`.

QC checks: flag review vs primary research; flag candidates with no sequence/accession; flag annotation-only support; flag duplicate proteins across papers; flag conflicting names for the same accession; flag papers where “holin” appears only in introduction/background; flag review-citation-only mentions.

**Scientific warning (must be in documentation):** This literature-mining stage is a curation aid, not an authority. Automated extraction can miss evidence, misread papers, or confuse annotation with experimental validation. The gold-positive set must be manually checked before use in HMM training or benchmarking.

## Code structure

```
holin-signature-benchmark/
  README.md  environment.yml  pyproject.toml
  config/config.yaml
  data/example/{gold_holins.csv,weak_annotated_holins.csv,hard_negatives.csv,proteins.faa,context.tsv,topology_predictions.tsv}
  results/.gitkeep
  workflow/Snakefile (optional)
  src/holinbench/{__init__,cli,io,validate,features,clustering,alignment,hmmer,topology,context,scoring,benchmark,motif,synteny,plots,report}.py
  tests/{test_io,test_features,test_scoring,test_context}.py
  docs/{method_notes,data_curation_guidelines,interpretation_guide}.md
  notebooks/exploratory_analysis.ipynb (optional)
```

## Command-line interface

`holinbench {validate,cluster,align,build-hmms,scan,features,context,score,benchmark,motifs,synteny,report,run-all}` plus Stage 0: `{literature-search,extract-evidence,map-sequences,prepare-review,export-curated}`, each `--config config/config.yaml`.

## Configuration

`config/config.yaml` holds: input file paths; output directory; clustering thresholds; minimum sequences per HMM; HMMER E-value cutoffs; bit score cutoffs; protein length range; TMD count range; gene neighborhood window; scoring weights; validation settings; plotting settings; full Stage 0 literature block. Scoring must be configurable, not hardcoded.

## Plots

1. Protein length distribution by dataset category. 2. TMD count distribution by category. 3. Hydrophobic fraction distribution. 4. Gold holin similarity/clustering plot. 5. HMM model performance comparison. 6. ROC/PR curves if appropriate. 7. Confusion matrices. 8. FP/FN summary. 9. Sequence logos. 10. Lysis cassette/synteny maps. 11. Final candidate score distribution.

## Tables

1. Cleaned gold holin dataset. 2. Cluster assignments. 3. HMM model metadata. 4. HMM scan results. 5. Feature table. 6. Context table. 7. Benchmark metrics. 8. Candidate ranking table. 9. False positives/negatives. 10. Motif/conservation summary.

## Documentation

* `README.md`: what the project tests; why a universal holin motif may not exist; installation; required external tools; example run; input formats; output explanation; interpretation warnings; how to add new experimentally characterized holins; how to run benchmarking; how to cite/use output.
* `docs/data_curation_guidelines.md`: what counts as a gold positive; what does not; how to record experimental evidence; how to avoid circular annotation; how to classify weak positives and hard negatives.
* `docs/interpretation_guide.md`: a high score does not prove holin function; HMM hits may reflect family similarity, not universal function; TMDs and hydrophobicity are not sufficient; genomic context is supportive but not definitive; experimental validation is required.

## Testing

Tests check: FASTA and CSV parsing; duplicate removal; feature calculations; context scoring from toy gene neighborhoods; candidate score calculation; benchmark metric calculation; that proteins appearing in both positive and negative sets are flagged; that missing citations in gold positives are flagged.

## Example toy data

Include a very small fake dataset for testing only. Example sequences are placeholders and not for scientific interpretation unless real curated data are supplied.

## Development strategy

1. Scaffold the repository and implement data parsing, validation, feature calculation, and scoring using toy/example data.
2. Implement HMMER wrappers and make them fail gracefully if HMMER is not installed.
3. Implement clustering and alignment wrappers.
4. Implement benchmark metrics.
5. Implement genomic context scoring.
6. Implement plotting and reporting.
7. Add tests and documentation.

Be analytical throughout. Do not write promotional or overconfident claims. This project should be designed to reveal whether the universal motif idea works or fails.
