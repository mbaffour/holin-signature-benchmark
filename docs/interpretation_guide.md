# Interpretation guide

This pipeline produces **rankings and benchmarks, not biological facts**. Read this before drawing
conclusions from any output.

## A high score does not prove holin function

`final_holin_score` is a weighted combination of HMM, architecture, and context evidence. It is a
prioritization for experiments. A high-confidence candidate is a *good thing to test in the lab*,
not a confirmed holin. Always read the per-candidate `explanation` field — it states exactly which
evidence drove the score.

## HMM hits may reflect family similarity, not universal function

A hit to a **family-specific** HMM means the protein resembles members of that family. It does not
mean a universal holin signature exists. A hit to the **universal** HMM is weak evidence: the
universal model is built from very diverse sequences and tends to capture generic membrane-protein
features. Treat universal-only hits with suspicion.

## TMDs and hydrophobicity are not sufficient

Many unrelated proteins are small, hydrophobic, and have 1–4 TMDs (transporters, tail proteins,
spanins, toxin–antitoxin toxins). That is exactly why the benchmark uses **hard negatives** with
these properties. If a candidate's support is only "small + hydrophobic + has TMDs", it is weak.

## Conserved columns may be hydrophobicity, not a motif

Stage 10 explicitly distinguishes a holin-specific linear motif from generic transmembrane
hydrophobicity. If `hydrophobic_fraction_of_conserved` is high, the "conservation" is
architecture-level (the residues are conserved because they sit in a membrane), **not** a universal
holin motif. Do not report such columns as a discovered motif.

## Genomic context is supportive, not definitive

Proximity to an endolysin or spanin raises plausibility, but lysis cassettes contain many genes,
and holins are not the only membrane proteins near lysis genes. Context is a tie-breaker and a
sanity check, never proof.

## The naive benchmark is circular — trust leave-one-family-out

The HMM models are built from the gold positives. Evaluating them on those same gold positives
(the **naive** regime) inflates performance and is labeled with a `circularity_warning`. The honest
question — *does a universal model generalize to a holin family it has never seen?* — is answered by
the **leave-one-family-out** table (`benchmark_lofo_universal.tsv`). If LOFO recall is low (or LOFO
false-hits on hard negatives are high), the universal model does **not** generalize, regardless of
its naive AUC.

## Small sample sizes

Experimentally validated holins are scarce. Per-family and per-topology metrics can rest on one or
two sequences. The pipeline warns when class sizes are below `benchmark.validation.min_class_size_warn`
and reports bootstrap confidence intervals — wide intervals mean "we don't really know yet".

## Experimental validation is required

No output here substitutes for experiments. To confirm a candidate is a holin, test it directly:
lysis-timing assays, membrane depolarization/permeabilization, complementation of a holin-defective
phage, and endolysin-dependent lysis.
