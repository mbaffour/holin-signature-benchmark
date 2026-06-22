# Data curation guidelines

The scientific value of this project depends entirely on the quality of the **gold positive** set.
Automated literature mining (Stage 0) is a *curation aid*, not an authority. Every candidate must
be checked by a human before it is used to train or benchmark anything.

## What counts as a gold positive

A protein is a **gold positive** only if the primary literature provides **direct experimental
evidence** that the protein itself functions in lysis. At least one of the following must be
demonstrated *for that protein* (not merely asserted, and not for a homolog):

- A loss-of-function mutation (amber/nonsense, deletion, knockout) **disrupts or retimes lysis**.
- **Complementation** by the protein restores lysis to a defective system.
- Expression of the protein causes **membrane depolarization or permeabilization**.
- Expression of the protein **plus an endolysin** causes lysis (endolysin-dependent lysis).
- A mutation in the protein **changes lysis timing**.
- The protein is experimentally shown to act as a **pinholin** with a SAR endolysin.
- The protein is directly tested in a **lysis assay, membrane assay, or phage infection assay**.

Record the evidence with `evidence_type` (which assay) and `citation` (PMID/DOI). Stage 1 will
flag any gold record missing these fields.

## What does NOT count as a gold positive

- The paper calls the protein a "holin" based only on **annotation**.
- The gene is **near an endolysin** but was never experimentally tested.
- The protein merely has **transmembrane domains** or is "small and hydrophobic".
- A **database product name** says "holin", or **BLAST/Pfam** says "holin-like".
- A **review** discusses it but presents no primary experimental data.
- The paper says **"putative holin"** / **"predicted holin"**.
- No **protein sequence or accession** can be linked with confidence.

These belong in the **weak positives** set, not the gold set.

## How to record experimental evidence

Use the `evidence_score` scale (Stage 0):

| Score | Meaning | Eligible for gold? |
| --- | --- | --- |
| 5 | direct genetic **and** functional evidence | yes (after manual check) |
| 4 | direct functional lysis or membrane evidence | yes (after manual check) |
| 3 | strong experimental evidence, sequence mapping uncertain | no (resolve mapping first) |
| 2 | literature claim, limited/indirect evidence | no → weak positive |
| 1 | annotation / gene-neighborhood / TMD prediction only | no → weak positive |
| 0 | rejected or irrelevant | no |

By default only `evidence_score >= 4` **and** `manually_verified = true` are exported to gold
(`literature.require_manual_verification_for_gold: true`).

## How to avoid circular annotation

- **Never** train on generic "annotated holin" labels as if they were experimentally confirmed.
- **Never** evaluate a model trained on database annotations using database annotations — that
  measures self-consistency, not biology.
- Keep gold (experimental) strictly separate from weak (annotation) throughout.
- For benchmarking, prefer the **leave-one-family-out** regime so a model is never tested on close
  homologs of its own training data. The naive regime is reported but labeled as inflated.

## Classifying weak positives and hard negatives

**Weak positives** — annotated/putative holins without direct experimental evidence. Useful for
exploratory searches; never gold-standard labels.

**Hard negatives** — proteins that *look* like holins but are not known holins. Good hard negatives
make the benchmark meaningful (random soluble proteins are too easy):

- Small phage proteins with 1–4 TMDs **not** near lysis genes.
- Phage membrane proteins with alternative annotations.
- Inner-membrane **spanin** components (when separable from holins).
- **Tail** membrane proteins.
- Small **bacterial** membrane proteins.
- **Toxin–antitoxin** membrane toxins.
- **Transporter** fragments.
- Hypothetical phage membrane proteins with **no lysis context**.
- Proteins with **enzymatic/structural domains** inconsistent with holin function.

## Adding a curated holin

1. Verify the protein meets the gold bar above.
2. Add a row to your gold CSV with `evidence_type` and `citation` filled.
3. Re-run `holinbench validate` and resolve any QC flags.
