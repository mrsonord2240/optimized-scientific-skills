---
name: ca-bio-microbiome-taxonomy-assignment
category: Data Analysis
description: Use when assigning taxonomy to 16S, ITS or 18S amplicon ASVs against a reference database.
tool_type: mixed
primary_tool: DADA2
license: MIT
author: GPTomics
---

# Taxonomy Assignment

Pick the one row that matches your setup. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your setup | Read |
|------------|------|
| R, a DADA2 `seqtab_nochim.rds`, SILVA-style reference | `routes/dada2.md` |
| R, conservative and novelty-aware calls, or a marker with no pre-trained set | `routes/idtaxa.md` |
| QIIME2, primers known, want a region-matched naive-Bayes classifier | `routes/qiime2-region-nb.md` |
| QIIME2, scikit-learn version error or no retraining wanted | `routes/qiime2-vsearch.md` |
| Which database for 16S, ITS or 18S, or which classifier for a habitat | `routes/choose-database.md` |
| How deep a rank the data supports; setting a confidence threshold | `routes/rank-and-confidence.md` |
| Host or plant sample with Mitochondria or Chloroplast labels | `routes/organelle.md` |
| ASV inference, diversity, differential abundance, shotgun reads | `routes/other-skills.md` |
| An error, or a result that looks wrong | `routes/errors.md` |

Rules on every route:

- The reference must match the marker and the primer region. A full-length-trained classifier on a V4 read both fabricates and erases calls; extract the region first.
- A label is classifier + database + region, not an identification. Confidence is certainty within the database; it cannot detect that the true taxon is absent.
- A short 16S read licenses genus at best. Species only from an exact match (`addSpecies`) or IDTAXA; never lower the threshold to force species. ITS may reach species.
- Seed (`set.seed`) before every `assignTaxonomy`, `IdTaxa` and `LearnTaxa` call.
- State the database release, the training region and the confidence threshold with the labels. Do not merge labels from different databases (SILVA, GTDB) without a crosswalk.
- Keep unassigned ranks as NA; do not drop or force-fill them.

Tested with DADA2 1.30+, DECIPHER 2.30+ (script checked on 3.2.0), QIIME2 2024.10+ (q2-feature-classifier), vsearch 2.22+. Database releases drift: SILVA 138.1/138.2, GTDB r220, Greengenes2 2024.09, UNITE, PR2 5.x. A DADA2 training FASTA, a DECIPHER `trainingSet` and a QIIME2 `.qza` classifier are not interchangeable.
