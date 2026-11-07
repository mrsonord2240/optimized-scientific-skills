---
name: ca-bio-microbiome-amplicon-processing
category: Data Analysis
description: Use when turning demultiplexed 16S rRNA or ITS amplicon FASTQ reads into an ASV feature table with DADA2.
tool_type: mixed
primary_tool: DADA2
license: MIT
author: GPTomics
---

# Amplicon Processing with DADA2

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data | Read |
|-----------|------|
| Raw paired 16S reads that still carry primers | `routes/remove-primers.md` |
| Primer-trimmed paired 16S reads, one or several runs | `routes/dada2-16s.md` |
| Choosing truncation lengths, or a near-zero merge rate | `routes/truncation.md` |
| Fungal ITS (variable length) | `routes/its.md` |
| Negative controls and a low-biomass sample set, to remove reagent contaminants | `routes/decontam.md` |
| DADA2 or Deblur inside QIIME2, or PacBio CCS or 454/Ion Torrent reads | `routes/qiime2.md` |
| DADA2, Deblur or OTU; which fits the design | `routes/choose-method.md` |
| Taxonomy, diversity, differential abundance, shotgun reads | `routes/other-skills.md` |
| An error, or a result that looks wrong | `routes/errors.md` |

Rules on every route:

- Primers come off first (cutadapt `--discard-untrimmed`), then filter, then learn errors. Leftover primers corrupt the error model and fake chimeras.
- The error model is per run. Never concatenate runs before `learnErrors`; infer each run, `mergeSequenceTables`, then one `removeBimeraDenovo`.
- `maxN = 0`; DADA2 cannot model ambiguous bases.
- Never fixed-truncate ITS.
- An ASV is an exact inferred sequence from one run, not an organism: one genome with several 16S copies can give several ASVs, so ASV count is an upper bound on richness.
- Do not merge a DADA2 table with a Deblur table; the features are defined differently.

Tested with DADA2 1.30+, cutadapt 4.6+, ITSxpress 2.0+, QIIME2 2024.2+. `dada()` defaults (`OMEGA_A` 1e-40, `mergePairs` `minOverlap` 12) and QIIME2 flag spellings drift between releases; confirm with `?dada` and `qiime dada2 --help`.
