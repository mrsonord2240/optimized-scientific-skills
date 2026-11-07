---
name: ca-bio-vcf-statistics
category: Data Analysis
description: Use when judging whether a VCF callset is trustworthy, or screening a cohort VCF for sample swaps, contamination, relatedness and sex mismatches.
tool_type: cli
primary_tool: bcftools
license: MIT
author: GPTomics
---

# VCF Statistics

Pick the one row that matches your task. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your task | Read |
|-----------|------|
| Counts, Ti/Tv, per-sample het/hom, QUAL for a VCF (CLI or Python) | `routes/callset-stats.md` |
| What a number means: low Ti/Tv, outlier het/hom, novel fraction, what to do about it | `routes/interpret.md` |
| Per-sample or per-site missingness, call rate, HWE | `routes/missingness-hwe.md` |
| Suspected contamination (het allele balance) | `routes/contamination.md` |
| Sample swaps, duplicates, relatedness, sex check | `routes/identity.md` |
| Accuracy in difficult regions (homopolymers, repeats, segdups) | `routes/stratified.md` |
| Apply filters, normalize, annotate, or compare call sets | `routes/other-skills.md` |
| An error message | `routes/errors.md` |

Rules on every route:

- No metric has a universal pass value. Compare a sample with a matched cohort (same assay, same inferred ancestry), never with an absolute number.
- Set low-quality genotypes to `./.` (GQ, DP, allele balance) before computing cohort missingness, call rate or HWE.
- Run identity QC (swaps, relatedness, sex) before any association or burden analysis.
- Only `./.` is missing; `0/0` is a confident hom-ref call.

Tested with bcftools 1.19+, vcftools 0.1.16+, somalier 0.2.19+, peddy 0.4.8+, cyvcf2 0.30+.
