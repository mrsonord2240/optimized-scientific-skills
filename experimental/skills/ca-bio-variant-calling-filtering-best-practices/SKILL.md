---
name: ca-bio-variant-calling-filtering-best-practices
category: Data Analysis
description: Use when filtering a germline or somatic VCF at the site or genotype level and choosing between VQSR, hard filters and ML recalibration.
tool_type: mixed
primary_tool: bcftools
license: MIT
author: GPTomics
---

# Variant Filtering

Pick the one row that matches your task. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your task | Read |
|-----------|------|
| Which method fits my cohort (hard filters, VQSR, VETS, NVScoreVariants, DeepVariant/DRAGEN output) | `routes/choose-method.md` |
| Germline hard filters, SNPs and indels | `routes/hard-filter.md` |
| Set low-quality genotypes to no-call (GQ, DP, allele balance) | `routes/genotype-filter.md` |
| VQSR on a large human cohort | `routes/vqsr.md` |
| Somatic calls from Mutect2 | `routes/somatic.md` |
| The same hard filters in Python | `routes/cyvcf2.md` |
| Check what filtering did (Ti/Tv, FILTER counts) | `routes/validate.md` |
| Depth, allele-frequency, missingness, soft-filter or region-exclusion filters | `routes/cohort-and-regions.md` |
| An error message | `routes/errors.md` |

Rules on every route:

- Site filter first, genotype filter second, cohort QC (missingness, HWE, allele frequency) last.
- SNPs and indels are filtered separately with their own thresholds, then merged.
- A missing annotation passes. MQRankSum and ReadPosRankSum are undefined at hom-alt sites; failing them deletes every hom-alt variant.
- Somatic calls use FilterMutectCalls, never VQSR or the germline hard filters. DeepVariant and DRAGEN output is filtered on the caller's own QUAL and FILTER.
- A VCF stays valid when a filter is wrong. Check the counts and Ti/Tv after every run (`routes/validate.md`).

Tested with bcftools 1.24, cyvcf2 0.34.0, GATK 4.6.1.0.
