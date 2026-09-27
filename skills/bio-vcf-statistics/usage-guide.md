# VCF Statistics Usage Guide

## Overview

Compute and interpret VCF quality-control metrics for a callset (Ti/Tv, het/hom, novel/known,
missingness, HWE, contamination) and cohort identity (sample swaps, relatedness, sex). No metric
has a universal pass value -- read each against a matched cohort, not a fixed threshold.

## Prerequisites

- bcftools (`conda install -c bioconda bcftools`) with plot-vcfstats (requires python + matplotlib)
- vcftools (`conda install -c bioconda vcftools`) for per-sample missingness, exact HWE, and KING-robust kinship
- somalier (`conda install -c bioconda somalier`) and peddy (`pip install peddy`) for scalable identity/relatedness/sex/ancestry QC
- cyvcf2 (`pip install cyvcf2`) and numpy for custom per-record statistics in Python

## Quick Start

Tell your AI agent what you want to do:
- "Generate comprehensive statistics for my VCF including Ti/Tv ratio"
- "My exome Ti/Tv is only 2.1 -- is the callset over-permissive?"
- "Flag het/hom outliers in my cohort accounting for ancestry"
- "Decide whether this HWE deviation is a genotyping error or real biology"
- "Screen my cohort for sample swaps, contamination, and wrong-sex annotations"
- "Compare variant statistics before and after filtering"

## Custom Statistics in Python (cyvcf2)

`examples/vcf_stats.py` computes counts, Ti/Tv, and mean QUAL in one pass. Two extensions not in
that script:

### Per-sample genotype distribution

```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
samples = vcf.samples
hom_ref = [0] * len(samples); het = [0] * len(samples)
hom_alt = [0] * len(samples); missing = [0] * len(samples)

for variant in vcf:
    for i, gt in enumerate(variant.gt_types):   # cyvcf2 codes: 0 hom-ref, 1 het, 2 unknown, 3 hom-alt
        if gt == 0:
            hom_ref[i] += 1
        elif gt == 1:
            het[i] += 1
        elif gt == 3:
            hom_alt[i] += 1
        else:
            missing[i] += 1

for i, s in enumerate(samples):
    ratio = het[i] / hom_alt[i] if hom_alt[i] else 0
    print(f'{s}: het/hom={ratio:.2f} HET={het[i]} HOM_ALT={hom_alt[i]} MISS={missing[i]}')
```

### Allele-frequency spectrum

```python
from cyvcf2 import VCF

bins = {'rare(<1%)': 0, 'low(1-5%)': 0, 'common(5-50%)': 0, 'frequent(>50%)': 0}
for variant in VCF('input.vcf.gz'):
    af = variant.INFO.get('AF')
    if af is None:
        continue
    af = af[0] if isinstance(af, tuple) else af
    if af < 0.01:
        bins['rare(<1%)'] += 1
    elif af < 0.05:
        bins['low(1-5%)'] += 1
    elif af < 0.5:
        bins['common(5-50%)'] += 1
    else:
        bins['frequent(>50%)'] += 1
bins
```

The metric decision table, bcftools/vcftools/somalier/peddy/KING commands, HWE and contamination
mechanics, and stratified evaluation are in `SKILL.md`.

## What the Agent Will Do

1. Identify the VCF/BCF and check indexing, then run `bcftools stats` (per-sample, region-based, or comparison mode as needed)
2. Extract key metrics -- counts, Ti/Tv (overall and novel), per-sample het/hom, depth, missingness
3. Contextualize each metric against the assay and, where relevant, inferred ancestry rather than an absolute threshold
4. For cohorts, run identity QC (gtcheck/peddy/somalier/KING) and reconcile relatedness against the manifest
5. Flag issues (low Ti/Tv, ancestry-adjusted het/hom outliers, excess-het sites, contamination signatures, unexpected relatedness) with the recommended action

## Tips

- Ti/Tv ~2.0-2.1 (WGS) and ~3.0-3.3 (WES) indicate good calls; below range signals false positives, since random errors have Ti/Tv ~0.5
- Never apply a single global het/hom cutoff across ancestries; stratify first
- Apply genotype-level filters before computing cohort missingness or HWE, or garbage genotypes drive both
- Filter HWE on excess heterozygosity only, within ancestry, in controls; heterozygote deficit is often real
- The het allele-balance distribution shifting off 0.5 is a contamination signal; confirm with VerifyBamID2/CHARR on the BAM
- Run identity QC (swap/relatedness/sex) early -- one undetected swap can create or erase a significant hit
- plot-vcfstats requires matplotlib for the PNGs and pdflatex or tectonic for `summary.pdf`; without a LaTeX engine it exits 2 at the PDF step but the PNGs are still written

## Related Skills

- variant-calling/filtering-best-practices - Apply the filters these metrics motivate
- variant-calling/gatk-variant-calling - Feed VerifyBamID2 contamination into calling
- variant-calling/variant-normalization - Normalize before comparing or annotating call sets
- variant-calling/vcf-basics - Query and understand VCF fields
- variant-calling/joint-calling - Cohort genotyping where population QC applies
