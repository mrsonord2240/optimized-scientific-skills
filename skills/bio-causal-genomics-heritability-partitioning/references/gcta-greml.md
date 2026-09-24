# GCTA-GREML (reference for heritability-partitioning SKILL.md)

## GCTA-GREML Pipeline

**Goal:** Individual-level SNP-heritability via GRM-based REML (Yang 2011 AJHG 88:76). Use when you
have individual-level genotypes (PLINK bed/bim/fam) and phenotypes, not summary statistics.

**Approach:** Build the genetic relationship matrix (GRM) from autosomal SNPs -> REML
variance-components estimation against the GRM + phenotype.

Run it as `bash examples/gcta_greml.sh <bfile_prefix> <pheno_file> <out_prefix>`: step 1
`gcta64 --autosome --make-grm`, step 2 `gcta64 --reml`. Output `<out_prefix>_h2.hsq` reports
`V(G)/Vp` (the h2 estimate) with SE and a likelihood-ratio p-value against h2 = 0.

Case-control data below 20% prevalence needs `--reml-pcgc` instead of plain `--reml` (PCGC
correction, not otherwise covered here).

## Install

```bash
# GCTA 1.94+ (bioconda, or the binary distribution)
micromamba install -c bioconda -c conda-forge gcta
# or: wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
```

Verified 2026-09-21: `gcta64` 1.94.1 (bioconda) built a GRM from a real 957-individual, 14389-SNP
multi-chromosome genotype panel and recovered h2 = 0.538 (SE 0.135, p = 2.3e-07) from a synthetic
phenotype planted at h2 = 0.5 (200 causal SNPs, standardized genotype effects) -- within 1 SE of the
planted value.

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `N <= 5000` gives very wide SE | GREML SE scales with sample size; not enough power | Report as exploratory; N > 5000 typical for a usable SE |
| REML fails to converge / negative h2 | Population stratification not captured in the GRM, or GRM built from too few SNPs | Include PCs as covariates via `--qcovar`; check GRM diagonal for outliers with `--grm-cutoff` |
