# Heritability Partitioning - Usage Guide

## Overview

Estimate SNP heritability (`h2_SNP`) and decompose it across functional categories, tissues, and individual loci using LD-score regression (LDSC), LDAK SumHer, HESS local h2, HDL high-definition likelihood, BOLT-REML, and GCTA-GREML. Choose between summary-statistic methods (LDSC, LDAK, HDL, HESS) and individual-level methods (BOLT-REML, GCTA), reconcile LDSC vs LDAK model-dependent enrichment estimates, prioritize trait-relevant cell types via Finucane 2018 chromatin partitioning, compute trans-ancestry rg with Popcorn, and report calibrated h2 with intercept and ratio diagnostics.

## Prerequisites

Install notes for LDSC, BOLT-LMM, GCTA, and Popcorn (which fork, which flags need a patch, and the
reference-resource downloads) are in SKILL.md's "Tool Install Notes"; LDAK, HESS and HDL installs are in
their `references/` files (SKILL.md, "Per-Method Reference Files").

Inputs: GWAS summary statistics with columns SNP, A1, A2, BETA (or Z), SE, P, N (per-SNP or column-supplied). Allele frequency column EAF strongly recommended. For case-control GWAS, also supply sample case fraction (`--samp-prev`) and population lifetime prevalence (`--pop-prev`).

## Quick Start

Tell your AI agent what you want to do:
- "Compute SNP heritability from this GWAS summary statistic file with LDSC, EUR ancestry"
- "Partition h2 across functional annotations using the baseline-LD v2.2 model"
- "Prioritize trait-relevant tissues from ENCODE chromatin marks via Finucane 2018 cell-type S-LDSC"
- "Compute genetic correlation between trait1 and trait2 from sumstats; use cross-trait LDSC because samples overlap"
- "Run HDL.rg for genetic correlation since cohorts are non-overlapping; want lower variance than LDSC"
- "Compute local heritability per LDetect locus with HESS and identify high-h2 loci for fine-mapping"
- "Reconcile LDSC and LDAK SumHer enrichment estimates for the conserved-region annotation"
- "Convert this case-control LDSC h2 from observed to liability scale; population prevalence is 0.005"
- "Run Popcorn for trans-ancestry rg between this EUR GWAS and the matched EAS GWAS"

## Example Prompts

### Total h2 from EUR Sumstats
> "Munge `gwas_T2D.tsv` with LDSC's `munge_sumstats.py` against the HapMap3 SNP list, then run `ldsc.py --h2` with `eur_w_ld_chr/` for the univariate estimate. Report intercept, mean chi-square, ratio, and h2 with SE on the liability scale (samp-prev 0.08, pop-prev 0.10)."

### Partitioned h2 with Baseline-LD
> "Partition T2D h2 across the baseline-LD v2.2 annotations. Use `--overlap-annot --print-coefficients`. Apply Bonferroni at `0.05 / N_annotations` for per-annotation enrichment claims. Report the top 5 categories with enrichment > 5x and joint p < 2e-3."

### Cell-Type / Tissue Prioritization (Finucane 2018)
> "Prioritize trait-relevant tissues using the Multi_tissue_chromatin_1000Gv3 ldcts file with `ldsc.py --h2-cts`. Apply Bonferroni at `0.05 / N_tissue` (~ 2.5e-4 for 200 tissues). Cross-reference top tissues against the disease's known biology."

### LDSC vs LDAK Reconciliation
> "Functional enrichment claim depends on the per-SNP heritability model. Run BOTH LDSC baseline-LD AND LDAK SumHer with LDAK-Thin tagging. If they disagree by > 2x, report enrichment as model-dependent (the unresolved S-LDSC vs LDAK enrichment debate: Gazal 2019 Nat Genet 51:1202 favours baseline-LD S-LDSC; the LDAK developers defend SumHer, Speed 2020 Nat Genet 52:458)."

### Cross-Trait Genetic Correlation
> "Estimate rg between trait1 and trait2 from sumstats. If sample overlap > 5%, use cross-trait LDSC (`ldsc.py --rg`) because HDL is biased under overlap. If non-overlapping, use HDL for ~60% lower variance. Report rg, SE, p, and cross-trait intercept."

### Local Heritability
> "Compute per-locus h2 across all 22 autosomes using HESS with the LDetect EUR partition (Berisa & Pickrell 2016). Require each locus to have >= 1000 SNPs. Identify loci with h2 > 0.001 for fine-mapping prioritization downstream."

### Liability-Scale Case-Control h2
> "Case-control GWAS for schizophrenia. Supply `--samp-prev <case_fraction>` and `--pop-prev 0.01` so LDSC reports h2 on the liability scale. Without these flags, observed-scale h2 is incomparable across studies."

### Non-EUR Ancestry h2
> "EAS GWAS for type 2 diabetes. Use ancestry-matched LD scores from the EAS folder on alkesgroup.broadinstitute.org/LDSCORE, NOT the default EUR scores. Mismatched LD biases h2 downward and inflates the intercept."

### Trans-Ancestry rg
> "Cross-population genetic correlation between EUR and EAS GWAS for the same trait. Use Popcorn with per-population LD scores. Require effective N > 5000 per population for stable estimation."

### Single-Cell ATAC Annotation
> "Build cell-type-specific .ldcts from per-cluster ATAC peaks (cross-reference atac-seq/single-cell-atac). Compute per-cluster LD scores via `ldsc.py --l2`. Then run `--h2-cts` to identify which cluster's open-chromatin landscape is most heritability-enriched for the trait."

## See Also

Method selection, the Finucane 2018 cell-type workflow, LDSC vs LDAK guidance, computational
footprint, and intercept/ratio interpretation are all in SKILL.md (Decision Tree by Scenario, Cell-Type
Prioritization, LDSC vs LDAK Reconciliation, Computational Footprint, and LDSC Intercept
Interpretation) -- kept in one place so a fix lands once.

## Tips

- Pre-compute LD scores once and reuse; downloading the EUR reference is ~3 GB but covers all standard LDSC analyses.

## Related Skills

causal-genomics/mendelian-randomization - h2 / rg-aware instrument selection and sample-overlap decisions
causal-genomics/colocalization-analysis - Per-locus shared-causal evidence complementary to HESS local h2
causal-genomics/fine-mapping - Credible-set construction at high-h2 HESS loci
causal-genomics/pleiotropy-detection - Cross-trait pleiotropy via LCV / LHC-MR using LDSC outputs
causal-genomics/genomic-sem - Genomic SEM extends LDSC rg to multivariate structural models
causal-genomics/transcriptome-wide-association - TWAS uses partitioned-h2 weights for gene-level testing
atac-seq/differential-accessibility - Per-cell-type chromatin annotations as S-LDSC input
atac-seq/single-cell-atac - scATAC peaks per cluster as .ldcts annotations
chip-seq/peak-calling - ENCODE / Roadmap chromatin marks for cell-type prioritization
population-genetics/association-testing - GWAS source summary statistics for LDSC munging
population-genetics/linkage-disequilibrium - LD reference panels for HESS / coloc.susie
workflows/gwas-pipeline - Upstream GWAS pipeline feeding sumstats to LDSC
