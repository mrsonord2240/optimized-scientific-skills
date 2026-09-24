# Differential Abundance - Usage Guide

## Overview
Identify proteins with significantly different abundance between experimental conditions. The central decisions are how to handle missing values (model the left-censored MNAR dropout, do not impute it), how to moderate per-protein variance at the small sample sizes proteomics uses, and at what level to test (protein summary vs feature/peptide). The skill covers limma, DEqMS, proDA, msqrob2, MSstats, and a Python Welch+BH fallback, plus minimum-fold-change testing and fold-change shrinkage. Method choice, workflows, thresholds and failure modes are in `SKILL.md`; the peptide-level (msqrob2/MSstats) material is in `references/feature_level.md`.

## Prerequisites
```bash
pip install numpy pandas scipy statsmodels
```
```r
BiocManager::install(c("limma", "DEqMS", "proDA", "msqrob2", "QFeatures", "MSstats", "ashr"))
```

## Quick Start
Tell your AI agent what you want to do:
- "Find differentially abundant proteins between treatment and control in my intensity matrix"
- "Run limma with empirical-Bayes moderation on my small-n protein data"
- "Use DEqMS because I have PSM counts per protein from a TMT experiment"
- "My label-free data has 30% missing values and some on/off proteins -- test it without imputing"
- "Test for at least a 1.5-fold change instead of just nonzero, without inflating FDR"
- "Test at the peptide level from my evidence.txt with msqrob2 instead of a protein matrix" (`examples/msqrob2_peptide_level.R`)

## Example Prompts

### Choosing a Method
> "I have 4 replicates per group of protein-level LFQ intensities. Pick and run the right differential test."

> "Analyze my TMT proteomics data for differential abundance. I have PSM counts per protein, so use DEqMS, and remember it is multi-batch."

> "My label-free data has many proteins detected in one group but missing in the other. Test these honestly instead of imputing."

### Missingness and Batch
> "Set up a limma model with condition and batch as covariates -- do not remove the batch effect before testing."

> "Some proteins are missing because they are below the detection limit. Use a method that models this dropout."

### Effect Size and Thresholds
> "Test whether fold changes exceed 1.5-fold using treat and topTreat, not a post-hoc filter."

> "Report raw fold changes for GSEA and shrunken estimates for the figure."

## Related Skills

- quantification - peptide-to-protein summarization, normalization, and IRS that produce the matrix this skill tests
- proteomics-qc - quality control and batch-effect assessment before testing
- protein-inference - razor/shared-peptide ambiguity that drives which protein group gets the quantity
- ptm-analysis - site-level differential testing for modified peptides
- differential-expression/de-results - analogous empirical-Bayes interpretation for RNA-seq DE
- data-visualization/volcano-and-ma-plots - volcano and MA plots of the result table
- pathway-analysis/go-enrichment - functional enrichment of the significant protein hit list
- machine-learning/biomarker-discovery - building predictive panels from differential proteins
- workflows/proteomics-pipeline - end-to-end pipeline that calls this skill as the testing stage
