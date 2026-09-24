---
name: bio-differential-expression-deseq2-basics
description: Performs differential expression on bulk RNA-seq or single-cell pseudobulk counts (aggregated per sample x cell type) with DESeq2's negative-binomial GLM, Wald and LRT testing, apeglm/ashr/normal LFC shrinkage, independent filtering, Cook's outlier handling, VST/rlog transforms, and design formulas including paired, batch, and interaction terms. Use when running bulk or pseudobulk DE, choosing DESeq2 over edgeR or limma-voom, building a paired or interaction design, applying LFC shrinkage for ranking or GSEA, choosing Wald vs LRT, troubleshooting padj=NA, picking VST vs rlog, importing salmon/kallisto via tximport, or analyzing prokaryotic RNA-seq.
tool_type: r
primary_tool: DESeq2
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, apeglm 1.28+, ashr 2.2+, IHW 1.34+, tximport 1.30+, edgeR 4.0+ (for cross-comparison), PyDESeq2 0.5+
Re-checked 2026-09-21 on DESeq2 1.46.0, apeglm 1.28.0, ashr 2.2.63, PyDESeq2 0.5.4 (R 4.4.3).

Install: `BiocManager::install(c('DESeq2', 'apeglm', 'ashr', 'IHW', 'tximport'))` (`BiocManager` from CRAN first if absent); Python: `pip install pydeseq2`.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# DESeq2 Basics

**"Find differentially expressed genes between conditions"** -> Fit a negative-binomial GLM per gene with shared dispersion shrinkage, test the coefficient of interest (Wald) or the joint effect of a factor (LRT), and report a shrunken effect-size estimate for ranking.

## The Single Most Important Modern Insight -- Shrunken LFC and the Wald p-value come from different models

`lfcShrink()` returns LFCs from a Bayesian posterior with apeglm/ashr/normal priors, BUT its `pvalue` column is still the **unshrunken Wald p-value** from `results()`. This is a deliberate design choice (Zhu, Ibrahim, Love 2019 *Bioinformatics* 35:2084) -- the shrunken estimate is for ranking and visualization; the p-value is for inference. Reporting "shrunken LFC = 0.4, padj = 1e-8" mixes two models, which is fine because both are correct for their stated purpose. What is NOT fine: using the shrunken LFC in a downstream filter and then claiming FDR control on that filter (it has none). For threshold-based FDR claims, use `lfcThreshold=` or TREAT (`glmTreat` in edgeR).

`padj` is the exception: `lfcShrink()` recomputes it by re-running independent filtering at its own default `alpha=0.1`, so it differs from `results(alpha = 0.05)` (measured: 1,669 of 1,669 non-NA genes differed, max |diff| 0.154; `metadata(shrunk)$alpha` is 0.1). Pass `res = res` to `lfcShrink()` to carry the `padj` from `results()` through unchanged (verified: identical), or take `padj` from `results()` and the LFC from `lfcShrink()`.

A second consequence: `results(dds)` with no `name=` or `contrast=` argument silently returns the **last coefficient in `resultsNames(dds)`** -- which depends on factor level order and design formula order. Always specify the contrast explicitly. Tutorials that hard-code `results(dds)` are setting an example that breaks the moment another factor is added.

## Algorithmic Taxonomy

| Test / estimator | What it tests | When mandatory | Failure mode |
|------------------|---------------|----------------|--------------|
| Wald | One coefficient = 0 | Two-level factor or single contrast | Anti-conservative with many low-count outliers |
| LRT (`test='LRT'`, `reduced=`) | Joint effect of dropped terms (>=1 df) | Multi-level factor, omnibus, interaction with >1 df | Reports LFC of the LAST coefficient, not omnibus -- read p-value but never report the LFC as "the effect" |
| `lfcShrink(type='apeglm')` (Zhu 2019) | Posterior LFC under heavy-tailed Cauchy prior | DEFAULT for ranking and visualization | Requires `coef=`; cannot use `contrast=` or numeric vectors |
| `lfcShrink(type='ashr')` (Stephens 2017) | Posterior LFC under unimodal prior; reports `lfsr`/`svalue` with `svalue=TRUE` | Arbitrary contrasts via `contrast=` | Slightly different inferential frame (sign-error rather than null-FDR) |
| `lfcShrink(type='normal')` | Posterior LFC under zero-centered normal; accepts `coef=` or `contrast=` (pass `res=` for numeric/list contrasts) | Quasi-deprecated since v1.16; the `pvalue` it returns is still the unshrunken Wald p-value (shrunken p-values need `DESeq(betaPrior = TRUE)`) | Errors on formulas containing interaction terms ("not implemented for designs with interactions") |
| TREAT / `lfcThreshold=` (McCarthy & Smyth 2009) | LFC magnitude exceeds threshold tau | Want FDR control for "|LFC| > 1.5x" claims | Conservative; use only when threshold is biologically pre-specified |

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Two-group bulk RNA-seq, n>=3/group | `DESeq()` + `results(name=...)` + `lfcShrink(coef=..., type='apeglm')` | Modern default; apeglm is the right prior (`references/lfc-shrinkage.md`) |
| Factor with 3+ levels, "any change" question | `DESeq(test='LRT', reduced=~1)`; read padj only | Wald + p-value combining is wrong (`references/lrt.md`) |
| Interaction `~ genotype * treatment` | Build combined factor `group = paste(genotype, treatment)` and design `~ 0 + group`; contrast pairs of interest | Avoids the resultsNames trap; works with apeglm via relevel (`references/interaction-designs.md`) |
| Paired design (tumor/normal same patient) | `~ patient + tissue`; pairing variable FIRST | Absorbs subject variability; n_paired effective sample size |
| Salmon/kallisto input | `tximport()` -> `DESeqDataSetFromTximport()` | Carries length offsets automatically; `DESeqDataSetFromMatrix(round(...))` loses length correction (`references/tximport.md`) |
| n=2/group, no choice | Continue but report results as exploratory; consider edgeR QL F-test as sensitivity | Schurch 2016 *RNA* 22:839: all tools miss 20-40% of true positives at n=3 |
| Single-cell pseudobulk (counts aggregated per donor) | DESeq2 standard pipeline on pseudobulk matrix | Crowell 2020 *Nat Commun* 11:6077: pseudobulk avoids the FDR inflation of cell-level DE |
| Many DE genes expected (>50% of genome) | `estimateSizeFactors(controlGenes=stable)` or spike-in normalization | Median-of-ratios assumes most genes unchanged (`references/size-factors.md`, `references/prokaryotic-rna-seq.md`) |
| GSEA preranked input | Shrunken LFC OR `stat` (Wald Z) as the rank | Unshrunken LFC dominated by low-count noise |
| Cross-sample heatmap, PCA, ML feature | `vst(dds, blind=FALSE)` (or `rlog` if n<30 and library sizes vary >4x) | Raw counts make PC1 = library size (`references/transforms.md`) |

## Reference Files

Read one only when the request needs it.

| File | Read when |
|------|-----------|
| `references/tximport.md` | Salmon / kallisto / RSEM input, or fractional counts |
| `references/interaction-designs.md` | `~ a * b` designs, "effect of B within level of A", combined-factor `~ 0 + group` |
| `references/lrt.md` | Multi-level factor, omnibus test, or an LRT LFC that looks wrong |
| `references/lfc-shrinkage.md` | Choosing apeglm / ashr / normal, or apeglm rejecting a contrast |
| `references/padj-na-remedies.md` | A gene of interest has `padj = NA` (Cook's outlier or independent filtering) |
| `references/transforms.md` | VST vs rlog for PCA, heatmaps, ML features |
| `references/size-factors.md` | Many zeros, spike-ins, or majority-DE biology where median-of-ratios fails |
| `references/dispersion-trend.md` | `plotDispEsts()` shows the fitted trend missing the cloud |
| `references/prokaryotic-rna-seq.md` | Bacterial / archaeal RNA-seq |
| `references/pydeseq2.md` | Python-only environment |
| `references/betaprior-history.md` | Old tutorials using `betaPrior` or contrast= vs name= differences |

## Standard Workflow

**Goal:** Take a raw integer count matrix and a sample table to a ranked, shrunken DE result table.

**Approach:** Construct DESeqDataSet with the design formula, set reference levels explicitly, run the pipeline, extract by explicit contrast, shrink for downstream use.

```r
library(DESeq2)
library(apeglm)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = 'control')

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name = 'condition_treated_vs_control', alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = 'condition_treated_vs_control', res = res, type = 'apeglm')  # res= keeps padj from results()

summary(res)
sig <- subset(res, padj < 0.05)
```

The reference level fix is non-cosmetic: DESeq2 picks alphabetically if not told otherwise, so `c('Treated','Untreated')` makes 'Treated' the reference and the LFC reads inverted. Set it BEFORE `DESeq()`.

`rownames(coldata)` must equal `colnames(counts)` in the same order, or `DESeqDataSetFromMatrix()` errors.

## What to Report

State the named coefficient or contrast, reference level, design formula, FDR alpha, pre-filter threshold, and the number of samples per condition. Report significant-gene and `padj=NA` counts separately; label the NA cause (all-zero, Cook's, or independent filtering). For effect sizes, identify whether the table uses the named Wald estimate or the shrunken LFC and which shrinkage prior supplied it.

## Design Formulas and the resultsNames Trap

**Goal:** Encode batch, paired, and interaction structure correctly and extract the intended contrast.

**Approach:** Put the variable of interest LAST for readability, but never trust the default `results(dds)` -- inspect `resultsNames(dds)` and pass `name=` or `contrast=` explicitly.

```r
design(dds) <- ~ batch + condition
dds <- DESeq(dds)
resultsNames(dds)
# "Intercept" "batch_B_vs_A" "condition_treated_vs_control"

res <- results(dds, name = 'condition_treated_vs_control')
```

Interaction designs (the `genotype:treatment` coefficients, the resultsNames trap and the combined-factor alternative): see `references/interaction-designs.md`.

## Independent Filtering, Cook's Outliers, padj=NA

`padj = NA` has three distinct causes (independent filtering, Cook's outlier, all-zero across every sample), each with a different remediation -- see `references/padj-na-remedies.md` for the Cook's and independent-filtering fixes, and `de-results` for the full diagnostic table, IHW alternative, and recovery code.

Two DESeq2-specific points worth knowing at this layer:

- A gene that is all-zero in ONE group but expressed in the other is testable and does NOT get `padj = NA`: 0,0,0,0,300,300,300,300 gave LFC 10.88, padj 5e-22 (DESeq2 1.46.0). Only `baseMean = 0` (all-zero across every sample) is NA; expect huge LFCs and check such genes by eye.
- Cook's distance filtering is NOT computed for continuous covariates -- a continuous-covariate analysis has effectively no automatic outlier filtering. Disable Cook's only when the outlier IS the signal (`results(dds, cooksCutoff = FALSE)`). At n>=7 per group, `DESeq()` REPLACES outliers via `replaceOutliers()` and refits (`minReplicatesForReplace = 7`).
- Pre-filtering (`rowSums(counts(dds)) >= 10`) is for memory and speed ONLY. It does NOT replace independent filtering, which operates downstream at `results()` time. Independent filtering can be swapped for IHW via `results(dds, filterFun = ihw)`.

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| `design matrix not full rank` | Confounded covariates | `alias(model.matrix(design, coldata))$Complete` to find the redundant column |
| `counts matrix should be integers` | Salmon/kallisto/RSEM counts are fractional | Use `DESeqDataSetFromTximport()` -- it rounds AND carries length offsets (`references/tximport.md`) |
| Wrong sign of LFC vs expected | Reference level set alphabetically | `relevel()` BEFORE `DESeq()` |
| `padj = NA` for biologically meaningful gene | Independent filtering or Cook's outlier | See padj=NA section |
| LRT LFC doesn't match Wald LFC for the same comparison | LRT reports last coefficient; Wald reports the named coefficient | Extract specific Wald per level for the effect size |
| `summary()` shows a different DE count than the `results(alpha = 0.05)` you ran | `summary(res)` uses the alpha stored by `results()`, but an `lfcShrink()` object stores the default `alpha=0.1` unless `res = res` was passed | Pass `res = res` to `lfcShrink()`, or `summary(x, alpha = 0.05)` |

## References

- Anders S, Huber W. 2010. Differential expression analysis for sequence count data. *Genome Biol* 11(10):R106. doi:10.1186/gb-2010-11-10-r106
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15(12):550. doi:10.1186/s13059-014-0550-8
- Zhu A, Ibrahim JG, Love MI. 2019. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35(12):2084-2092. doi:10.1093/bioinformatics/bty895
- Stephens M. 2017. False discovery rates: a new deal. *Biostatistics* 18(2):275-294. doi:10.1093/biostatistics/kxw041
- Ignatiadis N, Klaus B, Zaugg JB, Huber W. 2016. Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nat Methods* 13(7):577-580. doi:10.1038/nmeth.3885
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. *PNAS* 107(21):9546-9551. doi:10.1073/pnas.0914005107
- McCarthy DJ, Smyth GK. 2009. Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics* 25(6):765-771. doi:10.1093/bioinformatics/btp053
- Soneson C, Love MI, Robinson MD. 2015. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Res* 4:1521. doi:10.12688/f1000research.7563.2
- Schurch NJ et al. 2016. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22(6):839-851. doi:10.1261/rna.053959.115
- Crowell HL et al. 2020. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nat Commun* 11:6077. doi:10.1038/s41467-020-19894-4
- Jiang L et al. 2011. Synthetic spike-in standards for RNA-seq experiments. *Genome Res* 21(9):1543-1551. doi:10.1101/gr.121095.111

## Related Skills

- edger-basics - Cross-check or use when n<5/group; QL F-test framework
- de-results - padj=NA handling, IHW, TREAT, GSEA input preparation, gene annotation
- de-visualization - MA, volcano (with shrunken LFC), PCA, heatmap, dispersion plot
- batch-correction - Include batch in design vs Nygaard 2016 cardinal sin
- timeseries-de - DESeq2 LRT with splines for time-course
- expression-matrix/counts-ingest - tximport, featureCounts, STAR output decisions
- expression-matrix/normalization - RLE/TMM/VST/rlog mechanics and failure modes
- expression-matrix/metadata-joins - Reference level, paired design, interaction parameterization
- expression-matrix/gene-id-mapping - Annotating DE results with symbols
- rna-quantification/tximport-workflow - Detailed tximport mechanics
- pathway-analysis/gsea - Ranked-list input from DE
- pathway-analysis/go-enrichment - ORA with proper background
- data-visualization/volcano-and-ma-plots - Custom volcano with apeglm-shrunken LFC
