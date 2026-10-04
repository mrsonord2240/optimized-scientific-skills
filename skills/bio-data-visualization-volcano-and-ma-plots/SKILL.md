---
name: bio-data-visualization-volcano-and-ma-plots
description: Use when visualizing differential-expression results (RNA-seq, ChIP-seq, ATAC-seq, proteomics) or any per-feature effect-size and p-value table as a volcano or MA plot.
tool_type: mixed
primary_tool: ggplot2
license: MIT
author: GPTomics
category: Data Analysis
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+ (re-checked on 1.46.0), EnhancedVolcano 1.20+ (1.24.0), ggplot2 3.5+ (4.0.3), ggrepel 0.9.5+ (0.9.8), matplotlib 3.8+ (3.11.2), numpy 1.26+ (2.5.3), adjustText 1.1+ (1.4.0), apeglm 1.28+, ashr 2.2+ (2.2.63).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. EnhancedVolcano 1.24 still emits ggplot2 `size`-for-lines deprecation warnings under ggplot2 4.x; the plot renders.

# Volcano and MA Plots

**"Plot differential-expression results"** -> Place per-feature shrunken effect estimate on the x-axis and a significance measure (-log10 padj or -log10 p) on the y-axis. The decision space spans which effect estimate (raw vs shrunken), which significance measure (raw p vs adjusted vs s-value), how to encode categories (color by direction, not by gradient), how to label (top-N is rarely informative), and how to handle the tail (extreme p compresses the plot).

- R: `EnhancedVolcano::EnhancedVolcano()`, `ggplot2 + ggrepel`, `DESeq2::plotMA()`
- Python: `matplotlib.scatter` with `adjustText`, `sanbomics.plots.volcano` (`from sanbomics.plots import volcano`; `data` needs a `symbol` column, `pvalue='padj'` by default, one grey class for all DE points, no Up/Down colours), custom seaborn

## The Single Most Important Modern Insight -- Plot Shrunken LFC

Raw log2 fold change from DESeq2 / edgeR is the maximum-likelihood estimate and **inflates wildly at low counts**. A gene with 2 vs 0 reads gets log2FC = Inf; one with 4 vs 1 gets log2FC = 2 with a huge standard error. A naive volcano labels these as "top hits" purely because they have extreme estimates, not because they have a real signal.

`DESeq2::lfcShrink()` applies an empirical-Bayes prior that, on the airway data (baseMean < 5), cut the median |LFC| from 0.76 to 0.01 (apeglm) and 0.02 (ashr). Shrinkage is not strictly monotone: apeglm returns a MAP estimate that exceeded the MLE in absolute value for 108 of 29,391 airway genes (max 9.51 to 10.99; largest single increase 1.48, median baseMean 1444); ashr raised none. The default `type` in DESeq2 1.46.0 is `'apeglm'` (Zhu, Ibrahim, Love 2019 *Bioinformatics* 35:2084) which uses a Cauchy prior — heavy enough to preserve large real effects, sharp enough at zero to deflate noise. Plot the shrunken LFC. The unshrunken LFC is a misleading effect estimate for ranking, labeling, or thresholding.

A complementary modern alternative is the **s-value** (Stephens 2017 *Biostatistics* 18:275): the local false sign rate — probability that the sign of the effect is wrong. s-values rank genes by *how confident the sign is*, which is what a volcano plot is implicitly trying to communicate. Where padj answers "is the effect non-zero," s answers "do we know which direction."

## Shrinkage Method Selection

| Method | Prior | Best for | Fails when |
|--------|-------|----------|------------|
| `apeglm` (Zhu 2019) | Cauchy | Default; preserves large LFCs, deflates noise | Requires `coef=`; no support for `contrast=` |
| `ashr` (Stephens 2017) | Mixture of normals | Comparisons requiring `contrast=` (also accepts `coef=`) | Slightly more aggressive shrinkage of medium effects |
| `normal` (DESeq2 original) | Zero-centered normal | Legacy reproducibility only | Over-shrinks large effects; no `svalue=TRUE` (errors in 1.46.0) |
| Unshrunken MLE | None | NEVER for volcano/MA plots | Low-count genes dominate the tails with no real signal |

```r
library(DESeq2)
dds <- DESeq(dds)
res_apeglm <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
res_ashr <- lfcShrink(dds, contrast = c('condition', 'treated', 'control'), type = 'ashr')
# s-values (Stephens 2017 local false sign rate): add svalue = TRUE to apeglm or ashr.
# Neither returns a svalue column by default, and svalue = TRUE replaces pvalue/padj
# with svalue (columns baseMean, log2FoldChange, lfcSE, svalue), so keep results() for padj.
res_s <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm', svalue = TRUE)
```

For edgeR users: `topTags()` already provides moderated p-values but does not shrink LFC. Use `glmTreat()` for a moderated test against a non-zero LFC threshold; this is the edgeR equivalent of the shrunken-LFC philosophy.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Bulk RNA-seq with DESeq2 | `lfcShrink(type='apeglm')` + plot on padj threshold | apeglm is the default since DESeq2 v1.28 |
| Non-default contrast required | `lfcShrink(type='ashr')` | apeglm requires `coef=`, ashr accepts `contrast=` |
| Single-cell pseudobulk DE | `lfcShrink(type='apeglm')` + optionally rasterize the points | Rasterizing shrank a 17,994-point PDF from 540 KB to 38 KB (ggrastr 1.0.2) |
| Proteomics (limma/MSstats) | Already moderated; plot logFC vs adj.P.Val directly | limma's empirical-Bayes already shrinks |
| Microarray (limma) | `topTable()` adj.P.Val + logFC | Same — limma is the original shrunken-LFC method |
| ATAC/ChIP differential peaks | `lfcShrink(type='apeglm')` on DESeq2/DiffBind | Treat peaks as features identically to genes |
| Want to rank by sign-confidence | `lfcShrink(..., svalue = TRUE)` (apeglm or ashr), NOT padj | Stephens 2017 — sign-aware ranking |
| Many comparisons in one figure | Faceted MA plot with shared y-axis | MA scales better than volcano for >6 panels |

## Volcano with ggplot2 + ggrepel

**Goal:** Plot shrunken LFC vs -log10 adjusted p-value, color by significance class, label genes of interest with non-overlapping repulsion.

**Approach:** Compute a categorical significance variable from padj AND |LFC| thresholds; pre-select labels (genes of interest OR top-N by combined rank) before plotting; use `ggrepel::geom_text_repel` with `max.overlaps = Inf` to guarantee every selected label appears.

Function: `scripts/volcano_plot.R` (`volcano_plot(res, fdr, lfc_threshold, label_genes, top_n, y_cap)`; `source()` it, then call it on the shrunken results).

Key design choices encoded in `scripts/volcano_plot.R`:
- Colors from Okabe-Ito (Wong 2011 *Nat Methods* 8:441) — CVD-safe categorical palette
- `max.overlaps = Inf` because the ggrepel default of 10 silently drops labels with no error
- y is `-log10(padj)`, so the dashed `-log10(fdr)` line is exactly the colour boundary; a raw-p axis with an FDR line cannot match it. Genes with `padj = NA` (independent filtering) are not drawn and the count is reported
- No y cap by default. Extreme values (airway reaches 131) are better left visible; `y_cap` squishes higher genes onto the cap as triangles with a legend entry, so none are hidden (airway, cap 30: 118 triangles, 115 significant)
- `rank_score = -log10(padj) * abs(log2FoldChange)` selects labels that are both significant AND have non-trivial effect; pure top-N-by-p selects high-count genes with tiny effects

## EnhancedVolcano -- Production Use and Its Gotchas

```r
library(EnhancedVolcano)
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',                    # use padj NOT pvalue for the threshold line
    pCutoff = 0.05,
    FCcutoff = 1,
    selectLab = c('TP53', 'MYC', 'BRCA1'),
    drawConnectors = TRUE,
    widthConnectors = 0.3,
    maxoverlapsConnectors = Inf,
    colAlpha = 0.6,
    pointSize = 1.5,
    labSize = 3,
    colCustom = ev_col,
    ylab = bquote(-log[10]~'adjusted'~italic(P)),
    legendPosition = 'right')
```

`ev_col` is one colour per row of `res` (drop `padj = NA` rows first), named by class so the names become the legend entries: `setNames(okabe_ito[sig], sig)` with `sig` = Up/Down/NS from padj and LFC, built as in `scripts/volcano_phd.R`. The `col =` argument cannot colour by direction: its fourth colour is shared by Up and Down.

**Gotcha 1: `y = 'pvalue'` vs `y = 'padj'`.** EnhancedVolcano's default axis label is '-Log10 P' whichever column you pass, so set `ylab` to name the quantity plotted. With `y = 'pvalue'` the `pCutoff` line is a raw-p line, not an FDR; use `y = 'padj'`.

**Gotcha 2: Asymmetric x-limits hide the diverging null distribution.** Use `xlim = c(-max(abs(LFC)), max(abs(LFC)))` so the plot is symmetric around zero.

**`selectLab` behaviour (EnhancedVolcano 1.24.0, checked):** a listed gene that fails `pCutoff`/`FCcutoff` is still labelled, and a name absent from `lab` is silently ignored, so verify spelling and that `lab` holds the same identifiers.

## MA Plot -- The Underused Diagnostic

The MA plot (log2-mean vs log2-fold-change) is the original RNA-seq diagnostic (Dudoit 2002 *JASA*). It exposes the abundance-dependent variance structure that the volcano hides:

- A fan-shaped MA plot with extreme LFCs concentrated at low baseMean indicates **inadequate shrinkage**
- A horizontal "stripe" of significant genes at one LFC value indicates **batch confound with treatment**
- An asymmetric distribution (more Up than Down) at the low-count end indicates **library-size normalization failure**

```r
library(DESeq2)
plotMA(res_apeglm, alpha = 0.05, ylim = c(-5, 5))
# alpha colors significant points; ylim clips for readability without losing the gene
```

Python: `scripts/ma_plot.py` (`ma_plot(res, fdr, ax)`; matplotlib scatter with `rasterized=True`).

`rasterized=True` shrinks the PDF: for all 29,391 airway rows the MA scatter was 23 KB rasterized vs 444 KB vector. Vector scatter at this size stays under 1 MB, so rasterize for file size or editing speed, not out of necessity.

## Resources

- `scripts/volcano_plot.R` — ggplot2 + ggrepel volcano function (y = -log10 padj, optional marked `y_cap`, combined-rank labels).
- `scripts/volcano_phd.R` — end-to-end DESeq2 `lfcShrink` -> volcano -> MA -> `ggsave(cairo_pdf)` -> EnhancedVolcano (`colCustom`); expects a fitted `dds` and `condition` coefficient in the session.
- `scripts/ma_plot.py` — matplotlib MA plot.
- `references/failure-modes.md` — per-method failure modes (trigger/mechanism/symptom/fix) and error table.
- `references/reconciliation-thresholds-pushback.md` — when methods disagree, quantitative thresholds, reviewer responses.
- `usage-guide.md` — prerequisites and example prompts.

**Operational rule:** read a volcano in this order — (1) is LFC shrunken? (2) is the y-axis padj or raw p? (3) does the threshold line match the axis? (4) are labels selected by combined rank? If any is wrong, the plot is misleading.

## References

- Anders S, Huber W. 2010. Differential expression analysis for sequence count data. *Genome Biol* 11:R106.
- Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS-B* 57:289-300.
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. *PNAS* 107:9546-9551.
- Dudoit S, Yang YH, Callow MJ, Speed TP. 2002. Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments. *Stat Sin* 12:111-139.
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550.
- Stephens M. 2017. False discovery rates: a new deal. *Biostatistics* 18(2):275-294. doi:10.1093/biostatistics/kxw041
- Wong B. 2011. Points of view: Color blindness. *Nat Methods* 8(6):441. doi:10.1038/nmeth.1618
- Zhu A, Ibrahim JG, Love MI. 2019. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35(12):2084-2092. doi:10.1093/bioinformatics/bty895

## Related Skills

- differential-expression/de-results - Filter and rank DE result tables before plotting
- differential-expression/deseq2-basics - Run DESeq2 to produce the input results object
- differential-expression/de-visualization - DESeq2 / edgeR built-in plot helpers
- data-visualization/distribution-plots - Boxplot / raincloud follow-up for specific gene panels
- data-visualization/color-palettes - Okabe-Ito and CVD-safe palette selection
- data-visualization/ggplot2-fundamentals - Underlying grammar of graphics
- pathway-analysis/go-enrichment - Functional enrichment from the gene lists produced
