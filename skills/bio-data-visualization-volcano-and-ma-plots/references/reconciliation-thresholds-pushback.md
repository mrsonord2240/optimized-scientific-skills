# Reconciliation, thresholds, reviewer pushback

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| apeglm and ashr give different "top hits" | Different prior shapes; medium-effect, medium-count genes are most sensitive | Both are valid; pick one and document. apeglm is the DESeq2 default and the published recommendation |
| EnhancedVolcano shows fewer points than ggplot | EnhancedVolcano drops `padj = NA` (DESeq2 independent filtering) | Confirm by counting `is.na(res$padj)`; to include them, set NA padj to 1 before plotting |
| Volcano has many "significant" genes but MA plot shows them all at low baseMean | Unshrunken LFC; the volcano is showing fold-change noise | Re-plot with shrunken LFC; the MA-plot fan is the diagnostic |
| Forest of horizontal stripes in MA at integer LFC | Pseudocount-induced quantization in low-count genes | Increase normalization-method aggressiveness OR filter low-count genes upstream |
| Half the genes have padj = NA | DESeq2 independent filtering (Bourgon-Gentleman-Huber 2010 *PNAS*) excluded them as low-mean | This is correct behavior; do NOT set `independentFiltering = FALSE` to hide it. Report the NA count |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Default LFC cutoff for "biologically relevant" | \|log2FC\| > 1 (2-fold) | Convention; sensitive analyses use 0.58 (1.5-fold) for subtle effects |
| Default FDR cutoff | padj < 0.05 | Benjamini-Hochberg 1995 *JRSS-B* 57:289 |
| Stricter cutoff for unbiased screens | padj < 0.01 | Reduces false positives in unbiased genome-wide analyses |
| Relaxed cutoff for exploratory / hypothesis-generating | padj < 0.10 or 0.20 | Acceptable for follow-up enrichment, NOT for "hits" |
| s-value cutoff (Stephens 2017) | s < 0.005 is only a rough match to padj < 0.05 (airway: 4,684 vs 3,994 genes, 3,842 shared) | Stephens 2017 *Biostatistics* 18:275 |
| Raster threshold | None required | 17,994 vector points gave 540 KB; rasterize the points for size or editing speed, keep axes vector |
| ggrepel max.overlaps | Set to Inf for publication | Default 10 silently drops labels |
| Volcano y-axis cap | Optional; airway reaches -log10(padj) 131 and is readable uncapped | If used, mark capped genes (triangles) and state the cap; clipping hid 39 significant genes at 50 |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why is this LFC shrunken? Show me the unshrunken." | Shrunken LFC is the recommended estimate for ranking and visualization (Zhu 2019). Unshrunken LFC inflates at low counts and gives misleading rank. Unshrunken is available in the supplementary table |
| "Why padj < 0.05 not p < 0.05?" | padj controls FDR via Benjamini-Hochberg. Raw p < 0.05 across 20000 genes yields ~1000 false positives by chance; padj < 0.05 caps the expected false-positive rate at 5% of called hits |
| "Why are X gene and Y gene not labeled?" | Labels selected by combined rank (-log10(p) * \|LFC\|) or pre-specified gene list. List of all significant genes is in supplementary table T1 |
| "The volcano looks too clean / too sparse." | Color encodes 3 categories (Up/Down/NS), not a gradient. Gradient encoding implies a continuous interpretation of significance which is invalid — significance is a threshold decision |
| "Why is the x-axis asymmetric?" | Asymmetric x-axis reflects the asymmetry of the data. If symmetry is preferred for visual interpretation, use `xlim = c(-X, X)` with X = max(\|LFC\|) |
