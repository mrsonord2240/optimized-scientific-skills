# Volcano / MA failure modes

## Unshrunken LFC plotted as volcano

**Trigger:** Calling `results(dds)` and plotting `log2FoldChange` directly, without `lfcShrink()`.

**Mechanism:** ML estimate has infinite-variance tails at low counts; one read difference produces log2FC = Inf.

**Symptom:** "Top hits" by |LFC| are all genes with baseMean < 5; biologically interesting genes with moderate LFC are hidden in the noise cloud.

**Fix:** `lfcShrink(dds, coef=..., type='apeglm')` for the default case; `type='ashr'` if `contrast=` is needed.

## Raw p threshold line drawn on adjusted axis

**Trigger:** Drawing `geom_hline(yintercept = -log10(0.05))` on a plot whose y-axis is `-log10(padj)`.

**Mechanism:** padj < 0.05 corresponds to FDR < 5% control, NOT to raw p < 0.05. The drawn line is at the wrong y-value relative to the data.

**Symptom:** Visible "significant" points sit below the FDR line; the legend says "FDR < 0.05" but the line doesn't separate them correctly.

**Fix:** Be explicit: if y-axis is padj, use `-log10(fdr)` as the threshold line value AND label it "FDR threshold." If y-axis is raw p, an FDR threshold cannot be drawn as a horizontal line — the FDR threshold moves per gene.

## Top-N-by-p selects low-effect-size hits

**Trigger:** `head(arrange(res, pvalue), 20)` to choose labels.

**Mechanism:** With large N, the smallest p-values belong to high-count, low-variance, biologically-boring genes (housekeeping). Effect size and statistical confidence are not the same thing.

**Symptom:** Labels are GAPDH, ACTB, B2M — never the gene that drives the biology.

**Fix:** Rank by `-log10(p) * abs(log2FoldChange)` (geometric average of the two axes), OR pre-specify labels of interest from prior knowledge.

## ggrepel `max.overlaps` silently drops labels

**Trigger:** Default `max.overlaps = 10`; 30 genes labeled in code; only 10 render.

**Mechanism:** ggrepel emits a warning ("18 unlabeled data points (too many overlaps)") but no error. In a Quarto/Rmd render the warning is buried in the log.

**Symptom:** Reviewer asks "where is gene X?"; the label was specified in code but did not render.

**Fix:** `geom_text_repel(..., max.overlaps = Inf)` or `options(ggrepel.max.overlaps = Inf)` at the top of the script.

## Extreme p-values compress the upper axis

**Trigger:** Genes with p = 1e-200 or smaller (common in cancer datasets) push the y-axis maximum to 200; all biologically meaningful genes pile up at the bottom.

**Mechanism:** -log10 expands the tail; one ultra-significant gene visually dominates.

**Symptom:** Volcano looks like an Eiffel Tower with most genes squished near y = 0-20.

**Fix:** Cap only when a few extreme values really compress the rest, and mark what is capped: `volcano_plot(res, y_cap = 30)` squishes higher genes onto the cap as triangles. `coord_cartesian(ylim = ...)` alone clips points out of view; on the airway data a cap of 50 hid 39 significant genes. Alternatives: `sqrt(-log10(p))` to compress the tail, or `ggbreak::scale_y_break()` (ggbreak 0.1.7 runs on a volcano; the break hides the values inside it, so state the range).

## `lfcShrink(type='normal')` on a modern DESeq2

**Trigger:** Following old tutorials that pre-date DESeq2 v1.28 when apeglm became the default.

**Mechanism:** The `'normal'` prior over-shrinks large real effects toward zero.

**Symptom:** Volcano looks "too clean" — genuine 8-fold changes appear as 2-3 fold.

**Fix:** Use `type='apeglm'` (the `lfcShrink` default in DESeq2 1.46.0) or `type='ashr'`.

## EnhancedVolcano `selectLab` genes that do not appear

**Trigger:** A gene listed in `selectLab` is not labelled.

**Mechanism:** On EnhancedVolcano 1.24.0 a listed gene that fails `pCutoff`/`FCcutoff` is still labelled (checked). A name that is absent from `lab` (misspelled, or symbol vs Ensembl ID) is silently ignored, and so is a gene whose `padj` is `NA`.

**Fix:** Compare `selectLab` with `lab` (`setdiff(selectLab, rownames(res))`) and check `is.na(res[gene, 'padj'])`.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| All "top hits" have baseMean < 5 | Unshrunken LFC | `lfcShrink(type='apeglm')` |
| Threshold line doesn't separate colored from grey points | y = pvalue but threshold drawn at FDR | Switch y to padj OR redraw line at FDR-equivalent p |
| Labeled gene does not appear in EnhancedVolcano | Name not in `lab`, or `padj = NA` | `setdiff(selectLab, rownames(res))`; check `is.na(padj)` |
| Volcano renders as a flat horizontal cloud | Extreme p (e.g., 1e-200) dominates y-axis | `volcano_plot(res, y_cap = ...)` (capped genes drawn as triangles) or sqrt transform |
| PDF is large or slow to edit | Vector scatter of many thousand points (17,994 points: 540 KB vector, 38 KB rasterized) | `ggrastr::rasterise(geom_point())` or matplotlib `rasterized=True` |
| ggrepel labels 10 of 30 selected genes | Default `max.overlaps = 10` | `geom_text_repel(max.overlaps = Inf)` |
| Volcano "significant" gene count differs from DESeq2 summary | EnhancedVolcano drops `padj = NA` | Set NA padj to 1 or document the discrepancy |
| Up and Down counts asymmetric for a balanced experiment | Library-size normalization failure | Re-run DESeq2 with `estimateSizeFactors(type='poscounts')` for sparse data |
| EnhancedVolcano Up and Down share one colour | `col =` has one colour for both directions | `colCustom` named per gene (see SKILL.md) |
| `lfcShrink` result has no `svalue` | `svalue = TRUE` not passed (apeglm and ashr both support it; `normal` errors) | `lfcShrink(..., type = 'apeglm', svalue = TRUE)`; the result then has no `padj` |
