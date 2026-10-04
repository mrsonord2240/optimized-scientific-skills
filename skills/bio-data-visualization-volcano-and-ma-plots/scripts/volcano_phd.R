# Reference: DESeq2 1.42+ (run clean on 1.46.0), EnhancedVolcano 1.20+, ggplot2 3.5+ (4.0.3), ggrepel 0.9.5+ | Verify API if version differs

# PhD-level volcano + MA plot from DESeq2 results.
# Encodes the four correctness traps: (1) shrunken LFC, (2) padj on y-axis with the
# threshold line on the same quantity, (3) label by combined rank, (4) max.overlaps = Inf.

library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)

# 1. SHRINKAGE -- apeglm by default (Zhu et al 2019)
# coef = name from resultsNames(dds); ashr is the alternative when contrast= is needed
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
# alt: res <- lfcShrink(dds, contrast = c('condition', 'treated', 'control'), type = 'ashr')

# 2. CATEGORICAL SIGNIFICANCE -- threshold on padj NOT pvalue
fdr <- 0.05
lfc_threshold <- 1
res_df <- as.data.frame(res) %>%
    rownames_to_column('gene') %>%
    mutate(
        significance = case_when(
            is.na(padj) ~ 'NS',
            padj < fdr & log2FoldChange >  lfc_threshold ~ 'Up',
            padj < fdr & log2FoldChange < -lfc_threshold ~ 'Down',
            TRUE ~ 'NS'),
        neg_log10_padj = -log10(padj))

# 3. LABEL SELECTION -- combined rank, not pure top-N-by-p
# genes_of_interest must match rownames(res); with Ensembl IDs as rownames, map symbols first
genes_of_interest <- c('TP53', 'MYC', 'BRCA1')
top_by_rank <- res_df %>%
    filter(significance != 'NS') %>%
    mutate(rank_score = neg_log10_padj * abs(log2FoldChange)) %>%
    arrange(desc(rank_score)) %>% head(10) %>% pull(gene)
labels <- union(genes_of_interest, top_by_rank)
missing_genes <- setdiff(genes_of_interest, res_df$gene)
if (length(missing_genes) > 0) warning('genes_of_interest not in rownames(res): ', paste(missing_genes, collapse = ', '))
res_df$label <- ifelse(res_df$gene %in% labels, res_df$gene, '')

# 4. PLOT -- Okabe-Ito categorical palette (Wong 2011)
okabe_ito <- c(Up = '#D55E00', Down = '#0072B2', NS = '#999999')

# y is -log10(padj); padj = NA genes (independent filtering) are not drawn.
# Optional cap: set y_cap only when a few extreme values compress the rest (e.g. 50);
# capped genes are squished onto the cap and drawn as triangles so none are hidden.
y_cap <- NULL
res_df <- filter(res_df, !is.na(padj))
res_df$capped <- !is.null(y_cap) & res_df$neg_log10_padj > (if (is.null(y_cap)) Inf else y_cap)
res_df$y <- if (is.null(y_cap)) res_df$neg_log10_padj else pmin(res_df$neg_log10_padj, y_cap)

p_volcano <- ggplot(res_df, aes(log2FoldChange, y, color = significance, shape = capped)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_color_manual(values = okabe_ito, name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = 'none') +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = 'dashed', color = 'grey40', linewidth = 0.3) +
    # y is -log10(padj), so this line is exactly the padj < fdr colour boundary
    geom_hline(yintercept = -log10(fdr), linetype = 'dashed', color = 'grey40', linewidth = 0.3) +
    geom_text_repel(aes(label = label), color = 'black', size = 3,
                    max.overlaps = Inf,                    # critical -- default 10 drops labels
                    box.padding = 0.4, segment.size = 0.2, min.segment.length = 0.3) +
    labs(x = expression(log[2]~'fold change (shrunken)'),
         y = expression(-log[10]~'adjusted'~italic(p))) +
    theme_classic(base_size = 10) +
    theme(panel.grid = element_blank())

# 5. MA PLOT as the shrinkage diagnostic
# fan shape with extreme |LFC| at low baseMean indicates inadequate shrinkage
p_ma <- ggplot(res_df, aes(log10(baseMean), log2FoldChange,
                            color = significance != 'NS')) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c(`FALSE` = '#999999', `TRUE` = '#D55E00'),
                       guide = 'none') +
    geom_hline(yintercept = 0, color = 'black', linewidth = 0.4) +
    labs(x = expression(log[10]~'mean normalized count'),
         y = expression(log[2]~'fold change (shrunken)')) +
    theme_classic(base_size = 10) +
    theme(panel.grid = element_blank())

# 6. SAVE -- cairo_pdf embeds TrueType (17,994 drawn vector points gave 0.56 and 0.69 MB PDFs on the airway data;
# wrap the point layer in ggrastr::rasterise() if a PDF editor struggles)
ggsave('volcano.pdf', p_volcano, width = 89, height = 90, units = 'mm',
       device = cairo_pdf)
ggsave('ma_plot.pdf', p_ma, width = 89, height = 70, units = 'mm',
       device = cairo_pdf)

# 7. ENHANCEDVOLCANO equivalent
# colCustom (named by gene) is the only way to colour Up and Down differently:
# col[4] is shared by both directions. selectLab genes that fail pCutoff/FCcutoff
# are still labelled on EnhancedVolcano 1.24.0 (checked), so no manual workaround is needed.
ev_col <- setNames(unname(okabe_ito[as.character(res_df$significance)]), as.character(res_df$significance))
ev_res <- res[res_df$gene, ]
EnhancedVolcano(ev_res,
    lab = rownames(ev_res),
    x = 'log2FoldChange', y = 'padj',                       # y on padj NOT pvalue
    pCutoff = fdr, FCcutoff = lfc_threshold,
    selectLab = labels,
    colCustom = ev_col,
    ylab = bquote(-log[10]~'adjusted'~italic(P)),
    title = NULL, subtitle = NULL, caption = NULL,
    drawConnectors = TRUE, widthConnectors = 0.3,
    maxoverlapsConnectors = Inf,                            # match ggrepel guidance
    pointSize = 1.5, labSize = 3, colAlpha = 0.6,
    legendPosition = 'right')
