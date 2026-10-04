# Reference: ggplot2 3.5+ (run clean on 4.0.3), ggrepel 0.9.5+ | Verify API if version differs
# volcano_plot(res, fdr, lfc_threshold, label_genes, top_n, y_cap): res is a DESeq2 results object (shrunken LFC) or data frame with rownames, padj, log2FoldChange.
# y is -log10(padj), so the dashed FDR line is exactly the Up/Down/NS colour boundary. Genes with padj = NA (independent filtering) are not drawn.
# y_cap = NULL draws every gene. A numeric y_cap squishes higher points onto the cap and draws them as triangles.
library(ggplot2)
library(ggrepel)
library(dplyr)

volcano_plot <- function(res, fdr = 0.05, lfc_threshold = 1, label_genes = NULL, top_n = 10, y_cap = NULL) {
    res <- as.data.frame(res) %>%
        tibble::rownames_to_column('gene')
    n_na <- sum(is.na(res$padj))
    if (n_na > 0) message(n_na, ' genes with padj = NA are not drawn')
    res <- res %>%
        filter(!is.na(padj)) %>%
        mutate(
            significance = case_when(
                padj < fdr & log2FoldChange > lfc_threshold ~ 'Up',
                padj < fdr & log2FoldChange < -lfc_threshold ~ 'Down',
                TRUE ~ 'NS'
            ),
            neg_log10_padj = -log10(padj)
        )

    if (is.null(label_genes)) {
        label_genes <- res %>%
            filter(significance != 'NS') %>%
            mutate(rank_score = neg_log10_padj * abs(log2FoldChange)) %>%
            arrange(desc(rank_score)) %>%
            head(top_n) %>%
            pull(gene)
    } else if (!all(label_genes %in% res$gene)) {
        warning('label_genes not found in rownames (or padj is NA): ',
                paste(setdiff(label_genes, res$gene), collapse = ', '))
    }
    res$label <- ifelse(res$gene %in% label_genes, res$gene, '')

    res$capped <- !is.null(y_cap) & res$neg_log10_padj > (if (is.null(y_cap)) Inf else y_cap)
    res$y <- if (is.null(y_cap)) res$neg_log10_padj else pmin(res$neg_log10_padj, y_cap)

    okabe_ito <- c(Up = '#D55E00', Down = '#0072B2', NS = '#999999')

    p <- ggplot(res, aes(log2FoldChange, y, color = significance, shape = capped)) +
        geom_point(alpha = 0.6, size = 1.3) +
        scale_color_manual(values = okabe_ito, name = NULL) +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                           labels = c(`FALSE` = 'plotted value', `TRUE` = paste('capped at', y_cap)),
                           name = NULL, drop = FALSE) +
        geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                   linetype = 'dashed', color = 'grey40', linewidth = 0.3) +
        geom_hline(yintercept = -log10(fdr), linetype = 'dashed',
                   color = 'grey40', linewidth = 0.3) +
        geom_text_repel(aes(label = label), color = 'black', size = 3,
                        max.overlaps = Inf, box.padding = 0.4, segment.size = 0.2,
                        min.segment.length = 0.3) +
        labs(x = expression(log[2]~'fold change (shrunken)'),
             y = expression(-log[10]~'adjusted'~italic(p))) +
        theme_classic(base_size = 10) +
        theme(panel.grid = element_blank())
    if (is.null(y_cap)) p <- p + guides(shape = 'none')
    p
}
