# Reference: ggplot2 3.5+ (run clean on 4.0.3) | Verify API if version differs
# Every helper returns a ggplot already wearing theme_publication(); multi-panel composes them as-is.
library(ggplot2)
library(patchwork)
library(ggrepel)
library(dplyr)

# Okabe-Ito (CVD-safe); the one categorical palette for every helper
okabe_ito <- c('#0072B2', '#D55E00', '#009E73', '#E69F00', '#56B4E9', '#CC79A7', '#F0E442', '#000000')

theme_publication <- function(base_size = 10) {
    theme_classic(base_size = base_size) +
    theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black', linewidth = 0.3),
        axis.line = element_line(color = 'black', linewidth = 0.3),
        legend.position = 'right',
        legend.key.size = unit(0.4, 'cm'),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold', size = 9),
        plot.title = element_text(face = 'bold', size = 11),
        plot.tag = element_text(face = 'bold', size = 11)
    )
}

scale_okabe <- function(n, aesthetics = 'color') {
    if (n > length(okabe_ito)) stop('Okabe-Ito has ', length(okabe_ito), ' colours; got ', n, ' groups. Collapse groups or facet.')
    scale_color_manual(values = okabe_ito[seq_len(n)], aesthetics = aesthetics)
}

# res needs columns: log2FoldChange, padj, and a label column (default 'gene').
# The y axis is -log10(padj), so the horizontal threshold is exactly the FDR cutoff.
# top_n = NULL (default) labels the 10 smallest padj when those labels are short (<= 8 characters, gene symbols)
# and only the 3 smallest when they are long (Ensembl IDs); pass top_n to override.
create_volcano <- function(res, fdr_threshold = 0.05, lfc_threshold = 1, top_n = NULL,
                           label_col = 'gene') {
    res <- res %>% filter(!is.na(padj))
    if (is.null(top_n)) {
        lead <- as.character(res[[label_col]][order(res$padj)][1:min(10, nrow(res))])
        top_n <- if (max(nchar(lead)) <= 8) 10 else 3
    }
    res <- res %>%
        mutate(
            significance = case_when(
                padj < fdr_threshold & log2FoldChange > lfc_threshold ~ 'Up',
                padj < fdr_threshold & log2FoldChange < -lfc_threshold ~ 'Down',
                TRUE ~ 'NS'
            ),
            label = ifelse(rank(padj, ties.method = 'first') <= top_n & padj < fdr_threshold,
                           as.character(.data[[label_col]]), '')
        )

    ggplot(res, aes(log2FoldChange, -log10(padj), color = significance)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_text_repel(aes(label = label), color = 'black', size = 2.5, max.overlaps = Inf,
                        min.segment.length = 0, box.padding = 0.5, force = 4, max.iter = 20000, max.time = 5, seed = 1) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
        scale_color_manual(values = c(Up = '#D55E00', Down = '#0072B2', NS = 'grey60')) +
        geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                   linetype = 'dashed', color = 'grey40') +
        geom_hline(yintercept = -log10(fdr_threshold),
                   linetype = 'dashed', color = 'grey40') +
        labs(x = expression(log[2]~fold~change),
             y = expression(-log[10]~(adjusted~italic(P))),
             color = 'Significance') +
        theme_publication()
}

create_boxplot <- function(df, x_var, y_var, fill_var = NULL) {
    p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]]))

    if (!is.null(fill_var)) {
        p <- p + aes(fill = .data[[fill_var]]) +
            scale_okabe(nlevels(factor(df[[fill_var]])), 'fill')
    }

    p +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(position = position_jitter(width = 0.2, seed = 1), alpha = 0.5, size = 1.5) +
        labs(x = NULL) +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# pca_df needs PC1, PC2 and var_explained in PERCENT (0-100, one value per PC in PC order; row i = PCi).
# Fractions (0-1) are not rescaled: they would print as "PC1 (0.6%)", so a first-two sum <= 1 raises a warning.
create_pca_plot <- function(pca_df, color_var, shape_var = NULL) {
    if (!'var_explained' %in% names(pca_df))
        stop("pca_df needs a 'var_explained' column (percent variance, first two rows = PC1, PC2)")
    if (sum(pca_df$var_explained[1:2]) <= 1)
        warning("var_explained for PC1+PC2 sums to <= 1: it is read as percent, so axis labels show values under 1%. Multiply fractions by 100.")
    p <- ggplot(pca_df, aes(PC1, PC2, color = .data[[color_var]]))

    if (!is.null(shape_var)) {
        p <- p + aes(shape = .data[[shape_var]])
    }

    p +
        geom_point(size = 3, alpha = 0.8) +
        scale_okabe(nlevels(factor(pca_df[[color_var]]))) +
        labs(x = paste0('PC1 (', round(pca_df$var_explained[1], 1), '%)'),
             y = paste0('PC2 (', round(pca_df$var_explained[2], 1), '%)')) +
        theme_publication()
}

# Width in mm: 89 = single column, 183 = double column.
save_publication_figure <- function(plot, filename, width = 89, height = 70) {
    ggsave(paste0(filename, '.pdf'), plot, width = width, height = height, units = 'mm', device = cairo_pdf)
    ggsave(paste0(filename, '.png'), plot, width = width, height = height, units = 'mm', dpi = 300)
}

create_multi_panel <- function(p1, p2, p3, p4 = NULL) {
    if (is.null(p4)) {
        combined <- (p1 | p2) / p3 +
            plot_annotation(tag_levels = 'A') +
            plot_layout(heights = c(1, 1))
    } else {
        combined <- (p1 | p2) / (p3 | p4) +
            plot_annotation(tag_levels = 'A')
    }
    combined
}
