args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(2, 3))) stop('usage: volcano_csv.R results.csv output.pdf [y_cap]')
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'volcano_plot.R'))
df <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c('gene', 'log2FoldChange', 'padj') %in% names(df))) stop('missing required columns')
if (anyNA(df$gene) || any(!nzchar(df$gene)) || anyDuplicated(df$gene)) stop('gene identifiers must be unique and nonempty')
if (any(!is.finite(df$log2FoldChange))) stop('effects must be finite')
valid_p <- df$padj[!is.na(df$padj)]
if (any(!is.finite(valid_p)) || any(valid_p <= 0 | valid_p > 1)) stop('nonmissing padj must be in (0, 1]')
rownames(df) <- df$gene
df$gene <- NULL
cap <- if (length(args) == 3) as.numeric(args[3]) else NULL
if (!is.null(cap) && (!is.finite(cap) || cap <= 0)) stop('y_cap must be positive and finite')
p <- volcano_plot(df, y_cap = cap)
ggplot2::ggsave(args[2], p, width = 89, height = 90, units = 'mm', device = cairo_pdf)
