args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop('usage: plot_csv.R box|pca|volcano data.csv output-prefix [columns]')
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'publication_figures.R'))
df <- read.csv(args[2], check.names = FALSE)
if (args[1] == 'box' && length(args) == 5) {
    if (!all(args[4:5] %in% names(df))) stop('missing box columns')
    p <- create_boxplot(df, args[4], args[5])
} else if (args[1] == 'pca' && length(args) == 4) {
    if (!all(c('PC1', 'PC2', 'var_explained', args[4]) %in% names(df))) stop('missing PCA columns')
    p <- create_pca_plot(df, args[4])
} else if (args[1] == 'volcano' && length(args) == 3) {
    if (!all(c('log2FoldChange', 'padj', 'gene') %in% names(df))) stop('missing volcano columns')
    p <- create_volcano(df)
} else stop('invalid mode or arguments')
save_publication_figure(p, args[3])
