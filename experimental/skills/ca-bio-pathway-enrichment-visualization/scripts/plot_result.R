args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop('usage: plot_result.R result.rds dotplot|nes|running|emap output.pdf count_or_index')
library(enrichplot)
library(ggplot2)
result <- readRDS(args[1])
n <- nrow(as.data.frame(result))
if (n == 0) stop('result contains zero rows; report the empty result')
number <- as.integer(args[4])
if (is.na(number) || number < 1) stop('count or index must be positive')
if (args[2] == 'dotplot') {
    plot <- dotplot(result, showCategory = number)
} else if (args[2] == 'nes') {
    if (!inherits(result, 'gseaResult')) stop('NES requires gseaResult')
    plot <- dotplot(result, x = 'NES', showCategory = number, color = 'p.adjust')
} else if (args[2] == 'running') {
    if (!inherits(result, 'gseaResult') || number > n) stop('invalid GSEA set index')
    plot <- gseaplot2(result, geneSetID = number)
} else if (args[2] == 'emap') {
    plot <- emapplot(pairwise_termsim(result), showCategory = number)
} else stop('unknown mode')
pdf(args[3], width = 10, height = 8)
print(plot)
dev.off()
message('Result rows: ', n, '; requested count or index: ', number)
