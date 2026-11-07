args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: Rscript adjust_pvalues.R input.csv p_column method output.csv")
d <- read.csv(args[1], check.names = FALSE)
if (!args[2] %in% names(d)) stop("p-value column missing")
p <- d[[args[2]]]
if (!is.numeric(p) || any(!is.na(p) & (!is.finite(p) | p < 0 | p > 1))) stop("p-values must be numeric, missing or in [0,1]")
if (!args[3] %in% c("BH", "BY", "holm", "bonferroni")) stop("method must be BH, BY, holm or bonferroni")
d$adjusted_p <- p.adjust(p, method = args[3])
write.csv(d, args[4], row.names = FALSE)
