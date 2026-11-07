args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript proteomics_power.R scenarios.csv output.csv")
library(pwr)
d <- read.csv(args[1], check.names = FALSE)
if (!all(c("d", "power", "alpha", "n_proteins") %in% names(d))) stop("Missing scenario columns")
for (i in seq_len(nrow(d))) {
  x <- d[i, ]
  if (any(!is.finite(unlist(x[c("d", "power", "alpha", "n_proteins")]))) ||
      x$d <= 0 || x$power <= 0 || x$power >= 1 || x$alpha <= 0 ||
      x$alpha >= 1 || x$n_proteins < 1 || x$n_proteins != floor(x$n_proteins)) stop("Invalid scenario")
  d$raw_n[i] <- ceiling(pwr.t.test(d = x$d, sig.level = x$alpha, power = x$power, type = "two.sample")$n)
  d$bonferroni_n[i] <- ceiling(pwr.t.test(d = x$d, sig.level = x$alpha / x$n_proteins, power = x$power, type = "two.sample")$n)
}
write.csv(d, args[2], row.names = FALSE)
