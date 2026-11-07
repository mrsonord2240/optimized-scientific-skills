args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 4) stop("Usage: Rscript qvalue_summary.R input.csv p_column output.rds [alpha]")
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "qvalue_safe.R"))
d <- read.csv(args[1], check.names = FALSE)
if (!args[2] %in% names(d)) stop("p-value column missing")
p <- d[[args[2]]]
if (!is.numeric(p) || any(!is.finite(p) | p < 0 | p > 1)) stop("p-values must be finite and in [0,1]")
alpha <- if (length(args) == 4) as.numeric(args[4]) else 0.05
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
saveRDS(qvalue_safe(p, alpha), args[3])
