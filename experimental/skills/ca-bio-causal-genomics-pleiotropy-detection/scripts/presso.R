args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) stop("Usage: Rscript presso.R harmonized.tsv output.rds [simulations]")
library(MRPRESSO)
d <- read.delim(args[1], check.names = FALSE)
if ("mr_keep" %in% names(d)) d <- d[!is.na(d$mr_keep) & d$mr_keep, ]
n <- if (length(args) == 3) as.integer(args[3]) else 5000L
if (nrow(d) < 4 || is.na(n) || n <= nrow(d)) stop("Need at least four instruments and more simulations than instruments")
set.seed(42)
result <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE, data = d, NbDistribution = n, SignifThreshold = 0.05)
saveRDS(result, args[2])
