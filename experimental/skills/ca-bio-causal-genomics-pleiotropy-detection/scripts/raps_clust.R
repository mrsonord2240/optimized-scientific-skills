args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript raps_clust.R harmonized.tsv output.rds")
library(TwoSampleMR)
library(mrclust)
d <- read.delim(args[1], check.names = FALSE)
if ("mr_keep" %in% names(d)) d <- d[!is.na(d$mr_keep) & d$mr_keep, ]
raps <- mr_raps(b_exp = d$beta.exposure, b_out = d$beta.outcome,
  se_exp = d$se.exposure, se_out = d$se.outcome,
  parameters = list(over.dispersion = TRUE, loss.function = "huber", shrinkage = FALSE))
clusters <- mr_clust_em(theta = d$beta.outcome / d$beta.exposure,
  theta_se = abs(d$se.outcome / d$beta.exposure), bx = d$beta.exposure,
  by = d$beta.outcome, bxse = d$se.exposure, byse = d$se.outcome, obs_names = d$SNP)
saveRDS(list(raps = raps, clusters = clusters), args[2])
