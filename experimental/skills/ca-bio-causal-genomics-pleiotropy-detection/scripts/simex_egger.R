args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript simex_egger.R harmonized.tsv output.rds")
library(simex)
library(TwoSampleMR)
d <- read.delim(args[1], check.names = FALSE)
if ("mr_keep" %in% names(d)) d <- d[!is.na(d$mr_keep) & d$mr_keep, ]
isq <- Isq(d$beta.exposure, d$se.exposure)
if (isq < 0.6) stop("I-squared below 0.6: SIMEX unreliable; consider another method")
if (isq >= 0.9) stop("I-squared at least 0.9: source does not apply SIMEX")
w <- 1 / d$se.outcome^2
fit <- lm(beta.outcome ~ beta.exposure, weights = w, data = d, x = TRUE, y = TRUE)
set.seed(42)
corrected <- simex(fit, SIMEXvariable = "beta.exposure", measurement.error = d$se.exposure,
  lambda = seq(0.5, 2, 0.5), B = 1000, fitting.method = "quadratic", asymptotic = FALSE)
saveRDS(list(isq = isq, egger = fit, simex = corrected), args[2])
