args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript sensitivity_battery.R harmonized.tsv output.rds")
library(TwoSampleMR)
d <- read.delim(args[1], check.names = FALSE)
if ("mr_keep" %in% names(d)) d <- d[!is.na(d$mr_keep) & d$mr_keep, ]
set.seed(42)
result <- list(
  estimates = mr(d, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")),
  heterogeneity = mr_heterogeneity(d),
  egger_intercept = mr_pleiotropy_test(d),
  directionality = directionality_test(d),
  steiger = steiger_filtering(d),
  leave_one_out = mr_leaveoneout(d)
)
saveRDS(result, args[2])
