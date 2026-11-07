args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript cis_panel.R harmonized.tsv vep_annotations.tsv output.rds")
library(TwoSampleMR)
d <- read.delim(args[1], check.names = FALSE)
vep <- read.delim(args[2], check.names = FALSE)
if (!all(c("SNP", "Consequence") %in% names(vep))) stop("VEP table needs SNP and Consequence")
if ("mr_keep" %in% names(d)) d <- d[!is.na(d$mr_keep) & d$mr_keep, ]
if (!nrow(d)) stop("No retained instruments")
pav <- c("missense_variant", "stop_gained", "stop_lost", "frameshift_variant",
         "splice_acceptor_variant", "splice_donor_variant", "start_lost", "protein_altering_variant")
if (!all(d$SNP %in% vep$SNP) || any(is.na(vep$Consequence))) stop("Complete consequence annotation is required for every instrument")
is_pav <- vapply(strsplit(vep$Consequence, "[,&]"), function(x) any(trimws(x) %in% pav), logical(1))
d$is_pav <- d$SNP %in% vep$SNP[is_pav]
panel <- function(x) {
  if (!nrow(x)) return(NULL)
  methods <- if (nrow(x) == 1) "mr_wald_ratio" else if (nrow(x) >= 10)
    c("mr_ivw", "mr_egger_regression", "mr_weighted_median") else c("mr_ivw", "mr_weighted_median")
  mr(x, method_list = methods)
}
saveRDS(list(primary = panel(d), pav_excluded = panel(d[!d$is_pav, ]),
             n_instruments = nrow(d), n_after_pav_exclusion = sum(!d$is_pav)), args[3])
