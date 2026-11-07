# Matrix factor, extraction recovery and carryover for a targeted LC-MS/MS assay (base R only).
#
# Modes (key=value arguments; column names are configurable):
#   matrix_factor  file=<csv> [analyte_col=analyte_area] [istd_col=istd_area]
#                  [neat_analyte_col=neat_analyte_area] [neat_istd_col=neat_istd_area]
#                  [lot_col=lot] [filter=<col>:<value>] [tol_cv=15]
#       One row per matrix lot: post-extraction-spiked analyte/IS areas plus the neat-solvent
#       analyte/IS areas at the same concentration. MF = matrix response / neat response;
#       IS-normalized MF = MF_analyte / MF_IS. ICH M10 / Matuszewski 2003: CV of the
#       IS-normalized MF across >= 6 lots <= tol_cv (15%).
#   recovery       file=<csv> [pre_col=pre_area] [post_col=post_area] [level_col=level]
#       One row per replicate: analyte area from a pre-extraction spike and from a post-extraction
#       spike at the same level. Recovery % = pre / post * 100, reported per level with its CV.
#       ICH M10 sets no fixed limit: recovery must be consistent and reproducible across levels.
#   carryover      blank_ratio=<analyte/IS ratio of the blank injected after the ULOQ>
#                  slope=<calibration slope> intercept=<calibration intercept> lloq=<conc>
#                  [blank_istd=<IS area in a blank injected WITHOUT IS> ref_istd=<mean IS area>]
#       Carryover % of LLOQ = blank response / LLOQ response * 100 (limit 20%); the LLOQ response
#       is intercept + slope * lloq from the same response-ratio calibration. IS carryover (limit
#       5%) needs a blank injected without IS: blank_istd / ref_istd * 100.
#
# Usage: Rscript matrix_recovery_carryover.R matrix_factor file=lots.csv analyte_col=quantifier_area filter=sample_type:QC_mid
#        Rscript matrix_recovery_carryover.R recovery file=recovery.csv
#        Rscript matrix_recovery_carryover.R carryover blank_ratio=0.00852 slope=0.0010010 intercept=-0.0000366 lloq=2
# Exit status is 0 when the check passes (or, for recovery, always) and 1 when it fails.

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
kv <- strsplit(args[-1], "=", fixed = TRUE)
opt <- setNames(lapply(kv, function(x) paste(x[-1], collapse = "=")), sapply(kv, `[`, 1))
get <- function(name, default = NULL) if (!is.null(opt[[name]])) opt[[name]] else default
need <- function(name) if (is.null(opt[[name]])) stop("missing argument: ", name) else opt[[name]]
cv_pct <- function(x) sd(x) / mean(x) * 100

if (mode == "matrix_factor") {
  d <- read.csv(need("file"))
  if (!is.null(get("filter"))) {
    f <- strsplit(get("filter"), ":", fixed = TRUE)[[1]]
    d <- d[d[[f[1]]] == f[2], ]
  }
  a <- d[[get("analyte_col", "analyte_area")]]
  i <- d[[get("istd_col", "istd_area")]]
  na <- d[[get("neat_analyte_col", "neat_analyte_area")]]
  ni <- d[[get("neat_istd_col", "neat_istd_area")]]
  d$mf_analyte <- a / na
  d$mf_istd <- i / ni
  d$mf_is_norm <- d$mf_analyte / d$mf_istd
  lot <- get("lot_col", "lot")
  out <- d[, c(if (lot %in% names(d)) lot, "mf_analyte", "mf_istd", "mf_is_norm")]
  print(format(out, digits = 4), row.names = FALSE)
  cv <- cv_pct(d$mf_is_norm)
  tol <- as.numeric(get("tol_cv", 15))
  cat(sprintf("lots: %d | mean IS-normalized MF: %.3f | CV: %.2f%% (limit %g%%, >= 6 lots)\n",
              nrow(d), mean(d$mf_is_norm), cv, tol))
  if (nrow(d) < 6) cat("WARNING: fewer than 6 lots, ICH M10 requires >= 6\n")
  pass <- cv <= tol && nrow(d) >= 6
  cat("MATRIX FACTOR PASS:", pass, "\n")
  quit(status = if (pass) 0 else 1)
}

if (mode == "recovery") {
  d <- read.csv(need("file"))
  d$recovery_pct <- d[[get("pre_col", "pre_area")]] / d[[get("post_col", "post_area")]] * 100
  lv <- d[[get("level_col", "level")]]
  by_level <- split(d$recovery_pct, factor(lv, levels = unique(lv)))
  res <- data.frame(level = names(by_level), n = lengths(by_level),
                    mean_recovery_pct = sapply(by_level, mean), cv_pct = sapply(by_level, cv_pct))
  print(format(res, digits = 4), row.names = FALSE)
  cat(sprintf("range of level means: %.1f-%.1f%% (spread %.1f points; judge consistency, not closeness to 100%%)\n",
              min(res$mean_recovery_pct), max(res$mean_recovery_pct),
              max(res$mean_recovery_pct) - min(res$mean_recovery_pct)))
  quit(status = 0)
}

if (mode == "carryover") {
  blank <- as.numeric(need("blank_ratio"))
  lloq_resp <- as.numeric(need("intercept")) + as.numeric(need("slope")) * as.numeric(need("lloq"))
  co <- blank / lloq_resp * 100
  cat(sprintf("blank response %.5f | LLOQ response %.5f | analyte carryover %.1f%% of LLOQ (limit 20%%)\n",
              blank, lloq_resp, co))
  pass <- co <= 20
  if (!is.null(get("blank_istd"))) {
    co_is <- as.numeric(get("blank_istd")) / as.numeric(need("ref_istd")) * 100
    cat(sprintf("IS carryover %.1f%% of the reference IS response (limit 5%%)\n", co_is))
    pass <- pass && co_is <= 5
  } else {
    cat("IS carryover not assessed (needs a blank injected without IS: blank_istd, ref_istd)\n")
  }
  cat("CARRYOVER PASS:", pass, "\n")
  quit(status = if (pass) 0 else 1)
}

stop("first argument must be matrix_factor, recovery or carryover")
