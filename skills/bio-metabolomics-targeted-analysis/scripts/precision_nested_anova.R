# Intra-day and inter-day precision of QC replicates, inter-day by one-way ANOVA variance
# components (nested design), not by pooling every replicate into one SD.
#
# Input CSV: one row per QC replicate with a day column and a measured-concentration column.
# Arguments (key=value): file=<csv> [day_col=day] [value_col=measured_conc]
#                        [level_col=qc_level] [level=<one level>] [tol=15]
#   With level_col present, every level is reported (tolerance 20 at a level named LLOQ,
#   otherwise `tol`); level=<name> restricts to one level. Without a level column, all rows
#   are treated as one level.
# Method: intra-day CV = SD/mean within each day. Inter-day CV% = sqrt(within + between) / grand
#   mean * 100, within = MSwithin of aov(conc ~ factor(day)), between = max(0, (MSbetween -
#   MSwithin) / replicates_per_day). Needs >= 2 days with equal replicate counts.
# Usage: Rscript precision_nested_anova.R file=qc.csv [level=LLOQ]
# Exit status is 0 when every intra-day and inter-day CV is within tolerance, 1 otherwise.

args <- commandArgs(trailingOnly = TRUE)
kv <- strsplit(args, "=", fixed = TRUE)
opt <- setNames(lapply(kv, function(x) paste(x[-1], collapse = "=")), sapply(kv, `[`, 1))
get <- function(name, default = NULL) if (!is.null(opt[[name]])) opt[[name]] else default
if (is.null(opt$file)) stop("missing argument: file")

d <- read.csv(opt$file)
day_col <- get("day_col", "day"); val_col <- get("value_col", "measured_conc")
level_col <- get("level_col", "qc_level")
if (!level_col %in% names(d)) { d[[level_col]] <- "all"; }
if (!is.null(get("level"))) d <- d[d[[level_col]] == get("level"), ]

all_pass <- TRUE
for (lv in unique(d[[level_col]])) {
  qc <- d[d[[level_col]] == lv, ]
  qc$day <- qc[[day_col]]; qc$measured_conc <- qc[[val_col]]
  tol <- if (lv == "LLOQ") 20 else as.numeric(get("tol", 15))
  reps <- table(qc$day)
  if (length(reps) < 2 || length(unique(reps)) != 1) stop("level ", lv, ": need >= 2 days with equal replicate counts")

  intra <- aggregate(measured_conc ~ day, qc, function(x) sd(x) / mean(x) * 100)
  naive_pooled_cv <- sd(qc$measured_conc) / mean(qc$measured_conc) * 100  # WRONG for inter-day, shown for contrast

  fit <- aov(measured_conc ~ factor(day), data = qc)
  ms <- summary(fit)[[1]][["Mean Sq"]]
  n_per_day <- nrow(qc) / length(unique(qc$day))
  var_within <- ms[2]
  var_between <- max(0, (ms[1] - ms[2]) / n_per_day)
  inter_day_cv <- sqrt(var_within + var_between) / mean(qc$measured_conc) * 100

  cat(sprintf("%s (tolerance %g%%): intra-day CV by day: %s | naive pooled CV %.1f%% (not valid) | nested-ANOVA inter-day CV %.1f%% -> %s\n",
              lv, tol, paste(sprintf("d%s %.1f%%", intra$day, intra$measured_conc), collapse = ", "),
              naive_pooled_cv, inter_day_cv,
              if (all(intra$measured_conc <= tol) && inter_day_cv <= tol) "PASS" else "FAIL"))
  all_pass <- all_pass && all(intra$measured_conc <= tol) && inter_day_cv <= tol
}
cat("ALL LEVELS PASS:", all_pass, "\n")
quit(status = if (all_pass) 0 else 1)
