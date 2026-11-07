# Purpose: validate subject pairing, compute a paired Wilcoxon test, and annotate a paired plot.
# Inputs: CSV with subject_id,time,value; output image; p.adjust method.
# Usage: Rscript scripts/annotate_paired.R paired.csv paired.png holm

annotate_paired <- function(input_csv, output_plot, adjust_method = "holm") {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggpubr)
    library(rstatix)
  })

  df <- read.csv(input_csv, check.names = FALSE)
  stopifnot(all(c("subject_id", "time", "value") %in% names(df)))
  if (anyNA(df[, c("subject_id", "time", "value")])) stop("paired fields must not be missing")
  time_levels <- unique(df$time)
  if (length(time_levels) != 2L) stop("paired workflow requires exactly two time levels")
  if (any(duplicated(df[, c("subject_id", "time")]))) stop("duplicate subject/time rows")

  id_sets <- lapply(time_levels, function(level) sort(as.character(df$subject_id[df$time == level])))
  if (!identical(id_sets[[1]], id_sets[[2]])) stop("time levels do not contain identical subject IDs")
  df <- df |>
    mutate(time = factor(time, levels = time_levels)) |>
    arrange(subject_id, time)
  wide <- df |>
    select(subject_id, time, value) |>
    pivot_wider(names_from = time, values_from = value) |>
    arrange(subject_id)
  if (anyNA(wide)) stop("incomplete subject pairs")

  tests <- df |>
    pairwise_wilcox_test(value ~ time, paired = TRUE, p.adjust.method = adjust_method) |>
    add_xy_position(x = "time")
  expected <- wilcox.test(wide[[time_levels[[1]]]], wide[[time_levels[[2]]]], paired = TRUE)$p.value
  if (!isTRUE(all.equal(tests$p, expected, tolerance = 1e-8))) stop("paired p-value mismatch")

  plot <- ggpaired(df, x = "time", y = "value", id = "subject_id",
                   color = "time", line.color = "#888888", line.size = 0.25) +
    stat_pvalue_manual(tests, label = "p.adj", tip.length = 0.01) +
    labs(caption = sprintf("Wilcoxon signed-rank; %s-adjusted.", adjust_method)) +
    theme(legend.position = "none")
  ggsave(output_plot, plot, width = 5.5, height = 5, dpi = 160)
  results_csv <- sub("(\\.[^.]+)?$", ".results.csv", output_plot)
  write.csv(tests[vapply(tests, function(column) !is.list(column), logical(1))],
            results_csv, row.names = FALSE)
  if (!file.exists(output_plot) || file.info(output_plot)$size < 1000) stop("plot was not written")
  invisible(list(data = df, tests = tests, plot = plot, results_csv = results_csv))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L) stop("usage: annotate_paired.R INPUT.csv OUTPUT.png [adjust_method]")
  result <- annotate_paired(args[[1]], args[[2]], if (length(args) >= 3L) args[[3]] else "holm")
  cat("wrote", args[[2]], "and", result$results_csv, "\n")
  invisible(NULL)
}
