# Purpose: fit a subject-level random-intercept model and draw model-based adjusted contrasts.
# Inputs: CSV with group,subject_id,value; output image; p.adjust method.
# Usage: Rscript scripts/annotate_nested.R nested.csv nested.png holm

annotate_nested <- function(input_csv, output_plot, adjust_method = "holm") {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(lme4)
    library(emmeans)
    library(rstatix)
  })

  df <- read.csv(input_csv, check.names = FALSE)
  stopifnot(all(c("group", "subject_id", "value") %in% names(df)))
  if (anyNA(df[, c("group", "subject_id", "value")])) stop("nested fields must not be missing")
  memberships <- df |> distinct(subject_id, group) |> count(subject_id)
  if (any(memberships$n != 1L)) stop("each subject_id must belong to exactly one group")
  df$group <- factor(df$group, levels = unique(df$group))
  if (nlevels(df$group) < 2L) stop("at least two groups are required")

  model <- lmer(value ~ group + (1 | subject_id), data = df)
  emm <- emmeans(model, ~ group)
  contrast_table <- as.data.frame(pairs(emm, adjust = adjust_method))
  split_contrast <- strsplit(as.character(contrast_table$contrast), " - ", fixed = TRUE)
  if (any(lengths(split_contrast) != 2L)) stop("could not parse emmeans contrast labels")
  tests <- data.frame(
    group1 = vapply(split_contrast, `[[`, character(1), 1L),
    group2 = vapply(split_contrast, `[[`, character(1), 2L),
    p.adj = contrast_table$p.value,
    stringsAsFactors = FALSE
  ) |>
    add_significance(p.col = "p.adj", output.col = "p.adj.signif") |>
    add_y_position(formula = value ~ group, data = df, step.increase = 0.12)

  plot <- ggplot(df, aes(group, value, color = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.25) +
    stat_pvalue_manual(tests, label = "p.adj", tip.length = 0.01) +
    labs(x = NULL, caption = sprintf("Random-intercept LMM; emmeans %s-adjusted contrasts.",
                                     adjust_method)) +
    theme_classic() +
    theme(legend.position = "none")
  ggsave(output_plot, plot, width = 5.5, height = 5, dpi = 160)
  results_csv <- sub("(\\.[^.]+)?$", ".results.csv", output_plot)
  write.csv(tests[vapply(tests, function(column) !is.list(column), logical(1))],
            results_csv, row.names = FALSE)
  if (!file.exists(output_plot) || file.info(output_plot)$size < 1000) stop("plot was not written")
  invisible(list(data = df, model = model, emmeans = emm, tests = tests,
                 plot = plot, results_csv = results_csv))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L) stop("usage: annotate_nested.R INPUT.csv OUTPUT.png [adjust_method]")
  result <- annotate_nested(args[[1]], args[[2]], if (length(args) >= 3L) args[[3]] else "holm")
  cat("wrote", args[[2]], "and", result$results_csv, "\n")
  invisible(NULL)
}
