# Purpose: compute a declared pairwise family, adjust it, and draw adjusted brackets.
# Inputs: CSV with group,value; output image; mode wilcox|dunn|tukey; p.adjust method;
# optional bracket spacing as a fraction of the observed value range. Use the literal label
# tukey in Tukey mode; Tukey HSD performs its own simultaneous adjustment.
# Usage: Rscript scripts/annotate_pairwise.R data.csv annotated.png wilcox holm [bracket_spacing]

annotate_pairwise <- function(input_csv, output_plot, mode = "wilcox", adjust_method = "holm",
                              bracket_spacing = NULL) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(rstatix)
  })

  df <- read.csv(input_csv, check.names = FALSE)
  stopifnot(all(c("group", "value") %in% names(df)))
  if (anyNA(df[, c("group", "value")])) stop("group and value must not contain missing values")
  df$group <- factor(df$group, levels = unique(df$group))
  if (nlevels(df$group) < 2L) stop("at least two groups are required")

  pairs <- combn(levels(df$group), 2L, simplify = FALSE)
  mode <- match.arg(mode, c("wilcox", "dunn", "tukey"))
  if (mode == "wilcox") {
    tests <- df |>
      pairwise_wilcox_test(value ~ group, comparisons = pairs,
                           p.adjust.method = adjust_method)
    expected <- p.adjust(tests$p, method = adjust_method)
    if (!isTRUE(all.equal(tests$p.adj, expected, tolerance = 1e-8))) {
      stop("rstatix adjusted p-values do not match stats::p.adjust")
    }
    overall_method <- "kruskal.test"
    adjustment_label <- adjust_method
  } else if (mode == "dunn") {
    tests <- df |> dunn_test(value ~ group, p.adjust.method = adjust_method)
    overall_method <- "kruskal.test"
    adjustment_label <- adjust_method
  } else {
    if (tolower(adjust_method) != "tukey") stop("Tukey mode requires adjustment label 'tukey'")
    tests <- df |> tukey_hsd(value ~ group)
    overall_method <- "anova"
    adjustment_label <- "Tukey HSD"
  }

  span <- diff(range(df$value))
  if (!is.finite(span) || span == 0) span <- 1
  if (is.null(bracket_spacing)) {
    # A three-comparison family stays compact; larger families receive progressively
    # more vertical separation, capped to avoid an excessively tall bracket pyramid.
    bracket_spacing <- min(0.22, max(0.10, 0.06 + 0.02 * nrow(tests)))
  }
  if (length(bracket_spacing) != 1L || !is.finite(bracket_spacing) || bracket_spacing <= 0) {
    stop("bracket_spacing must be one positive finite number")
  }
  tests <- tests |> add_xy_position(x = "group", step.increase = bracket_spacing)

  rank_biserial <- function(left, right) {
    x <- df$value[df$group == left]
    y <- df$value[df$group == right]
    u_left <- unname(wilcox.test(x, y, exact = FALSE, correct = FALSE)$statistic)
    2 * u_left / (length(x) * length(y)) - 1
  }
  hedges_g <- function(left, right) {
    x <- df$value[df$group == left]
    y <- df$value[df$group == right]
    degrees_freedom <- length(x) + length(y) - 2
    pooled_sd <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) /
                      degrees_freedom)
    if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
    correction <- 1 - 3 / (4 * degrees_freedom - 1)
    correction * (mean(x) - mean(y)) / pooled_sd
  }
  tests$test <- if (mode == "tukey") "Tukey HSD" else if (mode == "dunn") "Dunn" else "Wilcoxon rank-sum"
  tests$adjust_method <- adjustment_label
  tests$family_size <- nrow(tests)
  tests$effect_type <- if (mode == "tukey") "hedges_g" else "rank_biserial_r"
  tests$effect_size <- vapply(seq_len(nrow(tests)), function(i) {
    if (mode == "tukey") hedges_g(tests$group1[[i]], tests$group2[[i]]) else
      rank_biserial(tests$group1[[i]], tests$group2[[i]])
  }, numeric(1))

  overall_y <- max(tests$y.position) + 0.12 * span
  plot_upper <- overall_y + max(0.12, bracket_spacing) * span
  plot_lower <- min(df$value) - 0.05 * span

  plot <- ggboxplot(df, x = "group", y = "value", color = "group",
                    add = "jitter", outlier.shape = NA) +
    stat_pvalue_manual(tests, label = "p.adj.signif", tip.length = 0.01) +
    stat_compare_means(method = overall_method, label = "p.format", label.y = overall_y) +
    labs(x = NULL,
         caption = sprintf("%s post-hoc; %s family of %d comparisons.",
                           mode, adjustment_label, nrow(tests))) +
    coord_cartesian(ylim = c(plot_lower, plot_upper), clip = "off") +
    theme_classic() +
    theme(legend.position = "none")

  ggsave(output_plot, plot, width = 6.5, height = 5.5, dpi = 160)
  results_csv <- sub("(\\.[^.]+)?$", ".results.csv", output_plot)
  write.csv(tests[vapply(tests, function(column) !is.list(column), logical(1))],
            results_csv, row.names = FALSE)
  if (!file.exists(output_plot) || file.info(output_plot)$size < 1000) stop("plot was not written")
  invisible(list(data = df, pairs = pairs, tests = tests, plot = plot,
                 results_csv = results_csv, bracket_spacing = bracket_spacing,
                 overall_y = overall_y, plot_upper = plot_upper))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L) {
    stop("usage: annotate_pairwise.R INPUT.csv OUTPUT.png [wilcox|dunn|tukey] [adjust_method] [bracket_spacing]")
  }
  result <- annotate_pairwise(args[[1]], args[[2]],
                              if (length(args) >= 3L) args[[3]] else "wilcox",
                              if (length(args) >= 4L) args[[4]] else "holm",
                              if (length(args) >= 5L) as.numeric(args[[5]]) else NULL)
  cat("wrote", args[[2]], "and", result$results_csv, "\n")
  invisible(NULL)
}
