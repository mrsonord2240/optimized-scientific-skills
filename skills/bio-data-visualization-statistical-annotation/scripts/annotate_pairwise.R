# Purpose: compute a declared pairwise family, adjust it, and draw adjusted brackets.
# Inputs: CSV with group,value; output image; mode wilcox|dunn|tukey; p.adjust method
# (use the literal label tukey in Tukey mode; Tukey HSD performs its own simultaneous adjustment).
# Usage: Rscript scripts/annotate_pairwise.R data.csv annotated.png wilcox holm

annotate_pairwise <- function(input_csv, output_plot, mode = "wilcox", adjust_method = "holm") {
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

  tests <- tests |> add_xy_position(x = "group", step.increase = 0.1)
  span <- diff(range(df$value))
  if (!is.finite(span) || span == 0) span <- 1
  overall_y <- max(tests$y.position) + 0.12 * span

  plot <- ggboxplot(df, x = "group", y = "value", color = "group",
                    add = "jitter", outlier.shape = NA) +
    stat_pvalue_manual(tests, label = "p.adj.signif", tip.length = 0.01) +
    stat_compare_means(method = overall_method, label = "p.format", label.y = overall_y) +
    labs(x = NULL,
         caption = sprintf("%s post-hoc; %s family of %d comparisons.",
                           mode, adjustment_label, nrow(tests))) +
    theme_classic() +
    theme(legend.position = "none")

  ggsave(output_plot, plot, width = 6.5, height = 5.5, dpi = 160)
  results_csv <- sub("(\\.[^.]+)?$", ".results.csv", output_plot)
  write.csv(tests[vapply(tests, function(column) !is.list(column), logical(1))],
            results_csv, row.names = FALSE)
  if (!file.exists(output_plot) || file.info(output_plot)$size < 1000) stop("plot was not written")
  invisible(list(data = df, pairs = pairs, tests = tests, plot = plot,
                 results_csv = results_csv))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L) {
    stop("usage: annotate_pairwise.R INPUT.csv OUTPUT.png [wilcox|dunn|tukey] [adjust_method]")
  }
  result <- annotate_pairwise(args[[1]], args[[2]],
                              if (length(args) >= 3L) args[[3]] else "wilcox",
                              if (length(args) >= 4L) args[[4]] else "holm")
  cat("wrote", args[[2]], "and", result$results_csv, "\n")
  invisible(NULL)
}
