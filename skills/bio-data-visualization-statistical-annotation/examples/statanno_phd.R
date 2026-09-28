# Purpose: run the adjusted pairwise, ID-checked paired, manual ggsignif, and nested workflows.
# Inputs: optional three_group.csv paired.csv nested.csv output_dir; shipped fixtures are defaults.
# Usage: Rscript examples/statanno_phd.R [THREE_GROUP.csv PAIRED.csv NESTED.csv OUTPUT_DIR]

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("could not determine example path")
example_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
skill_dir <- dirname(example_dir)

source(file.path(skill_dir, "scripts", "annotate_pairwise.R"))
source(file.path(skill_dir, "scripts", "annotate_paired.R"))
source(file.path(skill_dir, "scripts", "annotate_nested.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsignif)
})

args <- commandArgs(trailingOnly = TRUE)
three_group_csv <- if (length(args) >= 1L) args[[1]] else file.path(example_dir, "data", "three_group.csv")
paired_csv <- if (length(args) >= 2L) args[[2]] else file.path(example_dir, "data", "paired.csv")
nested_csv <- if (length(args) >= 3L) args[[3]] else file.path(example_dir, "data", "nested.csv")
output_dir <- if (length(args) >= 4L) args[[4]] else file.path(example_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Adjusted Wilcoxon family. Place the omnibus label above the highest bracket.
pairwise_result <- annotate_pairwise(
  three_group_csv, file.path(output_dir, "pairwise-adjusted.png"), "wilcox", "holm"
)
pairwise_result$plot <- pairwise_result$plot +
  labs(caption = sprintf(
    "Wilcoxon pairwise, Holm-adjusted; rank-biserial r range %.2f-%.2f.",
    min(pairwise_result$tests$effect_size), max(pairwise_result$tests$effect_size)
  ))
ggsave(file.path(output_dir, "pairwise-adjusted.png"), pairwise_result$plot,
       width = 6.5, height = 5.5, dpi = 160)

# ggsignif computes raw p-values when test= is used. Supply adjusted labels manually.
manual_signif <- ggplot(pairwise_result$data, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.45) +
  geom_signif(comparisons = pairwise_result$pairs,
              annotations = pairwise_result$tests$p.adj.signif,
              y_position = pairwise_result$tests$y.position,
              tip_length = 0.01) +
  labs(caption = "Manual ggsignif labels from Holm-adjusted rstatix results.") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(output_dir, "ggsignif-adjusted.png"), manual_signif,
       width = 6.5, height = 5.5, dpi = 160)

# Both time levels must contain exactly the same subject IDs before row-order pairing.
annotate_paired(paired_csv, file.path(output_dir, "paired-adjusted.png"), "holm")

# Annotate the mixed-model contrast, never a cell-level rank test.
annotate_nested(nested_csv, file.path(output_dir, "nested-lmm.png"), "holm")

expected <- file.path(output_dir, c(
  "pairwise-adjusted.png", "ggsignif-adjusted.png", "paired-adjusted.png", "nested-lmm.png"
))
stopifnot(all(file.exists(expected)), all(file.info(expected)$size > 1000))
print(expected)
