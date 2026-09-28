# Focused regression checks for adjusted labels, identity-safe pairing, and model-based annotation.
# Usage: Rscript tests/regression.R

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
test_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
skill_dir <- dirname(test_dir)
source(file.path(skill_dir, "scripts", "annotate_pairwise.R"))
source(file.path(skill_dir, "scripts", "annotate_paired.R"))
source(file.path(skill_dir, "scripts", "annotate_nested.R"))

scratch <- file.path(tempdir(), "statanno-regression")
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
data_dir <- file.path(skill_dir, "examples", "data")

pairwise <- annotate_pairwise(file.path(data_dir, "three_group.csv"),
                              file.path(scratch, "pairwise.png"), "wilcox", "holm")
stopifnot(isTRUE(all.equal(pairwise$tests$p.adj,
                          p.adjust(pairwise$tests$p, method = "holm"), tolerance = 1e-8)))
stopifnot(identical(as.character(pairwise$tests$p.adj.signif), c("ns", "ns", "*")))
stopifnot(all(c("test", "adjust_method", "family_size", "effect_type", "effect_size") %in%
              names(pairwise$tests)))
stopifnot(all(pairwise$tests$effect_type == "rank_biserial_r"),
          all(is.finite(pairwise$tests$effect_size)),
          all(abs(pairwise$tests$effect_size) <= 1))
persisted_pairwise <- read.csv(pairwise$results_csv, check.names = FALSE)
stopifnot(all(c("effect_type", "effect_size") %in% names(persisted_pairwise)))

set.seed(12)
dense <- data.frame(
  group = rep(LETTERS[1:4], each = 12),
  value = c(rnorm(12, 0), rnorm(12, 0.2), rnorm(12, 0.7), rnorm(12, 1.1))
)
dense_path <- file.path(scratch, "dense.csv")
write.csv(dense, dense_path, row.names = FALSE)
dense_result <- annotate_pairwise(dense_path, file.path(scratch, "dense.png"),
                                  "wilcox", "holm")
stopifnot(nrow(dense_result$tests) == 6L,
          dense_result$bracket_spacing > pairwise$bracket_spacing,
          dense_result$overall_y > max(dense_result$tests$y.position),
          dense_result$plot_upper > dense_result$overall_y)
custom_spacing <- annotate_pairwise(dense_path, file.path(scratch, "dense-custom.png"),
                                    "wilcox", "holm", 0.2)
stopifnot(identical(custom_spacing$bracket_spacing, 0.2))
dunn_result <- annotate_pairwise(dense_path, file.path(scratch, "dense-dunn.png"),
                                 "dunn", "holm")
stopifnot(all(dunn_result$tests$effect_type == "rank_biserial_r"),
          all(is.finite(dunn_result$tests$effect_size)))
tukey_result <- annotate_pairwise(dense_path, file.path(scratch, "dense-tukey.png"),
                                  "tukey", "tukey")
stopifnot(all(tukey_result$tests$effect_type == "hedges_g"),
          all(is.finite(tukey_result$tests$effect_size)))

paired <- annotate_paired(file.path(data_dir, "paired.csv"),
                          file.path(scratch, "paired.png"), "holm")
stopifnot(paired$tests$p.adj[[1]] < 0.01)

incomplete <- read.csv(file.path(data_dir, "paired.csv"))[-1, ]
incomplete_path <- file.path(scratch, "paired-incomplete.csv")
write.csv(incomplete, incomplete_path, row.names = FALSE)
rejected <- inherits(try(annotate_paired(incomplete_path,
                                         file.path(scratch, "should-not-exist.png")), silent = TRUE),
                     "try-error")
stopifnot(rejected)

nested <- annotate_nested(file.path(data_dir, "nested.csv"),
                          file.path(scratch, "nested.png"), "holm")
stopifnot(nrow(nested$tests) == 1L, is.finite(nested$tests$p.adj[[1]]))
message("R regression checks passed")
