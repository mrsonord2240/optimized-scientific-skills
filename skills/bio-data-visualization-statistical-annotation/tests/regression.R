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
