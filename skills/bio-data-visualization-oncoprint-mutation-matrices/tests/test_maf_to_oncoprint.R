# Regression test for scripts/maf_to_oncoprint.R using the audit's planted cohort.
# Usage: r.sh tests/test_maf_to_oncoprint.R synth_edge.maf synth_edge_clin.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: test_maf_to_oncoprint.R <maf> <clinical.tsv>")

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
skill_dir <- normalizePath(file.path(dirname(script_arg), ".."), mustWork = TRUE)
source(file.path(skill_dir, "scripts", "maf_to_oncoprint.R"))

maf <- read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE)
clinical <- read.delim(args[2], stringsAsFactors = FALSE, check.names = FALSE)
mat <- maf_to_oncoprint(maf, clinical$Tumor_Sample_Barcode, top = 6)
aligned <- align_oncoprint_clinical(
  clinical[sample(nrow(clinical)), , drop = FALSE], colnames(mat)
)

stopifnot(ncol(mat) == nrow(clinical))
stopifnot(identical(aligned$Tumor_Sample_Barcode, colnames(mat)))
stopifnot(any(colSums(mat != "") == 0L))
stopifnot(all(diff(rowSums(mat != "")) <= 0L))
stopifnot(all(unlist(strsplit(mat[mat != ""], ";", fixed = TRUE)) %in% ONCOPRINT_CLASSES))

empty_error <- try(maf_to_oncoprint(maf[0, ], clinical$Tumor_Sample_Barcode,
                                    genes = "TP53"), silent = TRUE)
stopifnot(inherits(empty_error, "try-error"))
cat("PASS: cohort-complete matrix, class map, alignment, frequency order, and empty guard\n")
