# Regression test for the documented maftools quick path and its denominator caveat.
# Usage: r.sh tests/test_maftools_oncoplot.R synth_edge.maf synth_edge_clin.tsv output.png
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: test_maftools_oncoplot.R <maf> <clinical.tsv> <output.png>")

suppressPackageStartupMessages(library(maftools))
clinical <- read.delim(args[2], stringsAsFactors = FALSE, check.names = FALSE)
maf <- read.maf(args[1], clinicalData = clinical, verbose = FALSE)
stopifnot(nrow(getClinicalData(maf)) < nrow(clinical))

png(args[3], width = 1400, height = 700, res = 120)
oncoplot(maf, top = 20,
         clinicalFeatures = c("Subtype", "Stage"),
         annotationColor = list(
           Subtype = c(Luminal = "#0072B2", Basal = "#D55E00", HER2 = "#009E73"),
           Stage = c(I = "#FFFFCC", II = "#FED976", III = "#FD8D3C", IV = "#BD0026")),
         sortByAnnotation = TRUE, removeNonMutated = FALSE)
dev.off()
stopifnot(file.info(args[3])$size > 1000)
cat("PASS: complete annotation colors; maftools kept",
    nrow(getClinicalData(maf)), "of", nrow(clinical), "cohort samples\n")
