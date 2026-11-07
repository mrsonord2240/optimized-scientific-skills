# Purpose: render a cohort-complete ComplexHeatmap OncoPrint and run pairwise interactions.
# Inputs: MAF, clinical TSV with Tumor_Sample_Barcode/Subtype/Stage, and output PDF or PNG.
# Usage: Rscript examples/oncoprint_phd.R cohort.maf clinical.tsv oncoprint.pdf

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: oncoprint_phd.R <cohort.maf> <clinical.tsv> <output.pdf|output.png>")
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(maftools)
})

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
skill_dir <- normalizePath(file.path(dirname(script_arg), ".."), mustWork = TRUE)
source(file.path(skill_dir, "scripts", "maf_to_oncoprint.R"))

maf_path <- args[1]
clinical <- read.delim(args[2], stringsAsFactors = FALSE, check.names = FALSE)
.require_columns(clinical, c("Tumor_Sample_Barcode", "Subtype", "Stage"), "clinical")
maf_df <- read.delim(maf_path, comment.char = "#", stringsAsFactors = FALSE,
                     check.names = FALSE, quote = "")

mat <- maf_to_oncoprint(maf_df, clinical$Tumor_Sample_Barcode, top = 20)
clinical <- align_oncoprint_clinical(clinical, colnames(mat))

if (!"tmb" %in% names(clinical)) {
  tmb <- table(maf_df$Tumor_Sample_Barcode)
  clinical$tmb <- as.numeric(tmb[match(clinical$Tumor_Sample_Barcode, names(tmb))])
  clinical$tmb[is.na(clinical$tmb)] <- 0
}
stopifnot(identical(clinical$Tumor_Sample_Barcode, colnames(mat)))

col <- c(
  Missense = "#56B4E9", Truncating = "#000000", Splice = "#CC79A7",
  Amp = "#D55E00", HomDel = "#0072B2", Fusion = "#009E73"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = "#EEEEEE", col = NA))
  },
  Amp = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["Amp"], col = NA))
  },
  HomDel = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["HomDel"], col = NA))
  },
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.5,
              gp = gpar(fill = col["Missense"], col = NA))
  },
  Truncating = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33,
              gp = gpar(fill = col["Truncating"], col = NA))
  },
  Splice = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.25,
              gp = gpar(fill = col["Splice"], col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid.points(x, y, pch = 17, size = unit(2, "mm"),
                gp = gpar(col = col["Fusion"]))
  }
)

subtype_col <- c(Luminal = "#0072B2", Basal = "#D55E00", HER2 = "#009E73")
stage_col <- c(I = "#FFFFCC", II = "#FED976", III = "#FD8D3C", IV = "#BD0026")
missing_subtype <- setdiff(unique(clinical$Subtype), names(subtype_col))
missing_stage <- setdiff(unique(clinical$Stage), names(stage_col))
if (length(missing_subtype) || length(missing_stage)) {
  stop("Add colors for every Subtype and Stage level before plotting.")
}

ha_top <- HeatmapAnnotation(
  TMB = anno_barplot(
    log10(clinical$tmb + 1), gp = gpar(fill = "#56B4E9", col = NA),
    axis_param = list(at = log10(c(1, 10, 100, 1000) + 1),
                      labels = c("1", "10", "100", "1000"))
  ),
  Subtype = clinical$Subtype,
  Stage = clinical$Stage,
  col = list(Subtype = subtype_col, Stage = stage_col),
  annotation_name_gp = gpar(fontsize = 8),
  show_legend = TRUE
)

ht <- oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col,
  top_annotation = ha_top,
  row_order = order(-rowSums(mat != "")),
  column_title = sprintf("Mutation landscape (N=%d)", ncol(mat)),
  row_names_gp = gpar(fontsize = 8),
  pct_gp = gpar(fontsize = 7),
  show_pct = TRUE,
  remove_empty_columns = FALSE,
  remove_empty_rows = FALSE,
  heatmap_legend_param = list(title = "Alteration", at = names(col), labels = names(col))
)

output_ext <- tolower(tools::file_ext(args[3]))
if (output_ext == "pdf") {
  pdf(args[3], width = 12, height = 7)
} else if (output_ext == "png") {
  png(args[3], width = 3600, height = 2100, res = 300)
} else {
  stop("Output extension must be .pdf or .png.")
}
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# maftools uses only samples represented by MAF rows; this result is not the
# denominator source for the cohort-complete plot above.
maf_object <- read.maf(maf = maf_path, clinicalData = clinical, verbose = FALSE)
interaction_output <- file.path(
  dirname(args[3]),
  paste0(tools::file_path_sans_ext(basename(args[3])), "-interactions.pdf")
)
pdf(interaction_output, width = 8, height = 7)
si <- somaticInteractions(maf = maf_object, top = min(20, nrow(mat)),
                          pvalue = c(0.05, 0.01), fontSize = 0.7)
dev.off()
required_si <- c("gene1", "gene2", "pValue", "oddsRatio", "00", "01",
                 "11", "10", "pAdj", "Event")
stopifnot(all(required_si %in% names(si)))
message("somaticInteractions returned ", nrow(si),
        " pairwise rows; plot: ", interaction_output)
