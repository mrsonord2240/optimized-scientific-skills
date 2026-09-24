# Runnable lollipop-map demonstration.
# Usage: Rscript lollipop_phd.R [optional-input.maf] [optional-clinical.tsv]
# With no inputs it creates a small deterministic MAF and clinical table in the
# working directory, then writes three PDFs and TP53_lollipop.html.

suppressPackageStartupMessages({
  library(data.table)
  library(maftools)
  library(trackViewer)
  library(GenomicRanges)
  library(g3viz)
  library(htmlwidgets)
})

args <- commandArgs(trailingOnly = TRUE)
maf_file <- if (length(args) >= 1L) args[[1L]] else "cohort.maf"
clinical_file <- if (length(args) >= 2L) args[[2L]] else "clinical.tsv"

make_demo_input <- function(maf_file, clinical_file) {
  positions <- c(rep("p.R175H", 8), rep("p.R248Q", 6), rep("p.R273H", 4), "p.Arg175His", "p.*394Wext*?", "p.M1?", "p.=")
  classes <- c(rep("Missense_Mutation", 19), "Nonstop_Mutation", "Translation_Start_Site", "Silent")
  n <- length(positions)
  maf <- data.table(Hugo_Symbol = "TP53", Entrez_Gene_Id = 7157L, Center = "demo", NCBI_Build = "GRCh38", Chromosome = "17", Start_Position = seq_len(n), End_Position = seq_len(n), Strand = "+", Variant_Classification = classes, Variant_Type = "SNP", Reference_Allele = "C", Tumor_Seq_Allele1 = "C", Tumor_Seq_Allele2 = "T", Tumor_Sample_Barcode = sprintf("DEMO-%03d", seq_len(n)), HGVSp_Short = positions)
  fwrite(maf, maf_file, sep = intToUtf8(9L))
  fwrite(data.table(Tumor_Sample_Barcode = maf$Tumor_Sample_Barcode, Subtype = rep(c("Luminal", "Basal"), length.out = n)), clinical_file, sep = intToUtf8(9L))
}

if (!file.exists(maf_file)) make_demo_input(maf_file, clinical_file)
if (!file.exists(clinical_file)) stop("Clinical table not found: ", clinical_file)
clinical <- fread(clinical_file)
maf <- read.maf(maf = maf_file, clinicalData = clinical, verbose = FALSE)
change_col <- intersect(c("HGVSp_Short", "Protein_Change", "AAChange"), names(maf@data))
if (!length(change_col)) stop("No supported protein-change column in MAF")
change_col <- change_col[[1L]]

class_col <- c(Missense_Mutation = "#D55E00", Nonsense_Mutation = "#000000", Frame_Shift_Del = "#0072B2", Frame_Shift_Ins = "#56B4E9", Splice_Site = "#CC79A7", In_Frame_Del = "#009E73", In_Frame_Ins = "#F0E442", Translation_Start_Site = "#999999", Nonstop_Mutation = "#E69F00", Other = "#666666")

# maftools uses its bundled CDD-domain table. refSeqID documents the intended TP53 transcript; printCount prints its summary table to the console.
pdf("TP53_lollipop.pdf", width = 8, height = 4)
lollipopPlot(maf = maf, gene = "TP53", AACol = change_col, refSeqID = "NM_000546", labelPos = c(175, 248, 273), labPosSize = 1, showMutationRate = TRUE, domainLabelSize = 0.8, printCount = TRUE, colors = class_col[names(class_col) != "Other"])
dev.off()

maf_luminal <- subsetMaf(maf, clinQuery = 'Subtype == "Luminal"')
maf_basal <- subsetMaf(maf, clinQuery = 'Subtype == "Basal"')
pdf("TP53_lollipop_subtype.pdf", width = 10, height = 5)
lollipopPlot2(m1 = maf_luminal, m2 = maf_basal, gene = "TP53", m1_name = "Luminal", m2_name = "Basal", AACol1 = change_col, AACol2 = change_col)
dev.off()

protein_position <- function(change) {
  x <- trimws(as.character(change))
  take <- function(pattern, x) { hits <- regmatches(x, regexec(pattern, x)); vapply(hits, function(h) if (length(h) >= 2L) as.integer(h[[2L]]) else NA_integer_, integer(1)) }
  out <- take("^p\\.[A-Z*](\\d+)", x)
  out[grepl("^p\\.[A-Z][0-9]+\\?$", x)] <- NA_integer_ # p.M1? is uncertain, unlike p.*394Wext*?
  miss <- is.na(out)
  out[miss] <- take("^p\\.[A-Z][a-z]{2}(\\d+)", x[miss])
  out
}

calls <- as.data.table(maf@data)[Hugo_Symbol == "TP53"]
calls[, aa_pos := protein_position(get(change_col))]
bad <- calls[is.na(aa_pos)]
if (nrow(bad)) warning("Dropping ", nrow(bad), " call(s) without a parseable protein position: ", paste(unique(bad[[change_col]]), collapse = ", "))
calls <- calls[!is.na(aa_pos) & aa_pos >= 1L & aa_pos <= 393L]
calls[, plot_class := fifelse(as.character(Variant_Classification) %in% names(class_col), as.character(Variant_Classification), "Other")]
summary <- calls[, .(count = .N, class = plot_class[1L], residue = unique(sub("^p\\.", "", get(change_col)))[1L]), by = aa_pos]
snps <- GRanges("TP53", IRanges(summary$aa_pos, width = 1L), color = unname(class_col[as.character(summary$class)]), score = summary$count)
names(snps) <- ifelse(summary$count >= 2L, summary$residue, "")
# The 102--292 DNA-binding core is a conventional supplied range; the other
# three ranges are current UniProt P04637 features.
features <- GRanges("TP53", IRanges(c(1, 102, 325, 368), c(44, 292, 356, 387), names = c("Transactivation", "DNA binding core (conventional)", "Oligomerization", "Basic")), fill = c("#56B4E9", "#0072B2", "#009E73", "#CC79A7"), height = 0.04)
pdf("TP53_trackviewer.pdf", width = 10, height = 4)
lolliplot(snps, features, ylab = "Mutation-row count", xaxis = TRUE, yaxis = TRUE, legend = list(labels = names(class_col), col = unname(class_col)))
dev.off()

# g3viz reads the MAF file rather than a maftools object. g3viz is archived on CRAN; this route requires an already-installed compatible copy.
g3_mutations <- readMAF(maf_file, protein.change.col = change_col)
widget <- g3Lollipop(g3_mutations, gene.symbol = "TP53", protein.change.col = change_col, plot.options = g3Lollipop.theme(theme.name = "nature"), output.filename = "TP53_lollipop")
saveWidget(widget, "TP53_lollipop.html", selfcontained = TRUE)

required_outputs <- c("TP53_lollipop.pdf", "TP53_lollipop_subtype.pdf", "TP53_trackviewer.pdf", "TP53_lollipop.html")
stopifnot(all(file.exists(required_outputs)), all(file.info(required_outputs)$size > 0L))
message("Wrote: ", paste(required_outputs, collapse = ", "))
