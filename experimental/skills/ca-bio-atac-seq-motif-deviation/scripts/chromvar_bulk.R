#!/usr/bin/env Rscript
# Purpose: chromVAR on a bulk peak-by-sample count matrix: GC-bias, matched-background z-scores, variability,
#   and (with a sample table) limma differential motifs, in one run.
# Usage:   Rscript chromvar_bulk.R peaks.bed counts.tsv depth.tsv outdir [key=value ...]
#            samples=<tsv>     header line, then sample<TAB>condition; enables the limma contrast
#            ref=<level>       REQUIRED with samples=: the reference (first) condition level
#            genome=BSgenome.Hsapiens.UCSC.hg38   BSgenome package matching the peaks
#            pfm_rds=<rds>     your own PFMatrixList (CIS-BP, plant, ...) instead of JASPAR2024 vertebrates CORE
#            seed=2024  niterations=50  min_logfc=0.5  alpha=0.05
# Input:   peaks.bed (chr, start0, end); counts.tsv (peaks x samples, header = sample names, first column = peak id);
#          depth.tsv (header line, then sample<TAB>total reads or fragments in the library, same unit as the counts;
#          from samtools/featureCounts totals, NOT colSums(counts)).
# Output:  outdir/chromvar_deviations.csv (z-scores, motifs labelled "ID (name)"), chromvar_variability.csv,
#          chromvar_differential.csv (with samples=), chromvar_variability.pdf, chromvar_heatmap.pdf, run.log
#          (seed, niterations, genome, motif source, package versions).
# Stops:   below 5000 peaks after filtering.
# Checked: chromVAR 1.30.1, motifmatchr 1.30.0, JASPAR2024 0.99.6, TFBSTools 1.44.0, limma 3.62.1 (Bioconductor 3.20);
#          not yet executed after the recut.
suppressPackageStartupMessages({
  library(chromVAR); library(motifmatchr); library(SummarizedExperiment)
  library(GenomicRanges); library(limma); library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
opts <- list(samples = '', ref = '', genome = 'BSgenome.Hsapiens.UCSC.hg38', pfm_rds = '', seed = '2024',
             niterations = '50', min_logfc = '0.5', alpha = '0.05')
positional <- character(0)
for (a in args) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) >= 2 && kv[1] %in% names(opts)) opts[[kv[1]]] <- paste(kv[-1], collapse = '=')
  else if (grepl('=', a, fixed = TRUE)) stop('unknown argument ', a)
  else positional <- c(positional, a)
}
if (length(positional) != 4) stop('usage: chromvar_bulk.R peaks.bed counts.tsv depth.tsv outdir [key=value ...]')
if (nzchar(opts$samples) && !nzchar(opts$ref)) stop('samples= needs ref=<level>: the reference level is never guessed')
dir.create(positional[4], showWarnings = FALSE, recursive = TRUE)
outdir <- positional[4]

if (!requireNamespace(opts$genome, quietly = TRUE)) stop('Install ', opts$genome, ' first')
genome <- getExportedValue(opts$genome, opts$genome)

peaks <- read.table(positional[1], col.names = c('chr', 'start', 'end'))[, 1:3]
peak_ranges <- GRanges(seqnames = peaks$chr, ranges = IRanges(peaks$start + 1, peaks$end))
counts_matrix <- as.matrix(read.delim(positional[2], row.names = 1))
if (nrow(counts_matrix) != length(peak_ranges)) stop('peaks.bed has ', length(peak_ranges), ' rows but counts.tsv has ', nrow(counts_matrix))
se <- SummarizedExperiment(assays = list(counts = counts_matrix), rowRanges = peak_ranges)

depth_tab <- read.delim(positional[3], row.names = 1)
depth <- setNames(depth_tab[[1]], rownames(depth_tab))
stopifnot(all(colnames(se) %in% names(depth)))
colData(se)$depth <- depth[colnames(se)]  # filterSamples needs total library size per sample

if (nzchar(opts$samples)) {
  st <- read.delim(opts$samples, row.names = 1)
  if (!all(colnames(se) %in% rownames(st))) stop('samples= table lacks some count-matrix samples')
  cond <- as.character(st[colnames(se), 1])
  if (!opts$ref %in% cond) stop('ref="', opts$ref, '" is not a level of the condition column: ', paste(unique(cond), collapse = ', '))
  colData(se)$condition <- relevel(factor(cond), ref = opts$ref)
}

se <- addGCBias(se, genome = genome)
se <- filterSamples(se, min_depth = 1500, min_in_peaks = 0.15, shiny = FALSE)
se <- filterPeaks(se, non_overlapping = TRUE, min_fragments_per_peak = 10)
cat(sprintf('After filtering: %d peaks, %d samples\n', nrow(se), ncol(se)))
stopifnot(nrow(se) >= 5000)

if (nzchar(opts$pfm_rds)) {
  pfm <- readRDS(opts$pfm_rds); motif_source <- paste('custom', opts$pfm_rds)
} else {
  suppressPackageStartupMessages({library(JASPAR2024); library(TFBSTools); library(RSQLite)})
  jaspar2024 <- JASPAR2024::JASPAR2024()
  sq <- dbConnect(SQLite(), db(jaspar2024))  # TFBSTools issue #39: getMatrixSet does not dispatch on the JASPAR2024 object
  pfm <- getMatrixSet(sq, opts = list(collection = 'CORE', tax_group = 'vertebrates'))
  motif_source <- 'JASPAR2024 vertebrates CORE'
}
motif_ix <- matchMotifs(pfm, se, genome = genome, p.cutoff = 5e-5)

set.seed(as.integer(opts$seed))  # background sampling is random
bg <- getBackgroundPeaks(object = se, niterations = as.integer(opts$niterations))
dev <- computeDeviations(object = se, annotations = motif_ix, background_peaks = bg)
zscores <- deviationScores(dev)
variability <- computeVariability(dev)
variability <- variability[order(-variability$variability), ]
labels <- setNames(paste0(rownames(variability), ' (', variability$name, ')'), rownames(variability))
rownames(zscores) <- labels[rownames(zscores)]
cat('\nTop 10 variable motifs:\n'); print(head(variability, 10))

diff_motifs <- NULL
if (nzchar(opts$samples)) {
  design <- model.matrix(~ colData(se)$condition)
  fit <- eBayes(lmFit(zscores, design))
  diff_motifs <- topTable(fit, coef = 2, number = Inf, sort.by = 'P')  # adj.P.Val is BH; limma has no FDR column
  sig <- diff_motifs[diff_motifs$adj.P.Val < as.numeric(opts$alpha) & abs(diff_motifs$logFC) >= as.numeric(opts$min_logfc), ]
  cat(sprintf('\nDifferential motifs (adj.P.Val < %s, |logFC| >= %s): %d\n', opts$alpha, opts$min_logfc, nrow(sig)))
  write.csv(diff_motifs, file.path(outdir, 'chromvar_differential.csv'))
}

pdf(file.path(outdir, 'chromvar_variability.pdf'), 8, 6); print(plotVariability(variability, use_plotly = FALSE)); dev.off()
top50 <- labels[head(rownames(variability), 50)]
ann <- if ('condition' %in% names(colData(se))) data.frame(Condition = colData(se)$condition, row.names = colnames(zscores)) else NA
pdf(file.path(outdir, 'chromvar_heatmap.pdf'), 10, 12)
pheatmap(zscores[top50, ], annotation_col = ann, scale = 'row', clustering_method = 'ward.D2', fontsize_row = 8)
dev.off()

write.csv(as.data.frame(zscores), file.path(outdir, 'chromvar_deviations.csv'))
write.csv(variability, file.path(outdir, 'chromvar_variability.csv'))  # row names are motif IDs; the other files use "ID (name)"
writeLines(c(paste('seed', opts$seed), paste('niterations', opts$niterations), paste('genome', opts$genome),
             paste('motifs', motif_source), paste('chromVAR', as.character(packageVersion('chromVAR'))),
             paste('motifmatchr', as.character(packageVersion('motifmatchr'))),
             sprintf('%d peaks, %d samples after filtering', nrow(se), ncol(se))), file.path(outdir, 'run.log'))
cat('\nResults in', outdir, '\n')
