#!/usr/bin/env Rscript
# Reference: ATACseqQC 1.26+, GenomicAlignments 1.38+ | Verify API if version differs
# Usage: Rscript atac_qc_metrics.R <bam> [peaks.narrowPeak|-] [output_prefix] [TxDb_package]
# Paired-end, coordinate-sorted BAM. Counts only MAPQ >= 30, non-duplicate-flagged, non-chrM,
# properly paired alignments (mark duplicates beforehand). Writes <prefix>_fragsize.pdf and
# <prefix>_metrics.csv. TxDb_package defaults to the hg38 knownGene package; its seqlevel style
# must match the BAM (stops otherwise).

suppressPackageStartupMessages({
    library(ATACseqQC); library(GenomicAlignments); library(GenomicRanges); library(rtracklayer)
})

# Periodicity from the smoothed fragment-size density: a peak counts when its density is >= 1.15x
# the lowest density between it and the previous (shorter) peak window. Heuristic, not ENCODE-defined.
classify_periodicity <- function(frag_widths) {
    d <- density(frag_widths[frag_widths <= 800], from = 20, to = 800, bw = 8, n = 781)
    peak_in <- function(lo, hi) {
        i <- which(d$x >= lo & d$x <= hi); i[which.max(d$y[i])]
    }
    win <- list(nfr = c(30, 120), mono = c(150, 260), di = c(300, 480))
    idx <- vapply(win, function(w) peak_in(w[1], w[2]), numeric(1))
    present <- logical(3); names(present) <- names(win)
    present['nfr'] <- TRUE
    for (k in 2:3) {
        valley <- min(d$y[idx[k - 1]:idx[k]])
        present[k] <- d$y[idx[k]] >= 1.15 * valley
    }
    nfr_present <- d$y[idx['nfr']] >= 1.15 * min(d$y[idx['nfr']:idx['mono']])
    if (!nfr_present || !present['mono']) return('flat/single')
    if (present['di']) return('3+ peaks')
    'NFR+mono'
}

read_peaks <- function(peaks_file) {
    tryCatch(import(peaks_file, format = 'narrowPeak'),
             error = function(e) tryCatch(import(peaks_file, format = 'BED'),
                 error = function(e2) stop(sprintf('cannot import peaks %s as narrowPeak or BED: %s',
                                                   peaks_file, conditionMessage(e2)))))
}

calculate_atac_qc <- function(bam_file, peaks_file = NULL, output_prefix = 'atac_qc',
                              txdb_pkg = 'TxDb.Hsapiens.UCSC.hg38.knownGene') {
    cat('=== ENCODE 4 ATAC-seq QC ===\n')
    suppressPackageStartupMessages(library(txdb_pkg, character.only = TRUE))
    txs <- transcripts(get(txdb_pkg))

    pdf(sprintf('%s_fragsize.pdf', output_prefix))
    fragSizeDist(bam_file, output_prefix)
    dev.off()

    param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isProperPair = TRUE),
                          mapqFilter = 30)
    gal <- readGAlignmentPairs(bam_file, param = param)
    gal_se <- readGAlignments(bam_file, param = param)
    if (length(intersect(seqlevels(gal), seqlevels(txs))) == 0)
        stop(sprintf('no shared seqlevels between BAM (%s...) and %s (%s...)',
                     paste(head(seqlevels(gal), 3), collapse = ','), txdb_pkg,
                     paste(head(seqlevels(txs), 3), collapse = ',')))
    mt <- c('chrM', 'chrMT', 'MT', 'M')
    gal <- gal[!as.character(seqnames(gal)) %in% mt]
    gal_se <- gal_se[!as.character(seqnames(gal_se)) %in% mt]
    n_pairs <- length(gal)
    cat(sprintf('Nuclear pairs (MAPQ >= 30, non-duplicate, non-chrM): %d\n', n_pairs))

    tsse <- TSSEscore(gal_se, txs)
    cat(sprintf('ATACseqQC TSSEscore: %.2f\n', tsse$TSSEscore))
    cat('  (ATACseqQC uses 100bp center / 1000bp flanks; not comparable to ENCODE-style TSS enrichment)\n')

    frags <- width(granges(gal))
    nfr <- sum(frags < 100); mono <- sum(frags >= 180 & frags <= 247)
    di <- sum(frags >= 315 & frags <= 473); tri <- sum(frags >= 558 & frags <= 615)
    periodicity <- classify_periodicity(frags)
    cat(sprintf('NFR: %d (%.1f%%)  Mono: %d (%.1f%%)  Di: %d (%.1f%%)  Tri: %d\n',
                nfr, 100 * nfr / n_pairs, mono, 100 * mono / n_pairs, di, 100 * di / n_pairs, tri))
    cat(sprintf('Periodicity: %s\n', periodicity))

    frip <- NA_real_
    if (!is.null(peaks_file) && nzchar(peaks_file) && peaks_file != '-' && file.exists(peaks_file)) {
        frip <- sum(countOverlaps(gal, read_peaks(peaks_file)) > 0) / n_pairs
        cat(sprintf('FRiP (pairs overlapping a peak): %.3f  (ENCODE: ideal >= 0.3, accept >= 0.2)\n', frip))
    }

    grade <- function(value, accept, ideal) {
        if (is.na(value)) 'NA' else if (value < accept) 'FAIL' else if (value >= ideal) 'PASS' else 'WARN'
    }
    metrics <- data.frame(
        sample = output_prefix,
        nuclear_reads_M = round(2 * n_pairs / 1e6, 2),
        ATACseqQC_TSSEscore = round(tsse$TSSEscore, 2),
        nfr_fraction = round(nfr / n_pairs, 3),
        mono_fraction = round(mono / n_pairs, 3),
        periodicity = periodicity,
        FRiP = ifelse(is.na(frip), NA, round(frip, 3)),
        nuc_reads_grade = grade(2 * n_pairs / 1e6, 25, 50),
        FRiP_grade = grade(frip, 0.2, 0.3)
    )
    write.csv(metrics, sprintf('%s_metrics.csv', output_prefix), row.names = FALSE)
    cat('\nReport card:\n'); print(metrics)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) calculate_atac_qc(args[1], if (length(args) > 1) args[2] else NULL,
                                        if (length(args) > 2) args[3] else 'sample',
                                        if (length(args) > 3) args[4] else 'TxDb.Hsapiens.UCSC.hg38.knownGene')
