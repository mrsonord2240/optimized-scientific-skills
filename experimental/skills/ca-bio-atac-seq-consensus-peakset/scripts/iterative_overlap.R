#!/usr/bin/env Rscript
# Purpose: Corces-style iterative overlap consensus in R (GenomicRanges): pool narrowPeak files, re-center on
#   summits, fix width to 2*half_width+1, keep peaks greedily by descending signalValue when they overlap no
#   kept peak, drop blacklist overlaps. The R counterpart of iterative_overlap.sh.
# Usage:   Rscript iterative_overlap.R chrom.sizes blacklist.bed out.bed half_width peaks1.narrowPeak [peaks2 ...]
#            blacklist.bed may be .gz; pass NONE to skip the blacklist filter.
#            half_width 250 gives the 501 bp standard.
# Output:  out.bed (chrom, start, end), sorted. The greedy loop is quadratic: above ~500k pooled peaks use
#          iterative_overlap.sh.
# Needs:   GenomicRanges, rtracklayer
# Checked: not yet executed after the recut.
suppressPackageStartupMessages({library(GenomicRanges); library(rtracklayer)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop('usage: iterative_overlap.R chrom.sizes blacklist.bed|NONE out.bed half_width peaks1.narrowPeak [peaks2 ...]')
sizes_file <- args[1]; blacklist_file <- args[2]; out_file <- args[3]
half_width <- as.integer(args[4]); peak_files <- args[-(1:4)]

# Validate every record before arithmetic: narrowPeak has 10 columns and -1 in column 10 means no summit.
for (f in peak_files) {
  rows <- Filter(function(l) nzchar(trimws(l)) && !startsWith(l, '#'), readLines(f))
  if (!length(rows)) stop(f, ': no non-comment narrowPeak records found')
  for (i in seq_along(rows)) {
    p <- strsplit(rows[i], '\t', fixed = TRUE)[[1]]
    if (length(p) != 10) stop(f, ' record ', i, ': expected 10 tab-delimited narrowPeak columns, found ', length(p))
    s <- suppressWarnings(as.numeric(p[2])); e <- suppressWarnings(as.numeric(p[3])); o <- suppressWarnings(as.numeric(p[10]))
    if (is.na(s) || is.na(e) || e <= s) stop(f, ' record ', i, ': invalid start/end')
    if (isTRUE(o == -1)) stop(f, ' record ', i, ': summit offset is -1 (no summit was called)')
    if (is.na(o) || o < 0 || o >= e - s) stop(f, ' record ', i, ': summit offset must be within the interval')
  }
}

peaks <- unlist(GRangesList(lapply(peak_files, import, format = 'narrowPeak')))
chrom_sizes <- read.table(sizes_file, col.names = c('chrom', 'len'))
lens <- setNames(chrom_sizes$len, chrom_sizes$chrom)
peaks <- peaks[as.character(seqnames(peaks)) %in% names(lens)]
seqlevels(peaks, pruning.mode = 'coarse') <- intersect(seqlevels(peaks), names(lens))
seqlengths(peaks) <- lens[seqlevels(peaks)]

summit <- start(peaks) + peaks$peak          # narrowPeak start is 0-based; import() shifts to 1-based, offset is unchanged
centered <- peaks
ranges(centered) <- IRanges(start = summit, width = 1)
fixed <- trim(resize(centered, width = 2 * half_width + 1, fix = 'center'))

sorted <- fixed[order(-fixed$signalValue)]
kept <- logical(length(sorted))
used <- GRanges(seqinfo = seqinfo(sorted))
for (i in seq_along(sorted)) {
  if (!overlapsAny(sorted[i], used)) { kept[i] <- TRUE; used <- c(used, sorted[i]) }
}
consensus <- sort(sorted[kept])
mcols(consensus) <- NULL
if (blacklist_file != 'NONE') {
  consensus <- consensus[!overlapsAny(consensus, import(blacklist_file))]
}
export(consensus, out_file, format = 'bed')
message('wrote ', out_file, ' (', length(consensus), ' peaks at ', 2 * half_width + 1, ' bp)')
