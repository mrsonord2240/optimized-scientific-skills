#!/usr/bin/env Rscript
# gtf, rMATS fromGTF directory, STAR SJ directory, Seurat RDS, intron counts,
# gene features, TPM, output RDS, optional event type (default SE), optional read length (required RI).
args <- commandArgs(TRUE)
if (length(args) < 8L || length(args) > 10L) stop('usage: marvel_psi.R gtf rmats_dir star_dir cells.rds intron_counts.tsv gene_features.tsv tpm.tsv out.rds [event_type] [read_length]')
suppressPackageStartupMessages({library(MARVEL); library(data.table); library(Seurat)})
event <- if (length(args) >= 9L) args[9] else 'SE'
events <- c('SE','A5SS','A3SS','MXE','RI')
if (!event %in% events) stop('event type must be SE, A5SS, A3SS, MXE or RI')
read_length <- if (length(args) >= 10L) as.integer(args[10]) else NA_integer_
if (event == 'RI' && (is.na(read_length) || read_length < 1L)) stop('RI requires the actual read length')
gtf <- fread(args[1], header=FALSE, sep='\t', quote='', data.table=FALSE)
features <- lapply(events, function(ev) {
  path <- file.path(args[2], sprintf('fromGTF.%s.txt', ev))
  Preprocess_rMATS(read.table(path, header=TRUE, sep='\t'), GTF=gtf, EventType=ev)
})
names(features) <- events
sj_files <- list.files(args[3], pattern='SJ.out.tab$', full.names=TRUE)
if (!length(sj_files)) stop('no STAR SJ.out.tab files found')
sj_long <- rbindlist(lapply(sj_files, function(f) {
  x <- fread(f, sep='\t', header=FALSE,
    col.names=c('chr','start','end','strand','motif','annot','unique','multi','overhang'))
  x$coord.intron <- paste(x$chr, x$start, x$end, sep=':')
  x$sample <- sub('_SJ.out.tab$', '', basename(f))
  x[, .(coord.intron, sample, unique)]
}))
sj <- dcast(sj_long, coord.intron ~ sample, value.var='unique', fill=0)
pheno <- readRDS(args[4])@meta.data
pheno$sample.id <- rownames(pheno)
marvel <- CreateMarvelObject(SpliceJunction=sj, SplicePheno=pheno, SpliceFeature=features,
  IntronCounts=read.table(args[5], header=TRUE, sep='\t'),
  GeneFeature=read.table(args[6], header=TRUE, sep='\t'),
  Exp=read.table(args[7], header=TRUE, sep='\t'), GTF=gtf)
marvel <- CheckAlignment(marvel, level='SJ')
if (event == 'RI') {
  marvel <- ComputePSI(marvel, CoverageThreshold=10, EventType=event, thread=2, read.length=read_length)
} else {
  marvel <- ComputePSI(marvel, CoverageThreshold=10, EventType=event)
}
saveRDS(marvel, args[8])
