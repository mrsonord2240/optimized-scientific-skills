#!/usr/bin/env Rscript
# Reference: cicero 1.3.x (GitHub cole-trapnell-lab/cicero-release, branch monocle3), monocle3 1.3+, GenomicRanges 1.54+ | Verify API if version differs
# Cicero co-accessibility: scATAC peak matrix -> metacell aggregation -> graphical lasso ->
# peak-pair connection scores -> enhancer-gene candidate pairs.
# CLI: Rscript cicero_workflow.R peaks.mtx peak_metadata.tsv cell_metadata.tsv [tss.bed]
#   peak_metadata.tsv row names = peak names "chr_start_end" (e.g. chr1_10244_10510), same order as matrix rows;
#   cell_metadata.tsv row names = cell barcodes, same order as matrix columns;
#   tss.bed (optional) = BED with the gene name in column 4 (e.g. GENCODE protein-coding TSS); omit to skip enhancer-gene mapping.

suppressPackageStartupMessages({
    library(cicero); library(monocle3); library(GenomicRanges); library(rtracklayer)
})

run_cicero_pipeline <- function(peak_matrix, peak_metadata, cell_metadata,
                                output_prefix='cicero',
                                k_metacell=50, window=500000, coaccess_threshold=0.25,
                                tss_bed=NULL, tss_pad=2000, seed=1) {

    # 0. Validate inputs; coerce to a numeric dgCMatrix (readMM returns a logical ngTMatrix for a
    #    binary 'pattern' .mtx, which monocle3 rejects).
    peak_re <- '^[^_:]+_[0-9]+_[0-9]+$'
    stopifnot(nrow(peak_matrix) == nrow(peak_metadata), ncol(peak_matrix) == nrow(cell_metadata))
    if (!all(grepl(peak_re, rownames(peak_metadata))))
        stop("peak_metadata row names must be 'chr_start_end' (e.g. chr1_10244_10510)")
    peak_matrix <- as(as(as(peak_matrix, 'dMatrix'), 'generalMatrix'), 'CsparseMatrix')
    dimnames(peak_matrix) <- list(rownames(peak_metadata), rownames(cell_metadata))
    if (!'gene_short_name' %in% colnames(peak_metadata))
        peak_metadata$gene_short_name <- rownames(peak_metadata)
    if (!is.null(tss_bed)) {   # fail before the long run, not after it
        if (!file.exists(tss_bed)) stop(sprintf('tss_bed not found: %s', tss_bed))
        tss <- import(tss_bed)
        if (is.null(tss$name)) stop('tss_bed needs the gene name in column 4')
    }
    set.seed(seed)

    # 1. Build CDS from inputs (peaks x cells matrix; expects 0/1)
    input_cds <- new_cell_data_set(peak_matrix,
                                   cell_metadata=cell_metadata,
                                   gene_metadata=peak_metadata)
    cat(sprintf('Loaded: %d peaks, %d cells\n', nrow(input_cds), ncol(input_cds)))

    # 2. Dimensionality reduction (LSI for sparse binary data; then UMAP)
    input_cds <- detect_genes(input_cds)
    input_cds <- estimate_size_factors(input_cds)
    input_cds <- preprocess_cds(input_cds, method='LSI')
    input_cds <- reduce_dimension(input_cds, reduction_method='UMAP',
                                  preprocess_method='LSI')

    # 3. Build metacells via k-NN (default k=50). Smaller k -> more variability captured but slower.
    umap_coords <- reducedDims(input_cds)$UMAP
    cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates=umap_coords, k=k_metacell)

    # 4. Run Cicero on the chromosomes present in the peak set (genome-build independent; cis only).
    #    Runtime grows with the modelled span: chr1 1-30 Mb, 1,726 peaks x 3,277 cells took 8.5-14 min on 24 cores (varies with machine load).
    peak_gr <- GRanges(sub('_([0-9]+)_([0-9]+)$', ':\\1-\\2', rownames(peak_metadata)))
    chr_end <- tapply(end(peak_gr), as.character(seqnames(peak_gr)), max)
    genome_df <- data.frame(chr=names(chr_end), length=as.numeric(chr_end))

    cat('Running Cicero (this can take time on large datasets)...\n')
    conns <- run_cicero(cicero_cds, genomic_coords=genome_df,
                        window=window, sample_num=100)

    # Cicero returns every pair in both orientations (Peak1 is a character, Peak2 a factor):
    # keep one row per unordered pair and drop unscored (NA) pairs.
    conns$Peak1 <- as.character(conns$Peak1); conns$Peak2 <- as.character(conns$Peak2)
    conns <- conns[!is.na(conns$coaccess), ]
    lo <- pmin(conns$Peak1, conns$Peak2); hi <- pmax(conns$Peak1, conns$Peak2)
    conns$Peak1 <- lo; conns$Peak2 <- hi
    conns <- conns[!duplicated(paste(lo, hi)), ]

    cat(sprintf('Scored peak pairs (unordered): %d (%d negative)\n', nrow(conns), sum(conns$coaccess < 0)))
    strong <- conns[conns$coaccess > coaccess_threshold, ]
    cat(sprintf('Strong (coaccess > %.2f): %d\n', coaccess_threshold, nrow(strong)))

    # 5. Save strong pairs for visualization
    write.table(strong, sprintf('%s_connections.tsv', output_prefix),
                sep='\t', row.names=FALSE, quote=FALSE)

    # 6. Map to enhancer-gene pairs via TSS overlap
    if (is.null(tss_bed)) {
        cat('No tss_bed given; skipping enhancer-gene mapping\n')
    } else {
        tss <- import(tss_bed)
        tss_extended <- resize(tss, width=2 * tss_pad, fix='center')

        # Cicero peaks are chr_start_end; convert to chr:start-end
        peak1 <- GRanges(sub('_([0-9]+)_([0-9]+)$', ':\\1-\\2', strong$Peak1))
        peak2 <- GRanges(sub('_([0-9]+)_([0-9]+)$', ':\\1-\\2', strong$Peak2))
        ov1 <- findOverlaps(peak1, tss_extended)
        ov2 <- findOverlaps(peak2, tss_extended)

        eg_pairs <- data.frame(
            enhancer = c(strong$Peak2[queryHits(ov1)], strong$Peak1[queryHits(ov2)]),
            gene = c(tss_extended$name[subjectHits(ov1)],
                     tss_extended$name[subjectHits(ov2)]),
            coaccess = c(strong$coaccess[queryHits(ov1)], strong$coaccess[queryHits(ov2)]),
            stringsAsFactors=FALSE)
        stopifnot(all(eg_pairs$enhancer %in% rownames(peak_metadata)))
        # An 'enhancer' anchor that itself lies in a TSS window is a promoter-promoter pair.
        enh_gr <- GRanges(sub('_([0-9]+)_([0-9]+)$', ':\\1-\\2', eg_pairs$enhancer))
        eg_pairs$enhancer_is_promoter <- overlapsAny(enh_gr, tss_extended)
        eg_pairs <- eg_pairs[order(-eg_pairs$coaccess), ]
        eg_pairs <- eg_pairs[!duplicated(eg_pairs[c('enhancer', 'gene')]), ]
        write.csv(eg_pairs, sprintf('%s_enhancer_gene_pairs.csv', output_prefix),
                  row.names=FALSE)
        cat(sprintf('Enhancer-gene pairs: %d unique (%d promoter-promoter)\n',
                    nrow(eg_pairs), sum(eg_pairs$enhancer_is_promoter)))
    }

    invisible(list(conns=conns, strong=strong))
}

# Usage: provide peak_matrix.mtx, peak_metadata.tsv, cell_metadata.tsv (+ optional tss.bed) as args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) >= 3) {
    peak_matrix <- Matrix::readMM(args[1])
    peak_meta <- read.delim(args[2], row.names=1)
    cell_meta <- read.delim(args[3], row.names=1)
    run_cicero_pipeline(peak_matrix, peak_meta, cell_meta,
                        tss_bed=if (length(args) >= 4) args[4] else NULL)
}
