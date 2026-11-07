#!/usr/bin/env Rscript
# Reference: Signac 1.13+, Seurat 5.0+, EnsDb.Hsapiens.v86 2.99+, BSgenome.Hsapiens.UCSC.hg38 1.4+ | Verify API if version differs
# Standard Signac scATAC-seq pipeline: 10X Cell Ranger output -> per-cell QC -> TF-IDF/LSI ->
# UMAP/Leiden (skipping depth-correlated dim 1) -> gene-activity scores for annotation.

suppressPackageStartupMessages({
    library(Signac); library(Seurat); library(EnsDb.Hsapiens.v86)
    library(BSgenome.Hsapiens.UCSC.hg38); library(GenomicRanges); library(ggplot2)
})

# QC rule set = the per-cell QC table in SKILL.md. Doublet detection is NOT run here.
run_signac <- function(h5_file='outs/filtered_peak_bc_matrix.h5',
                       fragments_file='outs/fragments.tsv.gz',
                       metadata_file='outs/singlecell.csv',
                       min_tss=4, min_frags=1000, output_prefix='scatac') {

    stopifnot(requireNamespace('leidenbase', quietly=TRUE))  # FindClusters(algorithm=4)
    counts <- Read10X_h5(h5_file)
    if (is.list(counts)) counts <- counts[['Peaks']]  # cellranger-arc matrix holds 'Gene Expression' and 'Peaks'
    metadata <- read.csv(metadata_file, header=TRUE, row.names=1)
    arc <- 'atac_fragments' %in% colnames(metadata)    # cellranger-arc per_barcode_metrics.csv
    if (arc) {
        # ARC: first column = GEX barcode = barcode used by the matrix and by atac_fragments.tsv.gz
        need <- c('atac_fragments', 'atac_peak_region_fragments')
        if (!all(need %in% colnames(metadata))) stop('ARC metadata lacks columns: ', paste(setdiff(need, colnames(metadata)), collapse=', '))
        metadata <- metadata[colnames(counts), , drop=FALSE]
        if (anyNA(rownames(metadata))) stop('matrix barcodes are missing from the ARC metadata (wrong sample pairing?)')
        metadata$passed_filters <- metadata$atac_fragments
        metadata$peak_region_fragments <- metadata$atac_peak_region_fragments
        if (all(c('atac_mitochondrial_reads', 'atac_raw_reads') %in% colnames(metadata))) {
            metadata$mitochondrial <- metadata$atac_mitochondrial_reads; metadata$total <- metadata$atac_raw_reads
        }
        cat('Metadata mode: cellranger-arc (blacklist ratio from the peak matrix; mito = mito reads / raw reads)\n')
    } else {
        need <- c('passed_filters', 'peak_region_fragments', 'blacklist_region_fragments')
        if (!all(need %in% colnames(metadata)))
            stop('metadata file lacks columns: ', paste(setdiff(need, colnames(metadata)), collapse=', '))
    }

    # EnsDb is Ensembl-style (1, 2, X) with no genome; the object is hg38/UCSC (chr1)
    ann <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
    seqlevelsStyle(ann) <- 'UCSC'
    genome(ann) <- 'hg38'

    # Build Signac assay
    chrom_assay <- CreateChromatinAssay(
        counts=counts, sep=c(':', '-'),
        genome='hg38', fragments=fragments_file,
        annotation=ann,
        min.cells=10, min.features=200)
    obj <- CreateSeuratObject(counts=chrom_assay, assay='ATAC', meta.data=metadata)

    # Per-cell QC -- looser per-cell thresholds than bulk
    obj <- NucleosomeSignal(obj)
    obj <- TSSEnrichment(obj, fast=FALSE)
    obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
    obj$blacklist_ratio <- if (arc) FractionCountsInRegion(obj, assay='ATAC', regions=blacklist_hg38_unified)
                           else obj$blacklist_region_fragments / obj$peak_region_fragments
    if ('mitochondrial' %in% colnames(metadata)) obj$mito_fraction <- obj$mitochondrial / obj$total

    # QC plots
    pdf(sprintf('%s_qc.pdf', output_prefix), 14, 4)
    print(VlnPlot(obj, features=c('nCount_ATAC', 'TSS.enrichment',
                                   'pct_reads_in_peaks', 'nucleosome_signal',
                                   'blacklist_ratio'), pt.size=0.1, ncol=5))
    dev.off()

    # Re-filter cellranger output at sensible per-cell thresholds (cellranger is lenient)
    keep <- obj$passed_filters >= min_frags & obj$passed_filters <= 80000 &
            obj$pct_reads_in_peaks >= 15 & obj$blacklist_ratio < 0.05 &
            obj$nucleosome_signal <= 4 & obj$TSS.enrichment >= min_tss
    if ('mito_fraction' %in% colnames(obj[[]])) keep <- keep & obj$mito_fraction < 0.05
    keep[is.na(keep)] <- FALSE
    cat(sprintf('After QC filter: %d / %d cells\n', sum(keep), ncol(obj)))
    if (!any(keep)) stop('No cells pass QC; inspect ', output_prefix, '_qc.pdf and lower min_tss/min_frags (args 4, 5) for shallow or subset data')
    obj_filt <- obj[, keep]

    # Dimensionality reduction
    # CRITICAL: dims = 2:30, NOT 1:30. LSI component 1 is depth, not biology.
    obj_filt <- RunTFIDF(obj_filt)
    obj_filt <- FindTopFeatures(obj_filt, min.cutoff='q0')
    obj_filt <- RunSVD(obj_filt)

    # Confirm depth correlation: component 1 should correlate with nCount_ATAC
    depth_cor <- DepthCor(obj_filt)
    pdf(sprintf('%s_depth_cor.pdf', output_prefix), 6, 4)
    print(depth_cor)
    dev.off()
    cat('  (Component 1 should correlate strongly with depth, |r| near 1, sign arbitrary; that is why we skip it.)\n')

    obj_filt <- RunUMAP(obj_filt, reduction='lsi', dims=2:30)
    obj_filt <- FindNeighbors(obj_filt, reduction='lsi', dims=2:30)
    obj_filt <- FindClusters(obj_filt, algorithm=4, resolution=0.5)        # Leiden = algorithm 4

    pdf(sprintf('%s_umap.pdf', output_prefix), 7, 6)
    print(DimPlot(obj_filt, label=TRUE, label.size=4))
    dev.off()

    # Gene activity (approximation of expression from accessibility) for annotation
    cat('Computing gene activity scores...\n')
    gene_activities <- GeneActivity(obj_filt)
    obj_filt[['ACT']] <- CreateAssayObject(counts=gene_activities)
    DefaultAssay(obj_filt) <- 'ACT'
    obj_filt <- NormalizeData(obj_filt, normalization.method='LogNormalize',
                              scale.factor=median(obj_filt$nCount_ACT))

    saveRDS(obj_filt, sprintf('%s_signac.rds', output_prefix))
    cat(sprintf('Saved to %s_signac.rds\n', output_prefix))
    invisible(obj_filt)
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
    run_signac(args[1], args[2], args[3],
               min_tss=if (length(args) >= 4) as.numeric(args[4]) else 4,
               min_frags=if (length(args) >= 5) as.numeric(args[5]) else 1000)
}
