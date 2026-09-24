# Purpose: chromVAR motif deviations for a Signac object, then differential motif activity between two groups.
# Input:   RDS of a Seurat object with a ChromatinAssay named 'peaks' (hg38, chr-prefixed seqnames) and cluster identities.
# Output:  <out_prefix>_obj.rds (object with a 'chromvar' assay of background-normalized z-scores) and
#          <out_prefix>_diff_motifs.csv (FindMarkers result, motif names in a 'motif_name' column).
# Usage:   Rscript run_chromvar.R obj.rds out_prefix cluster1 cluster2 [group_by=seurat_clusters]
# Note:    Signac::RunChromVAR() was removed in Signac 1.17.0; this calls chromVAR's own API. Checked on Signac 1.17.1, chromVAR 1.28.0.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)
rds_in <- args[1]; out_prefix <- args[2]; ident1 <- args[3]; ident2 <- args[4]
group_by <- if (length(args) >= 5) args[5] else 'seurat_clusters'

suppressPackageStartupMessages({
  library(Signac); library(Seurat)
  library(JASPAR2020); library(TFBSTools); library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(chromVAR); library(SummarizedExperiment); library(BiocParallel)
})
register(SerialParam())   # chromVAR/motifmatchr default to a multicore backend unsupported on Windows

obj <- readRDS(rds_in)
Idents(obj) <- group_by

pfm <- getMatrixSet(JASPAR2020, opts = list(collection = 'CORE', tax_group = 'vertebrates', all_versions = FALSE))
obj <- AddMotifs(obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

# Signac::RunChromVAR() was removed in Signac 1.17.0 (chromVAR became unavailable in Bioconductor
# 3.23, per Signac's own NEWS.md) -- call chromVAR's own lower-level API directly instead; this is
# the same sequence RunChromVAR used to wrap, and runs on any Signac version.
se <- SummarizedExperiment(assays = list(counts = as.matrix(GetAssayData(obj, assay = 'peaks', layer = 'counts'))),
                            rowRanges = granges(obj[['peaks']]))
se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
motif_ix <- matchMotifs(pfm, se, genome = BSgenome.Hsapiens.UCSC.hg38)
set.seed(1)                                        # getBackgroundPeaks() samples background peaks at
                                                    # random and is NOT internally seeded -- omitting
                                                    # this makes chromVAR's differential-motif calls and
                                                    # rankings change from run to run on identical input
bg_peaks <- getBackgroundPeaks(se)                 # GC- and accessibility-matched background
dev <- computeDeviations(object = se, annotations = motif_ix, background_peaks = bg_peaks)
obj[['chromvar']] <- CreateAssayObject(data = deviationScores(dev))   # background-normalized z-scores

DefaultAssay(obj) <- 'chromvar'
diff_motifs <- FindMarkers(obj, ident.1 = ident1, ident.2 = ident2,
                           mean.fxn = rowMeans, fc.name = 'avg_diff')
diff_motifs$motif_name <- rownames(diff_motifs)

saveRDS(obj, paste0(out_prefix, '_obj.rds'))
write.csv(diff_motifs, paste0(out_prefix, '_diff_motifs.csv'), row.names = FALSE)
cat('Motifs scored:', nrow(deviationScores(dev)), ' cells:', ncol(obj), ' significant (p_val_adj<0.05):', sum(diff_motifs$p_val_adj < 0.05), '\n')
