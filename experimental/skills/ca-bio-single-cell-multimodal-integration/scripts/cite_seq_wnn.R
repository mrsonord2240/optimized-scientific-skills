#!/usr/bin/env Rscript
# Purpose: CITE-seq in Seurat: QC, DSB denoising of ADT from the empty droplets, PCA per modality, WNN joint
#   graph, clusters and UMAP; writes the per-cell modality weights.
# Usage:   Rscript cite_seq_wnn.R raw_feature_bc_matrix/ filtered_feature_bc_matrix/ out_prefix [key=value ...]
#            isotypes=<name1,name2>   ADT row names of the isotype controls. Default: rows matching 'isotype|IgG'
#                                     (case-insensitive on Isotype); if none match, DSB runs without isotype controls.
#            resolution=0.5  adt_dims=18  rna_dims=30
# Input:   10x unfiltered (raw) and filtered directories with 'Gene Expression' and 'Antibody Capture'.
#          The raw matrix must contain the empty droplets; the script stops if the empty-droplet median total ADT
#          count is at least 0.5x the cell median (DSB gives no error when given cells as "empties").
# Output:  out_prefix.rds (Seurat object with ADT assay = DSB values, pca, apca, wnn.umap, wsnn clusters),
#          out_prefix_weights.csv (RNA.weight, ADT.weight, cluster per cell), out_prefix_umap.pdf.
# Checked: Seurat 5.x, dsb 2.0.1 (not yet executed after the recut).
suppressPackageStartupMessages({library(Seurat); library(dsb)})

args <- commandArgs(trailingOnly = TRUE)
opts <- list(isotypes = '', resolution = '0.5', adt_dims = '18', rna_dims = '30')
positional <- character(0)
for (a in args) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) >= 2 && kv[1] %in% names(opts)) opts[[kv[1]]] <- paste(kv[-1], collapse = '=')
  else if (grepl('=', a, fixed = TRUE)) stop('unknown argument ', a)
  else positional <- c(positional, a)
}
if (length(positional) != 3) stop('usage: cite_seq_wnn.R raw_dir filtered_dir out_prefix [isotypes=...] [resolution=0.5]')
out <- positional[3]

raw <- Read10X(positional[1])        # unfiltered: contains empty droplets, needed for DSB
data <- Read10X(positional[2])       # called cells
rna_counts <- data$`Gene Expression`
adt_counts <- data$`Antibody Capture`
obj <- CreateSeuratObject(counts = rna_counts, assay = 'RNA', min.cells = 3)
adt_cells <- as.matrix(adt_counts[, colnames(obj)])
adt_empty <- as.matrix(raw$`Antibody Capture`[, setdiff(colnames(raw$`Antibody Capture`), colnames(adt_counts))])

obj <- PercentageFeatureSet(obj, pattern = '^MT-', col.name = 'percent.mt')
obj <- subset(obj, nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
adt_cells <- adt_cells[, colnames(obj)]

obj <- NormalizeData(obj, assay = 'RNA')
obj <- FindVariableFeatures(obj, assay = 'RNA')
obj <- ScaleData(obj, assay = 'RNA')
obj <- RunPCA(obj, assay = 'RNA', reduction.name = 'pca')

med_cells <- median(colSums(adt_cells)); med_empty <- median(colSums(adt_empty))
if (med_empty >= med_cells * 0.5)
  stop(sprintf('empty_drop_matrix does not look like empty droplets (median total ADT %.1f vs cells %.1f)', med_empty, med_cells))

iso <- if (nzchar(opts$isotypes)) strsplit(opts$isotypes, ',', fixed = TRUE)[[1]] else grep('[Ii]sotype|IgG', rownames(adt_cells), value = TRUE)
use_iso <- length(iso) > 0
message(if (use_iso) paste('isotype controls:', paste(iso, collapse = ', ')) else 'no isotype controls found: DSB runs with use.isotype.control = FALSE')
adt_dsb <- if (use_iso) {
  DSBNormalizeProtein(cell_protein_matrix = adt_cells, empty_drop_matrix = adt_empty, denoise.counts = TRUE,
                      use.isotype.control = TRUE, isotype.control.name.vec = iso)
} else {
  DSBNormalizeProtein(cell_protein_matrix = adt_cells, empty_drop_matrix = adt_empty, denoise.counts = TRUE,
                      use.isotype.control = FALSE)
}

obj[['ADT']] <- CreateAssayObject(counts = adt_cells)
obj <- SetAssayData(obj, assay = 'ADT', layer = 'data', new.data = adt_dsb)  # DSB output is already normalized data
obj <- ScaleData(obj, assay = 'ADT')
adt_dims <- min(as.integer(opts$adt_dims), nrow(obj[['ADT']]) - 1)
obj <- RunPCA(obj, assay = 'ADT', reduction.name = 'apca', features = rownames(obj[['ADT']]), npcs = adt_dims)

obj <- FindMultiModalNeighbors(obj, reduction.list = list('pca', 'apca'),
                               dims.list = list(1:as.integer(opts$rna_dims), 1:min(adt_dims, ncol(obj[['apca']]))))
obj <- FindClusters(obj, graph.name = 'wsnn', resolution = as.numeric(opts$resolution))
obj <- RunUMAP(obj, nn.name = 'weighted.nn', reduction.name = 'wnn.umap')

write.csv(data.frame(cell = colnames(obj), RNA.weight = obj$RNA.weight, ADT.weight = obj$ADT.weight,
                     cluster = obj$seurat_clusters), paste0(out, '_weights.csv'), row.names = FALSE)
cat('Mean modality weights: RNA', mean(obj$RNA.weight), ' ADT', mean(obj$ADT.weight), '\n')
pdf(paste0(out, '_umap.pdf'), width = 10, height = 8); print(DimPlot(obj, reduction = 'wnn.umap', label = TRUE)); dev.off()
saveRDS(obj, paste0(out, '.rds'))
message('wrote ', out, '.rds')
