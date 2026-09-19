# Reference: Seurat 5.0+ | Verify API if version differs
library(Seurat)
library(dsb)

raw <- Read10X('raw_feature_bc_matrix/')           # unfiltered: contains empty droplets, needed for DSB
data <- Read10X('filtered_feature_bc_matrix/')      # called cells
rna_counts <- data$`Gene Expression`
adt_counts <- data$`Antibody Capture`

obj <- CreateSeuratObject(counts = rna_counts, assay = 'RNA', min.cells = 3)
adt_cells <- as.matrix(adt_counts[, colnames(obj)])
adt_empty <- as.matrix(raw$`Antibody Capture`[, setdiff(colnames(raw$`Antibody Capture`), colnames(adt_cells))])

obj <- PercentageFeatureSet(obj, pattern = '^MT-', col.name = 'percent.mt')
obj <- subset(obj, nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
adt_cells <- adt_cells[, colnames(obj)]

obj <- NormalizeData(obj, assay = 'RNA')
obj <- FindVariableFeatures(obj, assay = 'RNA')
obj <- ScaleData(obj, assay = 'RNA')
obj <- RunPCA(obj, assay = 'RNA', reduction.name = 'pca')

# Denoise ADT with DSB before joint embedding (SKILL.md's "CITE-seq: Denoise ADT" section):
# CLR alone rescales but does not remove ambient/technical background, and WNN itself does
# not denoise protein either. Guard against DSB's no-empty-droplets failure mode first --
# see SKILL.md's Common Errors table.
med_cells <- median(colSums(adt_cells))
med_empty <- median(colSums(adt_empty))
if (med_empty >= med_cells * 0.5) {
    stop(sprintf(
        "empty_drop_matrix does not look like empty droplets (median total ADT %.1f vs cells %.1f)",
        med_empty, med_cells))
}
adt_dsb <- DSBNormalizeProtein(
    cell_protein_matrix = adt_cells,
    empty_drop_matrix = adt_empty,
    denoise.counts = TRUE,
    use.isotype.control = TRUE,
    isotype.control.name.vec = grep('[Ii]sotype|IgG', rownames(adt_cells), value = TRUE)
)

obj[['ADT']] <- CreateAssayObject(counts = adt_cells)
obj <- SetAssayData(obj, assay = 'ADT', layer = 'data', new.data = adt_dsb)
obj <- ScaleData(obj, assay = 'ADT')
obj <- RunPCA(obj, assay = 'ADT', reduction.name = 'apca',
              features = rownames(obj[['ADT']]), npcs = min(18, nrow(obj[['ADT']]) - 1))

obj <- FindMultiModalNeighbors(obj, reduction.list = list('pca', 'apca'),
                                dims.list = list(1:30, 1:min(18, ncol(obj[['apca']]))))
obj <- FindClusters(obj, graph.name = 'wsnn', resolution = 0.5)
obj <- RunUMAP(obj, nn.name = 'weighted.nn', reduction.name = 'wnn.umap')

pdf('cite_seq_wnn_umap.pdf', width = 10, height = 8)
print(DimPlot(obj, reduction = 'wnn.umap', label = TRUE))
dev.off()

saveRDS(obj, 'cite_seq_analyzed.rds')
cat('Analysis complete. Saved to cite_seq_analyzed.rds\n')
