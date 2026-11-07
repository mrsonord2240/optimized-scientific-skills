args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript seurat_markers.R input.rds output.csv")
library(Seurat)
obj <- readRDS(args[1])
if (length(unique(Idents(obj))) < 2) stop("Need at least two clusters")
markers <- FindAllMarkers(obj, only.pos = TRUE, test.use = "wilcox",
  logfc.threshold = 0.25, min.pct = 0.1)
write.csv(markers, args[2], row.names = FALSE)
