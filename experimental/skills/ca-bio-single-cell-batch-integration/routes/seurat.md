# Seurat v5 integration
~~~sh
Rscript scripts/seurat_v5_integration.R merged.rds integrated.rds batch RPCAIntegration
~~~
Provide a Seurat object, output path, batch metadata column and available method. Script uses seeded PCA, split layers, integration, neighbors, clustering and UMAP. Optional wrapper methods need their packages loaded and configured; inspect available functions and memory before use. Keep original counts and review [Seurat notes](../references/seurat-v5.md).
**Done when:** output and uncorrected comparisons have been evaluated for mixing and preservation.
