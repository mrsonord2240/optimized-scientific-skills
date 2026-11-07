# Multi-sample handoffs
~~~bash
Rscript scripts/seurat_workflow.R filtered_feature_bc_matrix results_seurat
~~~
Done when samples have been processed independently, integrated before clustering, annotated with marker evidence and assessed for biological conservation; this command performs only the single-sample portion. Demultiplex pools first; ambient correction requires raw droplets and precedes per-capture doublets and normalization. Multi-sample integration, pseudobulk and differential abundance require the separate component workflows; the source snippets are not a complete file-based pipeline. Aggregate raw counts by sample and cell type for condition DE and separately test differential abundance. Do not use integrated values or cells as independent condition replicates.
