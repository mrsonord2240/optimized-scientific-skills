## Made-once commitments (decided at `cellranger count`, inherited downstream)

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Reference build + Ensembl vintage (CR v10 = Ensembl v110) | Every gene ID/marker lookup/cross-dataset merge; join on `gene_ids` (stable), not symbols (lossy across releases) |
| `--include-introns` (default TRUE since CR 7) | Total UMI/cell and snRNA-vs-scRNA comparability; nuclei need introns; mixing intron-on/off samples is a hidden batch axis |
| Cell-vs-empty-droplet calling | The filtered matrix is cell-CALLING only, NOT ambient-corrected; SoupX/CellBender and low-RNA-cell rescue need the RAW matrix (filtered-only is irreversible loss) |
| Feature/barcode namespace (`gex_only`) | `gex_only=True` silently drops ADT/HTO/CRISPR features; set False and split by `feature_types` for CITE-seq/hashed pools |

## Pipeline orchestration: the ordering decisions that make or break the result

Stage order is not arbitrary; getting it wrong silently corrupts every downstream step. The cross-cutting decisions this pipeline must get right:
- Correct-then-detect-then-normalize, never reorder. Ambient-RNA removal (SoupX/CellBender) reads the RAW droplet matrix and MUST precede doublet detection and the final QC thresholds (ambient inflation distorts both doublet scores and per-cell QC); doublet detection precedes integration (doublets seed fake intermediate clusters). See bio-single-cell-preprocessing.
- Multiplexed (hashed) pools are demultiplexed FIRST. When samples are pooled with cell hashing (HTO/MULTI-seq) or genotype multiplexing, assign each cell to its sample of origin and remove cross-sample doublets before any per-sample step; the per-sample QC and doublet operations below assume the pool is already split by sample. See bio-single-cell-hashing-demultiplexing.
- Per-sample QC and doublet detection come BEFORE any merge or integration. Ambient-RNA cleanup (SoupX/CellBender), mito/gene-count filtering, and doublet calling (expected rate ~0.8% per 1,000 recovered cells) are per-capture operations; running them on a pooled object lets one lane's artifacts contaminate the shared null and lets surviving doublets form fake intermediate clusters. See bio-single-cell-preprocessing and bio-single-cell-doublet-detection.
- Integrate, THEN cluster, THEN annotate - never reorder. For multi-sample designs, batch-correct on a shared embedding (Harmony/scVI/RPCA) and cluster on the corrected graph; clustering before correction makes clusters track samples or lanes instead of biology, and annotation has nothing to label until clusters exist. See bio-single-cell-batch-integration, bio-single-cell-clustering, bio-single-cell-cell-annotation.
- Over-integration erases biology. Batch-mixing metrics (kBET, iLISI) are maximized by destroying structure, so a method that flattens real cell states scores perfectly; pair every batch metric with a bio-conservation metric and re-inspect rare populations after correction. Batch-corrected expression is for embedding and clustering only - never feed it to differential expression. See bio-single-cell-batch-integration.
- Condition-level DE must use pseudobulk, not cells-as-replicates. Testing thousands of cells per donor as independent replicates is pseudoreplication and inflates false positives by orders of magnitude (Squair 2021); aggregate RAW counts per sample x cell type and hand them to DESeq2/edgeR/limma-voom. Step 7 marker p-values are descriptive ranking for labeling, not condition inference. See bio-single-cell-markers-annotation and bio-differential-expression-deseq2-basics.
- Composition shifts are tested separately and masquerade as DE. A cluster-level expression change between conditions can be a pure proportion shift (a mixed cluster's substates re-balance with no gene changing per cell), invisible if only DE is run; always pair condition DE with a differential-abundance test (Milo cluster-free, or scCODA/sccomp/propeller cluster-based). See bio-single-cell-differential-abundance.

The Seurat and Scanpy paths below are written single-sample for clarity. A multiplexed design demultiplexes the pool first (bio-single-cell-hashing-demultiplexing); a multi-sample design then inserts per-sample QC and doublet removal, a merge, and an integration step between normalization (Step 4) and dimensionality reduction (Step 5); the rest of the order is unchanged.


## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| QC | min.features | 200-500 |
| QC | max.features | 2500-5000 (depends on data) |
| QC | percent.mt | <10-20% (tissue-dependent) |
| SCTransform | vars.to.regress | none by default; only a validated non-confounded covariate |
| PCA | npcs | 30-50 |
| UMAP | dims | 15-30 (check elbow plot) |
| Clustering | resolution | 0.4-0.8 (start with 0.5) |

QC cutoffs above are illustrative starting points, not universal thresholds: set min/max features and mito % per dataset from the data (MAD-adaptive, tissue-aware) rather than porting flat values, since healthy mito fraction varies by tissue. See bio-single-cell-preprocessing.

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| All cells filtered | QC too strict | Relax thresholds |
| Poor UMAP separation | Too few HVGs or PCs | Increase nfeatures, check n_pcs |
| Too many/few clusters | Wrong resolution | Adjust resolution parameter |
| Unknown cell types | Missing markers | Check known marker genes manually |

## Multi-sample completion and checkpoints

For a multi-sample study, apply ambient removal, QC, and doublet removal to each
sample independently, then merge and integrate before neighbors/clustering. Emit
these four checkpoints in every run: after loading, after QC, after normalization,
and after clustering.

~~~r
stopifnot(!any(grepl("^ENSG", rownames(seu)))) # map Ensembl IDs before ^MT- QC
cat("[after_loading] cells=", ncol(seu), "\n")
seu <- subset(seu, subset=nFeature_RNA > median(nFeature_RNA)-3*mad(nFeature_RNA) &
                         nCount_RNA > median(nCount_RNA)-3*mad(nCount_RNA) &
                         percent.mt < min(25, median(percent.mt)+3*mad(percent.mt)))
cat("[after_qc] cells=", ncol(seu), "\n")
seu <- SCTransform(seu, verbose=FALSE)
cat("[after_normalization] cells=", ncol(seu), "\n")
merged <- merge(samples[[1]], y=samples[-1], add.cell.ids=names(samples))
merged <- RunHarmony(merged, group.by.vars="sample_id")
merged <- FindNeighbors(merged, reduction="harmony", dims=1:30)
merged <- FindClusters(merged)
cat("[after_clustering] clusters=", length(unique(merged$seurat_clusters)), "\n")
~~~

Use raw counts for sample-aware inference, never cell-level DE:

~~~r
pb <- AggregateExpression(merged, assays="RNA", group.by=c("sample_id","cell_type"),
                          return.seurat=FALSE)$RNA
# DESeq2: design = ~ condition; run separately per cell type after preserving sample IDs.
# Differential abundance: run propeller(clusters=merged$cell_type,
# sample=merged$sample_id, group=merged$condition) on sample-by-cell-type proportions.
~~~

If ambient correction is intended, start from the raw matrix; filtered matrices
cannot recover empty droplets. Warn before QC if ncol(seu) is implausibly close
to the raw barcode count, and inspect rownames(seu) before applying the MT pattern.


## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Ambient markers everywhere (e.g. Hb across all PBMCs) | Ran SoupX/CellBender on the FILTERED matrix (soup already removed) | Run ambient removal on the RAW droplet matrix, before QC |
| Fake intermediate cell states | Doublets removed AFTER clustering/integration | Call doublets per sample before merge/normalize |
| Clusters track sample/lane, not biology | Clustered before integration | Integrate -> cluster -> annotate |
| Inflated DE, thousands of false positives | Tested cells as replicates (on integrated values) | Pseudobulk RAW counts per sample x cell-type -> DESeq2/edgeR/limma (Squair 2021) |
| "DE change" is really a proportion shift | Only ran DE, not abundance | Pair with differential abundance (Milo/scCODA/propeller) |
| Biology collapses to one blob | Reflexively regressed total_counts/cell-cycle | Regress only a validated, non-confounded covariate |
| CITE-seq/hashtag features vanish | `gex_only=True` default | Set False; split by `feature_types` |
