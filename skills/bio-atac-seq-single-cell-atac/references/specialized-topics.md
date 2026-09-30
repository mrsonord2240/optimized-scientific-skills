# Specialized topics

## cellranger-atac vs cellranger-arc

| Pipeline | Use for | Output |
|----------|---------|--------|
| cellranger-atac | Single-modality 10X scATAC (chemistry v1, v2) | per-cell ATAC barcodes; fragments.tsv.gz; per-barcode metadata |
| cellranger-arc | 10X Multiome (paired RNA + ATAC, same cell) | joint barcodes for RNA + ATAC; separate fragment / count files; one barcode whitelist |

**Trigger:** Loading 10X output without checking which pipeline produced it.

**Mechanism:** The two produce different output structures. cellranger-arc fragments carry barcodes paired with the RNA matrix; cellranger-atac fragments are ATAC-only with their own barcode universe.

**Fix:** Verify the chemistry on the cellranger summary (look for "Multiome" in the run config). Use `Read10X_h5` for cellranger-arc Multiome RNA output and `CreateChromatinAssay` with the matched fragments for the ATAC. In `per_barcode_metrics.csv` the first column `barcode` equals `gex_barcode` and is the barcode used by the matrix, the BAM `CB` tag and `atac_fragments.tsv.gz`; `atac_barcode` is the raw ATAC barcode and matches none of them. Mixing barcodes across pipelines fails silently.

## Cell-cycle correction

**Trigger:** Proliferating cell types; cells spread across cell-cycle phases.

**Mechanism:** Replication-associated chromatin changes vary with cell-cycle phase and can confound clustering and DA (effect size not established here).

**Detection:** In Multiome data, score S phase from the paired RNA with Seurat `CellCycleScoring` (`cc.genes.updated.2019`). ATAC-only data needs a chromatin proxy (for example a Repli-seq peak overlap score); none is validated here.

**Fix:** `ScaleData(obj, vars.to.regress=...)` cannot change the LSI: `RunSVD` reads the TF-IDF `data` layer, not `scale.data`. Regress the score out of the LSI embedding instead (component 1 is depth and is skipped anyway; on 10x Multiome PBMC data this took max |cor(LSI dim, S score)| from 0.26 to 0):

```r
lsi <- Embeddings(obj, 'lsi')
lsi[, -1] <- resid(lm(lsi[, -1] ~ obj$S.Score))
obj[['lsi_sreg']] <- CreateDimReducObject(lsi, key='LSIS_', assay='ATAC')
obj <- RunUMAP(obj, reduction='lsi_sreg', dims=2:30)
```

For DA between cell-cycle-mismatched conditions, add S-phase as a covariate in pseudobulk DESeq2.

## Sex-chromosome QC

**Trigger:** Mixed-sex donors; sample-swap detection.

**Mechanism:** XIST locus accessibility is high in female cells (X-inactivation); chrY peak count is essentially zero in females. The XIST/chrY ratio identifies sex per cell or sample.

**Detection:** In a per-cell counts matrix, compute the fraction of fragments at the XIST locus (chrX:73820651-73852753 hg38) and the chrY peak count; classify cells; flag cells or samples whose assignment disagrees with metadata.

**XCI escapees:** KDM6A, DDX3X, and EIF1AX escape X inactivation, so females have two accessible copies; males carry Y gametologs (UTY, DDX3Y), so these loci are not a clean sex readout. Prefer XIST and chrY accessibility, and treat escapee accessibility as supporting evidence only.

## scArches reference mapping

**Trigger:** Projecting query scATAC onto a reference atlas; cross-study integration without batch effects.

**Mechanism:** scArches (Lotfollahi 2022) projects a new dataset onto an existing reference's latent space without retraining the reference. For ATAC the relevant model is PEAKVI (Ashuach 2022) in `scvi-tools` (`scvi.model.PEAKVI`); its `load_query_data` classmethod implements the scArches algorithm directly, so `import scarches` is not required.

```python
import scvi

# Pre-trained reference model (e.g. PBMC scATAC atlas)
# The query must be quantified on the reference peak set (same var_names, same order). A different
# feature count is rejected; same count with different names only warns and degrades the mapping (shuffled-peak-order test: median query-to-reference distance 0.099 vs 0.034 aligned).
adata_query = adata_query[:, adata_ref.var_names].copy()   # adata_ref: the reference training AnnData
assert list(adata_query.var_names) == list(adata_ref.var_names)
query_model = scvi.model.PEAKVI.load_query_data(adata_query, reference_path)
query_model.train(max_epochs=200)
adata_query.obsm['X_emb'] = query_model.get_latent_representation()
```

Reference atlases for ATAC are still emerging; the most developed are PBMC (Granja 2021) and brain (BRAIN Initiative).

## chromBPNet per-cluster pseudobulk

**Trigger:** Cell-type-specific variant effect prediction; per-cluster bias-corrected calling.

**Mechanism:** Train chromBPNet (`atac-seq/deep-learning-atac`) per pseudobulk cluster; outputs are bias-corrected per-base profiles and cell-type-specific variant effect predictions.

**Workflow:** Aggregate fragments per cluster into pseudobulk BAMs; run the chromBPNet pipeline per cluster (~24 h GPU each); use the model for in silico variant scoring at GWAS / rare-variant SNPs. Reserve this for the top 5-10 priority clusters; running it on >20 cell types is computationally heavy.

## Reconciliation across tools

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Signac UMAP tight clusters; ArchR UMAP loose | Different LSI implementation; ArchR uses iterative LSI | Both valid; cluster-label biology should match |
| Tools call different doublet rates | Different algorithms (collision vs simulation) | Cells flagged by 2+ tools are a high-precision set (lower sensitivity) |
| Cluster boundaries differ | Different clustering algorithm or resolution | Standardize on Leiden (algorithm 4) with the same resolution |
| Per-cluster peak count differs | Different peak callers or pseudobulk depth | Use the same MACS3 parameters; pool small clusters |

For high-confidence cell-type annotation, agree across two ecosystems (e.g. Signac + ArchR clusters) and report agreement metrics.

## Common errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| UMAP shows depth gradient | LSI component 1 included | `dims=2:30` instead of `1:30` |
| Many low-quality cells in Cell Ranger output | Lenient cell calling | Apply the filter rule in SKILL.md |
| ArchR "TileMatrix not found" | Missed `addTileMat=TRUE` in createArrowFiles | Re-create Arrow files with the flag |
| Signac CreateChromatinAssay fails on fragments | Wrong path to fragments.tsv.gz, or missing tabix index | Provide the full path; run `tabix -p bed fragments.tsv.gz` |
| Weak or unstable MACS3 peaks from a small cluster | Too few reads (a 67-cell cluster still ran without error) | Aggregate small clusters (~200 cells is a rule of thumb) |
| AMULET reports 100% doublets | Threshold mis-set, or input is technical replicates | Check the fragment-count distribution; AMULET needs depth (recall about 0.85-0.90 near 25K valid read pairs/cell) |
| AMULET `AttributeError: module 'numpy' has no attribute 'object'` | numpy >= 1.24 | Run AMULET in an env with `numpy<1.24` |
| Signac `Annotation genome does not match genome of the object` | EnsDb seqlevels are Ensembl style and have no genome | `seqlevelsStyle(ann) <- 'UCSC'; genome(ann) <- 'hg38'` before `CreateChromatinAssay` |
| Multiome WNN clusters dominated by ATAC noise | Equal modality weighting | Inspect modality weights; adjust if needed |
| chromVAR / motif assay all NA | Run before the peakset was finalized | Re-run AddMotifs / RunChromVAR after peaks are stable |
| EnsDb / BSgenome version mismatch | hg38 BSgenome with wrong-build EnsDb | Match builds; `EnsDb.Hsapiens.v86` is GRCh38 (Ensembl v86); use `EnsDb.Hsapiens.v75` for hg19. Newer hg38 EnsDb releases (v98+) reflect newer GENCODE annotations |
