# Single-cell ATAC-seq usage guide

## Prerequisites

```r
# Signac ecosystem
BiocManager::install(c('GenomicRanges', 'EnsDb.Hsapiens.v86', 'BSgenome.Hsapiens.UCSC.hg38'))
install.packages(c('Signac', 'Seurat', 'leidenbase', 'magrittr'))   # leidenbase: FindClusters(algorithm=4); magrittr: %>%

# ArchR ecosystem
remotes::install_github('GreenleafLab/ArchR', ref='master', repos=BiocManager::repositories())
ArchR::installExtraPackages()
```

```bash
# SnapATAC2 (Python)
pip install snapatac2 macs3   # macs3 also serves Signac CallPeaks and ArchR: pass its path (`which macs3`) as macs2.path / pathToMacs2

# Doublet detection: AMULET is a GitHub release (github.com/UcarLab/AMULET; needs numpy<1.24, pandas, scipy, statsmodels, Java 8+ for BAM input);
# environment and commands are in ecosystem-workflows.md. Alternative: amulet() in the scDblFinder Bioconductor package
```

```r
BiocManager::install('scDblFinder')
```

Inputs (from 10X Cell Ranger ATAC or Multiome):
- `outs/fragments.tsv.gz` (and `.tbi` index)
- `outs/filtered_peak_bc_matrix.h5` (Signac path)
- `outs/singlecell.csv` (ATAC 1.x/2.x) or `per_barcode_metrics.csv` (ARC) (per-barcode metadata)

## Example requests

- "Process 10X scATAC with Signac: NucleosomeSignal and TSSEnrichment; filter with the SKILL.md filter rule; then RunTFIDF -> RunSVD -> RunUMAP/FindNeighbors (dims=2:30) -> FindClusters(algorithm=4, resolution=0.5)."
- "150K-cell dataset: use ArchR (createArrowFiles minTSS=4, minFrags=1000; addDoubletScores; filterDoublets; addIterativeLSI on TileMatrix; addClusters; addUMAP; addReproduciblePeakSet per cluster)."
- "Use SnapATAC2: import_fragments with hg38 chrom_sizes, tsse, filter min_counts=1000 and min_tsse=4, add_tile_matrix bin_size=500, select_features 250k, spectral, leiden, umap."
- "Measure median valid read pairs per cell, pick the doublet tool by depth, and report the flagging rate per cluster."
- "Aggregate cells per cluster into pseudobulk BAMs, run MACS3 per cluster, build an iterative-overlap consensus, then DESeq2 between clusters (atac-seq/differential-accessibility)."
- "Multiome: RNA workflow (LogNormalize -> PCA, dims 1:30) and ATAC workflow (TF-IDF -> SVD, dims 2:30), then FindMultiModalNeighbors and RunUMAP nn.name='weighted.nn'."
- "Compute gene activity with Signac::GeneActivity (or ArchR getGeneScore); transfer labels from a reference scRNA-seq atlas via Seurat::FindTransferAnchors."

## Agent sequence

1. Verify input is Cell Ranger ATAC or Multiome output (chemistry version matters).
2. Re-filter cellranger cell calls with the SKILL.md filter rule.
3. Compute per-cell QC (NucleosomeSignal, TSS enrichment, FRiP, mt fraction, blacklist ratio).
4. Run doublet detection with the depth-appropriate tool (SKILL.md); use the intersection of two tools only when precision matters more than sensitivity.
5. Choose ecosystem, run TF-IDF + LSI / spectral skipping the depth-correlated first component, then UMAP + Leiden.
6. Per-cluster pseudobulk peak calling with MACS3, then iterative-overlap consensus.
7. Gene-activity scores for annotation; WNN integration if Multiome; optional trajectory analysis (ArchR getTrajectory, or Cicero for cis-regulatory inference).

## Tips

- Cell Ranger ATAC's cell calling is lenient; always re-filter.
- The population aggregate is what matters statistically, so per-cell thresholds are looser than bulk.
- Gene-activity scores (Signac::GeneActivity / ArchR getGeneScore) are approximate but useful for cross-modality cell-type transfer.
- chromVAR / AddMotifs must run after the peakset is finalized, or motif annotations become stale.
- Plant / non-model organisms need custom EnsDb / TxDb / BSgenome objects; build with `txdbmaker` or `AnnotationDbi`.
