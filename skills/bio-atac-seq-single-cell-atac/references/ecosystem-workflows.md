# Ecosystem workflows

Code sketches for the non-default ecosystems. The standard Signac pipeline is [`../scripts/signac_workflow.R`](../scripts/signac_workflow.R). Check installed versions and signatures before running (see the core Skill).

## Per-cluster pseudobulk peak calling (Signac + MACS3)

```r
# Signac wrapper around MACS3
peaks <- CallPeaks(obj, group.by='seurat_clusters',
                   macs2.path='/path/to/macs3', cleanup=FALSE,
                   format='BED', shift=-75, extsize=150,   # Tn5 cut-site recipe; shift only applies to format='BED'
                   additional.args='-p 0.01')
# Then iterative-overlap consensus across clusters (atac-seq/consensus-peakset)
```

## ArchR

**Goal:** Build an ArchR project from fragment files, filter doublets, cluster, and call per-cluster reproducible peaks.

```r
library(ArchR)
addArchRGenome('hg38')

# 1. Create Arrow files from fragment files
ArrowFiles <- createArrowFiles(
    inputFiles=c('fragments_rep1.tsv.gz', 'fragments_rep2.tsv.gz'),
    sampleNames=c('rep1', 'rep2'),
    minTSS=4, minFrags=1000,
    addTileMat=TRUE, addGeneScoreMat=TRUE)
# Failure is quiet: if no cell passes, ArchR logs "has encountered an error" and returns an empty vector
stopifnot(length(ArrowFiles) > 0)   # else lower minTSS/minFrags (shallow data) and read the ArchR log

# 2. Doublet detection (built-in)
doubletScores <- addDoubletScores(input=ArrowFiles, k=10, knnMethod='UMAP')

# 3. Project + filter
proj <- ArchRProject(ArrowFiles=ArrowFiles, outputDirectory='ArchR_out')
proj <- filterDoublets(proj)

# 4. LSI + UMAP + clustering
proj <- addIterativeLSI(proj, useMatrix='TileMatrix', name='IterativeLSI')
proj <- addClusters(proj, reducedDims='IterativeLSI', method='Seurat', resolution=0.5)
proj <- addUMAP(proj, reducedDims='IterativeLSI')

# 5. Reproducible peakset (per cluster)
proj <- addGroupCoverages(proj, groupBy='Clusters')
proj <- addReproduciblePeakSet(proj, groupBy='Clusters', pathToMacs2='/path/macs3')
proj <- addPeakMatrix(proj)
```

## SnapATAC2 (Python)

**Goal:** Run a Python-native scATAC pipeline from fragments through clusters, per-cluster peaks, and gene activity.

```python
import snapatac2 as snap

# 1. Load 10X fragments. SnapATAC2 uses snap.pp.import_fragments (NOT snap.read_10x, which doesn't exist).
data = snap.pp.import_fragments(
    fragment_file='outs/fragments.tsv.gz',
    chrom_sizes=snap.genome.hg38,
    file='out.h5ad',                        # backed AnnData; backend handled by file path
    sorted_by_barcode=False)

# 2. Per-cell QC
snap.metrics.tsse(data, gene_anno=snap.genome.hg38)
snap.metrics.frag_size_distr(data)
snap.pp.filter_cells(data, min_counts=1000, min_tsse=4)

# 3. Tile matrix + spectral
snap.pp.add_tile_matrix(data, bin_size=500)
snap.pp.select_features(data, n_features=250000)
snap.tl.spectral(data)
snap.pp.knn(data)                            # tl.leiden reads the kNN graph (obsp['distances']) that only pp.knn builds
snap.tl.umap(data)
snap.tl.leiden(data)

# 4. Per-cluster peak calling (uses MACS3)
# n_jobs=1 in a plain script: default multiprocessing re-imports __main__ (else guard it with
# `if __name__ == '__main__':`). On network/drvfs mounts set HDF5_USE_FILE_LOCKING=FALSE.
snap.tl.macs3(data, groupby='leiden', n_jobs=1)

# 5. Gene activity (gene score) for annotation
gene_mat = snap.pp.make_gene_matrix(data, gene_anno=snap.genome.hg38)
```

## Multiome WNN integration (Signac)

**Goal:** Build a joint RNA + ATAC embedding from a 10X Multiome dataset using Weighted Nearest Neighbors.

Run per-modality embeddings (PCA on RNA, TF-IDF + SVD on ATAC skipping LSI-1), then `FindMultiModalNeighbors` to learn per-cell modality weights and project a joint UMAP.

```r
library(Signac); library(Seurat); library(magrittr)   # magrittr supplies %>%; Signac/Seurat do not attach it
# Assume `obj` has both 'RNA' and 'ATAC' assays from same Multiome experiment
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(obj) <- 'ATAC'
obj <- RunTFIDF(obj) %>% FindTopFeatures(min.cutoff='q0') %>% RunSVD()

# Joint embedding
obj <- FindMultiModalNeighbors(obj, reduction.list=list('pca', 'lsi'),
                               dims.list=list(1:30, 2:30))
obj <- RunUMAP(obj, nn.name='weighted.nn', reduction.name='wnn.umap')
```

Default WNN weights modalities per cell; ATAC signal is much sparser than RNA and can swamp the joint embedding. Inspect per-modality weights; if ATAC noise dominates, adjust `FindMultiModalNeighbors` settings deliberately. Multiome ATAC peaks must be called from the Multiome ATAC fragments, not transferred from a separate scATAC dataset.

## Doublets: AMULET (fragment or BAM input)

Choose AMULET only at sufficient depth (see the doublet section of the core Skill). AMULET v1.1 fails on numpy >= 1.24 (`np.object` removed), so give it its own environment:

```bash
micromamba create -n amulet python=3.10 "numpy<1.24" pandas scipy statsmodels openjdk=17
# github.com/UcarLab/AMULET releases: AMULET-v1.1.zip (AMULET.sh, jar, human_autosomes.txt) and
# RestrictionRepeatLists.zip (restrictionlist_repeats_segdups_rmsk_hg38.bed / _hg19.bed)
# Fragment input (the csv must have columns named 'barcode' and 'is__cell_barcode'; the --*idx flags are ignored):
AMULET.sh fragments.tsv.gz singlecell.csv human_autosomes.txt restrictionlist_repeats_segdups_rmsk_hg38.bed outdir /path/to/AMULET
# BAM input (cell barcode in the CB tag; --forcesorted is required by the jar even for a coordinate-sorted BAM):
AMULET.sh --forcesorted possorted_bam.bam singlecell.csv human_autosomes.txt restrictionlist_repeats_segdups_rmsk_hg38.bed outdir /path/to/AMULET
```

The BAM route reads the csv by column index (defaults `--bcidx 0 --cellidx 0 --iscellidx 9`), so the flags must match the csv layout; the wrong index fails late (`IndexError` in `AMULET.py`) or silently gives 0 overlaps:

| Input | Barcode in the data | Flags |
|-------|--------------------|-------|
| ATAC 1.x / 2.x `singlecell.csv` (`is__cell_barcode` is column 9) | `CB` tag or fragments barcode = `barcode` | defaults |
| ARC `per_barcode_metrics.csv` (columns `barcode`=GEX, `gex_barcode`, `atac_barcode`, `is_cell`) BAM | `CB` tag = GEX barcode (column 0) | `--forcesorted --bcidx 0 --cellidx 0 --iscellidx 3` |
| ARC `atac_fragments.tsv.gz` (barcode = `barcode` column) | fragments barcode = `barcode` | fragment mode reads columns by name: write a two-column csv `barcode,is__cell_barcode` from `barcode` and `is_cell` |

Multiplet calls (q < 0.01) are in `outdir/MultipletBarcodes_01.txt`, per-barcode probabilities in `MultipletProbabilities.txt`. `--forcesorted` applies to BAM input only.
