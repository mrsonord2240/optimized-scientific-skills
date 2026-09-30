---
name: bio-atac-seq-single-cell-atac
category: Data Analysis
description: Process and analyze single-cell ATAC-seq data with Signac, ArchR, SnapATAC2, or Cell Ranger ATAC. Use when handling 10X scATAC or 10X Multiome (paired RNA+ATAC) data, performing per-cell QC, choosing between ArchR/Signac/SnapATAC2 ecosystems, building per-cluster consensus peaksets, integrating with paired scRNA-seq, doublet detection (AMULET vs ArchR vs scDblFinder), or running pseudobulk differential accessibility per cluster.
tool_type: mixed
primary_tool: Signac
license: MIT
author: GPTomics
---

# Single-Cell ATAC-seq

Build a per-cell fragment matrix, compute per-cell QC, reduce dimensionality (TF-IDF + LSI, spectral, or autoencoder), cluster, call cluster-level pseudobulk peaks, annotate cell types with gene-activity scores, and integrate with paired scRNA-seq for Multiome.

- R: `Signac::CreateChromatinAssay()` then the Seurat workflow (TF-IDF + SVD + UMAP + Leiden)
- R: `ArchR::createArrowFiles()` then an ArchR project (TileMatrix + LSI + UMAP)
- Python: `snapatac2.pp.import_fragments()` then SnapATAC2 (spectral clustering)
- CLI preprocessing: `cellranger-atac count` (10X) or `chromap` (alignment-only fragment files)

## Version compatibility

Reference examples were written against: Cell Ranger ATAC 2.1+, Signac 1.13+, Seurat 5.0+, ArchR 1.0.2+, SnapATAC2 2.8+, AMULET 1.1+, scDblFinder 1.16+, scater 1.30+, scvi-tools 1.1+, GenomicRanges 1.54+, JASPAR2024 0.99+, BSgenome.Hsapiens.UCSC.hg38 1.4+, EnsDb.Hsapiens.v86 2.99+, MACS3 3.0+. SnapATAC2 2.8+ uses `pp.import_fragments`; 2.5-2.7 used `pp.import_data` (2.8.0 ships both; 2.9.0 and 2.10.0 have only `import_fragments`). Executed 2026-09-30 with Signac 1.17.1, Seurat 5.5.1, ArchR 1.0.3, SnapATAC2 2.10.0, scvi-tools 1.5.1, scDblFinder 1.23.4, MACS3 3.0.4, AMULET 1.1; Cell Ranger ATAC/ARC were not run; their outputs were read from 10x public datasets (ATAC 1.0.1, ATAC 2.1.0, ARC 2.0.0).

Before use, check installed versions and signatures: `pip show <package>` and `help(module.function)` (Python), `packageVersion('<pkg>')` and `?function_name` (R), `<tool> --version` and `--help` (CLI). If code errors unexpectedly, introspect the installed package and adapt rather than retrying.

## Ecosystem choice (the most important decision)

| Ecosystem | Language | Strength | Fails when | Best for |
|-----------|---------|---------|------------|----------|
| Signac (Stuart 2021) | R, Seurat-based | Tightest scRNA-seq integration; mature | Memory hungry on >100K cells; slower than ArchR | Multiome RNA+ATAC; small-to-medium datasets; Seurat users |
| ArchR (Granja 2021) | R, Arrow/HDF5 | Memory-efficient; fast on 100K-1M cells; built-in trajectory + doublet | Less RNA integration; ArchR-specific format | Large scATAC cohorts; trajectories; ATAC-only |
| SnapATAC2 (Zhang 2024) | Python, AnnData | Memory-efficient; performant spectral clustering | Newer; smaller ecosystem | Python-first labs; very large datasets (>1M cells) |
| Cell Ranger ATAC | CLI (10X) | Official 10X preprocessing | Closed, fixed pipeline | Preprocessing only |
| scATAC-pro | CLI pipeline | Alternative preprocessing | Less maintained | Legacy; not recommended for new projects |

Methodology evolves; verify against Granja 2021, Stuart 2021, Zhang 2024, and Heumos 2023 before locking a pipeline.

Goal-to-pipeline routing:

| Goal | Pipeline |
|------|----------|
| Standard scATAC (R user) | Signac: CreateChromatinAssay -> RunTFIDF -> RunSVD -> RunUMAP (dims 2:30) -> FindClusters |
| Standard scATAC (Python user) | SnapATAC2: pp.import_fragments -> add_tile_matrix -> spectral -> pp.knn -> UMAP -> leiden |
| >100K cells | ArchR (Arrow files) |
| 10X Multiome | Signac + Seurat: per-modality embedding, then WNN |
| Trajectory / pseudotime | ArchR getTrajectory; or Signac + Cicero |
| Differential accessibility per cluster | Pseudobulk per cluster -> consensus peakset -> DESeq2 |
| Cell-type annotation | Gene-activity scores (ArchR or Signac), then `single-cell/markers-annotation` |
| Multimodal trajectories | MOFA+, ArchR + scRNA integration, or SCENIC+ |
| Plant / non-model | Signac with custom EnsDb / TxDb; ArchR with custom annotations |

## 10X Multiome caveat (paired RNA + ATAC)

Multiome chemistry profiles RNA and ATAC from the same cell, joined by a shared barcode; use Signac for ATAC and Seurat for RNA in one object, integrated by WNN (Hao 2021). Single-modality 10X scATAC (chemistry v1, v2) produces no paired RNA; verify the chemistry on the cellranger summary before assuming Multiome. See "cellranger-atac vs cellranger-arc" in [`references/specialized-topics.md`](references/specialized-topics.md).

## Per-cell QC thresholds

| Metric | Definition | Pass | Caution | Reject | Source |
|--------|-----------|------|---------|--------|--------|
| Fragment count per cell | n_fragments after dedup | 3000-50000 | 1000-3000 | < 1000 or > 80000 | 10X recommendation; high = doublet |
| TSS enrichment per cell | Signal at TSS / flanks | >= 4 | 2-4 | < 2 | ArchR / Signac default; lower than bulk |
| Nucleosome signal | Mono / NFR fragment ratio | <= 4 | 4-10 | > 10 | Signac default; high = poor library |
| % reads in peaks (per cell) | FRiP per cell at consensus | >= 0.15 | 0.10-0.15 | < 0.05 | ArchR / Signac defaults |
| Mitochondrial fraction | chrM / total per cell | < 0.05 | 0.05-0.15 | > 0.20 | Lower than bulk; per-cell more sensitive |
| Doublet score | AMULET / ArchR doublet | < 0.5 | 0.5-0.7 | > 0.7 | Tool-dependent threshold |
| Blacklist ratio | Blacklist-region fragments / peak-region fragments (Signac convention) | < 0.05 | 0.05-0.10 | > 0.10 | Standard |

Per-cell thresholds are looser than bulk because individual cells carry orders of magnitude less signal; the population aggregate is what matters.

**Filter rule** (single definition; [`scripts/signac_workflow.R`](scripts/signac_workflow.R) implements it, and other pages point here): keep cells with 1000-80000 fragments (the 1000-3000 caution band is kept, flagged for inspection), TSS enrichment >= 4, nucleosome signal <= 4, FRiP >= 15%, blacklist ratio < 0.05 and, when the metadata has it, mitochondrial fraction < 0.05. Doublet removal is a separate step. Cell Ranger ATAC cell calling is lenient (fragment-based heuristic, no UMIs): on the 10x PBMC 5k dataset 7.2% of called cells have < 1000 fragments. Shallow or subset data can leave zero cells at these values; lower TSS and fragment minimums deliberately and report the change.

## Doublet detection

| Tool | Method | Strength | Fails when |
|------|--------|---------|------------|
| AMULET (Thibodeau 2021) | Collision-based: too many fragments at one position (impossible from one diploid cell) | ATAC-specific; orthogonal to clustering | Shallow libraries: recall is about 0.85-0.90 near 25K valid read pairs/cell and ArchR outperforms it at lower depth |
| ArchR addDoubletScores | Synthetic doublets + LSI projection | Built into ArchR; auto-thresholds | Tied to ArchR's LSI; not portable |
| scDblFinder (Germain 2021) | Synthetic doublets + classifier | Works on Signac and SCE objects; well-benchmarked | RNA-developed; ATAC use needs careful settings |

Choose the primary tool by depth: measure median valid read pairs per cell first (in the 10x PBMC 5k data the median called cell has 10.7K passed fragments and 70% are below 15K). Use AMULET at about 25K or more; below that use ArchR or scDblFinder. Cells flagged by two tools are a high-precision doublet set at the cost of sensitivity, not the default. AMULET command and environment: [`references/ecosystem-workflows.md`](references/ecosystem-workflows.md).

## Workflows

- Signac (standard R pipeline): run [`scripts/signac_workflow.R`](scripts/signac_workflow.R) (`Rscript scripts/signac_workflow.R <h5> <fragments.tsv.gz> <singlecell.csv | per_barcode_metrics.csv>`). It loads 10X output, computes per-cell QC, filters, runs TF-IDF + SVD, UMAP and Leiden clustering with LSI component 1 skipped (`dims=2:30`, because component 1 tracks sequencing depth), and builds a gene-activity assay. Optional args 4 and 5 set the TSS and fragment minimums (defaults 4 and 1000). It assumes hg38 with `EnsDb.Hsapiens.v86` and converts that annotation to UCSC style (`seqlevelsStyle(ann) <- 'UCSC'`, `genome(ann) <- 'hg38'`) because `CreateChromatinAssay` rejects Ensembl-style annotation on an hg38/UCSC object; apply the same conversion to any assay you build from `GetGRangesFromEnsDb`. Edit the genome and annotation for other assemblies. Metadata is auto-detected: ATAC 1.x/2.x `singlecell.csv` (`passed_filters`, `peak_region_fragments`, `blacklist_region_fragments`) or cellranger-arc `per_barcode_metrics.csv` (has `atac_fragments`). For ARC it reads the `Peaks` matrix from the combined h5 (matrix, fragments and first csv column all use the GEX-style `barcode`), takes `atac_fragments` and `atac_peak_region_fragments` as the fragment and peak counts, takes the blacklist ratio from the peak matrix (`FractionCountsInRegion` on `blacklist_hg38_unified`, since ARC has no blacklist column; few ARC peaks overlap the blacklist, so this ratio stays near 0 (max 0.004 on PBMC 3k, versus up to 0.18 from the ATAC `singlecell.csv` column) and the 0.05 cut removes no cells) and mitochondrial fraction as `atac_mitochondrial_reads / atac_raw_reads`. It stops with a clear message if required columns are missing.
- ArchR, SnapATAC2, per-cluster pseudobulk peak calling (MACS3), and Multiome WNN: see [`references/ecosystem-workflows.md`](references/ecosystem-workflows.md).

Per-cluster pseudobulk peaks need enough reads: peak calls from small clusters are unreliable (rule of thumb, not a tested cutoff: aggregate clusters under ~200 cells into a "rare" group or drop them), and use the union of larger-cluster peaks for rare-cell analysis. Build the across-cluster consensus with `atac-seq/consensus-peakset`. Run chromVAR / AddMotifs only after the peakset is final.

## Per-tool failure modes

- Signac LSI component 1 is depth (|r| with log fragment count near 1, sign arbitrary): using `dims=1:30` adds depth-associated structure to the UMAP (on the 10x Multiome PBMC data max |cor(UMAP, depth)| was 0.10 vs 0.05 with `2:30`). Use `dims=2:30` for `RunUMAP` and `FindNeighbors` (ArchR and SnapATAC2 handle this automatically).
- ArchR TileMatrix vs PeakMatrix: TileMatrix (fixed 500 bp bins) is for embedding and clustering; use PeakMatrix (after `addReproduciblePeakSet`) for differential, motif, and gene-activity analysis.
- SnapATAC2 expects raw integer (Int32) fragment counts; convert float or negative-valued matrices before loading.
- Multiome WNN weights RNA and ATAC per cell; sparse ATAC can inject noise. Inspect per-modality weights and adjust deliberately.

Further failure modes, an error lookup table, and reconciliation between ecosystems are in [`references/specialized-topics.md`](references/specialized-topics.md).

## Routed material

- [`references/ecosystem-workflows.md`](references/ecosystem-workflows.md): ArchR, SnapATAC2, per-cluster MACS3 peak calling, Multiome WNN code.
- [`references/specialized-topics.md`](references/specialized-topics.md): cellranger-atac vs cellranger-arc, cell-cycle correction, sex-chromosome QC, scArches/PEAKVI mapping, chromBPNet per cluster, ecosystem reconciliation, common errors.
- [`references/usage-guide.md`](references/usage-guide.md): installation, expected inputs, example requests, working tips.

## Related skills

`atac-seq/atac-qc`, `atac-seq/atac-peak-calling`, `atac-seq/consensus-peakset`, `atac-seq/differential-accessibility`, `atac-seq/motif-deviation`, `atac-seq/footprinting`, `atac-seq/co-accessibility`, `atac-seq/deep-learning-atac`, `atac-seq/enhancer-gene-linking`, `atac-seq/allele-specific-accessibility`, `single-cell/preprocessing`, `single-cell/clustering`, `single-cell/cell-annotation`, `single-cell/multimodal-integration`, `single-cell/scatac-analysis`, `single-cell/batch-integration`.

## References

- Stuart T et al 2021 Nat Methods 18:1333 (Signac)
- Granja JM et al 2021 Nat Genet 53:403 (ArchR)
- Zhang K et al 2024 Nat Methods 21:217 (SnapATAC2)
- Hao Y et al 2021 Cell 184:3573 (Seurat WNN)
- Thibodeau A et al 2021 Genome Biol 22:252 (AMULET)
- Germain PL et al 2021 F1000Res 10:979 (scDblFinder)
- Cusanovich DA et al 2015 Science 348:910 (sciATAC; LSI for sc data)
- Chen H et al 2019 Genome Biol 20:241 (scATAC analysis benchmark)
- Heumos L et al 2023 Nat Rev Genet 24:550 (single-cell best practices)
- 10X Genomics Cell Ranger ATAC documentation

## License and provenance

This derived Skill retains GPTomics/bioSkills material by Domen Jemec, from commit `d91ed3d563019e649dc854c56ccd62551359488a` (`atac-seq/single-cell-atac`). See [`LICENSE`](LICENSE) for the preserved MIT notice.
