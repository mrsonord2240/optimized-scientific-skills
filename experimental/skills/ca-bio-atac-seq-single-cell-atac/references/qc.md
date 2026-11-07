# QC and doublets



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



- Signac (standard R pipeline): run [`scripts/signac_workflow.R`](scripts/signac_workflow.R) (`Rscript scripts/signac_workflow.R <h5> <fragments.tsv.gz> <singlecell.csv | per_barcode_metrics.csv>`). It loads 10X output, computes per-cell QC, filters, runs TF-IDF + SVD, UMAP and Leiden clustering with LSI component 1 skipped (`dims=2:30`, because component 1 tracks sequencing depth), and builds a gene-activity assay. Optional args 4 and 5 set the TSS and fragment minimums (defaults 4 and 1000). It assumes hg38 with `EnsDb.Hsapiens.v86` and converts that annotation to UCSC style (`seqlevelsStyle(ann) <- 'UCSC'`, `genome(ann) <- 'hg38'`) because `CreateChromatinAssay` rejects Ensembl-style annotation on an hg38/UCSC object; apply the same conversion to any assay you build from `GetGRangesFromEnsDb`. Edit the genome and annotation for other assemblies. Metadata is auto-detected: ATAC 1.x/2.x `singlecell.csv` (`passed_filters`, `peak_region_fragments`, `blacklist_region_fragments`) or cellranger-arc `per_barcode_metrics.csv` (has `atac_fragments`). For ARC it reads the `Peaks` matrix from the combined h5 (matrix, fragments and first csv column all use the GEX-style `barcode`), takes `atac_fragments` and `atac_peak_region_fragments` as the fragment and peak counts, takes the blacklist ratio from the peak matrix (`FractionCountsInRegion` on `blacklist_hg38_unified`, since ARC has no blacklist column; few ARC peaks overlap the blacklist, so this ratio stays near 0 (max 0.004 on PBMC 3k, versus up to 0.18 from the ATAC `singlecell.csv` column) and the 0.05 cut removes no cells) and mitochondrial fraction as `atac_mitochondrial_reads / atac_raw_reads`. It stops with a clear message if required columns are missing.
- ArchR, SnapATAC2, per-cluster pseudobulk peak calling (MACS3), and Multiome WNN: see [`references/ecosystem-workflows.md`](references/ecosystem-workflows.md).

Per-cluster pseudobulk peaks need enough reads: peak calls from small clusters are unreliable (rule of thumb, not a tested cutoff: aggregate clusters under ~200 cells into a "rare" group or drop them), and use the union of larger-cluster peaks for rare-cell analysis. Build the across-cluster consensus with `bio-atac-seq-consensus-peakset`. Run chromVAR / AddMotifs only after the peakset is final.
