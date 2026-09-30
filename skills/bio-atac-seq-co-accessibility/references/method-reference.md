# Co-Accessibility Method Reference

## Version Compatibility

Reference examples tested with: cicero 1.3.x (GitHub `cole-trapnell-lab/cicero-release`, branch `monocle3`; the Bioconductor `cicero` is the monocle2 package and lacks the API used here), monocle3 1.3+, ArchR 1.0.3, Signac 1.17, GenomicRanges 1.54+, GenomicInteractions 1.36+. SCENIC+ (1.0a2 alpha, Python <=3.11.8) and pycisTopic were not executed.

Verify before use:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws unexpected errors, introspect the installed package and adapt rather than retrying.

## What Co-accessibility Captures vs What It Doesn't

| Captures | Misses |
|----------|--------|
| Peak pairs that vary together across cell states | 3D physical contacts that don't vary in accessibility |
| Cis-regulatory grammar within a cell type | Trans-chromosomal interactions |
| Active enhancer-promoter pairs | Constitutive structural contacts |
| Lineage-specific regulation | Developmental contacts that opened before scATAC sample |
| Distance-decay biology of enhancer-promoter | Hub enhancers that contact many distal targets |

For physical contact, use Hi-C, Micro-C, or PCHi-C. Co-accessibility is the chromatin-only proxy.

## Algorithmic Taxonomy

| Tool | Method | Input | Output | Strength | Fails when |
|------|--------|-------|--------|----------|------------|
| Cicero (Pliner 2018) | Graphical lasso on aggregated cell metacells | scATAC peak-cell matrix + cell trajectory | Peak-pair connection score (-1 to 1) | Original, well-validated; integrates with Monocle3 | Slow on large datasets (">50K cells" is an unverified rule of thumb); sensitive to alpha tuning |
| ArchR getCoAccessibility | Cicero-based; uses ArchR's metacell aggregation | ArchR project | Same as Cicero | Built-in to ArchR pipeline; faster on large datasets | Tied to ArchR; same biology as Cicero |
| SCENIC+ (Bravo 2023) | Multi-step: co-accessibility + motif scoring + RNA correlation | Multiome (ATAC + RNA) or paired | TF-driven enhancer-gene networks | Most comprehensive; multi-modal | Multiome data required; computationally heavy |
| LinkPeaks (Signac) | Pearson correlation of accessibility with paired gene expression | Multiome | Peak-gene linkage score | Direct enhancer-gene from RNA correlation | Multiome-only; not pure ATAC |
| GeneHancer / FANTOM5 / EpiMap | Bulk-derived enhancer-gene reference | None (database lookup) | Pre-computed enhancer-gene pairs | Comprehensive; published references | Cell-type-agnostic; may not match the biology of interest |

Methodology evolves; verify against Pliner 2018 (Cicero), Bravo 2023 (SCENIC+), Nasser 2021 (ABC model alternative for enhancer-gene), and current Hi-C concordance benchmarks.

## How Cicero Works (Conceptually)

Cell-to-cell variability is too sparse for direct correlation. Cicero solves this via metacells:

1. Reduce dimensionality (UMAP from input).
2. Build k-NN graph of cells.
3. Aggregate k cells into metacells (default k = 50).
4. Compute correlation in accessibility across metacells, restricted to peak pairs within `window` (`run_cicero` argument, default 500 kb cis).
5. Apply graphical lasso with regularization `alpha` to sparsify the correlation matrix.
6. Output: per-pair connection score in -1..1 (observed -0.88..0.93 on PBMC chr1; about a third of scored pairs were negative); positive = co-variation, negative = anti-co-variation. Pairs Cicero cannot score are NA. Every pair appears in both orientations (Peak1/Peak2 swapped), so de-duplicate before counting.

Connection thresholds typically 0.05-0.5; > 0.25 is the Cicero-vignette working default.

## Per-Tool Failure Modes

### Cicero -- alpha tuning shifts results

**Trigger:** Default alpha (sometimes computed automatically from data); custom alpha < 0.5 or > 5.

**Mechanism:** Alpha controls graphical lasso regularization. Too low: dense graph with many spurious connections; too high: sparse with biology missing.

**Symptom:** Connection count varies 10-100x across alpha sweeps.

**Fix:** Use Cicero's `estimate_distance_parameter()` to get data-driven alpha; verify the connection count against a permutation control (below). The fraction of peaks with at least one strong connection is dataset dependent (60% on the PBMC chr1 slice at > 0.25), so it is not a validity test.

### Cicero -- metacell aggregation hides cell-type-specific connections

**Trigger:** Running Cicero on heterogeneous dataset spanning multiple cell types.

**Mechanism:** Metacells aggregate across cell types; connections that exist only in one cell type get diluted.

**Fix:** Run Cicero per-cluster separately; combine results with cluster annotations. Cell-type-specific connections often differ.

### Cicero -- distance assumption

**Trigger:** Default `window=500000` (500 kb cis only).

**Mechanism:** Distal connections beyond 500 kb cis are excluded; trans-chromosomal entirely missed.

**Fix:** For specific use cases (e.g., gene desertless TADs), increase `window` to 1 Mb or more (`run_cicero` keeps `distance_constraint`, used only for the alpha search, fixed at 250 kb). Trans connections require Hi-C, not co-accessibility.

**Executed check:** chr1 1-10 Mb, 547 peaks x 3,277 cells: `window=5e5` gave 16,158 scored pairs (max start-to-start span 498 kb); `window=1e6` gave 28,363 (7,804 beyond 500 kb, max 989 kb).

### SCENIC+ -- RNA scaling

**Trigger:** RNA-side dropouts in Multiome data.

**Mechanism:** SCENIC+ requires reasonable RNA quantification per cell. Sparse Multiome RNA with many zero genes causes correlation degradation.

**Fix:** Filter cells with insufficient RNA; aggregate cells if necessary. Multiome RNA should look comparable to standalone scRNA-seq.

### LinkPeaks (Signac) -- Distance default

**Trigger:** Default `LinkPeaks(..., distance=5e+05)`.

**Mechanism:** Same as Cicero; 500 kb cis only by default.

**Fix:** Same; widen if needed but trans not supported.

## Cicero Standard Workflow

**Goal:** Infer cis-regulatory peak-peak connections from a scATAC peak-cell matrix.

**Approach:** Build a Monocle3 CellDataSet, reduce dimensions via LSI + UMAP, aggregate cells into metacells, then run Cicero's graphical-lasso correlation across the cis window and threshold on connection score.

```r
library(cicero); library(monocle3); library(GenomicRanges)

# Input: peak-cell binary matrix from Signac/ArchR (rows = peaks, cols = cells)
# Convert peaks to "chrN_start_end" format
peak_names <- paste0(seqnames(peaks), '_', start(peaks), '_', end(peaks))
input_cds <- new_cell_data_set(peak_matrix, cell_metadata=metadata,
                               gene_metadata=peak_metadata)

# Reduce dimensionality (UMAP from input)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method='LSI')
input_cds <- reduce_dimension(input_cds, reduction_method='UMAP',
                              preprocess_method='LSI')

# Build metacell-aggregated CDS
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates=umap_coords, k=50)

# Model only the chromosomes present in the peak set (genome-build independent).
# Do not list every hg38 seqname: alt/random contigs and empty chromosomes only add runtime.
peak_gr <- GRanges(sub('_([0-9]+)_([0-9]+)$', ':\\1-\\2', peak_names))
chr_end <- tapply(end(peak_gr), as.character(seqnames(peak_gr)), max)
genome_df <- data.frame(chr=names(chr_end), length=as.numeric(chr_end))
conns <- run_cicero(cicero_cds, genomic_coords=genome_df,
                    window=500000, sample_num=100)

# One row per unordered pair; drop unscored (NA) pairs
conns$Peak1 <- as.character(conns$Peak1); conns$Peak2 <- as.character(conns$Peak2)
conns <- conns[!is.na(conns$coaccess), ]
lo <- pmin(conns$Peak1, conns$Peak2); hi <- pmax(conns$Peak1, conns$Peak2)
conns$Peak1 <- lo; conns$Peak2 <- hi
conns <- conns[!duplicated(paste(lo, hi)), ]

# Filter to high-confidence connections.
# Threshold 0.25 is a Cicero-documentation working default; the optimal cutoff
# is dataset-dependent and is best calibrated against orthogonal Hi-C / HiChIP.
strong <- conns[conns$coaccess > 0.25, ]
cat(sprintf('Total conns: %d; strong (>0.25): %d\n', nrow(conns), nrow(strong)))
```

## ArchR getCoAccessibility

```r
library(ArchR)
proj <- loadArchRProject('ArchR_out')
proj <- addCoAccessibility(proj, reducedDims='IterativeLSI',
                          k=100, knnIteration=500,
                          maxDist=250000)               # 250 kb cis (wider than the 100 kb default)
co_acc <- getCoAccessibility(proj, corCutOff=0.5,       # Default 0.5 in ArchR; lower for more (calibrate against Hi-C/HiChIP)
                             returnLoops=TRUE)           # TRUE (default) -> SimpleList, loops GRanges in element 1; FALSE -> DataFrame of peak-pair correlations
loops <- co_acc[[1]]                                    # GRanges of loops; ArchR's own help calls the return a GRanges, it is a SimpleList
```

With `returnLoops=TRUE` (the default) ArchR returns a SimpleList whose first element is the loops GRanges (compatible with `GenomicInteractions` for direct overlap with Hi-C loops); `returnLoops=FALSE` instead returns a DataFrame of peak-pair correlations.

## Visualizing Connections

```r
# Arc plot at ONE locus with a stringent cutoff: plotTracks(track) alone errors ("Unable to automatically
# determine plotting ranges"), and all strong arcs genome-wide draw an unreadable hairball.
library(Gviz); library(GenomicInteractions)
# Cicero Peak1/Peak2 are chr_start_end strings; convert to chr:start-end for GRanges()
to_gr <- function(x) GRanges(sub('_(\\d+)_(\\d+)$', ':\\1-\\2', x))
sel <- strong[strong$coaccess > 0.5, ]
gi <- GenomicInteractions(anchor1=to_gr(sel$Peak1), anchor2=to_gr(sel$Peak2),
                          counts=as.integer(sel$coaccess * 100))
track <- InteractionTrack(gi, name='co-accessibility > 0.5', chromosome='chr1', start=1e6, end=1.4e6)
plotTracks(list(GenomeAxisTrack(), track), chromosome='chr1', from=1e6, to=1.4e6)
```

For genome-browser visualization with ArchR: `plotPeak2GeneHeatmap()` shows the peak-gene linkage matrix; `plotBrowserTrack()` overlays connections on tracks.

## SCENIC+ TF-Driven Networks

SCENIC+ (1.0a2) runs as a Snakemake pipeline (CLI), not a single monolithic Python call; this route was reviewed statically only (not installable in the audit environment). Prepare the inputs first (a pycisTopic cisTopic object, motif-enrichment results, and paired RNA AnnData), then scaffold and run the workflow:

```bash
# Scaffold the pipeline, then edit its config.yaml to point at the cisTopic object,
# motif-enrichment results, and GEX AnnData
scenicplus init_snakemake --out_dir scplus_pipeline/
snakemake --cores 16 --snakefile scplus_pipeline/Snakemake/workflow/Snakefile
# eRegulons (TF + target genes + linked enhancers) are written to the output MuData (scplusmdata.h5mu)
```

SCENIC+ is significantly more complex than Cicero; budget 1-2 days for setup. The benefit is that outputs are TF -> enhancer -> gene triples, not just peak-peak co-accessibility.

## Cicero Alpha Mathematics

**Trigger:** Tuning Cicero's regularization parameter for the graphical lasso step.

**Mechanism:** `estimate_distance_parameter()` searches for the smallest distance-penalty scaling (Cicero's `distance_parameter`, called "alpha" here) such that, across random genomic windows, no more than ~5% of peak pairs beyond `distance_constraint` retain non-zero graphical-lasso entries and fewer than 80% of all entries are non-zero. This penalizes long-range co-accessibility so the graph sparsifies at biologically appropriate distance scales -- it is not a correlation-vs-distance regression slope.

**Implementation:** Cicero calls `estimate_distance_parameter(cicero_cds, window=window, maxit=100, sample_num=100, genomic_coords=genome_df)` over `sample_num` random windows and returns one `distance_parameter` per window; take the mean and pass it to `generate_cicero_models(cicero_cds, distance_parameter=mean(...))`. Supply `genomic_coords` explicitly -- its default is `cicero::human.hg19.genome`, wrong for an hg38 analysis; derive it from the peak set as in the workflow above.

**When manual tuning helps:** Very dense peaksets (>200k peaks) may need a higher `distance_parameter` to control false positives; very sparse (<10k peaks) may need a lower one to recover signal. Verify by running on a permutation / cell-label-shuffle negative control -- the expected outcome is ~0 strong connections (technical replicates should instead reproduce connections).

## ABC Model Cross-Reference

For enhancer-to-gene linking with paired Hi-C/Micro-C, the canonical method is the ABC model (Fulco 2019, Nasser 2021), not Cicero. ABC computes ABC = (Activity_E * Contact_E,G) / sum_e(Activity_e * Contact_e,G); standardizes on combined ATAC + H3K27ac activity and Hi-C contact frequencies. ENCODE-rE2G (Gschwind et al 2023, bioRxiv) is the modern logistic-regression enhancer-gene link predictor.

See atac-seq/enhancer-gene-linking for full ABC and ENCODE-rE2G coverage. Cicero is the ATAC-only fallback when no Hi-C is available.

## HiChIP H3K27ac as Orthogonal Anchor

| Decision | Action |
|----------|--------|
| Have Hi-C / Micro-C | Use ABC (atac-seq/enhancer-gene-linking) primary; Cicero as ATAC-only sanity check |
| Have HiChIP H3K27ac | FitHiChIP loops (FDR < 0.05, count >= 5) primary; ABC + HiChIP intersection is high-confidence |
| Have ATAC + H3K27ac, no 3D | ABC with average HiC fallback (Fulco 2019); document degraded performance |
| Have only ATAC | Cicero (this skill); Hi-C concordance is partial |

Cicero is appropriate when no 3D data exists; do not use Cicero in lieu of ABC when Hi-C/Micro-C are available.

## Hi-C / Micro-C Concordance

Bands below are unsourced heuristics, not benchmarks; calibrate on your own Hi-C.

| Hi-C concordance | Action |
|-----------------|--------|
| > 50% of strong Cicero connections overlap Hi-C loops | High-confidence; Cicero captures real 3D structure |
| 30-50% | Common; some 3D contacts don't vary in accessibility |
| < 20% | Co-accessibility may not reflect contacts; lineage-specific contacts may be missing |

**Goal:** Quantify what fraction of strong Cicero connections are supported by Hi-C loop calls.

**Approach:** Import HiCCUPS loops as GenomicInteractions, build a parallel object from Cicero connections, then count anchor-anchor overlaps and report the percentage.

```r
# Compare Cicero against published Hi-C loops
library(GenomicInteractions)
hic_loops <- makeGenomicInteractionsFromFile('hiccups_loops.bedpe', type='bedpe',
                                             experiment_name='hiccups', description='HiCCUPS loops')
ci <- GenomicInteractions(anchor1=GRanges(sub('_(\\d+)_(\\d+)$', ':\\1-\\2', strong$Peak1)),
                          anchor2=GRanges(sub('_(\\d+)_(\\d+)$', ':\\1-\\2', strong$Peak2)))
overlap <- countOverlaps(ci, hic_loops) > 0   # anchor-anchor 'any' overlap; 'equal' is too stringent at loop bin resolution
cat(sprintf('Cicero connections overlapping HiCCUPS loops: %.1f%%\n',
            100 * mean(overlap)))
```

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Cicero many weak connections; ArchR few strong | Different alpha or aggregation | Standardize parameters |
| LinkPeaks (Multiome) finds connections Cicero misses | LinkPeaks uses RNA expression as the anchor; Cicero is ATAC-only | Both valid; report intersection as high-confidence |
| Co-accessibility doesn't match Hi-C in heterochromatin | Heterochromatic contacts are constitutive; co-accessibility needs variation | Expected; co-accessibility complements Hi-C |
| SCENIC+ network has ENCODE-validated TFs but missing some | Motif database limited or RNA imputation missed | Expand motif database; integrate paired ChIP-seq if available |

**Operational rule:** Co-accessibility is a hypothesis generator. Validate with Hi-C, ChIP-seq, or experimental enhancer-promoter interaction (CRISPRi-FlowFISH).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Cicero `make_cicero_cds` slow / crashes | k too high or cell count too large | Reduce k or subsample cells |
| All connections near zero | alpha set too high | Use `estimate_distance_parameter()` |
| Connection score > 1 reported | Bug in older Cicero versions | Update; check `as.numeric(coaccess)` for outliers |
| ArchR getCoAccessibility "TileMatrix" error | Need PeakMatrix not TileMatrix | `addPeakMatrix()` first |
| SCENIC+ install fails | Many heavy dependencies; pybedtools wheel build failed here | Use a Python <=3.11 env with conda-forge `pybedtools`; untested here; the upstream README documents no Docker image |
| Connection count varies wildly per run | Stochastic metacell aggregation | Set seed; or aggregate at higher k for stability |
| LinkPeaks all NaN | RNA expression has too many zeros | Re-filter cells with sufficient RNA |
| Peak names not matching | format mismatch (chr_start_end vs chr:start-end); the script accepts only `chr_start_end` | Standardize naming convention |

## References

- Pliner HA et al 2018 Mol Cell 71:858 (Cicero)
- Granja JM et al 2021 Nat Genet 53:403 (ArchR getCoAccessibility)
- Bravo Gonzalez-Blas C et al 2023 Nat Methods 20:1355 (SCENIC+)
- Stuart T et al 2021 Nat Methods 18:1333 (Signac LinkPeaks)
- Nasser J et al 2021 Nature 593:238 (ABC model; alternative enhancer-gene)
- Fulco CP et al 2019 Nat Genet 51:1664 (CRISPRi-FlowFISH; gold-standard validation)
- Mumbach MR et al 2017 Nat Genet 49:1602 (HiChIP H3K27ac for enhancer-promoter)
- Boix CA et al 2021 Nature 590:300 (EpiMap; bulk enhancer-gene reference)
