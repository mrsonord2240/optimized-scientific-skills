# Cicero

The source script directly converts `chr_start_end` to GRanges coordinates. GRanges uses one-based closed coordinates, whereas BED commonly uses zero-based half-open coordinates; verify and convert the peak metadata to the convention expected by the script before gene mapping. Its `resize(..., width=4000, fix='center')` is centered on the supplied BED interval, so provide TSS intervals rather than arbitrary gene spans.



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



**Trigger:** Tuning Cicero's regularization parameter for the graphical lasso step.

**Mechanism:** `estimate_distance_parameter()` searches for the smallest distance-penalty scaling (Cicero's `distance_parameter`, called "alpha" here) such that, across random genomic windows, no more than ~5% of peak pairs beyond `distance_constraint` retain non-zero graphical-lasso entries and fewer than 80% of all entries are non-zero. This penalizes long-range co-accessibility so the graph sparsifies at biologically appropriate distance scales -- it is not a correlation-vs-distance regression slope.

**Implementation:** Cicero calls `estimate_distance_parameter(cicero_cds, window=window, maxit=100, sample_num=100, genomic_coords=genome_df)` over `sample_num` random windows and returns one `distance_parameter` per window; take the mean and pass it to `generate_cicero_models(cicero_cds, distance_parameter=mean(...))`. Supply `genomic_coords` explicitly -- its default is `cicero::human.hg19.genome`, wrong for an hg38 analysis; derive it from the peak set as in the workflow above.

**When manual tuning helps:** Very dense peaksets (>200k peaks) may need a higher `distance_parameter` to control false positives; very sparse (<10k peaks) may need a lower one to recover signal. Verify by running on a permutation / cell-label-shuffle negative control -- the expected outcome is ~0 strong connections (technical replicates should instead reproduce connections).
