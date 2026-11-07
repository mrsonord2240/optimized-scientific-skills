# Normalization

The retained source spike-in snippet uses inverse spike counts as DESeq2 size factors and calls arithmetic-mean scaling a geometric mean. The factor direction and scaling need method-specific verification before use; it is not an executable validated workflow in this draft.



DiffBind 3.x conflates two orthogonal choices: the normalization method (`normalize=`) and the library-size definition (`library=`). DiffBind's own defaults are `normalize=DBA_NORM_DEFAULT`, `library=DBA_LIBSIZE_DEFAULT` (native RLE/TMM on reads-in-peaks, verified on 3.16); the script's `lib` mode sets `library=DBA_LIBSIZE_FULL` (full mapped-read total) instead.

| Choice | DiffBind argument | What it does | When to use |
|--------|-------------------|--------------|-------------|
| Normalize by library size only | `normalize=DBA_NORM_LIB` (script default) | Scale counts by the chosen library size | Standard; pairs with full or RiP library |
| Reads-in-peaks library size | `library=DBA_LIBSIZE_PEAKREADS` | Library size = reads in consensus peaks (RiP) | When background varies independently of biology (protects against background drift) |
| Full mapped library size | `library=DBA_LIBSIZE_FULL` (script default) | Library size = total mapped reads | When global accessibility shifts must remain visible (e.g., chromatin compaction) |
| Native per-tool default | `normalize=DBA_NORM_NATIVE` | DESeq2 RLE or edgeR TMM, depending on backend | Use DESeq2/edgeR conventions directly |
| TMM (edgeR) | `normalize=DBA_NORM_TMM` | Trimmed mean of M-values | Robust to a few highly-DA peaks dominating |
| RLE (DESeq2) | `normalize=DBA_NORM_RLE` | DESeq2 geometric-mean size factors | DESeq2-conventional analysis |
| Spike-in / external | not built-in; pre-scale counts | Exogenous reference (e.g., spike-in chromatin) | Required when global scaling is biological |

### Global Accessibility Shift -- The Normalization Trap

**Trigger:** Treatment causes global chromatin compaction (e.g., HDAC inhibitor, DNMT inhibitor).

**Mechanism:** Full library-size normalization is robust to background but the default still scales background reads in; under uniform global compaction the magnitudes can collapse. RiP scaling (`library=DBA_LIBSIZE_PEAKREADS`) makes the opposite assumption (peak signal is stable, background absorbs the shift) and so erases the very biology of interest.

**Symptom:** Volcano plot is symmetric about zero; PCA shows treatment effect that vanishes after normalization.

**Fix:** Use spike-in normalization (add exogenous chromatin pre-Tn5; scale by spike-in reads), or keep `library=DBA_LIBSIZE_FULL` (script default) but interpret with the global shift in mind. Reske 2020 documented that the normalization choice materially changes which peaks are called differential under such a global shift.

## Per-Tool Failure Modes

### DiffBind -- Library-size choice confounds global change

**Trigger:** Treatment causes whole-genome accessibility shift; cell-cycle synchronized samples; differentiation timecourse.

**Mechanism:** Setting `library=DBA_LIBSIZE_PEAKREADS` (RiP-based) assumes total reads-in-peaks is comparable across samples. Global accessibility shifts break this assumption.

**Symptom:** Conditions clearly differ in PCA before normalization; after normalization PC1 is nearly noise.

**Fix:** Keep `library=DBA_LIBSIZE_FULL` (script default; or use spike-in scaling) and re-run `dba.contrast` and `dba.analyze`. Re-inspect PCA; if treatment now drives PC1, the global-shift biology is preserved.

### DiffBind summits parameter -- Width-driven differential

**Trigger:** Per-rep peaks have very different widths; consensus uses union.

**Mechanism:** Without `summits=250`, DiffBind counts reads in the original peak intervals. A peak called as 200 bp in one rep and 800 bp in another inflates the count for the wider rep.

**Symptom:** Top differential peaks track peak width, not signal intensity.

**Fix:** Always set `summits=250` (or 100, depending on resolution). This re-centers all peaks on the summit and uses identical 501 (or 201) bp windows.

### csaw -- Window size and filter choice dominates results

**Trigger:** Default `width=spacing=50` bp windows; default `filter=10` count cutoff.

**Mechanism:** Narrow windows have very low counts and inflated dispersion; the global background filter discards too many windows. Results are extremely sensitive to these.

**Symptom:** Number of significant windows ranges from 200 to 200,000 across reasonable parameter sweeps.

**Fix:** Use `width=150` for ATAC (matches typical NFR fragment); threshold with `filterWindowsGlobal(data, background)$filter > log2(3)` to discard low-signal windows. Validate by running on technical replicates -- ~0 differential windows is the expected outcome.

### DESeq2 -- Apeglm shrinkage with too few reps

**Trigger:** n=2 per condition; using `lfcShrink(type='apeglm')`.

**Mechanism:** Apeglm shrinks log2FC toward zero based on dispersion estimate; at n=2 dispersion is unreliable, shrinkage is over-aggressive, and biology is masked.

**Symptom:** All log2FC values cluster near zero post-shrinkage; FDR list has high p-values across the board.

**Fix:** Skip shrinkage at n=2 OR switch to edgeR QL test. If n=2 is unavoidable, report unshrunken log2FC alongside FDR; do not use shrunken FCs as the effect size.

### edgeR QL -- Filter must be aggressive enough

**Trigger:** Including peaks with mean count < 5 across all samples.

**Mechanism:** The QL F-test calibrates dispersion across all features. Including very-low-count peaks pulls dispersion estimates and inflates FDR.

**Fix:** `filterByExpr(y, group=group)` removes low-count peaks; restore peaks one at a time only if they are biologically critical and supported by at least one rep at depth.

## Reconciliation: When Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| DiffBind + DESeq2 differ wildly | Different normalization (DiffBind default = `DBA_NORM_DEFAULT`, native RLE/TMM on reads-in-peaks; script default = full library; plain DESeq2 = RLE on all counts) | Force same normalization; differences should shrink |
| DiffBind + csaw differ | csaw catches diffuse changes peaks miss; DiffBind catches narrow peaks csaw smooths | Both can be correct; report intersection as high-confidence |
| Top hits in DiffBind have FDR > 0.5 in DESeq2 | DiffBind's blacklist filter or width re-centering changes the per-region count | Re-run DESeq2 on the exact DiffBind consensus matrix (`dba.peakset` extract) |
| Effect-size ranking differs across reps | One rep is an outlier -- check PCA | Drop or block as covariate; never silently include |
| No significant peaks despite obvious browser-track differences | Library-size normalization eaten the global shift | Switch to spike-in or full-library normalization |

**Operational rule:** For high-confidence reporting, require concordant detection in two methods from different families (DiffBind/DESeq2-style on consensus peaks AND csaw-style sliding windows agreeing within +/- 500 bp). Report the intersection as primary; the union as exploratory.

## Effect Size and Threshold Selection

| Question | Threshold | Rationale |
|----------|-----------|-----------|
| Statistical significance | FDR < 0.05 | Standard BH FDR (DESeq2 / edgeR / DiffBind default) |
| Stringent biological change | abs(log2FC) >= 1 (= 2-fold) | Within-noise effects below 2-fold are unreliable in chromatin |
| Conservative reporting | FDR < 0.01 AND abs(log2FC) >= 1 | Per ENCODE differential reporting guidance |
| Exploratory / discovery | FDR < 0.1 OR shrunken log2FC >= 0.585 (1.5x) | For follow-up validation, not final claim |
| Proper effect-size reporting | Use shrunken log2FC (apeglm or DESeq2 lfcShrink) when n >= 3 | Raw log2FC at low counts is volatile |

abs(log2FC) >= 1 is *not* universal. ATAC effects in primary cells (immune subsets, neurons) often max at 1.5-fold; require log2FC >= 0.585 with FDR < 0.05 for those settings.

## Hidden Batch with SVA / RUVseq

**Goal:** Recover differential accessibility signal when unknown batch effects swamp the contrast.

**Approach:** Estimate surrogate variables via svaseq, append them to the design, refit, and extract the contrast. `scripts/diff_accessibility.R --sva=N` does this in DESeq2 directly on the DiffBind counts (limits: SKILL.md step 5); the block below is the standalone DESeq2 recipe. `n.sv` must not exceed `n_samples - ncol(mod) - 1` (4 samples: 1), or svaseq returns NaN.

```r
library(DESeq2); library(sva)

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized=TRUE)
dat <- dat[rowMeans(dat) > 1, ]

mod  <- model.matrix(~condition, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
nsv <- min(2, ncol(dat) - ncol(mod) - 1)
svobj <- svaseq(dat, mod, mod0, n.sv=nsv)

colData(dds) <- cbind(colData(dds), setNames(as.data.frame(svobj$sv), paste0('SV', seq_len(nsv))))
design(dds) <- as.formula(paste('~', paste0('SV', seq_len(nsv), collapse=' + '), '+ condition'))
dds <- DESeq(dds)
res <- results(dds, contrast=c('condition', 'treated', 'control'))
```

RUVseq is the alternative (recipe only, no shipped script) when negative-control regions (ChrM peaks NOT changing) or technical replicates are available. SVA is preferred when no controls exist.

## Spike-in Normalization

**Recipe only, no shipped script.**

**Trigger:** Treatment causes whole-genome accessibility shift (HDAC inhibitor, DNMT inhibitor); RPM/CPM/RiP normalization erases the global biology.

**Mechanism:** Exogenous chromatin spike-in (Drosophila S2 nuclei) is added at constant cell number ratio pre-Tn5; reads aligning to dm6 quantify the constant exogenous baseline. Sample-level scaling factor = inverse of dm6 reads per sample, applied to human-aligned counts.

**Goal:** Preserve global accessibility shifts that RiP / library-size normalization would erase.

**Approach:** Compute per-sample size factors from inverse spike-in read counts, override DESeq2's default size factors, then run the standard DESeq2 fit and contrast.

```r
library(DESeq2)

# spike_counts: per-sample dm6 read counts (one column per sample)
sf_spike <- 1 / spike_counts
sf_spike <- sf_spike / mean(sf_spike)                 # Geometric mean = 1 for stability

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
sizeFactors(dds) <- sf_spike                          # Override library-size factors
dds <- DESeq(dds)
res <- results(dds, contrast=c('condition', 'treated', 'control'))
```

After spike-in normalization, log2FC reflects absolute accessibility change (not just relative redistribution). Spike-in is the most direct control for global-shift biology.
