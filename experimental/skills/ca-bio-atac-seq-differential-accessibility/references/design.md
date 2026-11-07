# Design



Differential accessibility workflows differ by input type, statistical model, and handling of peak structure. Choose based on your replicate count, whether peaks are consistent across conditions, and whether global accessibility shifts are biological.

| Tool | Model | Input | Min reps | Strength | Fails when |
|------|-------|-------|----------|----------|------------|
| DiffBind 3.x (default DESeq2) | NB GLM via DESeq2 on consensus peaks | BAM + peak files | 2-3 per group | ATAC-aware defaults; built-in QC; blocking factors. DiffBind's own default is `DBA_NORM_DEFAULT` (RLE/TMM on a reads-in-peaks library, verified on 3.16); the shipped script defaults to full library size (`lib` + `DBA_LIBSIZE_FULL`, background-included) | Peaks differ dramatically between conditions (closed -> open shifts width); fewer than 2 reps per group |
| DiffBind with edgeR backend | NB GLM via edgeR-QL on consensus peaks | Same | 2-3 per group | Robust at low replicates (n=2 OK); QL test calibrates dispersion better than DESeq2 at small n | When global accessibility shifts dominate, switch to spike-in or full-library (library=DBA_LIBSIZE_FULL), never reads-in-peaks |
| DESeq2 directly on peak counts | NB GLM with shrinkage | featureCounts SAF | 3+ | Maximum control; integrates with apeglm shrinkage; modern interface | Need to manually build consensus peakset; per-region pre-filter required (low counts inflate dispersion) |
| edgeR QL F-test on peak counts | NB QL (quasi-likelihood) | featureCounts | 2 | Calibrated FDR at low n (n=2 viable); robust to outlier reps | Manual consensus peakset; small library bias unless normalization explicit |
| csaw (windows) | edgeR-QL on sliding windows | BAM only | 2 | No peak set required; detects diffuse changes peaks miss; merges adjacent windows | Computationally heavy; window size choice biases results; harder to annotate downstream |
| limma-voom | linear model with mean-variance trend | log2(CPM+offset) | 3 | Fast; good calibration at moderate count | Mis-calibrated at very low counts (atac peaks often have dropouts); needs explicit voom normalization |

Methodology evolves; verify the current consensus practice (Gontarz 2020 DA-strategy benchmark; Reske 2020 normalization comparison) before locking pipelines.

## Decision Tree by Experimental Scenario

| Scenario | Recommended workflow | Why |
|----------|---------------------|-----|
| 3+ reps, similar peak structure across conditions | DiffBind (DESeq2 backend), `summits=250`, normalization per the table below (script default `DBA_NORM_LIB`) | Standard pattern; peak-level inference is interpretable |
| 2 reps per condition | DiffBind with edgeR backend OR raw edgeR QL | DESeq2 underpowered at n=2; QL is robust |
| Peak structure differs dramatically (e.g., differentiation, KO of pioneer TF) | csaw windows OR rebuild consensus peakset post-hoc per condition then take union | Stable consensus peakset is invalid when chromatin landscape shifts |
| Multi-factor design (batch, sex, time) | DiffBind with `design='~Tissue + Condition'` (batch stored in the `Tissue` sample-sheet column; DiffBind rejects other names) | Standard linear model adjustment |
| Hidden batch / unknown variance | DESeq2 + SVA OR RUVseq before fitting | Empirical surrogate variables capture unknown nuisance |
| Long timecourse (5+ time points) | DESeq2 LRT (likelihood ratio test) on `~time + condition + time:condition` | Captures temporal interaction; use `time-course differential expression workflow (no matching shelf skill)` patterns |
| Diffuse / broad accessibility change (super-enhancers) | csaw with merged windows OR call broad peaks first | Narrow peaks fragment broad domains -> inflated peak count, deflated effect |
| Single-cell ATAC pseudobulk | DESeq2 on aggregated counts OR Signac::FindMarkers | See bio-atac-seq-single-cell-atac |
| Allele-specific accessibility | csaw on heterozygous SNPs OR HOMER tagDir | Peak-level invalid because alleles share peaks |
| Plant / non-model organism | DiffBind works; run the script with `--txdb=none` (or another organism's TxDb package) | Annotation defaults to hg38 UCSC TxDb; `none` writes an unannotated CSV |

## Consensus Peak Set Strategy

The consensus peakset choice drives FDR calibration. DiffBind defaults rarely match what a chromatin biologist wants.

| Strategy | Implementation | When to use |
|----------|---------------|-------------|
| Intersection (peak in all reps) | `dba.count(minOverlap=N)` with N = total reps | Strict; for high-confidence reproducible analysis (matches IDR philosophy) |
| Union (peak in any rep) | `minOverlap=1` | Maximum sensitivity; risks single-rep artefact peaks |
| Majority rule (peak in >= half reps) | `minOverlap=ceiling(N/2)` | DiffBind default-ish; balance |
| Per-condition union, then union of unions | Compute consensus per group, then merge | Best when conditions have very different peak counts |
| Iterative overlap removal (Corces 2018) | Sort peaks by significance; greedily keep non-overlapping; fixed-width 501 bp | Standard for fixed-width consensus; required for peak-count matrices used in machine learning |

Refer to bio-atac-seq-consensus-peakset for full coverage of fixed-width re-centering and the iterative overlap algorithm. For DiffBind, the key parameter is `summits=250` (re-center peaks on summit +/- 250 bp = 501 bp fixed width).



```csv
SampleID,Condition,Replicate,Tissue,bamReads,Peaks,PeakCaller
ctrl_rep1,control,1,A,ctrl_rep1.bam,ctrl_rep1_peaks.narrowPeak,macs
ctrl_rep2,control,2,B,ctrl_rep2.bam,ctrl_rep2_peaks.narrowPeak,macs
treat_rep1,treated,1,A,treat_rep1.bam,treat_rep1_peaks.narrowPeak,macs
treat_rep2,treated,2,B,treat_rep2.bam,treat_rep2_peaks.narrowPeak,macs
```

`Tissue` carries batch here. Covariates must occupy `Tissue`, `Factor`, `Treatment` or `Replicate`; other column names cannot enter a DiffBind design.




**Goal:** Assign each differentially accessible peak to its nearest gene and feature class for downstream interpretation.

**Approach:** Pull DiffBind / DESeq2 results as GRanges, annotate via ChIPseeker against a TxDb with a custom promoter window, then plot annotation distribution and extract gene IDs for enrichment.

```r
library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene)

diff_peaks <- dba.report(dba)
peakAnno <- annotatePeak(diff_peaks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                         tssRegion=c(-2000, 500), level='gene')
plotAnnoPie(peakAnno); plotDistToTSS(peakAnno)
genes <- as.data.frame(peakAnno)$geneId          # for GO enrichment via bio-pathway-go-enrichment
```

`tssRegion=c(-2000, 500)` defines promoter as TSS-2kb to TSS+500bp; ChIPseeker default (`-3000, 3000`) over-counts promoter assignments. Adjust per cell type / organism.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| DiffBind very slow | Counting all peaks across all BAMs sequentially | `dba.count(..., bParallel=TRUE)` and provide `BPPARAM` |
| `unable to use the provided design matrix` (DESeq2) | Confounded design (e.g., batch perfectly aligns with condition) | Replicate in a way that breaks the confound, or drop the batch term |
| FDR list empty despite obvious differences | RiP scaling (`library=DBA_LIBSIZE_PEAKREADS`) removed global biology | Use spike-in or keep `library=DBA_LIBSIZE_FULL` (script default); verify with browser tracks |
| Top peaks all on chrM | chrM not removed from BAM before counting | Always strip chrM upstream |
| dispersion estimate failure (DESeq2) | Too few peaks pass filter; too few reps | `filterByExpr` less aggressively; check rep count |
| `Error in if (any(out))`(csaw) | Window count below threshold | Reduce `bin.size`; check BAM is paired-end |
| ChIPseeker error on non-human TxDb | Wrong organism db loaded | Use `make_org_db` from biomartr or AnnotationDbi for non-model |
| Volcano plot symmetric about zero with no significant peaks | Hidden batch swamping signal | Run SVA/RUVseq |
