# ATAC-seq QC: Methods, Failure Modes, and Diagnostics

## Per-Metric Failure Modes

### Mitochondrial fraction > 50%

**Trigger:** Standard ATAC-seq protocol on intact cells (no nuclear isolation), or insufficient detergent in lysis.

**Mechanism:** Mitochondrial DNA is naked (no histones), so Tn5 hyperactively cuts it. Without nuclear-isolation steps (Omni-ATAC pre-spin, OR digitonin lysis with mt removal), chrM dominates the library.

**Symptom:** `samtools idxstats sample.bam | awk '$1=="chrM"'` shows >50% of mapped reads on chrM.

**Fix:** Re-prep with Omni-ATAC (Corces 2017) or fast-ATAC. Re-running QC on chrM-stripped BAM hides the underlying problem; the wasted sequencing remains. If chrM fraction is 30-50%, the library may still be salvageable via chrM removal but yield is reduced.

### NRF / PBC1 / PBC2 below threshold

**Trigger:** Over-amplified library; low input cell count combined with high PCR cycles.

**Mechanism:** Each PCR cycle doubles starting fragments. With low complexity input (<5000 cells) and >12 cycles, distinct fragments saturate and reads pile up at identical positions. NRF measures unique fragments / total; PBC2 specifically detects multi-copy duplication.

**Symptom:** NRF < 0.7; PBC2 < 1.0; massive duplicate-removal loss in `samtools markdup`.

**Fix:** No fix post-hoc. Re-prep with more starting cells and fewer PCR cycles. Note: ATAC has *legitimate* duplicates at hyperaccessible sites (Tn5 cuts identically there), so NRF < 0.9 is not by itself fatal. The combined PBC1 < 0.7 + PBC2 < 1.0 + visual coverage pile-ups confirm true bottlenecking.

### TSS enrichment < 5

**Trigger:** Generic chromatin opening throughout the genome (over-transposition), OR genome build mismatch between TSS BED and BAM, OR strand-flip in TSS file.

**Mechanism:** TSS enrichment requires that signal at TSSs is >> signal in genomic flanks. Over-transposition flattens the signal landscape. Strand-flipped TSSs subtract real signal because TSSs on - strand are calculated from the wrong direction.

**Symptom:** TSS profile is flat or shows a slight dip at TSS center. Genome browser shows accessibility everywhere, not concentrated at promoters.

**Fix:** Verify genome build (mm10 vs mm39 differ in TSS positions); verify GTF strand column; confirm signal track was generated post-deduplication. If TSS profile is genuinely flat, library is over-transposed and not recoverable; lower transposition time / Tn5 concentration in next prep.

### FRiP < 0.2

**Trigger:** Signal too diffuse to call peaks (over-transposition), low TSS enrichment, OR peak set is too narrow / restrictive.

**Mechanism:** FRiP correlates with TSS enrichment because both measure how concentrated the signal is. A diffuse library will have low FRiP regardless of peak count.

**Symptom:** Peak count looks normal but FRiP < 0.15.

**Fix:** Check TSS enrichment first. If TSS is also low, the library is over-transposed. If TSS is OK but FRiP is low, the peak caller may be undercalling -- try `-p 0.01` (looser) and recalculate FRiP.

### Replicate correlation < 0.85

**Trigger:** Batch effect, technical artefact, or cell-state drift between replicate biological collections.

**Mechanism:** Spearman correlation of binned counts (deepTools `multiBamSummary bins -bs 10000`, `--skipZeros`) tracks coverage similarity. Below 0.85 indicates non-trivial divergence; >= 0.9 is a common working target for biological replicates (not an ENCODE ATAC threshold).

**Fix:** Check PCA; if reps cluster apart from condition, drop the outlier or rerun. If the divergence aligns with batch, add batch as a covariate downstream (DiffBind `~Batch + Condition`). Do not silently merge with bad correlation.

## Library Complexity (NRF, PBC1, PBC2)

**Goal:** Detect over-amplification or low-input bottlenecks.

**Approach:** Key each paired-end fragment by (chrom, start, end) (single-end: 5' position and strand), tally how many keys have 1, 2, or more observations, and compute the three metrics with `scripts/library_complexity.py`. Per-mate read-start counting undercounts complexity (unfiltered GM12878 slice: 0.35 by mate start vs 0.62 by fragment).

**Computation Details:**
- **NRF (Non-Redundant Fraction):** Distinct fragments / total fragments
- **PBC1 (PCR Bottlenecking Coefficient 1):** Fragments seen exactly once / distinct fragments
- **PBC2:** Fragments seen once / fragments seen twice (undefined when none is seen twice; reported as null)

**Important:** Compute on the raw mapped BAM (MAPQ >= 30 and chrM excluded, both applied by the script; before deduplication). Computing post-dedup gives NRF = 1.0 trivially because every position becomes unique by definition.

## Sex-Chromosome QC

**Trigger:** Clinical-grade ATAC; biobank-scale studies; sample-mix-up detection.

**Mechanism:** chrY has minimal coverage in female samples; XIST locus (chrX) is highly accessible only in female cells (X-inactivation). Sample-swap or sex-misassignment detectable from these two loci.

```bash
# chrY read fraction
samtools idxstats sample.bam | awk '$1=="chrY"{print $3 / $2}'   # reads per bp

# XIST locus accessibility (chrX:73820651-73852753 in hg38)
samtools view -c sample.bam chrX:73820651-73852753
```

Female: chrY reads/bp ~0; XIST count high. Male: chrY reads/bp ~male coverage; XIST count low. Discrepancy with sample metadata flags swap.

## Cell-Cycle Effect on Accessibility

**Trigger:** Proliferating cell lines (K562, HEK293, HeLa); samples with high S/G2M signature.

**Mechanism:** Replication-associated chromatin opening adds 5-15% global accessibility shift in proliferating cells; without correction, condition-specific cell-cycle differences confound differential analysis.

**Detection:** Score cells/samples for S-phase signature (Macosko 2015 cell cycle gene set adapted for chromatin: regulated origin loci, replication-stress-response genes); for bulk ATAC, compute per-sample peak intersection with replication-origin atlas (Repli-seq peaks).

**Fix for differential:** Add S-phase score as covariate in DESeq2 design (`~Sphase + Condition`); for scATAC, regress on TF-IDF residuals analogous to Seurat CellCycleScoring.

## Spike-in QC (Drosophila or E. coli Chromatin)

**Trigger:** Studies where global accessibility shift is biological (HDAC inhibitor, DNMT inhibitor, differentiation).

**Mechanism:** Per-library normalization (RPM, CPM) erases global accessibility shifts because total reads are nominally constant. Exogenous chromatin spike-in (Drosophila S2 or E. coli Tn5-naive chromatin added pre-Tn5) provides an external scaling reference.

**Pipeline:** Align reads to a concatenated human + Drosophila reference; count spike-in reads per sample; normalize by spike-in (not by total reads). Reske 2020 Epigenetics Chromatin shows that normalization-method choice materially changes differential-accessibility results when a global accessibility shift is expected (ARID1A/PIK3CA endometrial-epithelium case study), motivating an external reference such as a chromatin spike-in.

**QC threshold:** spike-in fraction 0.5-5% of total reads is the workable range. Below 0.1% spike-in is unreliable; above 10% suggests too much spike-in (loss of cellular reads).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| TSS enrichment off by 3x from expected | Wrong implementation (ENCODE vs ATACseqQC) | State the formula; convert by recomputing |
| NRF = 1.0 exactly | BAM was already deduplicated -> all positions distinct | Compute NRF on raw mapped BAM (pre-dedup) |
| PBC2 = null | No fragment seen twice | Library is too sparse; PBC2 unreliable below ~5M reads |
| Mt fraction reported but BAM has no `chrM` | Mitochondrial chromosome named `MT`, `Mt`, or `chromosome:MT` | Match `samtools idxstats` chromosome name to the filter |
| Insert size distribution flat after Picard | Sample is single-end | Insert size only valid for paired-end; switch to deeptools fragmentSize |
| Replicates correlate poorly but PCA looks fine | High background dominates correlation | Use `--skipZeros`; or compute correlation on peak counts only |
| FRiP differs by 2x between identical pipeline runs | Peak set differs (q-value cutoff drift) | Pin caller version + cutoff; FRiP is peak-set-dependent |
| TSS enrichment lower than expected on Omni-ATAC | Used standard TSS BED on FFPE-prepped sample | FFPE TSSs are degraded; use peak-based metric instead |
