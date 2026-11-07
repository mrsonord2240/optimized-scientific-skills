## IsoformSwitchAnalyzeR Workflow

**Goal:** Identify isoform switches with functional consequences in one integrated workflow.

**Approach:** Import Salmon (raw counts), join the design to the quantification by sample name, pre-filter, run the statistical test, then annotate consequences (next section).

```r
library(IsoformSwitchAnalyzeR)
set.seed(1)  # importRdata() may call stochastic sva when it detects unwanted effects

salmonQuant <- importIsoformExpression(
    parentDir = 'salmon_quant/',              # one sub-directory per sample, each holding quant.sf
    calculateCountsFromAbundance = FALSE,     # raw NumReads as counts; see "Count route"
    addIsofomIdAsColumn = TRUE                # (sic) the package's own spelling
)

# Join the design to the quantification BY SAMPLE NAME. importIsoformExpression sorts samples
# alphabetically, so a hand-typed condition vector silently mislabels samples (on SRR-style IDs it
# recovered 0/20 planted switches with no warning; joined by name, 20/20).
meta <- read.delim('sample_metadata.tsv')     # columns: sample_id (= quant sub-directory name), condition, [batch]
ids <- setdiff(colnames(salmonQuant$counts), 'isoform_id')
if (!setequal(ids, meta$sample_id))
    stop('sample IDs differ between quantification and metadata.\n  only in quant: ', paste(setdiff(ids, meta$sample_id), collapse = ', '),
         '\n  only in metadata: ', paste(setdiff(meta$sample_id, ids), collapse = ', '))
stopifnot('duplicated sample_id in metadata' = !anyDuplicated(meta$sample_id))
design <- data.frame(sampleID = ids, condition = meta$condition[match(ids, meta$sample_id)])
# design$batch <- meta$batch[match(ids, meta$sample_id)]   # extra columns = covariates, see below
print(table(design$condition))

aSwitchList <- importRdata(
    isoformCountMatrix = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix = design,
    isoformExonAnnoation = 'annotation.gtf',   # (sic); a GTF with CDS lines also supplies annotated ORFs
    isoformNtFasta = 'transcripts.fa',
    addAnnotatedORFs = TRUE,
    showProgress = FALSE
)

aSwitchList <- preFilter(
    aSwitchList,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    IFcutoff = 0.01,
    removeSingleIsoformGenes = TRUE,
    keepIsoformInAllConditions = TRUE
)

# DEXSeq for <=5 replicates per condition, satuRn above (see Tool Selection).
# reduceToSwitchingGenes = FALSE: with TRUE (the function default) a run with no switches stops with an error
aSwitchList <- if (max(table(design$condition)) > 5) {
    isoformSwitchTestSatuRn(aSwitchList, reduceToSwitchingGenes = FALSE, alpha = 0.05, dIFcutoff = 0.1, diagplots = FALSE)
} else {
    isoformSwitchTestDEXSeq(aSwitchList, reduceToSwitchingGenes = FALSE, alpha = 0.05, dIFcutoff = 0.1)
}
f <- aSwitchList$isoformFeatures
sum(f$isoform_switch_q_value < 0.05 & abs(f$dIF) > 0.1, na.rm = TRUE)   # 0 -> see "Count route"
```

`scripts/isoform_switch_analysis.R` runs this whole workflow, including consequences and the plot, on a built-in synthetic dataset (`Rscript isoform_switch_analysis.R`) or on your own Salmon directory, GTF, FASTA and sample-metadata table.

`preFilter` parameters:
- `geneExpressionCutoff = 1` — minimum TPM for gene to be tested (raise for stricter)
- `isoformExpressionCutoff = 0` — minimum TPM per isoform (set to 1 for stricter)
- `IFcutoff = 0.01` — minimum isoform fraction; below = noise
- `removeSingleIsoformGenes = TRUE` — drop genes with only one detectable isoform (cannot have DTU)
- `keepIsoformInAllConditions = TRUE` — require expression across all conditions

**Replicates and covariates.** IsoformSwitchAnalyzeR stops on 1 replicate per condition ("A statistical test cannot be performed without replicates"). 2 per condition runs but is not trustworthy on the raw-count route (see "Count route": label-shuffled 2 v 2 splits of real data are called almost as often as the true split); plan >=3, ideally >=5. Any extra `designMatrix` column is treated as a covariate: on a batch-confounded synthetic set a batch column removed the batch artefacts (19/30 artefact genes called without it, 0 with it) at no cost to the 20/20 planted switches; a covariate identical to condition stops with "not full rank", and a constant one with "Contain constant information". `importRdata(detectUnwantedEffects = TRUE)` (the default) also runs stochastic `sva` and adds any surrogate variables it finds to the model. Set a seed before `importRdata()` (the workflow uses `set.seed(1)`) and record it with the results: on one label-shuffled chrX 2 v 2 split, `sv1` made the call count vary from 0 to 11 across seeds. Check `colnames(aSwitchList$designMatrix)` and treat such small-n lists as exploratory. Per `?importRdata` and the package NEWS, abundance, IF and dIF are then corrected for these covariates (counts stay unmodified but the test includes them), so with a covariate the reported dIF is no longer the plain mean of per-sample isoform fractions (up to 0.08 from a hand computation on the synthetic 3 v 3 set with a batch column; 7e-5 without).

For long-read data pass the long-read transcript count matrix to `importRdata` directly (a count-only import recovered 20/20 planted switches on the synthetic set) — no Salmon EM uncertainty.

### Count route (`calculateCountsFromAbundance`)

`importIsoformExpression()` derives counts from abundance by default (`calculateCountsFromAbundance = TRUE` = tximport `scaledTPM`); the package vignette recommends this because it carries Salmon's bias correction into the counts. The tests use those counts, while dIF comes from the abundance in either route. **On real chrX data this default called nothing in every label split (below), so it cannot be validated there. This Skill uses raw NumReads (`FALSE`), which does find switches but is not validated either: label-shuffled splits of the same data are called almost as often as the true one. Read "Null check" before trusting a call.**

Real chrX RNA-seq, 2 GBR v 2 YRI (nf-core rnasplice, Salmon 2.7.0, GRCh37), `isoformSwitchTestDEXSeq`, switching isoform = q < 0.05 and |dIF| > 0.1:

| Count route | Switching isoforms / genes | RPL10 ENST00000406022 (dIF -0.33) q |
|-------------|----------------------------|-------------------------------------|
| `calculateCountsFromAbundance = FALSE` (raw NumReads) | 26 / 20 | 4.6e-18 |
| tximport `lengthScaledTPM` | 28 / 22 | 1.9e-18 |
| default (`scaledTPM`) | 0 / 0 | 1 |
| tximport `dtuScaledTPM` | 0 / 0 | 1 |

- **Cause (consistent with, not proven by, the runs below):** scaledTPM counts are proportional to TPM, i.e. each isoform's counts are divided by its length. The factor median(length)/length ranges 0.019-576 across the 5,895 isoforms with `EffectiveLength` (0.023-18 with `Length`), so short isoforms are inflated and stop behaving like counts, and DEXSeq's isoform-versus-rest test degrades. Multiplying the raw counts by either factor (`Length` or `EffectiveLength`) reproduced the result (0 switches, RPL10 q = 1); multiplying by the same factors shuffled across isoforms did not (34 switches / 27 genes with `Length`, 31 / 23 with `EffectiveLength`). satuRn was not compared across routes on real data; the manual-pipeline route numbers are under "Manual DTU Pipeline".
- The synthetic planted set has near-uniform isoform lengths and recovered 20/20 with either route, so a clean simulation will not reveal the problem.

### Null check: shuffle the condition labels

The four chrX samples are 2 GBR + 2 YRI. Splits that put one GBR and one YRI in each group carry no population signal, so any calls there are individual-level differences that the n = 2 dispersion estimate does not absorb. Genes with a switching isoform (q < 0.05, |dIF| > 0.1), raw counts unless stated:

| Split | ISAR DEXSeq | ISAR default route | manual DEXSeq gene q < 0.05 | manual DRIMSeq | manual both | ISAR calls also DEXSeq-confirmed |
|-------|-------------|--------------------|-----------------------------|----------------|-------------|----------------------------------|
| true GBR v YRI | 20 | 0 | 13 | 11 | 7 | 12 |
| mixed 1 | 11 | 0 | 0 | 4 | 0 | 0 |
| mixed 2 | 15 | 0 | 4 | 3 | 1 | 4 |

- **The raw route is not validated at n = 2.** The mixed splits are called at 55% and 75% of the true split's count, and only 1 of the 20 true-split genes recurs in either. No cutoff separates them cleanly (genes, true / mixed 1 / mixed 2: q < 0.01: 12 / 3 / 7; q < 1e-4: 8 / 2 / 5; q < 1e-6: 5 / 2 / 1; |dIF| > 0.3: 17 / 7 / 10; a higher gene-TPM or isoform-TPM filter changed nothing). Report a 2 v 2 list as exploratory. The default route's zeros show only that it has no power on this data, not that it is the safe choice.
- **What reduced it here:** keep only calls that the manual gene-level DEXSeq test above also makes (12 of 20 true-split genes, 0 and 4 of the mixed-split calls; DRIMSeq as well: 6, 0, 1). It cuts the shuffled calls but does not remove them (4 of 15 remain), and the counts are from one 4-sample data set.
- **Run the permutation check on every real data set** (this block reuses `salmonQuant` and `design` from the workflow; each call repeats the import, about 20-40 s on chrX). Trust an observed list only when it is far above every permuted count. Group sizes of 2 give only two alternative splits (no p-value); 3 v 3 gives nine, 6 v 6 hundreds. Splits that share half of the true split keep part of the signal, so they overstate the null (conservative).

```r
# Same import, filter and test as the workflow, on any design; returns the genes called
call_genes <- function(design) {
    sl <- importRdata(isoformCountMatrix = salmonQuant$counts, isoformRepExpression = salmonQuant$abundance,
                      designMatrix = design, isoformExonAnnoation = 'annotation.gtf', isoformNtFasta = 'transcripts.fa',
                      addAnnotatedORFs = FALSE, showProgress = FALSE, quiet = TRUE)
    sl <- preFilter(sl, geneExpressionCutoff = 1, isoformExpressionCutoff = 0, IFcutoff = 0.01,
                    removeSingleIsoformGenes = TRUE, keepIsoformInAllConditions = TRUE, quiet = TRUE)
    sl <- if (max(table(design$condition)) > 5) {
        isoformSwitchTestSatuRn(sl, reduceToSwitchingGenes = FALSE, alpha = 0.05, dIFcutoff = 0.1, diagplots = FALSE, quiet = TRUE)
    } else {
        isoformSwitchTestDEXSeq(sl, reduceToSwitchingGenes = FALSE, alpha = 0.05, dIFcutoff = 0.1, quiet = TRUE)
    }
    f <- sl$isoformFeatures
    unique(f$gene_id[!is.na(f$isoform_switch_q_value) & f$isoform_switch_q_value < 0.05 & abs(f$dIF) > 0.1])
}

# Shuffled labels that keep the group sizes and share only the chance overlap with the true split
lab <- design$condition; n <- length(lab); in1 <- lab == lab[1]; k <- sum(in1)
canon <- function(v) paste(v == v[1], collapse = '')          # a split and its mirror image are the same split
seen <- canon(in1); perms <- list(); set.seed(1)
for (i in 1:5000) {
    if (length(perms) >= 10) break
    p <- sample(lab); key <- canon(p == lab[1])
    if (abs(sum(p == lab[1] & in1) - k * k / n) <= 0.5 && !key %in% seen) { seen <- c(seen, key); perms[[length(perms) + 1]] <- p }
}
observed <- length(call_genes(design))
permuted <- sapply(perms, function(p) { d <- design; d$condition <- p; tryCatch(length(call_genes(d)), error = function(e) NA) })
cat('observed', observed, 'genes; label-permuted:', permuted, '\n')
```

Checked with this block: real chrX 2 v 2 observed 20, permuted 15 and 11. Synthetic 6 v 6 (satuRn branch, batch covariate kept): observed 24 (20/20 planted, 0 null-type genes), permuted 0 0 0 1 0 0 0 0 0 0 over 10 splits. Synthetic 3 v 3 (DEXSeq): observed 25 (20/20 planted, 2 null-type genes), permuted 0 2 0 0 0 0 0 0 0 over the 9 other splits.

**Recommended route:** raw NumReads (`FALSE`); at least 3 replicates per condition, ideally 5 or more; run the permutation check; confirm calls with the manual gene-level test; if you also try `TRUE`, compare the two switch lists and inspect the top switches (`switchPlot`). A large gap between routes means the count route, not biology, drove the result.
