# limma_pipeline.R -- MaxQuant proteinGroups.txt -> differentially abundant proteins (LFQ, limma treat()).
# Purpose : the default route of bio-workflows-proteomics-pipeline: strip contaminant/reverse rows, log2,
#           inspect RAW per-sample distributions (drop failed loads), median-normalize, per-group
#           completeness filter, QC, then limma treat() with batch as a covariate. Nothing is imputed.
# Inputs  : proteinGroups.txt (needs Intensity.* and LFQ.intensity.* columns);
#           sample_annotation.csv with columns sample, condition, replicate, batch
#           (`sample` must equal the intensity column suffix; first condition level = reference).
# Usage   : Rscript limma_pipeline.R [proteinGroups.txt] [sample_annotation.csv] [results.csv] [min_frac=0.6] [min_fold_change=1.5]
# Output  : results.csv (protein x contrast: logFC, AveExpr, t, P.Value, adj.P.Val, significant) and
#           results_raw_boxplot.pdf next to it.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
input_file      <- arg(1, 'proteinGroups.txt')
annotation_file <- arg(2, 'sample_annotation.csv')
out_file        <- arg(3, 'proteomics_results.csv')
min_frac        <- as.numeric(arg(4, 0.6))
min_fc          <- as.numeric(arg(5, 1.5))

library(limma)
library(ggplot2)
library(pheatmap)

# === 1. DATA IMPORT ===
# quote = '' / comment.char = '': MaxQuant text fields contain quotes and '#'; the defaults
# silently truncate the table with only an 'EOF within quoted string' warning. Check the count.
proteins <- read.delim(input_file, stringsAsFactors = FALSE, quote = '', comment.char = '')
stopifnot(nrow(proteins) == length(readLines(input_file)) - 1)
cat('Loaded', nrow(proteins), 'protein groups\n')

# Filter contaminants, reverse, only-by-site. Use %in%, NOT `!= '+'`: when a run flags none of a
# category the column is entirely empty, read.delim types it as logical NA, `NA != '+'` is NA, and
# subsetting by NA replaces every row with NAs while nrow() still looks right.
proteins <- proteins[!(proteins$Potential.contaminant %in% '+') &
                      !(proteins$Reverse %in% '+') &
                      !(proteins$Only.identified.by.site %in% '+'), ]
stopifnot(!all(is.na(proteins$Majority.protein.IDs)))
cat('After filtering:', nrow(proteins), 'proteins\n')

# Extract LFQ intensities
lfq_cols <- grep('^LFQ\\.intensity\\.', colnames(proteins), value = TRUE)
intensities <- proteins[, lfq_cols]
rownames(intensities) <- proteins$Majority.protein.IDs
colnames(intensities) <- gsub('LFQ\\.intensity\\.', '', colnames(intensities))

# === 2. LOG2 TRANSFORM, THEN INSPECT RAW DISTRIBUTIONS ===
intensities[intensities == 0] <- NA
log2_int <- log2(intensities)

# Inspect BEFORE normalizing (rule 4). Median-centering rescales every sample onto a common median,
# so it mathematically erases the 3x-low load that marks a failed injection -- after this point the
# failure is invisible. Identify and drop failures HERE.
# Inspect the RAW `Intensity.` columns, NOT `LFQ intensity`: MaxLFQ has ALREADY renormalized the LFQ
# columns onto a common scale, so a low load is largely gone from them before you ever look.
raw_cols <- grep('^Intensity\\.', colnames(proteins), value = TRUE)
raw_log2 <- log2(replace(proteins[, raw_cols], proteins[, raw_cols] == 0, NA))
colnames(raw_log2) <- gsub('^Intensity\\.', '', colnames(raw_log2))
raw_log2 <- raw_log2[, colnames(log2_int), drop = FALSE]
raw_median <- apply(raw_log2, 2, median, na.rm = TRUE)
id_counts <- colSums(!is.na(log2_int))
pdf(sub('\\.csv$', '_raw_boxplot.pdf', out_file))   # Rscript has no screen device
boxplot(raw_log2, las = 2, main = 'RAW log2 Intensity (pre-normalization)', ylab = 'log2 intensity')
dev.off()
print(data.frame(id_count = id_counts, raw_median_log2 = round(raw_median, 2),
                 load_shift_log2 = round(raw_median - median(raw_median), 2)))

# Flag on BOTH axes and as an OUTLIER, not a fixed fraction: a >=2x (1 log2) drop in raw median
# signal, OR an ID count more than 3 MADs below the cohort median. Each axis alone misses cases --
# a 3x-low injection loses its dimmest proteins entirely, so the SURVIVING median understates the
# deficit (a 3.0x low load can read as only -0.8 log2), while a 50%-of-median ID-count rule is so
# loose that a run at 82% of the median IDs, visibly low on both axes, passes it.
mad_ids <- mad(id_counts)
failed <- colnames(raw_log2)[(raw_median - median(raw_median)) <= -1 |
                             (mad_ids > 0 & id_counts < median(id_counts) - 3 * mad_ids)]
if (length(failed) > 0) {
    message('Dropping failed samples: ', paste(failed, collapse = ', '))
    log2_int <- log2_int[, !colnames(log2_int) %in% failed, drop = FALSE]
}

# === 3. NORMALIZE (only after the raw inspection above) ===
# ASSUMPTION: median centering assumes most proteins are unchanged AND that the changes that do
# occur are roughly SYMMETRIC up/down. When they are not -- an enrichment pulldown, a secretome, a
# strong one-sided drug response -- the shift it removes is real signal, and the whole unchanged
# proteome acquires the opposite offset. Measured on a set with 44 up vs 26 down among 296
# proteins: every true null moved -0.19 log2 (t vs 0, p = 2e-55) and 21 of 179 true nulls were
# called at BH 5%, all negative. Normalize on a set expected to be unchanged (spike-in standards,
# an internal-reference protein set, or the housekeeping bulk) when the design is one-sided,
# and always check the null centre afterwards -- see the Common Errors row.
sample_medians <- apply(log2_int, 2, median, na.rm = TRUE)
global_median <- median(sample_medians)
normalized <- sweep(log2_int, 2, sample_medians - global_median)

# === 4. FILTER ON PER-GROUP COMPLETENESS (do NOT impute by default) ===
# Filter FIRST on completeness PER GROUP: keep a protein if it is valid in >= ~50-70%
# of replicates in AT LEAST ONE condition. A protein missing in every group fails QC.
sample_info <- read.csv(annotation_file)
# Re-align the annotation to the samples that SURVIVED the raw-distribution QC above; otherwise the
# column indexing below requests a dropped sample and errors (or silently misaligns the design).
# match() also puts the rows in COLUMN order: lmFit pairs design rows with matrix columns
# POSITIONALLY, so an annotation sorted differently from the intensity columns fits a wrong design
# in silence. Never subset the annotation with `%in%` alone here.
sample_info <- sample_info[match(colnames(normalized), sample_info$sample), ]
stopifnot(!any(is.na(sample_info$sample)))
sample_info$condition <- droplevels(factor(sample_info$condition))
# min_frac (default 0.6: >= 60% present within at least one group; tune 0.5-0.7 per design) is an argument
group_complete <- sapply(levels(sample_info$condition), function(g) {
    cols <- sample_info$sample[sample_info$condition == g]
    rowSums(!is.na(normalized[, cols, drop = FALSE])) >= ceiling(length(cols) * min_frac)
})
valid_rows <- rowSums(group_complete) > 0
filtered <- normalized[valid_rows, ]
cat('Proteins after per-group completeness filter:', nrow(filtered), '\n')

# Missingness in label-free DDA is left-censored MNAR (missing BECAUSE low). The modern,
# correct approach is to MODEL the missingness in the likelihood, NOT impute it. See
# proteomics/differential-abundance for the decision (proDA / msqrob2 / MSstats-AFT).
# NOTHING BELOW IMPUTES: limma fits each protein on its OBSERVED values, so the completeness
# filter above is the only missing-value handling in the executed path.

# --- UPGRADE: model the dropout explicitly with proDA, then use `da` in place of `results` ---
# library(proDA)
# fit <- proDA(as.matrix(filtered), design = ~ condition, col_data = sample_info,
#              reference_level = 'Control')            # add batch: design = ~ batch + condition
# da <- test_diff(fit, contrast = 'conditionTreatment')   # columns: diff (log2FC), pval, adj_pval

# --- DO NOT: left-censored downshift imputation, then limma ---
# Downshift MANUFACTURES systematic false positives for on/off proteins near the detection limit
# (the volcano "anchor arms"): it pins missing values ~1.8 SD below the mean with an artificially
# tight 0.3 SD spread, inflating the t-numerator and deflating its denominator. The honest report
# for a protein fully missing in one group is "undetected in group B", NOT a fold change. If a
# reviewer demands a complete matrix anyway, seed it (set.seed) so the call list is reproducible,
# and report the on/off proteins separately -- see proteomics/differential-abundance.

# === 5. QC ===
# PCA needs a complete matrix; use the proteins observed everywhere rather than inventing values.
# The complete-case set can be EMPTY -- 12 samples at realistic MNAR dropout is enough -- and
# prcomp then aborts the whole workflow at step 5 of 7 with `a dimension is zero`, which says
# nothing about missingness. Guard it, say what happened, and fall back to the pairwise-complete
# sample correlation, which answers the same clustering question without a complete matrix.
# QC is diagnostic: it must not gate the statistics below, and it is never a reason to impute.
complete <- filtered[stats::complete.cases(filtered), , drop = FALSE]
cat('Proteins used for PCA (complete cases):', nrow(complete), '\n')
pca_df <- NULL
if (nrow(complete) >= 3) {
    pca <- prcomp(t(complete), scale. = TRUE)
    pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Sample = rownames(pca$x))
} else {
    message('PCA skipped: only ', nrow(complete), ' protein(s) observed in EVERY sample. ',
            'That is missingness, not a corrupt matrix -- lower min_frac, drop the sparsest ',
            'samples, or read the correlation below. Do NOT impute to fill the matrix.')
    sample_cor <- cor(as.matrix(filtered), use = 'pairwise.complete.obs', method = 'spearman')
    print(round(sample_cor, 2))
}

# === 6. DIFFERENTIAL ANALYSIS (limma on the observed values) ===
# sample_info is already loaded and factored in step 4. Batch is a COVARIATE in the same model
# (rule 4); removeBatchEffect() is visualization-only, never an input to lmFit.
has_batch <- 'batch' %in% colnames(sample_info) && length(unique(sample_info$batch)) > 1
design <- if (has_batch) model.matrix(~ 0 + condition + factor(batch), data = sample_info) else
                         model.matrix(~ 0 + condition, data = sample_info)
colnames(design)[seq_along(levels(sample_info$condition))] <- levels(sample_info$condition)
colnames(design) <- make.names(colnames(design))
cat('design columns:', colnames(design), '| batch in design:', has_batch, '\n')

fit <- lmFit(as.matrix(filtered), design)
# Build the contrasts FROM the factor levels. A hard-coded `makeContrasts(Treatment - Control)`
# handles exactly one design: a three-condition dose series or a time course dies with
# `object 'Treatment' not found` even though the design matrix built one line above is correct.
# The reference is the FIRST level -- set it explicitly with relevel() if alphabetical order puts
# the wrong condition first. Two conditions give the single Treatment-Control contrast as before.
ref <- levels(sample_info$condition)[1]
others <- setdiff(levels(sample_info$condition), ref)
contrast <- makeContrasts(contrasts = paste(make.names(others), '-', make.names(ref)), levels = design)
colnames(contrast) <- paste0(others, '_vs_', ref)
cat('contrasts:', colnames(contrast), '\n')
fit2 <- contrasts.fit(fit, contrast)

# Select on FDR ALONE. A post-hoc fold-change + significance double filter inflates FDR
# (a collider/selection effect; realized FDR can exceed 50%). To require a minimum effect,
# use the moderated minimum-fold-change test treat()/topTreat() instead of filtering after.
# lfc = log2(1.5) is a per-contrast MINIMUM effect. On a dose series or time course the intermediate
# levels carry a SMALLER true effect than the extreme one, so the same floor can return zero calls
# there while the top level is significant. Measured (4 replicates per level, 120 proteins with a
# real 0.7 log2 low-dose effect): High_vs_Ctl / Low_vs_Ctl calls were 20 / 0 at lfc = log2(1.5),
# 52 / 0 at lfc = 0.3 and 78 / 8 at lfc = 0 -- the floor is not the only loss (low dose is also
# underpowered), so a zero on an intermediate level is "not detected at this n", never "no effect".
# For such designs lower lfc, and screen with the eBayes/topTable F-test below (75 hits on the same
# data) reporting the pairwise contrasts for direction only.
fit2_treat <- treat(fit2, lfc = log2(min_fc), trend = TRUE, robust = TRUE)   # moderated min-FC test; trend+robust ~mandatory for label-free LFQ
# topTreat adjusts WITHIN one contrast. With more than one contrast the family of tests is the
# whole set, so adjust ACROSS them: decideTests(method = 'global') applies one BH over every
# protein x contrast cell, which is what the significance calls below use. With a single contrast
# it is identical to the per-contrast BH, so the two-condition result is unchanged.
results <- do.call(rbind, lapply(colnames(contrast), function(cn) {
    tt <- topTreat(fit2_treat, coef = cn, number = Inf, sort.by = 'none')
    data.frame(protein = rownames(tt), contrast = cn, tt, row.names = NULL)
}))
# unclass(): decideTests returns an S4 TestResults, which as.matrix() does NOT demote, and whose
# `[` refuses the two-column character index used below ('Two subscripts required').
dt <- unclass(decideTests(fit2_treat, method = 'global', adjust.method = 'BH', p.value = 0.05))
# A protein observed in only one group (or one batch level) has no estimable contrast: limma
# returns NA, and `adj.P.Val < 0.05` would make every downstream sum() NA. Report those separately
# as "undetected in group B" -- that is the honest statement, not a fold change.
results$significant <- !is.na(results$adj.P.Val) &
                       dt[cbind(results$protein, results$contrast)] != 0
# To ask instead "does this protein change ANYWHERE across the conditions" (an omnibus/ANOVA
# question rather than a set of pairwise ones), F-test all contrast columns at once:
#   fit2_eb <- eBayes(fit2, trend = TRUE, robust = TRUE)
#   anova_res <- topTable(fit2_eb, number = Inf)   # F and its BH-adjusted p across all contrasts
# Use it to screen, then report the pairwise contrasts for direction. treat() has no F-test.

# === 7. OUTPUT ===
# `results` has one row per protein x contrast, columns: protein, contrast, logFC, AveExpr, t, P.Value,
# adj.P.Val, significant. `significant` (logical) is the call to use: it comes from the GLOBAL
# decideTests over the whole contrast family, not from adj.P.Val < 0.05 alone (adj.P.Val is
# per-contrast); NA logFC = not estimable, i.e. undetected in one group.
cat('\nResults:\n')
cat('  Contrast not estimable (report as undetected-in-group):', sum(is.na(results$logFC)), '\n')
cat('  Significant proteins:', sum(results$significant), '\n')
cat('  Up-regulated:', sum(results$significant & results$logFC > 0), '\n')
cat('  Down-regulated:', sum(results$significant & results$logFC < 0), '\n')
for (cn in colnames(contrast)) {                       # per-contrast breakdown (>=2 contrasts)
    k <- results$contrast == cn & results$significant
    cat('   ', cn, ':', sum(k), '(Up', sum(k & results$logFC > 0), 'Down', sum(k & results$logFC < 0), ')\n')
}

write.csv(results, out_file, row.names = FALSE)
cat('wrote', out_file, '\n')
