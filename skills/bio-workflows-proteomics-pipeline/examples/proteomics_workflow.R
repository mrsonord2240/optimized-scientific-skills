# Reference: limma 3.62.2, ggplot2 4.0.3, pheatmap 1.0.13 (checked 2026-09-15, R 4.4.3) | Verify API if version differs
# Complete proteomics workflow: MaxQuant to differential proteins. Runs as shipped: with no
# proteinGroups.txt in the working directory it simulates a synthetic one (see below).
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# === CONFIGURATION ===
input_file <- 'proteinGroups.txt'
output_prefix <- 'proteomics_results'
fdr_threshold <- 0.05
lfc_threshold <- 1

# Sample groups (modify for your experiment)
sample_groups <- c(rep('Control', 4), rep('Treatment', 4))

# Self-contained demo: with no proteinGroups.txt present, simulate one so the script runs as shipped.
# The last sample is given a 3x-low load so the raw-distribution QC below has something to catch.
if (!file.exists(input_file)) {
    cat('No', input_file, 'found -- simulating a SYNTHETIC one for the demo.\n')
    set.seed(1); n <- 1500; s <- paste0('S', seq_along(sample_groups))
    base <- rnorm(n, 25, 1.2)
    eff <- c(rep(0, n - 100), rep(c(-2, 2), each = 50))   # 100 true changers, 1400 true nulls
    raw <- sapply(seq_along(sample_groups), function(j)
        base + (if (sample_groups[j] == 'Treatment') eff else 0) + rnorm(n, 0, 0.3) +
            c(0, 0.1, -0.1, 0.05, -0.05, 0.08, -0.08, -1.585)[j])
    raw[raw < quantile(raw, 0.08)] <- NA                    # left-censored dropout
    colnames(raw) <- s
    pg <- data.frame(Majority.protein.IDs = sprintf('P%04d', seq_len(n)),
                     Potential.contaminant = '', Reverse = '', Only.identified.by.site = '')
    pg[paste0('Intensity.', s)] <- round(ifelse(is.na(raw), 0, 2^raw))
    lfq <- sweep(raw, 2, apply(raw, 2, median, na.rm = TRUE) - median(raw, na.rm = TRUE))
    pg[paste0('LFQ.intensity.', s)] <- round(ifelse(is.na(lfq), 0, 2^lfq))
    write.table(pg, input_file, sep = '\t', quote = FALSE, row.names = FALSE)
}

# === 1. DATA IMPORT ===
cat('=== Data Import ===\n')
# quote = '' / comment.char = '': MaxQuant text fields contain apostrophes and '#', and the
# read.delim defaults silently truncate the table with only an 'EOF within quoted string' warning.
proteins <- read.delim(input_file, stringsAsFactors = FALSE, quote = '', comment.char = '')
stopifnot(nrow(proteins) == length(readLines(input_file)) - 1)
cat('Loaded', nrow(proteins), 'protein groups\n')

# %in%, NOT `!= '+'`: an all-empty flag column is typed logical NA, `NA != '+'` is NA, and
# subsetting by NA silently replaces every row with NAs while nrow() still looks right.
proteins <- proteins[!(proteins$Potential.contaminant %in% '+') &
                      !(proteins$Reverse %in% '+') &
                      !(proteins$Only.identified.by.site %in% '+'), ]
stopifnot(!all(is.na(proteins$Majority.protein.IDs)))
cat('After filtering:', nrow(proteins), 'proteins\n')

lfq_cols <- grep('^LFQ\\.intensity\\.', colnames(proteins), value = TRUE)
intensities <- proteins[, lfq_cols]
rownames(intensities) <- proteins$Majority.protein.IDs
colnames(intensities) <- gsub('LFQ\\.intensity\\.', '', colnames(intensities))
cat('Samples:', ncol(intensities), '\n')

# === 2. LOG2, THEN INSPECT RAW DISTRIBUTIONS (before any normalization) ===
cat('\n=== Raw distribution QC ===\n')
intensities[intensities == 0] <- NA
log2_int <- log2(intensities)

# Median-centering rescales every sample onto a common median, so it mathematically erases the 3x-low
# load that marks a failed injection. A failed sample must be found and dropped HERE, not afterwards.
names(sample_groups) <- colnames(log2_int)   # key by column so a drop cannot misalign the design below
# Inspect the RAW `Intensity.` columns, NOT `LFQ intensity`: MaxLFQ has already renormalized the LFQ
# columns onto a common scale, so a low load is largely gone from them before you ever look.
raw_cols <- grep('^Intensity\\.', colnames(proteins), value = TRUE)
raw_log2 <- log2(replace(proteins[, raw_cols], proteins[, raw_cols] == 0, NA))
colnames(raw_log2) <- gsub('^Intensity\\.', '', colnames(raw_log2))
raw_log2 <- raw_log2[, colnames(log2_int), drop = FALSE]
raw_median <- apply(raw_log2, 2, median, na.rm = TRUE)
id_counts <- colSums(!is.na(log2_int))
print(data.frame(id_count = id_counts, raw_median_log2 = round(raw_median, 2),
                 load_shift_log2 = round(raw_median - median(raw_median), 2)))
pdf(paste0(output_prefix, '_raw_boxplot.pdf'), width = 8, height = 5)
boxplot(raw_log2, las = 2, main = 'RAW log2 Intensity (pre-normalization)', ylab = 'log2 intensity')
dev.off()

# Flag on BOTH axes and as an OUTLIER: a >=2x (1 log2) drop in raw median signal, or an ID count
# more than 3 MADs below the cohort median. A 50%-of-median ID rule alone sees almost nothing.
mad_ids <- mad(id_counts)
failed <- colnames(raw_log2)[(raw_median - median(raw_median)) <= -1 |
                             (mad_ids > 0 & id_counts < median(id_counts) - 3 * mad_ids)]
if (length(failed) > 0) {
    cat('Dropping failed samples:', paste(failed, collapse = ', '), '\n')
    log2_int <- log2_int[, !colnames(log2_int) %in% failed, drop = FALSE]
    sample_groups <- sample_groups[colnames(log2_int)]   # keep annotation aligned to surviving columns
}

# === 3. NORMALIZE (only after the raw inspection above) ===
cat('\n=== Normalization ===\n')
sample_medians <- apply(log2_int, 2, median, na.rm = TRUE)
cat('Sample medians before:', round(sample_medians, 2), '\n')
normalized <- sweep(log2_int, 2, sample_medians - median(sample_medians))
cat('Sample medians after:', round(apply(normalized, 2, median, na.rm = TRUE), 2), '\n')

# === 4. FILTER ON PER-GROUP COMPLETENESS (nothing is imputed) ===
# Per-group completeness: keep a protein valid in >= 60% of replicates in AT LEAST ONE condition,
# not a blanket across-cohort rule. Downshift imputation is deliberately NOT used -- it manufactures
# false positives for on/off proteins (the volcano-wing artifact); limma fits each protein on its
# observed values instead. Model the dropout explicitly with proDA when that matters (SKILL.md).
cat('\n=== Filtering ===\n')
group_complete <- sapply(unique(sample_groups), function(g)
    rowSums(!is.na(normalized[, names(sample_groups)[sample_groups == g], drop = FALSE])) >=
        ceiling(sum(sample_groups == g) * 0.6))
filtered <- normalized[rowSums(group_complete) > 0, ]
cat('Proteins after per-group completeness filter:', nrow(filtered), '\n')

# === 5. QC ===
cat('\n=== Quality Control ===\n')
# PCA needs a complete matrix; use the proteins observed everywhere rather than inventing values.
# The complete-case set can be EMPTY on a larger cohort with realistic MNAR dropout, and prcomp
# then aborts the whole script with `a dimension is zero`. Guard it and fall back to the
# pairwise-complete sample correlation; QC is diagnostic and must not gate the statistics below.
complete <- filtered[complete.cases(filtered), , drop = FALSE]
cat('Proteins used for PCA (complete cases):', nrow(complete), '\n')
if (nrow(complete) >= 3) {
    pca <- prcomp(t(complete), scale. = TRUE)
    pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Sample = rownames(pca$x), Group = sample_groups)
    var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

    p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
        geom_point(size = 4) + theme_minimal() +
        labs(x = paste0('PC1 (', var_exp[1], '%)'), y = paste0('PC2 (', var_exp[2], '%)'), title = 'PCA of Protein Abundances')
    ggsave(paste0(output_prefix, '_pca.pdf'), p_pca, width = 7, height = 6)
} else {
    cat('PCA skipped: only', nrow(complete), 'protein(s) observed in EVERY sample -- that is\n',
        'missingness, not a corrupt matrix. Sample correlation (pairwise complete) instead:\n')
    print(round(cor(as.matrix(filtered), use = 'pairwise.complete.obs', method = 'spearman'), 2))
}

# === 6. DIFFERENTIAL ANALYSIS ===
cat('\n=== Differential Analysis ===\n')
sample_info <- data.frame(sample = colnames(filtered), condition = factor(sample_groups, levels = c('Control', 'Treatment')))
design <- model.matrix(~ 0 + condition, data = sample_info)
colnames(design) <- levels(sample_info$condition)

fit <- lmFit(as.matrix(filtered), design)
# This demo is deliberately TWO conditions, so the contrast is written out. For three or more
# (dose series, time course) do NOT edit sample_groups and keep this line: it will fail with
# `object 'Treatment' not found`. Use the Complete R Workflow block in SKILL.md, which builds the
# contrasts from the condition levels and adjusts across them with decideTests(method='global').
contrast <- makeContrasts(Treatment - Control, levels = design)
# treat() folds the minimum fold-change INTO the test (a proper hypothesis against |lfc| > threshold),
# instead of a post-hoc logFC AND adj.P double filter -- the double filter is a collider/selection
# effect whose realized FDR can exceed the nominal rate (SKILL.md). Significance is then adj.P alone.
fit2 <- treat(contrasts.fit(fit, contrast), lfc = lfc_threshold, trend = TRUE, robust = TRUE)

results <- topTreat(fit2, coef = 1, number = Inf, adjust.method = 'BH')
results$protein <- rownames(results)
# A protein observed in only one group has no estimable contrast: limma returns NA, and a bare
# `adj.P.Val < t` would poison every downstream sum(). Report those as undetected-in-group.
results$significant <- !is.na(results$adj.P.Val) & results$adj.P.Val < fdr_threshold

cat('Total proteins tested:', nrow(results), '\n')
cat('Contrast not estimable (undetected in one group):', sum(is.na(results$logFC)), '\n')
cat('Significant:', sum(results$significant), '\n')
cat('  Up-regulated:', sum(results$significant & results$logFC > 0), '\n')
cat('  Down-regulated:', sum(results$significant & results$logFC < 0), '\n')

# === 7. VISUALIZATION ===
p_volcano <- ggplot(results, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(color = significant), alpha = 0.6) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = 'dashed') +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = 'dashed') +
    scale_color_manual(values = c('grey60', 'firebrick')) +
    theme_minimal() + labs(title = 'Volcano Plot', x = 'Log2 Fold Change', y = '-Log10 Adjusted P-value')
ggsave(paste0(output_prefix, '_volcano.pdf'), p_volcano, width = 8, height = 6)

# Heatmap of significant proteins
if (sum(results$significant) > 1) {
    sig_proteins <- rownames(results)[results$significant]
    mat <- as.matrix(filtered[sig_proteins, ])
    mat <- mat[complete.cases(mat), , drop = FALSE]   # pheatmap cannot cluster rows with NA
    mat_scaled <- t(scale(t(mat)))
    annotation_col <- data.frame(Group = sample_groups, row.names = colnames(mat))
    pheatmap(mat_scaled, annotation_col = annotation_col, show_rownames = nrow(mat_scaled) < 50,
             filename = paste0(output_prefix, '_heatmap.pdf'), width = 8, height = 10)
}

# === 8. EXPORT ===
write.csv(results, paste0(output_prefix, '.csv'), row.names = FALSE)
cat('\n=== Output Files ===\n')
cat(paste0(output_prefix, '.csv\n'))
cat(paste0(output_prefix, '_raw_boxplot.pdf\n'))
cat(paste0(output_prefix, '_pca.pdf\n'))
cat(paste0(output_prefix, '_volcano.pdf\n'))
cat(paste0(output_prefix, '_heatmap.pdf\n'))
