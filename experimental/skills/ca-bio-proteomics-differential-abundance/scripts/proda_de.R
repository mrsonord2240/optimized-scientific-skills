#!/usr/bin/env Rscript
# Purpose: proDA differential abundance with heavy intensity-dependent missingness: missing values enter the
#   likelihood as left-censored observations; nothing is imputed.
# Inputs:  same protein table, samples CSV (sample, condition, optionally batch), contrast and output as limma_de.R.
# Usage:   Rscript proda_de.R proteins.tsv samples.csv Treatment-Control results.csv [key=value ...]
#            normalize=check|median|none (as limma_de.R), prefix, id, format
# Output:  results CSV (protein, pval, adj_pval, diff, t_statistic, se, n_<group>, undetected_in); diff,
#          t_statistic and se are blank for proteins with no observed value in a group.
#          <results>_dropped.csv lists flagged rows and rows with no observed value.
# Checked: not executed. Written against proDA 1.20.0 (R 4.4.3 / Bioconductor 3.20).
.here <- local({
  f <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub('^--file=', '', f[1]))) else getwd()
})
source(file.path(.here, 'protein_input.R'))
suppressPackageStartupMessages(library(proDA))

OFFSET_MAX <- 0.1  # |median log2FC| trip-wire, as in limma_de.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop('Usage: Rscript proda_de.R proteins.tsv samples.csv Treatment-Control results.csv [key=value ...]')
kv <- parse_kv(args[-(1:4)])
normalize <- if (!is.null(kv$normalize)) kv$normalize else 'check'
samples <- read_samples(args[2])
inp <- read_protein_matrix(args[1], samples, kv)
parts <- strsplit(args[3], '-', fixed = TRUE)[[1]]
if (length(parts) != 2 || !all(parts %in% samples$condition)) stop('Contrast must be Test-Reference with both names in the condition column')
test_level <- parts[1]; ref_level <- parts[2]

m <- inp$matrix
if (normalize == 'median') m <- center_medians(m)
none <- rowSums(!is.na(m)) == 0
dropped <- rbind(inp$removed, data.frame(protein = rownames(m)[none], reason = rep('no observed value', sum(none))))
m <- m[!none, , drop = FALSE]

col_data <- samples
rownames(col_data) <- samples$sample
col_data$condition <- factor(col_data$condition)
form <- if ('batch' %in% names(col_data)) { col_data$batch <- factor(col_data$batch); ~condition + batch } else ~condition
fit <- proDA(m, design = form, col_data = col_data, reference_level = ref_level)
res <- test_diff(fit, paste0('condition', test_level))

n_obs <- sapply(c(test_level, ref_level), function(g) rowSums(!is.na(m[, col_data$condition == g, drop = FALSE])))
colnames(n_obs) <- paste0('n_', c(test_level, ref_level))
res <- cbind(protein = res$name, res[, c('pval', 'adj_pval', 'diff', 't_statistic', 'se')], n_obs[match(res$name, rownames(n_obs)), , drop = FALSE])
res$undetected_in <- ifelse(res[[paste0('n_', test_level)]] == 0, test_level, ifelse(res[[paste0('n_', ref_level)]] == 0, ref_level, ''))
both <- res$undetected_in == ''
offset <- median(res$diff[both], na.rm = TRUE)
if (normalize == 'check' && abs(offset) > OFFSET_MAX) stop(sprintf(paste(
  'median log2FC = %+.3f: the matrix does not look normalized. Rerun with normalize=median, or normalize=none',
  'if a global shift is the expected biology.'), offset))
res[!both, c('diff', 't_statistic', 'se')] <- NA  # the location prior sets these; the sign can be wrong
res <- res[order(res$pval), ]
write.csv(res, args[4], row.names = FALSE)
write.csv(dropped, sub('\\.csv$', '_dropped.csv', args[4]), row.names = FALSE)
cat('median log2FC (proteins seen in both groups):', sprintf('%+.3f', offset), '| tested:', nrow(res),
    '| undetected in one group:', sum(!both), '\n')
cat('significant (adj_pval < 0.05):', sum(res$adj_pval < 0.05, na.rm = TRUE), '\n')
