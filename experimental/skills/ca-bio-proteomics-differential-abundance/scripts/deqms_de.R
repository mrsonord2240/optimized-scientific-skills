#!/usr/bin/env Rscript
# Purpose: DEqMS differential abundance: the limma fit of limma_de.R plus a prior variance tied to each
#   protein's PSM or peptide count.
# Inputs:  same protein table, samples CSV, contrast and output as limma_de.R, plus the per-protein count:
#            count_col=<column>   raw TSV: column holding the count (e.g. Peptides)
#            count=<file.csv>     CSV with columns protein, count (use this for TMT PSM counts)
#          TMT: PSM count; label-free: peptide count; multi-batch TMT: the MINIMUM count across batches.
# Usage:   Rscript deqms_de.R proteins.tsv samples.csv Treatment-Control results.csv count_col=Peptides [key=value ...]
#          other key=value options are those of limma_de.R (normalize, min_valid, prefix, id, format).
# Output:  results CSV (protein, logFC, AveExpr, t, P.Value, adj.P.Val, count, sca.t, sca.P.Value, sca.adj.pval);
#          read the sca.* columns. <results>_dropped.csv lists proteins not tested.
# Checked: not executed. Written against DEqMS 1.24.0 (R 4.4.3 / Bioconductor 3.20).
.here <- local({
  f <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub('^--file=', '', f[1]))) else getwd()
})
source(file.path(.here, 'limma_de.R'))
suppressPackageStartupMessages(library(DEqMS))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop('Usage: Rscript deqms_de.R proteins.tsv samples.csv Treatment-Control results.csv count_col=<column>|count=<file.csv> [key=value ...]')
kv <- parse_kv(args[-(1:4)])
samples <- read_samples(args[2])
inp <- read_protein_matrix(args[1], samples, kv)
counts <- inp$counts
if (!is.null(kv$count)) {
  cf <- read.csv(kv$count, stringsAsFactors = FALSE)
  counts <- setNames(as.numeric(cf$count), cf$protein)
}
if (is.null(counts)) stop('Give count_col=<column> (raw TSV) or count=<file.csv> with columns protein, count')

out <- run_limma_de(inp$matrix, samples, args[3], if (!is.null(kv$min_valid)) as.integer(kv$min_valid) else 2L,
                    if (!is.null(kv$normalize)) kv$normalize else 'check')
fit2 <- out$fit2
stopifnot(all(fit2$df.residual > 0))  # NA-sigma rows make spectraCounteBayes recycle loess predictions
fit2$count <- counts[rownames(fit2$coefficients)]
if (anyNA(fit2$count) || any(fit2$count < 1)) stop('Every tested protein needs a count of at least 1; missing or zero for: ',
  paste(head(rownames(fit2$coefficients)[is.na(fit2$count) | fit2$count < 1], 5), collapse = ', '))
fit3 <- spectraCounteBayes(fit2)
res <- outputResult(fit3, coef_col = 1)
if (!'protein' %in% names(res)) res <- cbind(protein = rownames(res), res)
res <- res[order(res$sca.P.Value), ]
write.csv(res, args[4], row.names = FALSE)
write.csv(rbind(inp$removed, out$dropped), sub('\\.csv$', '_dropped.csv', args[4]), row.names = FALSE)
cat('median log2FC:', sprintf('%+.3f', out$offset), '| tested:', nrow(res), '\n')
cat('significant (sca.adj.pval < 0.05):', sum(res$sca.adj.pval < 0.05), '| limma adj.P.Val < 0.05:', sum(res$adj.P.Val < 0.05), '\n')
