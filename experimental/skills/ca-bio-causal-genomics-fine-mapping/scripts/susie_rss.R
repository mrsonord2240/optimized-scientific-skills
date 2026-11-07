#!/usr/bin/env Rscript
# Inputs: sumstats TSV with SNP, BETA, SE; signed named LD TSV with SNP names as first column and header.
# Usage: Rscript susie_rss.R sumstats.tsv ld.tsv N output_prefix [L=10]
suppressPackageStartupMessages(library(susieR))
args <- commandArgs(TRUE)
if (length(args) < 4 || length(args) > 5) stop('usage: susie_rss.R sumstats.tsv ld.tsv N output_prefix [L]')
x <- read.delim(args[1], check.names=FALSE)
R <- as.matrix(read.delim(args[2], row.names=1, check.names=FALSE))
n <- as.numeric(args[3]); out <- args[4]; L <- if (length(args)==5) as.integer(args[5]) else 10L
stopifnot(all(c('SNP','BETA','SE') %in% names(x)), identical(rownames(R), x$SNP),
          identical(colnames(R), x$SNP), is.finite(n), n > 0, L > 0)
z <- x$BETA/x$SE
s <- as.numeric(estimate_s_rss(z=z, R=R, n=n))
if (s > 0.10) stop('LD mismatch: lambda > 0.10')
k <- kriging_rss(z=z, R=R, n=n)
fit <- susie_rss(z=z, R=R, n=n, L=L, estimate_residual_variance=TRUE)
x$PIP <- fit$pip
write.table(x, paste0(out, '_pips.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
saveRDS(fit, paste0(out, '_fit.rds'))
cat('lambda:', s, 'kriging outliers:', sum(abs(k$conditional_dist$z_std_diff)>3, na.rm=TRUE), 'converged:', fit$converged, '\n')
