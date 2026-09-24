# ihw_safe.R -- covariate-weighted FDR (IHW) that survives the IHW solver crash.
# Purpose: ihw() can SEGFAULT the R process (a crash tryCatch cannot catch; 50-75% of runs at
#   m=18,000). ihw_safe() runs ihw() in a child Rscript process, retries, and falls back to plain BH.
# Inputs:  p (numeric p-values), covariate (independent of null p, same length), alpha, nbins, tries.
# Output:  list(padj = adjusted p-values, method = 'IHW' or the BH fallback label, attempts = n tries).
# Usage (function):  source('scripts/ihw_safe.R'); res <- ihw_safe(de_table$pvalue, de_table$mean_expression)
# Usage (CLI):       Rscript scripts/ihw_safe.R <in.csv> <p_col> <covariate_col> <out.csv> [alpha=0.05]
# Checked on IHW 1.34.0 / R 4.4.3.

# ihw_safe(): IHW in a child process (a segfault kills the child, not your session).
ihw_safe <- function(p, covariate, alpha = 0.05, nbins = 5, tries = 3) {
  inp <- tempfile(fileext = '.rds'); out <- tempfile(fileext = '.rds')
  on.exit(unlink(c(inp, out)))
  saveRDS(list(p = p, cov = covariate, alpha = alpha, nbins = nbins), inp)
  code <- sprintf("d <- readRDS('%s'); library(IHW); r <- ihw(d$p, d$cov, alpha = d$alpha, nbins = d$nbins); saveRDS(adj_pvalues(r), '%s')",
                  normalizePath(inp, winslash = '/'), normalizePath(out, winslash = '/', mustWork = FALSE))
  for (i in seq_len(tries)) {
    st <- system2(file.path(R.home('bin'), 'Rscript'), c('-e', shQuote(code)), stdout = FALSE, stderr = FALSE)
    if (identical(st, 0L) && file.exists(out)) return(list(padj = readRDS(out), method = 'IHW', attempts = i))
  }
  list(padj = p.adjust(p, 'BH'), method = 'BH (IHW solver crashed; fallback)', attempts = tries)
}

# CLI entry: only when run as a script with arguments, not when source()d.
if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) >= 4) {
    d <- read.csv(a[1])
    alpha <- if (length(a) >= 5) as.numeric(a[5]) else 0.05
    res <- ihw_safe(d[[a[2]]], d[[a[3]]], alpha = alpha)
    d$padj_ihw <- res$padj
    write.csv(d, a[4], row.names = FALSE)
    cat(sprintf('method: %s (attempts %d); discoveries at %.2f: %d\n',
                res$method, res$attempts, alpha, sum(res$padj < alpha)))
  }
}
