# qvalue_safe.R -- run qvalue in an isolated R process and preserve a completed result.
# Purpose: qvalue can terminate its R process during teardown in this Windows R environment.
# Inputs:  p (numeric p-values) and alpha (discovery threshold).
# Output:  list(pi0, discoveries, method, worker_status); falls back to BH only when no qvalue output exists.
# Usage:   source('scripts/qvalue_safe.R'); res <- qvalue_safe(pvalues)
# Checked on qvalue 2.38.0 / R 4.4.3.

qvalue_safe <- function(p, alpha = 0.05) {
  inp <- tempfile(fileext = '.rds'); out <- tempfile(fileext = '.rds')
  on.exit(unlink(c(inp, out)))
  saveRDS(list(p = p, alpha = alpha), inp)
  code <- sprintf("d <- readRDS('%s'); library(qvalue); q <- qvalue(d$p); saveRDS(list(pi0 = q$pi0, discoveries = sum(q$qvalues < d$alpha)), '%s')",
                  normalizePath(inp, winslash = '/'), normalizePath(out, winslash = '/', mustWork = FALSE))
  status <- system2(file.path(R.home('bin'), 'Rscript'), c('-e', shQuote(code)), stdout = FALSE, stderr = FALSE)
  if (file.exists(out)) {
    result <- tryCatch(readRDS(out), error = function(e) NULL)
    if (is.list(result) && is.finite(result$pi0) && is.finite(result$discoveries)) {
      return(c(result, list(method = 'qvalue', worker_status = status)))
    }
  }
  list(pi0 = NA_real_, discoveries = sum(p.adjust(p, 'BH') < alpha),
       method = 'BH fallback (qvalue worker did not produce output)', worker_status = status)
}
