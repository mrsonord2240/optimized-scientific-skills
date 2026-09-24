# check_balance.R -- verify a designit layout: confounding (hard fail) and avoidable imbalance (warning).
# Purpose: optimize_design() can converge to an unbalanced local optimum; always check the
#          covariate x batch table, not just the achieved score.
# Usage:   source('scripts/check_balance.R')
#          check_balance(assignment, batch_col = 'batch', vars = c('condition', 'sex'))
#   assignment : data.frame from bc$get_samples()
#   batch_col  : name of the batch/plex/lane column
#   vars       : every balanced covariate (not only the primary condition)
check_balance <- function(assignment, batch_col, vars) {
  for (v in vars) {
    tab <- table(assignment[[v]], assignment[[batch_col]])
    print(tab)
    # Hard fail: a whole level missing from a batch means `v` and batch are confounded in this
    # layout, not merely unbalanced -- the design must not be used as-is.
    if (!all(tab > 0))
      stop(sprintf('%s is confounded with %s in this layout (a batch has zero samples of some level)',
                   v, batch_col))
    # Soft check: flag an avoidable imbalance so the agent iterates instead of shipping a
    # sub-optimal split (7/7/7/9 was accepted when 7/8/7/8 was achievable).
    imbalance <- max(tab) - min(tab)
    if (imbalance > 1)
      warning(sprintf('%s x %s: largest cell minus smallest cell = %d; re-run optimize_design() with a higher max_iter or a different random seed before accepting this layout.',
                      v, batch_col, imbalance))
  }
  invisible(TRUE)
}
