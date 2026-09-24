# Purpose : High-dimensional mediation (HIMA2, DBlasso) of an exposure -> outcome effect through
#           thousands of candidate mediators (CpGs, genes), with BH-FDR control.
# Inputs  : pheno.csv     one row per sample; columns named in the formula (outcome, exposure, covariates)
#           mediators.csv one row per sample (same order as pheno.csv), one column per mediator
#           formula       e.g. "outcome ~ exposure + age + sex + cell_pc1 + cell_pc2"
#                         (LHS Surv(time, status) routes to Cox). Character/factor covariates are
#                         dummy-coded with model.matrix() automatically (HIMA 2.3.4 rejects factors).
# Output  : <out.csv> significant mediators (ID, alpha, beta, alpha*beta, rimp, p-value); prints count
# Usage   : Rscript scripts/hima_ewas.R pheno.csv mediators.csv "outcome ~ exposure + age + sex" out.csv \
#                  [mediator.type=gaussian] [penalty=DBlasso] [sigcut=0.05] [ncore=1]
#           (run R through your env wrapper, e.g. r.sh; checked on HIMA 2.3.4, R 4.4.3)
suppressPackageStartupMessages(library(HIMA))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 4) stop('usage: hima_ewas.R pheno.csv mediators.csv "formula" out.csv [mediator.type] [penalty] [sigcut] [ncore]')
pheno_csv <- a[1]; med_csv <- a[2]; fml <- a[3]; out_csv <- a[4]
mediator.type <- if (length(a) >= 5) a[5] else 'gaussian'   # 'negbin' for count, 'compositional' for microbiome
penalty       <- if (length(a) >= 6) a[6] else 'DBlasso'    # alternatives 'MCP', 'SCAD', 'lasso'
sigcut        <- if (length(a) >= 7) as.numeric(a[7]) else 0.05
ncore         <- if (length(a) >= 8) as.integer(a[8]) else 1L

dat <- read.csv(pheno_csv, stringsAsFactors = FALSE, check.names = FALSE)
M_matrix <- as.matrix(read.csv(med_csv, check.names = FALSE))
stopifnot(nrow(dat) == nrow(M_matrix))

# Drop NA rows for the variables in the formula; keep M_matrix rows aligned with dat
vars <- intersect(all.vars(as.formula(fml)), names(dat))
keep <- complete.cases(dat[, vars, drop = FALSE])
dat <- dat[keep, , drop = FALSE]
M_matrix <- M_matrix[keep, , drop = FALSE]

# Dummy-code categorical RHS covariates (model.matrix, drop intercept) and rewrite the formula terms
f <- as.formula(fml)
lhs <- paste(deparse(f[[2]]), collapse = '')
terms_rhs <- attr(terms(f), 'term.labels')
new_terms <- character(0)
for (tm in terms_rhs) {
  if (tm %in% names(dat) && (is.character(dat[[tm]]) || is.factor(dat[[tm]]))) {
    d <- model.matrix(~ f, data = data.frame(f = factor(dat[[tm]])))[, -1, drop = FALSE]
    colnames(d) <- paste0(tm, sub('^f', '', colnames(d)))
    dat <- cbind(dat[, setdiff(names(dat), tm), drop = FALSE], d)
    new_terms <- c(new_terms, colnames(d))
  } else new_terms <- c(new_terms, tm)
}
fml <- paste(lhs, '~', paste(new_terms, collapse = ' + '))

result <- hima(
  as.formula(fml),
  data.pheno = dat,
  data.M = M_matrix,
  mediator.type = mediator.type,
  penalty = penalty,
  scale = TRUE,
  sigcut = sigcut,
  parallel = ncore > 1, ncore = ncore, verbose = FALSE
)
# result is a LIST of class "hima" ($ID, $alpha, $beta, $`alpha*beta`, $rimp, $`p-value`),
# NOT a data.frame: nrow(result)/rownames(result) return NULL. Use result$ID / length(result$ID).
sig_mediators <- result$ID
n_sig <- length(sig_mediators)
cat('Significant mediators (FDR <', sigcut, '):', n_sig, '\n')
if (n_sig > 0) {
  out <- data.frame(ID = result$ID, alpha = result$alpha, beta = result$beta,
                    alpha_beta = result$`alpha*beta`, rimp = result$rimp, p = result$`p-value`)
  write.csv(out, out_csv, row.names = FALSE)
  cat('Wrote', out_csv, '\n')
} else {
  write.csv(data.frame(ID = character(0)), out_csv, row.names = FALSE)
}
