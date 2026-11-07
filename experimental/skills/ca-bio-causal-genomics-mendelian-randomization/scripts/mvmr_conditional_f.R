# Purpose : Multivariable MR with per-exposure conditional F: format_mvmr -> strength_mvmr (guard at
#           conditional F < 1) -> ivw_mvmr -> pleiotropy_mvmr (Q_A). Not the qhet_mvmr fallback
#           (see references/mvmr-conditional-f.md for its caveats).
# Inputs  : a TSV of harmonised instruments with columns SNP, beta.<x>, se.<x> for each exposure <x>
#           (default x1 and x2) and beta.y, se.y.
# Usage   : r.sh scripts/mvmr_conditional_f.R --dat mvmr_input.tsv [--exposures x1,x2] [--gencov 0]
#           gencov = 0 is valid ONLY if the exposure GWAS samples do not overlap; for overlapping
#           exposures pass the bivariate LDSC intercept.
# Output  : prints conditional F, MVMR-IVW and Q_A; stops with an error if any conditional F < 1.

suppressMessages(library(MVMR))

args <- commandArgs(TRUE)
opt <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
dat_file <- opt('--dat'); exposures <- strsplit(opt('--exposures', 'x1,x2'), ',')[[1]]
gencov <- as.numeric(opt('--gencov', '0'))
if (is.null(dat_file)) stop('usage: --dat mvmr_input.tsv [--exposures x1,x2] [--gencov 0]')
dat <- read.delim(dat_file, stringsAsFactors = FALSE)

mvmr_dat <- format_mvmr(
    BXGs = as.matrix(dat[paste0('beta.', exposures)]),
    BYG = dat$beta.y,
    seBXGs = as.matrix(dat[paste0('se.', exposures)]),
    seBYG = dat$se.y,
    RSID = dat$SNP
)

condF <- strength_mvmr(r_input = mvmr_dat, gencov = gencov)  # per-exposure conditional F
# condF must be > 10 for EACH exposure (Sanderson 2019); total F is misleading
print(condF)
if (any(condF < 1)) {
    stop("Conditional F < 1 for at least one exposure -- qhet_mvmr's own estimate is unreliable ",
         "at this floor (can flip an exposure's sign; see references/mvmr-conditional-f.md). Report MVMR-IVW ",
         "with a weak-instrument caveat instead, or acquire stronger/less-correlated instruments.")
}

mv_ivw <- ivw_mvmr(r_input = mvmr_dat)
mv_qa <- pleiotropy_mvmr(r_input = mvmr_dat, gencov = gencov)  # Q_A heterogeneity test
print(mv_ivw); print(mv_qa)
