# Purpose : SIMEX correction of MR-Egger for NOME violation (I^2_GX < 0.9). Prints the naive and the
#           SIMEX-corrected Egger slope.
# Inputs  : harmonised.tsv from twosample_workflow.R (columns beta.exposure, se.exposure,
#           beta.outcome, se.outcome).
# Usage   : r.sh scripts/simex_egger.R --dat mr_out/harmonised.tsv [--B 1000] [--seed 1]
#           --seed is optional (simex is Monte-Carlo; set it for a reproducible corrected slope).

suppressMessages(library(simex))

args <- commandArgs(TRUE)
opt <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
dat_file <- opt('--dat'); B <- as.integer(opt('--B', '1000')); seed <- opt('--seed')
if (is.null(dat_file)) stop('usage: --dat harmonised.tsv [--B 1000] [--seed S]')
if (!is.null(seed)) set.seed(as.integer(seed))
dat <- read.delim(dat_file, stringsAsFactors = FALSE)
if ('mr_keep' %in% names(dat)) dat <- dat[dat$mr_keep, ]

# Precompute the weights vector rather than dividing a data-frame column in-formula:
# simex() refits the model internally on perturbed data and cannot re-evaluate
# `1 / se.outcome^2` against its own working frame, which has no se.outcome column --
# that in-formula form crashes with "object 'se.outcome' not found" inside simex()'s
# refit (simex 1.8). Precomputing the vector and passing a fully-qualified `data = dat`
# avoids it.
w <- 1 / dat$se.outcome^2
egger_lm <- lm(beta.outcome ~ beta.exposure, weights = w, data = dat, x = TRUE, y = TRUE)

egger_simex <- simex(model = egger_lm, SIMEXvariable = 'beta.exposure',
                      measurement.error = dat$se.exposure,
                      lambda = seq(0.5, 2, 0.5), B = B,
                      fitting.method = 'quadratic', asymptotic = FALSE)

cat('SIMEX-corrected slope:', round(coef(egger_simex)['beta.exposure'], 4), '\n')
cat('Naive Egger slope:   ', round(coef(egger_lm)['beta.exposure'], 4), '\n')
