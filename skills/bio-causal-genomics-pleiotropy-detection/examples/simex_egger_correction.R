# Reference: simex 1.8+, TwoSampleMR 0.5.11+ | Verify API if version differs
## SIMEX correction for MR-Egger under NOME violation
##
## When I^2_GX (Bowden 2016 IJE 45:1961) falls below 0.9, exposure-side
## measurement error regression-dilutes the Egger slope toward the null and
## inflates the intercept. SIMEX (Cook & Stefanski 1994 JASA 89:1314) re-fits
## Egger at simulated, amplified measurement-error levels and extrapolates
## back to zero error. Below I^2_GX 0.6, Egger should be dropped entirely
## (MR-RAPS or CAUSE preferred). Between 0.6 and 0.9, SIMEX is the rescue, but
## its extrapolation is noisy in small samples and can land further from the truth
## than the naive slope: report both slopes and compare against IVW / median / mode.

library(simex)
library(TwoSampleMR)

isq <- Isq(dat$beta.exposure, dat$se.exposure)
cat('I^2_GX:', round(isq, 3), '\n')

if (isq < 0.6) {
    cat('NOME severely violated; drop Egger; use MR-RAPS or CAUSE\n')
} else if (isq < 0.9) {
    cat('NOME partially violated; applying SIMEX correction\n')

    # Precompute the weights vector rather than dividing a data-frame column
    # in-formula: simex() refits the model internally on perturbed data and
    # cannot re-evaluate `1 / se.outcome^2` against its own working frame,
    # which has no se.outcome column -- that in-formula form crashes with
    # "object 'se.outcome' not found" inside simex()'s refit (simex 1.8).
    w <- 1 / dat$se.outcome^2
    egger_lm <- lm(beta.outcome ~ beta.exposure,
                   weights = w, data = dat,
                   x = TRUE, y = TRUE)

    egger_simex <- simex(model = egger_lm,
                         SIMEXvariable = 'beta.exposure',
                         measurement.error = dat$se.exposure,
                         lambda = seq(0.5, 2, 0.5),
                         B = 1000,
                         fitting.method = 'quadratic',
                         asymptotic = FALSE)

    print(summary(egger_simex))
    cat('SIMEX-corrected slope:', round(coef(egger_simex)['beta.exposure'], 4), '\n')
    cat('Naive Egger slope:   ', round(coef(egger_lm)['beta.exposure'], 4), '\n')
} else {
    cat('NOME holds; report standard Egger slope\n')
}
