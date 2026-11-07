# GPTomics source examples; see LICENSE.
library(gsDesign)

design_obf <- gsDesign(
    k = 3,                  # 2 interim + 1 final
    test.type = 1,          # 1-sided efficacy
    alpha = 0.025,
    beta = 0.10,            # power = 0.90
    sfu = sfLDOF,           # Lan-DeMets approximation of OBF
    timing = c(0.33, 0.67, 1.0)
)
print(design_obf)
# Boundaries: very conservative early (~0.0001 nominal at 33% info),
# near-nominal at final (~0.0234 of 0.025)

plot(design_obf)


# ----------------------------------------------------------------------
# 2. Sample size for survival group-sequential
# ----------------------------------------------------------------------
n_gs <- gsSurv(
    k = 3,
    test.type = 2,           # 2-sided
    alpha = 0.025,           # one-sided alpha for each side
    beta = 0.10,             # 90% power
    sfu = sfLDOF,
    lambdaC = 0.04,          # control hazard per month
    hr = 0.65,               # treatment HR
    eta = 0.005,             # dropout hazard per month
    T = 30,                  # total study duration in months
    minfup = 18,             # minimum follow-up
    ratio = 1                # 1:1 allocation
)
print(n_gs)


# ----------------------------------------------------------------------
# 3. Blinded SSR for continuous endpoint (Friede-Kieser 2006)
# ----------------------------------------------------------------------
library(rpact)

# Pre-specify design with uncertain variance
design_initial <- getDesignGroupSequential(
    kMax = 1,                # fixed final-analysis design
    alpha = 0.025,
    beta = 0.20,             # 80% power
    sided = 1
)

# Initial sample size assuming SD = 12
ss_initial <- getSampleSizeMeans(
    design = design_initial,
    alternative = 5,         # detect mean diff of 5
    stDev = 12,
    groups = 2
)
print(ss_initial)

# At internal pilot, re-estimate SD from blinded data
# (manual: pool all observations, compute SD ignoring arm)
# If observed SD = 14 (higher than assumed), recompute n
ss_recomputed <- getSampleSizeMeans(
    design = design_initial,
    alternative = 5,
    stDev = 14,
    groups = 2
)
print(ss_recomputed)
stopifnot(
    is.finite(ss_initial$maxNumberOfSubjects),
    is.finite(ss_recomputed$maxNumberOfSubjects),
    ss_recomputed$maxNumberOfSubjects > ss_initial$maxNumberOfSubjects
)
# A prespecified blinded nuisance-only recalculation can preserve Type-I error
# when it contains no treatment-effect information and the valid final test is
# unchanged. This example does not establish that property for another rule.
