# Adaptive clinical trial design examples.
#
# Reference: rpact 4.2+ | Verify API if version differs
# Reference: gsDesign 3.6+ | Verify API if version differs
#
# Executable examples cover group-sequential designs, blinded SSR, BOIN and
# CRM dose finding, and a historical-data MAP prior. Promising-zone SSR,
# enrichment, platform graduation, and RAR remain documented-only scaffolds.

simulation_seeds <- list(
    boin = 20260928L,
    crm = 20260928L,
    map = 20260928L
)
simulation_replicates <- list(boin = 1000L, crm = 1000L)

# ----------------------------------------------------------------------
# 1. O'Brien-Fleming group-sequential design for survival
# ----------------------------------------------------------------------
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


# ----------------------------------------------------------------------
# 4. Promising-zone SSR [DOCUMENTED_ONLY]
# ----------------------------------------------------------------------
# No executable promising-zone design is supplied here. A valid implementation
# must pre-specify the conditional-power estimator and zones, the adaptation
# function, n_max, the original stage weights used by the final combination
# test, IDMC-only output, and null/alternative calibration with Monte Carlo
# uncertainty. Do not treat an inverse-normal object alone as that design.


# ----------------------------------------------------------------------
# 5. BOIN Phase 1 dose-finding (FDA Fit-for-Purpose 2021)
# ----------------------------------------------------------------------
library(BOIN)

# Generate pre-tabulated escalation decisions
boin_table <- get.boundary(
    target = 0.30,           # target DLT rate
    ncohort = 10,            # 10 cohorts of size 3 -> max 30 patients
    cohortsize = 3,
    n.earlystop = 12,        # stop early at lowest dose if 12 patients show futility/safety
    p.saf = 0.6 * 0.30,      # boundary for "safe" (escalate)
    p.tox = 1.4 * 0.30       # boundary for "toxic" (de-escalate)
)
print(boin_table)
# This table is printed in the protocol; investigator looks up at bedside
# No real-time Bayesian software is needed. FDA issued a Fit-for-Purpose
# determination for the BOIN method in December 2021; that is not a general
# preference or approval of any proposed BOIN trial.

# Simulate operating characteristics
boin_oc <- get.oc(
    target = 0.30,
    p.true = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.55),  # true DLT rates per dose
    ncohort = 10,
    cohortsize = 3,
    ntrial = simulation_replicates$boin,
    n.earlystop = 12,
    seed = simulation_seeds$boin
)
print(boin_oc)
stopifnot(
    is.numeric(boin_oc$selpercent),
    length(boin_oc$selpercent) == 6L,
    abs(sum(boin_oc$selpercent) - 100) < 1e-8,
    is.numeric(boin_oc$totaln),
    length(boin_oc$totaln) == 1L,
    boin_oc$totaln > 0,
    boin_oc$totaln <= 30
)
# Reports correct MTD selection rate, overdose risk, average sample size


# ----------------------------------------------------------------------
# 6. Continual Reassessment Method (CRM) -- alternative to BOIN
# ----------------------------------------------------------------------
library(dfcrm)

# CRM with a calibrated six-dose logistic skeleton
true_dlt_rates <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.55)
prior_skeleton <- c(0.05, 0.10, 0.20, 0.30, 0.50, 0.65)
dose_labels <- seq_along(true_dlt_rates)
target <- 0.30

validate_crm_grid <- function(truth, skeleton, labels) {
    if (length(truth) != length(skeleton) || length(truth) != length(labels)) {
        stop("CRM truth, skeleton, and dose labels must have identical lengths")
    }
    if (length(truth) < 2L || any(!is.finite(truth)) ||
        any(!is.finite(skeleton)) || any(truth <= 0 | truth >= 1) ||
        any(skeleton <= 0 | skeleton >= 1) || any(diff(truth) <= 0) ||
        any(diff(skeleton) <= 0)) {
        stop("CRM truth and skeleton must be finite, strictly increasing probabilities")
    }
    invisible(TRUE)
}

validate_crm_grid(true_dlt_rates, prior_skeleton, dose_labels)

# Simulate CRM operating characteristics
crm_sim <- crmsim(
    PI = true_dlt_rates,
    prior = prior_skeleton,
    target = target,
    n = 30,
    x0 = 1,                  # starting dose
    nsim = simulation_replicates$crm,
    mcohort = 1,             # cohort size
    count = FALSE,
    method = 'bayes',
    model = 'logistic',
    seed = simulation_seeds$crm
)
print(crm_sim)
stopifnot(
    inherits(crm_sim, "sim"),
    identical(length(crm_sim$PI), length(prior_skeleton)),
    all(c("PI", "prior", "MTD", "level", "tox", "seed") %in% names(crm_sim)),
    identical(as.integer(crm_sim$seed), simulation_seeds$crm)
)
# Compare to BOIN OCs; CRM more efficient under correct skeleton, sensitive to skeleton mis-spec


# ----------------------------------------------------------------------
# 7. Adaptive enrichment (population selection)
# ----------------------------------------------------------------------
# rpact has limited native enrichment support; use adaptr or custom implementation

# Conceptual: at interim, evaluate conditional power in:
# - Full population
# - Biomarker-positive subgroup
# If full CP < threshold and subgroup CP > threshold, drop biomarker-negative
# Closed-testing across full and subgroup preserves FWER

# library(adaptr)
# adapt_enrich_design <- create_adaptr_design(...)


# ----------------------------------------------------------------------
# 8. Historical-data MAP prior demonstration (not EXNEX)
# ----------------------------------------------------------------------
library(RBesT)

# This fits one meta-analytic predictive prior. It does not implement
# stratum-level EX/NEX components or automatic detachment for a conflicting
# basket. See the routed guidance before implementing EXNEX or robust MAP.
historical_data <- data.frame(
    study = c('study1', 'study2', 'study3'),
    n = c(40, 35, 50),
    r = c(8, 6, 12)
)

set.seed(simulation_seeds$map)
map_fit <- gMAP(
    cbind(r, n - r) ~ 1 | study,
    data = historical_data,
    family = binomial,
    tau.dist = 'HalfNormal',
    tau.prior = 0.5,
    beta.prior = cbind(0, 2),
    iter = 2400,
    warmup = 800,
    thin = 2,
    chains = 2,
    cores = 1
)
print(map_fit)
map_prior <- automixfit(map_fit)
map_ess <- ess(map_prior)
print(map_prior)
print(map_ess)
stopifnot(
    is.matrix(map_prior),
    ncol(map_prior) > 0L,
    is.numeric(map_ess),
    length(map_ess) == 1L,
    is.finite(map_ess),
    map_ess > 0
)

execution_metadata <- list(
    seeds = simulation_seeds,
    replicates = simulation_replicates,
    package_versions = vapply(
        c("gsDesign", "rpact", "BOIN", "dfcrm", "RBesT"),
        function(package) as.character(packageVersion(package)),
        character(1)
    )
)
print(execution_metadata)


# ----------------------------------------------------------------------
# 9. I-SPY 2 style graduation criterion (conceptual)
# ----------------------------------------------------------------------
# Posterior predictive probability of success in Phase 3
# graduation = PP(Phase 3 trial of size 300 succeeds | current data) >= 0.85
# Implementation requires Stan/JAGS or FACTS (Berry Consultants commercial)

# Conceptual scaffold (requires custom Bayesian implementation):
# 1. Fit hierarchical model to current platform data
# 2. Simulate forward: draw treatment effect from posterior
# 3. For each draw, simulate Phase 3 trial of size 300
# 4. Compute proportion of simulations achieving Phase 3 success
# 5. If proportion >= 0.85, arm graduates


# ----------------------------------------------------------------------
# 10. RAR with time-trend adjustment (Robertson 2023 consensus)
# ----------------------------------------------------------------------
# RAR appropriate for multi-arm (>=3 arms) and rare disease
# Pre-specify time-trend covariate (e.g., enrollment quarter) in primary analysis
# Use proper analysis weights (e.g., inverse-probability weighting)

# Conceptual: in 4-arm trial, update allocation probabilities based on posterior
# P(superior | data) every quarter; floor at minimum allocation (e.g., 10% per arm)
# Primary analysis includes enrollment quarter as covariate in Cox or logistic
