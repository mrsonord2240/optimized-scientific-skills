# GPTomics source examples; see LICENSE.
simulation_seeds <- list(
    boin = 20260928L,
    crm = 20260928L,
    map = 20260928L
)
simulation_replicates <- list(boin = 1000L, crm = 1000L)

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
