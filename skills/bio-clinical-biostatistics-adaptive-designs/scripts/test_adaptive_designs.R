# Focused regression harness for scripts/adaptive_designs.R.
# Usage: Rscript scripts/test_adaptive_designs.R [path/to/adaptive_designs.R]

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) > 0L) args[[1L]] else "scripts/adaptive_designs.R"
if (!file.exists(target)) {
    stop("adaptive-design example script not found: ", target)
}

example <- new.env(parent = globalenv())
source(target, local = example, echo = FALSE)

boin_repeat <- BOIN::get.oc(
    target = 0.30,
    p.true = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.55),
    ncohort = 10,
    cohortsize = 3,
    ntrial = example$simulation_replicates$boin,
    n.earlystop = 12,
    seed = example$simulation_seeds$boin
)

crm_repeat <- dfcrm::crmsim(
    PI = example$true_dlt_rates,
    prior = example$prior_skeleton,
    target = example$target,
    n = 30,
    x0 = 1,
    nsim = example$simulation_replicates$crm,
    mcohort = 1,
    count = FALSE,
    method = "bayes",
    model = "logistic",
    seed = example$simulation_seeds$crm
)

shape_error <- tryCatch(
    {
        example$validate_crm_grid(
            example$true_dlt_rates,
            example$prior_skeleton[-1L],
            example$dose_labels
        )
        NULL
    },
    error = conditionMessage
)

stopifnot(
    example$ss_recomputed$maxNumberOfSubjects >
        example$ss_initial$maxNumberOfSubjects,
    identical(boin_repeat$selpercent, example$boin_oc$selpercent),
    identical(boin_repeat$totaln, example$boin_oc$totaln),
    length(example$boin_oc$selpercent) == 6L,
    abs(sum(example$boin_oc$selpercent) - 100) < 1e-8,
    inherits(example$crm_sim, "sim"),
    identical(crm_repeat, example$crm_sim),
    length(example$crm_sim$PI) == length(example$prior_skeleton),
    grepl("identical lengths", shape_error, fixed = TRUE),
    is.finite(example$map_ess),
    example$map_ess > 0,
    !exists("design_promising", envir = example, inherits = FALSE),
    identical(example$simulation_seeds$boin, 20260928L),
    identical(example$simulation_seeds$crm, 20260928L),
    identical(example$simulation_seeds$map, 20260928L)
)

cat("adaptive_designs regression: PASS\n")
