# GPTomics source examples; see LICENSE.
simulation_seeds <- list(
    boin = 20260928L,
    crm = 20260928L,
    map = 20260928L
)
simulation_replicates <- list(boin = 1000L, crm = 1000L)

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
