library(phytools)
library(corHMM)

args <- commandArgs(trailingOnly = TRUE)
seed_arg <- grep("^--seed=", args, value = TRUE)
if (length(seed_arg) != 1L) {
    stop("Pass exactly one reproducibility seed as --seed=<integer>.", call. = FALSE)
}
seed_text <- sub("^--seed=", "", seed_arg)
if (!grepl("^[0-9]+$", seed_text)) {
    stop("The --seed value must be a non-negative integer.", call. = FALSE)
}
seed <- as.integer(seed_text)
if (is.na(seed)) stop("The --seed value is outside R's integer range.", call. = FALSE)
set.seed(seed)

tree <- read.tree('species_tree.nwk')
traits <- read.csv('traits.csv', row.names = 1)
x <- setNames(traits$state, rownames(traits))

if (anyDuplicated(tree$tip.label)) {
    stop("Tree has duplicate tip labels; provide one unique label per taxon.", call. = FALSE)
}
if (anyDuplicated(names(x))) {
    stop("Trait table has duplicate taxon row names; keep one row per tree tip.", call. = FALSE)
}
missing_traits <- setdiff(tree$tip.label, names(x))
extra_traits <- setdiff(names(x), tree$tip.label)
if (length(missing_traits) || length(extra_traits)) {
    show_taxa <- function(values) if (length(values)) paste(values, collapse = ", ") else "none"
    stop(sprintf(
        "Tree and trait taxa must match exactly. Missing traits for tree tips: %s. Trait-only taxa: %s. Align row names to tree tip labels and rerun.",
        show_taxa(missing_traits), show_taxa(extra_traits)
    ), call. = FALSE)
}
if (anyNA(x) || any(!nzchar(trimws(as.character(x))))) {
    stop("Trait states must be non-missing and non-blank for every tree tip.", call. = FALSE)
}
observed_states <- unique(as.character(x))
if (length(observed_states) < 2L) {
    stop("At least two observed trait states are required to fit a transition model.", call. = FALSE)
}

# Fit Mk model and pick ARD vs SYM vs ER by AIC
fit_er <- fitMk(tree, x, model = 'ER')
fit_sym <- fitMk(tree, x, model = 'SYM')
fit_ard <- fitMk(tree, x, model = 'ARD')
aic <- sapply(list(fit_er, fit_sym, fit_ard), AIC)
best_model <- c('ER', 'SYM', 'ARD')[which.min(aic)]

# Stochastic mapping under best model, nsim >= 1000 (raise for asymmetric rates)
smaps <- make.simmap(tree, x, model = best_model, nsim = 1000, pi = 'fitzjohn')
node_pp <- summary(smaps)$ace          # posterior probabilities per state per node

# Hidden-rate alternative for clade-rate heterogeneity
hmm_fit <- corHMM(phy = tree, data = data.frame(species = names(x), trait = x),
                   rate.cat = 2, model = 'ARD', node.states = 'marginal')
hmm_fit$states                          # marginal ancestral state matrix (rows: nodes)

result_metadata <- list(seed = seed, selected_model = best_model,
                        taxa = length(tree$tip.label), states = sort(observed_states),
                        stochastic_maps = length(smaps))
print(result_metadata)
