# Inputs expected from the calling analysis: a phylo object `tree` and a
# named numeric vector `traits` indexed by the tree's tip labels.
library(phytools)
library(geiger)

if (!inherits(tree, "phylo")) {
    stop("Provide a phylo object `tree`.", call. = FALSE)
}
trait_values <- if (is.data.frame(traits) && "mass" %in% names(traits)) {
    setNames(traits$mass, rownames(traits))
} else {
    traits
}
if (is.null(names(trait_values))) {
    stop("Provide `traits` as a named numeric vector or a data frame with numeric `mass` and taxon row names.", call. = FALSE)
}
if (anyDuplicated(tree$tip.label) || anyDuplicated(names(trait_values))) {
    stop("Tree tips and trait names must each be unique.", call. = FALSE)
}
missing_traits <- setdiff(tree$tip.label, names(trait_values))
extra_traits <- setdiff(names(trait_values), tree$tip.label)
if (length(missing_traits) || length(extra_traits)) {
    show_taxa <- function(values) if (length(values)) paste(values, collapse = ", ") else "none"
    stop(sprintf("Tree and trait taxa must match exactly. Missing traits: %s. Trait-only taxa: %s.",
                 show_taxa(missing_traits), show_taxa(extra_traits)), call. = FALSE)
}
if (!is.numeric(trait_values) || any(!is.finite(trait_values))) {
    stop("Continuous trait values must be finite numeric values.", call. = FALSE)
}
traits <- trait_values[tree$tip.label]
tree_p <- multi2di(tree)

models <- c("BM", "OU", "EB", "lambda", "kappa", "delta")
fits <- setNames(lapply(models, function(model) fitContinuous(tree_p, traits, model = model)), models)
aic_tbl <- data.frame(model = models, AIC = vapply(fits, function(fit) fit$opt$aic, numeric(1)))
if (any(!is.finite(aic_tbl$AIC))) stop("At least one continuous model returned non-finite AIC.", call. = FALSE)
best <- aic_tbl$model[which.min(aic_tbl$AIC)]

lambda_fit <- phylosig(tree_p, traits, method = "lambda", test = TRUE)
K_fit <- phylosig(tree_p, traits, method = "K", test = TRUE)

if (best == "BM") {
    reconstruction <- fastAnc(tree_p, traits, CI = TRUE)
    result <- list(status = "reconstructed", selected_model = "BM",
                   estimates = reconstruction$ace, uncertainty = reconstruction$CI,
                   uncertainty_type = "fastAnc 95% confidence intervals")
} else if (best == "OU") {
    # OUwie.anc accepts a fitted OUwie object (not a `data=` argument). Its
    # provider documentation says results are point estimates without
    # uncertainty and should be used for visualization/model intuition only.
    if (!ape::is.ultrametric(tree_p)) {
        stop("OUwie.anc requires an ultrametric tree; no OU reconstruction emitted.", call. = FALSE)
    }
    library(OUwie)
    ou_data <- data.frame(species = names(traits), regime = rep("1", length(traits)),
                          trait = unname(traits), stringsAsFactors = FALSE)
    # OUwie 3.0.3's identifiability diagnostic returns a missing-value error
    # for this validated single-regime OU1 setup; regimes are explicit and the
    # fitted object's finite likelihood/AICc are checked below.
    ou_fit <- OUwie(tree_p, ou_data, model = "OU1", simmap.tree = FALSE,
                    quiet = TRUE, check.identify = FALSE)
    if (!inherits(ou_fit, "OUwie") || !is.finite(ou_fit$AICc)) {
        stop("OUwie did not return a valid OU1 fitted object; no reconstruction emitted.", call. = FALSE)
    }
    reconstruction <- OUwie.anc(ou_fit, knowledge = TRUE)
    result <- list(status = "exploratory_point_estimates_only", selected_model = "OU",
                   estimates = reconstruction, uncertainty = NULL,
                   uncertainty_type = "not provided by OUwie.anc")
    warning("OUwie.anc estimates have no uncertainty and are intended for visualization/model intuition; they do not satisfy an uncertainty-bearing ASR deliverable.",
            call. = FALSE)
} else {
    stop(sprintf(paste0("Selected %s model is not implemented as an ancestral reconstruction here. ",
                        "No original-tree fastAnc/contMap result was substituted. Use a method that reconstructs ",
                        "under this fitted transformation and reports suitable uncertainty, then compare it independently."), best),
         call. = FALSE)
}

result$model_comparison <- aic_tbl
result$phylogenetic_signal <- list(lambda = lambda_fit, K = K_fit)
print(result)
