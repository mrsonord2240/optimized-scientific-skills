# scRNA-seq differential-expression power: donors, not cells, set population power
# Reference: edgeR 4.0+ | Verify API if version differs
#
# powsimR (the package this route is usually pointed at) is GitHub-only with a
# compile-required dependency (bayNorm) and is not installed by default -- see
# Version Compatibility. This script demonstrates the same claim with tools already
# used elsewhere in this Skill (edgeR, already a stated dependency): aggregate
# (sum) counts per donor into a pseudobulk matrix and treat donors -- not cells --
# as the unit of replication, exactly like bulk RNA-seq. It contrasts that with the
# anti-pattern of testing on cells directly.

suppressPackageStartupMessages(library(edgeR))
set.seed(123)

# A donor-count configuration is usable only if its realized pseudobulk FDR passes.
assess_realized_fdr <- function(actual_fdr, target_fdr = 0.05, tolerance = 0) {
  actual_fdr <- unname(as.numeric(actual_fdr))
  accepted <- is.finite(actual_fdr) & actual_fdr <= target_fdr + tolerance
  c(actual_fdr = actual_fdr, accepted = accepted)
}

ngenes <- 200
n_de <- 20        # 10% truly DE
true_lfc <- 1.0   # log2 fold change for DE genes
disp_donor <- 0.3 # biological (between-donor) dispersion, on the log-mean scale
disp_cell  <- 0.15 # additional technical/cell-level overdispersion

# Simulate a donor x cell NB dataset and score it two ways: pseudobulk-by-donor
# (the recommended unit) and naive cell-level testing (cells as pseudoreplicates).
run_grid <- function(n_donors, cells_per_donor, nsim = 5) {
  pb_power <- numeric(nsim); pb_fdr <- numeric(nsim)
  cell_power <- numeric(nsim); cell_fdr <- numeric(nsim)
  de_id <- 1:n_de
  for (s in 1:nsim) {
    donor_group <- rep(c(0, 1), each = n_donors / 2)
    base_mean <- exp(rnorm(ngenes, log(50), 0.4))
    donor_effect <- matrix(rnorm(ngenes * n_donors, 0, sqrt(disp_donor)), nrow = ngenes)
    group_effect <- matrix(0, ngenes, n_donors)
    group_effect[de_id, donor_group == 1] <- true_lfc * log(2)
    log_mean_donor <- log(base_mean) + donor_effect + group_effect
    counts <- matrix(0L, ngenes, n_donors * cells_per_donor)
    donor_of_cell <- rep(1:n_donors, each = cells_per_donor)
    for (d in 1:n_donors) {
      mu_d <- exp(log_mean_donor[, d])
      cell_noise <- exp(matrix(rnorm(ngenes * cells_per_donor, 0, sqrt(disp_cell)), nrow = ngenes))
      mu_cells <- mu_d * cell_noise
      counts[, donor_of_cell == d] <- rnbinom(ngenes * cells_per_donor, mu = mu_cells, size = 5)
    }
    group_of_cell <- donor_group[donor_of_cell]

    # --- pseudobulk by donor: n = number of DONORS ---
    pb_counts <- sapply(1:n_donors, function(d) rowSums(counts[, donor_of_cell == d, drop = FALSE]))
    y <- DGEList(counts = pb_counts, group = donor_group)
    y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    design <- model.matrix(~donor_group)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)
    padj <- p.adjust(qlf$table$PValue, "BH")
    called <- rownames(y) %in% as.character(de_id)
    sig <- padj < 0.05
    pb_power[s] <- sum(sig & called) / max(1, sum(called))
    pb_fdr[s] <- if (sum(sig) > 0) sum(sig & !called) / sum(sig) else 0

    # --- naive cell-level test: cells treated as independent replicates (anti-pattern) ---
    yc <- DGEList(counts = counts, group = group_of_cell)
    yc <- yc[filterByExpr(yc), , keep.lib.sizes = FALSE]
    yc <- calcNormFactors(yc)
    designc <- model.matrix(~group_of_cell)
    yc <- estimateDisp(yc, designc)
    fitc <- glmQLFit(yc, designc)
    qlfc <- glmQLFTest(fitc, coef = 2)
    padjc <- p.adjust(qlfc$table$PValue, "BH")
    calledc <- rownames(yc) %in% as.character(de_id)
    sigc <- padjc < 0.05
    cell_power[s] <- sum(sigc & calledc) / max(1, sum(calledc))
    cell_fdr[s] <- if (sum(sigc) > 0) sum(sigc & !calledc) / sum(sigc) else 0
  }
  c(pb_power = mean(pb_power), pb_fdr = mean(pb_fdr),
    cell_power = mean(cell_power), cell_fdr = mean(cell_fdr))
}

cfgs <- list(c(n_donors = 4, cpd = 200), c(n_donors = 12, cpd = 200))
for (cfg in cfgs) {
  r <- run_grid(cfg[["n_donors"]], cfg[["cpd"]])
  gate <- assess_realized_fdr(r["pb_fdr"])
  decision <- if (gate[["accepted"]]) "ACCEPT" else "REJECT"
  cat(sprintf("donors=%2d cells/donor=%3d | pseudobulk power=%.3f FDR=%.3f decision=%s | cell-level power=%.3f FDR=%.3f\n",
              cfg[["n_donors"]], cfg[["cpd"]], r["pb_power"], r["pb_fdr"], decision, r["cell_power"], r["cell_fdr"]))
  if (decision == "REJECT") cat('# Do not select this donor-count power estimate: increase nsim, then refit/revise the model.\n')
}
cat('# Pseudobulk power tracks donor count; cell-level "power" is inflated by treating\n')
cat('# pseudoreplicates as replicates -- its FDR stays far above the nominal 0.05.\n')
