# Reference: ropls 1.34+ | Verify API if version differs
# Unsupervised PCA with explicit scaling + Hotelling T2 multivariate outlier check.
# ropls exposes no ready accessor for T2 (pca@suppLs$outlierDF is NULL for a plain PCA fit),
# so it is computed directly from the scores.
#
# Input : CSV, features in rows x samples in columns (first column = feature IDs), NA-free
#         (impute upstream, see metabolomics/normalization-qc; opls() does not tolerate NA).
# Output: R2X per component, the T2 critical value, flagged samples; optional CSV of per-sample T2.
# Usage : Rscript scripts/pca_hotelling.R <peaks.csv> [scaleC=pareto] [alpha=0.05] [out_t2.csv]
#         source('scripts/pca_hotelling.R')   # defines hotelling_t2(scores, alpha) only
library(ropls)

# scores: samples x components (getScoreMN(pca)). Returns per-sample T2, the F-based critical
# value at 1 - alpha, and the names of samples beyond it.
hotelling_t2 <- function(scores, alpha = 0.05) {
    A <- ncol(scores); N <- nrow(scores)
    lambda <- apply(scores, 2, function(col) sum(col^2) / (N - 1))     # per-component variance
    t2 <- rowSums(sweep(scores^2, 2, lambda, '/'))                     # per-sample Hotelling T2
    crit <- A * (N - 1) / (N - A) * qf(1 - alpha, A, N - A)            # F-based critical value
    list(t2 = t2, crit = crit, outliers = rownames(scores)[t2 > crit])
}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 1) stop('usage: Rscript pca_hotelling.R <peaks.csv> [scaleC=pareto] [alpha=0.05] [out_t2.csv]')
    peaks_csv <- args[1]
    scaleC <- if (length(args) >= 2) args[2] else 'pareto'   # ropls default is "standard" (UV), NOT Pareto
    alpha <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
    out_csv <- if (length(args) >= 4) args[4] else NA

    peaks <- read.csv(peaks_csv, row.names = 1, check.names = FALSE)
    x <- t(as.matrix(peaks))                                 # samples x features, as opls() wants
    if (anyNA(x)) stop(sum(is.na(x)), ' NA values: impute upstream (metabolomics/normalization-qc)')

    pca <- opls(x, scaleC = scaleC, fig.pdfC = 'none', info.txtC = 'none')
    print(getSummaryDF(pca)[, 'R2X(cum)', drop = FALSE])     # R2X(cum) per component
    res <- hotelling_t2(getScoreMN(pca), alpha)
    cat(sprintf('N samples=%d, A components=%d, T2 crit (%g%%)=%.2f\n',
                nrow(x), ncol(getScoreMN(pca)), 100 * (1 - alpha), res$crit))
    cat('Outliers flagged:', length(res$outliers), 'of', nrow(x), '\n')
    if (length(res$outliers) > 0) print(res$outliers)
    if (!is.na(out_csv)) write.csv(data.frame(sample = names(res$t2), T2 = res$t2), out_csv, row.names = FALSE)
    # Tight pooled-QC clustering in the score plot = trustworthy run; QC scatter = analytical
    # variance dominates. For a plotted version: plot(pca, typeVc = 'outlier').
}
