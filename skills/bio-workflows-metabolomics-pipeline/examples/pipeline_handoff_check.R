# Reference: pmp 1.18.0, ropls 1.38.0, imputeLCMD 2.1 | checked 2026-09-16
#
# Second bundled example -- chains REAL pmp/ropls calls (not base-R stand-ins) across the exact
# Stage1->Stage2->Stage4 hand-offs eval_report_bio-workflows-metabolomics-pipeline_result.json's
# Inputs 1-2 found broken:
#   - Input 1: Stage1->Stage2 `fm <- t(feat)` inverted the orientation its own comment claimed
#   - Input 2: Stage2->Stage4 crashed on real MTBLS79 data (55.95% NA) because imputation was
#     named only in a comment, never called
# metabolomics_workflow.R (the first bundled example) only exercises Stages 2+4 in base R and
# never touches pmp/ropls, so it could not have caught either regression. This one does, on a
# small synthetic multi-batch table shaped like the real failure (one batch with fewer QCs than
# QCRSC's own minQC default).

suppressMessages({library(pmp); library(imputeLCMD); library(ropls)})
set.seed(42)

n_features <- 150
b1_qc <- 6; b1_bio <- 10     # batch 1: enough QCs
b2_qc <- 3; b2_bio <- 10     # batch 2: below minQC=5 -- the real MTBLS79 failure mode (5 of 8
                             # batches there had only 4 QCs)

sample_names <- c(paste0('B1QC', 1:b1_qc), paste0('B1S', 1:b1_bio),
                   paste0('B2QC', 1:b2_qc), paste0('B2S', 1:b2_bio))
sample_class <- c(rep('QC', b1_qc), rep(c('Control', 'Treatment'), each = b1_bio / 2),
                   rep('QC', b2_qc), rep(c('Control', 'Treatment'), each = b2_bio / 2))
batch_id <- c(rep(1, b1_qc + b1_bio), rep(2, b2_qc + b2_bio))
injection_order <- seq_along(sample_names)
n_samples <- length(sample_names)

base_intensity <- 2^runif(n_features, 8, 20)
feat <- base_intensity * matrix(rlnorm(n_features * n_samples, 0, 0.15), nrow = n_features)
true_hits <- 1:15
trt_idx <- sample_class == 'Treatment'
feat[true_hits, trt_idx] <- feat[true_hits, trt_idx] * 3   # planted 3-fold effect (|log2FC|~1.6)
drift <- 1 - 0.2 * (injection_order / max(injection_order))
feat <- sweep(feat, 2, drift, '*')
rownames(feat) <- paste0('FT', sprintf('%03d', 1:n_features))
colnames(feat) <- sample_names
defs <- data.frame(mzmed = runif(n_features, 80, 900), row.names = rownames(feat))

cat('=== Stage 1 (stand-in): featureValues()-shaped output ===\n')
cat(nrow(feat), 'features x', ncol(feat), 'samples\n\n')

cat('=== Stage 1 -> Stage 2 hand-off (Input 1 regression guard) ===\n')
# SKILL.md: `fm <- feat` -- NOT `t(feat)`. featureValues() is already features x samples.
fm <- feat
stopifnot(nrow(fm) == nrow(defs), ncol(fm) == length(sample_class))
cat('No transpose applied; dimension check passed (', nrow(fm), 'x', ncol(fm), ').\n\n')

cat('=== Stage 2: real pmp calls ===\n')
filtered <- filter_peaks_by_fraction(fm, classes = sample_class, min_frac = 0.5, qc_label = 'QC')
cat('After filter_peaks_by_fraction:', nrow(filtered), 'of', n_features, 'features kept\n')

corrected <- QCRSC(df = filtered, order = injection_order, batch = batch_id, classes = sample_class,
                    spar = 0, minQC = 5, qc_label = 'QC')
cat('After QCRSC:', nrow(corrected), 'x', ncol(corrected), '\n')

rsd_filtered <- filter_peaks_by_rsd(corrected, max_rsd = 30, classes = sample_class, qc_label = 'QC')
cat('After filter_peaks_by_rsd:', nrow(rsd_filtered), 'of', nrow(corrected), 'features kept\n')

normalized <- pqn_normalisation(rsd_filtered, classes = sample_class, qc_label = 'QC')
nm <- as.matrix(normalized)
cat('After pqn_normalisation:', nrow(nm), 'x', ncol(nm), '| NA cells:', sum(is.na(nm)),
    sprintf('(%.1f%%)', 100 * mean(is.na(nm))), '\n\n')

cat('=== Input 2 regression guard: drop QCRSC-wiped samples, then impute the sparse rest ===\n')
wiped <- colMeans(is.na(nm)) == 1
cat(sum(wiped), 'of', ncol(nm), 'samples came back all-NA (batch 2, minQC=5 > its', b2_qc, 'QCs)\n')
nm_keep <- nm[, !wiped, drop = FALSE]
class_keep <- sample_class[!wiped]

log_mat <- log2(nm_keep)
imputed <- 2^(impute.QRILC(log_mat, tune.sigma = 1)[[1]])
cat('After impute.QRILC (log2/2^x round-trip):', sum(is.na(imputed)), 'NA cells remain; min value:',
    round(min(imputed, na.rm = TRUE), 3), '\n')
stopifnot(!anyNA(imputed), min(imputed, na.rm = TRUE) >= 0)
cat('This is exactly the hand-off that crashed opls() before the fix -- now NA-free.\n\n')

cat('=== Stage 4: real ropls, guarded fit ===\n')
# Verbatim copy of metabolomics/statistical-analysis's fit_discriminant_guarded, kept in sync
# with that Skill (see its fix log if this diverges) so this example stays a faithful hand-off test.
fit_discriminant_guarded <- function(x, y, scaleC, permI = 1000, crossvalI = 7) {
    m <- opls(x, y, predI = 1, orthoI = NA, scaleC = scaleC, permI = permI,
              crossvalI = crossvalI, fig.pdfC = 'none', info.txtC = 'none')
    if (nrow(getSummaryDF(m)) > 0) return(list(model = m, type = 'OPLS-DA'))
    m2 <- opls(x, y, predI = 1, orthoI = 0, scaleC = scaleC, permI = permI,
               crossvalI = crossvalI, fig.pdfC = 'none', info.txtC = 'none')
    if (nrow(getSummaryDF(m2)) > 0) return(list(model = m2, type = 'PLS-DA (OPLS-DA fallback)'))
    stop('Neither OPLS-DA nor the PLS-DA fallback produced a usable model.')
}

study <- class_keep != 'QC'
group <- factor(class_keep[study])
x <- t(imputed)[study, ]
stopifnot(!anyNA(x))
cat('Stage2->Stage4 hand-off shape:', nrow(x), 'x', ncol(x), '| groups:',
    paste(levels(group), collapse = ' vs '), '\n')

fit <- fit_discriminant_guarded(x, group, scaleC = 'pareto', permI = 1000)
summ <- getSummaryDF(fit$model)
stopifnot(nrow(summ) > 0)
cat('Model type actually fit:', fit$type, '| getSummaryDF rows:', nrow(summ), '\n')
print(summ[, c('R2X(cum)', 'R2Y(cum)', 'Q2(cum)', 'pR2Y', 'pQ2')])

cat('\nEnd-to-end Stage1->Stage2->Stage4 hand-off completed: correct orientation, no NAs reached\n')
cat('opls(), and a usable model was produced -- both regressions from the audit stay fixed.\n')
