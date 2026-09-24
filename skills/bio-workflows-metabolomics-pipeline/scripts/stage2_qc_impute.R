# Stage 2 -- QC, drift correction, PQN normalization and mechanism-aware imputation (pmp + imputeLCMD).
# Inputs (define before source()): feat (features x samples matrix), defs (data.frame, one row per
#   feature), sample_class (per-sample, pooled QCs labelled 'QC'), injection_order, batch_id (per-sample).
# Outputs (left in the calling environment): imputed (NA-free, features x study+QC samples),
#   sample_class (OVERWRITTEN: re-synced after wiped samples are dropped, so re-define it before
#   re-running this script), normalized, nm.
# Usage: source('scripts/stage2_qc_impute.R')
# Checked on pmp 1.18.0, imputeLCMD 2.1 (R 4.4.3).
stopifnot(exists('feat'), exists('defs'), exists('sample_class'), exists('injection_order'), exists('batch_id'))
library(pmp)
library(imputeLCMD)
# feat (Stage 1 featureValues() output) is already features x samples -- the pmp convention.
# No transpose needed; assert it instead of assuming, so an API change is caught, not silently
# masked (pmp's own check_peak_matrix re-transposes with no warning when it can, which is what
# hid this exact bug before: a wrong `t(feat)` here still "worked").
fm <- feat
stopifnot(nrow(fm) == nrow(defs), ncol(fm) == length(sample_class))

filtered <- filter_peaks_by_fraction(fm, classes = sample_class, min_frac = 0.5, qc_label = 'QC')
corrected <- QCRSC(df = filtered, order = injection_order, batch = batch_id,
                   classes = sample_class, spar = 0, minQC = 5, qc_label = 'QC')  # CV-selected spline
rsd_filtered <- filter_peaks_by_rsd(corrected, max_rsd = 30, classes = sample_class, qc_label = 'QC')
normalized <- pqn_normalisation(rsd_filtered, classes = sample_class, qc_label = 'QC')

# QCRSC silently returns an ENTIRE batch as all-NA when that batch has fewer QCs than minQC (see
# normalization-qc's Common Errors) -- verified on real MTBLS79 data: 5 of 8 batches had only 4
# QCs (< minQC=5), wiping 90 of 172 samples. A wholly-missing sample carries NO measured
# information to impute -- report and drop it; do not fabricate a profile for it.
nm <- as.matrix(normalized)
wiped <- colMeans(is.na(nm)) == 1
if (any(wiped)) {
  cat(sum(wiped), 'of', ncol(nm), 'samples came back all-NA (QCRSC: their batch had < minQC QCs)',
      '-- dropping, not imputing:\n')
  print(table(batch_id[wiped]))
  nm <- nm[, !wiped, drop = FALSE]
  sample_class <- sample_class[!wiped]   # keep every per-sample vector in sync for Stage 4
}

# The remaining holes are the sparse, mechanism-driven kind imputation is actually valid for.
# QRILC (MNAR / left-censored) is only valid on log-scale intensities; on raw intensities it
# silently draws negative (impossible) values. Round-trip through log2/2^x, per normalization-qc.
log_mat <- log2(nm)
set.seed(123)   # impute.QRILC draws random values; seed it so the imputed matrix is reproducible
imputed <- 2^(impute.QRILC(log_mat, tune.sigma = 1)[[1]])
stopifnot(min(imputed, na.rm = TRUE) >= 0, !anyNA(imputed))  # no NAs may reach Stage 4 (opls() cannot tolerate them)
