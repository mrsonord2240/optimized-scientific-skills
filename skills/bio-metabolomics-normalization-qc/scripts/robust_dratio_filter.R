# Purpose: keep only reproducible features, using the robust (MAD-based) D-ratio (technical SD /
#          biological SD) and QC RSD. MAD because MS intensities are right-skewed.
# Input:   CSV with samples in ROWS and features in COLUMNS, plus one sample-type column
#          (QC rows carry the QC label). Feature columns must be numeric.
# Output:  the same CSV restricted to the kept feature columns (sample-type column retained).
# Usage:   Rscript robust_dratio_filter.R peaks.csv out.csv [type_col=sample_type] [qc_label=QC] \
#                  [dratio_max=0.5] [rsd_max=0.3]
#          or in R: source('scripts/robust_dratio_filter.R'); kept <- robust_dratio_filter(data, is_qc)

library(matrixStats)

robust_dratio_filter <- function(data, is_qc, dratio_max = 0.5, rsd_max = 0.3) {
    qc <- as.matrix(data[is_qc, ])
    bio <- as.matrix(data[!is_qc, ])
    # MAD-based (robust) form, because MS intensities are right-skewed
    sd_qc <- colMads(qc, na.rm = TRUE)
    sd_bio <- colMads(bio, na.rm = TRUE)
    dratio <- sd_qc / sd_bio
    rsd <- colSds(qc, na.rm = TRUE) / colMeans(qc, na.rm = TRUE)
    keep <- dratio <= dratio_max & rsd <= rsd_max
    keep[is.na(keep)] <- FALSE
    message(sprintf('D-ratio<=%.2f & RSD<=%.0f%%: kept %d / %d features',
                    dratio_max, rsd_max * 100, sum(keep), ncol(data)))
    data[, keep]
}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) stop('usage: Rscript robust_dratio_filter.R peaks.csv out.csv [type_col] [qc_label] [dratio_max] [rsd_max]')
    type_col <- if (length(args) >= 3) args[3] else 'sample_type'
    qc_label <- if (length(args) >= 4) args[4] else 'QC'
    dratio_max <- if (length(args) >= 5) as.numeric(args[5]) else 0.5
    rsd_max <- if (length(args) >= 6) as.numeric(args[6]) else 0.3
    tbl <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
    is_qc <- tbl[[type_col]] == qc_label
    feats <- tbl[, setdiff(colnames(tbl), type_col), drop = FALSE]
    kept <- robust_dratio_filter(feats, is_qc, dratio_max, rsd_max)
    write.csv(cbind(tbl[type_col], kept), args[2], row.names = FALSE)
}
