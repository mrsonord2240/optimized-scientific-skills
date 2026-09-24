# Purpose : class-based internal-standard normalization (normalize_istd) with two guards that
#           normalize_istd() lacks: lipids with Class = NA are dropped loudly, and any class with
#           zero recognized internal standards stops the run (else it passes through uncorrected, factor = 1).
# Inputs  : Skyline export CSV(s) + a sample-annotation CSV, or none to use lipidr's bundled raw export
#           (normalize_istd() refuses already-normalized data such as lipidr's data_normalized).
# Usage   : Rscript istd_normalize.R out_dir [annotation.csv skyline1.csv [skyline2.csv ...]]
# Outputs : out_dir/istd_normalized.csv (log2 Area per lipid x sample) and out_dir/lipidclass_sd.png
# Checked : lipidr 2.20.0
suppressMessages(library(lipidr))

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 1) stop('usage: Rscript istd_normalize.R out_dir [annotation.csv skyline1.csv ...]')
out_dir <- a[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# normalize_istd() refuses already-normalized data, so start from a raw Skyline export.
if (length(a) >= 3) {
  d_raw <- add_sample_annotation(read_skyline(a[3:length(a)]), a[2])
} else {
  datadir <- system.file('extdata', package = 'lipidr')
  d_raw <- add_sample_annotation(
    read_skyline(list.files(datadir, 'A1_data.csv|F1_data.csv|F2_data.csv', full.names = TRUE)),
    file.path(datadir, 'clin.csv')
  )
}

# table() drops NA groups, so an unclassified lipid (Class = NA) would reach normalize_istd() unflagged.
# Drop it loudly (lipidr's own bundled export has one such row, 'PI 34:1p'); fix the name instead if it is real.
unclassified <- is.na(rowData(d_raw)$Class)
if (any(unclassified)) {
  warning(sprintf('Dropping %d lipid(s) with Class = NA (unparsed names): %s',
                  sum(unclassified), paste(rowData(d_raw)$Molecule[unclassified], collapse = ', ')))
  d_raw <- d_raw[!unclassified, ]
}

# GUARD (non-negotiable, not optional): normalize_istd() does NOT enforce this. Any class with
# zero recognized standards passes through with a silent correction factor of 1 -- i.e. completely
# uncorrected data reported as if it had been normalized (verified against lipidr 2.20.0's
# internal normalize_istd(): `if (length(istd_list[[i]]) == 0) f <- 1`). Fail loud instead:
istd_coverage <- table(rowData(d_raw)$Class, rowData(d_raw)$istd)
uncovered <- rownames(istd_coverage)[
  !('TRUE' %in% colnames(istd_coverage)) | istd_coverage[, 'TRUE'] == 0
]
if (length(uncovered) > 0) {
  stop(sprintf(
    'No recognized internal standard for class(es): %s -- normalize_istd() would pass these through uncorrected (factor=1), not normalized. Add a labeled standard for this class or exclude it from ISTD-normalized reporting.',
    paste(uncovered, collapse = ', ')
  ))
}

d_istd <- normalize_istd(d_raw, measure = 'Area', exclude = 'blank', log = TRUE)

# Class-level summary is only valid WITHIN a class unless per-class response factors were calibrated:
# cross-class molar ratios (e.g. 'PE is 3x PC') carry head-group response bias and are not licensed here.
ggplot2::ggsave(file.path(out_dir, 'lipidclass_sd.png'), plot_lipidclass(d_istd, 'sd'), width = 8, height = 5)

m <- assay(d_istd, 'Area')
write.csv(cbind(Molecule = rownames(m), as.data.frame(m, check.names = FALSE)),
          file.path(out_dir, 'istd_normalized.csv'), row.names = FALSE)
cat(sprintf('ISTD-normalized %d lipids x %d samples -> %s\n', nrow(d_istd), ncol(d_istd), out_dir))
