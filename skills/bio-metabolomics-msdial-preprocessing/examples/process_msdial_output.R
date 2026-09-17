# Reference: R 4.3+ (base only) | Verify API if version differs
# Parses an MS-DIAL alignment-result export, splits metadata from per-sample
# intensities, and filters honestly on Fill% and MS/MS support. Generates a
# synthetic MS-DIAL-style table matching the real MS-DIAL 5.x console's actual
# output conventions (verified against a real MSDIALCUI.exe 5.5.260820 run on
# real mzML LC-MS/MS data, 2026-09-16):
#   - 4 header rows above the column header (class / file type / injection
#     order / batch), each written as "<tabs>\t<row label>\t<value 1>...<value N>"
#     -- ONE label cell immediately followed by exactly N value cells (N =
#     number of samples), not a per-column label on every sample column.
#   - "Fill %" is a 0-1 fraction (e.g. 1.00 = 100%), not a 0-100 percentage.
#   - "MS/MS assigned" is title-case text ("True"/"False"), not "TRUE"/"FALSE".
#   - Extra annotation/QC columns (Formula, INCHIKEY, dot-product scores,
#     spectra text, ...) sit between the known metadata columns and the real
#     per-sample intensity columns -- exactly the layout that breaks a naive
#     "everything after the last known metadata column is a sample column" cut.
# All output goes to tempdir().

make_synthetic_export <- function(path, n_features = 60, n_samples = 8) {
    set.seed(42)
    tags <- sample(c('Metabolite', 'Lipid', 'Suggested', 'Unknown'), n_features,
                   replace = TRUE, prob = c(0.25, 0.15, 0.2, 0.4))
    has_msms <- ifelse(tags %in% c('Metabolite', 'Lipid'),
                       sample(c('True', 'False'), n_features, replace = TRUE, prob = c(0.7, 0.3)),
                       'False')
    # Fill % is 0-1 in real MS-DIAL 5.x output (e.g. 1.00 = 100%), not 0-100.
    fill_frac <- round(sample(20:100, n_features, replace = TRUE) / 100, 2)
    meta <- data.frame(
        'Alignment ID' = seq_len(n_features),
        'Average Rt(min)' = round(runif(n_features, 0.5, 15), 3),
        'Average Mz' = round(runif(n_features, 80, 900), 4),
        'Metabolite name' = ifelse(tags == 'Unknown', 'Unknown', paste0(tags, '_', seq_len(n_features))),
        'Adduct type' = sample(c('[M+H]+', '[M+Na]+', '[M+NH4]+', '[M-H]-'), n_features, replace = TRUE),
        'Fill %' = fill_frac,
        'MS/MS assigned' = has_msms,
        'Annotation tag (VS1.0)' = tags,
        check.names = FALSE
    )
    # Real exports carry extra annotation/QC columns (Formula, Ontology, INCHIKEY,
    # SMILES, RT/m/z/MS-MS matched flags, dot-product scores, spectra text, ...)
    # between the last known metadata column and the first real sample column.
    # Reproduce a few so the sample-column detector below is actually exercised.
    extra_cols <- data.frame(
        'Formula' = 'null', 'Ontology' = 'null', 'INCHIKEY' = 'null',
        'Total score' = round(runif(n_features, 0, 100), 2),
        'MS/MS spectrum' = 'null',
        check.names = FALSE
    )
    sample_names <- paste0('Sample_', seq_len(n_samples))
    intensities <- as.data.frame(matrix(round(rlnorm(n_features * n_samples, 11, 1.6)),
                                        nrow = n_features, dimnames = list(NULL, sample_names)),
                                 check.names = FALSE)
    body <- cbind(meta, extra_cols, intensities)

    # Real MS-DIAL exports prepend 4 metadata rows (class / file type / injection
    # order / batch) above the column header, each row written as one label cell
    # followed immediately by N value cells (one per sample) -- not a label
    # repeated under every sample column. Reproduce that exact layout.
    n_annot <- ncol(meta) + ncol(extra_cols)
    row_of <- function(label, values) c(rep('', n_annot - 1), label, values)
    header_block <- rbind(
        row_of('Class', rep('0', n_samples)),
        row_of('File type', rep('Sample', n_samples)),
        row_of('Injection order', as.character(seq_len(n_samples))),
        row_of('Batch ID', rep('1', n_samples))
    )
    writeLines(apply(header_block, 1, paste, collapse = '\t'), path)
    suppressWarnings(write.table(body, path, sep = '\t', append = TRUE,
                                 row.names = FALSE, col.names = TRUE, quote = FALSE))
}

export_path <- file.path(tempdir(), 'AlignResult-synthetic.mdalign')
make_synthetic_export(export_path)

# Real MS-DIAL exports mark where the per-sample columns start with a single
# "Class" label cell in the first header row, immediately followed by one value
# cell per sample -- NOT a label repeated under every sample column, and NOT
# simply "everything after the last metadata column named in this script".
# Real exports interleave extra annotation/QC columns (Formula, INCHIKEY,
# dot-product scores, spectra text, ...) before the real per-sample columns,
# so anchor on the "Class" cell itself rather than guessing a fixed cut point.
header_row1 <- as.character(read.table(export_path, sep = '\t', header = FALSE,
                                        nrows = 1, check.names = FALSE,
                                        colClasses = 'character', comment.char = '')[1, ])
class_idx <- which(header_row1 == 'Class')
stopifnot(length(class_idx) == 1)

msdial <- read.csv(export_path, sep = '\t', skip = 4, check.names = FALSE)
cat('Parsed', nrow(msdial), 'features x', ncol(msdial), 'columns\n')

meta_cols <- intersect(c('Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name',
                         'Adduct type', 'Fill %', 'MS/MS assigned', 'Annotation tag (VS1.0)'),
                       colnames(msdial))
sample_cols <- colnames(msdial)[(class_idx + 1):ncol(msdial)]

feature_info <- msdial[, meta_cols]
intensity <- as.matrix(msdial[, sample_cols])
rownames(intensity) <- msdial[['Alignment ID']]
cat('Split into', length(meta_cols), 'metadata cols and', length(sample_cols), 'sample cols\n')
storage.mode(intensity) <- 'numeric'  # confirms every detected sample col is truly numeric

# Fill% floor: below this a feature is mostly gap-filled noise-floor integrals, not
# measurements, so it would fabricate intensity for truly below-detection samples.
# MS-DIAL 5.x reports Fill % as a 0-1 fraction (confirmed on a real console run,
# MSDIALCUI.exe 5.5.260820, 2026-09-16) -- 0.70, not 70, is the floor to compare against.
fill_floor <- 0.70
keep_fill <- feature_info[['Fill %']] >= fill_floor

# A named hit without MS/MS is at best a putative (MSI Level 3) ID; require MS/MS for
# identity. Compare case-insensitively: MS-DIAL 5.x's real "MS/MS assigned" column is
# title-case text ("True"/"False"), not all-caps "TRUE"/"FALSE" -- an exact-case
# comparison silently keeps zero features on real output.
has_msms <- tolower(trimws(feature_info[['MS/MS assigned']])) == 'true'
feature_info$msi_level <- ifelse(feature_info[['Annotation tag (VS1.0)']] %in% c('Metabolite', 'Lipid') & has_msms, 2L,
                          ifelse(feature_info[['Annotation tag (VS1.0)']] == 'Suggested', 3L, NA_integer_))

filtered_intensity <- intensity[keep_fill, ]
filtered_info <- feature_info[keep_fill, ]
cat('After Fill% >=', fill_floor, ':', nrow(filtered_intensity), '/', nrow(intensity), 'features\n')

cat('\nMSI confidence levels (survivors):\n')
print(table(filtered_info$msi_level, useNA = 'ifany'))

out_path <- file.path(tempdir(), 'msdial_filtered.tsv')
write.table(cbind(filtered_info, log2(filtered_intensity + 1)), out_path,
            sep = '\t', row.names = FALSE, quote = FALSE)
cat('\nWrote filtered table to', out_path, '\n')

unlink(c(export_path, out_path))
cat('Cleaned up temp files.\n')
