# MS-DIAL alignment export -> Stage 2 inputs (feat, defs, sample_class, injection_order, batch_id).
# Input (define before source(), else the first AlignResult-*.mdalign in the working directory is used):
#   export_file (path to an MSDIALCUI AlignResult-<timestamp>.mdalign file).
# Outputs (left in the calling environment): msdial, feat, defs, sample_class, injection_order, batch_id.
# Usage: source('scripts/msdial_import.R')   then run Stage 2 from `fm <- feat`
# Checked on a real MSDIALCUI 5.5.260820 export (R 4.4.3).
if (!exists('export_file')) export_file <- Sys.glob('AlignResult-*.mdalign')[1]
stopifnot(!is.na(export_file), file.exists(export_file))
hdr <- strsplit(readLines(export_file, n = 4), '\t', fixed = TRUE)   # class / file type / order / batch
msdial <- read.csv(export_file, sep = '\t', skip = 4, check.names = FALSE)
s_idx <- (which(hdr[[1]] == 'Class') + 1):length(hdr[[1]])          # sample columns follow the "Class" cell
file_type <- hdr[[2]][s_idx]
use <- file_type %in% c('Sample', 'QC')      # Blank/Standard injections leave the matrix (blank filter: normalization-qc)

feat <- as.matrix(msdial[, colnames(msdial)[s_idx][use]]); storage.mode(feat) <- 'numeric'
feat[feat == 0] <- NA                        # MS-DIAL writes not-detected as 0; Stage 2's detection filter counts NA
rownames(feat) <- paste0('FT', msdial[['Alignment ID']])
defs <- data.frame(mzmed = msdial[['Average Mz']], rtmed = msdial[['Average Rt(min)']] * 60,  # seconds, as xcms
                   row.names = rownames(feat))
sample_class <- ifelse(file_type[use] == 'QC', 'QC', hdr[[1]][s_idx][use])   # pooled QCs must be labelled 'QC'
injection_order <- as.integer(hdr[[3]][s_idx][use])
batch_id <- as.integer(hdr[[4]][s_idx][use])
stopifnot(ncol(feat) == length(sample_class), !anyNA(injection_order), !anyNA(batch_id))
