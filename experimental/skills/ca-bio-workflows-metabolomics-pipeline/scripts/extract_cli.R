args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop('usage: Rscript scripts/extract_cli.R metadata.csv positive|negative output.rds')
if (file.exists(args[3])) stop('output already exists')
pd <- read.csv(args[1], stringsAsFactors = FALSE)
stopifnot(all(c('sample', 'sample_group') %in% names(pd)), args[2] %in% c('positive', 'negative'))
mzml_files <- pd$sample
ionization_mode <- args[2]
source('scripts/stage1_xcms_extract.R')
defs$mode <- ionization_mode
saveRDS(list(feat = feat, defs = defs, pd = pd, xdata = xdata), args[3])
