args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('usage: Rscript scripts/msdial_cli.R alignment.mdalign output.rds')
if (file.exists(args[2])) stop('output already exists')
export_file <- args[1]
source('scripts/msdial_import.R')
saveRDS(list(feat=feat, defs=defs, sample_class=sample_class, injection_order=injection_order, batch_id=batch_id), args[2])
