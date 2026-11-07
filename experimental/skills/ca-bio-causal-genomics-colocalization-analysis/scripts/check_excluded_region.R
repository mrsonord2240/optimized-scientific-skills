#!/usr/bin/env Rscript
# Usage: Rscript scripts/check_excluded_region.R chromosome hg38_position result.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop('usage: check_excluded_region.R chromosome hg38_position result.tsv')
script_dir <- dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, 'flag_excluded_region.R'))
pos <- suppressWarnings(as.numeric(args[2]))
if (!is.finite(pos) || pos < 1 || pos != floor(pos)) stop('position must be a positive integer in hg38')
flag <- flag_excluded_region(args[1], pos, build = 'hg38')
result <- data.frame(chromosome = args[1], position_hg38 = pos,
                     exclusion = if (is.na(flag)) 'none' else flag,
                     action = if (is.na(flag)) 'outside_named_exclusions' else 'stop_and_use_specialized_method_or_exclude')
write.table(result, args[3], sep = '\t', quote = FALSE, row.names = FALSE)
print(result)
