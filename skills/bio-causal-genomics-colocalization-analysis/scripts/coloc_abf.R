# coloc_abf.R -- single-causal coloc.abf for one locus, with the MHC / chr 8 gate and p12 sensitivity.
#
# Purpose: test whether a GWAS signal and an eQTL share a causal variant at one locus.
# Inputs:  two tab-separated files with the SAME SNP set and allele coding (run scripts/harmonise.R
#          first). Columns: SNP, CHR, POS, BETA, SE (hg38 positions; MAF also needed if sdY is estimated).
#          --gwas-type cc|quant  --gwas-s <case fraction, cc>  --gwas-sdy <quant>  --gwas-n <N>
#          --eqtl-type cc|quant  --eqtl-s ...                  --eqtl-sdy <1 if standardised> --eqtl-n <N>
#          --p1 1e-4 --p2 1e-4 --p12 5e-6   --out <prefix>
# Output:  <prefix>_summary.tsv (PP.H0-H4), <prefix>_snps.tsv (per-SNP PP.H4), <prefix>_sensitivity.pdf;
#          PP table on stdout.
# Usage:   Rscript scripts/coloc_abf.R gwas.tsv eqtl.tsv --gwas-type cc --gwas-s 0.30 --gwas-n 50000 \
#              --eqtl-type quant --eqtl-sdy 1 --eqtl-n 500 --p12 5e-6 --out coloc_out

library(coloc)

script_dir <- dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, 'flag_excluded_region.R'))

argv <- commandArgs(TRUE)
files <- argv[!grepl('^--', argv) & !c(FALSE, grepl('^--', head(argv, -1)))]
opt <- list(`gwas-type` = 'cc', `gwas-s` = '0.30', `gwas-sdy` = NA, `gwas-n` = NA,
            `eqtl-type` = 'quant', `eqtl-s` = NA, `eqtl-sdy` = '1', `eqtl-n` = NA,
            p1 = '1e-4', p2 = '1e-4', p12 = '5e-6', out = 'coloc_out')
for (i in which(grepl('^--', argv))) opt[[sub('^--', '', argv[i])]] <- argv[i + 1]
stopifnot(length(files) == 2, !is.na(opt$`gwas-n`), !is.na(opt$`eqtl-n`))
num <- function(x) if (is.na(x)) NULL else as.numeric(x)

gwas_df <- read.delim(files[1], stringsAsFactors = FALSE)
eqtl_df <- read.delim(files[2], stringsAsFactors = FALSE)
stopifnot(identical(gwas_df$SNP, eqtl_df$SNP))   # harmonised, same SNP set and order

# Programmatic MHC / chr 8 inversion gate -- run BEFORE coloc.abf
lead <- which.min(if ('P' %in% names(gwas_df)) gwas_df$P else -abs(gwas_df$BETA / gwas_df$SE))
stop_if_excluded_region(chr = gwas_df$CHR[1], pos_bp = gwas_df$POS[lead])

dataset <- function(df, pfx) {
  d <- list(beta = df$BETA, varbeta = df$SE^2, snp = df$SNP, position = df$POS,
            type = opt[[paste0(pfx, '-type')]], N = num(opt[[paste0(pfx, '-n')]]))
  if (d$type == 'cc') d$s <- num(opt[[paste0(pfx, '-s')]])                 # case fraction
  else if (!is.na(opt[[paste0(pfx, '-sdy')]])) d$sdY <- num(opt[[paste0(pfx, '-sdy')]])
  d
}
gwas_input <- dataset(gwas_df, 'gwas')
eqtl_input <- dataset(eqtl_df, 'eqtl')

res <- coloc.abf(dataset1 = gwas_input, dataset2 = eqtl_input,
                 p1 = num(opt$p1), p2 = num(opt$p2), p12 = num(opt$p12))

print(res$summary)
write.table(t(res$summary), paste0(opt$out, '_summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(res$results[order(-res$results$SNP.PP.H4), c('snp', 'SNP.PP.H4')],
            paste0(opt$out, '_snps.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
pdf(paste0(opt$out, '_sensitivity.pdf'))
sens <- coloc::sensitivity(res, rule = 'H4 > 0.75')   # generates plot + table
invisible(dev.off())
