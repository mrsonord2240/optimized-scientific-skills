# coloc_susie.R -- multi-causal colocalization (SuSiE per trait, then coloc per credible-set pair).
#
# Purpose: test colocalization at a locus with several independent signals (allelic heterogeneity).
# Inputs:  gwas.tsv, eqtl.tsv: columns SNP, CHR, POS, BETA, SE on the SAME SNP set and order (harmonised).
#          ld.tsv: signed Pearson r matrix, tab-separated, first column = SNP ids, header = SNP ids,
#          in the SAME order as the beta vectors and ancestry-matched (see references/coloc-susie.md).
#          --gwas-type cc|quant --gwas-s <case fraction> --gwas-sdy <quant> --gwas-n <N>
#          --eqtl-type cc|quant --eqtl-s ...             --eqtl-sdy <1 if standardised> --eqtl-n <N>
#          --L 10  --lambda-max 0.05  --out <prefix>
# Output:  <prefix>_susie_summary.tsv (one row per credible-set pair: hit1, hit2, PP.H0-H4); stdout summary.
#          Stops if the locus is in the MHC / chr 8 inversion, or if estimate_s_rss lambda > --lambda-max.
# Usage:   Rscript scripts/coloc_susie.R gwas.tsv eqtl.tsv ld.tsv --gwas-type cc --gwas-s 0.3 --gwas-n 50000 \
#              --eqtl-type quant --eqtl-sdy 1 --eqtl-n 500 --L 10 --out coloc_susie_out

library(coloc); library(susieR)

script_dir <- dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, 'flag_excluded_region.R'))

argv <- commandArgs(TRUE)
files <- argv[!grepl('^--', argv) & !c(FALSE, grepl('^--', head(argv, -1)))]
opt <- list(`gwas-type` = 'cc', `gwas-s` = '0.3', `gwas-sdy` = NA, `gwas-n` = NA,
            `eqtl-type` = 'quant', `eqtl-s` = NA, `eqtl-sdy` = '1', `eqtl-n` = NA,
            L = '10', `lambda-max` = '0.05', out = 'coloc_susie_out')
for (i in which(grepl('^--', argv))) opt[[sub('^--', '', argv[i])]] <- argv[i + 1]
stopifnot(length(files) == 3, !is.na(opt$`gwas-n`), !is.na(opt$`eqtl-n`))
num <- function(x) if (is.na(x)) NULL else as.numeric(x)

gwas_df <- read.delim(files[1], stringsAsFactors = FALSE)
eqtl_df <- read.delim(files[2], stringsAsFactors = FALSE)
ld_matrix <- as.matrix(read.delim(files[3], row.names = 1, check.names = FALSE))
gwas_n <- num(opt$`gwas-n`)

# Same MHC / chr 8 inversion gate as the coloc.abf pipeline
lead <- which.min(if ('P' %in% names(gwas_df)) gwas_df$P else -abs(gwas_df$BETA / gwas_df$SE))
stop_if_excluded_region(chr = gwas_df$CHR[1], pos_bp = gwas_df$POS[lead])

# LD matrix MUST be in the same SNP order as the beta vector; mis-ordering silently produces nonsense
stopifnot(identical(rownames(ld_matrix), gwas_df$SNP), identical(colnames(ld_matrix), gwas_df$SNP),
          identical(gwas_df$SNP, eqtl_df$SNP))

# Diagnostic: z-score vs LD consistency MUST be checked independently for each trait.
z_gwas <- gwas_df$BETA / gwas_df$SE
lam_gwas <- susieR::estimate_s_rss(z=z_gwas, R=ld_matrix, n=gwas_n)
if (lam_gwas > num(opt$`lambda-max`)) stop('LD reference mismatched to z-scores; lambda=', lam_gwas)
z_eqtl <- eqtl_df$BETA / eqtl_df$SE
lam_eqtl <- susieR::estimate_s_rss(z=z_eqtl, R=ld_matrix, n=num(opt$`eqtl-n`))
if (lam_eqtl > num(opt$`lambda-max`)) stop('LD reference mismatched to eQTL z-scores; lambda=', lam_eqtl)

dataset <- function(df, pfx) {
  d <- list(beta = df$BETA, varbeta = df$SE^2, snp = df$SNP, position = df$POS,
            type = opt[[paste0(pfx, '-type')]], N = num(opt[[paste0(pfx, '-n')]]), LD = ld_matrix)
  if (d$type == 'cc') d$s <- num(opt[[paste0(pfx, '-s')]])
  else if (!is.na(opt[[paste0(pfx, '-sdy')]])) d$sdY <- num(opt[[paste0(pfx, '-sdy')]])
  d
}
s1 <- runsusie(dataset(gwas_df, 'gwas'), L = num(opt$L))
s2 <- runsusie(dataset(eqtl_df, 'eqtl'), L = num(opt$L))

res_susie <- coloc.susie(s1, s2)   # NULL if no overlapping CS
# res_susie$summary rows: each (hit1, hit2) pair of credible sets
if (is.null(res_susie$summary)) {
  cat('No overlapping credible sets between the traits\n')
} else {
  print(res_susie$summary[, c('hit1', 'hit2', 'PP.H3.abf', 'PP.H4.abf')])
  write.table(res_susie$summary, paste0(opt$out, '_susie_summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
}
