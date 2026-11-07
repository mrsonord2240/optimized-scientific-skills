# Purpose : TwoSampleMR standard workflow. Genome-wide-significant instruments -> F filter -> LD clump
#           -> outcome -> harmonise -> IVW/Egger/weighted median/mode -> Cochran Q, Egger intercept,
#           leave-one-out -> Steiger directionality. Writes the tables and the harmonised data
#           (harmonised.tsv) that mr_presso_outliers.R and simex_egger.R read.
# Inputs  : two tab-separated GWAS files with columns SNP, BETA, SE, A1 (effect allele), A2, EAF, P, N
#           (rename via the col_* variables below if yours differ).
# Usage   : r.sh scripts/twosample_workflow.R --exposure exposure_gwas.tsv --outcome outcome_gwas.tsv \
#               --outdir mr_out --bfile 1kg_EUR/EUR [--plink /path/to/plink] [--no-clump]
#           --bfile: LD reference (1KG EUR or matched ancestry). --plink defaults to
#           genetics.binaRies::get_plink_binary(). --no-clump skips clumping: use ONLY for instruments
#           that are already LD-independent (e.g. synthetic data); the default is to clump.
# Output  : <outdir>/{primary,heterogeneity,pleiotropy,leaveoneout,steiger,harmonised}.tsv

suppressMessages({library(TwoSampleMR); library(ieugwasr)})

# ---- inputs -------------------------------------------------------------------------------
args <- commandArgs(TRUE)
opt <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
exposure_file <- opt('--exposure'); outcome_file <- opt('--outcome'); outdir <- opt('--outdir', 'mr_out')
bfile <- opt('--bfile'); plink_bin <- opt('--plink'); do_clump <- !('--no-clump' %in% args)
if (is.null(exposure_file) || is.null(outcome_file)) stop('usage: --exposure FILE --outcome FILE [--outdir DIR] [--bfile PREFIX] [--plink BIN] [--no-clump]')
if (do_clump && is.null(bfile)) stop('LD clumping needs --bfile <1KG reference prefix> (or pass --no-clump for pre-clumped instruments)')
col_snp <- 'SNP'; col_beta <- 'BETA'; col_se <- 'SE'; col_a1 <- 'A1'; col_a2 <- 'A2'
col_eaf <- 'EAF'; col_p <- 'P'; col_n <- 'N'
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

exposure_raw <- read_exposure_data(
    filename = exposure_file, sep = '\t',
    snp_col = col_snp, beta_col = col_beta, se_col = col_se,
    effect_allele_col = col_a1, other_allele_col = col_a2,
    eaf_col = col_eaf, pval_col = col_p, samplesize_col = col_n  # required for directionality_test() below
)

exposure_sig <- subset(exposure_raw, pval.exposure < 5e-08)  # genome-wide significance

# F-statistic computed from EXPOSURE (Burgess 2011); ratio of squared effect to its variance
exposure_sig$f_stat <- (exposure_sig$beta.exposure / exposure_sig$se.exposure)^2
exposure_sig <- subset(exposure_sig, f_stat >= 10)  # Staiger-Stock 1997 weak-IV heuristic

if (do_clump) {
    clumped <- ld_clump(
        data.frame(rsid = exposure_sig$SNP, pval = exposure_sig$pval.exposure),
        clump_r2 = 0.001, clump_kb = 10000,  # polygenic MR convention
        plink_bin = if (is.null(plink_bin)) genetics.binaRies::get_plink_binary() else plink_bin,
        bfile = bfile
    )
    exposure_dat <- subset(exposure_sig, SNP %in% clumped$rsid)
} else {
    message('--no-clump: instruments assumed LD-independent')
    exposure_dat <- exposure_sig
}

outcome_dat <- read_outcome_data(
    filename = outcome_file, snps = exposure_dat$SNP, sep = '\t',
    snp_col = col_snp, beta_col = col_beta, se_col = col_se,
    effect_allele_col = col_a1, other_allele_col = col_a2,
    eaf_col = col_eaf, pval_col = col_p, samplesize_col = col_n  # required for directionality_test() below
)

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)  # infer from EAF; drops MAF~0.5 palindromes

primary <- mr(dat, method_list = c('mr_ivw', 'mr_egger_regression',
                                    'mr_weighted_median', 'mr_weighted_mode'))

heterogeneity <- mr_heterogeneity(dat)         # Cochran Q
pleiotropy <- mr_pleiotropy_test(dat)          # Egger intercept
loo <- mr_leaveoneout(dat)                     # influential-SNP check

steiger <- directionality_test(dat)            # variance-explained direction
if (is.null(steiger)) {
    stop("directionality_test() returned NULL -- dat is missing pval.exposure/pval.outcome/",
         "samplesize.exposure/samplesize.outcome (or supply pre-computed r.exposure/r.outcome, ",
         "e.g. via get_r_from_lor() for binary traits). It fails silently, not loudly, so check ",
         "for NULL rather than trusting a downstream NULL$correct_causal_direction.")
}

# ---- outputs ------------------------------------------------------------------------------
wt <- function(x, f) write.table(x, file.path(outdir, f), sep = '\t', row.names = FALSE, quote = FALSE)
wt(primary, 'primary.tsv'); wt(heterogeneity, 'heterogeneity.tsv'); wt(pleiotropy, 'pleiotropy.tsv')
wt(loo, 'leaveoneout.tsv'); wt(steiger, 'steiger.tsv'); wt(dat, 'harmonised.tsv')
cat('instruments after F/clump:', nrow(exposure_dat), ' harmonised:', nrow(dat),
    ' mr_keep:', sum(dat$mr_keep), '\n')
print(primary[, c('method', 'nsnp', 'b', 'se', 'pval')])
print(steiger[, c('correct_causal_direction', 'steiger_pval')])
