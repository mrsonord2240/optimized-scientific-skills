# Purpose : MR-PRESSO global / outlier / distortion tests and the outlier-SNP list, using the
#           corrected rule (adjusted P <= SignifThreshold; MRPRESSO's Outlier Test P is already
#           Bonferroni-adjusted, so it is NOT divided by nrow(dat) again).
# Inputs  : harmonised.tsv from twosample_workflow.R (columns beta.exposure, se.exposure,
#           beta.outcome, se.outcome, SNP, mr_keep).
# Usage   : r.sh scripts/mr_presso_outliers.R --dat mr_out/harmonised.tsv [--nb 10000] [--seed 42]
#               [--signif 0.05] [--out mr_out/presso_outliers.txt]
#           --nb must exceed nrow(dat)/signif (100 SNPs need > 2000); 10000 draws on ~100 SNPs can
#           exceed 40 minutes on a shared CPU, so run it in the background and use 3000-5000 to explore.
# Output  : prints Global and Distortion tests; writes one outlier SNP per line to --out.

suppressMessages(library(MRPRESSO))

args <- commandArgs(TRUE)
opt <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
dat_file <- opt('--dat'); nb <- as.integer(opt('--nb', '10000')); seed <- as.integer(opt('--seed', '42'))
signif_threshold <- as.numeric(opt('--signif', '0.05')); out_file <- opt('--out', 'presso_outliers.txt')
if (is.null(dat_file)) stop('usage: --dat harmonised.tsv [--nb N] [--seed S] [--signif 0.05] [--out FILE]')
dat <- read.delim(dat_file, stringsAsFactors = FALSE)

dat_p <- dat[dat$mr_keep, ]  # harmonise_data() keeps dropped palindromes as mr_keep = FALSE rows; MR-PRESSO would fit them
set.seed(seed)  # mr_presso()'s global/outlier tests are Monte-Carlo; seed for a reproducible p-value
presso <- mr_presso(
    BetaOutcome = 'beta.outcome', BetaExposure = 'beta.exposure',
    SdOutcome = 'se.outcome', SdExposure = 'se.exposure',
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
    data = dat_p, NbDistribution = nb,  # >= 10000 for publication-grade p-value precision
    SignifThreshold = signif_threshold
)

print(presso$`MR-PRESSO results`$`Global Test`)         # any pleiotropy
print(presso$`MR-PRESSO results`$`Distortion Test`)     # change after outlier removal

# Outlier SNPs. `Outlier Test` is NULL unless the global test is significant. Its Pvalue is
# ALREADY Bonferroni-adjusted inside MRPRESSO (raw p x nrow(dat)) and may be a string such as
# "<3e-04", so do not divide the threshold by nrow(dat) again (that double correction flagged 0
# outliers where MRPRESSO's own rule flagged 11, 9 of them planted; checked on MRPRESSO 1.0).
# MRPRESSO's own rule is adjusted P <= SignifThreshold.
ot <- presso$`MR-PRESSO results`$`Outlier Test`
outlier_snps <- character(0)
if (!is.null(ot)) {
    p_adj <- suppressWarnings(as.numeric(sub('^<', '', ot$Pvalue)))
    outlier_snps <- dat_p[rownames(ot)[which(p_adj <= signif_threshold)], 'SNP']
}
writeLines(outlier_snps, out_file)
cat('outlier SNPs (', length(outlier_snps), ') written to ', out_file, '\n', sep = '')
