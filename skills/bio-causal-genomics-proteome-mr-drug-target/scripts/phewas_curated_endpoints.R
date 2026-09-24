# Phenome-wide cis-MR for one drug target over a curated endpoint list (on-target adverse-effect scan).
# Inputs:  cis-pQTL instrument table (TSV; columns SNP BETA SE A1 A2 EAF P), curated endpoint table
#          (TSV; column id = OpenGWAS study ids, e.g. the FinnGen DF12 endpoints).
# Output:  all IVW results and a Bonferroni-significant table (TSV), plus the top hits printed.
# Needs:   TwoSampleMR, ieugwasr; OpenGWAS access (JWT token) for available_outcomes()/extract_outcome_data().
# Usage:   Rscript scripts/phewas_curated_endpoints.R <pqtl.tsv> <endpoints.tsv> [out.tsv] [min_sample_size=50000] [population=European]
# Example: Rscript scripts/phewas_curated_endpoints.R pcsk9_cis_pqtls.tsv finngen_DF12_endpoints.tsv pcsk9_phewas_mr.tsv
library(TwoSampleMR); library(ieugwasr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('usage: phewas_curated_endpoints.R <pqtl.tsv> <endpoints.tsv> [out.tsv] [min_sample_size] [population]')
pqtl_file <- args[1]
endpoints_file <- args[2]
out_file <- if (length(args) >= 3) args[3] else 'phewas_mr.tsv'
min_sample_size <- if (length(args) >= 4) as.numeric(args[4]) else 50000
pop_filter <- if (length(args) >= 5) args[5] else 'European'

target_pqtl <- read.table(pqtl_file, header = TRUE)
exposure_dat <- format_data(target_pqtl, type = 'exposure',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2', eaf_col = 'EAF', pval_col = 'P')

curated_endpoints <- read.table(endpoints_file, header = TRUE)  # ~3000 curated endpoints from finngen.fi
outcomes <- available_outcomes()
outcomes_filt <- subset(outcomes, id %in% curated_endpoints$id & sample_size >= min_sample_size & population == pop_filter)

results <- lapply(outcomes_filt$id, function(out_id) {
    outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = out_id)
    if (is.null(outcome_dat) || nrow(outcome_dat) < 2) return(NULL)
    dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    mr(dat, method_list = 'mr_ivw')
})

results_df <- do.call(rbind, Filter(Negate(is.null), results))
n_tests <- nrow(curated_endpoints)
results_df$p_bonf <- pmin(results_df$pval * n_tests, 1)
top_hits <- subset(results_df, pval < 0.05 / n_tests)   # 0.05 / 3000 = 1.7e-5

cat('Outcomes scanned:', length(results), ' with results:', nrow(results_df), ' Bonferroni-significant:', nrow(top_hits), '\n')
print(top_hits)
write.table(results_df, file = out_file, sep = '\t', row.names = FALSE, quote = FALSE)
