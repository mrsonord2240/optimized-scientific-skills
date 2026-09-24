# Common-factor GWAS with Q_SNP classification (Grotzinger 2019), DWLS plus an ML cross-check.
# Tested: GenomicSEM 0.0.5 + lavaan 0.6.19 (lavaan >= 0.7.0 crashes commonfactorGWAS(); see SKILL.md).
#
# Inputs (positional):
#   1. covstruc.rds  saveRDS() of the ldsc() output list (needs S, V, I)
#   2. snps.rds|csv  the sumstats() data frame (SNP, A1, A2, MAF, beta.<trait>, se.<trait> ...)
#   3. out.tsv       per-SNP table: commonfactorGWAS columns + factor_sig, qsnp_sig,
#                    factor_only_dwls, Q_pval_ML, factor_only
#   4. cores         optional, default 1 (parallel = cores > 1)
# Usage: Rscript scripts/commonfactor_gwas_qsnp.R ldsc_results.rds sumstats_snps.rds cfgwas_qsnp.tsv 8
#
# GenomicSEM internally detects the OS via Sys.info()[['sysname']] and chooses
# PSOCK (Windows) vs FORK (Linux/Mac) clusters automatically; there is no user
# `Operating=` argument. MPI=TRUE switches to an mpirun-based strategy for
# cluster job submission. parallel=TRUE is the default.
#
# Q_pval IS the per-SNP heterogeneity test (Q_SNP in the literature). factor_only needs
# factor p < 5e-8, DWLS Q_pval and ML Q_pval both above 0.05 / (number of factor-significant SNPs).
# DWLS Q_pval is provisional: never call a SNP factor-only on it alone (SKILL.md).
library(GenomicSEM)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop('usage: commonfactor_gwas_qsnp.R covstruc.rds snps.rds|csv out.tsv [cores]')
cores <- if (length(args) >= 4) as.integer(args[4]) else 1L
ldsc_results <- readRDS(args[1])
ss <- if (endsWith(tolower(args[2]), '.rds')) readRDS(args[2]) else read.csv(args[2], stringsAsFactors = FALSE)

cfgwas <- commonfactorGWAS(
    covstruc = ldsc_results,
    SNPs = ss,
    estimation = 'DWLS',
    parallel = cores > 1,
    cores = cores,
    MPI = FALSE
)

cfgwas$factor_sig <- cfgwas$Pval_Estimate < 5e-08
cfgwas$qsnp_sig <- cfgwas$Q_pval < (0.05 / sum(cfgwas$factor_sig))
cfgwas$factor_only_dwls <- cfgwas$factor_sig & !cfgwas$qsnp_sig

cfgwas_ml <- commonfactorGWAS(covstruc = ldsc_results, SNPs = ss, estimation = 'ML', parallel = cores > 1, cores = cores)
cfgwas$Q_pval_ML <- cfgwas_ml$Q_pval[match(cfgwas$SNP, cfgwas_ml$SNP)]
cfgwas$factor_only <- cfgwas$factor_sig & !cfgwas$qsnp_sig & cfgwas$Q_pval_ML > (0.05 / sum(cfgwas$factor_sig))

# lavaan warnings are multi-line strings; flatten so the TSV keeps one row per SNP
for (col in c('fail', 'warning')) if (col %in% names(cfgwas)) cfgwas[[col]] <- gsub('[\r\n\t]+', ' ', as.character(cfgwas[[col]]))
write.table(cfgwas, args[3], sep = '\t', quote = FALSE, row.names = FALSE)
cat(sprintf('SNPs %d | factor-significant %d | factor-only (DWLS) %d | factor-only (DWLS + ML) %d\n',
            nrow(cfgwas), sum(cfgwas$factor_sig), sum(cfgwas$factor_only_dwls), sum(cfgwas$factor_only)))
