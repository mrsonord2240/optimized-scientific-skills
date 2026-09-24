# test_coloc_susie_eqtl_ld_guard.R -- regression test for the eQTL LD-consistency stop.
# Inputs: path to scripts/coloc_susie.R.
# Usage:  r.sh scripts/test_coloc_susie_eqtl_ld_guard.R scripts/coloc_susie.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop('usage: test_coloc_susie_eqtl_ld_guard.R <coloc_susie.R>')
script <- normalizePath(args[[1]], mustWork = TRUE)
if (!requireNamespace('susieR', quietly = TRUE)) stop('susieR is required')

set.seed(101026)
nind <- 1800L
nsnp <- 100L
causal <- 45L
snp <- paste0('rs', seq_len(nsnp))
ld_target <- 0.75 ^ abs(outer(seq_len(nsnp), seq_len(nsnp), '-'))
gwas_genotypes <- scale(matrix(rnorm(nind * nsnp), nind, nsnp) %*% chol(ld_target))
ld <- cor(gwas_genotypes)
gwas_y <- 0.48 * gwas_genotypes[, causal] + rnorm(nind)
gwas_fits <- do.call(rbind, lapply(seq_len(nsnp), function(j) summary(lm(gwas_y ~ gwas_genotypes[, j]))$coefficients[2, 1:2]))

# These eQTL z-scores come from an independent near-identity-LD panel, intentionally
# paired with the GWAS LD matrix above. The script must reject this mismatch.
eqtl_genotypes <- scale(matrix(rnorm(nind * nsnp), nind, nsnp))
eqtl_y <- 0.48 * eqtl_genotypes[, causal] + rnorm(nind)
eqtl_fits <- do.call(rbind, lapply(seq_len(nsnp), function(j) summary(lm(eqtl_y ~ eqtl_genotypes[, j]))$coefficients[2, 1:2]))
lambda_eqtl <- susieR::estimate_s_rss(z = eqtl_fits[, 1] / eqtl_fits[, 2], R = ld, n = nind)
stopifnot(lambda_eqtl > 0.05)

work <- tempfile('coloc-susie-eqtl-ld-guard-')
dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
gwas <- data.frame(SNP = snp, CHR = 1L, POS = seq_len(nsnp) * 10000L,
                   BETA = gwas_fits[, 1], SE = gwas_fits[, 2])
eqtl <- gwas
eqtl$BETA <- eqtl_fits[, 1]
eqtl$SE <- eqtl_fits[, 2]
ld_frame <- data.frame(SNP = snp, ld, check.names = FALSE)
names(ld_frame)[-1] <- snp
gwas_path <- file.path(work, 'gwas.tsv')
eqtl_path <- file.path(work, 'eqtl.tsv')
ld_path <- file.path(work, 'ld.tsv')
out_path <- file.path(work, 'coloc')
write.table(gwas, gwas_path, sep = '\t', quote = FALSE, row.names = FALSE)
write.table(eqtl, eqtl_path, sep = '\t', quote = FALSE, row.names = FALSE)
write.table(ld_frame, ld_path, sep = '\t', quote = FALSE, row.names = FALSE)

rscript <- file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript')
captured <- file.path(work, 'stderr.txt')
status <- suppressWarnings(system2(
  rscript,
  c(script, gwas_path, eqtl_path, ld_path,
    '--gwas-type', 'quant', '--gwas-sdy', '1', '--gwas-n', nind,
    '--eqtl-type', 'quant', '--eqtl-sdy', '1', '--eqtl-n', nind,
    '--L', '10', '--out', out_path),
  stdout = captured, stderr = captured
))
message_text <- paste(readLines(captured, warn = FALSE), collapse = '\n')
stopifnot(status != 0L,
          grepl('LD reference mismatched to eQTL z-scores', message_text, fixed = TRUE),
          !file.exists(paste0(out_path, '_susie_summary.tsv')))
cat(sprintf('PASS: eQTL LD mismatch rejected (lambda=%.6f) before coloc output\n', lambda_eqtl))
quit(save = 'no', status = 0L, runLast = FALSE)
