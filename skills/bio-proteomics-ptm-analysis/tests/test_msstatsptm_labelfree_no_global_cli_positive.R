#!/usr/bin/env Rscript
# On the affected Windows cli 3.6.6 runtime, invoke this test with
# R_PROFILE_USER=/absolute/path/to/scripts/cli_windows_cleanup_guard.Rprofile.
# Positive literal-CLI regression. It intentionally requires natural exit 0; do not weaken this
# assertion to accommodate a host/runtime teardown fault.
# Run with an R 4.4 runtime carrying MSstatsPTM 2.8.1 and its private library:
# Rscript tests/test_msstatsptm_labelfree_no_global_cli_positive.R

this_file <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))][1])
skill_dir <- normalizePath(file.path(dirname(this_file), '..'))
analysis_script <- normalizePath(file.path(skill_dir, 'scripts', 'msstatsptm_labelfree.R'))
runtime <- Sys.getenv('PTM_RSCRIPT', file.path(R.home('bin'), 'Rscript.exe'))
if (!file.exists(runtime)) stop('Set PTM_RSCRIPT to the isolated Rscript executable.')

# Build a deterministic valid proxy fixture from MSstatsPTM's real MaxQuant/FASTA-mappable P09938 rows.
data('maxq_lf_evidence', package = 'MSstatsPTM')
fasta <- system.file('extdata', 'maxq_lf_fasta.fasta', package = 'MSstatsPTM')
source_evidence <- maxq_lf_evidence[maxq_lf_evidence$Proteins == 'P09938', ]
source_mod <- source_evidence[grepl('Phospho \\(STY\\)', source_evidence$Modified.sequence), ]
source_mod <- source_mod[!duplicated(source_mod$Modified.sequence), ]
source_unmod <- source_evidence[!grepl('Phospho \\(STY\\)', source_evidence$Modified.sequence), ]
source_unmod <- source_unmod[!duplicated(source_unmod$Sequence), ][1:2, ]
stopifnot(nrow(source_mod) >= 2L, nrow(source_unmod) == 2L)

annotation <- data.frame(
  Raw.file = c(paste0('Control_', 1:4), paste0('Treatment_', 1:4)),
  Condition = c(rep('Control', 4), rep('Treatment', 4)),
  BioReplicate = rep(1:4, 2)
)
# Four biological replicates per condition and a feature-by-run interaction avoid the degenerate
# all-proportional synthetic design that leaves lmer without estimable comparison columns.
run_factor <- c(0.93, 1.01, 1.08, 1.16, 1.39, 1.51, 1.63, 1.77)
synthetic <- do.call(rbind, lapply(seq_len(nrow(annotation)), function(i) {
  rows <- rbind(source_mod, source_unmod)
  rows$Raw.file <- annotation$Raw.file[i]
  rows$Experiment <- annotation$Raw.file[i]
  feature_factor <- 1 + 0.012 * seq_len(nrow(rows)) * ((i %% 3L) - 1) +
    0.004 * (seq_len(nrow(rows)) %% 2L) * i
  rows$Intensity <- as.numeric(rows$Intensity) * run_factor[i] * feature_factor
  rows
}))

test_root <- Sys.getenv('PTM_TEST_ROOT', tempdir())
if (!dir.exists(test_root)) dir.create(test_root, recursive = TRUE)
root <- file.path(normalizePath(test_root, winslash = '/', mustWork = TRUE), 'msstatsptm_no_global_positive')
input_dir <- file.path(root, 'input')
out_dir <- file.path(root, 'out')
unlink(root, recursive = TRUE, force = TRUE)
dir.create(input_dir, recursive = TRUE)
write.table(synthetic, file.path(input_dir, 'evidence_phospho.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(annotation, file.path(input_dir, 'annotation_ptm.csv'), row.names = FALSE)
file.copy(fasta, file.path(input_dir, 'uniprot_human.fasta'))
stopifnot(!any(file.exists(file.path(input_dir,
  c('evidence_global.txt', 'proteinGroups_global.txt', 'annotation_protein.csv')))))

log_file <- file.path(root, 'literal_no_global.log')
status <- system2(runtime,
  c(analysis_script, paste0('dir=', normalizePath(input_dir, winslash = '/')),
    paste0('out=', normalizePath(out_dir, winslash = '/', mustWork = FALSE)), 'use_unmod=TRUE'),
  stdout = log_file, stderr = log_file)
log_text <- paste(readLines(log_file, warn = FALSE), collapse = '\n')
if (status != 0L) stop('Positive no-global literal CLI must exit 0; got ', status, '.\n', log_text)

artifacts <- c('proxy_adjusted_sites.csv', 'proxy_ptm_model.csv', 'proxy_adjusted_qc.csv')
manifest_path <- file.path(out_dir, 'proxy_adjusted_manifest.csv')
stopifnot(all(file.exists(file.path(out_dir, artifacts))), file.exists(manifest_path))
proxy <- read.csv(file.path(out_dir, 'proxy_adjusted_sites.csv'), check.names = FALSE)
qc <- read.csv(file.path(out_dir, 'proxy_adjusted_qc.csv'), check.names = FALSE)
manifest <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c('raw_ptm_log2FC', 'raw_ptm_adj.pvalue', 'adjustment_source', 'interpretation') %in% names(proxy)),
          all(proxy$adjustment_source == 'co-enriched unmodified peptides; no paired global proteome'),
          all(proxy$interpretation == 'proxy-adjusted candidate; not a regulation call'),
          all(c('Protein', 'Condition', 'BioReplicate', 'unique_unmodified_peptides') %in% names(qc)),
          identical(manifest$schema_version, rep(1L, 3L)),
          identical(manifest$mode, rep('no-global-proxy', 3L)),
          identical(manifest$complete, rep(TRUE, 3L)),
          identical(manifest$artifact, artifacts),
          identical(as.numeric(manifest$bytes), as.numeric(file.info(file.path(out_dir, artifacts))$size)),
          identical(manifest$md5, unname(tools::md5sum(file.path(out_dir, artifacts)))))
message('PASS positive no-global literal CLI: ', out_dir)
