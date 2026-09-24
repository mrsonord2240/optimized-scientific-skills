#!/usr/bin/env Rscript
# Focused contract regression for scripts/msstatsptm_labelfree.R.
# On the affected Windows cli 3.6.6 runtime, run with the opt-in process-local profile, for example:
# R_PROFILE_USER="$(pwd)/scripts/cli_windows_cleanup_guard.Rprofile" \
#   /f/OpenScience/audit-envs/mass-spec-proteomics-analyst/r.sh \
#   tests/test_msstatsptm_labelfree_no_global.R

this_file <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))][1])
skill_dir <- normalizePath(file.path(dirname(this_file), '..'))
analysis_script <- normalizePath(file.path(skill_dir, 'scripts', 'msstatsptm_labelfree.R'))
Sys.setenv(MSSTATSPTM_LF_NO_AUTORUN = '1')
source(analysis_script)
Sys.unsetenv('MSSTATSPTM_LF_NO_AUTORUN')

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  stopifnot(!is.null(error), grepl(pattern, error, fixed = TRUE))
}

# The effect-size threshold must fail closed before package loading or modeling.
stopifnot(identical(validate_lfc('0'), 0), identical(validate_lfc('1.5'), 1.5))
expect_error(validate_lfc('not-a-number'), 'finite, non-negative')
expect_error(validate_lfc('-0.01'), 'finite, non-negative')
expect_error(validate_lfc('Inf'), 'finite, non-negative')

root <- file.path(tempdir(), 'msstatsptm_labelfree_contract')
unlink(root, recursive = TRUE, force = TRUE)
dir.create(root, recursive = TRUE)
opt <- list(use_unmod = 'TRUE', evidence_prot = 'evidence_global.txt',
            proteinGroups = 'proteinGroups_global.txt', annotation_protein = 'annotation_protein.csv')
inp <- function(name) file.path(root, opt[[name]])

# A no-global request is valid only with the explicit proxy flag; FALSE, invalid flags, partial global
# inputs, and mixed proxy/global requests all fail before any converter can read a global file.
mode <- validate_input_mode(opt, inp)
stopifnot(identical(mode, list(use_unmod = TRUE, has_paired_global = FALSE)))
opt$use_unmod <- 'FALSE'
expect_error(validate_input_mode(opt, inp), 'No paired-global input found')
opt$use_unmod <- 'yes'
expect_error(validate_input_mode(opt, inp), 'must be exactly TRUE or FALSE')
opt$use_unmod <- 'TRUE'
file.create(inp('evidence_prot'))
expect_error(validate_input_mode(opt, inp), 'Incomplete paired-global input')
file.create(inp('proteinGroups'))
file.create(inp('annotation_protein'))
expect_error(validate_input_mode(opt, inp), 'mutually exclusive')
opt$use_unmod <- 'FALSE'
paired_a <- validate_input_mode(opt, inp)
paired_b <- validate_input_mode(opt, inp)
stopifnot(identical(paired_a, list(use_unmod = FALSE, has_paired_global = TRUE)), identical(paired_a, paired_b))

# The installed MSstatsPTM 2.8.1 converter accepts the no-global input contract with fasta_path
# (not the stale fasta alias) and returns both PTM and PROTEIN proxy inputs.
data('maxq_lf_evidence', package = 'MSstatsPTM')
data('maxq_lf_annotation', package = 'MSstatsPTM')
fasta <- system.file('extdata', 'maxq_lf_fasta.fasta', package = 'MSstatsPTM')
proxy_input <- MSstatsPTM::MaxQtoMSstatsPTMFormat(
  evidence = maxq_lf_evidence, annotation = maxq_lf_annotation, fasta_path = fasta,
  fasta_protein_name = 'uniprot_ac', mod_id = '\\(Phospho \\(STY\\)\\)',
  use_unmod_peptides = TRUE, labeling_type = 'LF', which_proteinid_ptm = 'Proteins'
)
stopifnot(identical(sort(names(proxy_input)), c('PROTEIN', 'PTM')),
          nrow(proxy_input$PTM) > 0L, nrow(proxy_input$PROTEIN) > 0L)

# Execute the literal no-global CLI with only enriched evidence, annotation, and FASTA. The bundled
# fixture intentionally fails the stricter proxy QC, but must get that far rather than reading a
# missing global file; this guards the previously unreachable documented route.
cli_dir <- file.path(root, 'literal_no_global')
cli_out <- file.path(root, 'literal_no_global_out')
dir.create(cli_dir)
write.table(maxq_lf_evidence, file.path(cli_dir, 'evidence_phospho.txt'), sep = '\t',
            row.names = FALSE, quote = FALSE)
write.csv(maxq_lf_annotation, file.path(cli_dir, 'annotation_ptm.csv'), row.names = FALSE)
file.copy(fasta, file.path(cli_dir, 'uniprot_human.fasta'))
stopifnot(!any(file.exists(file.path(cli_dir,
  c('evidence_global.txt', 'proteinGroups_global.txt', 'annotation_protein.csv')))))
cli_log <- file.path(root, 'literal_no_global.log')
cli_status <- system2(file.path(R.home('bin'), 'Rscript.exe'),
  c(analysis_script, paste0('dir=', normalizePath(cli_dir, winslash = '/')),
    paste0('out=', normalizePath(cli_out, winslash = '/', mustWork = FALSE)), 'use_unmod=TRUE'),
  stdout = cli_log, stderr = cli_log)
cli_text <- paste(readLines(cli_log, warn = FALSE), collapse = '\n')
stopifnot(cli_status != 0L, grepl('No-global proxy QC', cli_text, fixed = TRUE),
          !grepl('No paired-global input found', cli_text, fixed = TRUE),
          !grepl('evidence_global', cli_text, fixed = TRUE))

# Proxy QC is deterministic and rejects both ambiguous protein mappings and missing
# protein x Condition x BioReplicate coverage. This minimal fixture meets the declared
# two-unmodified-peptide and two-biological-replicate minima for one PTM protein.
qc_annotation <- data.frame(Raw.file = c('A1', 'A2', 'B1', 'B2'),
                            Condition = c('A', 'A', 'B', 'B'),
                            BioReplicate = c(1, 2, 1, 2))
qc_evidence <- data.frame(
  Sequence = rep(c('MOD1', 'UNMOD1', 'UNMOD2'), 4),
  Proteins = rep('P1', 12), Raw.file = rep(qc_annotation$Raw.file, each = 3),
  Modified.sequence = rep(c('PEPTIDE(Phospho (STY))', 'UNMODIFIED', 'UNMODIFIED'), 4),
  stringsAsFactors = FALSE
)
is_mod <- grepl('Phospho \\(STY\\)', qc_evidence$Modified.sequence)
qc_one <- proxy_qc(qc_evidence, qc_annotation, is_mod)
qc_two <- proxy_qc(qc_evidence, qc_annotation, is_mod)
stopifnot(identical(qc_one, qc_two), nrow(qc_one) > 0L,
          all(qc_one$unique_unmodified_peptides > 0L))
ambiguous <- qc_evidence
first_unmod <- which(!is_mod)[1]
ambiguous$Proteins[first_unmod] <- paste0(ambiguous$Proteins[first_unmod], ';P99999')
expect_error(proxy_qc(ambiguous, qc_annotation, is_mod), 'unambiguous single-protein')
missing_coverage <- qc_evidence[is_mod | qc_evidence$Raw.file != 'A1', ]
expect_error(proxy_qc(missing_coverage, qc_annotation,
                      grepl('Phospho \\(STY\\)', missing_coverage$Modified.sequence)),
             'every PTM protein/Condition/BioReplicate')
unannotated_raw <- qc_evidence
unannotated_raw$Raw.file[1] <- 'UNANNOTATED'
expect_error(proxy_qc(unannotated_raw, qc_annotation,
                      grepl('Phospho \\(STY\\)', unannotated_raw$Modified.sequence)),
             'absent from annotation')
conflicting_annotation <- rbind(qc_annotation, data.frame(Raw.file = 'A1', Condition = 'B', BioReplicate = 9))
expect_error(proxy_qc(qc_evidence, conflicting_annotation, is_mod),
             'one-to-one Raw.file to Condition/BioReplicate')

# Publication is deterministic and refuses both an existing file and a dangling symlink.
out <- file.path(root, 'out')
dir.create(out)
destination <- file.path(out, 'proxy_adjusted_sites.csv')
stage_dir <- create_output_staging_dir(out)
staged <- stage_csv_no_clobber(data.frame(value = 1), destination, stage_dir)
publish_staged_no_clobber(list(staged))
stopifnot(identical(read.csv(destination)$value, 1L))
expect_error(stage_csv_no_clobber(data.frame(value = 2), destination), 'Refusing to overwrite')
link_path <- file.path(out, 'dangling.csv')
# Windows developer mode is not enabled in the audit runtime, so model a dangling link at the
# filesystem seam: file.exists is false while Sys.readlink is non-empty.
original_readlink <- Sys.readlink
assign('Sys.readlink', function(path) if (identical(path, link_path)) 'missing-target.csv' else original_readlink(path),
       envir = .GlobalEnv)
on.exit(assign('Sys.readlink', original_readlink, envir = .GlobalEnv), add = TRUE)
stopifnot(!file.exists(link_path), nzchar(Sys.readlink(link_path)))
expect_error(stage_csv_no_clobber(data.frame(value = 3), link_path), 'Refusing to overwrite')

# The manifest has exactly three artifact/hash rows and is linked last, after all three artifacts.
manifest_out <- file.path(root, 'manifest_out')
dir.create(manifest_out)
manifest_stage_dir <- create_output_staging_dir(manifest_out)
artifact_names <- c('proxy_adjusted_sites.csv', 'proxy_ptm_model.csv', 'proxy_adjusted_qc.csv')
artifact_staged <- lapply(seq_along(artifact_names), function(i) {
  stage_csv_no_clobber(data.frame(value = i), file.path(manifest_out, artifact_names[i]), manifest_stage_dir)
})
manifest <- make_proxy_manifest(artifact_staged)
stopifnot(identical(manifest$schema_version, rep(1L, 3L)),
          identical(manifest$mode, rep('no-global-proxy', 3L)),
          identical(manifest$complete, rep(TRUE, 3L)),
          identical(manifest$artifact, artifact_names),
          identical(manifest$bytes, file.info(vapply(artifact_staged, `[[`, character(1), 'stage'))$size),
          identical(manifest$md5, unname(tools::md5sum(vapply(artifact_staged, `[[`, character(1), 'stage')))))
marker <- stage_csv_no_clobber(manifest, file.path(manifest_out, 'proxy_adjusted_manifest.csv'), manifest_stage_dir)
link_order <- character()
recording_link <- function(from, to) { link_order <<- c(link_order, basename(to)); file.link(from, to) }
publish_staged_no_clobber(c(artifact_staged, list(marker)), link_fn = recording_link)
published_manifest <- read.csv(file.path(manifest_out, 'proxy_adjusted_manifest.csv'), check.names = FALSE,
                               stringsAsFactors = FALSE)
stopifnot(identical(link_order, c(artifact_names, 'proxy_adjusted_manifest.csv')),
          identical(names(published_manifest), names(manifest)),
          identical(as.integer(published_manifest$schema_version), manifest$schema_version),
          identical(as.character(published_manifest$mode), manifest$mode),
          identical(as.logical(published_manifest$complete), manifest$complete),
          identical(as.character(published_manifest$artifact), manifest$artifact),
          identical(as.numeric(published_manifest$bytes), as.numeric(manifest$bytes)),
          identical(as.character(published_manifest$md5), manifest$md5))

# A simulated marker-link failure rolls back the three artifacts and leaves no completion marker.
marker_fail_out <- file.path(root, 'marker_fail_out')
dir.create(marker_fail_out)
marker_fail_stage <- create_output_staging_dir(marker_fail_out)
marker_fail_artifacts <- lapply(seq_along(artifact_names), function(i) {
  stage_csv_no_clobber(data.frame(value = i), file.path(marker_fail_out, artifact_names[i]), marker_fail_stage)
})
marker_fail <- stage_csv_no_clobber(make_proxy_manifest(marker_fail_artifacts),
                                   file.path(marker_fail_out, 'proxy_adjusted_manifest.csv'), marker_fail_stage)
fail_marker_link <- function(from, to) {
  if (identical(basename(to), 'proxy_adjusted_manifest.csv')) return(FALSE)
  file.link(from, to)
}
expect_error(publish_staged_no_clobber(c(marker_fail_artifacts, list(marker_fail)), link_fn = fail_marker_link),
             'rolled back this invocation')
stopifnot(!any(file.exists(file.path(marker_fail_out, c(artifact_names, 'proxy_adjusted_manifest.csv')))))

# A simulated second atomic-link failure rolls back only the first file published by this invocation.
rollback_dir <- create_output_staging_dir(out)
first_destination <- file.path(out, 'first.csv')
second_destination <- file.path(out, 'second.csv')
rollback_staged <- list(stage_csv_no_clobber(data.frame(value = 1), first_destination, rollback_dir),
                        stage_csv_no_clobber(data.frame(value = 2), second_destination, rollback_dir))
link_calls <- 0L
fail_second_link <- function(from, to) {
  link_calls <<- link_calls + 1L
  if (link_calls == 2L) return(FALSE)
  file.link(from, to)
}
expect_error(publish_staged_no_clobber(rollback_staged, link_fn = fail_second_link), 'rolled back this invocation')
stopifnot(!file.exists(first_destination), !file.exists(second_destination))

message('PASS no-global proxy contract, QC, and no-clobber publication')
