#!/usr/bin/env Rscript
# Protein-adjusted phosphosite testing with MSstatsPTM, label-free (DDA/DIA) route.
#
# Class-I pre-filter on the enriched MaxQuant evidence -> MaxQtoMSstatsPTMFormat with either a paired
# global run or the explicit co-enriched-unmodified-peptide proxy -> dataSummarizationPTM ->
# groupComparisonPTM (explicit Treatment-vs-Control contrast) -> site rows of ADJUSTED.Model ->
# TREAT-style test against |log2FC| > lfc.
#
# Inputs (key=value; files are looked up inside dir=):
#   dir=<folder>            folder holding the input files (default .)
#   out=<folder>            where adjusted_sites.csv or proxy_adjusted_sites.csv is written (default .)
#   evidence=evidence_phospho.txt          enriched-run MaxQuant evidence.txt
#   annotation=annotation_ptm.csv          Raw.file, Condition, BioReplicate, IsotopeLabelType
#   fasta=uniprot_human.fasta
#   evidence_prot=evidence_global.txt      global (unenriched) run evidence.txt
#   proteinGroups=proteinGroups_global.txt
#   annotation_protein=annotation_protein.csv
#   use_unmod=FALSE         TRUE permits the no-global co-enriched-unmodified-peptide proxy
#   lfc=1                   log2 fold-change threshold inside the test
# Output: <out>/adjusted_sites.csv for a paired-global analysis, or
# <out>/proxy_adjusted_sites.csv for the no-global proxy. Both contain site rows of ADJUSTED.Model plus
# pvalue_lfc and adj.pvalue_lfc; proxy-adjusted rows are candidates, not regulation calls.
#   MSstatsPTM also writes its own log files into the working directory.
# Usage: Rscript scripts/msstatsptm_labelfree.R dir=<data_dir> out=<out_dir> [use_unmod=FALSE] [lfc=1]

output_exists <- function(path) file.exists(path) || nzchar(Sys.readlink(path))

validate_lfc <- function(value) {
  lfc <- suppressWarnings(as.numeric(value))
  if (length(lfc) != 1L || !is.finite(lfc) || lfc < 0) {
    stop('lfc must be a finite, non-negative number')
  }
  lfc
}

create_output_staging_dir <- function(output_dir) {
  stage_dir <- tempfile('.msstatsptm-staging-', tmpdir = output_dir)
  if (!dir.create(stage_dir, recursive = TRUE)) stop('Could not create private output staging directory: ', stage_dir)
  stage_dir
}

stage_csv_no_clobber <- function(data, destination, stage_dir = dirname(destination)) {
  if (output_exists(destination)) stop('Refusing to overwrite existing output: ', destination)
  stage <- file.path(stage_dir, basename(destination))
  if (output_exists(stage)) stop('Refusing to reuse staged output: ', stage)
  tryCatch(write.csv(data, stage, row.names = FALSE), error = function(e) {
    unlink(stage)
    stop(e)
  })
  list(stage = stage, destination = destination)
}

publish_staged_no_clobber <- function(staged, link_fn = file.link) {
  stages <- vapply(staged, `[[`, character(1), 'stage')
  on.exit(unlink(stages), add = TRUE)
  destinations <- vapply(staged, `[[`, character(1), 'destination')
  if (anyDuplicated(destinations)) stop('Each staged output must have a distinct destination.')
  if (any(vapply(destinations, output_exists, logical(1)))) {
    stop('Refusing to overwrite existing output: ', paste(destinations[vapply(destinations, output_exists, logical(1))], collapse = ', '))
  }
  published <- character()
  for (item in staged) {
    # A same-filesystem hard link is an atomic no-clobber publication: it fails if another process
    # creates the destination after the preflight, including a dangling symlink. Stages live under
    # opt$out so linking stays on the same filesystem.
    if (!link_fn(item$stage, item$destination)) {
      # Roll back only destinations successfully created by this invocation; pre-existing and
      # dangling-symlink destinations were rejected before the first link.
      unlink(published, force = TRUE)
      stop('Atomic no-clobber publication failed for: ', item$destination, '; rolled back this invocation.')
    }
    unlink(item$stage)
    published <- c(published, item$destination)
  }
}

make_proxy_manifest <- function(staged) {
  expected <- c('proxy_adjusted_sites.csv', 'proxy_ptm_model.csv', 'proxy_adjusted_qc.csv')
  artifacts <- basename(vapply(staged, `[[`, character(1), 'destination'))
  if (!identical(artifacts, expected)) {
    stop('Proxy manifest requires exactly these staged artifacts in order: ', paste(expected, collapse = ', '))
  }
  stages <- vapply(staged, `[[`, character(1), 'stage')
  data.frame(
    schema_version = rep(1L, length(expected)),
    mode = rep('no-global-proxy', length(expected)),
    complete = rep(TRUE, length(expected)),
    artifact = artifacts,
    bytes = file.info(stages)$size,
    md5 = unname(tools::md5sum(stages)),
    stringsAsFactors = FALSE
  )
}

validate_input_mode <- function(opt, inp) {
  if (!(opt$use_unmod %in% c('TRUE', 'FALSE'))) stop('use_unmod must be exactly TRUE or FALSE')
  use_unmod <- identical(opt$use_unmod, 'TRUE')
  global_inputs <- c(evidence_prot = inp('evidence_prot'),
                     proteinGroups = inp('proteinGroups'),
                     annotation_protein = inp('annotation_protein'))
  global_present <- file.exists(global_inputs)
  if (any(global_present) && !all(global_present)) {
    stop('Incomplete paired-global input: provide all of evidence_prot, proteinGroups, and annotation_protein, or provide none of them for use_unmod=TRUE. Missing: ',
         paste(names(global_inputs)[!global_present], collapse = ', '))
  }
  has_paired_global <- all(global_present)
  if (has_paired_global && use_unmod) {
    stop('use_unmod=TRUE is the no-global proxy route and is mutually exclusive with paired-global files; use use_unmod=FALSE for paired-global protein-adjusted calls.')
  }
  if (!has_paired_global && !use_unmod) {
    stop('No paired-global input found. Provide evidence_prot, proteinGroups, and annotation_protein for protein-adjusted calls, or set use_unmod=TRUE for the explicitly proxy-adjusted no-global route.')
  }
  list(use_unmod = use_unmod, has_paired_global = has_paired_global)
}

proxy_qc <- function(ev, annotation, is_mod) {
  required_ev <- c('Sequence', 'Proteins', 'Raw.file')
  required_annotation <- c('Raw.file', 'Condition', 'BioReplicate')
  if (!all(required_ev %in% names(ev))) stop('No-global proxy QC requires evidence columns: ', paste(required_ev, collapse = ', '))
  if (!all(required_annotation %in% names(annotation))) stop('No-global proxy QC requires annotation columns: ', paste(required_annotation, collapse = ', '))
  unmod <- ev[!is_mod, required_ev, drop = FALSE]
  if (!nrow(unmod)) stop('use_unmod=TRUE requires unmodified rows in the enriched evidence after the class-I filter; none were retained.')
  proteins <- trimws(as.character(unmod$Proteins))
  if (any(is.na(unmod$Sequence) | is.na(unmod$Proteins) | !nzchar(proteins) | proteins == 'NA' | grepl(';', proteins, fixed = TRUE))) {
    stop('No-global proxy QC requires unambiguous single-protein assignments for every retained unmodified peptide.')
  }
  peptide_proteins <- split(proteins, as.character(unmod$Sequence))
  ambiguous <- names(Filter(function(x) length(unique(x)) != 1L, peptide_proteins))
  if (length(ambiguous)) {
    stop('No-global proxy QC requires each retained unmodified peptide sequence to map to one protein; ambiguous sequences include: ',
         paste(head(ambiguous, 5L), collapse = ', '))
  }
  target_proteins <- trimws(as.character(ev$Proteins[is_mod]))
  if (any(is.na(ev$Proteins[is_mod]) | !nzchar(target_proteins) | target_proteins == 'NA' | grepl(';', target_proteins, fixed = TRUE))) {
    stop('No-global proxy QC requires unambiguous single-protein assignments for every retained modified peptide.')
  }
  target_proteins <- unique(target_proteins)
  if (any(is.na(ev$Raw.file) | !nzchar(as.character(ev$Raw.file)))) {
    stop('No-global proxy QC requires non-missing Raw.file values in retained evidence.')
  }
  annotation_rows <- annotation[, required_annotation, drop = FALSE]
  # Exact duplicate annotation rows are harmless and are deliberately collapsed; conflicting rows
  # are rejected below so one Raw.file cannot receive two experimental identities.
  ann <- unique(annotation_rows)
  evidence_raw_files <- unique(as.character(ev$Raw.file))
  absent_annotation <- setdiff(evidence_raw_files, as.character(ann$Raw.file))
  if (length(absent_annotation)) {
    stop('No-global proxy QC found retained evidence Raw.file values absent from annotation: ',
         paste(absent_annotation, collapse = ', '))
  }
  if (any(is.na(ann$Raw.file) | !nzchar(as.character(ann$Raw.file)) | is.na(ann$Condition) |
          !nzchar(as.character(ann$Condition)) | is.na(ann$BioReplicate) | !nzchar(as.character(ann$BioReplicate)))) {
    stop('No-global proxy QC requires non-missing Raw.file, Condition, and BioReplicate annotation values.')
  }
  raw_assignments <- split(paste(ann$Condition, ann$BioReplicate, sep = '\r'), as.character(ann$Raw.file))
  if (any(vapply(raw_assignments, function(x) length(unique(x)) != 1L, logical(1)))) {
    stop('No-global proxy QC requires a one-to-one Raw.file to Condition/BioReplicate annotation mapping.')
  }
  coverage <- do.call(rbind, unlist(lapply(target_proteins, function(protein) lapply(seq_len(nrow(ann)), function(i) {
    rows <- unmod[proteins == protein & as.character(unmod$Raw.file) == as.character(ann$Raw.file[i]), , drop = FALSE]
    data.frame(Protein = protein, Condition = as.character(ann$Condition[i]), BioReplicate = as.character(ann$BioReplicate[i]),
               Raw.file = as.character(ann$Raw.file[i]), unmodified_rows = nrow(rows),
               unique_unmodified_peptides = length(unique(rows$Sequence)), stringsAsFactors = FALSE)
  })), recursive = FALSE))
  missing <- coverage[coverage$unique_unmodified_peptides == 0L, , drop = FALSE]
  if (nrow(missing)) {
    stop('No-global proxy QC requires unmodified-peptide coverage for every PTM protein/Condition/BioReplicate; missing: ',
         paste(paste(missing$Protein, missing$Condition, missing$BioReplicate, missing$Raw.file, sep = '/'), collapse = ', '))
  }
  peptide_counts <- vapply(target_proteins, function(protein) {
    length(unique(unmod$Sequence[proteins == protein]))
  }, integer(1))
  too_few_peptides <- names(peptide_counts)[peptide_counts < 2L]
  if (length(too_few_peptides)) {
    stop('No-global proxy QC requires at least 2 unique unmodified peptides per PTM protein; insufficient: ',
         paste(too_few_peptides, collapse = ', '))
  }
  reps <- unique(coverage[, c('Protein', 'Condition', 'BioReplicate'), drop = FALSE])
  replicate_counts <- vapply(split(reps$BioReplicate, paste(reps$Protein, reps$Condition, sep = '\r')),
                            function(x) length(unique(x)), integer(1))
  too_few_reps <- names(replicate_counts)[replicate_counts < 2L]
  if (length(too_few_reps)) {
    stop('No-global proxy QC requires at least 2 biological replicates per condition for every PTM protein; insufficient: ',
         paste(too_few_reps, collapse = ', '))
  }
  coverage
}

if (!identical(Sys.getenv('MSSTATSPTM_LF_NO_AUTORUN'), '1')) {
opt <- list(dir = '.', out = '.', evidence = 'evidence_phospho.txt', annotation = 'annotation_ptm.csv',
            fasta = 'uniprot_human.fasta', evidence_prot = 'evidence_global.txt',
            proteinGroups = 'proteinGroups_global.txt', annotation_protein = 'annotation_protein.csv',
            use_unmod = 'FALSE', lfc = '1')
for (a in strsplit(commandArgs(trailingOnly = TRUE), '=', fixed = TRUE)) opt[[a[1]]] <- paste(a[-1], collapse = '=')
inp <- function(n) file.path(opt$dir, opt[[n]])
lfc <- validate_lfc(opt$lfc)

library(MSstatsPTM)
rd <- function(f) read.table(f, sep = '\t', header = TRUE, quote = '')

mode <- validate_input_mode(opt, inp)
use_unmod <- mode$use_unmod
has_paired_global <- mode$has_paired_global
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

# The converter uses the best-localized sequence as is and keeps unmodified peptides from the
# enriched runs, so apply the class-I rule to the enriched evidence FIRST: keep rows carrying the
# modification whose best site probability (from 'Phospho (STY) Probabilities') is >= 0.75.
# The use_unmod proxy needs the unmodified rows, so keep those too when it is on.
ev <- rd(inp('evidence'))
annotation <- read.csv(inp('annotation'))
site_prob <- vapply(regmatches(ev$Phospho..STY..Probabilities,
                               gregexpr('(?<=\\()[0-9.]+(?=\\))', ev$Phospho..STY..Probabilities, perl = TRUE)),
                    function(p) if (length(p)) max(as.numeric(p)) else NA_real_, numeric(1))
is_mod <- grepl('Phospho \\(STY\\)', ev$Modified.sequence)
keep <- is_mod & !is.na(site_prob) & site_prob >= 0.75
if (use_unmod) keep <- keep | !is_mod
ev <- ev[keep, ]
is_mod <- is_mod[keep]
proxy_coverage <- if (!has_paired_global) proxy_qc(ev, annotation, is_mod) else NULL

# Converters are <Tool>toMSstatsPTMFormat and return a list with $PTM and $PROTEIN.
# MaxQtoMSstatsPTMFormat reads the MaxQuant 'evidence.txt' (NOT the Phospho (STY)Sites
# table -- the pandas multiplicity-expansion is a SEPARATE workflow); the FASTA maps
# peptides back to site coordinates. With a paired global run, $PROTEIN comes from that run. With
# use_unmod=TRUE and no global files, MSstatsPTM builds $PROTEIN from co-enriched unmodified peptides;
# this is a weaker proxy and is never labelled as a regulation result below.
converter_args <- list(
  evidence = ev,
  annotation = annotation,
  fasta_path = inp('fasta'),
  mod_id = '\\(Phospho \\(STY\\)\\)',
  which_proteinid_ptm = 'Proteins',
  use_unmod_peptides = use_unmod
)
if (has_paired_global) {
  converter_args <- c(converter_args, list(
    evidence_prot = rd(inp('evidence_prot')),
    proteinGroups = rd(inp('proteinGroups')),
    annotation_protein = read.csv(inp('annotation_protein')),
    which_proteinid_protein = 'Proteins'
  ))
}
input <- do.call(MaxQtoMSstatsPTMFormat, converter_args)
if (!('PROTEIN' %in% names(input))) {
  stop(if (has_paired_global) {
    'MSstatsPTM did not construct PROTEIN from the supplied paired-global input; cannot make protein-adjusted calls.'
  } else {
    'MSstatsPTM did not construct PROTEIN from co-enriched unmodified peptides; cannot make proxy-adjusted calls.'
  })
}

# append defaults to TRUE and requires a log file, so set it FALSE when use_log_file = FALSE
summarized <- dataSummarizationPTM(input, use_log_file = FALSE, append = FALSE)

# data.type is 'LabelFree' (DDA/DIA label-free) or 'TMT' -- NOT the converter's labeling_type 'LF'.
# Pass an explicit contrast: the default pairwise Label is 'Control vs Treatment' (log2FC = Control - Treatment).
# Columns must follow the sorted Condition levels.
contrast <- matrix(c(-1, 1), nrow = 1, dimnames = list('Treatment vs Control', c('Control', 'Treatment')))
result <- groupComparisonPTM(summarized, data.type = 'LabelFree', contrast.matrix = contrast)

# Three models; the adjusted one is the paired-global deliverable or the explicitly labelled no-global
# proxy. Keep site rows only (Protein_<residue><position>).
adjusted <- result$ADJUSTED.Model
adjusted <- adjusted[grepl('_[STY][0-9]+', adjusted$Protein), ]

# A claim of "changed more than 2-fold" needs the threshold INSIDE the test (TREAT-style), not
# adj.pvalue < 0.05 & |log2FC| > 1 as a post-hoc double filter, whose FDR refers to log2FC != 0.
adjusted$pvalue_lfc <- pt((abs(adjusted$log2FC) - lfc) / adjusted$SE, adjusted$DF, lower.tail = FALSE)
adjusted$adj.pvalue_lfc <- p.adjust(adjusted$pvalue_lfc, method = 'BH')
if (has_paired_global) {
  regulated <- adjusted[!is.na(adjusted$adj.pvalue_lfc) & adjusted$adj.pvalue_lfc < 0.05, ]
  output_file <- 'adjusted_sites.csv'
  result_kind <- 'paired-global protein-adjusted regulated sites (TREAT)'
  result_count <- nrow(regulated)
} else {
  ptm_model <- result$PTM.Model
  ptm_model <- ptm_model[grepl('_[STY][0-9]+', ptm_model$Protein), c('Protein', 'Label', 'log2FC', 'adj.pvalue')]
  names(ptm_model)[3:4] <- c('raw_ptm_log2FC', 'raw_ptm_adj.pvalue')
  adjusted <- merge(adjusted, ptm_model, by = c('Protein', 'Label'), all.x = TRUE, sort = FALSE)
  adjusted$adjustment_source <- 'co-enriched unmodified peptides; no paired global proteome'
  adjusted$interpretation <- 'proxy-adjusted candidate; not a regulation call'
  output_file <- 'proxy_adjusted_sites.csv'
  result_kind <- 'proxy-adjusted candidates (not regulation calls)'
  result_count <- sum(!is.na(adjusted$adj.pvalue_lfc) & adjusted$adj.pvalue_lfc < 0.05)
}

# How much of each call was protein-driven: compare PTM.Model vs ADJUSTED.Model. All outputs are first
# staged, then published without replacing an existing file (including a dangling symlink).
stage_dir <- create_output_staging_dir(opt$out)
staged <- list(stage_csv_no_clobber(adjusted, file.path(opt$out, output_file), stage_dir))
if (!has_paired_global) {
  staged <- c(staged, list(
    stage_csv_no_clobber(result$PTM.Model, file.path(opt$out, 'proxy_ptm_model.csv'), stage_dir),
    stage_csv_no_clobber(proxy_coverage, file.path(opt$out, 'proxy_adjusted_qc.csv'), stage_dir)
  ))
  # Consumers accept the proxy set only if this manifest is present and validates all three artifacts.
  # Stage it before the worker-local artifacts and append it last within this private output stage.
  # This is only the internal integrity manifest; run_proxy_checked.py writes the public outer
  # proxy_checked_complete.json marker after the R process exits cleanly.
  completion <- make_proxy_manifest(staged)
  staged <- c(staged, list(stage_csv_no_clobber(
    completion, file.path(opt$out, 'proxy_adjusted_manifest.csv'), stage_dir
  )))
}
tryCatch({
  publish_staged_no_clobber(staged)
}, finally = unlink(stage_dir, recursive = TRUE, force = TRUE))
cat('mode:', if (has_paired_global) 'paired-global' else 'no-global-proxy',
    '| names(input):', names(input), '| ADJUSTED site rows:', nrow(adjusted), '|', result_kind, ':', result_count,
    '| Label:', unique(as.character(adjusted$Label)), '\n')
}
