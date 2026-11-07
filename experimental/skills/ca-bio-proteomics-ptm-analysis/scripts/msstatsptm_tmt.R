#!/usr/bin/env Rscript
# Protein-adjusted phosphosite testing with MSstatsPTM, TMT / isobaric route.
#
# Same adjustment as msstatsptm_labelfree.R; only three calls differ: labeling_type = 'TMT' on the
# converter, dataSummarizationPTM_TMT, and data.type = 'TMT'. Every plex needs a pooled reference
# channel carried as Condition = 'Norm' in the annotation.
#
# Inputs (key=value; files are looked up inside dir=):
#   dir=<folder>   out=<folder>   lfc=1 (non-negative log2 fold-change threshold inside the test)
#   evidence=evidence_phospho_tmt.txt   annotation=annotation_ptm_tmt.csv   fasta=uniprot_human.fasta
#   evidence_prot=evidence_global_tmt.txt   proteinGroups=proteinGroups_global.txt
#   annotation_protein=annotation_protein_tmt.csv
#   Channel names in the annotations follow MaxQuant's 0-indexed reporter suffixes (channel.0 .. channel.9).
# Output: <out>/adjusted_sites_tmt.csv: site rows of ADJUSTED.Model with the declared columns
#   Protein, Label, log2FC, SE, Tvalue, DF, pvalue, adj.pvalue, GlobalProtein, Adjusted,
#   pvalue_lfc, adj.pvalue_lfc. MSstatsPTM/MSstatsTMT write log files into the working directory.
# Usage: Rscript scripts/msstatsptm_tmt.R dir=<data_dir> out=<out_dir> [lfc=1]
# A bare Rscript invocation is diagnostic-only. Unattended completion requires scripts/run_checked.py
# with a private stage directory and this output as its single validated, no-clobber publication.
opt <- list(dir = '.', out = '.', evidence = 'evidence_phospho_tmt.txt', annotation = 'annotation_ptm_tmt.csv',
            fasta = 'uniprot_human.fasta', evidence_prot = 'evidence_global_tmt.txt',
            proteinGroups = 'proteinGroups_global.txt', annotation_protein = 'annotation_protein_tmt.csv', lfc = '1')
for (a in strsplit(commandArgs(trailingOnly = TRUE), '=', fixed = TRUE)) opt[[a[1]]] <- paste(a[-1], collapse = '=')
lfc <- suppressWarnings(as.numeric(opt$lfc))
if (length(lfc) != 1L || !is.finite(lfc) || lfc < 0) stop('lfc must be a finite, non-negative number')
inp <- function(n) file.path(opt$dir, opt[[n]])
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(opt$out, 'adjusted_sites_tmt.csv')
output_exists <- function(path) file.exists(path) || nzchar(Sys.readlink(path))
if (output_exists(output_file)) stop('Refusing to overwrite existing output: ', output_file)

publish_csv_no_clobber <- function(data, destination) {
  # Create the stage under opt$out so file.link is same-filesystem and is an atomic create-only
  # operation. Do not fall back to copy/rename: those can replace an existing destination.
  stage_dir <- tempfile('.msstatsptm-tmt-stage-', tmpdir = dirname(destination))
  if (!dir.create(stage_dir, recursive = TRUE)) stop('Could not create private output staging directory: ', stage_dir)
  on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  stage_file <- file.path(stage_dir, basename(destination))
  tryCatch(write.csv(data, stage_file, row.names = FALSE), error = function(e) stop(e))
  if (output_exists(destination)) stop('Refusing to overwrite existing output: ', destination)
  if (!file.link(stage_file, destination)) {
    stop('Atomic no-clobber publication failed for: ', destination)
  }
}

library(MSstatsPTM)
rd <- function(f) read.table(f, sep = '\t', header = TRUE, quote = '')

# Class-I pre-filter on the ENRICHED evidence, exactly as in the label-free route.
ev <- rd(inp('evidence'))
site_prob <- vapply(regmatches(ev$Phospho..STY..Probabilities,
                               gregexpr('(?<=\\()[0-9.]+(?=\\))', ev$Phospho..STY..Probabilities, perl = TRUE)),
                    function(p) if (length(p)) max(as.numeric(p)) else NA_real_, numeric(1))
ev <- ev[grepl('Phospho \\(STY\\)', ev$Modified.sequence) & !is.na(site_prob) & site_prob >= 0.75, ]

# TMT annotation: Run, Fraction, TechRepMixture, Channel, Condition, Mixture, BioReplicate.
# Channel names follow the reporter-column suffixes MaxQuant wrote, and those are 0-indexed: a
# 10-plex has 'Reporter intensity corrected 0' .. '9', so the annotation needs 'channel.0' .. 'channel.9'
# (the CORRECTED columns, not the raw reporters). Read the suffixes off your own evidence header;
# 'channel.1' .. 'channel.10' is rejected with 'the channel name must be matched with that in input
# data', which never mentions the off-by-one.
# Give the pooled reference channel Condition = 'Norm' in EVERY plex.
input <- MaxQtoMSstatsPTMFormat(
  evidence = ev,
  annotation = read.csv(inp('annotation')),
  fasta_path = inp('fasta'),
  evidence_prot = rd(inp('evidence_prot')),
  proteinGroups = rd(inp('proteinGroups')),
  annotation_protein = read.csv(inp('annotation_protein')),
  labeling_type = 'TMT',          # the single converter switch; default is 'LF'
  mod_id = '\\(Phospho \\(STY\\)\\)',
  which_proteinid_ptm = 'Proteins',
  which_proteinid_protein = 'Proteins',
  use_unmod_peptides = FALSE)
stopifnot('PROTEIN' %in% names(input))

# TMT summarization is its own function. reference_norm / reference_norm.PTM (default TRUE) apply
# the 'Norm'-channel bridge; remove_norm_channel (default TRUE) drops that channel afterwards, so
# the contrast below names only the biological conditions.
summarized <- dataSummarizationPTM_TMT(input, use_log_file = FALSE, append = FALSE)

contrast <- matrix(c(-1, 1), nrow = 1, dimnames = list('Treatment vs Control', c('Control', 'Treatment')))
result <- groupComparisonPTM(summarized, data.type = 'TMT', contrast.matrix = contrast)

adjusted <- as.data.frame(result$ADJUSTED.Model, stringsAsFactors = FALSE)
adjusted <- adjusted[grepl('_[STY][0-9]+', adjusted$Protein), ]

# Keep the effect-size threshold inside the hypothesis test, as in the label-free route. A post-hoc
# |log2FC| filter would leave the adjusted p-value referring to the null log2FC == 0.
adjusted$pvalue_lfc <- pt((abs(adjusted$log2FC) - lfc) / adjusted$SE, adjusted$DF, lower.tail = FALSE)
adjusted$adj.pvalue_lfc <- p.adjust(adjusted$pvalue_lfc, method = 'BH')
regulated <- adjusted[!is.na(adjusted$adj.pvalue_lfc) & adjusted$adj.pvalue_lfc < 0.05, ]

# Publish only the site rows and this stable, declared schema; fail rather than silently producing a
# partial table if a supported MSstatsPTM version changes a required model column.
output_columns <- c('Protein', 'Label', 'log2FC', 'SE', 'Tvalue', 'DF', 'pvalue', 'adj.pvalue',
                    'GlobalProtein', 'Adjusted', 'pvalue_lfc', 'adj.pvalue_lfc')
missing_columns <- setdiff(output_columns, names(adjusted))
if (length(missing_columns)) stop('ADJUSTED.Model is missing required output columns: ', paste(missing_columns, collapse = ', '))
adjusted <- adjusted[, output_columns, drop = FALSE]
publish_csv_no_clobber(adjusted, output_file)
cat('names(input):', names(input), '| ADJUSTED site rows:', nrow(adjusted), '| regulated (TREAT):', nrow(regulated),
    '| Label:', unique(as.character(adjusted$Label)), '\n')
