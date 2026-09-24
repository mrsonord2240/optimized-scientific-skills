#!/usr/bin/env Rscript
# Static regression for the direct-TMT invocation boundary. This does not invoke MSstatsPTM;
# the positive workflow test owns package/runtime execution.

this_file <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))][1])
script <- normalizePath(file.path(dirname(this_file), '..', 'scripts', 'msstatsptm_tmt.R'))
lines <- readLines(script, warn = FALSE)
line_number <- function(pattern) {
  found <- grep(pattern, lines, fixed = TRUE)
  stopifnot(length(found) >= 1L)
  found[1]
}

# An invalid threshold and an existing (including dangling-link) destination must fail before the
# expensive converter/modeling branch. Publication must remain a staged hard-link create, never a
# write or replacement at the final destination.
stopifnot(
  line_number("lfc <- suppressWarnings(as.numeric(opt$lfc))") < line_number('library(MSstatsPTM)'),
  line_number("if (output_exists(output_file)) stop('Refusing to overwrite existing output: ', output_file)") <
    line_number('library(MSstatsPTM)'),
  line_number("output_exists <- function(path) file.exists(path) || nzchar(Sys.readlink(path))") <
    line_number('library(MSstatsPTM)'),
  line_number('stage_dir <- tempfile') < line_number('write.csv(data, stage_file, row.names = FALSE)'),
  line_number('write.csv(data, stage_file, row.names = FALSE)') < line_number('file.link(stage_file, destination)'),
  line_number('adjusted <- as.data.frame(result$ADJUSTED.Model, stringsAsFactors = FALSE)') <
    line_number('adjusted <- adjusted[, output_columns, drop = FALSE]'),
  line_number('publish_csv_no_clobber(adjusted, output_file)') > line_number('output_columns <- c('),
  any(grepl('diagnostic-only', lines, fixed = TRUE)),
  any(grepl('run_checked.py', lines, fixed = TRUE)),
  !any(grepl('write.csv(adjusted, output_file', lines, fixed = TRUE))
)
message('PASS TMT direct-output contract: ', script)
