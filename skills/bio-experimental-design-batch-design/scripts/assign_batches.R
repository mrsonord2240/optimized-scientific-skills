# assign_batches.R -- constrained sample-to-batch assignment with designit, then verify the layout.
# Purpose: balance the biological factors and covariates across batches so batch is estimable.
# Usage:   Rscript scripts/assign_batches.R samples.csv n_batch n_position feature_vars out.csv [max_iter] [seed]
#   samples.csv  : one row per sample, an `id` column plus the columns named in feature_vars
#   n_batch      : number of batches/lanes/plates;  n_position: positions per batch
#   feature_vars : comma-separated columns to balance, e.g. condition,sex
#   out.csv      : layout written here (designit's get_samples() table, with batch and position)
#   max_iter     : default 10000 (raise it if check_balance() warns);  seed: default 1
# Example: Rscript scripts/assign_batches.R samples.csv 3 8 condition,sex layout.csv
args <- commandArgs(trailingOnly = TRUE)
stopifnot('usage: assign_batches.R samples.csv n_batch n_position feature_vars out.csv [max_iter] [seed]' =
            length(args) >= 5)
samples_file <- args[1]; n_batch <- as.integer(args[2]); n_position <- as.integer(args[3])
feature_vars <- strsplit(args[4], ',')[[1]]; out_file <- args[5]
max_iter <- if (length(args) >= 6) as.integer(args[6]) else 10000L
set.seed(if (length(args) >= 7) as.integer(args[7]) else 1L)

script_dir <- dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, 'check_balance.R'))
suppressPackageStartupMessages(library(designit))          # verify API against installed vignette

samples <- read.csv(samples_file, stringsAsFactors = FALSE)
stopifnot('samples.csv needs an id column and every feature_vars column' =
            all(c('id', feature_vars) %in% names(samples)),
          'more samples than n_batch x n_position slots' = nrow(samples) <= n_batch * n_position)
bc <- BatchContainer$new(dimensions = list(batch = n_batch, position = n_position))
bc <- assign_in_order(bc, samples = samples)
bc <- optimize_design(
  bc,
  scoring = osat_score_generator(batch_vars = 'batch', feature_vars = feature_vars),  # balance every factor
  max_iter = max_iter)
assignment <- bc$get_samples()                # R6 method on the container (no standalone get_samples())

check_balance(assignment, 'batch', feature_vars)   # every balanced covariate, not only condition
write.csv(assignment, out_file, row.names = FALSE)
cat('layout written to', out_file, '\n')
