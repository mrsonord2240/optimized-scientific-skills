# bridge_layout.R -- reserve one channel per plex for a pooled reference/bridge sample (TMT,
#                    multiplexed proteomics), then assign and verify the remaining channels.
# Purpose: make cross-plex normalization possible without reusing a biological channel.
# Usage:   Rscript scripts/bridge_layout.R samples.csv n_plex n_channel reserved_channel feature_vars out.csv [max_iter] [seed]
#   samples.csv     : one row per biological sample, an `id` column plus the feature_vars columns
#   n_plex          : number of plexes;  n_channel: channels per plex
#   reserved_channel: channel index held by the pooled bridge sample in every plex (e.g. 16)
#   feature_vars    : comma-separated columns to balance across plexes, e.g. condition,site
#   out.csv         : layout of the biological samples (reserved channel absent from it)
#   max_iter        : default 10000;  seed: default 1
# Example: Rscript scripts/bridge_layout.R samples.csv 4 16 16 condition,site layout.csv
# Randomize LC-MS injection order separately, independently of plex/channel.
args <- commandArgs(trailingOnly = TRUE)
stopifnot('usage: bridge_layout.R samples.csv n_plex n_channel reserved_channel feature_vars out.csv [max_iter] [seed]' =
            length(args) >= 6)
samples_file <- args[1]; n_plex <- as.integer(args[2]); n_channel <- as.integer(args[3])
reserved_channel <- as.integer(args[4]); feature_vars <- strsplit(args[5], ',')[[1]]; out_file <- args[6]
max_iter <- if (length(args) >= 7) as.integer(args[7]) else 10000L
set.seed(if (length(args) >= 8) as.integer(args[8]) else 1L)

script_dir <- dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, 'check_balance.R'))
suppressPackageStartupMessages(library(designit))

samples <- read.csv(samples_file, stringsAsFactors = FALSE)
stopifnot('samples.csv needs an id column and every feature_vars column' =
            all(c('id', feature_vars) %in% names(samples)),
          'more samples than non-reserved positions' = nrow(samples) <= n_plex * (n_channel - 1))

# `exclude` needs one row per excluded (plex, channel) combination, not just the channel column alone.
bc <- BatchContainer$new(
  dimensions = list(plex = n_plex, channel = n_channel),
  exclude = data.frame(plex = seq_len(n_plex), channel = reserved_channel))
bc <- assign_in_order(bc, samples = samples)     # only fills the non-reserved positions
bc <- optimize_design(
  bc,
  scoring = osat_score_generator(batch_vars = 'plex', feature_vars = feature_vars),
  max_iter = max_iter)
assignment <- bc$get_samples()

# Verify: the reserved channel must stay empty (it holds the pooled reference, added outside this
# table) and no plex may exceed its non-reserved positions.
stopifnot(
  'the reserved channel was assigned a biological sample' = !reserved_channel %in% assignment$channel,
  'a plex holds more biological samples than it has non-reserved channels' =
    all(table(assignment$plex) <= n_channel - 1)
)
check_balance(assignment, 'plex', feature_vars)
write.csv(assignment, out_file, row.names = FALSE)
cat('layout written to', out_file, '\n')
