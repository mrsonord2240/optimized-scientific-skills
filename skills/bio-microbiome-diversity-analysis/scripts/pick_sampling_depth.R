# Reference: vegan 2.6+ | Verify API if version differs
# Numeric stand-in for reading the alpha-rarefaction plateau by eye.
# Picks the smallest rarefaction depth at which 1000 more reads add less than `tol` (relative, default
# 1%) to mean expected richness, restricted to depths that keep >= `min_keep` of the samples.
# The gain at each step is measured on the samples that survive BOTH depths (paired), so samples
# falling out of the pool do not bias the curve. Reports the depth, survivors and dropped samples.
#
# Input : feature x sample integer count table (TSV, first column = feature IDs; QIIME2:
#         `qiime tools export --input-path table.qza --output-path e; biom convert -i e/feature-table.biom -o counts.tsv --to-tsv`,
#         then `sed -i 1d counts.tsv` to drop the "# Constructed from biom file" line; the "#OTU ID" header is read as-is).
# Usage : Rscript pick_sampling_depth.R counts.tsv [tol=0.01] [min_keep=0.85] [n_steps=40]
# Output: prints the curve, the chosen depth, dropped samples; writes depth_curve.tsv to the working directory.
suppressPackageStartupMessages(library(vegan))

args     <- commandArgs(trailingOnly = TRUE)
counts_f <- args[1]
tol      <- if (length(args) >= 2) as.numeric(args[2]) else 0.01
min_keep <- if (length(args) >= 3) as.numeric(args[3]) else 0.85
n_steps  <- if (length(args) >= 4) as.integer(args[4]) else 40
if (is.na(counts_f)) stop('usage: Rscript pick_sampling_depth.R counts.tsv [tol] [min_keep] [n_steps]')

tab <- read.delim(counts_f, row.names = 1, check.names = FALSE, comment.char = '', skip = 0)
tab <- tab[, !grepl('^taxonomy$', colnames(tab), ignore.case = TRUE), drop = FALSE]
x <- t(as.matrix(tab))                     # samples x features
storage.mode(x) <- 'integer'
tot <- rowSums(x)
cat('Samples:', nrow(x), '| depth range:', min(tot), '-', max(tot), '| median:', median(tot), '\n')

# depth grid: from ~1% of the median up to the deepest sample
grid <- unique(as.integer(round(seq(max(50, 0.01 * median(tot)), max(tot), length.out = n_steps))))
# expected (analytical, RNG-free) rarefied richness of each sample at each depth it can reach
rich <- sapply(grid, function(d) {
    ok <- tot >= d
    r  <- rep(NA_real_, nrow(x)); r[ok] <- rarefy(x[ok, , drop = FALSE], sample = d)
    r
})
rownames(rich) <- rownames(x)

curve <- data.frame(depth = grid, survivors = colSums(!is.na(rich)), mean_richness = colMeans(rich, na.rm = TRUE),
                    gain_per_1000 = NA_real_)
for (i in seq_len(length(grid) - 1)) {
    both <- !is.na(rich[, i + 1])          # survives the next depth, hence this one too
    if (any(both)) curve$gain_per_1000[i] <- mean(rich[both, i + 1] - rich[both, i]) / mean(rich[both, i]) /
        ((grid[i + 1] - grid[i]) / 1000)   # relative richness gain per 1000 added reads
}
curve$frac_kept <- curve$survivors / nrow(x)
write.table(curve, 'depth_curve.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
print(curve, digits = 3, row.names = FALSE)

allowed <- which(curve$frac_kept >= min_keep)
hit     <- allowed[!is.na(curve$gain_per_1000[allowed]) & curve$gain_per_1000[allowed] < tol]
if (length(hit)) {
    chosen <- curve$depth[hit[1]]; note <- sprintf('plateau reached (gain < %.1f%% per 1000 reads)', 100 * tol)
} else {
    chosen <- curve$depth[max(allowed)]; note <- sprintf('NO plateau within the >=%.0f%%-survivor range; deepest allowed depth used, richness is still rising - report that', 100 * min_keep)
}
dropped <- rownames(x)[tot < chosen]
cat(sprintf('\nChosen sampling depth: %d (%s)\n', chosen, note))
cat(sprintf('Survivors: %d of %d | dropped: %s\n', sum(tot >= chosen), nrow(x),
            if (length(dropped)) paste(dropped, collapse = ', ') else 'none'))
cat('Confirm the conclusion at a nearby depth (e.g. +/- 20%) before reporting.\n')
