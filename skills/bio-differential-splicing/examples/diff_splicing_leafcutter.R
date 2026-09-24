#!/usr/bin/env Rscript
# Reference: leafcutter 0.2.9, regtools 1.0.0 (checked 2026-09-20)
# Differential splicing analysis with leafcutter
# Intron cluster-based approach for novel junction discovery
#
# Run in a directory holding one <sample>.bam (+ .bai) per sample. regtools, python (3) and Rscript
# must be on PATH; LEAFCUTTER_DIR is a clone of https://github.com/davidaknowles/leafcutter
# (it holds clustering/ and scripts/; the R package does not put leafcutter_ds.R on PATH).

LEAFCUTTER_DIR <- Sys.getenv('LEAFCUTTER_DIR', 'leafcutter')
# Sample names = BAM basenames without .bam (they become the .junc names and the counts-table column names)
CONTROL <- c('control1', 'control2', 'control3')
TREATMENT <- c('treatment1', 'treatment2', 'treatment3')
NONSTANDARD_CONTIGS <- FALSE   # TRUE if contigs are not chr1..22,X,Y / 1..22,X,Y (clustering then drops every junction unless -k True)
EXON_FILE <- NULL              # optional, only labels clusters with gene names; make with leafcutter/scripts/gtf_to_exons.R gencode.gtf.gz gencode_exons.txt.gz

run <- function(cmd, args) {
    status <- system2(cmd, args)
    if (status != 0) stop(sprintf('%s failed (exit %d)', cmd, status))
}

# Step 1: BAMs -> junction files (strand from the XS tag: the BAM needs XS, e.g. STAR --outSAMstrandField intronMotif)
samples <- c(CONTROL, TREATMENT)
for (s in samples) {
    run('regtools', c('junctions', 'extract', '-a', '8', '-m', '50', '-s', 'XS', paste0(s, '.bam'), '-o', paste0(s, '.junc')))
}

# Step 2: junction file list + intron clustering (-m: min reads per cluster, lower it for shallow data; -l: max intron length)
writeLines(paste0(samples, '.junc'), 'juncfiles.txt')
run('python', c(file.path(LEAFCUTTER_DIR, 'clustering', 'leafcutter_cluster_regtools.py'),
                '-j', 'juncfiles.txt', '-o', 'leafcutter', '-m', '50', '-l', '500000',
                if (NONSTANDARD_CONTIGS) c('-k', 'True')))

# Step 3: groups file. Sample names must equal the counts-table column names (the .junc basenames)
counts_file <- 'leafcutter_perind_numers.counts.gz'
if (!file.exists(counts_file)) {
    stop('Clustering did not write ', counts_file, '; check LEAFCUTTER_DIR and the clustering log.')
}
counts_rows <- length(readLines(gzfile(counts_file), warn = FALSE)) - 1L
if (counts_rows < 1L) {
    stop('Clustering produced zero introns. Use -k True for nonstandard contigs, lower -m for shallow data, and verify BAM XS tags.')
}
counts_samples <- scan(gzfile(counts_file), what = '', nlines = 1, quiet = TRUE)
groups <- data.frame(sample = samples, group = rep(c('control', 'treatment'), c(length(CONTROL), length(TREATMENT))))
missing <- setdiff(groups$sample, counts_samples)
if (length(missing) > 0) stop('groups.txt names not in the counts table: ', paste(missing, collapse = ', '))
write.table(groups, 'groups.txt', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

# Step 4: differential splicing. leafcutter_ds.R defaults (-i 5 -g 3 -c 20) stop when a group has < 5 samples:
# set -i and -g to the smaller group size (capped at 5); -c is the read threshold per sample
n_min <- min(length(CONTROL), length(TREATMENT), 5)
ds_args <- c(file.path(LEAFCUTTER_DIR, 'scripts', 'leafcutter_ds.R'), '-i', n_min, '-g', n_min, '-c', 10,
             '-o', 'differential', if (!is.null(EXON_FILE)) c('-e', EXON_FILE), counts_file, 'groups.txt')
run('Rscript', ds_args)

# Step 5: load and analyze results
load_leafcutter_results <- function(cluster_sig_file, effect_sizes_file) {
    sig <- read.table(cluster_sig_file, header = TRUE, sep = '\t')
    effects <- read.table(effect_sizes_file, header = TRUE, sep = '\t')

    # Filter significant clusters
    # FDR < 0.05 (p.adjust column)
    significant <- sig[!is.na(sig$p.adjust) & sig$p.adjust < 0.05, ]

    cat(sprintf('Total clusters tested: %d\n', sum(sig$status == 'Success')))
    cat(sprintf('Significant clusters (FDR < 0.05): %d\n', nrow(significant)))

    # Merge with effect sizes. The effect-sizes file is per-intron (no `cluster`
    # column); its cluster id lives inside the `intron` string (chr:start:end:clu_N).
    # The significance file keys on `cluster` (chr:clu_N). Match on the run-unique clu_N.
    extract_clu <- function(x) sub('.*?(clu_[0-9]+).*', '\\1', x)
    significant$clu <- extract_clu(significant$cluster)
    effects$clu <- extract_clu(effects$intron)
    results <- merge(significant, effects, by = 'clu', all.x = TRUE)

    # Sort by significance
    results <- results[order(results$p.adjust), ]

    return(results)
}

# Parse the intron key (chr:start:end:clu_N) into coordinates
parse_intron_coordinates <- function(results) {
    parts <- strsplit(results$intron, ':')
    results$chr <- sapply(parts, '[', 1)
    results$start <- as.numeric(sapply(parts, '[', 2))
    results$end <- as.numeric(sapply(parts, '[', 3))
    return(results)
}

results <- load_leafcutter_results('differential_cluster_significance.txt', 'differential_effect_sizes.txt')
results <- parse_intron_coordinates(results)
# deltapsi = group2 (treatment) - group1 (control)
print(head(results[, c('intron', 'p.adjust', 'deltapsi', 'chr', 'start', 'end')]))
