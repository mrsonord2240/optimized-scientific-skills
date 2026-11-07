#!/usr/bin/env Rscript
# Checked on DADA2 1.30+
# Assign taxonomy to ASVs with DADA2 naive Bayes (assignTaxonomy) and, optionally, exact-match species (addSpecies).
# A short 16S read licenses genus at best; species comes only from an exact match.
#
# Usage: Rscript dada2_assign.R seqtab_nochim.rds train.fa.gz out.rds [species.fa.gz] [rerun]
#   train.fa.gz    DADA2-formatted (rank-labelled headers) reference, ideally trimmed to the amplicon region
#   species.fa.gz  optional species-level reference for the exact-match step
#   rerun          repeat the seeded classification once and require identical genus calls (doubles the run time)
#   MIN_BOOT       env, default 50 (DADA2 default and the RDP recommendation for reads <=250 nt; 80 is a stricter choice)
library(dada2)

args <- commandArgs(trailingOnly = TRUE)
rerun_check <- 'rerun' %in% args
args <- setdiff(args, 'rerun')
if (length(args) < 3) stop('Usage: Rscript dada2_assign.R seqtab_nochim.rds train.fa.gz out.rds [species.fa.gz] [rerun]')
seqtab_file <- args[1]
train_file <- args[2]
out_file <- args[3]
species_file <- if (length(args) >= 4) args[4] else NULL
minBoot <- as.numeric(Sys.getenv('MIN_BOOT', '50'))

for (f in c(seqtab_file, train_file)) if (!file.exists(f)) stop('not found: ', f)
seqtab_nochim <- readRDS(seqtab_file)
cat('ASVs to classify:', ncol(seqtab_nochim), '\n')

# assignTaxonomy() runs 100 stochastic bootstrap resamples per sequence; seed before EVERY call or two runs
# differ at the genus call for ~2-3% of ASVs. Seeded runs are reproducible at genus, not bitwise at every rank
# (about 2/770 cells at Kingdom/Order can still differ; multithread=FALSE does not remove that).
set.seed(100)
taxa <- assignTaxonomy(seqtab_nochim, train_file, minBoot = minBoot, tryRC = TRUE, multithread = TRUE)

# A wrong-region or wrong-marker reference, bad orientation or mis-formatted FASTA gives all-NA calls with no error.
genus_frac <- mean(!is.na(taxa[, 'Genus']))
if (genus_frac == 0) stop('assignTaxonomy returned no genus calls: check the reference format, region and marker')
cat(sprintf('Genus assigned: %.1f%% of ASVs\n', 100 * genus_frac))

if (rerun_check) {
    set.seed(100)
    again <- assignTaxonomy(seqtab_nochim, train_file, minBoot = minBoot, tryRC = TRUE, multithread = TRUE)
    n_diff <- sum(is.na(taxa[, 'Genus']) != is.na(again[, 'Genus']) |
                  (!is.na(taxa[, 'Genus']) & !is.na(again[, 'Genus']) & taxa[, 'Genus'] != again[, 'Genus']))
    cat(sprintf('Rerun check: %d genus calls differ between two seeded runs\n', n_diff))
    if (n_diff > 0) stop('genus calls are not reproducible between seeded runs')
}

# addSpecies: species ONLY by exact (100%) match against the species reference; most ASVs stay NA at species.
if (!is.null(species_file)) {
    if (!file.exists(species_file)) stop('species reference not found: ', species_file)
    taxa <- addSpecies(taxa, species_file)
}

cat('\nAssigned fraction per rank (NA = unassigned at threshold, kept honestly):\n')
for (rank in colnames(taxa)) {
    assigned <- sum(!is.na(taxa[, rank]))
    cat(sprintf('  %-8s %d/%d (%.1f%%)\n', rank, assigned, nrow(taxa), 100 * assigned / nrow(taxa)))
}

# State the conditioning choices alongside the labels (classifier, database file, confidence floor).
attr(taxa, 'classifier') <- 'DADA2 assignTaxonomy (RDP naive Bayes)'
attr(taxa, 'reference') <- basename(train_file)
attr(taxa, 'minBoot') <- minBoot

saveRDS(taxa, out_file)
cat('\nSaved', out_file, '(genus via naive Bayes, species via exact match only)\n')
