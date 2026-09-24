# Reference: DADA2 1.30+, DECIPHER 2.30+ | Verify API if version differs
# Assign taxonomy to 16S ASVs with DADA2 naive Bayes (genus) + addSpecies (exact-match species).
# Honest message: a short 16S read licenses genus at best; species comes only from an exact match.
# Usage: Rscript assign_silva.R [seqtab.rds] [train.fa.gz] [species.fa.gz] [out.rds] [rerun]
#   rerun = the literal word 'rerun' repeats the classification once and diffs the genus calls
#   (doubles the run time).
library(dada2)

args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
seqtab_file <- arg(1, 'seqtab_nochim.rds')
silva_train <- arg(2, 'silva_nr99_v138.1_train_set.fa.gz')        # DADA2-formatted, rank-labelled headers
silva_species <- arg(3, 'silva_species_assignment_v138.1.fa.gz')  # species-level reference for exact match
out_file <- arg(4, 'taxa.rds')
rerun_check <- identical(arg(5, ''), 'rerun')

seqtab_nochim <- readRDS(seqtab_file)
cat('ASVs to classify:', ncol(seqtab_nochim), '\n')

# The reference must match the marker AND ideally the primer region. A FULL-LENGTH SILVA training
# set applied to V4 (~250 bp) reads mismatches k-mer composition and degrades calls (Werner 2012;
# Bokulich 2018). For V4 data prefer a region-extracted reference (see assign_qiime2_region.sh).
minBoot <- 50  # DADA2 default and the RDP recommendation for reads <=250 nt; tutorials often use
               # 80 (a stricter CHOICE). Raising it truncates to shallower-but-reliable ranks;
               # ranks below the floor are returned as NA, not guessed.

# assignTaxonomy() runs 100 stochastic bootstrap resamples per sequence; without a seed, two runs
# on identical input produce different genus calls for a real fraction of ASVs. set.seed() before
# the call makes the GENUS calls reproducible (any fixed integer works; 100 is just a convention
# here). It is not bitwise identity at every rank: a residual of ~2/770 cells at Kingdom/Order can
# remain between seeded runs, with multithread=TRUE and also with multithread=FALSE.
set.seed(100)
taxa <- assignTaxonomy(seqtab_nochim, silva_train, minBoot = minBoot, tryRC = TRUE, multithread = TRUE)

# Sanity check before trusting the table: a wrong-region/wrong-marker reference, bad orientation or
# a mis-formatted training FASTA gives all-NA calls with no error or warning.
genus_frac <- mean(!is.na(taxa[, 'Genus']))
if (genus_frac == 0) stop('assignTaxonomy returned no genus calls: check the reference format, region and marker')
cat(sprintf('Genus assigned: %.1f%% of ASVs
', 100 * genus_frac))
if (rerun_check) {  # optional determinism check: same seed, same input, genus calls must match
    set.seed(100)
    again <- assignTaxonomy(seqtab_nochim, silva_train, minBoot = minBoot, tryRC = TRUE, multithread = TRUE)
    n_diff <- sum(is.na(taxa[, 'Genus']) != is.na(again[, 'Genus']) |
                  (!is.na(taxa[, 'Genus']) & !is.na(again[, 'Genus']) & taxa[, 'Genus'] != again[, 'Genus']))
    cat(sprintf('Rerun check: %d genus calls differ between two seeded runs
', n_diff))
    if (n_diff > 0) stop('genus calls are not reproducible between seeded runs')
}

# addSpecies assigns species ONLY by exact (100%) match against the species reference - it does
# not infer species from a noisy read. Most ASVs stay NA at species; that is the honest 16S result.
if (file.exists(silva_species)) {
    taxa <- addSpecies(taxa, silva_species)
}

cat('\nAssigned fraction per rank (NA = unassigned at threshold, kept honestly):\n')
for (rank in colnames(taxa)) {
    assigned <- sum(!is.na(taxa[, rank]))
    cat(sprintf('  %-8s %d/%d (%.1f%%)\n', rank, assigned, nrow(taxa), 100 * assigned / nrow(taxa)))
}

# State the conditioning choices alongside the labels (classifier + database release + region).
attr(taxa, 'classifier') <- 'DADA2 assignTaxonomy (RDP naive Bayes)'
attr(taxa, 'reference') <- 'SILVA 138.1'
attr(taxa, 'region') <- 'verify the reference region matches the amplicon primers'

saveRDS(taxa, out_file)
cat('\nSaved', out_file, '(genus via naive Bayes, species via exact match only)\n')
