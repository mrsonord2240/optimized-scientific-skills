# Reference: DECIPHER 2.30+ (checked on DECIPHER 3.2.0, R 4.4.3) | Verify API if version differs
# Classify ASVs with DECIPHER IDTAXA and flatten the result to an ASV x rank matrix (NA where unclassified).
# Inputs : seqtab.rds (DADA2 seqtab_nochim, ASV sequences as column names)
#          a trainingSet: a pre-trained .RData (provides the `trainingSet` object) or .rds, OR train one
#          with LearnTaxa() from a reference FASTA + "Root;domain;...;genus;" taxonomy strings.
# Output : out.rds (character matrix, rows = ASV sequences, columns = domain..species)
# Usage  : Rscript idtaxa_classify.R seqtab.rds trainingset.(RData|rds) out.rds [rerun]
#          Rscript idtaxa_classify.R seqtab.rds trainingset.rds out.rds refseqs.fasta reftax.txt [rerun]
#            (5-argument form: trains with LearnTaxa() and saves the trainingSet to trainingset.rds first)
#          rerun = repeat the seeded IdTaxa() call once and require identical() output.
library(DECIPHER)

args <- commandArgs(trailingOnly = TRUE)
rerun_check <- 'rerun' %in% args
args <- setdiff(args, 'rerun')
stopifnot(length(args) %in% c(3, 5))
seqtab_file <- args[1]
trainingset_file <- args[2]
out_file <- args[3]
refseqs_file <- if (length(args) == 5) args[4] else NULL
reftax_file <- if (length(args) == 5) args[5] else NULL

seqtab_nochim <- readRDS(seqtab_file)

if (is.null(refseqs_file)) {
    # A pre-trained trainingSet, if you have one: .RData provides the trainingSet object.
    if (grepl('\\.rds$', trainingset_file, ignore.case = TRUE)) {
        trainingSet <- readRDS(trainingset_file)
    } else {
        load(trainingset_file)
    }
} else {
    # No pre-trained .RData for your marker/region? Train one directly from a reference FASTA +
    # matching "Root;domain;phylum;...;genus;" taxonomy strings (one per sequence, same order).
    # LearnTaxa() tunes its tree-descent k-mer sampling with repeated random subsamples (its own
    # documentation: "this process is repeated with 100 random subsamples") -- inherently stochastic,
    # same class of bug as IdTaxa() below. Verified: two unseeded LearnTaxa() calls on identical input
    # produce non-identical trainingSet objects; set.seed() before EVERY LearnTaxa() call makes the
    # trainingSet object itself reproducible (identical() TRUE).
    refseqs <- readDNAStringSet(refseqs_file)
    reftax  <- readLines(reftax_file)  # e.g. "Root;Bacteria;Firmicutes;...;"
    set.seed(100)
    trainingSet <- LearnTaxa(refseqs, taxonomy = reftax)
    saveRDS(trainingSet, trainingset_file)
    # MEMORY: LearnTaxa() against a full, un-subsampled reference (400K+ sequences) needs tens of GB
    # of RAM and can crash on constrained hardware; subsample the reference (e.g. ~60,000 sequences)
    # if it does.
    # LearnTaxa's OPTIONAL rank= argument (a 5-column Index/Name/Parent/Level/Rank data.frame, rarely
    # available outside DECIPHER's own pre-built .RData sets) is not required to train or classify --
    # see the flattening note below for why it matters anyway.
}

dna <- DNAStringSet(colnames(seqtab_nochim))

# IdTaxa() descends its classification tree with an internal stochastic step -- inherently
# stochastic, same as assignTaxonomy(), and by a LARGER margin (verified: unseeded, two
# back-to-back calls on the identical trainingSet and identical query set differ at ~3-4% of
# genus calls). set.seed() before EVERY IdTaxa() call -- without it, repeated runs on the same
# input differ at the genus call for ~3-4% of ASVs. Verified: with set.seed() before each call,
# repeated runs are bit-identical at every rank including genus, in both the default
# multithreaded (processors=NULL) and single-threaded (processors=1) configurations; any fixed
# integer works, 100 is just a convention here (matches the assignTaxonomy() seed).
set.seed(100)

# threshold 60 = DECIPHER default confidence cutoff; raise for stricter calls. IDTAXA's
# tree-descent stops (leaves the rank unclassified) when the query likely belongs to a taxon
# absent from the reference - this is the intended anti-over-classification behaviour.
ids <- IdTaxa(dna, trainingSet, strand = 'both', threshold = 60, processors = NULL)

ranks <- c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')

# Flatten POSITIONALLY, not by name. x$rank is populated ONLY when trainingSet was built with
# LearnTaxa's rank= data.frame (see above) -- absent that, x$rank is NULL for every result, and
# match(ranks, x$rank) silently returns all-NA with no error or warning (confirmed empirically
# against a real LearnTaxa()-trained set: 100% NA at every rank, no exception raised). x$taxon[1]
# is always "Root"; the remaining entries are domain..genus/species in taxonomic order regardless
# of whether rank= was supplied, so index positionally instead.
taxa_idtaxa <- t(sapply(ids, function(x) {
    taxa <- x$taxon[-1]                    # drop "Root"
    taxa[startsWith(taxa, 'unclassified_')] <- NA
    length(taxa) <- length(ranks)          # pad/truncate to the fixed rank depth above
    taxa
}))
colnames(taxa_idtaxa) <- ranks
rownames(taxa_idtaxa) <- colnames(seqtab_nochim)

# Sanity check: a wrong-region or wrong-marker training set, or a flattening bug, gives all-NA
# calls with no error.
if (all(is.na(taxa_idtaxa[, 'genus']))) stop('IdTaxa returned no genus calls: check the trainingSet marker/region and the flattening')
cat(sprintf('Genus assigned: %d/%d ASVs\n', sum(!is.na(taxa_idtaxa[, 'genus'])), nrow(taxa_idtaxa)))

if (rerun_check) {  # optional determinism check: same seed, same input, output must be identical
    set.seed(100)
    ids_again <- IdTaxa(dna, trainingSet, strand = 'both', threshold = 60, processors = NULL)
    if (!identical(ids, ids_again)) stop('IdTaxa output is not reproducible between seeded runs')
    cat('Rerun check: identical() TRUE between two seeded IdTaxa() runs\n')
}

saveRDS(taxa_idtaxa, out_file)
cat('Saved', out_file, '\n')
