args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('Usage: Rscript saturation.R alignment.fasta pairs.tsv')
library(ape)
aln <- read.dna(args[1], format = 'fasta')
L <- ncol(aln)
p <- dist.dna(aln, model = 'raw')
d <- dist.dna(aln, model = 'TN93')
ts <- dist.dna(aln, model = 'TS') / L
tv <- dist.dna(aln, model = 'TV') / L
pairs <- data.frame(p = as.vector(p), tn93 = as.vector(d), transitions = as.vector(ts), transversions = as.vector(tv))
write.table(pairs, args[2], sep = '\t', quote = FALSE, row.names = FALSE)
cat('fraction p > 0.5:', mean(p > 0.5), '\n')
cat('transition vs distance slope:', coef(lm(transitions ~ tn93, data = pairs))[2], '\n')
