# harmonise.R -- allele harmonisation of two summary-statistic tables before coloc.
#
# Purpose: merge two datasets by SNP, resolve same / flip / strand-complement allele coding, negate
#          betas on flips, and drop unresolvable and high-MAF palindromic SNPs (MAF > 0.42).
# Inputs:  df1, df2 with columns SNP, A1, A2, BETA, MAF (extra columns are kept, suffixed .1/.2).
#          Genome builds must already match (lift over first).
# Output:  merged data.frame of kept SNPs (columns suffixed .1 / .2); in file mode, a TSV.
# Usage:   source('scripts/harmonise.R'); m <- harmonise(df1, df2)
#          Rscript scripts/harmonise.R gwas.tsv eqtl.tsv harmonised.tsv

harmonise <- function(df1, df2) {
    comp <- function(a) c(A='T', T='A', C='G', G='C')[a]
    m <- merge(df1, df2, by='SNP', suffixes=c('.1','.2'))
    palindromic <- (m$A1.1 %in% c('A','T') & m$A2.1 %in% c('A','T')) |
                   (m$A1.1 %in% c('C','G') & m$A2.1 %in% c('C','G'))
    # Non-palindromic strand mismatch: complement dataset 2's alleles, then treat as same/flip
    strand <- !palindromic & m$A1.1 == comp(m$A1.2) & m$A2.1 == comp(m$A2.2) |
              !palindromic & m$A1.1 == comp(m$A2.2) & m$A2.1 == comp(m$A1.2)
    strand[is.na(strand)] <- FALSE
    a1 <- ifelse(strand, comp(m$A1.2), m$A1.2); a2 <- ifelse(strand, comp(m$A2.2), m$A2.2)
    m$A1.2 <- unname(a1); m$A2.2 <- unname(a2)
    same <- m$A1.1 == m$A1.2 & m$A2.1 == m$A2.2
    flip <- m$A1.1 == m$A2.2 & m$A2.1 == m$A1.2
    m$BETA.2[flip] <- -m$BETA.2[flip]
    m$MAF.2[flip] <- 1 - m$MAF.2[flip]
    keep <- (same | flip) & !(palindromic & m$MAF.1 > 0.42)
    m[keep, ]
}

if (sys.nframe() == 0) {   # run as a script, not sourced
    args <- commandArgs(TRUE)
    stopifnot(length(args) == 3)
    out <- harmonise(read.delim(args[1], stringsAsFactors = FALSE),
                     read.delim(args[2], stringsAsFactors = FALSE))
    write.table(out, args[3], sep = '\t', quote = FALSE, row.names = FALSE)
    cat('Harmonised', nrow(out), 'SNPs ->', args[3], '\n')
}
