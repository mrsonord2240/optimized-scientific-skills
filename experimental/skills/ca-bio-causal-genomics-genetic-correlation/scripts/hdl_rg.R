#!/usr/bin/env Rscript
# Usage: Rscript hdl_rg.R trait1.tsv trait2.tsv ld_reference output.txt overlap_count
args <- commandArgs(TRUE)
if (length(args) != 5) stop('usage: hdl_rg.R trait1.tsv trait2.tsv ld_reference output.txt overlap_count')
suppressPackageStartupMessages(library(HDL))
x <- read.table(args[1], header=TRUE, stringsAsFactors=FALSE)
y <- read.table(args[2], header=TRUE, stringsAsFactors=FALSE)
overlap <- as.numeric(args[5])
if (!is.finite(overlap) || overlap < 0) stop('provide known nonnegative overlap count')
result <- HDL.rg(gwas1.df=x, gwas2.df=y, LD.path=args[3], N0=overlap, output.file=args[4])
print(result)
