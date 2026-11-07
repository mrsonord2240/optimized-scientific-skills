#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('usage: supplementary.R GSE_ACCESSION fresh_output_directory')
if (file.exists(args[2])) stop('output directory exists; choose a fresh path')
dir.create(args[2], recursive=TRUE)
suppressPackageStartupMessages(library(GEOquery))
files <- getGEOSuppFiles(args[1], baseDir=args[2])
write.csv(files, file.path(args[2], 'supplementary_manifest.csv'))
