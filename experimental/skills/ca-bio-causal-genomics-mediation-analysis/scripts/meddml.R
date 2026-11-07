#!/usr/bin/env Rscript
# data.csv outcome treatment mediator covariates(comma-separated) output.rds
args <- commandArgs(TRUE)
if (length(args) != 6L) stop('usage: meddml.R data.csv outcome treatment mediator covariates output.rds')
suppressPackageStartupMessages(library(causalweight))
dat <- read.csv(args[1], check.names=FALSE)
covars <- strsplit(args[5], ',', fixed=TRUE)[[1]]
if (!length(covars) || !all(c(args[2:4],covars) %in% names(dat))) stop('missing requested columns')
if (anyNA(dat[, c(args[2:4],covars), drop=FALSE])) stop('handle missing data before fitting')
x <- as.matrix(dat[,covars,drop=FALSE])
if (!is.numeric(x)) stop('covariates must be numeric and named')
result <- medDML(y=dat[[args[2]]], d=dat[[args[3]]], m=dat[[args[4]]], x=x, trim=0.05, order=1)
saveRDS(result,args[6])
