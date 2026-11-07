#!/usr/bin/env Rscript
# data.csv exposure mediator outcome covariates(comma-separated or -) output.rds
args <- commandArgs(TRUE)
if (length(args) != 6L) stop('usage: cmaverse_4way.R data.csv exposure mediator outcome covariates|- output.rds')
suppressPackageStartupMessages(library(CMAverse))
dat <- read.csv(args[1], check.names=FALSE)
covars <- if (args[5] == '-') character() else strsplit(args[5], ',', fixed=TRUE)[[1]]
if (!all(c(args[2:4],covars) %in% names(dat))) stop('missing requested column')
if (anyNA(dat[, c(args[2:4],covars), drop=FALSE])) stop('handle missing data before fitting')
result <- cmest(data=dat, model='rb', outcome=args[4], exposure=args[2], mediator=args[3],
  basec=covars, EMint=TRUE, mreg=list('linear'), yreg='logistic', astar=0, a=1, mval=list(0),
  estimation='paramfunc', inference='bootstrap', nboot=1000, boot.ci.type='per')
saveRDS(result,args[6])
