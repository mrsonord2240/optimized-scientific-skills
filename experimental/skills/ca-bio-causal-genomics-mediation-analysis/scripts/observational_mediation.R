#!/usr/bin/env Rscript
# data.csv exposure mediator outcome covariates(comma-separated or -) output.rds outcome_family
args <- commandArgs(TRUE)
if (length(args) != 7L) stop('usage: observational_mediation.R data.csv exposure mediator outcome covariates|- output.rds linear|probit')
suppressPackageStartupMessages(library(mediation))
dat <- read.csv(args[1], check.names=FALSE)
exposure <- args[2]; mediator <- args[3]; outcome <- args[4]
covars <- if (args[5] == '-') character() else strsplit(args[5], ',', fixed=TRUE)[[1]]
if (!all(c(exposure,mediator,outcome,covars) %in% names(dat))) stop('missing requested column')
if (anyNA(dat[, c(exposure,mediator,outcome,covars), drop=FALSE])) stop('handle missing data before fitting')
mf <- reformulate(c(exposure,covars), response=mediator)
of <- reformulate(c(exposure,mediator,covars), response=outcome)
med_model <- eval(bquote(lm(.(mf), data=.(dat))))
if (args[7] == 'linear') {
  out_model <- eval(bquote(lm(.(of), data=.(dat))))
} else if (args[7] == 'probit') {
  out_model <- eval(bquote(glm(.(of), data=.(dat), family=binomial(link='probit'))))
} else stop('outcome family must be linear or probit')
result <- mediate(med_model, out_model, treat=exposure, mediator=mediator, boot=TRUE, sims=1000)
sensitivity <- medsens(result, rho.by=0.05, effect.type='indirect', sims=1000)
saveRDS(list(mediation=result, sensitivity=sensitivity), args[6])
