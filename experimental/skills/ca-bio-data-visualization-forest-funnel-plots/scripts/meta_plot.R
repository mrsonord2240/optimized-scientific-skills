args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('usage: meta_plot.R studies.csv output-directory')
library(metafor)
studies <- read.csv(args[1], stringsAsFactors = FALSE)
if (!all(c('author', 'year', 'log_or', 'log_or_se') %in% names(studies))) stop('missing columns')
if (nrow(studies) < 3) stop('do not pool fewer than three studies')
if (anyNA(studies[, c('author', 'year', 'log_or', 'log_or_se')]) ||
    any(!is.finite(studies$log_or)) || any(!is.finite(studies$log_or_se)) ||
    any(studies$log_or_se <= 0)) stop('invalid or missing study values')
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
res <- rma(yi = log_or, vi = log_or_se^2, data = studies,
           slab = paste(author, year), method = 'REML',
           test = if (nrow(studies) < 5) 'knha' else 'z')
pred <- predict(res)
ticks <- log(c(0.25, 0.5, 1, 2, 4))
limits <- range(c(studies$log_or - qnorm(.975) * studies$log_or_se,
                  studies$log_or + qnorm(.975) * studies$log_or_se,
                  res$ci.lb, res$ci.ub, pred$pi.lb, pred$pi.ub, ticks), finite = TRUE)
label <- if (res$k < 5) sprintf('RE model (k=%d; HKSJ; heterogeneity withheld)', res$k) else
    sprintf('RE model (tau^2=%.3f; I^2=%.1f%%; Q p=%s)', res$tau2, res$I2, format.pval(res$QEp, digits=2))
pdf(file.path(args[2], 'forest.pdf'), width=8, height=6)
forest(res, atransf=exp, at=ticks, alim=limits, refline=0, xlab='Odds ratio (95% CI)', mlab=label, addpred=TRUE)
dev.off()
pdf(file.path(args[2], 'funnel.pdf'), width=7, height=5)
funnel(res, level=c(90,95,99), shade=c('white','gray55','gray75'), refline=0, legend=TRUE, xlab='log(OR)')
dev.off()
sink(file.path(args[2], 'summary.txt'))
cat('k:', res$k, ' pooled log OR:', as.numeric(res$b[1]), ' CI:', res$ci.lb, res$ci.ub, '\n')
if (res$k >= 5) cat('tau-squared:', res$tau2, ' I-squared:', res$I2, ' Q p:', res$QEp, '\n') else cat('Heterogeneity metrics withheld for k < 5\n')
print(pred)
if (res$k >= 10) print(regtest(res, model='lm', predictor='sei')) else cat('Egger withheld: k < 10\n')
tf <- trimfill(res)
cat('Trim-and-fill sensitivity only; imputed studies:', tf$k-res$k, '\n')
cat('Original OR:', exp(res$b[1]), ' sensitivity OR:', exp(tf$b[1]), '\n')
sink()
