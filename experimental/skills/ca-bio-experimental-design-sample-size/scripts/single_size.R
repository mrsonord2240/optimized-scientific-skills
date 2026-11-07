args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript single_size.R scenario.csv output.rds")
library(ssizeRNA)
d <- read.csv(args[1], check.names = FALSE)
cols <- c("nGenes", "pi0", "m", "mu", "disp", "fc", "fdr", "power", "maxN", "seed")
if (nrow(d) != 1 || !all(cols %in% names(d))) stop("Need one row and all required columns")
x <- as.list(d[1, cols])
if (any(!is.finite(unlist(x))) || x$nGenes < 1 || x$m < 1 || x$mu <= 0 || x$disp <= 0 ||
    x$fc <= 0 || x$fc == 1 || x$pi0 <= 0 || x$pi0 > 1 || x$fdr <= 0 || x$fdr >= 1 ||
    x$power <= 0 || x$power >= 1 || x$maxN < 2) stop("Invalid scenario parameters")
set.seed(x$seed)
result <- tryCatch(ssizeRNA_single(nGenes = x$nGenes, pi0 = x$pi0, m = x$m, mu = x$mu,
  disp = x$disp, fc = x$fc, fdr = x$fdr, power = x$power, maxN = x$maxN),
  error = function(e) {
    if (identical(conditionMessage(e), "argument is of length zero")) return(NULL)
    stop(e)
  })
if (is.null(result) || length(result$ssize[, "ssize"]) != 1 || is.na(result$ssize[, "ssize"]))
  stop(sprintf("No n <= %s returned; revise assumptions or search ceiling", x$maxN))
saveRDS(list(input = d, result = result), args[2])
