args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript count_power.R scenarios.csv output.csv")
library(RNASeqPower)
d <- read.csv(args[1], check.names = FALSE)
required <- c("depth", "cv", "effect", "n", "power", "alpha")
if (!all(required %in% names(d))) stop("Missing required scenario columns")
validate <- function(x) {
  if (any(!is.finite(unlist(x[c("depth", "cv", "alpha")])))) stop("depth, cv, alpha must be finite")
  if (x$depth <= 0 || x$cv <= 0 || x$alpha <= 0 || x$alpha >= 1) stop("invalid depth, cv or alpha")
  if (sum(is.na(unlist(x[c("effect", "n", "power")]))) != 1) stop("exactly one of effect, n and power must be missing")
  if (!is.na(x$effect) && (!is.finite(x$effect) || x$effect <= 0 || x$effect == 1)) stop("invalid effect")
  if (!is.na(x$n) && (!is.finite(x$n) || x$n <= 0)) stop("invalid n")
  if (!is.na(x$power) && (!is.finite(x$power) || x$power <= 0 || x$power >= 1)) stop("invalid power")
}
d$result <- vapply(seq_len(nrow(d)), function(i) {
  x <- as.list(d[i, required]); validate(x)
  a <- x[!vapply(x, is.na, logical(1))]
  as.numeric(do.call(RNASeqPower::rnapower, a))
}, numeric(1))
write.csv(d, args[2], row.names = FALSE)
