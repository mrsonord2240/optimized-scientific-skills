args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: Rscript randomize_blocks.R units.csv output.csv label1,label2 seed")
d <- read.csv(args[1], stringsAsFactors = FALSE)
labels <- strsplit(args[3], ",", fixed = TRUE)[[1]]
seed <- as.integer(args[4])
if (!all(c("id", "block") %in% names(d)) || anyNA(d$id) || anyNA(d$block) ||
    anyDuplicated(d$id) || length(labels) != 2 || any(!nzchar(labels)) || anyDuplicated(labels) || is.na(seed)) stop("Invalid IDs, blocks, labels or seed")
set.seed(seed)
d$treatment <- NA_character_
d$run_order <- NA_integer_
for (b in unique(d$block)) {
  i <- which(d$block == b)
  assignments <- rep(labels, length.out = length(i))
  d$treatment[i] <- assignments[sample.int(length(i))]
  d$run_order[i] <- sample.int(length(i))
}
print(table(d$block, d$treatment))
write.csv(d, args[2], row.names = FALSE)
