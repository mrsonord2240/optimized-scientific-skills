args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript summarize_units.R observations.csv output.csv")
d <- read.csv(args[1], stringsAsFactors = FALSE)
if (!all(c("unit", "condition", "measurement") %in% names(d)) ||
    anyNA(d$unit) || anyNA(d$condition) || !is.numeric(d$measurement) ||
    any(!is.finite(d$measurement))) stop("Invalid unit, condition or measurement")
if (any(vapply(split(d$condition, d$unit), function(x) length(unique(x)) != 1, logical(1)))) stop("A unit has more than one condition")
means <- aggregate(measurement ~ unit + condition, data = d, FUN = mean)
counts <- aggregate(measurement ~ unit + condition, data = d, FUN = length)
names(counts)[3] <- "observations"
write.csv(merge(means, counts, by = c("unit", "condition")), args[2], row.names = FALSE)
