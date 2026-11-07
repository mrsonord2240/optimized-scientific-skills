args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript sva_from_files.R expression.csv metadata.csv output.rds")
library(sva)
x <- read.csv(args[1], row.names = 1, check.names = FALSE)
meta <- read.csv(args[2], stringsAsFactors = FALSE)
if (!all(c("id", "condition") %in% names(meta)) || anyDuplicated(meta$id) ||
    !setequal(colnames(x), meta$id)) stop("Expression columns and unique metadata IDs must match")
meta <- meta[match(colnames(x), meta$id), , drop = FALSE]
expr <- as.matrix(x)
if (!is.numeric(expr)) stop("Expression matrix must be numeric")
good <- apply(expr, 1, function(row) all(is.finite(row)))
message(sprintf("Removed %d of %d non-finite features", sum(!good), length(good)))
expr <- expr[good, , drop = FALSE]
if (!nrow(expr)) stop("No complete features remain")
mod <- model.matrix(~ condition, data = meta)
mod0 <- model.matrix(~ 1, data = meta)
n_sv <- num.sv(expr, mod, method = "leek")
sv <- if (n_sv > 0) sva(expr, mod, mod0, n.sv = n_sv) else NULL
saveRDS(list(n_sv = n_sv, sva = sv, removed_features = sum(!good)), args[3])
