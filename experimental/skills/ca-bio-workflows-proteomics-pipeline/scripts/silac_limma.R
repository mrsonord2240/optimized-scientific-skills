# silac_limma.R -- MaxQuant SILAC H/L ratios -> moderated one-sample limma against log2 ratio 0.
# Purpose : the SILAC route of bio-workflows-proteomics-pipeline. Caveats (Arg->Pro conversion, labelling
#           efficiency) are in references/silac.md.
# Inputs  : proteinGroups.txt with `Ratio H/L normalized <sample>` columns
# Usage   : Rscript silac_limma.R [proteinGroups.txt] [out.csv]
# Output  : out.csv = topTable (protein, logFC, AveExpr, t, P.Value, adj.P.Val BH, B)
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
input_file <- arg(1, 'proteinGroups.txt')
out_file   <- arg(2, 'silac_results.csv')

library(limma)

# SILAC ratios from MaxQuant (quote/comment.char as in scripts/msstats_maxquant.R)
silac <- read.delim(input_file, quote = '', comment.char = '')
ratio_cols <- grep('Ratio.H.L.normalized', colnames(silac), value = TRUE)

# Log2 transform ratios. MaxQuant leaves NaN (and 0 for an absent channel) wherever it could not
# form a ratio, so log2 produces NaN/-Inf; coerce every non-finite cell to NA before testing.
silac_log2 <- log2(as.matrix(silac[, ratio_cols]))
silac_log2[!is.finite(silac_log2)] <- NA
rownames(silac_log2) <- silac$Majority.protein.IDs   # without this the result table has no identity

# Keep proteins with >= 2 finite ratios. apply(t.test) over the raw matrix STOPS the whole script
# with "not enough 'x' observations" on the first protein quantified in one replicate only.
keep <- rowSums(!is.na(silac_log2)) >= 2
cat('proteins tested:', sum(keep), 'of', nrow(silac_log2), '\n')

# One-sample moderated test against log2 ratio 0 (no change): an intercept-only limma fit borrows
# variance across proteins, which an unmoderated per-protein t-test at n = 3 cannot. Report the
# BH-adjusted p-value -- a raw p-value per protein controls nothing across thousands of tests.
fit <- eBayes(lmFit(silac_log2[keep, , drop = FALSE]), trend = TRUE, robust = TRUE)
results <- topTable(fit, coef = 1, number = Inf, adjust.method = 'BH')
# Proteins dropped by `keep` are an on/off presence table, not a fold change (quantification).
write.csv(cbind(protein = rownames(results), results), out_file, row.names = FALSE)
cat('wrote', out_file, '\n')
