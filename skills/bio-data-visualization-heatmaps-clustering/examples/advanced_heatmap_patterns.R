# Self-contained patterns for three decisions that are easy to get wrong:
# sample-correlation QC, shared row ordering, and single-cell pseudobulk.

library(ComplexHeatmap)
library(circlize)

set.seed(20260923)

# 1. SAMPLE-CORRELATION QC -----------------------------------------------
# Select features before looking at the QC clustering. Here the objective
# variance rule is deliberately independent of sample labels and dendrograms.
expression <- matrix(rnorm(36 * 9, sd = 0.5), nrow = 36,
                     dimnames = list(paste0("gene", 1:36), paste0("sample", 1:9)))
expression[1:12, 4:9] <- expression[1:12, 4:9] + 1.2
top_variable <- order(apply(expression, 1, var), decreasing = TRUE)[1:20]
sample_correlation <- cor(expression[top_variable, , drop = FALSE])
sample_distance <- as.dist(1 - sample_correlation)
sample_dendrogram <- as.dendrogram(hclust(sample_distance, method = "average"))
correlation_range <- range(sample_correlation, na.rm = TRUE)
correlation_col <- colorRamp2(correlation_range, c("#F7FBFF", "#2166AC"))
stopifnot(isTRUE(all.equal(as.matrix(sample_distance), 1 - sample_correlation)))

pdf("sample_correlation_qc.pdf", width = 5, height = 5)
draw(Heatmap(sample_correlation, name = "Pearson r", col = correlation_col,
             cluster_rows = sample_dendrogram, cluster_columns = sample_dendrogram,
             heatmap_legend_param = list(at = correlation_range)))
dev.off()

# 2. LINKED HEATMAPS SHARE THE EXPRESSION-DERIVED ROW ORDER --------------
expression_z <- t(scale(t(expression)))
expression_z[is.na(expression_z)] <- 0
expression_dendrogram <- as.dendrogram(hclust(dist(expression_z), method = "ward.D2"))
shared_row_order <- order.dendrogram(expression_dendrogram)
methylation <- matrix(runif(nrow(expression) * ncol(expression), 0.2, 0.8),
                      nrow = nrow(expression), dimnames = dimnames(expression))

expression_ht <- Heatmap(expression_z, name = "Expression z", cluster_rows = FALSE,
                         row_order = shared_row_order, cluster_columns = FALSE)
methylation_ht <- Heatmap(methylation, name = "Methylation beta", cluster_rows = FALSE,
                          row_order = shared_row_order, cluster_columns = FALSE)
stopifnot(identical(as.integer(expression_ht@row_order), as.integer(shared_row_order)))
stopifnot(identical(as.integer(methylation_ht@row_order), as.integer(shared_row_order)))

pdf("shared_row_order_heatmaps.pdf", width = 8, height = 7)
draw(expression_ht + methylation_ht)
dev.off()

# 3. PSEUDOBULK BEFORE GROUP-LEVEL HEATMAPPING ---------------------------
cell_counts <- matrix(rpois(24 * 18, lambda = 4), nrow = 24,
                      dimnames = list(paste0("gene", 1:24), paste0("cell", 1:18)))
cell_group <- factor(rep(c("B_cell", "T_cell", "Myeloid"), each = 6),
                     levels = c("B_cell", "T_cell", "Myeloid"))
pseudobulk <- t(rowsum(t(cell_counts), group = cell_group, reorder = FALSE))
stopifnot(identical(colnames(pseudobulk), levels(cell_group)))
stopifnot(ncol(pseudobulk) == nlevels(cell_group))

pseudobulk_z <- t(scale(t(log1p(pseudobulk))))
pseudobulk_z[is.na(pseudobulk_z)] <- 0
pdf("pseudobulk_heatmap.pdf", width = 5, height = 7)
draw(Heatmap(pseudobulk_z, name = "Pseudobulk z", clustering_method_rows = "ward.D2",
             cluster_columns = FALSE, column_split = cell_group[match(colnames(pseudobulk), levels(cell_group))]))
dev.off()

message("PASS: correlation QC, shared row order, and pseudobulk heatmap patterns")
