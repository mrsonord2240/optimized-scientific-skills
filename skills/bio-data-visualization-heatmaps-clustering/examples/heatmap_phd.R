# Reference: ComplexHeatmap 2.18+, circlize 0.4.16+, seriation 1.5+ | Verify API if version differs

# PhD-level annotated heatmap with the four correctness traps encoded:
# (1) ward.D2 NOT ward.D, (2) robust symmetric color bounds, (3) OLO for leaf ordering,
# (4) draw() not bare Heatmap() for non-interactive use. With an explicit
# OLO dendrogram, keep pathway labels in a row annotation rather than row_split:
# ComplexHeatmap 2.22.0 cannot combine a categorical row_split with a supplied dendrogram.

library(ComplexHeatmap)
library(circlize)
library(seriation)

# 1. INPUT -- self-contained features (rows) x samples (columns) demonstration.
# Replace this deterministic toy matrix with a normalized expression matrix in practice.
set.seed(20260923)
n_genes <- 36
n_samples <- 12
pathway_levels <- c('Metabolism', 'Signaling', 'CellCycle')
gene_info <- data.frame(
    pathway = factor(rep(pathway_levels, each = n_genes / length(pathway_levels)),
                     levels = pathway_levels),
    log2FC = c(rnorm(12, 0.8, 0.25), rnorm(12, 0, 0.25), rnorm(12, -0.4, 0.25)),
    row.names = paste0('Gene', seq_len(n_genes))
)
metadata <- data.frame(
    condition = factor(rep(c('Control', 'Treatment'), each = n_samples / 2),
                       levels = c('Control', 'Treatment')),
    batch = factor(rep(c('A', 'B', 'C'), length.out = n_samples), levels = c('A', 'B', 'C')),
    age = seq(35, 35 + 2 * (n_samples - 1), by = 2),
    row.names = paste0('Sample', seq_len(n_samples))
)
mat <- matrix(rnorm(n_genes * n_samples, sd = 0.35), nrow = n_genes,
              dimnames = list(rownames(gene_info), rownames(metadata)))
pathway_profiles <- rbind(
    Metabolism = seq(-0.6, 0.6, length.out = n_samples),
    Signaling = c(rep(-0.5, n_samples / 2), rep(0.8, n_samples / 2)),
    CellCycle = sin(seq(0, 2 * pi, length.out = n_samples))
)
for (pathway in pathway_levels) {
    rows <- gene_info$pathway == pathway
    mat[rows, ] <- mat[rows, ] + rep(pathway_profiles[pathway, ], each = sum(rows))
}

# Assume row z-score is requested (typical bulk RNA-seq case).
mat_scaled <- t(scale(t(mat)))                          # row z-score
mat_scaled[is.na(mat_scaled)] <- 0                      # rows with sd=0 become NaN

# 2. ROBUST SYMMETRIC COLOR BOUNDS -- 1st-99th percentile of |z|
bounds <- as.numeric(quantile(abs(mat_scaled), 0.99, na.rm = TRUE))
col_fun <- colorRamp2(c(-bounds, 0, bounds),
                      c('#0072B2', 'white', '#D55E00'))  # CVD-safe diverging

# 3. CLUSTERING -- ward.D2 not ward.D (Murtagh-Legendre 2014)
d_rows <- dist(mat_scaled, method = 'euclidean')
hc_rows <- hclust(d_rows, method = 'ward.D2')

# 4. OPTIMAL LEAF ORDERING (Bar-Joseph 2001) -- reveals block structure
olo_rows <- seriate(d_rows, method = 'OLO', control = list(hclust = hc_rows))
dend_rows <- as.dendrogram(olo_rows[[1]])

# 5. COLUMN ANNOTATION -- Okabe-Ito categorical palette (Wong 2011)
ha_col <- HeatmapAnnotation(
    Condition = metadata$condition,
    Batch     = metadata$batch,
    Age       = anno_barplot(metadata$age, gp = gpar(fill = '#56B4E9')),
    col = list(
        Condition = c(Control = '#56B4E9', Treatment = '#D55E00'),
        Batch     = c(A = '#009E73', B = '#0072B2', C = '#CC79A7')),
    annotation_name_gp = gpar(fontsize = 8),
    show_legend = TRUE)

# 6. ROW ANNOTATION
ha_row <- rowAnnotation(
    Pathway = gene_info$pathway,
    LogFC   = anno_barplot(gene_info$log2FC, baseline = 0,
                            gp = gpar(fill = ifelse(gene_info$log2FC > 0,
                                                     '#D55E00', '#0072B2'))),
    col = list(Pathway = setNames(c('#8491B4', '#91D1C2', '#F39B7F'), pathway_levels)))
stopifnot(identical(ha_row@anno_list$Pathway@color_mapping@levels, pathway_levels))

# 7. ASSEMBLE
ht <- Heatmap(mat_scaled,
              name = 'Z-score',
              col  = col_fun,
              cluster_rows    = dend_rows,                # OLO dendrogram
              cluster_columns = FALSE,                   # retain input order within column_split groups
              clustering_method_columns   = 'ward.D2',    # explicit; never 'ward'
              clustering_distance_columns = 'euclidean',
              top_annotation  = ha_col,
              left_annotation = ha_row,                    # Pathway labels remain visible
              column_split = metadata$condition,          # grouped column layout
              show_row_names = FALSE,
              show_column_names = TRUE,
              column_names_gp  = gpar(fontsize = 7),
              use_raster = TRUE,                          # raster the cell layer
              raster_quality = 5,                         # publication quality (default 1 is pixelated)
              heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 8, fontface = 'bold'),
                  labels_gp = gpar(fontsize = 7),
                  at = c(-bounds, 0, bounds)))

# 8. RENDER -- draw() NOT bare Heatmap() in a script
pdf('heatmap.pdf', width = 7, height = 9)
ht_drawn <- draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right',
     padding = unit(c(2, 2, 2, 2), 'mm'))
dev.off()

# 9. EXTRACT CLUSTER ASSIGNMENTS for downstream use; reuse the rendered object.
row_order_vector <- row_order(ht_drawn)                  # the OLO order for all rows
column_order_list <- column_order(ht_drawn)
stopifnot(identical(as.integer(row_order_vector), as.integer(order.dendrogram(dend_rows))))
# To get k=4 cuts of the row dendrogram:
row_clusters <- cutree(as.hclust(dend_rows), k = 4)
