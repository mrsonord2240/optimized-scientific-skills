# Reference: ComplexHeatmap 2.18+, circlize 0.4.16+, RColorBrewer 1.1+ | Verify API if version differs
library(ComplexHeatmap)
library(circlize)

# --- ALTERNATIVE: Use real Bioconductor datasets ---
# For realistic expression heatmaps, try these publicly available datasets:
#
# library(airway)
# data('airway')
# vst_data <- assay(vst(airway))
# top_var <- head(order(rowVars(vst_data), decreasing = TRUE), 100)
# mat <- t(scale(t(vst_data[top_var, ])))  # Z-score normalize
# metadata <- as.data.frame(colData(airway))
#
# Or use the ALL leukemia dataset:
# library(ALL)
# data('ALL')
# mat <- exprs(ALL)[1:100, ]
# metadata <- pData(ALL)

# Simulated data with realistic expression patterns
set.seed(42)
n_genes <- 100
n_samples <- 20

pathway_assignment <- factor(rep(c('Immune', 'Signaling', 'Metabolism'),
                                 c(30, 30, 40)),
                             levels = c('Immune', 'Signaling', 'Metabolism'))
pathway_profiles <- rbind(
    Immune = c(rep(-0.7, 10), rep(0.9, 10)),
    Signaling = sin(seq(0, 2 * pi, length.out = n_samples)),
    Metabolism = seq(0.8, -0.8, length.out = n_samples)
)

# Create expression with planted pathway-specific patterns.
mat <- matrix(rnorm(n_genes * n_samples, 0, 0.3), nrow = n_genes,
              dimnames = list(paste0('Gene', 1:n_genes), paste0('Sample', 1:n_samples)))
for (pathway in levels(pathway_assignment)) {
    rows <- pathway_assignment == pathway
    mat[rows, ] <- mat[rows, ] + rep(pathway_profiles[pathway, ], each = sum(rows))
}

metadata <- data.frame(
    sample = colnames(mat),
    condition = rep(c('Control', 'Treatment'), each = 10),
    batch = rep(c('A', 'B'), 10),
    row.names = colnames(mat)
)

gene_info <- data.frame(
    gene = rownames(mat),
    pathway = pathway_assignment,
    log2FC = c(rnorm(30, 1.2, 0.5), rnorm(70, 0, 0.8)),  # DE genes have positive log2FC
    row.names = rownames(mat)
)

# Row scaling makes color and distance represent within-gene expression patterns.
mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0
stopifnot(all(abs(rowMeans(mat_scaled)) < 1e-10))
stopifnot(all(abs(apply(mat_scaled, 1, sd) - 1) < 1e-10))
bounds <- as.numeric(quantile(abs(mat_scaled), 0.99, na.rm = TRUE))
col_fun <- colorRamp2(c(-bounds, 0, bounds), c('#4DBBD5', 'white', '#E64B35'))
pathway_levels <- levels(gene_info$pathway)
pathway_colors <- setNames(c('#F39B7F', '#91D1C2', '#8491B4'), pathway_levels)

ha_col <- HeatmapAnnotation(
    Condition = metadata$condition,
    Batch = metadata$batch,
    col = list(
        Condition = c(Control = '#4DBBD5', Treatment = '#E64B35'),
        Batch = c(A = '#00A087', B = '#3C5488')
    ),
    annotation_name_side = 'left'
)

ha_row <- rowAnnotation(
    Pathway = gene_info$pathway,
    LogFC = anno_barplot(gene_info$log2FC, baseline = 0,
                          gp = gpar(fill = ifelse(gene_info$log2FC > 0, '#E64B35', '#4DBBD5')),
                          width = unit(2, 'cm')),
    col = list(Pathway = pathway_colors)
)
stopifnot(identical(ha_row@anno_list$Pathway@color_mapping@levels, pathway_levels))

ht <- Heatmap(mat_scaled,
              name = 'Z-score',
              col = col_fun,
              top_annotation = ha_col,
              left_annotation = ha_row,
              row_split = gene_info$pathway,
              cluster_rows = TRUE,
              cluster_columns = FALSE,                   # preserve Control -> Treatment input order
              clustering_method_rows = 'ward.D2',
              clustering_method_columns = 'ward.D2',
              column_split = metadata$condition,
              cluster_row_slices = FALSE,
              cluster_column_slices = FALSE,
              show_row_names = FALSE,
              show_column_names = TRUE,
              column_names_rot = 45,
              row_title_rot = 0,
              column_title = 'Expression Heatmap',
              heatmap_legend_param = list(title = 'Z-score', direction = 'horizontal',
                                          at = c(-bounds, 0, bounds)))

pdf('expression_heatmap.pdf', width = 12, height = 10)
ht_drawn <- draw(ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'right')
dev.off()

message('Heatmap saved: expression_heatmap.pdf')
