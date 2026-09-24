# Deterministic regression test for ComplexHeatmap 2.22.0:
# an explicit OLO dendrogram renders with a Pathway row annotation, whereas a
# categorical row_split cannot be supplied alongside that dendrogram.

library(ComplexHeatmap)
library(circlize)
library(seriation)

mat <- matrix(c(
  0, 0, 1, 1,
  0, 1, 1, 2,
  8, 8, 9, 9,
  8, 9, 9, 10,
  3, 3, 4, 4,
  3, 4, 4, 5
), nrow = 6, byrow = TRUE,
dimnames = list(paste0("gene", 1:6), paste0("sample", 1:4)))

declared_pathway_levels <- c("Metabolism", "Signaling", "Immune")
pathway <- factor(c("Metabolism", "Metabolism", "Signaling",
                    "Signaling", "Immune", "Immune"),
                  levels = declared_pathway_levels)
pathway_levels <- unique(as.character(pathway))

d_rows <- dist(mat, method = "euclidean")
hc_rows <- hclust(d_rows, method = "ward.D2")
olo_rows <- seriate(d_rows, method = "OLO", control = list(hclust = hc_rows))
dend_rows <- as.dendrogram(olo_rows[[1]])
expected_order <- order.dendrogram(dend_rows)

pathway_colors <- setNames(grDevices::hcl.colors(length(pathway_levels), palette = "Dark 3"),
                           pathway_levels)
ha_row <- rowAnnotation(
  Pathway = pathway,
  col = list(Pathway = pathway_colors)
)
annotation_pathway <- ha_row@anno_list$Pathway@fun@var_env$value
stopifnot(identical(ha_row@anno_list$Pathway@color_mapping@levels, pathway_levels))
stopifnot(identical(sort(unique(annotation_pathway)), sort(pathway_levels)))

# ComplexHeatmap 2.22.0 rejects this incompatible combination. Keep this
# assertion so a future package change is explicit rather than silently assumed.
incompatible_error <- tryCatch(
  {
    Heatmap(
      mat,
      name = "value",
      cluster_rows = dend_rows,
      cluster_columns = FALSE,
      row_split = pathway
    )
    NULL
  },
  error = conditionMessage
)
stopifnot(is.character(incompatible_error))
stopifnot(grepl("row_split.*single", incompatible_error))

ht <- Heatmap(
  mat,
  name = "value",
  col = colorRamp2(c(0, 5, 10), c("#0072B2", "white", "#D55E00")),
  cluster_rows = dend_rows,
  cluster_columns = FALSE,
  left_annotation = ha_row,
  show_row_names = TRUE
)

pdf("heatmap_olo_row_annotation_test.pdf", width = 5, height = 4)
ht_drawn <- draw(ht)
dev.off()

stopifnot(file.exists("heatmap_olo_row_annotation_test.pdf"))
stopifnot(file.info("heatmap_olo_row_annotation_test.pdf")$size > 0)
stopifnot(identical(as.integer(row_order(ht_drawn)), as.integer(expected_order)))

message("PASS: explicit OLO dendrogram rendered with Pathway annotation and retained OLO order")
