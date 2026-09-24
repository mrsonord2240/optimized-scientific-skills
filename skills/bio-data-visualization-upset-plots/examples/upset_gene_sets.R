# Purpose: deterministic UpSetR examples with identifier preflight and Cairo PDF export.
# Usage: Rscript upset_gene_sets.R (run in a directory where PDFs may be written).
library(UpSetR)

prepare_sets <- function(sets) {
  stopifnot(is.list(sets), !is.null(names(sets)), all(nzchar(names(sets))))
  cleaned <- lapply(names(sets), function(nm) {
    x <- trimws(as.character(sets[[nm]]))
    bad <- is.na(x) | !nzchar(x)
    if (any(bad)) warning(sprintf("%s: removed %d NA/blank identifier(s)", nm, sum(bad)))
    x <- x[!bad]
    old_n <- length(x); x <- unique(x)
    if (length(x) != old_n) warning(sprintf("%s: removed %d duplicate identifier(s)", nm, old_n - length(x)))
    x
  })
  names(cleaned) <- names(sets)
  empty <- names(cleaned)[lengths(cleaned) == 0L]
  if (length(empty)) stop("empty set(s) after preflight: ", paste(empty, collapse = ", "))
  message("set sizes: ", paste(sprintf("%s=%d", names(cleaned), lengths(cleaned)), collapse = "; "),
          " | union=", length(unique(unlist(cleaned, use.names = FALSE))))
  cleaned
}

set.seed(42)
all_genes <- paste0("Gene", 1:500)
gene_sets <- list(
  Treatment_vs_Control = sample(all_genes, 150), Timepoint_6h = sample(all_genes, 120),
  Timepoint_24h = sample(all_genes, 180), Drug_A = sample(all_genes, 100),
  Drug_B = sample(all_genes, 90), Combined_Treatment = sample(all_genes, 200)
)
core_genes <- sample(all_genes, 30)
for (i in 1:4) gene_sets[[i]] <- c(gene_sets[[i]], core_genes)
gene_sets <- prepare_sets(gene_sets)
upset_data <- UpSetR::fromList(gene_sets)
colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")

cairo_pdf("upset_basic.pdf", width = 12, height = 7)
UpSetR::upset(upset_data, nsets = length(gene_sets), order.by = "freq",
              mainbar.y.label = "Genes in exclusive intersection", sets.x.label = "Total genes per set")
dev.off()

cairo_pdf("upset_customized.pdf", width = 12, height = 8)
UpSetR::upset(upset_data, nsets = length(gene_sets), nintersects = 30, order.by = "freq", decreasing = TRUE,
              mb.ratio = c(0.55, 0.45), point.size = 3.5, line.size = 1.2, sets.bar.color = colors,
              main.bar.color = "#7E6148", matrix.color = "#7E6148",
              text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1.2), set_size.show = TRUE)
dev.off()

cairo_pdf("upset_queries.pdf", width = 12, height = 8)
UpSetR::upset(upset_data, nsets = length(gene_sets), order.by = "freq",
  queries = list(
    list(query = intersects, params = list("Combined_Treatment"), color = "#E64B35", active = TRUE,
         query.name = "Combined treatment"),
    list(query = intersects, params = list("Timepoint_6h", "Timepoint_24h"), color = "#4DBBD5", active = TRUE,
         query.name = "Both timepoints")), query.legend = "bottom")
dev.off()
message("Saved upset_basic.pdf, upset_customized.pdf, upset_queries.pdf")
