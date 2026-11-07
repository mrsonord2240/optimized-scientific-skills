# Demonstrate hierarchical edge bundling on ggraph's documented flare hierarchy.
# Inputs: output PNG path. The bundled flare hierarchy supplies vertices and imports.
# Usage: Rscript scripts/edge_bundling.R out/edge-bundling.png

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/edge_bundling.R output.png")
}
output <- normalizePath(args[[1]], mustWork = FALSE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

edges <- flare$edges
vertices <- flare$vertices
imports <- flare$imports
graph <- graph_from_data_frame(edges, vertices = vertices)
from_index <- match(imports$from, vertices$name)
to_index <- match(imports$to, vertices$name)
if (anyNA(from_index) || anyNA(to_index)) stop("Flare imports do not match vertices")

plot <- ggraph(graph, layout = "dendrogram", circular = TRUE) +
  geom_conn_bundle(
    data = get_con(from = from_index, to = to_index),
    alpha = 0.4,
    tension = 0.8,
    edge_colour = "grey60"
  ) +
  geom_node_point(size = 0.5) +
  theme_void()
ggsave(output, plot, width = 7, height = 7, dpi = 160, bg = "white")
if (file.size(output) == 0) stop("Empty output: ", output)
message(length(from_index), " bundled connections; ", file.size(output), " bytes")
