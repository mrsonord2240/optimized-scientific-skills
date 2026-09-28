# Render reproducible ggraph layouts from GraphML.
# Inputs: GraphML path and output directory.
# Usage: Rscript scripts/ggraph_layouts.R network.graphml out/r-layouts

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript scripts/ggraph_layouts.R network.graphml output_dir")
}
input <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
graph <- read_graph(input, format = "graphml")
if (vcount(graph) == 0) stop("The input graph is empty")
V(graph)$node_degree <- degree(graph)

for (layout_name in c("fr", "kk", "circle")) {
  set.seed(42)
  plot <- ggraph(graph, layout = layout_name) +
    geom_edge_link(alpha = 0.3, colour = "grey60") +
    geom_node_point(aes(size = node_degree), colour = "#4DBBD5") +
    scale_size_continuous(range = c(2, 5), guide = "none") +
    theme_graph()
  output <- file.path(output_dir, paste0(layout_name, ".png"))
  ggsave(output, plot, width = 7, height = 6, dpi = 160)
  if (file.size(output) == 0) stop("Empty output: ", output)
  message(layout_name, ": ", vcount(graph), " nodes; ", file.size(output), " bytes")
}
