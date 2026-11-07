args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('usage: audit_palette.R colors.csv output.csv')
library(colorspace)
palette <- read.csv(args[1], stringsAsFactors = FALSE)$color
if (length(palette) < 2 || anyNA(palette)) stop('color column needs at least two colors')
sim <- list(deutan = deutan(palette), protan = protan(palette), tritan = tritan(palette))
distance <- vapply(sim, function(p) min(dist(coords(as(hex2RGB(p), 'LAB')))), numeric(1))
L <- coords(as(hex2RGB(palette), 'LAB'))[, 'L']
monotonic <- all(diff(L) >= 0) || all(diff(L) <= 0)
write.csv(data.frame(simulation = names(distance), min_cielab_distance = unname(distance),
                     luminance_monotonic = monotonic), args[2], row.names = FALSE)
