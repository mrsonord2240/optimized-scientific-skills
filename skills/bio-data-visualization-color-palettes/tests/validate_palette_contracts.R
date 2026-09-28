# Purpose: regression-test the R palette guidance and both shipped examples.
# Inputs: audit data directory and skill directory.
# Usage: r.sh validate_palette_contracts.R <audit-data-dir> <skill-dir>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
data_dir <- normalizePath(args[[1]], mustWork = TRUE)
skill_dir <- normalizePath(args[[2]], mustWork = TRUE)
script_path <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[[1]])
invisible(parse(file = script_path))
invisible(lapply(list.files(file.path(skill_dir, 'examples'), pattern = '\\.R$',
                            full.names = TRUE), parse))

suppressPackageStartupMessages({
    library(colorspace)
    library(scico)
})

umap <- read.csv(file.path(data_dir, 'synthetic_umap_celltypes.csv'))
cell_levels <- unique(as.character(umap$cell_type))
chromatic_levels <- setdiff(cell_levels, 'Unassigned')
stopifnot(length(chromatic_levels) == 7, 'Unassigned' %in% cell_levels)
okabe_ito <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442',
               '#0072B2', '#D55E00', '#CC79A7')
cell_colors <- c(setNames(okabe_ito, chromatic_levels), Unassigned = '#BBBBBB')
full_mapping <- cell_colors[cell_levels]
subset_levels <- setdiff(cell_levels, chromatic_levels[[1]])
stopifnot(identical(unname(full_mapping[subset_levels]),
                    unname(cell_colors[subset_levels])))

simulated <- list(deutan = deutan(cell_colors), protan = protan(cell_colors),
                  tritan = tritan(cell_colors))
min_distances <- sapply(simulated, function(p) {
    min(dist(coords(as(hex2RGB(p), 'LAB'))))
})
stopifnot(all(is.finite(min_distances)), all(min_distances > 0))

batlow <- scico(10, palette = 'batlow')
L <- coords(as(hex2RGB(batlow), 'LAB'))[, 'L']
stopifnot(all(diff(L) >= 0) || all(diff(L) <= 0))
stopifnot(scico(255, palette = 'vik')[[128]] == '#EBE5E0')
stopifnot(scico(255, palette = 'roma')[[128]] == '#C0E9C2')
stopifnot(scico(255, palette = 'bam')[[128]] == '#F5F0F0')

lfc <- read.csv(file.path(data_dir, 'synthetic_skewed_lfc.csv'))$lfc
vmax <- quantile(abs(lfc), 0.99, na.rm = TRUE)
stopifnot(is.finite(vmax), identical(unname(c(-vmax, vmax)), c(-unname(vmax), unname(vmax))))
custom <- colorRampPalette(c('#0072B2', 'white', '#D55E00'))(101)
stopifnot(custom[[51]] == '#FFFFFF')

old_wd <- getwd()
run_dir <- tempfile('color-palettes-r-')
dir.create(run_dir)
on.exit(setwd(old_wd), add = TRUE)
setwd(run_dir)
example_env <- new.env()
source(file.path(skill_dir, 'examples', 'palette_examples.R'), local = example_env)
stopifnot(file.exists('palette_examples.pdf'), file.info('palette_examples.pdf')$size > 1000)
ggplot2::ggsave('palette_examples.png', example_env$combined, width = 12, height = 10,
                dpi = 120)
stopifnot(file.info('palette_examples.png')$size > 1000)
pdf('palettes_phd.pdf')
source(file.path(skill_dir, 'examples', 'palettes_phd.R'), local = new.env())
invisible(dev.off())
stopifnot(file.exists('palettes_phd.pdf'), file.info('palettes_phd.pdf')$size > 1000)

cat('R palette contracts PASS\n')
cat('render directory:', run_dir, '\n')
cat('rows:', nrow(umap), 'q99:', unname(vmax), '\n')
cat('CVD minimum CIELAB distances:', paste(names(min_distances),
    round(min_distances, 2), collapse = '; '), '\n')
