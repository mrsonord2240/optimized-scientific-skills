# Customize a plot

```bash
Rscript scripts/plot_csv.R box data.csv output x_column y_column
```

Use this file-based plot as a starting point; `scripts/publication_figures.R` provides `theme_publication`, `okabe_ito`, box, PCA, volcano, save and patchwork helpers. Adapt its functions for facets, rich labels or other geoms using [recipes](../references/geoms-scales-facets.md); these alternatives do not have a standalone CLI. Put constants outside `aes()`, use `.data[[name]]` for string-based mapping and `linewidth` for lines. Fixed facet scales support comparison; justify free scales. PDF export uses Cairo and explicit millimeters; raster output uses 300 dpi. For dense points consider rasterization while retaining vector axes. See [failure modes](../references/failure-modes.md).

Done when rendered output has been opened and checked for clipping, labels, colors, fonts and intended dimensions.
