---
name: bio-data-visualization-ggplot2-fundamentals
description: Use when producing static publication-quality figures in R with ggplot2 for papers, presentations, or reports.
tool_type: r
primary_tool: ggplot2
goal_approach_exempt: true
license: MIT
author: GPTomics
category: Data Analysis
---

## Version Compatibility

Reference examples tested with: ggplot2 3.5+ (re-checked on ggplot2 4.0.3), scales 1.3+ (1.4.0), ggrepel 0.9.5+ (0.9.8), ggtext 0.1.2+ (0.2.0), viridis 0.6+ (0.6.5), scico 1.5+, patchwork 1.2+ (1.3.2), ggrastr 1.0+ (1.0.2).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws an error or deprecation warning, introspect the installed package and adapt the example to match the actual API rather than retrying. ggplot2 4.0 objects are S7 (`class(p)` is no longer just `gg`/`ggplot`); extension packages built against 3.5 may fail on 4.x.

# ggplot2 Fundamentals

**"Build a publication figure in R"** -> Express the figure as **data + aesthetic mappings + one or more geometries + scales + facets + theme**. The grammar of graphics (Wilkinson 2005; Wickham 2010 *J Comput Graph Stat* 19:3) makes each visual element separately addressable — change scales without rewriting geoms; swap geom_point for geom_violin without touching aesthetics.

- R: `ggplot(data, aes(x, y)) + geom_point() + scale_color_manual(...) + theme_classic()`
- Programmatic: `aes(x = .data[[var]])` for tidy-eval; `!!sym(var)` for dplyr-style symbol evaluation

## The Three Modern Defaults

1. **theme_classic() + remove panel grid + Okabe-Ito palette** as the publication baseline. `theme_minimal` adds light gridlines; `theme_bw` adds a panel border; both work but `theme_classic` is the cleanest for journals.

2. **cairo_pdf for every PDF export** — it embeds TrueType fonts (searchable PDFs); the default `pdf()` device leaves Helvetica/Symbol unembedded, which journals reject. See Saving.

3. **Tidy evaluation for programmatic aes** — `aes(x = .data[[var]])` is the modern idiom (ggplot2 3.0+); `aes_string(x = var)` is deprecated.

## Grammar in Layers

```r
library(ggplot2)

# data + aes + geom is the minimum
ggplot(df, aes(x = condition, y = expression)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = condition), width = 0.2, alpha = 0.5) +
    # scales (log10 axis with plain-number tick labels)
    scale_y_continuous(transform = 'log10', breaks = scales::breaks_log(), labels = scales::label_number()) +
    scale_color_manual(values = c('#0072B2', '#D55E00')) +
    # labels
    labs(x = NULL, y = 'Expression',
         title = 'Gene X across conditions',
         caption = 'Source: ...') +
    # facets
    facet_wrap(~ tissue, ncol = 3, scales = 'free_y') +
    # theme
    theme_classic(base_size = 10) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = 'bold'))
```

Geoms, aesthetic mappings (constant vs `aes()` mapping, color vs fill), scales, and facets: `references/geoms-scales-facets.md`.

## Theme

The baseline (theme_classic, no grid, black axes and ticks, bold strips) and the Okabe-Ito palette live in `scripts/publication_figures.R` as `theme_publication()` and `okabe_ito`; source it and use them rather than redefining.

```r
source('scripts/publication_figures.R')
ggplot(df, aes(x, y, color = group)) + geom_point() +
    scale_color_manual(values = okabe_ito[1:2]) + theme_publication()
```

## Programmatic Plots (Tidy Evaluation)

```r
# Pass variable name as a string
plot_var <- function(df, x_var, y_var) {
    ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]])) +
        geom_point()
}
plot_var(df, 'PC1', 'PC2')

# Alternative: bare names via embracing
plot_var2 <- function(df, x_var, y_var) {
    ggplot(df, aes(x = {{ x_var }}, y = {{ y_var }})) +
        geom_point()
}
plot_var2(df, PC1, PC2)   # bare names only; plot_var2(df, 'PC1', 'PC2') maps a constant string
```

## Labels with ggtext (rich-text)

```r
library(ggtext)
ggplot(df, aes(x, y)) + geom_point() +
    labs(x = 'log<sub>2</sub> fold change',
         y = '\u2212log<sub>10</sub>(*p*)') +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())
```

ggtext renders inline HTML / Markdown in titles, captions, axis labels — much better than `expression(...)` for italics + subscripts + special characters.

## Saving — TrueType Embedding

```r
# cairo_pdf for TrueType embedded; portable across systems
ggsave('figure.pdf', plot = p,
       width = 89, height = 70, units = 'mm',
       device = cairo_pdf)

# Vector + raster mix via ggrastr (for large scatter)
library(ggrastr)
ggplot(df, aes(x, y)) +
    rasterise(geom_point(alpha = 0.5), dpi = 300) +
    theme_publication()
ggsave('out.pdf', device = cairo_pdf)

# PNG for raster
ggsave('figure.png', p, width = 89, height = 70, units = 'mm', dpi = 300)

# TIFF for some journals
ggsave('figure.tiff', p, width = 89, height = 70, units = 'mm', dpi = 300,
       compression = 'lzw')
```

Always set `units = 'mm'` explicitly (default is inches). Nature single column = 89 mm; double column = 183 mm.

## Guardrails

- `geom_text_repel(..., max.overlaps = Inf)` — the default (10) silently drops labels that overlap more than 10 neighbours; count the labels actually drawn.
- `geom_boxplot(outlier.shape = NA)` when overlaying jitter.
- `linewidth` (not `size`) for lines; `facet_wrap(scales = 'fixed')` when panels are compared.
- Constants go outside `aes()`; `aes(color = 'red')` creates a one-level category.

Trigger/mechanism/symptom/fix for each: `references/failure-modes.md`.

## Resources

- `scripts/publication_figures.R` — reusable helpers: `theme_publication()`, `okabe_ito`, `create_volcano()` (y = -log10(padj); needs `log2FoldChange`, `padj`, `gene` or a `label_col`; labels the 10 smallest padj when they are short symbols, 3 when long IDs, or `top_n`; no label-label overlaps in the airway DESeq2 results (the only dataset measured) for symbols in a standalone figure down to 89 mm wide and Ensembl IDs down to 120 mm wide (heights not recorded), and in `create_multi_panel()` composites at 183 x 120 mm and 183 x 150 mm; a 183 x 90 mm composite had 1 overlapping pair, a re-ranked Ensembl set had 3 at 183 x 120 mm, and a 120 mm wide composite collides with either; the counts cover label boxes only, threshold lines and points still cross labels in rows counted clean, so keep composites at 183 mm wide and 120 mm or taller and open the figure rather than rely on this list), `create_boxplot()`, `create_pca_plot()` (needs a `var_explained` column in percent, 0-100; fractions are read as percent and warn), `save_publication_figure()` (mm, cairo_pdf), `create_multi_panel()` (patchwork). Source it, then call the functions.
- `usage-guide.md` — prerequisites and example prompts.

## References

- Wickham H. 2016. *ggplot2: Elegant Graphics for Data Analysis* (2nd ed). Springer.
- Wickham H. 2010. A layered grammar of graphics. *J Comput Graph Stat* 19(1):3-28.
- Wilkinson L. 2005. *The Grammar of Graphics* (2nd ed). Springer.

## Related Skills

- data-visualization/color-palettes - Scale_color/_fill palette selection
- data-visualization/multipanel-figures - patchwork composition
- data-visualization/distribution-plots - Box / violin / raincloud geoms
- data-visualization/volcano-and-ma-plots - ggplot2 volcano with ggrepel
- data-visualization/heatmaps-clustering - ComplexHeatmap and ggplot2 geom_tile
