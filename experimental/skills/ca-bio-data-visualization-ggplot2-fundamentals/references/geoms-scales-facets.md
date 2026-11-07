# ggplot2 geoms, aesthetics, scales, facets

Lookup for ggplot2 4.x (tested 4.0.3; patterns hold on 3.5+).

## Common geoms

```r
geom_point(alpha = 0.7, size = 1)                      # large N: wrap in ggrastr::rasterise(geom_point(...), dpi = 300); ggplot2 has no rasterize argument
geom_line(linewidth = 0.5)                             # linewidth replaces size for lines (ggplot2 3.4+)
geom_col()                                              # bar with y values (use this; geom_bar(stat='identity') is older)
geom_bar()                                              # bar with counts
geom_boxplot(outlier.shape = NA)                       # always suppress when overlaying jitter
geom_violin(bw = 'SJ', trim = FALSE)                   # Sheather-Jones bandwidth; show full tails
geom_histogram(bins = 30)                              # bins NOT binwidth for control
geom_density(alpha = 0.5)
geom_tile(aes(fill = z))                               # heatmap building block
geom_text(aes(label = label), check_overlap = TRUE)
geom_text_repel(aes(label = label), max.overlaps = Inf)   # ggrepel; max.overlaps = Inf prevents silent label drops
```

## Aesthetic Mappings

```r
aes(x, y, color, fill, shape, size, alpha, linetype, linewidth, group)

# Color vs fill: color = stroke; fill = interior (boxplot, bar, area, polygon)
# Use both when needed: geom_point(aes(color = group, fill = group), shape = 21)
```

**Constant inside vs mapping inside aes** is a common confusion:
```r
geom_point(color = 'red')             # constant: every point red
geom_point(aes(color = group))        # mapping: color varies with group
```

## Scales

```r
# Continuous
scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, 2),
                    labels = scales::label_number(scale = 1e-6, suffix = 'M'))
scale_y_log10()
scale_y_continuous(transform = 'sqrt')

# Discrete
scale_x_discrete(limits = c('Control', 'Treatment', 'Vehicle'))   # explicit order
scale_color_manual(values = c(Control = '#0072B2', Treatment = '#D55E00'))

# Colormap (sequential, diverging, cyclic) -- see color-palettes
scale_color_viridis_c(option = 'viridis')
scale_color_scico(palette = 'batlow')                              # Crameri
scale_fill_gradient2(low = '#0072B2', mid = 'white', high = '#D55E00', midpoint = 0)

# Date / time
scale_x_date(date_breaks = '1 year', date_labels = '%Y')
```

## Facets

```r
facet_wrap(~ var, ncol = 3, scales = 'free_y')
facet_grid(rows = vars(condition), cols = vars(timepoint), scales = 'free_x')
facet_grid(condition ~ timepoint)                                  # formula syntax
```

`scales = 'free_y'` lets each panel have its own y-range — appropriate when biological scales differ across facets. `scales = 'fixed'` (default) is the right choice when comparing across panels.
