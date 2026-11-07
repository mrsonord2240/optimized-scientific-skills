# Apply a palette

```bash
Rscript scripts/plot_palette.R data.csv x y value sequential figure.pdf
```

Supply a CSV and column names followed by `sequential`, `diverging`, `cyclic`, or `categorical`, and an output path. The script plots tiles with batlow, vik, or romaO; diverging bounds use the 99th percentile of absolute finite values, centered on zero and squished. Categorical plots require at most seven nonmissing groups; sorted names receive the source Okabe-Ito colors. Preserve a separately named mapping across panels and use a separate neutral encoding for missing values. For a nonzero reference or custom bounds, adapt the scale rather than using the zero-centered wrapper. Check output location before use; it is overwritten. Source demonstrations are in `examples/` and use simulated data.

Done when the figure, mapping, bounds, clipping, missingness and labels are reviewed; run the separate [audit](audit.md).
