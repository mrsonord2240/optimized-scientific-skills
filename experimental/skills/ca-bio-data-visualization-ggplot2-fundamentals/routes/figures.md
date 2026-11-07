# File-based figures

```bash
Rscript scripts/plot_csv.R box data.csv output x_column y_column
```

The box route reads CSV and uses source `create_boxplot` with jitter and suppressed duplicate outliers. For PCA use `Rscript scripts/plot_csv.R pca data.csv output group_column`: the file needs `PC1`, `PC2`, and `var_explained`, with the first two rows giving PC1 and PC2 variance in percent, not fractions. For volcano use `Rscript scripts/plot_csv.R volcano data.csv output`: columns `log2FoldChange`, `padj`, and `gene` are required. Source thresholds are adjusted p < 0.05 and absolute log2 fold change > 1; labels use 10 short or 3 long names by default. Output is `output.pdf` and `output.png`, overwritten if present. The source helper caps its palette at eight groups. Count displayed labels and inspect overlaps: published observations are dataset- and layout-specific, not a guarantee.

Done when data, dimensions, labels, palette and exported figures have been checked.
