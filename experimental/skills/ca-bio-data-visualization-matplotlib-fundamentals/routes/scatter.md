# Scatter from CSV

```bash
python scripts/scatter_csv.py data.csv x_column y_column figure.pdf
```

The input must contain the named numeric columns. The wrapper uses explicit Figure/Axes, constrained layout, Type-42 PDF fonts, a rasterized point layer and 89 by 70 mm page size. It overwrites the output. For categories or additional panels, adapt the source example and bind colors to category names rather than input order. Verify fonts and saved page size, and inspect the plot and color accessibility.

Done when the figure has correct axes, size, encoding and readable labels in its intended format.
