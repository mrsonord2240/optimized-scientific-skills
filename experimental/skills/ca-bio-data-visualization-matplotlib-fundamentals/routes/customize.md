# Charts and export

```bash
python scripts/matplotlib_phd.py
```

This source demonstration uses simulated data and writes five fixed PDF names to the current directory; run only in a new directory and adapt data before substantive use. It covers grids, heatmaps, seaborn Axes and seaborn.objects. Its objects legend uses a private seaborn attribute and a fixed page reserve that can overlap long legends; prefer the documented public `Plot.on(fig)` approach and inspect placement. Other chart and axis recipes require adaptation to input data; see [recipes](../references/chart-recipes.md). Set `pdf.fonttype=42`, use constrained layout and avoid tight bounding-box recropping when exact page size matters. Figure sizes are inches, so divide millimeters by 25.4. Axes-level seaborn functions return Axes; figure-level functions return FacetGrid. See [failure modes](../references/failure-modes.md).

Done when the actual export is inspected for page size, font embedding, labels, legends and accessibility.
