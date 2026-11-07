# MA diagnostic

```bash
python scripts/ma_csv.py results.csv ma.pdf
```

CSV requires `baseMean`, `log2FoldChange`, `padj`; this wrapper requires finite positive means and finite effects, with missing adjusted p-values permitted. The source helper plots log10 mean against LFC, colors adjusted p < 0.05, and rasterizes points; it labels the effect as shrunken, so supply documented shrunken estimates. Zero mean rows require a separately justified handling rule before plotting. Output can overwrite an existing file. An apparent pattern can suggest a diagnostic question, but it does not by itself establish a normalization or batch problem.

Done when effect scale, mean definition, excluded data and the displayed diagnostic are checked.
