# Fold changes for ranking or a figure

Raw `logFC` is the estimate to rank by and to pool in a meta-analysis. For shrunken fold changes add `shrink=ashr`:

```bash
Rscript scripts/limma_de.R proteinGroups.txt samples.csv Treatment-Control results.csv normalize=median shrink=ashr
```

The output gains `logFC_shrunk` and `lfsr`.

- GSEA or other ranked lists: use the raw `logFC` (or the `t` statistic), not `logFC_shrunk`.
- Use `logFC_shrunk` to say which proteins truly changed and by how much. Report it beside the raw `logFC`, not instead of it.
- Do not set shrunken fold changes to 0 at `adj.P.Val > 0.05`.
- Significance stays `adj.P.Val`.
- Needs the ashr package; the script stops if it is missing.

Done when the ranked list or figure uses the column named above.
