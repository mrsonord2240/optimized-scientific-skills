# Plot a saved result

```bash
Rscript scripts/plot_result.R result.rds dotplot figure.pdf 20
```

Input is a nonempty clusterProfiler enrichment result saved as RDS. The wrapper selects up to 20 terms and reports the object row count; this is the object's retained result count, not necessarily all tested terms. For a GSEA result, `Rscript scripts/plot_result.R result.rds nes figure.pdf 20` maps signed NES to x. `Rscript scripts/plot_result.R result.rds running figure.pdf 1` draws one GSEA set by index. Verify object class and installed enrichplot API; output may be overwritten. A signed NES is direction relative to the ranked contrast, not proof of activation.

Done when class, ranking direction, term selection, counts and rendered figure are reviewed.
