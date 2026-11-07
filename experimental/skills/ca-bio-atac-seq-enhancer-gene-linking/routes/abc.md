# ABC predictions

```bash
ABC_REPO=/path/to/ABC ACCESS_BAM=access.bam OUTDIR=abc_out bash scripts/run_abc.sh
```

Set ACCESS_TYPE to ATAC or DHS; optional H3K27AC_BAM and HIC_TYPE/HIC_FILE select inputs. Default HIC_TYPE=none uses powerlaw; report fallback or proxy. Use matching genome, gene and TSS annotations, blacklist and Hi-C. The wrapper includes TSS candidates, quantile normalization and input-calibrated threshold lookup; 0.02 is its fallback when no row matches. See `references/setup.md` and `references/method-reference.md` for configuration. For scATAC, pseudobulk by cluster with adequate unique fragments.

Done when thresholded links and input, contact source and thresholds are reported.
