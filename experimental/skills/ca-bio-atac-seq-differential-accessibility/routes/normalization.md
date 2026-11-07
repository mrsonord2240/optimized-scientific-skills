# Normalization or hidden batch

```bash
Rscript scripts/diff_accessibility.R samples.csv results lib 0.05 1 --case=treated --control=control --sva=1 --txdb=none
```

SVA fits DESeq2 directly on DiffBind counts, ignores `--method` and `--design`, caps surrogate count by sample size, omits the heatmap and treats fold threshold as a post-filter. Its results are not directly comparable to the default threshold test. No branch of the script applies a blacklist. For full-library versus reads-in-peaks, native scaling, spike-in and RUV alternatives read `references/normalization.md`; full-library scaling alone cannot establish absolute global changes.

Done when assumptions, surrogate count, filtering and contrasts are documented with the results.
