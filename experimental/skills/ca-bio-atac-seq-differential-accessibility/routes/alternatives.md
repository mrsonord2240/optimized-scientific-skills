# Alternative designs

```bash
Rscript scripts/diff_accessibility.R samples.csv results lib 0.05 1 --method=edger --case=treated --control=control --txdb=none
```

For unstable peak boundaries consider window-based csaw; for existing matrices use edgeR or DESeq2 directly. The source supplies partial recipes rather than a complete executable window pipeline; see `references/alternatives.md` and the relevant `bio-differential-expression-deseq2-basics` skill. Time-course LRT, low-replicate permutation and loop-anchored aggregation are also described there; treat the loop snippet as a pseudopattern, not a validated test. A single sample per group cannot support replicate-based inference.

Done when the requested method is fitted to the actual data and model, replicates and limitations are reported.
