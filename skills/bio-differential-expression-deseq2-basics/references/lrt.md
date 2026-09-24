# Wald vs LRT in DESeq2

## Wald vs LRT

**Goal:** Choose Wald for single-coefficient hypotheses, LRT for joint hypotheses involving more than 1 df.

**Approach:** Wald is default. LRT is mandatory for multi-level factors tested as "any change", interactions with >1 df, and ANOVA-style omnibus tests. With LRT, the reported LFC is for the LAST coefficient in `resultsNames(dds)` -- not the omnibus effect.

```r
dds <- DESeq(dds, test = 'LRT', reduced = ~ batch)
res_lrt <- results(dds)
# res_lrt$padj is the LRT joint p-value (correct).
# res_lrt$log2FoldChange is for the last coefficient in resultsNames -- NOT an omnibus summary.
```

When reporting LRT results, name what the LFC actually represents or extract specific Wald coefficients per level for the effect-size table.

### LRT reports the wrong LFC

**Trigger:** Multi-level factor analyzed with `DESeq(dds, test='LRT', reduced=~1)`; user reports the `log2FoldChange` column as "the effect".

**Mechanism:** The LRT p-value is the omnibus test of "any difference among levels". The LFC reported by `results()` after LRT is for the LAST coefficient in `resultsNames(dds)`, which is one specific level-vs-reference comparison.

**Symptom:** A 4-level factor produces a single LFC value per gene; reviewer asks "the effect of which condition?"

**Fix:** Treat LRT padj as a screen for "any change". For effect sizes, extract per-level Wald coefficients individually via `results(dds, name='<specific coefficient>')` for each non-reference level.
