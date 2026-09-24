# coloc.susie Integration

## Coloc.susie Integration

**Goal:** Test colocalization between two traits using credible sets, not single SNPs.

**Approach:** Fit susie_rss separately per trait; pass both `susie` objects to `coloc.susie`; per-credible-set colocalization probabilities are returned.

**Precondition:** `coloc.susie` matches SNPs between the two fits' `lbf_variable` matrices via `intersect(colnames(...))`. If `z`/`R` are unnamed, that intersect is empty and `coloc.susie` fails with a cryptic, unrelated `data.table` error (`Check that is.data.table(DT) == TRUE ... := is defined for use in j`) instead of a clear message. Name `z1`/`z2` and set matching `ld_matrix` dimnames to the same SNP IDs before fitting.

```r
library(coloc)

names(z1) <- names(z2) <- colnames(ld_matrix) <- rownames(ld_matrix) <- snp_ids

fit_trait1 <- susie_rss(z = z1, R = ld_matrix, n = N1, L = 10)
fit_trait2 <- susie_rss(z = z2, R = ld_matrix, n = N2, L = 10)

stopifnot("no shared SNP names between the two fits: name z and dimnames(R) first (see Precondition)" =
          length(intersect(colnames(fit_trait1$lbf_variable), colnames(fit_trait2$lbf_variable))) > 0)
coloc_res <- coloc.susie(fit_trait1, fit_trait2)
# coloc_res$summary: per-credible-set PP.H4 (shared causal probability)
print(coloc_res$summary)
```

PP.H4 > 0.8 per credible set is the conventional shared-causal threshold; weaker thresholds suggest distinct or conditional signals. See causal-genomics/colocalization-analysis.
