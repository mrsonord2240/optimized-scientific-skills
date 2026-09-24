# Allele Harmonization

### Allele Harmonization with the LD Reference

**Trigger:** Effect allele in GWAS sumstats differs from the coding/A1 allele in the LD reference panel; or palindromic SNPs (A/T, C/G) carried without strand resolution.

**Mechanism:** susie_rss treats `z` and `R` as defined on the same allele coding. If the effect allele is swapped relative to the LD-reference A1, the sign of z is wrong and the LD row/column for that SNP is implicitly flipped. SNPs matching by rsID can silently swap alleles between sumstats and reference, breaking the `z' R z` consistency the model relies on.

**Symptom:** `estimate_s_rss` lambda inflated despite ancestry-matched panel; `kriging_rss` flags many SNPs with `|z_obs - z_exp| > 3` clustered at SNPs where reference A1 != GWAS effect allele; credible sets pick up tag-only SNPs anti-correlated with the lead.

**Fix:** Harmonize before fitting:

```r
harmonize_z_to_ref <- function(z, gwas_a1, gwas_a2, ref_a1, ref_a2) {
    palindromic <- (gwas_a1 == 'A' & gwas_a2 == 'T') | (gwas_a1 == 'T' & gwas_a2 == 'A') |
                   (gwas_a1 == 'C' & gwas_a2 == 'G') | (gwas_a1 == 'G' & gwas_a2 == 'C')
    flip <- (gwas_a1 == ref_a2) & (gwas_a2 == ref_a1)
    z[flip] <- -z[flip]
    drop <- palindromic | !((gwas_a1 == ref_a1 & gwas_a2 == ref_a2) | flip)
    list(z = z, keep = !drop)
}
```

Drop palindromic SNPs at MAF > 0.42 (ambiguous strand); or resolve via external strand info (TopMed, 1000G strand files). `TwoSampleMR::harmonise_data()` offers an alternative implementation. See causal-genomics/colocalization-analysis for an equivalent harmonize helper used downstream.
